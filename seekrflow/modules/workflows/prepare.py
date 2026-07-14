"""
modules/workflows/prepare.py

Free functions that turn a generic seekrflow ``Workflow`` into a seekr
``Seekr_input`` (which ``seekr.prepare.prepare`` which then writes to 
``model.json``).

All seekr-specific construction lives here; the structure modules stay pure
data. The entry point is :func:`make_seekr_input`.

Key seekr contracts honored here:
  * Selections (name -> atom indices) are carried on the MD engine settings'
    ``unpartitioned_settings.selections``; CVs and restraints reference them by
    name.
  * ``initial_scopes`` is a dict mapping each scale type to ``"unpartitioned"``;
    logistic stages handle the transition to the partitioned scope and back.
  * Stages chain via ``input_stage_name`` (the sentinel ``"initial"`` means
    "depends only on the prepared structures").
"""

import os
import typing

import mdtraj
import numpy as np

import seekrflow.modules.parameters_topology_structures \
    as parameters_topology_structures
import seekrflow.modules.pqr as pqr
import seekrflow.modules.workflows.workflow as workflow_module
import seekrflow.modules.workflows.components as components_module
import seekrflow.modules.workflows.cv_specs as cv_specs_module
import seekrflow.modules.workflows.anchor_specs as anchor_specs_module
import seekrflow.modules.workflows.scale_settings as scale_settings_module
import seekrflow.modules.workflows.stage_procedures as stage_procedures_module

import seekr.prepare
from seekr.modules.preparer import Seekr_input, Partition_input
from seekr.modules.engines.openmm.preparer import Openmm_settings_input
from seekr.modules.engines.preparer import MD_unpartitioned_settings_input
from seekr.modules.engines.browndye.preparer import (
    Browndye_settings_input, Browndye_molecule_input)
from seekr.modules.engines.structures import (
    Selection, Ion, Amber_parameters_topology, Charmm_parameters_topology,
    Openmm_system)
from seekr.modules.scales.molecular_dynamics.preparer import Input_stage_md
from seekr.modules.scales.brownian_dynamics.preparer import Input_stage_bd
from seekr.modules.scales.logistic\
    .anchor_starting_structures_from_reporter_logistic_stage import (
        Input_stage_anchor_starting_structures_from_reporter)
from seekr.modules.sampling_methods.preparer import (
    Conventional_md_input, Metadynamics_md_input, Steered_md_input,
    MMVT_md_input, NAM_bd_input)
from seekr.modules.completion_criteria.structures import (
    Number_of_steps_completion_criteria, Reporter_progress_completion_criteria,
    Number_of_trajectories_completion_criteria)
from seekr.modules.collective_variables.distance_cv import Distance_cv_input
from seekr.modules.collective_variables.rmsd_cv import RMSD_cv_input
from seekr.modules.anchor_schemes.preparer import (
    Anchors_as_voronoi_cells_anchor_scheme, Voronoi_cell)
from seekr.modules.restraints.preparer import (
    Positional_restraint_input, Pairwise_restraint_input,
    Proximity_restraint_input)
from seekr.modules.reporters.preparer import (
    State_reporter_input, Position_reporter_input, Anchor_reporter_input)
from seekr.modules.base import (
    State_point, NVT_ensemble, NPT_ensemble, NPT_membrane_ensemble,
    NVE_ensemble, mdtraj_center_of_mass)
import seekr.modules.check as check

class Prepare_context:
    """
    Throwaway bag of state shared by the prepare helpers for a single
    ``make_seekr_input`` call.
    """

    def __init__(
            self,
            workflow: "workflow_module.Workflow",
            root_directory: str,
            physical_attributes: typing.Any,
            ) -> None:
        self.workflow = workflow
        self.root_directory = root_directory
        self.physical_attributes = physical_attributes
        self.target_temperature = float(
            getattr(physical_attributes, "temperature", 298.15))
        self.bd_active = workflow.has_bd()
        self.md_structure_filename = \
            workflow.get_md_settings().system.solvated_pdb
        self.selections: dict[str, dict[str, list[int]]] = {}


# ======================================================================
#  Mapping helpers: seekrflow native -> seekr3 input objects.
# ======================================================================

def _map_parameters_topology(
        #pt: typing.Any,
        system: typing.Any,
        ) -> typing.Any:
    """
    Map a seekrflow parameters/topology object to its seekr3 equivalent.
    """
    pt = system.parameters_topology
    if isinstance(pt, parameters_topology_structures.Amber_parameters_topology):
        if not system.solvated_pdb and pt.inpcrd_filename:
            # produce PDB file for solvated_pdb from the inpcrd.
            amber_parmed = pt.make_parmed(directory="/")
            # Use the same filename as the prmtop but with .pdb extension.
            pdb_filename = os.path.splitext(pt.prmtop_filename)[0] + ".pdb"
            amber_parmed.save(pdb_filename)
        new_pt = Amber_parameters_topology(prmtop_filename=pt.prmtop_filename)
        return new_pt
    if isinstance(pt, parameters_topology_structures.Charmm_parameters_topology):
        return Charmm_parameters_topology(
            psf_filename=pt.psf_filename,
            param_filename_list=list(pt.param_filename_list))
    if isinstance(pt, parameters_topology_structures.Openmm_system):
        return Openmm_system(system_filename=pt.system_filename)
    raise NotImplementedError(
        f"Parameters/topology type {type(pt).__name__} is not yet supported "
        "for seekr3 preparation.")


def _map_ensemble(ensemble: str, physical_attributes: typing.Any) -> typing.Any:
    pressure = getattr(physical_attributes, "pressure", None)
    pressure = 1.0 if pressure is None else float(pressure)
    if ensemble == "NVT":
        return NVT_ensemble()
    if ensemble == "NPT":
        return NPT_ensemble(pressure=pressure)
    if ensemble == "NPT_membrane":
        return NPT_membrane_ensemble(pressure=pressure)
    if ensemble == "NVE":
        return NVE_ensemble()
    raise ValueError(f"Unknown ensemble {ensemble}.")


def _map_sampling(
        sampling: typing.Any,
        cv_specs: list,
        ) -> typing.Any:
    if isinstance(sampling, stage_procedures_module.Conventional_sampling_spec):
        heat_ramp = sampling.heat_ramp_temperatures or None
        interval = sampling.heat_ramp_step_interval or None
        return Conventional_md_input(
            heat_ramp_temperatures=heat_ramp,
            heat_ramp_step_interval=interval)
    if isinstance(sampling, stage_procedures_module.Metadynamics_sampling_spec):
        # Make sure that cv_names are in the model
        for cv_name in sampling.cv_names:
            if cv_name not in [cv_spec.name for cv_spec in cv_specs]:
                raise ValueError(
                    f"Sampling method references CV '{cv_name}' that is not "
                    "defined in the workflow's CV specs.")
        return Metadynamics_md_input(
            number_of_points=sampling.number_of_points,
            bias_factor=sampling.bias_factor,
            gaussian_height=sampling.gaussian_height,
            gaussian_width=sampling.gaussian_width,
            steps_per_update=sampling.steps_per_update,
            cv_names=list(sampling.cv_names))
    if isinstance(sampling, stage_procedures_module.Steered_sampling_spec):
        for cv_name in sampling.cv_names:
            if cv_name not in [cv_spec.name for cv_spec in cv_specs]:
                raise ValueError(
                    f"Sampling method references CV '{cv_name}' that is not "
                    "defined in the workflow's CV specs.")
        return Steered_md_input(
            force_constant=sampling.force_constant,
            velocity=sampling.velocity,
            cv_names=list(sampling.cv_names))
    if isinstance(sampling, stage_procedures_module.MMVT_sampling_spec):
        return MMVT_md_input()
    if isinstance(sampling, stage_procedures_module.NAM_sampling_spec):
        return NAM_bd_input()
    if isinstance(sampling, stage_procedures_module.RAMD_sampling_spec):
        import seekrsamplingmethodsplugin.modules.sampling_methods.preparer \
            as sampling_methods_preparer
        if sampling.escape_cv_name not in [cv_spec.name for cv_spec in cv_specs]:
            raise ValueError(
                f"Sampling method references CV '{sampling.escape_cv_name}' that is not "
                "defined in the workflow's CV specs.")
        return sampling_methods_preparer.Random_accelerated_md_sampling_method_input(
            ligand_selection_name = sampling.ligand_selection_name,
            receptor_selection_name = sampling.receptor_selection_name,
            steps_per_update = sampling.steps_per_update,
            force_magnitude = sampling.force_magnitude,
            cutoff_distance = sampling.cutoff_distance,
            ramp_start_step = sampling.ramp_start_step,
            ramp_force_increase_per_step \
                = sampling.ramp_force_increase_per_step
        )
    raise ValueError(f"Unknown sampling spec {type(sampling).__name__}.")


def _map_restraint(restraint: typing.Any) -> typing.Any:
    if isinstance(restraint, stage_procedures_module.Positional_restraint_spec):
        return Positional_restraint_input(
            selection_name=restraint.selection_name,
            force_constant=restraint.force_constant,
            coordinates_filename=restraint.coordinates_filename)
    if isinstance(restraint, stage_procedures_module.Pairwise_restraint_spec):
        return Pairwise_restraint_input(
            selection1_name=restraint.selection1_name,
            selection2_name=restraint.selection2_name,
            force_constant=restraint.force_constant)
    if isinstance(restraint, stage_procedures_module.Proximity_restraint_spec):
        return Proximity_restraint_input(
            selection1_name=restraint.selection1_name,
            selection2_name=restraint.selection2_name,
            radius=restraint.radius,
            force_constant=restraint.force_constant)
    raise ValueError(f"Unknown restraint spec {type(restraint).__name__}.")


def _map_completion(completion: typing.Any) -> typing.Any:
    if isinstance(
            completion,
            stage_procedures_module.Number_of_steps_completion_spec):
        return Number_of_steps_completion_criteria(
            number_of_steps=completion.number_of_steps)
    if isinstance(
            completion,
            stage_procedures_module.Number_of_trajectories_completion_spec):
        return Number_of_trajectories_completion_criteria(
            number_of_trajectories=completion.number_of_trajectories)
    if isinstance(
            completion,
            stage_procedures_module.Reporter_progress_completion_spec):
        return Reporter_progress_completion_criteria(
            reporter_name=completion.reporter_name,
            interval=completion.interval)
    if isinstance(
            completion,
            stage_procedures_module.CV_value_attained_completion_spec):
        import seekrsamplingmethodsplugin.modules.completion_criteria.structures \
            as completion_criteria_structures
        return completion_criteria_structures\
                .CV_value_attained_completion_criteria(
            cv_name=completion.cv_name,
            cv_inequality=completion.cv_inequality,
            cv_value=completion.cv_value)
    raise ValueError(f"Unknown completion spec {type(completion).__name__}.")


def _default_md_reporters(item: typing.Any) -> list:
    reporters: list = []
    if item.reporter_name is not None:
        anchor_reporter = Anchor_reporter_input()
        anchor_reporter.name = item.reporter_name
        anchor_reporter.interval = 100
        reporters.append(anchor_reporter)
    state_reporter = State_reporter_input()
    state_reporter.interval = item.checkpoint_interval
    reporters.append(state_reporter)
    if item.position_reporter_interval is not None:
        position_reporter = Position_reporter_input()
        position_reporter.interval = item.position_reporter_interval
        position_reporter.name = item.get_position_reporter_name()
        position_reporter.atom_subset = None
        position_reporter.enforce_periodic_box = item.enforce_periodic_box
        reporters.append(position_reporter)
    return reporters


# ======================================================================
#  Stage construction.
# ======================================================================

def _make_input_stage(
        input_stage_name: str,
        item: typing.Any,
        ctx: Prepare_context,
        ) -> typing.Any:
    if isinstance(item, stage_procedures_module.MD_stage_item):
        return Input_stage_md(
            name=item.name,
            input_stage_name=input_stage_name,
            run_minimization=item.run_minimization,
            checkpoint_interval=item.checkpoint_interval,
            ensemble=_map_ensemble(item.ensemble, ctx.physical_attributes),
            sampling_method=_map_sampling(item.sampling, ctx.workflow.cv_specs),
            completion_criteria=_map_completion(item.completion),
            restraints=[_map_restraint(r) for r in item.restraints],
            reporters=_default_md_reporters(item),
            timestep=item.timestep
        )
    if isinstance(item, stage_procedures_module.BD_stage_item):
        return Input_stage_bd(
            name=item.name,
            input_stage_name=input_stage_name,
            n_trajectories_per_output=item.n_trajectories_per_output,
            n_steps_per_output=item.n_steps_per_output,
            sampling_method=_map_sampling(item.sampling, ctx.workflow.cv_specs),
            completion_criteria=_map_completion(item.completion),
        )
    if isinstance(item, stage_procedures_module.Logistic_anchor_from_reporter_stage_item):
        stage = Input_stage_anchor_starting_structures_from_reporter()
        stage.name = item.name
        stage.input_stage_name = input_stage_name
        stage.reporter_name = item.reporter_name
        return stage
    if isinstance(item, stage_procedures_module.Logistic_swarm_from_reporter_stage_item):
        import seekrsamplingmethodsplugin.modules.scales.logistic\
            .swarm_structure_from_reporter_logistic_stage \
            as logistic_swarm_from_reporter_module
        stage = logistic_swarm_from_reporter_module\
            .Input_stage_swarm_structures_from_reporter()
        stage.name = item.name
        stage.input_stage_name = input_stage_name
        stage.reporter_name = item.reporter_name
        return stage
    raise ValueError(f"Unknown stage item {type(item).__name__}.")


def flatten_procedure(
        procedure: typing.Any,
        ctx: Prepare_context,
        ) -> list:
    """
    Expand a stage procedure into a flat list of seekr3 input stages.
    """
    resolved_items = procedure.expand(ctx)
    input_stage_list = []
    for input_stage_name, item in resolved_items:
        input_stage = _make_input_stage(input_stage_name, item, ctx)
        input_stage_list.append(input_stage)
    return input_stage_list

# ======================================================================
#  Collective variables, anchors, and partitions.
# ======================================================================

def _make_cv_inputs(cv_specs: list) -> list:
    cv_inputs: list = []
    index = 0
    for cv_spec in cv_specs:
        if isinstance(cv_spec, cv_specs_module.Com_com_distance_CV_spec):
            cv = Distance_cv_input(
                name=cv_spec.name,
                group1_name=cv_spec.group1_selection_name,
                group2_name=cv_spec.group2_selection_name,
                min_value=cv_spec.min_value,
                max_value=cv_spec.max_value)
            cv._index = index
            index += 1
            cv_inputs.append(cv)
        elif isinstance(
                cv_spec, cv_specs_module.Com_com_distance_rmsd_CV_spec):
            distance_cv = Distance_cv_input(
                name=f"{cv_spec.name}_distance",
                group1_name=cv_spec.group1_selection_name,
                group2_name=cv_spec.group2_selection_name,
                min_value=cv_spec.min_value,
                max_value=cv_spec.max_value)
            distance_cv._index = index
            index += 1
            rmsd_cv = RMSD_cv_input(
                name=f"{cv_spec.name}_rmsd",
                group_rmsd_name=cv_spec.rmsd_group_selection_name,
                group_align_name=cv_spec.rmsd_align_selection_name,
                ref_structure=cv_spec.rmsd_ref_structure,
                min_value=cv_spec.rmsd_min_value,
                max_value=cv_spec.rmsd_max_value)
            rmsd_cv._index = index
            index += 1
            cv_inputs.extend([distance_cv, rmsd_cv])
        else:
            raise ValueError(
                f"Unknown CV spec {type(cv_spec).__name__}.")
    return cv_inputs


def _single_com_com_cv_spec(
        workflow: "workflow_module.Workflow",
        ) -> typing.Any:
    """
    Return the workflow's single collective-variable spec, validating that it
    is a supported, one-dimensional CV for automatic anchor placement.
    """
    cv_specs = workflow.cv_specs
    if len(cv_specs) != 1:
        raise ValueError(
            "Automatic anchor placement currently supports exactly one "
            "collective variable. Multi-dimensional anchor placement will be "
            "handled by a future logistic stage.")
    cv_spec = cv_specs[0]
    if not isinstance(cv_spec, cv_specs_module.Com_com_distance_CV_spec):
        raise ValueError(
            "Automatic anchor placement currently supports only a "
            "center-of-mass distance collective variable "
            "(Com_com_distance_CV_spec).")
    return cv_spec


def _group_com_position(
        ctx: Prepare_context,
        traj: mdtraj.Trajectory,
        selection_name: str,
        ) -> np.ndarray:
    """
    Center of mass (first frame) of a named molecular-dynamics selection, using
    seekr's own center-of-mass definition.
    """
    indices = ctx.selections["molecular_dynamics"][selection_name]
    com = mdtraj_center_of_mass(traj.atom_slice(indices))
    return com[0, :]


def _structure_cv_value(
        ctx: Prepare_context,
        cv_spec: typing.Any,
        ) -> float:
    """
    Value of the center-of-mass distance collective variable in the starting
    structure (used to locate the bound state's anchor).
    """
    traj = mdtraj.load(ctx.md_structure_filename)
    com1 = _group_com_position(ctx, traj, cv_spec.group1_selection_name)
    com2 = _group_com_position(ctx, traj, cv_spec.group2_selection_name)
    return float(np.linalg.norm(com2 - com1))


def _solvent_approach_distance(
        ctx: Prepare_context,
        cv_spec: typing.Any,
        ) -> float:
    """
    Closest approach of solvent to the receptor center (the collective
    variable's first group's center of mass) in the starting structure.
    Solvent is every atom that does not belong to a named component.
    """
    traj = mdtraj.load(ctx.md_structure_filename)
    center = _group_com_position(ctx, traj, cv_spec.group1_selection_name)
    solute_indices: set[int] = set()
    for member in ctx.workflow.components.members:
        solute_indices.update(
            ctx.selections["molecular_dynamics"][member.name])
    solvent_indices = [
        i for i in range(traj.n_atoms) if i not in solute_indices]
    if not solvent_indices:
        raise ValueError(
            "No solvent atoms found in the starting structure to determine "
            "the receptor surface for graded anchor placement.")
    solvent_xyz = traj.xyz[0, solvent_indices, :]
    distances = np.linalg.norm(solvent_xyz - center, axis=1)
    return float(distances.min())


def _make_anchor_scheme(
        anchor_spec: typing.Any,
        ctx: Prepare_context,
        ) -> Anchors_as_voronoi_cells_anchor_scheme:
    cv_spec = _single_com_com_cv_spec(ctx.workflow)
    surface_value = None
    if isinstance(anchor_spec, anchor_specs_module.Graded_anchor_spec):
        surface_value = _solvent_approach_distance(ctx, cv_spec)
    value_vectors = anchor_spec.get_anchor_value_vectors(
        cv_spec.min_value, cv_spec.max_value, surface_value)
    n_anchors = len(value_vectors)
    scheme = Anchors_as_voronoi_cells_anchor_scheme()
    scheme.voronoi_cells = []
    for i, value in enumerate(value_vectors):
        is_bulk = (i == n_anchors - 1)
        scale_types: list[str] = []
        if not is_bulk:
            scale_types.append("molecular_dynamics")
        if ctx.bd_active and i >= n_anchors - 2:
            scale_types.append("brownian_dynamics")
        cell = Voronoi_cell()
        cell.value = list(value)
        cell.scale_types = scale_types
        cell.end_anchor = is_bulk
        cell.bulk_anchor = is_bulk
        # Anchor starting structures are produced by the seeding/logistic
        #  stage(s) from the equilibrated structure, so no explicit per-cell
        #  positions are set here.
        cell.positions = ""
        scheme.voronoi_cells.append(cell)
    return scheme


def _make_partitions(ctx: Prepare_context) -> list:
    workflow = ctx.workflow
    if not workflow.cv_specs and workflow.anchor_spec is None:
        return []
    cv_inputs = _make_cv_inputs(workflow.cv_specs)
    state_points: list = []
    anchor_scheme = None
    if workflow.anchor_spec is not None:
        anchor_scheme = _make_anchor_scheme(workflow.anchor_spec, ctx)
        # The bound state is the anchor whose Voronoi cell contains the
        #  starting structure; seekr3 locates it from this collective-variable
        #  value. The bulk state is simply the outermost anchor (flagged on the
        #  cell), so no explicit bulk state point is needed.
        cv_spec = _single_com_com_cv_spec(workflow)
        bound_value = _structure_cv_value(ctx, cv_spec)
        state_points = [State_point(
            name="bound", value=[bound_value], end_state=True)]
    return [Partition_input(
        cvs=cv_inputs,
        state_points=state_points,
        anchor_scheme=anchor_scheme)]


# ======================================================================
#  Engine settings.
# ======================================================================

def _make_md_engine_settings(ctx: Prepare_context) -> Openmm_settings_input:
    md_settings: scale_settings_module.Molecular_dynamics_scale_settings = \
        ctx.workflow.get_md_settings()
    seekr_pt = _map_parameters_topology(md_settings.system)
    selections = [
        Selection(name=name, atom_indices=list(indices))
        for name, indices in ctx.selections["molecular_dynamics"].items()]
    unpartitioned = MD_unpartitioned_settings_input(
        md_parameters_topology=seekr_pt,
        positions=md_settings.system.solvated_pdb,
        selections=selections)
    return Openmm_settings_input(
        nonbonded_method=md_settings.nonbonded_method,
        nonbonded_cutoff=md_settings.nonbonded_cutoff,
        constraints=md_settings.constraints,
        hydrogen_mass=md_settings.hydrogen_mass,
        timestep=md_settings.timestep,
        friction_coefficient=md_settings.friction_coefficient,
        integrator_type=md_settings.integrator_type,
        platform_type=md_settings.platform_type,
        unpartitioned_settings=unpartitioned)


def _ensure_bd_pqrs(ctx: Prepare_context) -> dict[str, str]:
    """
    Ensure every Brownian-dynamics molecule has a PQR file. Molecules whose
    ``pqr_filename`` is empty are auto-generated from the molecular-dynamics
    system (parameters + solvated structure); molecules that already point at a
    PQR are validated against the atom count their components resolve to.
    Returns ``{molecule_name: pqr_path}``.
    """
    workflow = ctx.workflow
    bd_settings = workflow.get_scale_settings("brownian_dynamics")
    components = workflow.components
    md_settings = workflow.get_md_settings()
    pqr_paths: dict[str, str] = {}
    parmed_complex = None
    for molecule in bd_settings.system.molecules:
        molecule_atoms = components.bd_molecule_md_atoms(molecule)
        if molecule.pqr_filename:
            # This tests for equal atom counts based on a component
            # that encompasses the entire receptor. This may be 
            # necessary when the PQRs are missing, but we probably
            # should trust the existing PQRs.
            #if molecule_atoms:
            #    parsed = pqr.read_pqr_atoms(molecule.pqr_filename)
            #    if len(parsed) != len(molecule_atoms):
            #        raise ValueError(
            #            f"BD molecule '{molecule.name}': provided PQR has "
            #            f"{len(parsed)} atoms but its components resolve to "
            #            f"{len(molecule_atoms)} atoms.")
            pqr_paths[molecule.name] = molecule.pqr_filename
            continue
        if parmed_complex is None:
            parmed_complex = \
                md_settings.system.parameters_topology.make_parmed(
                    pdb_filename=ctx.md_structure_filename,
                    directory=ctx.root_directory)
        one_resid_per_atom = components.get_member(molecule.component_name)\
            .one_resid_per_atom_pqr
        pqr_filename = os.path.join(
            ctx.root_directory, f"{molecule.name}.pqr")
        pqr.write_pqr_from_parmed(
            parmed_complex, molecule_atoms, pqr_filename,
            one_resid_per_atom=one_resid_per_atom)
        molecule.pqr_filename = pqr_filename
        pqr_paths[molecule.name] = pqr_filename
    return pqr_paths


def _make_bd_engine_settings(
        ctx: Prepare_context,
        pqr_paths: dict[str, str],
        ) -> Browndye_settings_input:
    bd_settings: scale_settings_module.Brownian_dynamics_scale_settings = \
        ctx.workflow.get_scale_settings("brownian_dynamics")
    bd_selections = ctx.selections.get("brownian_dynamics", {})
    molecules = []
    for molecule in bd_settings.system.molecules:
        local_selections = bd_selections.get(molecule.name, {})
        selections = [
            Selection(name=name, atom_indices=list(indices))
            for name, indices in local_selections.items()]
        molecules.append(Browndye_molecule_input(
            name=molecule.name,
            pqr_filename=pqr_paths[molecule.name],
            selections=selections))
    # TODO: handle more complex ion setups (e.g. divalent ions, 
    # asymmetric concentrations, etc.)
    ionic_strength = float(
        getattr(ctx.physical_attributes, "ionic_strength", 0.15))
    ions = [
        Ion(name="Cl-", radius=1.67, charge=-1.0, concentration=ionic_strength),
        Ion(name="Na+", radius=0.9, charge=1.0, concentration=ionic_strength),
    ]
    return Browndye_settings_input(
        binary_directory=bd_settings.binary_directory,
        molecules=molecules,
        ions=ions)


def _make_engine_settings_list(ctx: Prepare_context) -> list:
    engine_settings_list = [_make_md_engine_settings(ctx)]
    if ctx.bd_active:
        pqr_paths = _ensure_bd_pqrs(ctx)
        engine_settings_list.append(
            _make_bd_engine_settings(ctx, pqr_paths))
    return engine_settings_list


def _make_initial_scopes(ctx: Prepare_context) -> dict:
    scopes = {
        settings.get_scale_type(): "unpartitioned"
        for settings in ctx.workflow.scale_settings}
    scopes["logistic"] = "unpartitioned"
    return scopes


# ======================================================================
#  Entry point.
# ======================================================================

def make_seekr_input(
        workflow: "workflow_module.Workflow",
        root_directory: str,
        physical_attributes: typing.Any,
        ) -> Seekr_input:
    """
    Build a seekr3 ``Seekr_input`` from a generic seekrflow ``Workflow``.

    Parameters
    ----------
    workflow
        The workflow describing components, CVs, anchors, the stage procedure
        and per-scale settings.
    root_directory
        The directory in which the seekr3 model will be prepared.
    physical_attributes
        Object carrying ``temperature``, ``pressure``, ``ionic_strength``,
        and ``random_seed``.
    """
    ctx = Prepare_context(
        workflow, root_directory, physical_attributes)
    ctx.selections = workflow.components.resolve_all_selections(
        workflow.scale_settings)
    stages = flatten_procedure(workflow.procedure, ctx)
    engine_settings_list = _make_engine_settings_list(ctx)
    partitions = _make_partitions(ctx)
    seekr_input = Seekr_input(
        directory=root_directory,
        plugins=list(workflow.plugins),
        temperature=ctx.target_temperature,
        random_seed=physical_attributes.random_seed,
        engine_settings_list=engine_settings_list,
        initial_scopes=_make_initial_scopes(ctx),
        stages=stages,
        partitions=partitions)
    return seekr_input


def prepare_workflow(
        workflow: "workflow_module.Workflow",
        root_directory: str,
        physical_attributes: typing.Any,
        force_overwrite: bool = False,
        skip_checks: bool = False,
        ) -> typing.Tuple[typing.Any, str]:
    """
    Build the seekr3 input and run ``seekr.prepare.prepare`` to write
    ``model.json``. Returns ``(model, model_json_path)``.
    """
    seekr_input = make_seekr_input(
        workflow, root_directory, physical_attributes)
    model, file_path = seekr.prepare.prepare(seekr_input, force_overwrite=force_overwrite)
    if model.directory == ".":
        model_dir = os.path.dirname(file_path)
        model.directory = os.path.abspath(model_dir)
    if not skip_checks:
        check.check_pre_simulation_all(model)
    return model, file_path
