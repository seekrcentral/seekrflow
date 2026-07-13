"""
modules/workflows/stage_procedures.py

Stage procedures are ORDERED, NESTABLE building blocks that expand into a flat
list of stage items. Each stage item is a seekrflow-native, JSON-serializable
description of one seekr3 input stage. ``modules/workflows/prepare.py`` converts
stage items into seekr3 input-stage objects and chains them via
``input_stage_name``.

Design notes:
  * Single-input / single-output chaining only (linear). Branching/joins are a
    future feature implemented as a dedicated logistic-only stage procedure.
  * The reserved sentinel ``"initial"`` means "depends only on the prepared
    structures".
  * Sampling methods, restraints and completion criteria are thin native specs
    here; they map 1:1 onto seekr3 objects in prepare.py. This keeps the JSON
    seekrflow-owned and stable across seekr3 versions.
"""

import os
import typing

import numpy as np
from attrs import define, field, validators, Factory
import mdtraj

import seekrflow.modules.workflows.utilities as utilities_module

INITIAL = "initial"


# ======================================================================
#  Sampling-method specs.
# ======================================================================

@define
class Sampling_spec_base:
    type: typing.Literal["base"] = "base"


@define
class Conventional_sampling_spec(Sampling_spec_base):
    """
    Plain MD, optionally with a temperature heat ramp.
    """
    type: typing.Literal["conventional"] = "conventional"
    heat_ramp_temperatures: list[float] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(float),
            iterable_validator=validators.instance_of(list),
        ))
    heat_ramp_step_interval: int = field(
        default=0,
        validator=validators.instance_of(int),
        )


@define
class Metadynamics_sampling_spec(Sampling_spec_base):
    type: typing.Literal["metadynamics"] = "metadynamics"
    cv_names: list[str] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(str),
            iterable_validator=validators.instance_of(list),
        ))
    bias_factor: float = field(
        default=10.0, validator=validators.instance_of(float))
    gaussian_height: float = field(
        default=0.5, validator=validators.instance_of(float))
    gaussian_width: float = field(
        default=0.1, validator=validators.instance_of(float))
    steps_per_update: int = field(
        default=250, validator=validators.instance_of(int))
    number_of_points: int = field(
        default=101, validator=validators.instance_of(int))


@define
class Steered_sampling_spec(Sampling_spec_base):
    type: typing.Literal["steered"] = "steered"
    cv_names: list[str] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(str),
            iterable_validator=validators.instance_of(list),
        ))
    force_constant: float = field(
        default=10000.0, validator=validators.instance_of(float))
    velocity: float = field(
        default=1.0, validator=validators.instance_of(float))


@define
class MMVT_sampling_spec(Sampling_spec_base):
    type: typing.Literal["mmvt"] = "mmvt"


@define
class NAM_sampling_spec(Sampling_spec_base):
    """
    Northrup-Allison-McCammon Brownian-dynamics sampling.
    """
    type: typing.Literal["nam"] = "nam"


@define
class RAMD_sampling_spec(Sampling_spec_base):
    """
    Random-acceleration molecular dynamics sampling.
    """
    type: typing.Literal["ramd"] = "ramd"
    swarm_size: int = field(
        default=20, validator=validators.instance_of(int))
    ligand_selection_name: str = field(
        default="ligand", validator=validators.instance_of(str))
    receptor_selection_name: str = field(
        default="receptor", validator=validators.instance_of(str))
    escape_cv_name: str = field(
        default="escape_cv", validator=validators.instance_of(str))
    steps_per_update: int = field(
        default=50, validator=validators.instance_of(int))
    force_magnitude: float = field(
        default=14.0, validator=validators.instance_of(float))
    cutoff_distance: float = field(
        default=0.0025, validator=validators.instance_of(float))
    ramp_start_step: int = field(
        default=0, validator=validators.instance_of(int))
    ramp_force_increase_per_step: float = field(
        default=0.0, validator=validators.instance_of(float))

Sampling_spec = (
    Conventional_sampling_spec | Metadynamics_sampling_spec
    | Steered_sampling_spec | MMVT_sampling_spec | NAM_sampling_spec
    # | RAMD_sampling_spec
)

# ======================================================================
#  Restraint specs (reference selections by name).
# ======================================================================

@define
class Restraint_spec_base:
    type: typing.Literal["base"] = "base"


@define
class Positional_restraint_spec(Restraint_spec_base):
    type: typing.Literal["positional"] = "positional"
    selection_name: str = field(
        default="", validator=validators.instance_of(str))
    force_constant: float = field(
        default=10.0, validator=validators.instance_of(float))
    coordinates_filename: str | None = field(
        default=None, validator=validators.optional(
            validators.instance_of(str)))
    


@define
class Pairwise_restraint_spec(Restraint_spec_base):
    type: typing.Literal["pairwise"] = "pairwise"
    selection1_name: str = field(
        default="", validator=validators.instance_of(str))
    selection2_name: str = field(
        default="", validator=validators.instance_of(str))
    force_constant: float = field(
        default=1000.0, validator=validators.instance_of(float))


@define
class Proximity_restraint_spec(Restraint_spec_base):
    type: typing.Literal["proximity"] = "proximity"
    selection1_name: str = field(
        default="", validator=validators.instance_of(str))
    selection2_name: str = field(
        default="", validator=validators.instance_of(str))
    radius: float = field(
        default=2.0, validator=validators.instance_of(float))
    force_constant: float = field(
        default=1000.0, validator=validators.instance_of(float))


Restraint_spec = (
    Positional_restraint_spec | Pairwise_restraint_spec
    | Proximity_restraint_spec
)


# ======================================================================
#  Completion-criteria specs.
# ======================================================================

@define
class Completion_spec_base:
    type: typing.Literal["base"] = "base"


@define
class Number_of_steps_completion_spec(Completion_spec_base):
    type: typing.Literal["number_of_steps"] = "number_of_steps"
    number_of_steps: int = field(
        default=100000, validator=validators.instance_of(int))


@define
class Number_of_trajectories_completion_spec(Completion_spec_base):
    type: typing.Literal["number_of_trajectories"] = "number_of_trajectories"
    number_of_trajectories: int = field(
        default=1000, validator=validators.instance_of(int))


@define
class Reporter_progress_completion_spec(Completion_spec_base):
    type: typing.Literal["reporter_progress"] = "reporter_progress"
    reporter_name: str = field(
        default="", validator=validators.instance_of(str))
    interval: int = field(
        default=100, validator=validators.instance_of(int))

@define 
class CV_value_attained_completion_spec(Completion_spec_base):
    type: typing.Literal["cv_value_attained"] = "cv_value_attained"
    eval_interval: int = field(
        default=100, validator=validators.instance_of(int))
    cv_name: str = field(
        default="", validator=validators.instance_of(str))
    cv_inequality: str = field(
        default="greater_than", validator=validators.in_(
            ["greater_than", "less_than", "greater_than_or_equal_to", 
             "less_than_or_equal_to"]))
    cv_value: float = field(
        default=0.0, validator=validators.instance_of(float))
    
Completion_spec = (
    Number_of_steps_completion_spec | Number_of_trajectories_completion_spec
    | Reporter_progress_completion_spec | CV_value_attained_completion_spec
)


# ======================================================================
#  Stage items (leaf nodes; one seekr3 input stage each).
# ======================================================================

@define
class Stage_item_base:
    type: typing.Literal["base"] = "base"
    name: str = field(
        default="stage", validator=validators.instance_of(str))
    scope: str = field(
        default="unpartitioned", validator=validators.instance_of(str))


@define
class MD_stage_item(Stage_item_base):
    type: typing.Literal["md"] = "md"
    run_minimization: bool = field(
        default=False, validator=validators.instance_of(bool))
    checkpoint_interval: int = field(
        default=1000, validator=validators.instance_of(int))
    ensemble: str = field(
        default="NVT", validator=validators.in_(
            ["NVT", "NPT", "NPT_membrane", "NVE"]))
    sampling: Sampling_spec = field(
        default=Factory(Conventional_sampling_spec),
        validator=validators.instance_of(Sampling_spec_base))
    restraints: list[Restraint_spec] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(Restraint_spec_base),
            iterable_validator=validators.instance_of(list),
        ))
    completion: Completion_spec = field(
        default=Factory(Number_of_steps_completion_spec),
        validator=validators.instance_of(Completion_spec_base))
    #NOTE: None defaults to the MD_settings default.
    timestep: float | None = field(
        default=None, validator=validators.optional(validators.instance_of(float)))
    reporter_name: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)))
    position_reporter_interval: int | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(int)))
    enforce_periodic_box: bool = field(
        default=True, validator=validators.instance_of(bool))

    def get_position_reporter_name(self) -> str:
        return f"{self.name}_position_reporter"

@define
class BD_stage_item(Stage_item_base):
    type: typing.Literal["bd"] = "bd"
    scope: str = field(
        default="unpartitioned", validator=validators.instance_of(str))
    sampling: Sampling_spec = field(
        default=Factory(NAM_sampling_spec),
        validator=validators.instance_of(Sampling_spec_base))
    n_trajectories_per_output: int = field(
        default=1000, validator=validators.instance_of(int))
    n_steps_per_output: int = field(
        default=100000, validator=validators.instance_of(int))
    completion: Completion_spec = field(
        default=Factory(Number_of_trajectories_completion_spec),
        validator=validators.instance_of(Completion_spec_base))


@define
class Logistic_anchor_from_reporter_stage_item(Stage_item_base):
    """
    A logistic stage that assigns anchor starting structures from a reporter.
    """
    type: typing.Literal["logistic_anchor_from_reporter"] = "logistic_anchor_from_reporter"
    scope: str = field(
        default="partitioned", validator=validators.instance_of(str))
    reporter_name: str = field(
        default="", validator=validators.instance_of(str))

@define
class Logistic_swarm_from_reporter_stage_item(Stage_item_base):
    """
    A logistic stage that assigns anchor starting structures from a reporter.
    """
    type: typing.Literal["logistic_swarm_from_reporter"] = "logistic_swarm_from_reporter"
    scope: str = field(
        default="partitioned", validator=validators.instance_of(str))
    reporter_name: str = field(
        default="", validator=validators.instance_of(str))


Stage_item = MD_stage_item | BD_stage_item \
    | Logistic_swarm_from_reporter_stage_item \
    | Logistic_anchor_from_reporter_stage_item


# ======================================================================
#  A reasonable default equilibration schedule.
#
#  NOTE: this is a SIMPLIFIED replacement for the old imperative 8-stage
#  backbone/sidechain/membrane restraint ladder. seekr3 positional restraints
#  reference a single named selection, so we relax a single restrained
#  selection across a heat ramp and progressively weaker holds. Force constants
#  are in seekr3 restraint units.
# ======================================================================

# Each row: (name, run_minimization, restraint_force_constant, number_of_steps,
#            heat_ramp_start_temperature_or_None)
DEFAULT_GLOBULAR_EQUIL_SCHEDULE: list[dict] = [
    {"name": "equil_stage_1", 
     "run_minimization": True,
     "backbone_force_constant": 41840.0, 
     "sidechain_force_constant": 41840.0, 
     "nonprotein_force_constant": 41840.0, 
     "number_of_steps": 10000,
     "heat_ramp_start": 100.0,
     "heat_ramp_stride": 10.0,
     "heat_ramp_step_interval": 1000,
     "timestep": 0.001,
     "ensemble": "NPT"},
    {"name": "equil_stage_2", 
     "run_minimization": False,
     "backbone_force_constant": 41840.0, 
     "sidechain_force_constant": 41840.0, 
     "nonprotein_force_constant": 41840.0,
     "number_of_steps": 100000,
     "heat_ramp_start": None,
     "heat_ramp_stride": None,
     "heat_ramp_step_interval": None,
     "timestep": 0.001,
     "ensemble": "NPT"},
    {"name": "equil_stage_3", 
     "run_minimization": False,
     "backbone_force_constant": 4184.0, 
     "sidechain_force_constant": 4184.0, 
     "nonprotein_force_constant": 4184.0, 
     "number_of_steps": 10000,
     "heat_ramp_start": 100.0,
     "heat_ramp_stride": 10.0,
     "heat_ramp_step_interval": 1000,
     "timestep": 0.001,
     "ensemble": "NPT"},
    {"name": "equil_stage_4",
     "run_minimization": False,
     "backbone_force_constant": 4184.0, 
     "sidechain_force_constant": 4184.0, 
     "nonprotein_force_constant": 4184.0,
     "number_of_steps": 100000,
     "heat_ramp_start": None,
     "heat_ramp_stride": None,
     "heat_ramp_step_interval": None,
     "timestep": 0.001,
     "ensemble": "NPT"},
    {"name": "equil_stage_5",
     "run_minimization": False,
     "backbone_force_constant": 418.4, 
     "sidechain_force_constant": 418.4, 
     "nonprotein_force_constant": 418.4,
     "number_of_steps": 200000,
     "heat_ramp_start": None,
     "heat_ramp_stride": None,
     "heat_ramp_step_interval": None,
     "timestep": 0.001,
     "ensemble": "NPT"},
    {"name": "equil_stage_6",
     "run_minimization": False,
     "backbone_force_constant": 418.4, 
     "sidechain_force_constant": 0.0, 
     "nonprotein_force_constant": 0.0,
     "number_of_steps": 200000,
     "heat_ramp_start": None,
     "heat_ramp_stride": None,
     "heat_ramp_step_interval": None,
     "timestep": 0.002,
     "ensemble": "NPT"},
    {"name": "equil_stage_7",
     "run_minimization": False,
     "backbone_force_constant": 0.0, 
     "sidechain_force_constant": 0.0, 
     "nonprotein_force_constant": 0.0,
     "number_of_steps": 500000,
     "heat_ramp_start": None,
     "heat_ramp_stride": None,
     "heat_ramp_step_interval": None,
     "timestep": 0.002,
     "ensemble": "NPT"}
]


# ======================================================================
#  Stage procedures (ordered, nestable; expand into stage items).
# ======================================================================

# A resolved stage item paired with the name of the stage it depends on. The
#  dependency is computed during expansion and consumed by prepare.py; it is
#  deliberately NOT stored on the (serializable) stage item, since it is fully
#  determined by procedure order rather than set by the user.
Resolved_stage_item = typing.Tuple[str, Stage_item_base]

# (address path relative to the owning procedure, role for policy defaults)
Stage_address_info = typing.Tuple[list[str], str | None]


# TODO: move to structures?
@define
class Dispatch:
    """
    Define how anchors/swarms within a stage are 'dispatched' to a resource,
    i.e. how the anchors or swarms are run in parallel or serially in an
    array (group size) or concurrently in separate processes. The ``dimensions`` 
    field defines whether parallelism is used across anchors or swarms or both.
    Group size is the number of anchors/swarms to run within a single array job. 
    Concurrency is the number of anchors/swarms to run in separate processes 
    simultaneously. If group_size is larger than concurrency, then the anchors/
    swarms will run serially within the same array job, only 'concurrency' at a 
    time.
    """
    type: typing.Literal["dispatch"] = "dispatch"
    dimensions: list[str] | None = field(
        default=None,
        validator=validators.optional(validators.deep_iterable(
            member_validator=validators.in_(["anchor", "swarm"]),
            iterable_validator=validators.instance_of(list),
        )))
    group_size: int | None = field(
        default=None,
        validator=validators.optional(
            [validators.instance_of(int), validators.ge(1)]))
    concurrency: int = field(
        default=1,
        validator=[validators.instance_of(int), validators.ge(1)])

    def merged_with(self, override: "Dispatch") -> "Dispatch":
        """
        Field-wise merge: non-None ``dimensions`` / ``group_size`` on
        ``override`` replace; ``concurrency`` is always taken from ``override``.
        """
        return Dispatch(
            dimensions=(
                override.dimensions
                if override.dimensions is not None else self.dimensions),
            group_size=(
                override.group_size
                if override.group_size is not None else self.group_size),
            concurrency=override.concurrency,
        )


@define
class Emitted_stage:
    """One stage a procedure emits, with intrinsic execution defaults."""
    name: str
    role: str | None = field(default=None)
    dispatch: Dispatch = field(default=Factory(Dispatch))
    co_schedule_with: str | None = field(
        default=None,
        validator=validators.optional(
            validators.in_(["predecessor", "successor"])))


@define
class Stage_procedure_base:
    type: typing.Literal["base"] = "base"
    name: str = field(
        default="procedure", validator=validators.instance_of(str))
    position_reporter_interval: int | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(int)))

    def expand(
            self,
            ctx: typing.Any,
            input_stage_name: str = INITIAL,
            ) -> list[Resolved_stage_item]:
        """
        Expand this procedure into a flat, internally-chained list of
        ``(input_stage_name, stage_item)`` pairs. The first item depends on
        ``input_stage_name`` (defaulting to the ``INITIAL`` sentinel, i.e. the
        prepared structures); subsequent items chain to their predecessor.
        """
        raise NotImplementedError

    def get_emitted_stages(self) -> list[Emitted_stage]:
        """
        Return, in order, each stage this procedure emits with its intrinsic
        execution defaults. ``role`` is the trailing address segment; ``None``
        means the procedure emits a single addressable stage (no sub-role
        segment). ctx-free, like :meth:`expand`.
        """
        raise NotImplementedError

    def get_stage_names(self) -> list[str]:
        return [stage.name for stage in self.get_emitted_stages()]


def _chain(items: list[Stage_item_base], input_stage_name: str
           ) -> list[Resolved_stage_item]:
    """
    Wire a list of items into a linear chain of ``(predecessor_name, item)``
    pairs.
    """
    resolved: list[Resolved_stage_item] = []
    previous = input_stage_name
    for item in items:
        resolved.append((previous, item))
        previous = item.name
    return resolved


@define
class Explicit_stage_procedure(Stage_procedure_base):
    """
    A procedure consisting of explicitly-listed stage items, used verbatim.
    """
    type: typing.Literal["explicit"] = "explicit"
    items: list[Stage_item_base] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(Stage_item_base),
            iterable_validator=validators.instance_of(list),
        ))

    def expand(
            self,
            ctx: typing.Any,
            input_stage_name: str = INITIAL,
            ) -> list[Resolved_stage_item]:
        return _chain(list(self.items), input_stage_name)

    def get_emitted_stages(self) -> list[Emitted_stage]:
        return [
            Emitted_stage(name=item.name, role=item.name)
            for item in self.items
        ]


@define
class Equilibration_globular_stage_procedure(Stage_procedure_base):
    """
    A heat-ramp + progressive restraint-relaxation equilibration, expanded from
    a schedule of positional-restraint holds on a single named selection.
    """
    type: typing.Literal["equilibration_globular"] = "equilibration_globular"
    reference_structure_filename: str | None = field(
        default=None, validator=validators.optional(validators.instance_of(str)))
    schedule: list[dict] = field(
        default=Factory(lambda: [dict(row) for row in DEFAULT_GLOBULAR_EQUIL_SCHEDULE]),
        validator=validators.instance_of(list))

    def expand(
            self,
            ctx: typing.Any,
            input_stage_name: str = INITIAL,
            ) -> list[Resolved_stage_item]:
        target_temperature = getattr(ctx, "target_temperature", 298.15)
        items: list[Stage_item_base] = []
        for row in self.schedule:
            sampling = Conventional_sampling_spec()
            if row.get("heat_ramp_start") is not None:
                heat_ramp_start = float(row["heat_ramp_start"])
                heat_ramp_end = target_temperature
                heat_ramp_stride = float(row["heat_ramp_stride"])
                heat_ramp_temperatures = np.arange(
                    heat_ramp_start, heat_ramp_end, heat_ramp_stride)
                sampling.heat_ramp_temperatures = heat_ramp_temperatures.tolist()
                sampling.heat_ramp_temperatures.append(target_temperature)
                sampling.heat_ramp_step_interval = int(row["heat_ramp_step_interval"])
                number_of_steps = int(row["heat_ramp_step_interval"]) * len(heat_ramp_temperatures)
            
            else:
                number_of_steps = int(row["number_of_steps"])
            
            if self.reference_structure_filename is not None:
                ref_mdtraj = mdtraj.load(self.reference_structure_filename)
            else:
                ref_mdtraj = mdtraj.load(ctx.md_structure_filename)
            md_mdtraj = mdtraj.load(ctx.md_structure_filename)
            ref_atom_indices, md_atom_indices = utilities_module\
                .obtain_md_ref_structure_map(ref_mdtraj, md_mdtraj)
            restraints: list[Restraint_spec] = []
            
            # Backbone restraints
            backbone_selection = "protein and backbone and not element H"
            protein_backbone_md_indices, protein_backbone_ref_indices \
                = utilities_module.select_both_from_mdtraj_and_indices(
                    md_mdtraj, md_atom_indices, ref_mdtraj, ref_atom_indices,
                    backbone_selection
                )
            BACKBONE_SELECTION_NAME = "backbone_for_equilibration"
            assert len(protein_backbone_ref_indices) == len(protein_backbone_md_indices), \
                "Number of protein backbone atoms in solvated and unsolvated structures "\
                "should be the same"
            backbone_coordinates_filename = os.path.join(
                ctx.root_directory, f"{self.name}_backbone_restraint_coordinates.npy")
            np.savetxt(backbone_coordinates_filename, 
                       ref_mdtraj.xyz[0, protein_backbone_ref_indices])
            if row.get("backbone_force_constant", 0.0) > 0.0:
                restraints.append(Positional_restraint_spec(
                    selection_name=BACKBONE_SELECTION_NAME,
                    force_constant=float(row["backbone_force_constant"]),
                    coordinates_filename=backbone_coordinates_filename,
                ))
            ctx.selections["molecular_dynamics"][BACKBONE_SELECTION_NAME] \
                = protein_backbone_md_indices

            # Sidechain restraints
            SIDECHAIN_SELECTION_NAME = "sidechain_for_equilibration"
            sidechain_selection = "protein and sidechain and not element H"
            protein_sidechain_md_indices, protein_sidechain_ref_indices \
                = utilities_module.select_both_from_mdtraj_and_indices(
                    md_mdtraj, md_atom_indices, ref_mdtraj, ref_atom_indices,
                    sidechain_selection)
            assert len(protein_sidechain_ref_indices) == len(protein_sidechain_md_indices), \
                "Number of protein sidechain atoms in solvated and unsolvated structures "\
                "should be the same"
            sidechain_coordinates_filename = os.path.join(
                ctx.root_directory, f"{self.name}_sidechain_restraint_coordinates.npy")
            np.savetxt(sidechain_coordinates_filename, 
                       ref_mdtraj.xyz[0, protein_sidechain_ref_indices])
            if row.get("sidechain_force_constant", 0.0) > 0.0:
                restraints.append(Positional_restraint_spec(
                    selection_name=SIDECHAIN_SELECTION_NAME,
                    force_constant=float(row["sidechain_force_constant"]),
                    coordinates_filename=sidechain_coordinates_filename,
                ))
            ctx.selections["molecular_dynamics"][SIDECHAIN_SELECTION_NAME] \
                = protein_sidechain_md_indices
            
            # Nonprotein restraints
            NONPROTEIN_SELECTION_NAME = "nonprotein_for_equilibration"
            nonprotein_selection = "not protein and not element H"
            nonprotein_md_indices, nonprotein_ref_indices \
                = utilities_module.select_both_from_mdtraj_and_indices(
                    md_mdtraj, md_atom_indices, ref_mdtraj, ref_atom_indices,
                    nonprotein_selection)
            assert len(nonprotein_ref_indices) == len(nonprotein_md_indices), \
                "Number of non-protein atoms in solvated and unsolvated structures "\
                "should be the same"
            nonprotein_coordinates_filename = os.path.join(
                ctx.root_directory, f"{self.name}_nonprotein_restraint_coordinates.npy")
            np.savetxt(nonprotein_coordinates_filename, 
                       ref_mdtraj.xyz[0, nonprotein_ref_indices])
            if row.get("nonprotein_force_constant", 0.0) > 0.0:
                restraints.append(Positional_restraint_spec(
                    selection_name=NONPROTEIN_SELECTION_NAME,
                    force_constant=float(row["nonprotein_force_constant"]),
                    coordinates_filename=nonprotein_coordinates_filename,
                ))
            ctx.selections["molecular_dynamics"][NONPROTEIN_SELECTION_NAME] \
                = nonprotein_md_indices

            items.append(MD_stage_item(
                name=row["name"],
                run_minimization=bool(row.get("run_minimization", False)),
                ensemble=row["ensemble"],
                sampling=sampling,
                restraints=restraints,
                completion=Number_of_steps_completion_spec(
                    number_of_steps=number_of_steps),
                timestep=float(row["timestep"]),
                position_reporter_interval=self.position_reporter_interval,
                ))
        return _chain(items, input_stage_name)

    def get_emitted_stages(self) -> list[Emitted_stage]:
        return [
            Emitted_stage(
                name=row["name"],
                role=row["name"],
            )
            for row in self.schedule
        ]

#@define
#class Seeding_method_input:
#    type: typing.Literal["base"] = "base"

@define
class Steered_seeding_method_input: #(Seeding_method_input):
    """
    A steered MD input specification.
    """
    type: typing.Literal["steered"] = "steered"
    force_constant: float = field(
        default=10000.0, validator=validators.instance_of(float))
    velocity: float = field(
        default=1.0, validator=validators.instance_of(float))

@define
class Metadynamics_seeding_method_input: #(Seeding_method_input):
    """
    A metadynamics input specification.
    """
    type: typing.Literal["metadynamics"] = "metadynamics"
    bias_factor: float = field(
        default=10.0, validator=validators.instance_of(float))
    gaussian_height: float = field(
        default=0.5, validator=validators.instance_of(float))
    gaussian_width: float = field(
        default=0.1, validator=validators.instance_of(float))
    steps_per_update: int = field(
        default=250, validator=validators.instance_of(int))
    number_of_points: int = field(
        default=101, validator=validators.instance_of(int))

Seeding_method_input = Steered_seeding_method_input \
                     | Metadynamics_seeding_method_input

@define
class Seeding_stage_procedure(Stage_procedure_base):
    """
    A seeding procedure: a biased MD stage (metadynamics or steered MD) that
    drives the system along the collective variable(s) and records an anchor
    reporter, followed by a logistic stage that assigns per-anchor starting
    structures from that reporter.
    """
    type: typing.Literal["seeding"] = "seeding"
    # TODO: implement additional per-method parameters - this will
    #  be objects like Sampling_spec.
    #method: str = field(
    #    default="metadynamics",
    #    validator=validators.in_(["metadynamics", "steered"]))
    method_input: Seeding_method_input = field(
        default=Factory(Metadynamics_seeding_method_input),
        validator=validators.instance_of(Seeding_method_input))
    cv_names: list[str] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(str),
            iterable_validator=validators.instance_of(list),
        ))
    ensemble: str = field(
        default="NVT", validator=validators.in_(
            ["NVT", "NPT", "NPT_membrane", "NVE"]))
    reporter_interval: int = field(
        default=100, validator=validators.instance_of(int))
    

    def get_emitted_stages(self) -> list[Emitted_stage]:
        return [
            Emitted_stage(name=f"{self.name}_seed", role="sampling"),
            Emitted_stage(
                name=f"{self.name}_assign_structures",
                role="logistic",
                co_schedule_with="predecessor"),
        ]

    def expand(
            self,
            ctx: typing.Any,
            input_stage_name: str = INITIAL,
            ) -> list[Resolved_stage_item]:
        seed_name, assign_name = self.get_stage_names()
        reporter_name = f"{self.name}_anchor_reporter"
        if self.method_input.type == "metadynamics":
            sampling: Sampling_spec = Metadynamics_sampling_spec(
                cv_names=list(self.cv_names), bias_factor=self.method_input.bias_factor,
                gaussian_height=self.method_input.gaussian_height,
                gaussian_width=self.method_input.gaussian_width,
                steps_per_update=self.method_input.steps_per_update,
                number_of_points=self.method_input.number_of_points)
        elif self.method_input.type == "steered":
            sampling = Steered_sampling_spec(
                cv_names=list(self.cv_names),
                force_constant=self.method_input.force_constant,
                velocity=self.method_input.velocity)
        else:
            raise ValueError(f"Unknown seeding method type: {self.method_input.type}")
        seed_stage = MD_stage_item(
            name=seed_name,
            scope="unpartitioned",
            ensemble=self.ensemble,
            sampling=sampling,
            completion=Reporter_progress_completion_spec(
                reporter_name=reporter_name, interval=self.reporter_interval),
            reporter_name=reporter_name,
            position_reporter_interval=self.position_reporter_interval,
        )
        logistic_stage = Logistic_anchor_from_reporter_stage_item(
            name=assign_name,
            scope="partitioned",
            reporter_name=reporter_name,
        )
        return _chain([seed_stage, logistic_stage], input_stage_name)


@define
class MMVT_stage_procedure(Stage_procedure_base):
    """
    A Markovian Milestoning with Voronoi Tessellations (MMVT) sampling stage,
    run in the partitioned scope (one milestoning simulation per anchor).
    """
    type: typing.Literal["mmvt"] = "mmvt"
    ensemble: str = field(
        default="NVT", validator=validators.in_(
            ["NVT", "NPT", "NPT_membrane", "NVE"]))
    number_of_steps: int = field(
        default=200000, validator=validators.instance_of(int))
    checkpoint_interval: int = field(
        default=5000, validator=validators.instance_of(int))

    def get_emitted_stages(self) -> list[Emitted_stage]:
        return [
            Emitted_stage(
                name=self.name,
                dispatch=Dispatch(dimensions=["anchor"], group_size=1)),
        ]

    def expand(
            self,
            ctx: typing.Any,
            input_stage_name: str = INITIAL,
            ) -> list[Resolved_stage_item]:
        mmvt_stage = MD_stage_item(
            name=self.name,
            scope="partitioned",
            ensemble=self.ensemble,
            checkpoint_interval=self.checkpoint_interval,
            sampling=MMVT_sampling_spec(),
            completion=Number_of_steps_completion_spec(
                number_of_steps=self.number_of_steps),
            position_reporter_interval=self.position_reporter_interval,
        )
        return _chain([mmvt_stage], input_stage_name)


@define
class BD_stage_procedure(Stage_procedure_base):
    """
    A Brownian-dynamics stage (rate of association into the outermost anchors).
    """
    type: typing.Literal["bd"] = "bd"
    n_trajectories_per_output: int = field(
        default=1000, validator=validators.instance_of(int))
    n_steps_per_output: int = field(
        default=100000, validator=validators.instance_of(int))
    number_of_trajectories: int = field(
        default=10000, validator=validators.instance_of(int))

    def get_emitted_stages(self) -> list[Emitted_stage]:
        return [Emitted_stage(name=self.name)]

    def expand(
            self,
            ctx: typing.Any,
            input_stage_name: str = INITIAL,
            ) -> list[Resolved_stage_item]:
        bd_stage = BD_stage_item(
            name=self.name,
            scope="unpartitioned",
            sampling=NAM_sampling_spec(),
            n_trajectories_per_output=self.n_trajectories_per_output,
            n_steps_per_output=self.n_steps_per_output,
            completion=Number_of_trajectories_completion_spec(
                number_of_trajectories=self.number_of_trajectories),
        )
        return _chain([bd_stage], input_stage_name)

@define
class RAMD_stage_procedure(Stage_procedure_base):
    """
    A RAMD procedure.
    """
    type: typing.Literal["ramd"] = "ramd"
    # TODO: implement additional per-method parameters - this will
    #  be objects like Sampling_spec.
    swarm_size: int = field(
        default=20, validator=validators.instance_of(int))
    starting_structure_interval: int = field(
        default=2000, validator=validators.instance_of(int))
    ligand_selection_name: str = field(
        default="ligand", validator=validators.instance_of(str))
    receptor_selection_name: str = field(
        default="receptor", validator=validators.instance_of(str))
    escape_cv_name: str = field(
        default="escape_cv", validator=validators.instance_of(str))
    escape_distance: float = field(
        default=3.0, validator=validators.instance_of(float))
    steps_per_update: int = field(
        default=50, validator=validators.instance_of(int))
    force_magnitude: float = field(
        default=14.0, validator=validators.instance_of(float))
    cutoff_distance: float = field(
        default=0.0025, validator=validators.instance_of(float))
    ramp_start_step: int = field(
        default=0, validator=validators.instance_of(int))
    ramp_force_increase_per_ns: float = field(
        default=0.0, validator=validators.instance_of(float))
    ensemble: str = field(
        default="NVT", validator=validators.in_(
            ["NVT", "NPT", "NPT_membrane", "NVE"]))

    def get_emitted_stages(self) -> list[Emitted_stage]:
        return [
            Emitted_stage(name=f"{self.name}_sampling", role="sampling"),
            Emitted_stage(
                name=f"{self.name}_assign_swarm",
                role="logistic",
                co_schedule_with="predecessor"),
            Emitted_stage(
                name=f"{self.name}_ramd",
                role="ramd",
                dispatch=Dispatch(dimensions=["swarm"], group_size=1)),
        ]

    def expand(
            self,
            ctx: typing.Any,
            input_stage_name: str = INITIAL,
            ) -> list[Resolved_stage_item]:
        sampling_name, assign_swarm_name, ramd_name = self.get_stage_names()
        equil_sampling = Conventional_sampling_spec()
        number_of_equil_steps = self.starting_structure_interval \
            * self.swarm_size
        step_size = ctx.workflow.get_md_settings().timestep
        NS_PER_STEP = step_size * 1e-3
        ramd_sampling = RAMD_sampling_spec(
            swarm_size=self.swarm_size,
            ligand_selection_name=self.ligand_selection_name,
            receptor_selection_name=self.receptor_selection_name,
            escape_cv_name=self.escape_cv_name,
            steps_per_update=self.steps_per_update,
            force_magnitude=self.force_magnitude,
            cutoff_distance=self.cutoff_distance,
            ramp_start_step=self.ramp_start_step,
            ramp_force_increase_per_step=NS_PER_STEP \
                * self.ramp_force_increase_per_ns,
        )
        #reporter_name = f"{self.name}_swarm_reporter"
        sampling_stage = MD_stage_item(
            name=sampling_name,
            run_minimization=False,
            ensemble=self.ensemble,
            sampling=equil_sampling,
            completion=Number_of_steps_completion_spec(
                number_of_steps=number_of_equil_steps),
            position_reporter_interval=self.starting_structure_interval,
            enforce_periodic_box=False,
            )
        reporter_name = sampling_stage.get_position_reporter_name()
        logistic_stage = Logistic_swarm_from_reporter_stage_item(
            name=assign_swarm_name,
            scope="partitioned",
            reporter_name=reporter_name,
        )
        ramd_stage = MD_stage_item(
            name=ramd_name,
            scope="unpartitioned",
            ensemble=self.ensemble,
            sampling=ramd_sampling,
            completion=CV_value_attained_completion_spec(
                cv_name=self.escape_cv_name, cv_inequality="greater_than", 
                cv_value=self.escape_distance,),
            position_reporter_interval=self.position_reporter_interval,
        )
        return _chain([sampling_stage, logistic_stage, ramd_stage], input_stage_name)

@define
class Composite_stage_procedure(Stage_procedure_base):
    """
    An ordered sequence of sub-procedures. Sub-procedures are expanded in order
    and chained end-to-start.
    """
    type: typing.Literal["composite"] = "composite"
    procedures: list[Stage_procedure_base] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(Stage_procedure_base),
            iterable_validator=validators.instance_of(list),
        ))

    def expand(
            self,
            ctx: typing.Any,
            input_stage_name: str = INITIAL,
            ) -> list[Resolved_stage_item]:
        all_items: list[Resolved_stage_item] = []
        previous_exit = input_stage_name
        for procedure in self.procedures:
            items = procedure.expand(ctx, previous_exit)
            all_items.extend(items)
            if items:
                previous_exit = items[-1][1].name
        return all_items

    def get_emitted_stages(self) -> list[Emitted_stage]:
        stages: list[Emitted_stage] = []
        for procedure in self.procedures:
            stages.extend(procedure.get_emitted_stages())
        return stages


def _insert_unique_stage_maps(
        address_result: dict[str, Stage_address_info],
        policy_result: dict[str, Emitted_stage],
        stage_name: str,
        address_info: Stage_address_info,
        emitted: Emitted_stage,
        ) -> None:
    if stage_name in address_result:
        raise ValueError(
            f"Duplicate stage name {stage_name!r}. Stage names must be unique "
            f"across all procedures.")
    address_result[stage_name] = address_info
    policy_result[stage_name] = emitted


def _build_stage_maps(
        procedure: Stage_procedure_base,
        prefix_path: list[str],
        ) -> tuple[dict[str, Stage_address_info], dict[str, Emitted_stage]]:
    address_result: dict[str, Stage_address_info] = {}
    policy_result: dict[str, Emitted_stage] = {}
    if isinstance(procedure, Composite_stage_procedure):
        for child in procedure.procedures:
            child_address, child_policy = _build_stage_maps(
                child, prefix_path + [child.name])
            for stage_name in child_address:
                _insert_unique_stage_maps(
                    address_result,
                    policy_result,
                    stage_name,
                    child_address[stage_name],
                    child_policy[stage_name],
                )
    else:
        for emitted in procedure.get_emitted_stages():
            role_path = [] if emitted.role is None else [emitted.role]
            _insert_unique_stage_maps(
                address_result,
                policy_result,
                emitted.name,
                (prefix_path + role_path, emitted.role),
                emitted,
            )
    return address_result, policy_result


def _initial_prefix_path(
        procedure: Stage_procedure_base,
        ) -> list[str]:
    if isinstance(procedure, Composite_stage_procedure):
        return []
    return [procedure.name]


def build_stage_address_map(
        procedure: Stage_procedure_base,
        ) -> dict[str, Stage_address_info]:
    """
    Map each produced stage name to ``(address_path, role)``.

    ``address_path`` is ``[procedure_name, (child...)?, (role)?]``. The root
    composite's own name is omitted; a lone non-composite top procedure is
    prefixed with its name.
    """
    address_map, _policy_map = _build_stage_maps(
        procedure, _initial_prefix_path(procedure))
    return address_map


def build_stage_policy_map(
        procedure: Stage_procedure_base,
        ) -> dict[str, Emitted_stage]:
    """
    Map each produced stage name to its intrinsic execution defaults as declared
    by the emitting procedure.
    """
    _address_map, policy_map = _build_stage_maps(
        procedure, _initial_prefix_path(procedure))
    return policy_map


def build_stage_to_procedure_map(
        top_procedure: Stage_procedure_base,
        ) -> dict[str, str]:
    """
    Map each generated stage name to the top-level procedure that authored it.

    Thin wrapper over :func:`build_stage_address_map`: for a composite it takes
    the first address segment (``path[0]``, i.e. the composite's direct child
    procedure), and for a lone non-composite it uses that procedure's name.

    Example (composite of ``equilibration``, ``metadynamics``, ``MMVT``, ``bd``)::

        {"sampling": "equilibration",
         "metadynamics_seed": "metadynamics",
         "metadynamics_assign_structures": "metadynamics",
         "MMVT": "MMVT",
         "bd": "bd"}

    Use it to answer "which authored procedure produced this stage?" without
    caring about the full role-path (see ``build_stage_address_map`` for that).
    """
    if isinstance(top_procedure, Composite_stage_procedure):
        return {
            name: path[0]
            for name, (path, _role) in
            build_stage_address_map(top_procedure).items()
            if len(path) > 0
        }
    return {
        name: top_procedure.name
        for name in build_stage_address_map(top_procedure)
    }


Stage_procedure = (
    Explicit_stage_procedure | Equilibration_globular_stage_procedure
    | Seeding_stage_procedure | MMVT_stage_procedure | BD_stage_procedure
    | Composite_stage_procedure
)
