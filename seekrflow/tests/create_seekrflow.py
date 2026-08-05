"""
tests/create_seekrflow.py

Produce test seekflow objects
"""

import os

import seekrflow.modules.structures as structures
import seekrflow.modules.parameterize_structures as parameterizer_structures
import seekrflow.modules.parameters_topology_structures \
    as parameters_topology_structures
import seekrflow.modules.workflows.workflow as workflow_module
import seekrflow.modules.workflows.components as components_module
import seekrflow.modules.workflows.cv_specs as cv_specs_module
import seekrflow.modules.workflows.anchor_specs as anchor_specs_module
import seekrflow.modules.workflows.scale_settings as scale_settings_module
import seekrflow.modules.workflows.stage_procedures as stage_procedures_module

TEST_DIRECTORY = os.path.dirname(__file__)


def _minimal_prepare_procedure(
        reference_structure_filename: str,
        align_selection_str: str = "protein and not type H",
        cv_name: str = "distance_receptor_ligand",
        ) -> stage_procedures_module.Composite_stage_procedure:
    """
    A minimal equilibration -> metadynamics seeding -> MMVT -> BD procedure
    sufficient for prepare tests.
    """
    equilibration = stage_procedures_module.Equilibration_globular_stage_procedure(
        name="equilibration",
        reference_structure_filename=reference_structure_filename,
        align_selection_str=align_selection_str,
    )
    metad_method_input = stage_procedures_module.Metadynamics_seeding_method_input(
        bias_factor=10.0,
        gaussian_height=0.1,
        gaussian_width=0.1,
        steps_per_update=250,
        number_of_points=101)
    seeding = stage_procedures_module.Seeding_stage_procedure(
        name="metadynamics",
        method_input=metad_method_input,
        cv_names=[cv_name])
    mmvt = stage_procedures_module.MMVT_stage_procedure(
        name="MMVT",
        ensemble="NVT",
        number_of_steps=1000)
    bd = stage_procedures_module.BD_stage_procedure(
        name="bd",
        n_trajectories_per_output=1000,
        n_steps_per_output=100000,
        number_of_trajectories=100)
    return stage_procedures_module.Composite_stage_procedure(
        name="full_procedure",
        procedures=[equilibration, seeding, mmvt, bd])


def _make_protein_ligand_components(
        ligand_resname: str,
        ligand_indices: list[int] | None = None,
        sdf_file: str = "",
        ) -> components_module.Components:
    if ligand_indices is not None and len(ligand_indices) > 0:
        ligand_selector = components_module.Selector_by_indices(
            indices=list(ligand_indices))
    else:
        ligand_selector = components_module.Selector_by_resname(
            resname=ligand_resname, include_hydrogens=True)
    protein = components_module.Protein_component(
        name="receptor",
        selector=components_module.Selector_mdtraj(selection_string="protein"),
    )
    ligand = components_module.Small_molecule_component(
        name="ligand",
        selector=ligand_selector,
        sdf_file=sdf_file,
        resname=ligand_resname,
    )
    return components_module.Components(members=[protein, ligand])


def create_unparameterized_seekrflow(
        starting_structure: str,
        ligand_resname: str,
        anchor_radius_list: list[float],
        ff: str = "amber",
        ) -> structures.Seekrflow:
    """
    Create a testing seekrflow object before the parameterization step.
    """
    seekrflow = structures.Seekrflow()
    seekrflow.name = "test_seekrflow"
    components = _make_protein_ligand_components(ligand_resname)
    n_anchors = len(anchor_radius_list)
    cv_spec = cv_specs_module.Com_com_distance_CV_spec(
        name="distance_receptor_ligand",
        group1_selection_name="receptor",
        group2_selection_name="ligand",
        min_value=min(anchor_radius_list) if anchor_radius_list else 0.0,
        max_value=max(anchor_radius_list) if anchor_radius_list else 1.0,
    )
    md_scale = scale_settings_module.Molecular_dynamics_scale_settings()
    # Blank system: parameterize will produce it.
    md_scale.system = None
    md_scale.hydrogen_mass = 1.5
    workflow = workflow_module.Workflow(
        components=components,
        cv_specs=[cv_spec],
        anchor_spec=anchor_specs_module.Uniform_anchor_spec(n_anchors=n_anchors),
        procedure=_minimal_prepare_procedure(
            reference_structure_filename=starting_structure),
        scale_settings=[md_scale],
    )
    seekrflow.workflow = workflow
    parameterizer = parameterizer_structures.Parameterizer()
    parameterizer.complex_pdb_filename = starting_structure
    parameterizer.water_model = "tip3p"
    if ff == "amber":
        parameterizer.small_molecule_forcefield = "gaff-2.11"
        parameterizer.forcefields = [
            "amber/ff14SB.xml",
            "amber/tip3p_standard.xml",
            "amber/tip3p_HFE_multivalent.xml",
        ]
    else:
        parameterizer.small_molecule_forcefield = ff
        parameterizer.forcefields = [
            "amber/ff14SB.xml",
            "amber/tip3p_standard.xml",
            "amber/tip3p_HFE_multivalent.xml",
        ]
    seekrflow.parameterizer = parameterizer
    seekrflow.run_settings = structures.Run_settings()
    seekrflow.run_settings.placements = []
    return seekrflow


def create_parameterized_seekrflow(
        solvated_pdb: str,
        parameters_topology: parameters_topology_structures.Parameters_topology,
        receptor_pqr_filename: str,
        ligand_pqr_filename: str,
        ligand_resname: str,
        anchor_radius_list: list[float],
        ff: str = "amber",
        ligand_indices: list[int] | None = None,
        ) -> structures.Seekrflow:
    """
    Create a parameterized seekrflow object ready for the prepare stage.
    """
    seekrflow = create_unparameterized_seekrflow(
        starting_structure=solvated_pdb,
        ligand_resname=ligand_resname,
        anchor_radius_list=anchor_radius_list,
        ff=ff,
    )
    # Already parameterized: clear parameterizer and set MD system explicitly.
    seekrflow.parameterizer = None
    if ligand_indices is not None:
        seekrflow.workflow.components = _make_protein_ligand_components(
            ligand_resname, ligand_indices=ligand_indices)
    md_scale = seekrflow.workflow.get_md_settings()
    md_scale.system = scale_settings_module.System_for_md(
        parameters_topology=parameters_topology,
        solvated_pdb=solvated_pdb,
    )
    bd_scale = scale_settings_module.Brownian_dynamics_scale_settings()
    bd_scale.system.molecules = [
        scale_settings_module.BD_molecule(
            name="receptor",
            component_name="receptor",
            pqr_filename=receptor_pqr_filename,
            role="receptor"),
        scale_settings_module.BD_molecule(
            name="ligand",
            component_name="ligand",
            pqr_filename=ligand_pqr_filename,
            role="ligand"),
    ]
    seekrflow.workflow.procedure = _minimal_prepare_procedure(
        reference_structure_filename=solvated_pdb)
    seekrflow.workflow.scale_settings = [md_scale, bd_scale]
    return seekrflow


def create_unparameterized_tryp_ben_seekrflow(
        ff: str = "amber",
        ) -> structures.Seekrflow:
    """
    Create a seekrflow object for the Tryp-Ben system without ff parameters.
    """
    starting_structure = os.path.join(
        TEST_DIRECTORY, "data", "trypsin_benzamidine.pdb")
    ligand_resname = "BEN"
    anchor_radius_list = [
        0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
    return create_unparameterized_seekrflow(
        starting_structure=starting_structure,
        ligand_resname=ligand_resname,
        anchor_radius_list=anchor_radius_list,
        ff=ff,
    )


def create_parameterized_tryp_ben_openmm_xml_seekrflow(
        ff: str = "amber",
        ) -> structures.Seekrflow:
    """
    Create a seekrflow object for the Tryp-Ben system with ff parameters.
    """
    starting_pdb = os.path.join(
        TEST_DIRECTORY, "data", "tryp_ben_system_bound.pdb")
    starting_system_xml = os.path.join(
        TEST_DIRECTORY, "data", "tryp_ben_system.xml")
    receptor_pqr_filename = os.path.join(
        TEST_DIRECTORY, "data", "tryp_ben_receptor.pqr")
    ligand_pqr_filename = os.path.join(
        TEST_DIRECTORY, "data", "tryp_ben_ligand.pqr")
    md_parameters_topology = parameters_topology_structures.Openmm_system()
    md_parameters_topology.system_filename = starting_system_xml
    ligand_resname = "BEN"
    anchor_radius_list = [
        0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
    return create_parameterized_seekrflow(
        solvated_pdb=starting_pdb,
        parameters_topology=md_parameters_topology,
        receptor_pqr_filename=receptor_pqr_filename,
        ligand_pqr_filename=ligand_pqr_filename,
        ligand_resname=ligand_resname,
        anchor_radius_list=anchor_radius_list,
        ff=ff,
    )


def create_parameterized_host_guest_amber_seekrflow(
        ff: str = "amber",
        ) -> structures.Seekrflow:
    """
    Create a seekrflow object for the host-guest system with Amber parameters.
    """
    starting_pdb = os.path.join(TEST_DIRECTORY, "data", "hostguest_at0.5.pdb")
    starting_system_parm7 = os.path.join(
        TEST_DIRECTORY, "data", "hostguest.parm7")
    receptor_pqr_filename = os.path.join(
        TEST_DIRECTORY, "data", "hostguest_receptor.pqr")
    ligand_pqr_filename = os.path.join(
        TEST_DIRECTORY, "data", "hostguest_ligand.pqr")
    md_parameters_topology = \
        parameters_topology_structures.Amber_parameters_topology()
    md_parameters_topology.prmtop_filename = starting_system_parm7
    ligand_resname = "APN"
    anchor_radius_list = [
        0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
    # Host-guest: both are small molecules in the full example; for this
    # fixture keep a protein-style receptor selector by indices for BD tests.
    flow = create_parameterized_seekrflow(
        solvated_pdb=starting_pdb,
        parameters_topology=md_parameters_topology,
        receptor_pqr_filename=receptor_pqr_filename,
        ligand_pqr_filename=ligand_pqr_filename,
        ligand_resname=ligand_resname,
        anchor_radius_list=anchor_radius_list,
        ff=ff,
        ligand_indices=list(range(147, 162)),
    )
    # Replace receptor with index-based host selection.
    host = components_module.Small_molecule_component(
        name="receptor",
        selector=components_module.Selector_by_indices(
            indices=list(range(147))),
        resname="MGO",
    )
    guest = components_module.Small_molecule_component(
        name="ligand",
        selector=components_module.Selector_by_indices(
            indices=list(range(147, 162))),
        resname="APN",
    )
    flow.workflow.components.members = [host, guest]
    flow.workflow.procedure = _minimal_prepare_procedure(
        reference_structure_filename=starting_pdb,
        align_selection_str="resname MGO and not type H",
        cv_name="distance_receptor_ligand",
    )
    flow.physical_attributes.ionic_strength = 0.0
    return flow
