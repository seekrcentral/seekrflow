"""
tools/make_input.py

Produce an example input object for a seekrflow calculation.
"""

import os
import argparse

import seekrflow.modules.base as base
import seekrflow.modules.structures as structures
import seekrflow.modules.workflows.structures as workflow_structures
import seekrflow.modules.parameters_topology_structures as parameters_topology_structures
import seekrflow.modules.parameterize_workflow as parameterize_workflow_module
import seekrflow.modules.workflows.workflow as workflow_module
import seekrflow.modules.workflows.components as components_module
import seekrflow.modules.workflows.cv_specs as cv_specs_module
import seekrflow.modules.workflows.anchor_specs as anchor_specs_module
import seekrflow.modules.workflows.scale_settings as scale_settings_module
import seekrflow.modules.workflows.stage_procedures as stage_procedures_module


def _hostguest_data_directory() -> str:
    """
    Locate the seekr3 ``hostguest_files`` data directory (parameter, structure
    and PQR files for the host-guest test system).
    """
    import seekr
    return os.path.join(
        os.path.dirname(seekr.__file__), "data", "hostguest_files")

def create_example_seekrflow(ff="amber") -> structures.Seekrflow:
    seekrflow = structures.Seekrflow()
    seekrflow.name = "example_seekrflow"
    workflow = parameterize_workflow_module.Parameterize_workflow()
    workflow.solvated_system_for_md = workflow_structures.Solvated_system_for_md()
    if ff == "amber":
        workflow.solvated_system_for_md.parameters_topology = parameters_topology_structures.Amber_parameters_topology()
        workflow.solvated_system_for_md.parameters_topology.prmtop_filename = ""
    elif ff == "gromacs":
        workflow.solvated_system_for_md.parameters_topology = parameters_topology_structures.Gromacs_parameters_topology()
        workflow.solvated_system_for_md.parameters_topology.gro_filename = ""
        workflow.solvated_system_for_md.parameters_topology.top_filename = ""
    elif ff == "charmm":
        workflow.solvated_system_for_md.parameters_topology = parameters_topology_structures.Charmm_parameters_topology()
        workflow.solvated_system_for_md.parameters_topology.psf_filename = ""
        workflow.solvated_system_for_md.parameters_topology.param_filename_list = []
    elif ff == "forcefield":
        workflow.solvated_system_for_md.parameters_topology = parameters_topology_structures.Forcefield_parameters()
        workflow.solvated_system_for_md.parameters_topology.built_in_forcefield_filenames = []
        workflow.solvated_system_for_md.parameters_topology.custom_forcefield_filenames = []
    elif ff == "openmm":
        workflow.solvated_system_for_md.parameters_topology = parameters_topology_structures.Openmm_system()
        workflow.solvated_system_for_md.parameters_topology.system_filename = ""
    else:
        raise ValueError(f"Unrecognized forcefield type: {ff}")

    workflow.solvated_system_for_md.solvated_pdb = ""
    workflow.ligand_indices = []
    workflow.receptor_indices = []
    workflow.receptor_pqr_filename_for_bd = ""
    workflow.ligand_pqr_filename_for_bd = ""
    workflow.parameterizer_information = parameterize_workflow_module.Parameterizer_information()
    workflow.md_settings = workflow_structures.MD_settings()
    workflow.bd_settings = workflow_structures.BD_settings()
    seekrflow.parameterize_workflow = workflow
    physical_attributes = base.Physical_attributes()
    physical_attributes.temperature = 298.15
    physical_attributes.hydrogen_mass = 3.0
    seekrflow.physical_attributes = physical_attributes
    seekrflow.run_settings = structures.Run_settings()
    delta_slurm_resource = structures.Resource_remote_slurm()
    delta_slurm_resource.name = "delta"
    delta_slurm_resource.remote_working_directory = "/scratch/kif/USERNAME/seekrflow_playground/"
    delta_slurm_resource.partition = "gpuA100x4,gpuA40x4"
    delta_slurm_resource.account = "kif-delta-gpu"
    delta_slurm_resource.constraint = "scratch"
    delta_slurm_resource.nodes_per_block = 1
    delta_slurm_resource.cpus_per_task = 32
    delta_slurm_resource.memory_per_node = 48000
    delta_slurm_resource.time_limit = "00:30:00"
    delta_slurm_resource.scheduler_options = "#SBATCH --gpus-per-node=1 --gpu-bind=closest"
    delta_slurm_resource.worker_init = "source $HOME/.bashrc; "\
                            "conda activate SEEKR2; "\
                            "export OPENMM_CUDA_COMPILER=`which nvcc`"
    delta_slurm_resource.remote_interface = structures.Remote_interface_globus_compute_sdk()
    delta_slurm_resource.remote_interface.endpoint_id = "MY_ENDPOINT_ID"
    delta_slurm_resource.transfer_settings = structures.Transfer_settings_globus()
    delta_slurm_resource.transfer_settings.local_collection_id = "MY_LOCAL_COLLECTION_ID"
    delta_slurm_resource.transfer_settings.remote_collection_id = "MY_REMOTE_COLLECTION_ID"
    seekrflow.run_settings.resources = [delta_slurm_resource]
    #anvil_slurm_resource = structures.Slurm_resource()
    seekrflow.run_settings.placements = [
        structures.Placement(target=[], resource="local"),
        structures.Placement(target=["equilibration"], resource="local"),
        structures.Placement(
            target=["metadynamics", "sampling"], resource="delta"),
        structures.Placement(
            target=["metadynamics", "logistic"], resource="local"),
        structures.Placement(target=["MMVT"], resource="delta"),
        structures.Placement(target=["bd"], resource="local"),
    ]
    
    return seekrflow

def create_host_guest_example_seekrflow(
        work_directory: str,
        provide_pqr_filenames: bool = True
        ) -> structures.Seekrflow:
    """
    Build a complete, runnable seekrflow input for the beta-cyclodextrin /
    1-butanol host-guest test system, using the composable ``Workflow``.

    Mirrors the seekr3 host-guest example (``seekr/tools/make_input.py``):
    a distance collective variable between the host (resname MGO) and guest
    (resname APN), spherical anchors from 0.05 to 1.35 nm, an
    equilibration -> metadynamics-seeding -> MMVT -> Brownian-dynamics
    procedure, and two Brownian-dynamics rigid bodies (one per molecule).
    """
    hostguest_directory = _hostguest_data_directory()
    prmtop_filename = os.path.join(hostguest_directory, "hostguest.parm7")
    solvated_pdb_filename = os.path.join(
        hostguest_directory, "hostguest_at0.5.pdb")
    receptor_pqr_filename = os.path.join(
        hostguest_directory, "hostguest_receptor.pqr")
    ligand_pqr_filename = os.path.join(
        hostguest_directory, "hostguest_ligand.pqr")

    # Components: the host and the guest, each a small molecule selected by
    #  residue name. Classifying them as small molecules means their
    #  Brownian-dynamics PQRs are automatically emitted with one residue per
    #  atom.
    receptor_component = components_module.Small_molecule_component(
        name="host",
        selector=components_module.Selector_by_resname(resname="MGO"))
    ligand_component = components_module.Small_molecule_component(
        name="guest",
        selector=components_module.Selector_by_resname(resname="APN"))
    components = components_module.Components(
        members=[receptor_component, ligand_component])

    # Collective variable: host-guest center-of-mass distance.
    cv_spec = cv_specs_module.Com_com_distance_CV_spec(
        name="distance_host_guest",
        group1_selection_name="host",
        group2_selection_name="guest",
        min_value=0.0,
        max_value=1.4)

    # Anchors: 14 evenly-spaced anchors spanning the collective-variable range
    #  (0.0 - 1.4 nm), i.e. cell centers at 0.05, 0.15, ..., 1.35 nm. The bound
    #  state (the anchor containing the starting structure) and the bulk state
    #  (the outermost anchor) are determined automatically during prepare.
    anchor_spec = anchor_specs_module.Uniform_anchor_spec(n_anchors=14)

    # Procedure: equilibrate -> seed anchors by metadynamics -> MMVT -> BD.
    """
    equilibration = stage_procedures_module.Explicit_stage_procedure(
        name="equilibration",
        items=[stage_procedures_module.MD_stage_item(
            name="equilibration_sampling",
            run_minimization=True,
            ensemble="NVT",
            sampling=stage_procedures_module.Conventional_sampling_spec(
                heat_ramp_temperatures=[100.0, 298.15],
                heat_ramp_step_interval=1000),
            completion=stage_procedures_module
                .Number_of_steps_completion_spec(number_of_steps=100000),
        )])
    """
    equilibration = stage_procedures_module.Equilibration_globular_stage_procedure(
        name="equilibration",
        reference_structure_filename=solvated_pdb_filename,
        align_selection_str="resname MGO and not type H"
    )
    metad_method_input = stage_procedures_module.Metadynamics_seeding_method_input(
        bias_factor=10.0,
        gaussian_height=0.1,
        gaussian_width=0.1,
        steps_per_update=250,
        number_of_points=101)
    seeding = stage_procedures_module.Seeding_stage_procedure(
        name="metadynamics",
        #method="metadynamics",
        method_input=metad_method_input,
        cv_names=["distance_host_guest"])
    mmvt = stage_procedures_module.MMVT_stage_procedure(
        name="MMVT",
        ensemble="NVT",
        number_of_steps=250000)
    bd = stage_procedures_module.BD_stage_procedure(
        name="bd",
        n_trajectories_per_output=1000,
        n_steps_per_output=100000,
        number_of_trajectories=10000)
    procedure = stage_procedures_module.Composite_stage_procedure(
        name="full_procedure",
        procedures=[equilibration, seeding, mmvt, bd])

    # Per-scale settings: molecular dynamics + Brownian dynamics.
    md_scale = scale_settings_module.Molecular_dynamics_scale_settings()
    md_scale.system.parameters_topology = \
        parameters_topology_structures.Amber_parameters_topology()
    md_scale.system.parameters_topology.prmtop_filename = prmtop_filename
    md_scale.system.solvated_pdb = solvated_pdb_filename
    md_scale.hydrogen_mass = 3.0
    md_scale.timestep = 0.004
    md_scale.platform_type = "cuda"

    bd_scale = scale_settings_module.Brownian_dynamics_scale_settings()
    bd_scale.engine = "browndye"
    bd_scale.binary_directory = ""
    bd_scale.system.molecules = [
        scale_settings_module.BD_molecule(
            name="beta-cyclodextrin",
            component_name="host",
            pqr_filename=receptor_pqr_filename if provide_pqr_filenames else "",
            role="receptor"),
        scale_settings_module.BD_molecule(
            name="1-butanol",
            component_name="guest",
            pqr_filename=ligand_pqr_filename if provide_pqr_filenames else "",
            role="ligand"),
    ]

    workflow = workflow_module.Workflow(
        components=components,
        cv_specs=[cv_spec],
        anchor_spec=anchor_spec,
        procedure=procedure,
        scale_settings=[md_scale, bd_scale])

    seekrflow = structures.Seekrflow()
    seekrflow.name = "host_guest_example"
    seekrflow.workflow = workflow
    # The host-guest system is pre-parameterized, so no parameterize step.
    seekrflow.parameterize_workflow = None
    physical_attributes = base.Physical_attributes()
    physical_attributes.temperature = 298.15
    physical_attributes.ionic_strength = 0.0
    physical_attributes.random_seed = 111
    seekrflow.physical_attributes = physical_attributes
    seekrflow.work_directory = work_directory
    seekrflow.run_settings = structures.Run_settings()
    mmvt_dispatch = stage_procedures_module.Dispatch(
        dimensions=["anchor"], group_size=1, concurrency=1)
    seekrflow.run_settings.placements = [
        structures.Placement(target=[], resource="local"),
        structures.Placement(target=["equilibration"], resource="panamint"),
        structures.Placement(
            target=["metadynamics", "sampling"], resource="panamint"),
        structures.Placement(
            target=["metadynamics", "logistic"], resource="panamint", 
            co_schedule_with="predecessor"),
        structures.Placement(target=["MMVT"], resource="panamint", dispatch=mmvt_dispatch),
        structures.Placement(target=["bd"], resource="local"),
    ]
    return seekrflow

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument(
        "input_file", metavar="INPUT_FILE", type=str, 
        help="The name of the input JSON file to generate.")
    argparser.add_argument(
        "-s", "--system", dest="system", type=str, default="parameterize",
        choices=["parameterize", "hostguest"],
        help="Which example input to generate: 'parameterize' (default) for "
             "the parameterize-side example, or 'hostguest' for the complete "
             "prepared host-guest workflow.")
    argspace = argparser.parse_args()
    args = vars(argspace)
    json_filename = args["input_file"]
    if args["system"] == "hostguest":
        provide_pqr_filenames = False
        work_directory = "~/tmp/seekrflow_hostguest_example/"
        seekrflow = create_host_guest_example_seekrflow(
            provide_pqr_filenames=provide_pqr_filenames,
            work_directory=work_directory)
    else:
        seekrflow = create_example_seekrflow()
    seekrflow.save(json_filename)