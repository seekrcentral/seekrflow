"""
tools/make_input.py

Produce an example input object for a seekrflow calculation.
"""

import json
import argparse

from cattrs import unstructure

import seekrflow.modules.base as base
import seekrflow.modules.structures as structures
import seekrflow.modules.workflows.structures as workflow_structures
import seekrflow.modules.parameters_topology_structures as parameters_topology_structures
import seekrflow.modules.workflows.protein_ligand_seekr2.structures as protein_ligand_seekr2_structures

def create_example_seekrflow(ff="amber") -> structures.Seekrflow:
    seekrflow = structures.Seekrflow()
    seekrflow.name = "example_seekrflow"
    workflow = protein_ligand_seekr2_structures.Protein_ligand_seekr2_workflow()
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
    workflow.parameterizer_information = protein_ligand_seekr2_structures.Parameterizer_information()
    workflow.mmvt_settings = protein_ligand_seekr2_structures.MMVT_seekr_settings()
    workflow.hidr_settings = protein_ligand_seekr2_structures.HIDR_settings_metaD()
    workflow.md_settings = workflow_structures.MD_settings()
    workflow.bd_settings = workflow_structures.BD_settings()
    seekrflow.workflow = workflow
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
    seekrflow.run_settings.bd_stage_resource_name = "local"
    seekrflow.run_settings.hidr_stage_resource_name = "delta"
    seekrflow.run_settings.seekr_stage_resource_name = "delta"
    
    return seekrflow

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument(
        "input_file", metavar="INPUT_FILE", type=str, 
        help="The name of the input JSON file to generate.")
    argspace = argparser.parse_args()
    args = vars(argspace)
    json_filename = args["input_file"]
    seekrflow = create_example_seekrflow()
    seekrflow_dict = unstructure(seekrflow)
    json_dump = json.dumps(seekrflow_dict, indent=4)
    with open(json_filename, "w") as file:
        file.write(json_dump)