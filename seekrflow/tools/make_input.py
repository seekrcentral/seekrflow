"""
tools/make_input.py

Produce an example input object for a seekrflow calculation.
"""

import json
import argparse

from cattrs import unstructure

import seekrflow.modules.structures as structures

def create_example_seekrflow() -> structures.Seekrflow:
    seekrflow = structures.Seekrflow()
    seekrflow.seekr_settings = structures.MMVT_seekr_settings()
    seekrflow.hidr_settings = structures.HIDR_settings_metaD()
    seekrflow.run_settings = structures.Run_settings()
    delta_slurm_resource = structures.Resource_remote_slurm()
    delta_slurm_resource.name = "delta"
    delta_slurm_resource.remote_working_directory = "/scratch/kif/lvotapka/seekrflow_playground/"
    delta_slurm_resource.partition = "gpuA100x4,gpuA40x4"
    delta_slurm_resource.account = "kif-delta-gpu"
    delta_slurm_resource.constraint = "scratch"
    delta_slurm_resource.cores_per_node = 32
    delta_slurm_resource.memory_per_node = "220"
    delta_slurm_resource.time_limit = "00:30:00"
    delta_slurm_resource.scheduler_options = "#SBATCH --gpus-per-node=1 --gpu-bind=closest"
    delta_slurm_resource.worker_init = "source $HOME/.bashrc; "\
                            "conda activate SEEKR2; "\
                            "export OPENMM_CUDA_COMPILER=`which nvcc`"
    delta_slurm_resource.globus_compute_endpoint_id = "29f14463-e6e1-4cd0-b25d-74f63aa10ca7"
    delta_slurm_resource.transfer_settings = structures.Transfer_settings_globus()
    delta_slurm_resource.transfer_settings.local_collection_id = "3896d182-5123-11f0-8180-0affcfc1d1e5"
    delta_slurm_resource.transfer_settings.remote_collection_id = "7e936164-de58-4e3d-85da-21aa23c07169"
    seekrflow.run_settings.resources = [delta_slurm_resource]
    #anvil_slurm_resource = structures.Slurm_resource()
    seekrflow.run_settings.bd_stage_resource_name = "local"
    seekrflow.run_settings.hidr_stage_resource_name = "delta"
    seekrflow.run_settings.seekr_stage_resource_name = "delta"
    seekrflow.run_settings.allow_parsl_usage_tracking = False
    
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