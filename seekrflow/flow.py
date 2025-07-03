"""
flow.py

Manage all stages of seekrflow (except for parametrization). This 
includes preparation and running.
"""

import os
import typing
import pathlib
import argparse

import seekrflow.modules.base as base
import seekrflow.modules.structures as structures
import seekrflow.modules.seekr_input as seekr_input
import seekrflow.modules.seekr_run as seekr_run

PREPARE_SEEKRFLOW_GLOB = "seekrflow_pre_prepare_*.json"
PREPARE_SEEKRFLOW_BASE = "seekrflow_pre_prepare_{}.json"
RUN_SEEKRFLOW_GLOB = "seekrflow_pre_run_*.json"
RUN_SEEKRFLOW_BASE = "seekrflow_pre_run_{}.json"

def flow(
        seekrflow: structures.Seekrflow, 
        instruction: str
        ) -> None | typing.Tuple[str, str]:
    """
    Execute the instructed seekrflow stage.
    """
    if instruction == "any":
        # Automatically perform the next step
        print("Preparing system...")
        seekr_input.prepare_model(seekrflow)
        seekrflow_glob = PREPARE_SEEKRFLOW_GLOB
        seekrflow_base = PREPARE_SEEKRFLOW_BASE
        # Run
        print("Running system...")
        seekr_run.run_model(seekrflow)
        
        return seekrflow_glob, seekrflow_base
        
    elif instruction == "prepare":
        # Prepare the system
        print("Preparing system...")
        seekr_input.prepare_model(seekrflow)
        return PREPARE_SEEKRFLOW_GLOB, PREPARE_SEEKRFLOW_BASE
        
    elif instruction == "run":
        # Run the system
        print("Running system...")
        seekr_run.run_model(seekrflow)
        return None, None
        
    else:
        raise ValueError(
            f"Invalid instruction '{instruction}'. Options are: 'any', 'prepare', 'run'.")
    

def main():
    argparser = argparse.ArgumentParser(
        description="Automates the preparation and running of SEEKR calculations"
        "for particular purposes, such as a ligand-receptor calculation.")
    argparser.add_argument(
        "input_json", metavar="INPUT_JSON", type=str, 
        help="Path to the input JSON file containing the parameters for seekrflow.")
    argparser.add_argument(
        "instruction", metavar="INSTRUCTION", type=str, 
        help="The instruction for what step/stage to perform. Options include: "\
            "'any', which will automatically perform the next step, 'prepare', "\
            "which will prepare the system, 'run', which will run the system.")
    argparser.add_argument(
        "-w", "--work_directory", dest="work_directory",
        metavar="WORK_DIRECTORY", type=str, default=None,
        help="Path to the work directory for the parametrization.")
    args = argparser.parse_args()
    args = vars(args)
    input_json = pathlib.Path(args["input_json"])
    instruction = args["instruction"]
    
    assert input_json.exists(), \
        f"Input JSON file {input_json} does not exist."
    seekrflow = structures.load_seekrflow(input_json)
    if args["work_directory"] is None:
        work_dir = pathlib.Path(seekrflow.work_directory)
    else:
        work_dir = pathlib.Path(args["work_directory"])
    seekrflow.make_work_directory(work_dir)
    curdir = os.getcwd()
    os.chdir(work_dir)
    if len(seekrflow.ligand_indices) == 0:
        if seekrflow.ligand_resname != "":
            seekrflow.ligand_indices = base.get_ligand_indices(seekrflow.starting_pdb_filename, 
                                                               seekrflow.ligand_resname)
        else:
            seekrflow.ligand_indices = []
    os.chdir(curdir)
    seekrflow_glob, seekrflow_base = flow(seekrflow, instruction)
    if seekrflow_glob is not None and seekrflow_base is not None:
        seekrflow.work_directory = str(work_dir)
        structures.save_new_seekrflow(seekrflow, seekrflow_glob, seekrflow_base,
                                    save_old_seekrflow=True, directory=seekrflow.work_directory)

if __name__ == "__main__":
    main()