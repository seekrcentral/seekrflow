"""
flow.py

Manage all stages of seekrflow (except for parameterization). This 
includes preparation and running.
"""

import os
import argparse

import seekrflow.modules.parameters_topology_structures as parameters_topology_structures
import seekrflow.modules.structures as structures
import seekrflow.modules.workflows.structures as workflow_structures
import seekrflow.modules.seekr_input as seekr_input
import seekrflow.modules.seekr_run as seekr_run

def assign_seekrflow_defaults(
        seekrflow: structures.Seekrflow,
        pdb_with_system: str | None = None
        ) -> None:
    """
    Assign default values to any missing seekrflow parameters.
    """
    if seekrflow.workflow.solvated_system_for_md is None:
        starting_pdb_basename = "complex-equil.pdb"
        seekrflow.workflow.solvated_system_for_md = \
            workflow_structures.Solvated_system_for_md()
        src_pdb_filename = os.path.join(structures.PARAMETERIZE, starting_pdb_basename)
        seekrflow.workflow.solvated_system_for_md.solvated_pdb = src_pdb_filename
        seekrflow.workflow.solvated_system_for_md.parameters_topology \
            = parameters_topology_structures.Openmm_system()
        starting_system_basename = "complex.xml"
        src_system_filename = os.path.join(
            structures.PARAMETERIZE, starting_system_basename)
        seekrflow.workflow.solvated_system_for_md.parameters_topology\
            .system_filename = src_system_filename
    elif seekrflow.workflow.solvated_system_for_md.solvated_pdb == "":
        starting_pdb_basename = "complex-equil.pdb"
        src_pdb_filename = os.path.join(structures.PARAMETERIZE, starting_pdb_basename)
        seekrflow.workflow.solvated_system_for_md.solvated_pdb = src_pdb_filename
        seekrflow.workflow.solvated_system_for_md.parameters_topology \
            = parameters_topology_structures.Openmm_system()
        starting_system_basename = "complex.xml"
        src_system_filename = os.path.join(
            structures.PARAMETERIZE, starting_system_basename)
        seekrflow.workflow.solvated_system_for_md.parameters_topology\
            .system_filename = src_system_filename
    else:
        if pdb_with_system is not None:
            seekrflow.workflow.solvated_system_for_md.solvated_pdb = pdb_with_system
        src_pdb_filename = seekrflow.workflow.solvated_system_for_md.solvated_pdb
    return src_pdb_filename

def assign_seekrflow_pdb_with_system(
        seekrflow: structures.Seekrflow,
        pdb_with_system: str
        ) -> None:
    """
    Assign the provided PDB file with the solvated system to the seekrflow
    structure.
    """
    seekrflow.workflow.solvated_system_for_md.solvated_pdb = pdb_with_system
    assert os.path.exists(pdb_with_system), \
        f"Expected PDB file {pdb_with_system} to exist."

def flow(
        seekrflow: structures.Seekrflow, 
        instruction: str,
        src_pdb_filename: str | None = None,
        transfer_before: str | None = None,
        transfer_from_host_only: str | None = None,
        force_rerun: list[str] | None = None,
        benchmark_mode: bool = False,
        bd_resource_name: str | None = None,
        hidr_resource_name: str | None = None,
        seekr_resource_name: str | None = None
        ) -> None:
    """
    Execute the instructed seekrflow stage.
    """
    if instruction == "any":
        # Automatically perform the next step
        print("Preparing system...")
        seekr_input.prepare_model(seekrflow, src_pdb_filename)
        # Run
        print("Running system...")
        seekr_run.run_model(
            seekrflow, transfer_before, transfer_from_host_only, force_rerun,
            benchmark_mode, bd_resource_name, hidr_resource_name, seekr_resource_name)

        return
        
    elif instruction == "prepare":
        # Prepare the system
        print("Preparing system...")
        seekr_input.prepare_model(seekrflow, src_pdb_filename)
        return
        
    elif instruction == "run":
        # Run the system
        print("Running system...")
        seekr_run.run_model(
            seekrflow, transfer_before, transfer_from_host_only, force_rerun,
            benchmark_mode, bd_resource_name, hidr_resource_name, seekr_resource_name)
        return
        
    else:
        raise ValueError(
            f"Invalid instruction '{instruction}'. Options are: 'any', "
            f"'prepare', 'run'.")
    
def main():
    argparser = argparse.ArgumentParser(
        description="Automates the preparation and running of seekr calculations"
        "for particular purposes, such as a ligand-receptor calculation.")
    argparser.add_argument(
        "instruction", metavar="INSTRUCTION", type=str, 
        help="The instruction for what step/stage to perform. Options include: "\
        "'any', which will automatically perform the next step, 'prepare', "\
        "which will prepare the system, 'run', which will run the system.")
    argparser.add_argument(
        "-i", "--input_json", dest="input_json",
        metavar="INPUT_JSON", type=str, default="",
        help="Path to the input JSON file containing the seekrflow "\
        "configuration file. If this is provided, then any other "\
        "command-line arguments will override the values in the JSON file.")
    argparser.add_argument(
        "-m", "--model_filename", dest="model_filename",
        metavar="MODEL_FILENAME", type=str, default="",
        help="Path to the model file to use for the simulation. This activates the "\
        "so-called 'hotshot mode', where an existing model.xml file can be run "\
        "easily without needing to prepare a full seekrflow structure. Note that "\
        "if this is provided, the seekrflow JSON file must not be provided. "\
        "Also, if HIDR will be run, then -p must be provided to specify the "\
        "PDB file with the solvated system.")
    argparser.add_argument(
        "-p", "--pdb_with_system", dest="pdb_with_system", 
        metavar="PDB_WITH_SYSTEM", type=str, default="",
        help="Path to the input PDB file that contains the solvated molecules.")
    #argparser.add_argument(
    #    "-P", "--parameter_topology_files", dest="parameter_topology_files", 
    #    metavar="PARAMETER_TOPOLOGY_FILES", type=str, nargs="+", default=None,
    #    help="List of parameter and topology files to use for the prepare stage.")
    argparser.add_argument(
        "-w", "--work_directory", dest="work_directory",
        metavar="WORK_DIRECTORY", type=str, default=None,
        help="Path to the work directory for the parameterization.")
    argparser.add_argument(
        "-t", "--transfer_before", dest="transfer_before",
        metavar="STAGE", type=str, default=None,
        help="If set, transfer files to remote host for a particular stage. " \
        "This is done automatically if the model.xml file does not already "\
        "exist on the remote host. However, if you have already transferred "\
        "files and want to re-transfer them, use this flag. Default: None, "\
        "which will not transfer files. Other values might be: 'bd', 'hidr', " \
        "'seekr', 'all'.")
    argparser.add_argument(
        "-T", "--transfer_from_host_only", dest="transfer_from_host_only", 
        metavar="STAGE", type=str, default=None,
        help="If set, transfer files from the remote host for a particular stage "\
        "without doing any preparation or running. Default: None, "\
        "which will not transfer files. Other values might be: 'bd', 'hidr', " \
        "'seekr', 'all'.")
    argparser.add_argument(
        "-f", "--force_rerun", dest="force_rerun", metavar="STAGE_LIST",
        type=str, nargs="+", default=None,
        help="If set, activate force rerun for particular stages, listed "\
        "as a space-separated list.")
    argparser.add_argument(
        "-b", "--benchmark_mode", dest="benchmark_mode", 
        help="If set, run in benchmark mode - quick seekr run on resource "\
        "to get approximate timings in ns/day. Default: False.", action="store_true",
        default=False)
    argparser.add_argument(
        "-B", "--Brownian_dynamics_resource_name", dest="bd_resource_name", 
        metavar="BD_RESOURCE_NAME", type=str, default=None,
        help="Name of the resource to run Brownian dynamics on. The configuration" \
        "for the resource must be provided in the seekrflow resources JSON file or" \
        "in a .seekrflow_resources.json file in the home, working, or root directory." \
        "If not provided, defaults to the values in the seekrflow JSON file or 'local'."\
        "Default: None.")
    argparser.add_argument(
        "-H", "--HIDR_resource_name", dest="hidr_resource_name", 
        metavar="HIDR_RESOURCE_NAME", type=str, default=None,
        help="Name of the resource to run HIDR on. The configuration" \
        "for the resource must be provided in the seekrflow resources JSON file or" \
        "in a .seekrflow_resources.json file in the home, working, or root directory." \
        "If not provided, defaults to the values in the seekrflow JSON file or 'local'."\
        "Default: None.")
    argparser.add_argument(
        "-S", "--SEEKR_resource_name", dest="seekr_resource_name", 
        metavar="SEEKR_RESOURCE_NAME", type=str, default=None,
        help="Name of the resource to run SEEKR on. The configuration" \
        "for the resource must be provided in the seekrflow resources JSON file or" \
        "in a .seekrflow_resources.json file in the home, working, or root directory." \
        "If not provided, defaults to the values in the seekrflow JSON file or 'local'."\
        "Default: None.")

    args = argparser.parse_args()
    args = vars(args)
    instruction = args["instruction"]
    input_json = args["input_json"]
    model_filename = args["model_filename"]
    pdb_with_system = args["pdb_with_system"]
    #parameter_topology_files = args["parameter_topology_files"]
    work_dir = args["work_directory"]
    transfer_before = args["transfer_before"]
    transfer_from_host_only = args["transfer_from_host_only"]
    force_rerun = args["force_rerun"]
    benchmark_mode = args["benchmark_mode"]
    bd_resource_name = args["bd_resource_name"]
    hidr_resource_name = args["hidr_resource_name"]
    seekr_resource_name = args["seekr_resource_name"]

    hotshot_mode = False
    src_pdb_filename = None
    if input_json == "":
        # Create default seekrflow structure
        seekrflow = structures.Seekrflow()
        if model_filename != "":
            # Hotshot mode - load existing model
            hotshot_mode = True
            assert instruction == "run", \
                "In hotshot mode, only 'run' instruction is allowed."
            assert os.path.exists(model_filename), \
                f"Model file {model_filename} does not exist."
            seekrflow.name = "hotshot"
            seekrflow.work_directory = None
            seekrflow.root_directory = os.path.dirname(
                os.path.abspath(model_filename))
            structures.try_to_load_resources_json(seekrflow)
            seekrflow.workflow.solvated_system_for_md = \
                workflow_structures.Solvated_system_for_md()
            seekrflow.workflow.solvated_system_for_md.solvated_pdb \
                = pdb_with_system
    else:
        assert model_filename == "", \
            "If input JSON file is provided, model filename must not be provided."
        assert os.path.exists(input_json), \
            f"Input JSON file {input_json} does not exist."
        seekrflow = structures.load_seekrflow(input_json)
        
    if not hotshot_mode:
        if work_dir is None:
            work_dir = seekrflow.work_directory
        seekrflow.make_work_directory(work_dir)
        curdir = os.getcwd()
        os.chdir(seekrflow.work_directory)
        src_pdb_filename = assign_seekrflow_defaults(
            seekrflow, pdb_with_system)
        os.chdir(curdir)
        
    if pdb_with_system != "":
        assign_seekrflow_pdb_with_system(
            seekrflow, pdb_with_system)
        
    flow(seekrflow, instruction, src_pdb_filename, transfer_before, 
         transfer_from_host_only, force_rerun, benchmark_mode, bd_resource_name,
         hidr_resource_name, seekr_resource_name)
    
if __name__ == "__main__":
    main()
    os._exit(0)  # Force exit to prevent Parsl threads from blocking
    # TODO: remove above line later when parsl is totally removed - it
    #  might no longer be necessary.