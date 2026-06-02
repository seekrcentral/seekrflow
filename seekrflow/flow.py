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
        pdb_with_system: str
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
        if pdb_with_system != "":
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
        transfer_from_remote_only: str | None = None,
        force_rerun: list[str] | None = None,
        benchmark_stage: str | None = None,
        semaphore_dict: dict[str, str] | None = None,
        resource_dict: dict[str, str] | None = None,
        keystrokes_enabled: bool = True,
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
            seekrflow, transfer_before, transfer_from_remote_only, force_rerun,
            benchmark_stage=benchmark_stage, stage_resource_dict=resource_dict, 
            semaphore_dict=semaphore_dict,
            keystrokes_enabled=keystrokes_enabled)
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
            seekrflow, transfer_before, transfer_from_remote_only, force_rerun,
            benchmark_stage=benchmark_stage, stage_resource_dict=resource_dict, 
            semaphore_dict=semaphore_dict,
            keystrokes_enabled=keystrokes_enabled)
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
        "if this is provided, the seekrflow JSON file must not be provided.")
    argparser.add_argument(
        "-p", "--pdb_with_system", dest="pdb_with_system", 
        metavar="PDB_WITH_SYSTEM", type=str, default="",
        help="Path to the input PDB file that contains the solvated molecules.")
    argparser.add_argument(
        "-n", "--name", dest="name",
        metavar="NAME", type=str, default="",
        help="Name for the simulation or calculation. This is particularly useful "\
            "in hotshot mode to give a more informative name to the run - otherwise, "\
            "it will be named 'hotshot_<inode>_<device>'.")
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
        "which will not transfer files.")
    argparser.add_argument(
        "-T", "--transfer_from_remote_only", dest="transfer_from_remote_only", 
        metavar="STAGE", type=str, default=None,
        help="If set, transfer files from the remote host for a particular stage "\
        "without doing any preparation or running. Default: None, "\
        "which will not transfer files.")
    argparser.add_argument(
        "-f", "--force_rerun", dest="force_rerun", metavar="STAGE",
        type=str, nargs="*", default=None,
        help="If set, force re-run of stages. Provide a space-separated list "
        "of stage names (e.g. -f bd mmvt) to force-rerun only those stages, "
        "or pass -f with no arguments to force-rerun all stages. Any running "
        "processes for the affected stages will be killed and the stages "
        "restarted with force_overwrite=True.")
    argparser.add_argument(
        "-b", "--benchmark", dest="benchmark",
        metavar="STAGE", type=str, default=None,
        help="If set, run the named stage in benchmark mode - a quick run on "\
        "the resource to get approximate timings. At most, a single "\
        "stage can be benchmarked per invocation, since dependent stages "\
        "would not produce the outputs needed downstream. Default: None.")
    argparser.add_argument(
        "--no-keystrokes", dest="no_keystrokes", action="store_true",
        default=False,
        help="Disable the interactive keystroke command watcher. Useful for "\
        "non-interactive runs, GUI integration, or when stdin is not a TTY. "\
        "Default: keystrokes enabled when stdin is a TTY.")
    argparser.add_argument(
        "--semaphore", dest="semaphore", 
        metavar="STAGE_CONTROL", type=str, default=None,
        help="Control stage execution. Format: 'stage:value,stage:value,...' " \
        "where stage is any configured stage name and value is go/wait/stop. " \
        "go: normal operation (default), wait: don't submit new jobs but let "\
        "running jobs finish, stop: don't submit new jobs AND kill any running "\
        "jobs. Example: --semaphore stage1:stop,stage2:wait. --force_overwrite "\
        "and --benchmark not allowed for 'stop' semaphore in a given stage, "\
        "if one desires to launch a stage with --force_overwrite or --benchmark "\
        "at a later time, set semaphore to 'wait' for that stage.")
    argparser.add_argument(
        "-r", "--resources", dest="resources", 
        metavar="STAGE_CONTROL", type=str, default=None,
        help="Control stage resource assignments. Format: 'stage:value,"\
        "stage:value,...' where stage is any configured stage name and value "\
        "is the name of the resource as defined in the seekrflow input JSON file "\
        "or in the user's .seekrflow_resources.json file. Example: --resources "\
        "stage1:local,stage2:supercomputerName.")
    
    args = argparser.parse_args()
    args = vars(args)
    instruction = args["instruction"]
    input_json = args["input_json"]
    model_filename = args["model_filename"]
    pdb_with_system = args["pdb_with_system"]
    name = args["name"]
    work_dir = args["work_directory"]
    transfer_before = args["transfer_before"]
    transfer_from_remote_only = args["transfer_from_remote_only"]
    force_rerun = args["force_rerun"]
    benchmark = args["benchmark"]
    keystrokes_enabled = not args["no_keystrokes"]


    if pdb_with_system != "":
        #pdb_with_system = os.path.abspath(pdb_with_system)
        assert os.path.exists(pdb_with_system), \
            f"PDB file with system {pdb_with_system} does not exist."
        
    # Parse semaphore argument using generic stage names.
    semaphore_dict = {}
    if args["semaphore"]:
        for item in args["semaphore"].split(","):
            parts = item.strip().split(":")
            if len(parts) != 2:
                raise ValueError(
                    f"Invalid semaphore format: {item}. Expected 'stage:value'")
            stage, value = parts[0].strip(), parts[1].strip()
            if value not in ["go", "wait", "stop"]:
                raise ValueError(
                    f"Invalid semaphore value: {value}. Must be go, wait, or stop")
            semaphore_dict[stage] = value
    
    # Check for conflicting force_rerun and semaphore="stop"
    if force_rerun is not None:
        stages_to_check = force_rerun if force_rerun \
            else list(semaphore_dict.keys())
        for stage in stages_to_check:
            if semaphore_dict.get(stage) == "stop":
                raise ValueError(
                    f"Conflicting options: --force-rerun {stage} and "
                    f"--semaphore {stage}:stop. "
                    "Cannot force rerun a stopped stage.")
    
    # Check for conflicting benchmark and semaphore="stop"
    if benchmark is not None:
        stages_to_check = benchmark if benchmark \
            else list(semaphore_dict.keys())
        for stage in stages_to_check:
            if semaphore_dict.get(stage) == "stop":
                raise ValueError(
                    f"Conflicting options: --benchmark {stage} and "
                    f"--semaphore {stage}:stop. "
                    "Cannot benchmark a stopped stage.")
            
    # Parse resource assignment argument using generic stage names.
    resource_dict = {}
    if args["resources"]:
        for item in args["resources"].split(","):
            parts = item.strip().split(":")
            if len(parts) != 2:
                raise ValueError(
                    f"Invalid resources format: {item}. Expected 'stage:value'")
            stage, value = parts[0].strip(), parts[1].strip()
            resource_dict[stage] = value

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
            # Generate a unique name based on the directory's inode and device
            model_dirname = os.path.dirname(model_filename)
            if model_dirname == "":
                model_dirname = os.path.abspath(os.curdir)
            st = os.stat(model_dirname)
            if name != "":
                seekrflow.name = name
            else:
                seekrflow.name = "hotshot_" + str(st.st_dev) + "_" + str(st.st_ino)
            print("Running in hotshot mode with seekrflow name:", seekrflow.name)
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
        if name != "":
            seekrflow.name = name
        
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
         transfer_from_remote_only, force_rerun, benchmark, 
         semaphore_dict=semaphore_dict,
         resource_dict=resource_dict,
         keystrokes_enabled=keystrokes_enabled)
    
if __name__ == "__main__":
    main()