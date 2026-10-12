"""
client.py

Orchestrate seekrflow single or batch runs by launching workflows,
monitoring progress and statuses, perform transfers, produce output
for reporting, and take input from semaphores or other commands to
direct the workflows.
"""

import os
import asyncio
import argparse
from typing import List, Dict

import seekr.modules.structures as seekr_structures

import seekrflow.modules.structures as structures
import seekrflow.modules.client.structures as client_structures
import seekrflow.modules.client.start as client_start
import seekrflow.modules.client.validation as client_validation
import seekrflow.modules.client.run as client_run
import seekrflow.modules.client.status_and_cancel as client_status_and_cancel

def client(
        batch: client_structures.Batch,
        instruction: str,
        output_file: str,
        control_file: str,
        transfer_targets: List[tuple[str, str]] | None,
        force_targets: List[tuple[str, str]] | None,
        benchmark_target: tuple[str, str] | None,
        semaphore_dict: Dict[tuple[str, str], str] | None,
        ) -> None:
    """
    Orchestrate a batch of seekrflow workflows.
    """
    curdir = os.getcwd()
    output_file = os.path.abspath(output_file)
    control_file = os.path.abspath(control_file)
    if batch.batch_directory is not None:
        os.chdir(batch.batch_directory)
    session = batch.create_session(output_file, control_file)
    for systemrun in session.systemrun_objects:
        seekrflow, model = systemrun.seekrflow, systemrun.model
        client_validation.validate_run_settings(seekrflow, model)
        this_system_semaphore_dict = {}
        if semaphore_dict is not None:
            this_system_semaphore_dict = client_start.handle_semaphore(
                model, semaphore_dict, seekrflow.name)
            systemrun.semaphore_dict = this_system_semaphore_dict
        if benchmark_target is not None:
            benchmark_stage = client_start.handle_benchmark_stage(
                model, benchmark_target, this_system_semaphore_dict,seekrflow.name)
            systemrun.benchmark_stage = benchmark_stage
        if force_targets is not None:
            force_rerun_stages_this_system = client_start.handle_force_rerun(
                model, force_targets, this_system_semaphore_dict, seekrflow.name)
            systemrun.force_rerun_stages = force_rerun_stages_this_system

        all_stage_names = [s.name for s in model.stages]
        if transfer_targets is not None:
            transfer_stages_this_system= client_start.handle_transfer(
                model, transfer_targets, seekrflow.name)
            transferred_resources = client_start.transfer_unique_stage_resources(
                seekrflow, transfer_stages_this_system, backwards=True)
            continue
        
    if transfer_targets is not None:
        return

    if instruction == "run":
        asyncio.run(client_run.launch_session_pipelines(session))
        perform_final_transfer = True
    elif instruction == "status":
        client_status_and_cancel.report_detached_status(session)
    elif instruction == "stop":
        client_status_and_cancel.stop_detached_jobs(session)
    else:
        raise ValueError(f"Unknown instruction: {instruction!r}")
    
    if perform_final_transfer:
        for systemrun in session.systemrun_objects:
            seekrflow, model = systemrun.seekrflow, systemrun.model
            all_stage_names = [s.name for s in model.stages]
            transfer_stages_this_system = all_stage_names
            transferred_resources = client_start.transfer_unique_stage_resources(
                seekrflow, transfer_stages_this_system, backwards=True)
            continue
    
    os.chdir(curdir)
    return

def parse_semaphore_token(token: str) -> tuple[str, str, str]:
    parts = token.split(":")
    if len(parts) == 2:
        stage, value = parts
        return "*", stage, value
    if len(parts) == 3:
        system, stage, value = parts
        return system, stage, value
    raise ValueError(f"Invalid --semaphore token: {token!r}")

def parse_force_benchmark_transfer_token(token: str) -> tuple[str, str]:
    parts = token.split(":")
    if len(parts) == 1:
        return "*", parts[0]          # stage only
    if len(parts) == 2:
        return parts[0], parts[1]     # system, stage
    raise ValueError(f"Invalid --force_rerun or --benchmark token: {token!r}")

def main():
    parser = argparse.ArgumentParser(
        description="Run a single or batch of seekrflow workflows")
    parser.add_argument(
        "instruction", metavar="INSTRUCTION", type=str,
        help="The instruction for client. Options include 'run', which will "
        "run the workflow normally, 'status' will produce a quick status report "
        "and then exit, 'stop' will stop any detached jobs and exit.")
    parser.add_argument(
        "-i", "--input-json", dest="input_json", metavar="INPUT_JSON",
        help="Input JSON file containing the job specification")
    parser.add_argument(
        "-m", "--model_filename", dest="model_filename",
        metavar="MODEL_FILENAME", type=str, default="",
        help="Path to the model file to use for the simulation. This activates "
        "the so-called 'hotshot mode', where an existing model.xml file can be "
        "run easily without needing to prepare a full seekrflow structure. Note "
        "that if this is provided, the seekrflow JSON file must not be provided.")
    parser.add_argument(
        "-n", "--name", dest="name",
        metavar="NAME", type=str, default="",
        help="Name for the simulation or calculation. This is particularly "
        "useful in hotshot mode to give a more informative name to the run - "
        "otherwise, it will be named 'hotshot_<inode>_<device>'.")
    parser.add_argument(
        "-o", "--output-file", dest="output_file", metavar="OUTPUT_FILE",
        type=str, default="", required=True,
        help="Output JSON file containing internal states of jobs")
    parser.add_argument(
        "-c", "--control-file", dest="control_file", metavar="CONTROL_FILE",
        type=str, default="", required=True,
        help="Control JSON file that can change semaphore states, force transfers "\
            "detach jobs, and other control actions.")
    parser.add_argument(
        "-T", "--transfer-from-remote-only", dest="transfer_from_remote_only", 
        metavar="SPEC", type=str, nargs="*", default=None,
        help="Pull files from remote for the given system/stage targets. Each SPEC is "
        "'SYSTEM:STAGE' (use '*' for all systems or all stages). "
        "Examples: -T lig1:mmvt  |  -T *:mmvt  |  -T lig1:*  |  -T lig1:mmvt lig2:mmvt lig3:*"
        "This argument will cause nothing to happen after the transfers - the "
        "client will simply exit.")
    parser.add_argument(
        "-f", "--force_rerun", dest="force_rerun", metavar="SPEC",
        type=str, nargs="*", default=None,
        help="Force-rerun the given system/stage targets. Each SPEC is "
        "'STAGE' (that stage on every system) or 'SYSTEM:STAGE'. "
        "Use '*' for all systems or all stages (e.g. lig1:* or *:mmvt). "
        "Pass -f with no SPEC to force-rerun everything. "
        "Examples: -f mmvt  |  -f lig1:bd lig2:mmvt  |  -f")
    parser.add_argument(
        "-b", "--benchmark", dest="benchmark",
        metavar="SPEC", type=str, default=None,
        help="Run one stage in benchmark mode (short timing run). "
        "SPEC is 'STAGE' or 'SYSTEM:STAGE'. Only one benchmark per "
        "invocation. Example: -b mmvt  |  -b lig1:mmvt")
    parser.add_argument(
        "--semaphore", dest="semaphore", metavar="SPEC",
        type=str, nargs="*", default=None,
        help="Set go/wait/stop for system/stage pairs. Each SPEC is "
        "'STAGE:VALUE' (every system) or 'SYSTEM:STAGE:VALUE'. "
        "VALUE is go (default), wait (no new submits), or stop "
        "(no new submits and kill running work). '*' means all "
        "systems or all stages. "
        "Examples: --semaphore mmvt:stop  |  "
        "--semaphore lig1:mmvt:stop lig2:*:wait")
    #parser.add_argument(
    #    "-O", "--output-socket", dest="output_socket", metavar="OUTPUT_SOCKET",
    #    help="Output socket number for internal state updates"
    #)
    args = vars(parser.parse_args())
    instruction = args.get("instruction")
    input_json = args.get("input_json")
    model_filename = args.get("model_filename")
    name = args.get("name")
    output_file = args.get("output_file")
    control_file = args.get("control_file")
    transfer_from_remote_only = args.get("transfer_from_remote_only")
    force_rerun = args.get("force_rerun")
    benchmark = args.get("benchmark")

    # TODO: implement 'status' instruction to produce quick status reports for
    #  display before launch

    semaphore_dict = {}
    if args["semaphore"]:
        for item in args["semaphore"]:
            system, stage, value = parse_semaphore_token(item)
            # TODO: move this check into the client code
            #if value not in ["go", "wait", "stop"]:
            #    raise ValueError(
            #        f"Invalid semaphore value: {value}. Must be go, wait, or stop")
            semaphore_dict[(system, stage)] = value
    
    benchmark_target = None
    if benchmark is not None:
        benchmark_target = parse_force_benchmark_transfer_token(benchmark)

    force_targets = None
    if force_rerun is not None:
        if not force_rerun:
            force_targets = [("*", "*")]
        else:
            force_targets = [parse_force_benchmark_transfer_token(t) \
                for t in force_rerun]

    transfer_targets = None
    if transfer_from_remote_only is not None:
        if not transfer_from_remote_only:
            transfer_targets = [("*", "*")]
        else:
            transfer_targets = [parse_force_benchmark_transfer_token(t) \
                for t in transfer_from_remote_only]

    hotshot_mode = False
    if input_json == "":
        seekrflow = structures.Seekrflow()
        assert model_filename != "", \
            "Model filename must be provided in hotshot mode."
        hotshot_mode = True
        assert os.path.exists(model_filename), \
            f"Model file {model_filename} does not exist."
        if name != "":
            seekrflow.name = name
        else:
            model_dirname = os.path.dirname(model_filename)
            if model_dirname == "":
                model_dirname = os.path.abspath(os.curdir)
            st = os.stat(model_dirname)
            seekrflow.name = "hotshot_" + str(st.st_dev) + "_" + str(st.st_ino)
        print("Running in hotshot mode with seekrflow name:", seekrflow.name)
        seekrflow.work_directory = None
        seekrflow.root_directory = os.path.dirname(
            os.path.abspath(model_filename))
        structures.try_to_load_resources_json(seekrflow)
        batch = client_start.make_batch_from_single_seekrflow(seekrflow)

    else:
        assert model_filename == "", \
            "If input JSON file is provided, model filename must not be provided."
        assert os.path.exists(input_json), \
            f"Input JSON file {input_json} does not exist."
        
        # TODO: load either the seekrflow JSON or the batch JSON

        if client_start.input_is_batch_file(input_json):
            #session = client_start.make_session_from_batch_file(input_json)
            batch = client_structures.Batch.load_batch_file(input_json)
        else:
            seekrflow = structures.load_seekrflow(input_json)
            if name != "":
                seekrflow.name = name
            work_dir = os.path.abspath(
                os.path.expanduser(seekrflow.work_directory))
            seekrflow.make_work_directory(work_dir)
            batch = client_start.make_batch_from_single_seekrflow(seekrflow)
    
    client(batch, instruction, output_file, control_file, transfer_targets, 
           force_targets, benchmark_target, semaphore_dict)

if __name__ == "__main__":
    main()