"""
modules/seekr_run.py

Using the seekr API, run stages of the seekr calculation.
"""

import os
import sys
import time
import json
import glob
import select
import signal
import pathlib
import multiprocessing

import seekr2.modules.common_base as seekr2_base

import seekrflow.modules.base as base
import seekrflow.modules.structures as structures
import seekrflow.modules.transfer.base as transfer_base
import seekrflow.modules.workload_managers.local_multiprocessing as workload_local_mp
import seekrflow.modules.workload_managers.remote as workload_remote

KEYSTROKE_CHECK_INTERVAL = 1.0 # seconds
MAIN_LOOP_INTERVAL = 30.0 # 600.0  # seconds

def determine_stage_manager_status(manager_status: dict | None) -> str:
    """
    Determine the manager status string based on job states.
    
    Returns "idle" if no jobs, or "running", "queued", or "running/queued"
    based on the State field of jobs in the manager_status["jobs"] list.
    """
    if manager_status is None:
        return "idle"
    
    jobs = manager_status.get("jobs", [])
    if len(jobs) == 0:
        return "idle"
    
    # Collect all unique job states (convert to lowercase)
    states = set(job.get("State", "UNKNOWN").lower() for job in jobs)
    
    # Determine combined status
    if "running" in states and (("queued" in states) or ("pending" in states)):
        return "running/queued"
    elif "running" in states:
        return "running"
    elif ("queued" in states) or ("pending" in states):
        return "queued"
    else:
        return "unknown"

def get_stage_pids(stage_name: str, stage_processes: dict, 
                   root_directory_path: pathlib.Path) -> list:
    """
    Get all PIDs for a stage from either process objects or state files.
    
    Returns:
        For BD/HIDR: [(pid, None)] or []
        For SEEKR: [(pid, anchor_index), ...] or []
    """
    pids = []
    process_or_dict = stage_processes.get(stage_name)
    
    if stage_name == "seekr" and isinstance(process_or_dict, dict):
        # Check direct process objects
        for anchor_idx, process in process_or_dict.items():
            if process is not None and process.is_alive():
                pids.append((process.pid, anchor_idx))
        
        # Check state files if no direct processes (reattached case)
        if len(pids) == 0:
            state_dir = workload_local_mp.get_local_state_dir(root_directory_path)
            for state_file in state_dir.glob("seekr_anchor_*_state_*.json"):
                try:
                    state = workload_local_mp.LocalProcessState.load(state_file)
                    if state.status == "running" and workload_local_mp.check_pid_exists(state.pid):
                        anchor_idx = state.progress_info.get("anchor_index") if state.progress_info else None
                        if anchor_idx is not None:
                            pids.append((state.pid, anchor_idx))
                except Exception as e:
                    print(f"Warning: Could not load state file: {e}")
    else:
        # Single process (BD/HIDR)
        if process_or_dict is not None and process_or_dict.is_alive():
            pids.append((process_or_dict.pid, None))
        else:
            # Check state file for reattached process
            state_file = workload_local_mp.find_latest_state_file(root_directory_path, stage_name)
            if state_file and state_file.exists():
                try:
                    state = workload_local_mp.LocalProcessState.load(state_file)
                    if state.status == "running" and workload_local_mp.check_pid_exists(state.pid):
                        pids.append((state.pid, None))
                except Exception as e:
                    print(f"Warning: Could not load state file: {e}")
    
    return pids

def kill_process_group(pid: int, stage_name: str, root_directory_path: pathlib.Path, 
                       anchor_index: int = None, timeout: int = 5) -> None:
    """
    Kill a process group with escalation from SIGTERM to SIGKILL.
    Updates state file if force-killed.
    """
    stage_label = f"{stage_name} anchor {anchor_index}" if anchor_index is not None else stage_name
    print(f"Terminating {stage_label} process group PID {pid}...")
    
    try:
        os.killpg(pid, signal.SIGTERM)
    except (ProcessLookupError, PermissionError) as e:
        print(f"Warning: {e}")
        return
    
    time.sleep(timeout)
    
    # Check if still alive and force kill
    if workload_local_mp.check_pid_exists(pid):
        print(f"Force killing {stage_label} process group PID {pid}...")
        try:
            os.killpg(pid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError) as e:
            print(f"Warning: {e}")
        
        time.sleep(1)
        
        # Update state file if successfully killed
        if not workload_local_mp.check_pid_exists(pid):
            try:
                if anchor_index is not None:
                    state_dir = workload_local_mp.get_local_state_dir(root_directory_path)
                    state_file = state_dir / f"seekr_anchor_{anchor_index}_state_{pid}.json"
                else:
                    state_file = workload_local_mp.find_latest_state_file(
                        root_directory_path, stage_name)
                
                if state_file and state_file.exists():
                    state = workload_local_mp.LocalProcessState.load(state_file)
                    if state.status == "running":
                        state.status = "killed"
                        state.ended_at = time.time()
                        state.error = "Process force-killed by user command"
                        state.save(state_file)
            except Exception as e:
                print(f"Warning: Could not update state file: {e}")

def check_and_raise_if_process_failed(stage_name: str, stage_processes: dict, 
                                      root_directory_path: pathlib.Path) -> None:
    """
    Check if a local process has failed by examining state files.
    Raises an exception if failure detected.
    """
    if stage_name == "seekr":
        # Check all anchor processes
        seekr_processes = stage_processes.get("seekr", {})
        if isinstance(seekr_processes, dict):
            for anchor_index, process in seekr_processes.items():
                if process is not None and not process.is_alive():
                    state_pattern = os.path.join(
                        root_directory_path, f"seekr_anchor_{anchor_index}_state_*.json")
                    state_files = glob.glob(state_pattern)
                    for state_file in state_files:
                        try:
                            with open(state_file, 'r') as f:
                                state = json.load(f)
                            if state.get("status") == "failed":
                                error_msg = f"SEEKR anchor {anchor_index} failed. Error: {state.get('error')}"
                                if state.get("traceback"):
                                    error_msg += f"\nTraceback: {state.get('traceback')}"
                                raise Exception(error_msg)
                        except Exception as e:
                            if "failed" in str(e):
                                raise
    else:
        # Check single process (BD/HIDR)
        process = stage_processes.get(stage_name)
        if process is not None and not process.is_alive():
            # Check exit code
            if hasattr(process, 'exitcode') and process.exitcode is not None and process.exitcode != 0:
                state_file = workload_local_mp.find_latest_state_file(root_directory_path, stage_name)
                if state_file:
                    try:
                        state = workload_local_mp.LocalProcessState.load(state_file)
                        error_msg = f"Local {stage_name.upper()} stage failed. Error: {state.error}"
                        if state.traceback:
                            error_msg += f"\nTraceback:\n{state.traceback}"
                        raise Exception(error_msg)
                    except json.JSONDecodeError:
                        raise Exception(f"Local {stage_name.upper()} stage failed with exit code {process.exitcode}")
                else:
                    raise Exception(f"Local {stage_name.upper()} stage failed with exit code {process.exitcode}")

def run_model(
        seekrflow: structures.Seekrflow,
        transfer_before: str | None = None,
        transfer_from_host_only: str | None = None,
        force_rerun: list[str] | None = None,
        benchmark_mode: bool = False,
        bd_resource_name: str | None = None,
        hidr_resource_name: str | None = None,
        seekr_resource_name: str | None = None
        ) -> None:
    """
    Run the SEEKR calculation using remote and local resources.
    """
    curdir = os.getcwd()
    if seekrflow.work_directory is not None:
        seekrflow.work_directory = os.path.abspath(seekrflow.work_directory)
        os.chdir(seekrflow.work_directory)

    bd_n_threads = 1
    if seekrflow.workflow.bd_settings is not None:
        bd_n_threads = seekrflow.workflow.bd_settings.num_threads
    if seekrflow.root_directory is None:
        root_directory = os.path.join(seekrflow.work_directory, structures.ROOT)
    else:
        root_directory = seekrflow.root_directory
    root_directory_path = pathlib.Path(root_directory)
    model_filename = os.path.join(root_directory, "model.xml")
    job_status_filename = os.path.join(root_directory, ".seekrflow_job_status.json")
    model = seekr2_base.load_model(model_filename)
    using_bd = model.using_bd()
    if benchmark_mode:
        print("Benchmark mode enabled: only short MD runs will be performed.")
        using_bd = False
    # Run through the stages, see which are completed, started, running, queued, 
    # unstarted, and skipping.
    stage_names = ["bd", "hidr", "seekr"]
    stage_dependencies = {
        "bd": [],
        "hidr": [],
        "seekr": ["hidr"]
    }
    stage_locations = {}
    stage_progress = {}
    stage_progress_bar = {}
    stage_should_perform_transfer = {}
    stage_manager_status = {}
    stage_run_results = {}
    stage_job_ids = {}
    # seekr below is dict of {anchor_index: process}
    stage_processes = {"bd": None, "hidr": None, "seekr": {}}  
    stage_statuses = {}
    stage_running_locally = {}
    
    # Override resource names if provided
    if bd_resource_name is not None:
        seekrflow.run_settings.bd_stage_resource_name = bd_resource_name
    if hidr_resource_name is not None:
        seekrflow.run_settings.hidr_stage_resource_name = hidr_resource_name
    if seekr_resource_name is not None:
        seekrflow.run_settings.seekr_stage_resource_name = seekr_resource_name

    # Check which stages run locally
    stage_running_locally["bd"] \
        = using_bd and seekrflow.run_settings.bd_stage_resource_name == "local"
    stage_running_locally["hidr"] \
        = seekrflow.run_settings.hidr_stage_resource_name == "local"
    stage_running_locally["seekr"] \
        = seekrflow.run_settings.seekr_stage_resource_name == "local"

    # Check for existing processes from previous session (reattachment after detach)
    print("Checking for existing processes from previous session...")
    for stage_name in ["bd", "hidr", "seekr"]:
        existing_state = workload_local_mp.check_for_existing_local_processes(
            root_directory_path, stage_name)
        if existing_state:
            # Note: We can't truly "reattach" to a multiprocessing.Process object,
            # but we can track the PID and monitor/kill it via the state file
            print(f"  Reattached to {stage_name} process (PID: {existing_state.pid})")
            # stage_processes[stage_name] remains None - we'll use the PID from state file

    # Handle force_rerun - kill any existing processes for those stages, then delete files
    if force_rerun is not None:
        for stage_name in force_rerun:
            if stage_name in stage_names:
                print(f"Forcing re-run of stage: {stage_name}")
                
                if stage_running_locally[stage_name]:
                    # Local force rerun
                    workload_local_mp.force_rerun_stage(
                        model, stage_name, root_directory_path)
                else:
                    # Remote force rerun
                    workload_remote.force_rerun_stage_remote(
                        seekrflow, stage_name, silent=False)
                    # Reload model after HIDR force rerun
                    if stage_name == "hidr":
                        model = seekr2_base.load_model(model_filename)
            else:
                raise ValueError(f"Stage name {stage_name} not recognized.")
    
    # Setup signal handler for graceful shutdown
    def signal_handler_local(sig, frame, curdir):
        print("\n\n=== Received interrupt signal ===")
        workload_local_mp.cleanup_local_processes(stage_processes, root_directory_path)
        print("Exiting seekrflow.")
        os.chdir(curdir)
        sys.exit(0)

    signal.signal(signal.SIGINT, signal_handler_local)
    signal.signal(signal.SIGTERM, signal_handler_local)

    iteration = 0
    silent = False
    keep_running = True
    perform_final_transfer = True
    do_process_cleanup = True
    instruction_str = f"""Options: 
    i [STAGE] - print (i)nformation for STAGE (bd, hidr, seekr, all),
    d - (d)etach to let jobs keep running and exit seekrflow,
    k - (k)ill all running/queued jobs,
    t [STAGE] - (t)ransfer remote files back for STAGE (bd, hidr, seekr, all).
    """
    print(instruction_str, flush=True)
    
    try:
        while keep_running:
            # TODO: check each stage based on whether they were submitted or queued in the
            #  previous step. If there's a sudden change (like a job finishing), there might
            #  be an error in configuration or something.
            print("ITERATION:", iteration)
            end_of_loop_time = time.time() + MAIN_LOOP_INTERVAL
            # Get stage and manager statuses
            for stage_name in stage_names:
                stage_progress[stage_name] = "unknown"
                stage_should_perform_transfer[stage_name] = False
                stage_job_ids[stage_name] = []
                if stage_name == "bd":
                    stage_locations[stage_name] = seekrflow.run_settings.bd_stage_resource_name
                    if not using_bd:
                        bd_status = None
                        stage_progress[stage_name] = "omitted"
                        stage_manager_status[stage_name] = "n/a"
                        continue
                    if seekrflow.run_settings.bd_stage_resource_name == "local":
                        bd_status = workload_local_mp.bd_status_local(
                            model, stage_processes, root_directory_path,
                            n_trajectories_for_completion=model.k_on_info\
                                .b_surface_num_trajectories)
                        # Check if process failed
                        check_and_raise_if_process_failed(stage_name, stage_processes, root_directory_path)
                    else:
                        bd_status = workload_local_mp.bd_status_remote(
                            seekrflow, model, silent=silent)
                        if not bd_status["stage_status"]["model_xml_found"]:
                            stage_should_perform_transfer[stage_name] = True
                    
                    silent = True
                    if bd_status is not None:
                        stage_statuses[stage_name] = bd_status
                        stage_progress[stage_name] = bd_status["stage_status"]["progress"]
                        stage_progress_bar[stage_name] = bd_status["stage_status"]["progress_bar"]
                        manager_status = bd_status["manager_status"]
                        stage_manager_status[stage_name] = determine_stage_manager_status(manager_status)
                        if manager_status and manager_status.get("jobs"):
                            for job in manager_status["jobs"]:
                                stage_job_ids[stage_name].append(job["JobID"])
            
                elif stage_name == "hidr":
                    stage_locations[stage_name] = seekrflow.run_settings.hidr_stage_resource_name
                    stage_job_ids[stage_name] = []
                    if seekrflow.run_settings.hidr_stage_resource_name == "local":
                        model = seekr2_base.load_model(model_filename)
                        hidr_status = workload_local_mp.hidr_status_local(
                            model, stage_processes, root_directory_path)
                        # Check if process failed
                        check_and_raise_if_process_failed(stage_name, stage_processes, root_directory_path)
                    else:
                        hidr_status = workload_remote.hidr_status_remote(
                            seekrflow, model, silent=silent)
                        stage_progress[stage_name] = hidr_status["stage_status"]["progress"]
                        if not hidr_status["stage_status"]["model_xml_found"] \
                                and (not stage_progress[stage_name] in ["completed", "unknown"]):
                            stage_should_perform_transfer[stage_name] = True

                    silent = True
                    stage_statuses[stage_name] = hidr_status
                    stage_progress[stage_name] = hidr_status["stage_status"]["progress"]
                    stage_progress_bar[stage_name] = hidr_status["stage_status"]["progress_bar"]
                    manager_status = hidr_status["manager_status"]
                    stage_manager_status[stage_name] = determine_stage_manager_status(manager_status)
                    if manager_status and manager_status.get("jobs"):
                        for job in manager_status["jobs"]:
                            stage_job_ids[stage_name].append(job["JobID"])

                elif stage_name == "seekr":
                    stage_locations[stage_name] = seekrflow.run_settings.seekr_stage_resource_name
                    stage_job_ids[stage_name] = []
                    if seekrflow.run_settings.seekr_stage_resource_name == "local":
                        seekr_status = workload_local_mp.seekr_anchors_status_local(
                            model, stage_processes, root_directory_path, benchmark_mode)
                        # Check if any anchor processes failed
                        check_and_raise_if_process_failed(stage_name, stage_processes, root_directory_path)
                    else:
                        seekr_status = workload_remote.seekr_status_remote(
                            seekrflow, model, silent=silent, benchmark_mode=benchmark_mode)
                        if (not seekr_status["stage_status"]["model_xml_found"]) \
                                and (not stage_progress[stage_name] in ["completed", "unknown"]):
                            stage_should_perform_transfer[stage_name] = True
                        
                    silent = True
                    stage_statuses[stage_name] = seekr_status
                    stage_progress[stage_name] = seekr_status["stage_status"]["progress"]
                    stage_progress_bar[stage_name] = seekr_status["stage_status"]["progress_bar"]
                    manager_status = seekr_status["manager_status"]
                    stage_manager_status[stage_name] = determine_stage_manager_status(manager_status)
                    if manager_status and manager_status.get("jobs"):
                        for job in manager_status["jobs"]:
                            stage_job_ids[stage_name].append(job["JobID"])

            # Write to the status file
            with open(job_status_filename, "w") as f:
                json.dump(stage_statuses, f, indent=4)

            stage_should_submit_new_run = {}
            number_of_stages_to_submit_new_run = 0
            number_of_running_stages = 0
            for stage_name in stage_names:
                if stage_name == "bd" and not using_bd:
                    continue
                stage_should_submit_new_run[stage_name] = False
                if stage_progress[stage_name] in ["unstarted", "started"]:
                    stage_should_submit_new_run[stage_name] = True
                else:
                    stage_should_submit_new_run[stage_name] = False
                # Check dependencies
                for dependency in stage_dependencies[stage_name]:
                    if stage_progress[dependency] != "completed":
                        stage_should_submit_new_run[stage_name] = False
                if stage_manager_status[stage_name] not in ["idle", "n/a"]:
                    stage_should_submit_new_run[stage_name] = False
                    number_of_running_stages += 1
                if stage_should_submit_new_run[stage_name]:
                    number_of_stages_to_submit_new_run += 1

            if number_of_stages_to_submit_new_run == 0 and number_of_running_stages == 0:
                print("No calculations remaining. Exiting seekrflow...")
                keep_running = False
            
            if transfer_from_host_only:
                for stage_name in stage_names:
                    resource_name = stage_locations[stage_name]
                    if transfer_from_host_only in [stage_name, "all"] \
                            and resource_name != "local":
                        resource: structures.Resource_remote_base \
                            = seekrflow.run_settings.get_resource_by_name(resource_name)
                        if stage_progress[stage_name] in ["started", "completed"]:
                            transfer_base.transfer_files_to_from_remote_resource(
                                seekrflow.name, resource, root_directory, backwards=True)
                return

            remotes_transferred_to = set()
            for stage_name in stage_names:
                if stage_should_perform_transfer[stage_name] \
                        or transfer_before in [stage_name, "all"]:
                    resource_name = stage_locations[stage_name]
                    if (resource_name != "local") \
                            and (resource_name not in remotes_transferred_to):
                        resource: structures.Resource_remote_base \
                            = seekrflow.run_settings.get_resource_by_name(resource_name)
                        transfer_base.transfer_files_to_from_remote_resource(
                            seekrflow.name, resource, root_directory)
                        remotes_transferred_to.add(resource_name)

            transfer_before = None  # only transfer once at beginning
            
            if using_bd and stage_should_submit_new_run["bd"]:
                if stage_locations["bd"] != "local":
                    resource: structures.Resource_remote_base \
                        = seekrflow.run_settings.get_resource_by_name(stage_locations["bd"])
                    destination_path = os.path.join(
                        resource.remote_working_directory, seekrflow.name)
                    destination_model_filename = os.path.join(destination_path, "model.xml")
                    output_filename = os.path.join(destination_path, "bd_run.out")
                    command_string = f"python {resource.remote_seekr2_directory}/run.py "\
                        f"any_bd model.xml -N {bd_n_threads} > {output_filename}" 
                    bd_run_result = workload_remote.submit_remote_run_workflow(
                                seekrflow, "bd", destination_path, resource,
                                command_string, destination_model_filename, workflow_type="bd",
                                silent=silent)

                    stage_run_results["bd"] = bd_run_result
                    stage_manager_status["bd"] = "running"
                    stage_progress["bd"] = "started"
                
            if stage_should_submit_new_run["hidr"]:
                if seekrflow.run_settings.hidr_stage_resource_name != "local":
                    resource: structures.Resource_remote_base \
                        = seekrflow.run_settings.get_resource_by_name(stage_locations["hidr"])
                    destination_path = os.path.join(
                        resource.remote_working_directory, seekrflow.name)
                    destination_model_filename = os.path.join(destination_path, "model.xml")
                    input_pdb_filename = os.path.basename(
                        seekrflow.workflow.solvated_system_for_md.solvated_pdb)
                    output_filename = os.path.join(destination_path, "hidr_run.out")
                    if seekrflow.workflow.hidr_settings.type == "hidr_metaD":
                        gaussian_height = seekrflow.workflow.hidr_settings.gaussian_height
                        gaussian_width = seekrflow.workflow.hidr_settings.gaussian_width
                        bias_factor = seekrflow.workflow.hidr_settings.bias_factor
                        command_string \
                            = f"python {resource.remote_seekrtools_directory}/hidr/hidr.py any "\
                            f"model.xml -M metadyn -p {input_pdb_filename} "\
                            f"-H {gaussian_height} -w {gaussian_width} -b {bias_factor} -c 0 "\
                            f" > {output_filename}"

                    elif seekrflow.workflow.hidr_settings.type == "hidr_SMD":
                        restraint_force_constant = seekrflow.workflow.hidr_settings.restraint_force_constant
                        translation_velocity = seekrflow.workflow.hidr_settings.translation_velocity
                        command_string \
                            = f"python {resource.remote_seekrtools_directory}/hidr/hidr.py any"\
                            f" model.xml -M SMD -p {input_pdb_filename} "\
                            f"-k {restraint_force_constant} -v {translation_velocity} -c 0 "\
                            f" > {output_filename}"

                    else:
                        raise NotImplementedError(
                            f"HIDR type {seekrflow.workflow.hidr_settings.type} is not implemented.")
                    hidr_run_result = workload_remote.submit_remote_run_workflow(
                        seekrflow, "hidr", destination_path, resource,
                        command_string, destination_model_filename, workflow_type="hidr",
                        silent=silent)

                    stage_run_results["hidr"] = hidr_run_result
                    stage_manager_status["hidr"] = "running"
                    stage_progress["hidr"] = "started"
                
            if stage_should_submit_new_run["seekr"]:
                # TODO: what if the resources are different between hidr and seekr stages??
                # gotta perform a transfer. Assert they are the same for now...
                assert seekrflow.run_settings.get_resource_by_name(stage_locations["seekr"]) \
                    == seekrflow.run_settings.get_resource_by_name(stage_locations["hidr"]), \
                    "For now, SEEKR and HIDR stages must be placed on the same resources."
                if seekrflow.run_settings.seekr_stage_resource_name != "local":
                    resource: structures.Resource_remote_base \
                        = seekrflow.run_settings.get_resource_by_name(stage_locations["seekr"])
                    destination_path = os.path.join(
                        resource.remote_working_directory, seekrflow.name)
                    destination_model_filename = os.path.join(destination_path, "model.xml")
                    output_filename = os.path.join(destination_path, "seekr_{index}_run.out")
                    if benchmark_mode:
                        benchmark_string = f"-t {base.BENCHMARK_MIN_SIMULATION_LENGTH} "
                    else:
                        benchmark_string = ""
                    command_string = f"python {resource.remote_seekr2_directory}/run.py "\
                        f"{{index}} model.xml -c 0 {benchmark_string}> {output_filename}"
                    indices = seekr_status["stage_status"]["incomplete_anchors"]
                    seekr_run_result = workload_remote.submit_remote_run_workflow(
                        seekrflow, "seekr", destination_path, resource,
                        command_string, destination_model_filename, workflow_type="seekr", 
                        indices=indices)
                    stage_run_results["seekr"] = seekr_run_result
                    stage_manager_status["seekr"] = "running"
                    stage_progress["seekr"] = "started"
                
                # Run BD - if needed
            if using_bd:
                if stage_running_locally["bd"] and stage_should_submit_new_run["bd"]:
                    # Start BD process
                    bd_process = multiprocessing.Process(
                        target=workload_local_mp.run_bd_locally,
                        args=(model_filename, bd_n_threads, root_directory_path)
                    )
                    bd_process.start()
                    stage_processes["bd"] = bd_process
                
                    stage_manager_status["bd"] = "running"
                    stage_progress["bd"] = "started"
                
            # Run HIDR if needed
            if stage_running_locally["hidr"] and stage_should_submit_new_run["hidr"]:
                # Get the input PDB file basename
                input_pdb_file = os.path.basename(
                    seekrflow.workflow.solvated_system_for_md.solvated_pdb)
                # Start HIDR process
                process = multiprocessing.Process(
                    target=workload_local_mp.run_hidr_locally,
                    args=(model_filename, input_pdb_file, root_directory_path),
                    daemon=False
                )
                process.start()
                stage_processes["hidr"] = process
                stage_manager_status["hidr"] = "running"
                stage_progress["hidr"] = "started"
            
            # Run SEEKR if needed and HIDR is finished
            if stage_running_locally["seekr"] and stage_should_submit_new_run["seekr"]:
                # Check which anchors need to be filled
                seekr_status = workload_local_mp.seekr_anchors_status_local(
                    model, stage_processes, root_directory_path, benchmark_mode)
                for unfilled_anchor_index in seekr_status["stage_status"]["incomplete_anchors"]:
                    # Start process for this anchor
                    process = multiprocessing.Process(
                        target=workload_local_mp.run_seekr_locally,
                        args=(model_filename, unfilled_anchor_index, root_directory_path, benchmark_mode),
                        daemon=False
                    )
                    process.start()
                    stage_processes["seekr"][unfilled_anchor_index] = process
                    break
            
                # NOTE: this probably must be done for local job where there's just one GPU
                # but for the cluster, we will want to spin them off all at once.
                stage_manager_status["seekr"] = "running"
                stage_progress["seekr"] = "started"
        
            # Now pause to accept user keystrokes
            while True:
                if not keep_running:
                    break
                rlist, _, _ = select.select([sys.stdin], [], [], KEYSTROKE_CHECK_INTERVAL)
                break_out_of_keystroke_loop = False
                if rlist:
                    user_input = sys.stdin.readline().strip().split()
                    command = user_input[0]
                    arg = user_input[1] if len(user_input) > 1 else None
                    if command == "i":
                        printed_nothing = True
                        for stage_name in stage_names:
                            if stage_name == "bd" and not using_bd:
                                continue
                            if arg is None or arg == stage_name or arg == "all":
                                print(f"Stage: {stage_name}")
                                print(f"  Location: {stage_locations[stage_name]}")
                                print(f"  Progress: {stage_progress[stage_name]}")
                                print(f"  Progress bar: {100.0*stage_progress_bar[stage_name]:.1f}%")
                                print(f"  Manager status: {stage_manager_status[stage_name]}")
                                print(f"  Should perform transfer: {stage_should_perform_transfer[stage_name]}")
                                printed_nothing = False
                    
                        if printed_nothing:
                            print(f"No stage found for argument: {arg}.")

                        print(instruction_str, flush=True)
                    
                    elif command == "d":
                        print(f"Detaching from job...")
                        # Allow detaching even with local processes running
                        # The processes will continue running and can be reattached later
                        keep_running = False
                        break_out_of_keystroke_loop = True
                        perform_final_transfer = False
                        
                        # Inform user about running local processes
                        running_local_stages = []
                        for stage_name in stage_names:
                            if stage_locations[stage_name] == "local" and stage_manager_status[stage_name] != "idle":
                                running_local_stages.append(stage_name)
                        
                        if running_local_stages:
                            print(f"Local processes for {', '.join(running_local_stages)} will continue running.")
                            print(f"Run this program again to reattach and monitor progress.")
                            do_process_cleanup = False
                        
                        break

                    elif command == "k":
                        for stage_name in stage_names:
                            if stage_manager_status[stage_name] != "idle":
                                print(f"Killing stage {stage_name}...")
                                if stage_locations[stage_name] == "local":
                                    # Get all PIDs (handles both direct and reattached processes)
                                    pids_to_kill = get_stage_pids(stage_name, stage_processes, root_directory_path)
                                    
                                    # Kill all found processes
                                    for pid, anchor_idx in pids_to_kill:
                                        kill_process_group(pid, stage_name, root_directory_path, anchor_idx)
                                    
                                    # Clear manager status
                                    if stage_statuses.get(stage_name):
                                        stage_statuses[stage_name]["manager_status"]["error"] = "killed"
                                        stage_statuses[stage_name]["manager_status"]["jobs"] = []
                                else:
                                    # Remote job killing
                                    job_ids = stage_job_ids.get(stage_name)
                                    for job_id in job_ids:
                                        print(f"killing job_id: {job_id}")
                                        if job_id is not None:
                                            workload_remote.submit_remote_cancel_workflow(
                                                seekrflow, stage_locations[stage_name], job_id, 
                                                silent=silent)
                        keep_running = False
                        break_out_of_keystroke_loop = True
                        perform_final_transfer = False
                        break

                    elif command == "t":
                        printed_nothing = True
                        for stage_name in stage_names:
                            if arg is None or arg == stage_name or arg == "all":
                                printed_nothing = False
                                if stage_locations[stage_name] != "local":
                                    print(f"Transferring files back for stage {stage_name}...")
                                    resource = seekrflow.run_settings.get_resource_by_name(
                                        stage_locations[stage_name])
                                    transfer_base.transfer_files_to_from_remote_resource(
                                        seekrflow.name, resource, 
                                        seekrflow.work_directory, backwards=True)
                                else:
                                    print(f"Stage {stage_name} is local; no transfer needed.")
                    
                        if printed_nothing:
                            print(f"No stage found for argument: {arg}.")

                        print(instruction_str, flush=True)
                    else:
                        print(f"Command not recognized: {user_input}.")
                        print(instruction_str, flush=True)

                if break_out_of_keystroke_loop or time.time() >= end_of_loop_time:
                    break
        
            iteration += 1     

    except KeyboardInterrupt:
        print("\nKeyboardInterrupt received, cleaning up...")
        raise
    finally:
        # Only cleanup processes if NOT detaching (perform_final_transfer is our indicator)
        # When detaching (perform_final_transfer=False from 'd' command), leave processes running
        if do_process_cleanup:
            # Cleanup processes on normal exit, error, or Ctrl+C (but not on detach)
            workload_local_mp.cleanup_local_processes(stage_processes, root_directory_path)
        else:
            print("Detaching - local processes will continue running...")
            # Explicitly don't wait for daemon processes
            for stage_name in ["bd", "hidr"]:
                process = stage_processes.get(stage_name)
                if process is not None and process.is_alive():
                    print(f"Leaving {stage_name} process running (PID: {process.pid})")
                    # Don't join - let it keep running
            
            # Handle SEEKR dict of processes
            seekr_processes = stage_processes.get("seekr")
            if isinstance(seekr_processes, dict):
                for anchor_index, process in seekr_processes.items():
                    if process is not None and process.is_alive():
                        print(f"Leaving SEEKR anchor {anchor_index} process running (PID: {process.pid})")
                        
    # Write to the status file
    with open(job_status_filename, "w") as f:
        json.dump(stage_statuses, f, indent=4)

    if perform_final_transfer:
        remotes_transferred_to = set()
        for stage_name in stage_names:
            if (stage_locations[stage_name] != "local") \
                    and (stage_locations[stage_name] not in remotes_transferred_to):
                resource: structures.Resource_remote_base \
                    = seekrflow.run_settings.get_resource_by_name(stage_locations[stage_name])
                print(f"Transferring files back from resource {resource.name}...")
                transfer_base.transfer_files_to_from_remote_resource(
                    seekrflow.name, resource, root_directory, backwards=True)
                remotes_transferred_to.add(resource.name)
    
    if benchmark_mode:
        print("Benchmark mode was enabled - approximate timings can be found in "\
            "hidr_run.out and seekr_run.out files.")

    os.chdir(curdir)
    