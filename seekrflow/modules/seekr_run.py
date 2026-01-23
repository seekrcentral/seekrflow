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
MAX_SUBSEQUENT_NONCOMPLETED_RUNS = 3

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
        seekr_resource_name: str | None = None,
        semaphore_dict: dict[str, str] | None = None,
        mps: int = 1
        ) -> None:
    """
    Run the SEEKR calculation using remote and local resources.
    
    Semaphore values:
        - "go": Normal operation (default)
        - "wait": Don't submit new jobs, but let running jobs finish
        - "stop": Don't submit new jobs AND kill any running jobs
    """
    curdir = os.getcwd()
    if seekrflow.work_directory is not None:
        seekrflow.work_directory = os.path.abspath(seekrflow.work_directory)
        os.chdir(seekrflow.work_directory)
    assert mps >= 1, "MPS must be at least 1"
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
        # Set remote submission times to small values 
        # TODO: move to beginning of seekr stage?
        #for resource in seekrflow.run_settings.resources:
        #    if resource.type == "slurm_remote":
        #        resource.time_limit = "00:30:00"
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
    stage_subsequent_noncompleted_runs = {}
    
    # Track consecutive empty job checks to avoid premature resubmission
    MAX_EMPTY_CHECKS_BEFORE_RESUBMIT = 10
    empty_jobs_check_count = {stage: MAX_EMPTY_CHECKS_BEFORE_RESUBMIT + 1 for stage in stage_names}
    
    # Track quick failures (jobs that start running but fail within 2 cycles)
    previous_stage_manager_status = {stage: "idle" for stage in stage_names}
    stage_running_start_time = {stage: None for stage in stage_names}
    stage_quick_failure_warned = {stage: False for stage in stage_names}
    MIN_RUNNING_TIME_BEFORE_IDLE = 2 * MAIN_LOOP_INTERVAL  # 2 check cycles (60 seconds)
    
    # Initialize semaphore dictionary with defaults
    if semaphore_dict is None:
        semaphore_dict = {"bd": "go", "hidr": "go", "seekr": "go"}
    
    stage_semaphore = semaphore_dict
    
    # Validate semaphore values
    for stage_name, value in stage_semaphore.items():
        if value not in ["go", "wait", "stop"]:
            raise ValueError(f"Invalid semaphore value for {stage_name}: {value}")
    
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
        stage_subsequent_noncompleted_runs[stage_name] = 0

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
    i [STAGE] - print (i)nformation for STAGE (bd, hidr, seekr, all), default: all,
    s [STAGE] - set semaphore to (s)top for STAGE (bd, hidr, seekr, all), kills and prevents jobs, default: all,
    g [STAGE] - set semaphore to (g)o for STAGE (bd, hidr, seekr, all), allows job submission, default: all,
    w [STAGE] - set semaphore to (w)ait for STAGE (bd, hidr, seekr, all), pause new submissions, default: all,
    t [STAGE] - (t)ransfer remote files back for STAGE (bd, hidr, seekr, all),
    q - (q)uit seekrflow, any running jobs will be detached.
    """
    
    try:
        while keep_running:
            # TODO: check each stage based on whether they were submitted or queued in the
            #  previous step. If there's a sudden change (like a job finishing), there might
            #  be an error in configuration or something.
            #print("ITERATION:", iteration)
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
                    
                    # Track consecutive empty job checks for remote stages
                    if stage_locations[stage_name] != "local":
                        jobs = manager_status.get("jobs", []) if manager_status else []
                        empty_jobs_check_count[stage_name] = (
                            empty_jobs_check_count[stage_name] + 1 if len(jobs) == 0 else 0
                        )
            
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
                    
                    # Track consecutive empty job checks for remote stages
                    if stage_locations[stage_name] != "local":
                        jobs = manager_status.get("jobs", []) if manager_status else []
                        empty_jobs_check_count[stage_name] = (
                            empty_jobs_check_count[stage_name] + 1 if len(jobs) == 0 else 0
                        )

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
                    
                    # Track consecutive empty job checks for remote stages
                    if stage_locations[stage_name] != "local":
                        jobs = manager_status.get("jobs", []) if manager_status else []
                        empty_jobs_check_count[stage_name] = (
                            empty_jobs_check_count[stage_name] + 1 if len(jobs) == 0 else 0
                        )
            
            # Detect jobs that start running but fail quickly (within 2 check cycles)
            for stage_name in stage_names:
                if stage_name == "bd" and not using_bd:
                    continue
                
                current_status = stage_manager_status.get(stage_name, "idle")
                previous_status = previous_stage_manager_status[stage_name]
                
                # Track when job enters running/queued state
                if current_status in ["running", "queued", "running/queued"] and previous_status == "idle":
                    stage_running_start_time[stage_name] = time.time()
                    stage_quick_failure_warned[stage_name] = False  # Reset warning flag
                
                # Detect quick failure: running/queued -> idle within MIN_RUNNING_TIME_BEFORE_IDLE
                elif current_status == "idle" \
                    and previous_status in ["running", "queued", "running/queued"] \
                    and stage_semaphore[stage_name] == "go":
                    start_time = stage_running_start_time[stage_name]
                    if start_time is not None:
                        running_duration = time.time() - start_time
                        if running_duration < MIN_RUNNING_TIME_BEFORE_IDLE and not stage_quick_failure_warned[stage_name]:
                            # Warn about quick failure
                            resource = seekrflow.run_settings.get_stage_resource(stage_name)
                            if resource is None:
                                work_dir = root_directory_path
                            else:
                                work_dir = resource.remote_working_directory
                            log_path = os.path.join(work_dir, "logs")
                            print(f"\n{'='*70}")
                            print(f"WARNING: Stage '{stage_name}' ran for only {running_duration:.0f} seconds before disappearing.")
                            print(f"This suggests the job may have failed immediately after starting.")
                            print(f"Common causes:")
                            print(f"  - Missing module or incorrect environment")
                            print(f"  - Syntax error in submission script")
                            print(f"  - Insufficient memory (killed by SLURM)")
                            print(f"  - Bad file permissions or missing files")
                            print(f"Check logs at: {log_path}")
                            print(f"Setting semaphore to 'wait' to prevent resubmission.")
                            print(f"Press 't {stage_name}' to transfer logs, fix the issue, then restart.")
                            print(f"{'='*70}\n")
                            
                            # Set semaphore to "wait" to prevent automatic resubmission
                            stage_semaphore[stage_name] = "wait"
                            stage_quick_failure_warned[stage_name] = True
                        stage_running_start_time[stage_name] = None
                
                # Reset tracking when job runs successfully for longer than threshold
                elif current_status in ["running", "queued", "running/queued"]:
                    start_time = stage_running_start_time[stage_name]
                    if start_time is not None:
                        running_duration = time.time() - start_time
                        if running_duration >= MIN_RUNNING_TIME_BEFORE_IDLE:
                            stage_quick_failure_warned[stage_name] = False  # Reset - this was a legitimate run
                
                # Update previous status
                previous_stage_manager_status[stage_name] = current_status
            
            # Handle semaphore controls
            for stage_name in stage_names:
                if stage_name == "bd" and not using_bd:
                    continue
                
                # Handle "stop" semaphore - kill running jobs immediately
                if stage_semaphore[stage_name] == "stop":
                    if stage_manager_status[stage_name] not in ["idle", "n/a", "waiting"]:
                        print(f"\nSemaphore 'stop' for {stage_name} - killing running jobs...")
                        
                        if stage_locations[stage_name] == "local":
                            # Kill local processes
                            pids_to_kill = get_stage_pids(stage_name, stage_processes, root_directory_path)
                            for pid, anchor_idx in pids_to_kill:
                                kill_process_group(pid, stage_name, root_directory_path, anchor_idx)
                            
                            if stage_statuses.get(stage_name):
                                if stage_statuses[stage_name].get("manager_status"):
                                    stage_statuses[stage_name]["manager_status"]["error"] = "stopped by semaphore"
                                    stage_statuses[stage_name]["manager_status"]["jobs"] = []
                        else:
                            # Kill remote jobs
                            job_ids = stage_job_ids.get(stage_name, [])
                            for job_id in job_ids:
                                if job_id is not None:
                                    print(f"  Canceling job {job_id}...")
                                    workload_remote.submit_remote_cancel_workflow(
                                        seekrflow, stage_locations[stage_name], job_id, silent=silent)
                        
                        stage_manager_status[stage_name] = "idle"

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
                
                # Apply semaphore logic FIRST
                if stage_semaphore[stage_name] in ["wait", "stop"]:
                    # Don't submit new jobs if waiting or stopped
                    stage_should_submit_new_run[stage_name] = False
                    # Count as running if jobs are still active (for "wait" mode)
                    if stage_manager_status[stage_name] not in ["idle", "n/a", "waiting"]:
                        number_of_running_stages += 1
                    continue
                
                # Normal "go" logic
                if stage_progress[stage_name] in ["unstarted", "started"]:
                    stage_should_submit_new_run[stage_name] = True
                else:
                    stage_should_submit_new_run[stage_name] = False
                    stage_manager_status[stage_name] = "idle"

                number_of_incomplete_stages = 0
                if stage_progress[stage_name] not in ["completed", "omitted"]:
                    number_of_incomplete_stages += 1

                # Check dependencies
                for dependency in stage_dependencies[stage_name]:
                    if stage_progress[dependency] != "completed":
                        stage_should_submit_new_run[stage_name] = False
                        stage_manager_status[stage_name] = "waiting"
                        
                if stage_manager_status[stage_name] not in ["idle", "n/a", "waiting"]:
                    stage_should_submit_new_run[stage_name] = False
                    number_of_running_stages += 1
                
                # For remote stages with empty jobs, require multiple checks before resubmit
                if (stage_locations[stage_name] != "local" 
                        and stage_manager_status[stage_name] == "idle"
                        and empty_jobs_check_count[stage_name] <= MAX_EMPTY_CHECKS_BEFORE_RESUBMIT):
                    stage_should_submit_new_run[stage_name] = False
                
                if stage_should_submit_new_run[stage_name]:
                    number_of_stages_to_submit_new_run += 1

            # TODO: remove to force manual exit?
            if number_of_stages_to_submit_new_run == 0 and number_of_running_stages == 0 \
                    and number_of_incomplete_stages == 0:
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
                stage_subsequent_noncompleted_runs[stage_name] += 1
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
                    empty_jobs_check_count["bd"] = 0  # Reset counter after submission
            else:
                if using_bd:
                    stage_subsequent_noncompleted_runs["bd"] = 0
                
            if stage_should_submit_new_run["hidr"]:
                stage_subsequent_noncompleted_runs["hidr"] += 1
                if seekrflow.run_settings.hidr_stage_resource_name != "local":
                    resource: structures.Resource_remote_base \
                        = seekrflow.run_settings.get_resource_by_name(stage_locations["hidr"])
                    destination_path = os.path.join(
                        resource.remote_working_directory, seekrflow.name)
                    destination_model_filename = os.path.join(destination_path, "model.xml")
                    input_pdb_filename = os.path.basename(
                        seekrflow.workflow.solvated_system_for_md.solvated_pdb)
                    #input_pdb_filename = seekrflow.workflow.solvated_system_for_md.solvated_pdb
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
                    empty_jobs_check_count["hidr"] = 0  # Reset counter after submission

            else:
                stage_subsequent_noncompleted_runs["hidr"] = 0
                
            if stage_should_submit_new_run["seekr"]:
                stage_subsequent_noncompleted_runs["seekr"] += 1
                # TODO: what if the resources are different between hidr and seekr stages??
                # gotta perform a transfer. Assert they are the same for now...
                if stage_locations["seekr"] != stage_locations["hidr"]:
                    seekr_resource = seekrflow.run_settings.get_resource_by_name(stage_locations["seekr"])
                    hidr_resource = seekrflow.run_settings.get_resource_by_name(stage_locations["hidr"])
                    transfer_base.transfer_files_to_from_remote_resource(
                        seekrflow.name, hidr_resource,
                        root_directory, backwards=True)
                    transfer_base.transfer_files_to_from_remote_resource(
                        seekrflow.name, seekr_resource, root_directory)

                if seekrflow.run_settings.seekr_stage_resource_name != "local":
                    resource: structures.Resource_remote_base \
                        = seekrflow.run_settings.get_resource_by_name(stage_locations["seekr"])
                    if benchmark_mode:
                        resource.time_limit = "00:30:00"
                    destination_path = os.path.join(
                        resource.remote_working_directory, seekrflow.name)
                    destination_model_filename = os.path.join(destination_path, "model.xml")
                    if benchmark_mode:
                        benchmark_string = f"-t {base.BENCHMARK_MIN_SIMULATION_LENGTH} "
                    else:
                        benchmark_string = ""
                    if mps == 1:
                        output_filename = os.path.join(destination_path, "seekr_{index}_run.out")
                        command_string = f"python {resource.remote_seekr2_directory}/run.py "\
                            f"{{index}} model.xml -c 0 {benchmark_string}> {output_filename}"
                    else:
                        command_list = []
                        for mps_index in range(mps):
                            output_filename = os.path.join(destination_path, f"seekr_{{index{mps_index}}}_run.out")
                            command_list.append(
                                f"if [ {mps_index} -lt {{len_array}} ]; then "
                                f"python {resource.remote_seekr2_directory}/run.py "\
                                f"{{index{mps_index}}} model.xml -c 0 {benchmark_string}> {output_filename} & "
                                f"fi"
                            )
                        command_list.append("wait")
                        command_string = ";".join(command_list)
                        
                    indices = seekr_status["stage_status"]["incomplete_anchors"]
                    seekr_run_result = workload_remote.submit_remote_run_workflow(
                        seekrflow, "seekr", destination_path, resource,
                        command_string, destination_model_filename, workflow_type="seekr", 
                        indices=indices, mps=mps)
                    stage_run_results["seekr"] = seekr_run_result
                    stage_manager_status["seekr"] = "running"
                    stage_progress["seekr"] = "started"
                    empty_jobs_check_count["seekr"] = 0  # Reset counter after submission

            else:
                stage_subsequent_noncompleted_runs["seekr"] = 0

            
            for stage_name in stage_names:
                if stage_subsequent_noncompleted_runs[stage_name] \
                        >= MAX_SUBSEQUENT_NONCOMPLETED_RUNS:
                    raise RuntimeError(
                        f"Stage {stage_name} has been submitted "
                        f"{stage_subsequent_noncompleted_runs[stage_name]} times without progress. "
                        f"Please check the logs for errors.")
                
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
                #input_pdb_file = seekrflow.workflow.solvated_system_for_md.solvated_pdb
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
                for mps_index, unfilled_anchor_index in enumerate(
                        seekr_status["stage_status"]["incomplete_anchors"]):
                    # Start process for this anchor
                    process = multiprocessing.Process(
                        target=workload_local_mp.run_seekr_locally,
                        args=(model_filename, unfilled_anchor_index, root_directory_path, benchmark_mode),
                        daemon=False
                    )
                    process.start()
                    stage_processes["seekr"][unfilled_anchor_index] = process
                    # In an MPS process, launch an 'mps' number of concurrent jobs.
                    if mps_index >= mps-1:
                        break
            
                # NOTE: this probably must be done for local job where there's just one GPU
                # but for the cluster, we will want to spin them off all at once.
                stage_manager_status["seekr"] = "running"
                stage_progress["seekr"] = "started"
        
            # Now pause to accept user keystrokes
            if iteration == 0:
                print(instruction_str, flush=True)
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
                                print(f"  Semaphore: {stage_semaphore[stage_name]}")
                                print(f"  Should perform transfer: {stage_should_perform_transfer[stage_name]}")
                                printed_nothing = False
                    
                        if printed_nothing:
                            print(f"No stage found for argument: {arg}.")

                        print(instruction_str, flush=True)
                    
                    elif command == "q":
                        print(f"Quitting seekrflow. Running jobs will be detached...")
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

                    elif command == "s":
                        # Determine which stages to stop
                        stages_to_stop = []
                        if arg is None or arg == "all":
                            # Stop all stages
                            for stage_name in stage_names:
                                if stage_name == "bd" and not using_bd:
                                    continue
                                stages_to_stop.append(stage_name)
                        elif arg in stage_names:
                            # Stop specific stage
                            if arg == "bd" and not using_bd:
                                print(f"Stage {arg} is not being used.")
                                print(instruction_str, flush=True)
                                continue
                            stages_to_stop.append(arg)
                        else:
                            print(f"Invalid stage name: {arg}. Must be bd, hidr, seekr, or all.")
                            print(instruction_str, flush=True)
                            continue
                        
                        if not stages_to_stop:
                            print("No stages to stop.")
                            print(instruction_str, flush=True)
                            continue
                        
                        # Set semaphores to "stop" for selected stages
                        for stage_name in stages_to_stop:
                            stage_semaphore[stage_name] = "stop"
                        
                        # Immediately execute the semaphore "stop" logic for selected stages
                        for stage_name in stages_to_stop:
                            if stage_manager_status[stage_name] not in ["idle", "n/a", "waiting"]:
                                print(f"Killing stage {stage_name}...")
                                
                                if stage_locations[stage_name] == "local":
                                    # Kill local processes
                                    pids_to_kill = get_stage_pids(stage_name, stage_processes, root_directory_path)
                                    for pid, anchor_idx in pids_to_kill:
                                        kill_process_group(pid, stage_name, root_directory_path, anchor_idx)
                                    
                                    if stage_statuses.get(stage_name):
                                        if stage_statuses[stage_name].get("manager_status"):
                                            stage_statuses[stage_name]["manager_status"]["error"] = "stopped by user"
                                            stage_statuses[stage_name]["manager_status"]["jobs"] = []
                                else:
                                    # Kill remote jobs
                                    job_ids = stage_job_ids.get(stage_name, [])
                                    for job_id in job_ids:
                                        if job_id is not None:
                                            print(f"  Canceling job {job_id}...")
                                            workload_remote.submit_remote_cancel_workflow(
                                                seekrflow, stage_locations[stage_name], job_id, silent=silent)
                                
                                stage_manager_status[stage_name] = "idle"
                            else:
                                print(f"Stage {stage_name} is already {stage_manager_status[stage_name]}.")
                        
                        # Only exit if all stages were stopped
                        #active_stage_count = sum(1 for s in stage_names if not (s == "bd" and not using_bd))
                        #if len(stages_to_stop) == active_stage_count:
                        #    keep_running = False
                        #    break_out_of_keystroke_loop = True
                        #    perform_final_transfer = False
                        #    break
                        #else:
                        print(instruction_str, flush=True)

                    elif command == "g":
                        # Determine which stages to set to "go"
                        stages_to_go = []
                        if arg is None or arg == "all":
                            for stage_name in stage_names:
                                if stage_name == "bd" and not using_bd:
                                    continue
                                stages_to_go.append(stage_name)
                        elif arg in stage_names:
                            if arg == "bd" and not using_bd:
                                print(f"Stage {arg} is not being used.")
                                print(instruction_str, flush=True)
                                continue
                            stages_to_go.append(arg)
                        else:
                            print(f"Invalid stage name: {arg}. Must be bd, hidr, seekr, or all.")
                            print(instruction_str, flush=True)
                            continue
                        
                        # Set semaphores to "go" for selected stages
                        for stage_name in stages_to_go:
                            stage_semaphore[stage_name] = "go"
                            print(f"Stage {stage_name} semaphore set to 'go'.")
                        
                        print(instruction_str, flush=True)

                    elif command == "w":
                        # Determine which stages to set to "wait"
                        stages_to_wait = []
                        if arg is None or arg == "all":
                            for stage_name in stage_names:
                                if stage_name == "bd" and not using_bd:
                                    continue
                                stages_to_wait.append(stage_name)
                        elif arg in stage_names:
                            if arg == "bd" and not using_bd:
                                print(f"Stage {arg} is not being used.")
                                print(instruction_str, flush=True)
                                continue
                            stages_to_wait.append(arg)
                        else:
                            print(f"Invalid stage name: {arg}. Must be bd, hidr, seekr, or all.")
                            print(instruction_str, flush=True)
                            continue
                        
                        # Set semaphores to "wait" for selected stages
                        for stage_name in stages_to_wait:
                            stage_semaphore[stage_name] = "wait"
                            print(f"Stage {stage_name} semaphore set to 'wait'.")
                        
                        print(instruction_str, flush=True)

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
                                        root_directory, backwards=True)
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
    