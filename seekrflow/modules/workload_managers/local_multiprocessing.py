"""
modules/workload_managers/local_multiprocessing.py

Local multiprocessing workload manager for seekrflow.
"""

import os
import glob
import time
import json
import signal
import multiprocessing
from dataclasses import dataclass, asdict
from typing import Optional, Dict

import seekr2.modules.common_base as seekr2_base
import seekr2.run as seekr2_run

import seekrflow.modules.base as base

# ============================================================================
# Multiprocessing State Management (similar to SLURM RunState pattern)
# ============================================================================

@dataclass
class LocalProcessState:
    """State information for a locally running multiprocessing process."""
    stage: str
    pid: int
    status: str  # "running", "completed", "failed", "killed", "unknown"
    started_at: float
    ended_at: Optional[float]  # Set when process completes, fails, or is killed
    error: Optional[str]
    traceback: Optional[str]
    model_file: str
    work_dir: str
    # Stage-specific progress information
    progress_info: Optional[Dict] = None
    
    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return asdict(self)
    
    @staticmethod
    def from_dict(data: dict) -> "LocalProcessState":
        """Create from dictionary loaded from JSON."""
        return LocalProcessState(**data)
    
    @staticmethod
    def load(path: str) -> "LocalProcessState":
        """Load state from file."""
        with open(path, "r") as f:
            data = json.load(f)
        return LocalProcessState.from_dict(data)
    
    def save(self, path: str) -> None:
        """Save state to file."""
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=4)


def get_local_state_dir(root_dir: str) -> str:
    """
    Get the directory for local multiprocessing state files.
    Located next to model.xml in the root directory.
    """
    state_dir = os.path.join(
        root_dir, ".multiprocessing")
    os.makedirs(state_dir, exist_ok=True)
    return state_dir

def get_local_state_file(
        root_dir: str, 
        stage: str, 
        pid: int,
        anchor: str = "any",
        swarm_id: int | None = None
        ) -> str:
    """
    Get the state file path for a local process.
    Includes PID in filename to avoid conflicts.
    """
    state_dir = get_local_state_dir(root_dir)
    anchor_str = ""
    if anchor != "any":
        anchor_str = f"anchor_{anchor}_"

    swarm_id_str = ""
    if swarm_id is not None:
        swarm_id_str = f"swarm_{swarm_id}_"

    return os.path.join(
        state_dir, f"{stage}_{anchor_str}{swarm_id_str}state_{pid}.json")


def check_pid_exists(pid: int) -> bool:
    """
    Check if a process with given PID is still running.
    
    Uses os.kill with signal 0, which doesn't actually send a signal
    but checks if the process exists.
    """
    if pid is None or pid <= 0:
        return False
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False

def find_latest_state_file(
        root_dir: str, 
        stage_name: str,
        anchor: str = "any",
        swarm_id: int | None = None
        ) -> Optional[str]:
    """
    Find the most recent state file for a given stage.
    Returns None if no state files exist for this stage.
    """
    state_dir = get_local_state_dir(root_dir)
    anchor_str = ""
    if anchor != "any":
        anchor_str = f"anchor_{anchor}_"
    swarm_id_str = ""
    if swarm_id is not None:
        swarm_id_str = f"swarm_{swarm_id}_"
    pattern = f"{stage_name}_{anchor_str}{swarm_id_str}state_*.json"
    state_glob = os.path.join(state_dir, pattern)
    state_files = list(glob.glob(state_glob))
    if not state_files:
        return None
    state_files.sort(key=lambda p: os.path.getmtime(p), reverse=True)
    return state_files[0]

def cleanup_old_state_files(
        root_dir: str, 
        stage_name: str, 
        anchor: str = "any",
        swarm_id: int | None = None,
        keep_latest: bool = True):
    """
    Clean up old state files for a stage, optionally keeping the latest one.
    """
    state_dir = get_local_state_dir(root_dir)
    anchor_str = ""
    if anchor != "any":
        anchor_str = f"anchor_{anchor}_"
    swarm_id_str = ""
    if swarm_id is not None:
        swarm_id_str = f"swarm_{swarm_id}_"
    pattern = f"{stage_name}_{anchor_str}{swarm_id_str}state_*.json"
    state_files = list(state_dir.glob(pattern))
    
    if not state_files:
        return
    
    # Sort by modification time
    state_files.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    
    # Determine which files to delete
    files_to_delete = state_files[1:] if keep_latest else state_files
    
    for state_file in files_to_delete:
        try:
            state_file.unlink()
        except Exception as e:
            print(f"Warning: Could not delete old state file {state_file}: {e}")

def check_for_existing_local_processes(
        root_dir: str, 
        stage_name: str,
        anchor: str = "any",
        swarm_id: int | None = None
        ) -> Optional[LocalProcessState]:
    """
    Check if there's a process from a previous run still active.
    Returns the process state if found and running, None otherwise.
    
    This is used for reattachment after detaching with 'd' command,
    and for monitoring.
    """
    state_file = find_latest_state_file(
        root_dir, stage_name, anchor, swarm_id)
    if not state_file or not os.path.exists(state_file):
        return None
    
    try:
        state = LocalProcessState.load(state_file)
        
        # Only consider it if it's marked as running
        if state.status == "running":
            # Verify the PID actually exists
            if check_pid_exists(state.pid):
                print(f"Found existing {stage_name} process from previous "\
                      f"session (PID: {state.pid})")
                return state
            else:
                # Process died without updating state file
                print(f"Found stale {stage_name} state file - " \
                      f"process {state.pid} no longer exists")
                state.status = "unknown"
                state.ended_at = time.time()
                state.error = "Process terminated unexpectedly"
                state.save(state_file)
                return None
        
    except Exception as e:
        print(f"Warning: Could not load state file for {stage_name}: {e}")
    
    return None

def status_local(
        root_dir: str,
        stage_index: int,
        stage_name: str,
        stage_process: multiprocessing.Process | None,
        anchor: str = "any",
        swarm_id: int | None = None,
        ) -> dict: 
    """
    Check if the BD stage has finished locally. This is done by:
    1. Checking the process state file
    2. Checking if the process is still alive
    3. Reading the results XML files to see how many BD steps have elapsed
    """
    import os
    import sys
    import traceback
    
    # Check for state file
    state_file = find_latest_state_file(root_dir, stage_name)
    process_state = None
    
    if state_file and os.path.exists(state_file):
        try:
            process_state = LocalProcessState.load(state_file)
        except Exception as e:
            print(f"Warning: Could not load BD state file: {e}")
    
    benchmark_mode = False
    if process_state is not None and process_state.progress_info is not None:
        benchmark_mode = bool(
            process_state.progress_info.get("benchmark_mode", False))
    
    # TODO: perhaps extract this from the state file instead?
    if anchor == "any":
        partitioned_arg = None
    else:
        partitioned_arg = anchor

    stage_status = {}
    stage_status["model_xml_found"] = True
    stage_status["notes"] = ""

    if benchmark_mode:
        # status.py doesn't know benchmark-mode criteria. For benchmark runs,
        # derive state from the process-state lifecycle and local process liveness.
        running_now = False
        if stage_process is not None and stage_process.is_alive():
            running_now = True
        elif process_state is not None \
                and process_state.status == "running" \
                and check_pid_exists(process_state.pid):
            running_now = True

        if process_state is not None and process_state.status == "completed":
            stage_status["finished"] = True
            stage_status["state"] = "completed"
            stage_status["progress"] = 1.0
        elif running_now:
            stage_status["finished"] = False
            stage_status["state"] = "started"
            stage_status["progress"] = 0.0
        else:
            stage_status["finished"] = False
            stage_status["state"] = "unstarted"
            stage_status["progress"] = 0.0
    else:
        # Check actual simulation progress from model files
        model_filename = os.path.join(root_dir, "model.json")

        try:
            # Import here so child process gets fresh imports
            import seekr.modules.structures as structures
            import seekr.status as seekr_status

            # Only care about stage progress, not analysis statistics
            model = structures.load_model(model_filename)
            instruction = "progress"
            message_dict = seekr_status.status(
                model,
                instruction,
                stage_arg=stage_index,
                # seekr.status/run.extract_anchor_index expects "any" or a
                # concrete anchor identifier; passing None raises TypeError.
                anchor_arg=anchor,
                swarm_id=swarm_id,
                print_json=True,
            )
            progress_by_stage = message_dict.get("progress", {})
            stage_progress = progress_by_stage.get(stage_index)
            if stage_progress is None:
                # Be tolerant of key-type differences from status() output.
                stage_progress = progress_by_stage.get(str(stage_index), {})

        except Exception as e:
            # Update state to failed (process itself updates this)
            print(f"\nStage {stage_index} status check failed: {e}")
            traceback.print_exc()
            sys.stdout.flush()
            raise  # Re-raise so process exits with error code

        if stage_progress.get("finished", False):
            stage_status["finished"] = True
            stage_status["state"] = "completed"
            stage_status["progress"] = 1.0
        else:
            progress_map = stage_progress.get("progress", {})
            if partitioned_arg is None:
                # "any" means stage-level view. seekr.status stores progress
                # per-partition (e.g., per-anchor) rather than under key None,
                # so aggregate all partition progress values.
                progress_values = []
                if isinstance(progress_map, dict):
                    for _, per_partition in progress_map.items():
                        if isinstance(per_partition, dict):
                            value = per_partition.get("progress")
                            if isinstance(value, (int, float)):
                                progress_values.append(float(value))
                if progress_values:
                    progress = sum(progress_values) / len(progress_values)
                else:
                    progress = 0.0
                partitioned_status = {}
            else:
                partitioned_status = progress_map.get(partitioned_arg, {})
                # Empty dict = no progress data written yet (run just started, or
                # the partitioned_arg isn't present in seekr's output yet). Treat
                # that as 0.0 progress rather than an invariant violation.
                if partitioned_status and "progress" not in partitioned_status:
                    raise RuntimeError(
                        f"seekr status output for stage {stage_index} "
                        f"(stage_name={stage_name!r}, partition={partitioned_arg!r}) "
                        f"is missing the required 'progress' key. Got keys: "
                        f"{sorted(partitioned_status.keys())}")
                progress = partitioned_status.get("progress", 0.0)

            if progress > 0.0:
                state = "started"
            else:
                state = "unstarted"
            stage_status["finished"] = False
            stage_status["state"] = state
            stage_status["progress"] = progress
            for key, value in partitioned_status.items():
                if key != "progress":
                    stage_status[key] = value

    # seekr.status() above reflects model.json, which may still show a
    # previous run as "finished" until the freshly-spawned child clears it
    # via force_overwrite. The orchestrator's own Process handle is the
    # authoritative source for whether the CURRENT run is still active:
    # if it is alive, the stage cannot possibly be completed. The state
    # file is unreliable here because, just after a force_rerun, the
    # latest-by-mtime file is the OLD "killed" record from the prior run
    # (the newly-spawned child hasn't written its state file yet due to
    # "spawn" import overhead).
    current_run_alive = False
    if stage_process is not None and stage_process.is_alive():
        current_run_alive = True
    elif process_state is not None \
            and process_state.status == "running" \
            and check_pid_exists(process_state.pid):
        current_run_alive = True

    if current_run_alive and stage_status["state"] == "completed":
        stage_status["finished"] = False
        stage_status["state"] = "started"
    elif current_run_alive and stage_status["state"] == "unstarted":
        # During partitioned runs, seekr.status can briefly return no progress
        # data even though the process is actively producing outputs.
        stage_status["state"] = "started"

    # If the child process finished cleanly (seekr_run.run returned without
    # raising), trust that terminal state only for benchmark runs. In BD
    # benchmark mode, do_run_instruction_bd intentionally wipes BD outputs
    # afterwards, so seekr.status can report finished=False/progress=0.0 and
    # the progress-derived state would otherwise pin to "unstarted" forever.
    if benchmark_mode \
            and not current_run_alive \
            and process_state is not None \
            and process_state.status == "completed" \
            and stage_status["state"] != "completed":
        stage_status["finished"] = True
        stage_status["state"] = "completed"
        # Preserve the underlying progress number so the UI can still show
        # what seekr currently sees (e.g. 0.0 after benchmark cleanup).
    
    # Build manager status based on process state
    stage_manager_status = {}
    stage_manager_status["tool"] = "multiprocessing"
    stage_manager_status["timestamp"] = time.time()
    stage_manager_status["error"] = ""
    stage_manager_status["jobs"] = []
    
    # Check if process is running and build job object with integrated state
    if stage_process is not None and stage_process.is_alive():
        job = {
            "JobID": f"local_pid_{stage_process.pid}",
            "State": "RUNNING",
            "PID": stage_process.pid
        }
        # Integrate process state info into job object
        if process_state:
            job["process_status"] = process_state.status
            job["started_at"] = process_state.started_at
            job["output_file"] = process_state.progress_info.get(
                "output_file", "N/A") if process_state.progress_info \
                else "N/A"
        stage_manager_status["jobs"].append(job)
    elif process_state:
        # Check state file
        if process_state.status == "running":
            # Process claims to be running, verify
            if check_pid_exists(process_state.pid):
                job = {
                    "JobID": f"local_pid_{process_state.pid}",
                    "State": "RUNNING",
                    "PID": process_state.pid,
                    "process_status": process_state.status,
                    "started_at": process_state.started_at,
                    "output_file": process_state.progress_info.get(
                        "output_file", "N/A") if process_state.progress_info \
                        else "N/A"
                }
                stage_manager_status["jobs"].append(job)
            else:
                # Process died without updating state file
                stage_status["notes"] = f"Process {process_state.pid} "\
                    "exited unexpectedly"
                stage_manager_status["error"] = "Process crashed"
        elif process_state.status in ["failed", "killed"]:
            stage_status["notes"] = f"Process {process_state.status}: "\
                f"{process_state.error}"
            if process_state.status == "failed":
                stage_manager_status["error"] = process_state.error
            # DON'T add job to jobs list - these processes are done and 
            # not running, this allows stage_manager_status to be "idle" 
            # so a new run can be submitted
    
    results = {
        "success": True,
        "error": None,
        "manager_status": stage_manager_status,
        "stage_status": stage_status
    }
    
    return results
     
# ============================================================================
# Local Execution Functions with State Management
# ============================================================================

def run_locally(
        root_dir: str,
        stage_name: str,
        anchor: str = "any",
        swarm_id: int | None = None,
        device_index: str | None = None,
        force_overwrite: bool = False,
        benchmark_mode: bool = False,
        #restart_attempts: int = 0,
        ):
    """
    Run seekr stage locally using multiprocessing.
    Writes state file with PID for monitoring and reattachment.
    Redirects stdout/stderr to {stage}_run.out.
    
    This function is executed in a separate process.
    """
    # TODO: handle benchmark mode?

    import os
    import sys
    import time
    import traceback
    import signal
    
    # Make this process a process group leader so all subprocesses 
    # (like openmm / nam_simulation) are in the same group and can 
    # be killed together
    os.setpgrp()
    
    pid = os.getpid()
    # TODO: check if this is a file or a dir
    state_file = get_local_state_file(
        root_dir, stage_name, pid, anchor, swarm_id)
    anchor_str = ""
    if anchor != "any":
        anchor_str = f"anchor_{anchor}_"

    swarm_id_str = ""
    if swarm_id is not None:
        swarm_id_str = f"swarm_id_{swarm_id}_"

    output_file = os.path.join(
        root_dir, f"{stage_name}_{anchor_str}{swarm_id_str}run.out")
    model_filename = os.path.join(root_dir, "model.json")
    
    # Create initial state BEFORE redirecting output (so any errors 
    # are visible)
    state = LocalProcessState(
        stage=stage_name,
        pid=pid,
        status="running",
        started_at=time.time(),
        ended_at=None,
        error=None,
        traceback=None,
        model_file=model_filename,
        work_dir=root_dir,
        progress_info={
            "output_file": str(output_file),
            "anchor": anchor,
            "swarm_id": swarm_id,
            "device_index": device_index,
            "force_overwrite": force_overwrite,
            "benchmark_mode": benchmark_mode,
            "imports_successful": False,
            "model_loaded": False,
            "finished": False
        }
    )
    
    try:
        # Write initial state file
        state.save(state_file)
        
        # Create/update symlink to latest state file
        latest_link = os.path.join(
            root_dir, ".multiprocessing", f"{stage_name}_{anchor_str}latest.json")
        if os.path.exists(latest_link) or os.path.islink(latest_link):
            os.unlink(latest_link)
        os.symlink(state_file, latest_link)
        
    except Exception as e:
        # If we can't write state file, fail loudly BEFORE redirecting output
        print(f"FATAL: Cannot create state file {state_file}: {e}", 
              file=sys.stderr)
        traceback.print_exc()
        raise
    
    # NOW redirect output at OS level using dup2 (so subprocesses via 
    # os.system() are also redirected)
    try:
        output_fd = os.open(
            str(output_file), os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
        os.dup2(output_fd, 1)  # Redirect stdout (fd 1)
        os.dup2(output_fd, 2)  # Redirect stderr (fd 2)
        os.close(output_fd)
        
        # Also redirect at Python level for consistency
        # Line buffered, append mode
        sys.stdout = open(output_file, 'a', buffering=1) 
        sys.stderr = sys.stdout
    except Exception as e:
        print(f"FATAL: Cannot redirect output to {output_file}: {e}", 
              file=sys.stderr)
        raise
    
    print(f"{stage_name} stage started at {time.ctime()}")
    print(f"PID: {pid}")
    print(f"Model file: {model_filename}")
    print(f"State file: {state_file}")
    print(f"Anchor(s): {anchor}")
    print(f"Device index: {device_index}")
    print(f"Force overwrite: {force_overwrite}")
    print(f"Benchmark mode: {benchmark_mode}")
    print(f"Swarm ID: {swarm_id}")
    print("-" * 60)
    sys.stdout.flush()
    
    # Handle termination signals to update state before exit
    def signal_handler(signum, frame):
        print(f"\nReceived signal {signum}, updating state to 'killed'...")
        sys.stdout.flush()
        state.status = "killed"
        state.ended_at = time.time()
        state.error = f"Process killed by signal {signum}"
        state.save(state_file)
        sys.exit(128 + signum)  # Standard exit code for signals
    
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)
    
    try:
        # Import here so child process gets fresh imports
        import seekr.modules.structures as structures
        import seekr.run as seekr_run
        
        # Update state - imports successful
        state.progress_info["imports_successful"] = True
        state.save(state_file)
        print("Imports successful")
        sys.stdout.flush()
        
        # Load model
        model = structures.load_model(model_filename)
        state.progress_info["model_loaded"] = True
        state.save(state_file)
        print("Model loaded")
        sys.stdout.flush()
        
        # Run stage
        print(f"Starting {stage_name}...")
        sys.stdout.flush()
        seekr_run.run(
            model, stage_name, anchor, device_index, force_overwrite, 
            swarm_id, benchmark=benchmark_mode)
        
        # Update state to completed (process itself updates this)
        state.status = "completed"
        state.ended_at = time.time()
        state.progress_info["finished"] = True
        state.save(state_file)
        print(f"\n{stage_name} completed at {time.ctime()}")
        sys.stdout.flush()
        
    except Exception as e:
        # Update state to failed (process itself updates this)
        state.status = "failed"
        state.error = str(e)
        state.traceback = traceback.format_exc()
        state.ended_at = time.time()
        state.save(state_file)
        print(f"\n{stage_name} failed: {e}")
        traceback.print_exc()
        sys.stdout.flush()
        raise  # Re-raise so process exits with error code

    return

def kill_existing_local_stage_processes(
        stage_name: str,
        root_directory_path: str,
        ) -> None:
    """
    Kill any running local processes for the given stage, identified via
    state files in the .multiprocessing directory. Idempotent.
    """
    state_dir = get_local_state_dir(root_directory_path)
    pattern = "seekr_anchor_*_state_*.json" if stage_name == "seekr" \
        else f"{stage_name}_state_*.json"
    for state_file in glob.glob(os.path.join(state_dir, pattern)):
        try:
            state = LocalProcessState.load(state_file)
        except Exception as e:
            print(f"  Warning: Could not load state file {state_file}: {e}")
            continue
        if state.status != "running" or not check_pid_exists(state.pid):
            continue
        pid = state.pid
        print(f"  Killing existing {stage_name} process (PID: {pid})...")
        try:
            os.killpg(pid, signal.SIGTERM)
            time.sleep(2)
            if check_pid_exists(pid):
                os.killpg(pid, signal.SIGKILL)
                time.sleep(1)
            state.status = "killed"
            state.ended_at = time.time()
            state.error = "Process killed for force re-run"
            state.save(state_file)
        except (ProcessLookupError, PermissionError) as e:
            print(f"  Warning: Could not kill process: {e}")
    return

def check_and_raise_if_process_failed(
        stage_name: str, 
        stage_process: multiprocessing.Process | None, 
        root_directory_path: str,
        anchor: str = "any",
        swarm_id: int | None = None
        ) -> None:
    """
    Check if a local process has failed by examining state files.
    Raises an exception if failure detected.
    """
    if stage_process is not None and not stage_process.is_alive():
        # Check exit code
        if hasattr(stage_process, 'exitcode') \
                and stage_process.exitcode is not None \
                and stage_process.exitcode != 0:
            state_file = find_latest_state_file(
                root_directory_path, stage_name, anchor, swarm_id)
            if state_file:
                try:
                    state = LocalProcessState.load(state_file)
                    # Trust the state file: if the child reported success,
                    # ignore a noisy nonzero exit code from interpreter
                    # shutdown (common with fork + asyncio/threads).
                    if state.status == "completed":
                        return
                    error_text = state.error if state.error is not None \
                        else f"process exited with code " \
                             f"{stage_process.exitcode} while state was " \
                             f"still '{state.status}'"
                    error_msg = f"Local {stage_name.upper()} stage failed. "\
                                f"Error: {error_text}"
                    if state.traceback:
                        error_msg += f"\nTraceback:\n{state.traceback}"
                    raise Exception(error_msg)
                except json.JSONDecodeError:
                    raise Exception(f"Local {stage_name.upper()} stage failed "\
                                    f"with exit code {stage_process.exitcode}")
            else:
                raise Exception(f"Local {stage_name.upper()} stage failed "\
                                f"with exit code {stage_process.exitcode}")
