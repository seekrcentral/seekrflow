"""
modules/workload_managers/local_multiprocessing.py

Local multiprocessing workload manager for seekrflow.
"""

import os
import sys
import glob
import time
import json
import signal
import pathlib
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
    stage: str  # "bd", "hidr", or "seekr"
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
    def load(path: pathlib.Path) -> "LocalProcessState":
        """Load state from file."""
        with open(path, "r") as f:
            data = json.load(f)
        return LocalProcessState.from_dict(data)
    
    def save(self, path: pathlib.Path) -> None:
        """Save state to file."""
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=4)


# TODO: send to local_mp
def get_local_state_dir(root_dir: pathlib.Path) -> pathlib.Path:
    """
    Get the directory for local multiprocessing state files.
    Located next to model.xml in the root directory.
    """
    state_dir = root_dir / ".multiprocessing"
    state_dir.mkdir(exist_ok=True)
    return state_dir

# TODO: send to local_mp
def get_local_state_file(
        root_dir: pathlib.Path, 
        stage: str, pid: int
        ) -> pathlib.Path:
    """
    Get the state file path for a local process.
    Includes PID in filename to avoid conflicts.
    """
    state_dir = get_local_state_dir(root_dir)
    return state_dir / f"{stage}_state_{pid}.json"


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

# TODO: send to local_mp
def find_latest_state_file(
        root_dir: pathlib.Path, 
        stage: str
        ) -> Optional[pathlib.Path]:
    """
    Find the most recent state file for a given stage.
    Returns None if no state files exist for this stage.
    """
    state_dir = get_local_state_dir(root_dir)
    if stage == "seekr":
        pattern = "seekr_anchor_*_state_*.json"
    else:
        pattern = f"{stage}_state_*.json"
    state_files = list(state_dir.glob(pattern))
    if not state_files:
        return None
    state_files.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return state_files[0]

# TODO: send to local_mp
def cleanup_old_state_files(
        root_dir: pathlib.Path, 
        stage: str, 
        keep_latest: bool = True):
    """
    Clean up old state files for a stage, optionally keeping the latest one.
    """
    state_dir = get_local_state_dir(root_dir)
    if stage == "seekr":
        pattern = "seekr_anchor_*_state_*.json"
    else:
        pattern = f"{stage}_state_*.json"
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

# TODO: send to local_mp
def check_for_existing_local_processes(
        root_dir: pathlib.Path, 
        stage: str
        ) -> Optional[LocalProcessState]:
    """
    Check if there's a process from a previous run still active.
    Returns the process state if found and running, None otherwise.
    
    This is used for reattachment after detaching with 'd' command.
    """
    state_file = find_latest_state_file(root_dir, stage)
    if not state_file or not state_file.exists():
        return None
    
    try:
        state = LocalProcessState.load(state_file)
        
        # Only consider it if it's marked as running
        if state.status == "running":
            # Verify the PID actually exists
            if check_pid_exists(state.pid):
                print(f"Found existing {stage} process from previous session (PID: {state.pid})")
                return state
            else:
                # Process died without updating state file
                print(f"Found stale {stage} state file - process {state.pid} no longer exists")
                state.status = "unknown"
                state.ended_at = time.time()
                state.error = "Process terminated unexpectedly"
                state.save(state_file)
                return None
        
    except Exception as e:
        print(f"Warning: Could not load state file for {stage}: {e}")
    
    return None

def bd_status_local(
        model: seekr2_base.Model,
        stage_processes: dict[str, multiprocessing.Process],
        work_dir: pathlib.Path,
        n_trajectories_for_completion: int
        ) -> dict: 
    """
    Check if the BD stage has finished locally. This is done by:
    1. Checking the process state file
    2. Checking if the process is still alive
    3. Reading the results XML files to see how many BD steps have elapsed
    """
    
    # Check for state file
    state_file = find_latest_state_file(work_dir, "bd")
    process_state = None
    
    if state_file and state_file.exists():
        try:
            process_state = LocalProcessState.load(state_file)
        except Exception as e:
            print(f"Warning: Could not load BD state file: {e}")
    
    # Check actual simulation progress from model files
    bd_milestone_info_to_run = seekr2_run.choose_next_simulation_browndye2(
            model, "any_bd", None, False, None)
    
    stage_status = {}
    stage_status["model_xml_found"] = True
    stage_status["notes"] = ""
    stage_status["n_trajectories_for_completion"] = n_trajectories_for_completion
    
    if len(bd_milestone_info_to_run) == 0:
        stage_status["finished"] = True
        stage_status["progress"] = "completed"
        stage_status["n_trajectories"] = n_trajectories_for_completion
        stage_status["progress_bar"] = 1.0
    else:
        steps_to_go, total_b_surface_encounters, surface_name, restart, \
            total_b_surface_sims = bd_milestone_info_to_run[0]
        if restart:
            progress = "started"
        else:
            progress = "unstarted"
        stage_status["finished"] = False
        stage_status["progress"] = progress
        stage_status["n_trajectories"] = n_trajectories_for_completion - steps_to_go
        stage_status["progress_bar"] = (n_trajectories_for_completion - steps_to_go) / n_trajectories_for_completion
    
    # Build manager status based on process state
    stage_manager_status = {}
    stage_manager_status["tool"] = "multiprocessing"
    stage_manager_status["timestamp"] = time.time()
    stage_manager_status["error"] = ""
    stage_manager_status["jobs"] = []
    
    # Check if process is running and build job object with integrated state
    process = stage_processes.get("bd")
    if process is not None and process.is_alive():
        job = {
            "JobID": f"local_pid_{process.pid}",
            "State": "RUNNING",
            "PID": process.pid
        }
        # Integrate process state info into job object
        if process_state:
            job["process_status"] = process_state.status
            job["started_at"] = process_state.started_at
            job["output_file"] = process_state.progress_info.get("output_file", "N/A") if process_state.progress_info else "N/A"
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
                    "output_file": process_state.progress_info.get("output_file", "N/A") if process_state.progress_info else "N/A"
                }
                stage_manager_status["jobs"].append(job)
            else:
                # Process died without updating state file
                stage_status["notes"] = f"Process {process_state.pid} exited unexpectedly"
                stage_manager_status["error"] = "Process crashed"
        elif process_state.status in ["failed", "killed"]:
            stage_status["notes"] = f"Process {process_state.status}: {process_state.error}"
            if process_state.status == "failed":
                stage_manager_status["error"] = process_state.error
            # DON'T add job to jobs list - these processes are done and not running
            # This allows stage_manager_status to be "idle" so a new run can be submitted
    
    bd_results = {
        "success": True,
        "error": None,
        "manager_status": stage_manager_status,
        "stage_status": stage_status
    }
    
    return bd_results
        
def hidr_status_local(
        model: seekr2_base.Model,
        stage_processes: dict[str, multiprocessing.Process],
        work_dir: pathlib.Path
        ) -> dict:
    """
    Check the Hidr stage status locally. This is done by:
    1. Checking the process state file
    2. Checking if the process is still alive
    3. Reading the model files to see if all anchors have starting PDB files
    """
    # Check for state file
    state_file = find_latest_state_file(work_dir, "hidr")
    process_state = None
    
    if state_file and state_file.exists():
        try:
            process_state = LocalProcessState.load(state_file)
        except Exception as e:
            print(f"Warning: Could not load HIDR state file: {e}")
    
    # Check actual HIDR progress from model files
    filled_anchors = []
    unfilled_anchors = []
    for alpha, anchor in enumerate(model.anchors):
        has_structure = False
        if anchor.bulkstate:
            continue
        if anchor.amber_params is not None:
            if anchor.amber_params.pdb_coordinates_filename != "":
                full_path = os.path.join(
                    model.anchor_rootdir, anchor.directory, anchor.building_directory,
                    anchor.amber_params.pdb_coordinates_filename)
                if os.path.exists(
                        full_path):
                    has_structure = True
        elif anchor.charmm_params is not None:
            if anchor.charmm_params.pdb_coordinates_filename != "":
                full_path = os.path.join(
                    model.anchor_rootdir, anchor.directory, anchor.building_directory,
                    anchor.charmm_params.pdb_coordinates_filename)
                if os.path.exists(
                        full_path):
                    has_structure = True
        elif anchor.forcefield_params is not None:
            if anchor.forcefield_params.pdb_coordinates_filename != "":
                full_path = os.path.join(
                    model.anchor_rootdir, anchor.directory, anchor.building_directory,
                    anchor.forcefield_params.pdb_coordinates_filename)
                if os.path.exists(
                        full_path):
                    has_structure = True
        if has_structure:
            filled_anchors.append(alpha)
        else:
            unfilled_anchors.append(alpha)
    if len(filled_anchors) == 0:
        progress = "unstarted"
    elif len(unfilled_anchors) == 0:
        progress = "completed"
    else:
        progress = "started"
    
    stage_status = {
        "filled_anchors": sorted(filled_anchors),
        "unfilled_anchors": sorted(unfilled_anchors),
        "total_anchors": len(filled_anchors) + len(unfilled_anchors),
        "progress": progress,
        "progress_bar": len(filled_anchors) / (len(filled_anchors) \
            + len(unfilled_anchors)) if (len(filled_anchors) \
            + len(unfilled_anchors)) > 0 else 0.0,
        "model_xml_found": True,
        "notes": ""
    }
    
    # Build manager status based on process state (similar to BD)
    stage_manager_status = {}
    stage_manager_status["tool"] = "multiprocessing"
    stage_manager_status["timestamp"] = time.time()
    stage_manager_status["error"] = ""
    stage_manager_status["jobs"] = []
    
    # Check if process is running
    process = stage_processes.get("hidr")
    if process is not None and process.is_alive():
        job = {
            "JobID": f"local_pid_{process.pid}",
            "State": "RUNNING",
            "PID": process.pid
        }
        if process_state:
            job["process_status"] = process_state.status
            job["started_at"] = process_state.started_at
            job["output_file"] = process_state.progress_info.get("output_file", "N/A") if process_state.progress_info else "N/A"
        stage_manager_status["jobs"].append(job)
    elif process_state:
        if process_state.status == "running":
            if check_pid_exists(process_state.pid):
                job = {
                    "JobID": f"local_pid_{process_state.pid}",
                    "State": "RUNNING",
                    "PID": process_state.pid,
                    "process_status": process_state.status,
                    "started_at": process_state.started_at,
                    "output_file": process_state.progress_info.get("output_file", "N/A") if process_state.progress_info else "N/A"
                }
                stage_manager_status["jobs"].append(job)
            else:
                stage_status["notes"] = f"Process {process_state.pid} exited unexpectedly"
                stage_manager_status["error"] = "Process crashed"
        elif process_state.status in ["failed", "killed"]:
            stage_status["notes"] = f"Process {process_state.status}: {process_state.error}"
            if process_state.status == "failed":
                stage_manager_status["error"] = process_state.error
    
    hidr_results = {"success": True, "error": None, "manager_status": stage_manager_status,
                  "stage_status": stage_status}
    return hidr_results

def get_anchor_simulation_time(anchor_dir: str) -> float:
        """
        Get the current simulation time for an anchor by reading the 
        mmvt.restart#.dcd files.
        
        Returns the simulation time in the same units as stored in the file,
        or 0.0 if no restart files are found or can't be parsed.
        """
        import re
        import pathlib
        anchor_dir = pathlib.Path(anchor_dir)
        prod_dir = anchor_dir / "prod"
        if not prod_dir.exists():
            return 0.0
        
        # Find all mmvt.restart*.dcd files
        restart_files = list(prod_dir.glob("mmvt.restart*.out"))
        if not restart_files:
            return 0.0
        
        # Extract the numbers from filenames and find the highest
        def extract_restart_number(filepath):
            """Extract number from mmvt.restart#.out filename"""
            match = re.search(r'mmvt\.restart(\d+)\.out', filepath.name)
            if match:
                return int(match.group(1))
            return -1
        
        restart_files.sort(key=extract_restart_number)
        latest_restart = restart_files[-1]
        
        try:
            # Read the last line of the file
            with open(latest_restart, 'r') as f:
                # Read all lines and get the last one
                lines = f.readlines()
                if not lines:
                    return 0.0
                
                last_line = lines[-1].strip()
                if last_line.startswith("#"):
                    return 0.0
                # Parse the last line
                # Format: "1,1493,550.786" or "CHECKPOINT,20.000"
                parts = last_line.split(',')
                if len(parts) >= 2:
                    # The last part is the simulation time
                    time_str = parts[-1].strip()
                    return float(time_str)
                
        except Exception:
            # If we can't read or parse the file, return 0.0
            return 0.0
        
        return 0.0

def seekr_anchors_status_local(
        model: seekr2_base.Model,
        stage_processes: dict[str, dict[int, multiprocessing.Process]],
        work_dir: pathlib.Path,
        benchmark_mode: bool = False
        ) -> dict:
    """
    Check if SEEKR has finished locally.
    Checks for multiple per-anchor processes and state files.
    """
    # Check for state files for each anchor
    state_dir = get_local_state_dir(work_dir)
    anchor_process_states = {}
    
    # Find all seekr anchor state files
    seekr_state_files = list(state_dir.glob("seekr_anchor_*_state_*.json"))
    for state_file in seekr_state_files:
        try:
            state = LocalProcessState.load(state_file)
            anchor_index = state.progress_info.get("anchor_index") if state.progress_info else None
            if anchor_index is not None:
                # Keep track of state for each anchor (latest one if multiple)
                if anchor_index not in anchor_process_states or state.started_at > anchor_process_states[anchor_index].started_at:
                    anchor_process_states[anchor_index] = state
        except Exception as e:
            print(f"Warning: Could not load SEEKR state file {state_file}: {e}")
    
    # Check actual SEEKR progress from model files
    total_sim_time = 0.0
    total_time_for_completion = 0.0
    anchor_times = {}
    completed_anchors = []
    incomplete_anchors = []
    if benchmark_mode:
        time_for_completion = base.BENCHMARK_MIN_SIMULATION_LENGTH * model.get_timestep()
    else:
        time_for_completion = model.calculation_settings.num_production_steps \
            * model.get_timestep()
    for alpha, anchor in enumerate(model.anchors):
        if anchor.bulkstate:
            continue
        if not anchor.md:
            continue
        
        anchor_dir = os.path.join(model.anchor_rootdir, anchor.directory)
        sim_time = get_anchor_simulation_time(anchor_dir)
        anchor_times[alpha] = sim_time
        total_sim_time += sim_time
        total_time_for_completion += time_for_completion
        
        if sim_time >= time_for_completion:
            completed_anchors.append(alpha)
        elif sim_time > 0.0:
            incomplete_anchors.append(alpha)
        else:
            incomplete_anchors.append(alpha)

    # Determine SEEKR stage
    if len(incomplete_anchors) == 0 and len(completed_anchors) > 0:
        progress = "completed"
    elif total_sim_time > 0.0:
        progress = "started"
    else:
        progress = "unstarted"
    
    stage_status = {
        "completed_anchors": sorted(completed_anchors),
        "incomplete_anchors": sorted(incomplete_anchors),
        "anchor_times": anchor_times,
        "anchor_time_for_completion": time_for_completion,
        "total_anchors": len(completed_anchors) + len(incomplete_anchors),
        "progress": progress,
        "progress_bar": total_sim_time / total_time_for_completion if total_time_for_completion > 0 else 0.0,
        "model_xml_found": True,
        "notes": ""
    }
    
    # Build manager status - check all anchor processes
    stage_manager_status = {}
    stage_manager_status["tool"] = "multiprocessing"
    stage_manager_status["timestamp"] = time.time()
    stage_manager_status["error"] = ""
    stage_manager_status["jobs"] = []
    
    # Check processes from stage_processes dict (keyed by anchor index)
    seekr_processes = stage_processes.get("seekr", {})
    for anchor_index, process in seekr_processes.items():
        if process is not None and process.is_alive():
            job = {
                "JobID": f"local_pid_{process.pid}_anchor_{anchor_index}",
                "State": "RUNNING",
                "PID": process.pid,
                "anchor_index": anchor_index
            }
            # Add state info if available
            if anchor_index in anchor_process_states:
                state = anchor_process_states[anchor_index]
                job["process_status"] = state.status
                job["started_at"] = state.started_at
                job["output_file"] = state.progress_info.get("output_file", "N/A") if state.progress_info else "N/A"
            stage_manager_status["jobs"].append(job)
    
    # Check for running processes from state files (for reattached processes)
    for anchor_index, state in anchor_process_states.items():
        # Skip if we already added this from stage_processes
        if anchor_index in seekr_processes and seekr_processes[anchor_index] is not None:
            continue
            
        if state.status == "running":
            if check_pid_exists(state.pid):
                job = {
                    "JobID": f"local_pid_{state.pid}_anchor_{anchor_index}",
                    "State": "RUNNING",
                    "PID": state.pid,
                    "anchor_index": anchor_index,
                    "process_status": state.status,
                    "started_at": state.started_at,
                    "output_file": state.progress_info.get("output_file", "N/A") if state.progress_info else "N/A"
                }
                stage_manager_status["jobs"].append(job)
            else:
                stage_status["notes"] += f" Anchor {anchor_index} process {state.pid} exited unexpectedly."
        elif state.status in ["failed", "killed"]:
            stage_status["notes"] += f" Anchor {anchor_index} {state.status}: {state.error}."
            if state.status == "failed":
                stage_manager_status["error"] += f" Anchor {anchor_index} failed."
    
    seekr_results = {"success": True, "error": None, "manager_status": stage_manager_status,
                  "stage_status": stage_status}
    return seekr_results

# ============================================================================
# Local Execution Functions with State Management
# ============================================================================

def run_bd_locally(
        model_filename: str, 
        n_threads: int, 
        root_dir: pathlib.Path):
    """
    Run BD stage locally using multiprocessing.
    Writes state file with PID for monitoring and reattachment.
    Redirects stdout/stderr to bd_run.out.
    
    This function is executed in a separate process.
    """
    import os
    import sys
    import time
    import traceback
    import signal
    
    # Make this process a process group leader so all subprocesses (like nam_simulation) 
    # are in the same group and can be killed together
    os.setpgrp()
    
    pid = os.getpid()
    state_file = get_local_state_file(root_dir, "bd", pid)
    output_file = root_dir / "bd_run.out"
    
    # Create initial state BEFORE redirecting output (so any errors are visible)
    state = LocalProcessState(
        stage="bd",
        pid=pid,
        status="running",
        started_at=time.time(),
        ended_at=None,
        error=None,
        traceback=None,
        model_file=model_filename,
        work_dir=str(root_dir),
        progress_info={"n_threads": n_threads, "output_file": str(output_file)}
    )
    
    try:
        # Write initial state file
        state.save(state_file)
        
        # Create/update symlink to latest state file
        latest_link = root_dir / ".multiprocessing" / "bd_latest.json"
        if latest_link.exists() or latest_link.is_symlink():
            latest_link.unlink()
        latest_link.symlink_to(state_file.name)
        
    except Exception as e:
        # If we can't write state file, fail loudly BEFORE redirecting output
        print(f"FATAL: Cannot create state file {state_file}: {e}", file=sys.stderr)
        traceback.print_exc()
        raise
    
    # NOW redirect output at OS level using dup2 (so subprocesses via os.system() are also redirected)
    try:
        output_fd = os.open(str(output_file), os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
        os.dup2(output_fd, 1)  # Redirect stdout (fd 1)
        os.dup2(output_fd, 2)  # Redirect stderr (fd 2)
        os.close(output_fd)
        
        # Also redirect at Python level for consistency
        sys.stdout = open(output_file, 'a', buffering=1)  # Line buffered, append mode
        sys.stderr = sys.stdout
    except Exception as e:
        print(f"FATAL: Cannot redirect output to {output_file}: {e}", file=sys.stderr)
        raise
    
    print(f"BD stage started at {time.ctime()}")
    print(f"PID: {pid}")
    print(f"Model file: {model_filename}")
    print(f"Threads: {n_threads}")
    print(f"State file: {state_file}")
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
        import seekr2.modules.common_base as seekr2_base
        import seekr2.run as seekr2_run
        
        # Update state - imports successful
        state.progress_info["imports_successful"] = True
        state.save(state_file)
        print("Imports successful")
        sys.stdout.flush()
        
        # Load model
        model = seekr2_base.load_model(model_filename)
        state.progress_info["model_loaded"] = True
        state.save(state_file)
        print("Model loaded")
        sys.stdout.flush()
        
        # Run BD simulation
        print("Starting BD simulation...")
        sys.stdout.flush()
        seekr2_run.run(model, "any_bd", n_threads=n_threads)
        
        # Update state to completed (process itself updates this)
        state.status = "completed"
        state.ended_at = time.time()
        state.progress_info["bd_finished"] = True
        state.save(state_file)
        print(f"\nBD simulation completed at {time.ctime()}")
        sys.stdout.flush()
        
    except Exception as e:
        # Update state to failed (process itself updates this)
        state.status = "failed"
        state.error = str(e)
        state.traceback = traceback.format_exc()
        state.ended_at = time.time()
        state.save(state_file)
        print(f"\nBD simulation failed: {e}")
        traceback.print_exc()
        sys.stdout.flush()
        raise  # Re-raise so process exits with error code

def run_hidr_locally(
        model_filename: str, 
        input_pdb_file: str, 
        root_dir: pathlib.Path):
    """
    Run HIDR stage locally using multiprocessing.
    Writes state file with PID for monitoring and reattachment.
    Redirects stdout/stderr to hidr_run.out.
    
    This function is executed in a separate process.
    """
    import os
    import sys
    import time
    import traceback
    import signal
    
    # Make this process a process group leader
    os.setpgrp()
    
    pid = os.getpid()
    state_file = get_local_state_file(root_dir, "hidr", pid)
    output_file = root_dir / "hidr_run.out"
    
    # Create initial state BEFORE redirecting output
    state = LocalProcessState(
        stage="hidr",
        pid=pid,
        status="running",
        started_at=time.time(),
        ended_at=None,
        error=None,
        traceback=None,
        model_file=model_filename,
        work_dir=str(root_dir),
        progress_info={"input_pdb_file": input_pdb_file, "output_file": str(output_file)}
    )
    
    try:
        # Write initial state file
        state.save(state_file)
        
        # Create/update symlink to latest state file
        latest_link = root_dir / ".multiprocessing" / "hidr_latest.json"
        if latest_link.exists() or latest_link.is_symlink():
            latest_link.unlink()
        latest_link.symlink_to(state_file.name)
        
    except Exception as e:
        print(f"FATAL: Cannot create state file {state_file}: {e}", file=sys.stderr)
        traceback.print_exc()
        raise
    
    # Redirect output at OS level
    try:
        output_fd = os.open(str(output_file), os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
        os.dup2(output_fd, 1)
        os.dup2(output_fd, 2)
        os.close(output_fd)
        
        sys.stdout = open(output_file, 'a', buffering=1)
        sys.stderr = sys.stdout
    except Exception as e:
        print(f"FATAL: Cannot redirect output to {output_file}: {e}", file=sys.stderr)
        raise
    
    print(f"HIDR stage started at {time.ctime()}")
    print(f"PID: {pid}")
    print(f"Model file: {model_filename}")
    print(f"Input PDB: {input_pdb_file}")
    print(f"State file: {state_file}")
    print("-" * 60)
    sys.stdout.flush()
    
    # Handle termination signals
    def signal_handler(signum, frame):
        print(f"\nReceived signal {signum}, updating state to 'killed'...")
        sys.stdout.flush()
        state.status = "killed"
        state.ended_at = time.time()
        state.error = f"Process killed by signal {signum}"
        state.save(state_file)
        sys.exit(128 + signum)
    
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)
    
    try:
        # Import here so child process gets fresh imports
        import seekr2.modules.common_base as seekr2_base
        import seekrtools.hidr.hidr as seekr2_hidr
        
        state.progress_info["imports_successful"] = True
        state.save(state_file)
        print("Imports successful")
        sys.stdout.flush()
        
        # Load model
        model = seekr2_base.load_model(model_filename)
        state.progress_info["model_loaded"] = True
        state.save(state_file)
        print("Model loaded")
        sys.stdout.flush()
        
        # Change to anchor root directory
        os.chdir(os.path.abspath(model.anchor_rootdir))
        print(f"Changed to directory: {os.getcwd()}")
        sys.stdout.flush()
        
        # Run HIDR simulation
        print("Starting HIDR simulation...")
        sys.stdout.flush()
        seekr2_hidr.hidr(model, "any", pdb_files=[input_pdb_file], mode="MetaD", skip_checks=True)
        
        # Update state to completed
        state.status = "completed"
        state.ended_at = time.time()
        state.progress_info["hidr_finished"] = True
        state.save(state_file)
        print(f"\nHIDR simulation completed at {time.ctime()}")
        sys.stdout.flush()
        
    except Exception as e:
        # Update state to failed
        state.status = "failed"
        state.error = str(e)
        state.traceback = traceback.format_exc()
        state.ended_at = time.time()
        state.save(state_file)
        print(f"\nHIDR simulation failed: {e}")
        traceback.print_exc()
        sys.stdout.flush()
        raise

def run_seekr_locally(
        model_filename: str, 
        anchor_index: int, 
        root_dir: pathlib.Path,
        benchmark_mode: bool = False):
    """
    Run SEEKR stage for a single anchor locally using multiprocessing.
    Writes state file with PID and anchor index for monitoring and reattachment.
    Redirects stdout/stderr to seekr_<index>_run.out.
    
    This function is executed in a separate process.
    """
    import os
    import sys
    import time
    import traceback
    import signal
    
    # Make this process a process group leader
    os.setpgrp()
    
    pid = os.getpid()
    # Include anchor index in state filename for per-anchor tracking
    state_dir = get_local_state_dir(root_dir)
    state_file = state_dir / f"seekr_anchor_{anchor_index}_state_{pid}.json"
    output_file = root_dir / f"seekr_{anchor_index}_run.out"
    
    # Create initial state BEFORE redirecting output
    state = LocalProcessState(
        stage="seekr",
        pid=pid,
        status="running",
        started_at=time.time(),
        ended_at=None,
        error=None,
        traceback=None,
        model_file=model_filename,
        work_dir=str(root_dir),
        progress_info={"anchor_index": anchor_index, "output_file": str(output_file)}
    )
    
    try:
        # Write initial state file
        state.save(state_file)
        
        # Create/update symlink to latest state file for this anchor
        latest_link = root_dir / ".multiprocessing" / f"seekr_anchor_{anchor_index}_latest.json"
        if latest_link.exists() or latest_link.is_symlink():
            latest_link.unlink()
        latest_link.symlink_to(state_file.name)
        
    except Exception as e:
        print(f"FATAL: Cannot create state file {state_file}: {e}", file=sys.stderr)
        traceback.print_exc()
        raise
    
    # Redirect output at OS level
    try:
        output_fd = os.open(str(output_file), os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
        os.dup2(output_fd, 1)
        os.dup2(output_fd, 2)
        os.close(output_fd)
        
        sys.stdout = open(output_file, 'a', buffering=1)
        sys.stderr = sys.stdout
    except Exception as e:
        print(f"FATAL: Cannot redirect output to {output_file}: {e}", file=sys.stderr)
        raise
    
    print(f"SEEKR stage started at {time.ctime()}")
    print(f"PID: {pid}")
    print(f"Model file: {model_filename}")
    print(f"Anchor index: {anchor_index}")
    print(f"State file: {state_file}")
    print("-" * 60)
    sys.stdout.flush()
    
    # Handle termination signals
    def signal_handler(signum, frame):
        print(f"\nReceived signal {signum}, updating state to 'killed'...")
        sys.stdout.flush()
        state.status = "killed"
        state.ended_at = time.time()
        state.error = f"Process killed by signal {signum}"
        state.save(state_file)
        sys.exit(128 + signum)
    
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)
    
    try:
        # Import here so child process gets fresh imports
        import seekr2.modules.common_base as seekr2_base
        import seekr2.run as seekr2_run
        
        state.progress_info["imports_successful"] = True
        state.save(state_file)
        print("Imports successful")
        sys.stdout.flush()
        
        # Load model
        model = seekr2_base.load_model(model_filename)
        state.progress_info["model_loaded"] = True
        state.save(state_file)
        print("Model loaded")
        sys.stdout.flush()
        
        # Run SEEKR simulation for this anchor
        print(f"Starting SEEKR simulation for anchor {anchor_index}...")
        sys.stdout.flush()
        if benchmark_mode:
            seekr2_run.run(model, f"{anchor_index}", 
                           min_total_simulation_length=base.BENCHMARK_MIN_SIMULATION_LENGTH)
        else:
            seekr2_run.run(model, f"{anchor_index}")
        
        # Update state to completed
        state.status = "completed"
        state.ended_at = time.time()
        state.progress_info["seekr_finished"] = True
        state.save(state_file)
        print(f"\nSEEKR simulation for anchor {anchor_index} completed at {time.ctime()}")
        sys.stdout.flush()
        
    except Exception as e:
        # Update state to failed
        state.status = "failed"
        state.error = str(e)
        state.traceback = traceback.format_exc()
        state.ended_at = time.time()
        state.save(state_file)
        print(f"\nSEEKR simulation for anchor {anchor_index} failed: {e}")
        traceback.print_exc()
        sys.stdout.flush()
        raise

def force_bd_rerun(model: seekr2_base.Model) -> None:
    """
    Force a re-run of the BD stage by deleting the requisite files.
    """
    bd_results_files = glob.glob(
        os.path.join(model.anchor_rootdir,  model.k_on_info.b_surface_directory, 
                     model.k_on_info.bd_output_glob))
    for filename in bd_results_files:
        os.remove(filename)
    traj_files = glob.glob(
        os.path.join(model.anchor_rootdir,  model.k_on_info.b_surface_directory, 
                     "traj*"))
    for filename in traj_files:
        os.remove(filename)
    simulation_xml_files = glob.glob(
        os.path.join(model.anchor_rootdir,  model.k_on_info.b_surface_directory, 
                     "*_simulation.xml"))
    for filename in simulation_xml_files:
        os.remove(filename)
    return

def force_hidr_rerun(model: seekr2_base.Model) -> None:
    """
    Force a re-run of the HIDR stage by deleting the requisite files.
    """
    os.chdir(os.path.abspath(model.anchor_rootdir))
    original_model_filename = "model_pre_hidr_0.xml"
    model_filename = "model.xml"
    if os.path.exists(original_model_filename):
        os.remove(model_filename)
        os.rename(original_model_filename, model_filename)
        other_pre_hidr_models = glob.glob("model_pre_hidr_*.xml")
        for other_model in other_pre_hidr_models:
            os.remove(other_model)
    return

def force_seekr_rerun(model: seekr2_base.Model) -> None:
    """
    Force a re-run of the SEEKR stage by deleting the requisite files.
    """
    for anchor in model.anchors:
        if anchor.bulkstate:
            continue
        prod_directory_all_files = glob.glob(
            os.path.join(model.anchor_rootdir, anchor.directory, 
                         anchor.production_directory, "*"))
        for filename in prod_directory_all_files:
            os.remove(filename)
    return

def force_rerun_stage(
        model: seekr2_base.Model,
        stage_name: str,
        root_directory_path: pathlib.Path,
        ) -> None:
    """
    Force a re-run of a given stage by deleting requisite files.
    Supported stages: "bd", "hidr", "seekr"
    """
    # First, kill any running process for this stage
    # TODO: for seekr, find_latest_state_file only returns one anchor's state file; remaining
    #  anchor processes are not killed. Fix by iterating all seekr_anchor_*_state_*.json files.
    state_file = find_latest_state_file(root_directory_path, stage_name)
    if state_file and state_file.exists():
        try:
            state = LocalProcessState.load(state_file)
            if state.status == "running" and check_pid_exists(state.pid):
                pid = state.pid
                print(f"  Killing existing {stage_name} process (PID: {pid})...")
                try:
                    os.killpg(pid, signal.SIGTERM)
                    time.sleep(2)
                    
                    # Force kill if still alive
                    if check_pid_exists(pid):
                        os.killpg(pid, signal.SIGKILL)
                        time.sleep(1)
                    
                    # Update state file
                    state.status = "killed"
                    state.ended_at = time.time()
                    state.error = "Process killed for force re-run"
                    state.save(state_file)
                    print(f"  Killed {stage_name} process (PID: {pid})")
                except (ProcessLookupError, PermissionError) as e:
                    print(f"  Warning: Could not kill process: {e}")
        except Exception as e:
            print(f"  Warning: Could not check for existing process: {e}")
    
    # Then delete the stage files
    if stage_name == "bd":
        force_bd_rerun(model)
    elif stage_name == "hidr":
        force_hidr_rerun(model)
        # Must reload model.
        model_filename = str(root_directory_path / "model.xml")
        model = seekr2_base.load_model(model_filename)
    elif stage_name == "seekr":
        force_seekr_rerun(model)
    return

# Setup cleanup function for multiprocessing
def cleanup_local_processes(
        stage_processes: dict[str, multiprocessing.Process | dict[int, multiprocessing.Process]],
        root_directory_path: pathlib.Path
        ) -> None:
    """
    Clean up all running local processes and their subprocesses using process groups.
    Note: The child process updates its own state file when it receives SIGTERM.
    Parent only updates state if force kill (SIGKILL) was needed.
    """
    print("\nCleaning up local processes...")
    for stage_name, process_or_dict in stage_processes.items():
        # Handle BD and HIDR (single process)
        if stage_name in ["bd", "hidr"]:
            process = process_or_dict
            if process is not None and process.is_alive():
                pid = process.pid
                print(f"Terminating {stage_name} process group (PID: {pid})...")
                try:
                    os.killpg(pid, signal.SIGTERM)
                except ProcessLookupError:
                    print(f"Process group {pid} already terminated")
                    continue
                
                process.join(timeout=5)
                
                if process.is_alive():
                    print(f"Force killing {stage_name} process group (PID: {pid})...")
                    try:
                        os.killpg(pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    process.join()
                    
                    # Update state file since child couldn't
                    try:
                        state_file = find_latest_state_file(root_directory_path, stage_name)
                        if state_file and state_file.exists():
                            state = LocalProcessState.load(state_file)
                            if state.status == "running":
                                state.status = "killed"
                                state.ended_at = time.time()
                                state.error = "Process force-killed (SIGKILL) by parent"
                                state.save(state_file)
                                print(f"Updated state file for {stage_name} to 'killed'")
                    except Exception as e:
                        print(f"Warning: Could not update state file for {stage_name}: {e}")
        
        # Handle SEEKR (dict of processes by anchor index)
        elif stage_name == "seekr":
            anchor_processes = process_or_dict
            for anchor_index, process in anchor_processes.items():
                if process is not None and process.is_alive():
                    pid = process.pid
                    print(f"Terminating seekr anchor {anchor_index} process group (PID: {pid})...")
                    try:
                        os.killpg(pid, signal.SIGTERM)
                    except ProcessLookupError:
                        print(f"Process group {pid} already terminated")
                        continue
                    
                    process.join(timeout=5)
                    
                    if process.is_alive():
                        print(f"Force killing seekr anchor {anchor_index} process group (PID: {pid})...")
                        try:
                            os.killpg(pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                        process.join()
                        
                        # Update state file since child couldn't
                        try:
                            state_dir = get_local_state_dir(root_directory_path)
                            state_file = state_dir / f"seekr_anchor_{anchor_index}_state_{pid}.json"
                            if state_file.exists():
                                state = LocalProcessState.load(state_file)
                                if state.status == "running":
                                    state.status = "killed"
                                    state.ended_at = time.time()
                                    state.error = "Process force-killed (SIGKILL) by parent"
                                    state.save(state_file)
                                    print(f"Updated state file for seekr anchor {anchor_index} to 'killed'")
                        except Exception as e:
                            print(f"Warning: Could not update state file for seekr anchor {anchor_index}: {e}")
    print("Process cleanup complete.")

