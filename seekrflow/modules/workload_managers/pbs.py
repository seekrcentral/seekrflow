"""
PBS/Torque workload manager functions for remote job submission and status checking.

This module provides PBS/Torque equivalents to the SLURM workload manager functions,
following the same patterns and conventions.
"""

def calculate_optimal_seekr_time_limit(
    job_status_file: str,
    incomplete_anchors: list[int],
    default_time_limit: str
) -> str:
    """
    Calculate optimal PBS walltime for SEEKR stage based on remaining work
    and benchmark data from previous runs.

    This function estimates how much wallclock time is needed to complete the
    remaining simulation work for incomplete anchors by:
    1. Reading anchor simulation times from .seekrflow_job_status.json
    2. Calculating remaining simulation time per anchor
    3. Estimating performance from previous job's elapsed time and progress
    4. Computing time needed using the slowest anchor
    5. Adding 20% safety buffer and rounding up to nearest hour
    6. Capping at default_time_limit

    Args:
        job_status_file: Path to .seekrflow_job_status.json file
        incomplete_anchors: List of anchor indices that need more simulation
        default_time_limit: Fallback time limit (HH:MM:SS format) if calculation fails

    Returns:
        Time limit string in HH:MM:SS format, capped at default_time_limit
        Returns default_time_limit if:
        - No benchmark data available (first run)
        - Previous run made no progress
        - Any calculation errors occur
        - Calculated time exceeds default_time_limit
    """
    import os
    import json
    import math
    import pathlib

    # Fallback to default if file doesn't exist
    if not os.path.exists(job_status_file):
        print("DEBUG: Job status file not found, using default time limit.")
        return default_time_limit

    try:
        # Load status file
        with open(job_status_file, "r") as f:
            status_data = json.load(f)

        # Extract SEEKR stage data
        seekr_data = status_data.get("seekr", {})
        stage_status = seekr_data.get("stage_status", {})
        manager_status = seekr_data.get("manager_status")

        # Get anchor times and target
        anchor_times = stage_status.get("anchor_times", {})
        target_time = stage_status.get("anchor_time_for_completion", 0.0)

        if target_time <= 0 or not anchor_times:
            print("DEBUG: Invalid target time or empty anchor times, using default time limit.")
            return default_time_limit

        # Get elapsed time - prioritize last_known_elapsed from status result
        elapsed_str = seekr_data.get("last_known_elapsed")

        # Fallback to manager_status if available
        if not elapsed_str and manager_status:
            jobs = manager_status.get("jobs", [])
            if jobs:
                elapsed_str = jobs[0].get("resources_used", {}).get("walltime")

        if not elapsed_str:
            print("DEBUG: No elapsed time found (no previous job data), using default time limit.")
            return default_time_limit

        def parse_pbs_time(time_str: str) -> float:
            """Convert PBS time string (HH:MM:SS) to seconds"""
            parts = time_str.split(":")
            if len(parts) == 3:
                hours, minutes, seconds = map(int, parts)
                return hours * 3600 + minutes * 60 + seconds
            return 0.0

        elapsed_seconds = parse_pbs_time(elapsed_str)
        print(f"DEBUG: Previous job elapsed time: {elapsed_seconds} seconds.")

        if elapsed_seconds <= 0:
            print("DEBUG: Elapsed time is zero, using default time limit.")
            return default_time_limit

        # Read anchor_times_at_submission from local state file for baseline
        job_status_path = pathlib.Path(job_status_file)
        root_dir = job_status_path.parent
        pbs_runner_dir = root_dir / ".pbs_runner"
        seekr_latest = pbs_runner_dir / "seekr_latest.json"

        anchor_times_at_submission = None
        if seekr_latest.exists():
            try:
                with open(seekr_latest, "r") as f:
                    state_data = json.load(f)
                    anchor_times_at_submission = state_data.get("anchor_times_at_submission")
            except Exception:
                pass

        # Calculate simulation progress made in the previous job
        if anchor_times_at_submission is not None:
            # We have baseline - calculate actual progress
            total_progress = 0.0
            for anchor_idx in incomplete_anchors:
                anchor_str = str(anchor_idx)
                current_time = anchor_times.get(anchor_str, 0.0)
                previous_time = anchor_times_at_submission.get(anchor_str, 0.0)
                progress = max(0.0, current_time - previous_time)
                total_progress += progress
        else:
            # No baseline - estimate using current simulation times for incomplete anchors
            total_progress = sum(
                anchor_times.get(str(idx), 0.0) for idx in incomplete_anchors
            )

        # If no progress was made, fall back to default
        print(f"DEBUG: Total progress made in previous job: {total_progress} ps.")
        if total_progress <= 0:
            print("DEBUG: No progress made in previous job, using default time limit.")
            return default_time_limit

        # Calculate benchmark: seconds per picosecond of simulation
        benchmark = elapsed_seconds / total_progress
        print(f"DEBUG: Benchmark calculated: {benchmark} seconds per unit simulation time.")

        # Calculate time needed for each incomplete anchor
        max_time_needed = 0.0
        for anchor_idx in incomplete_anchors:
            anchor_str = str(anchor_idx)
            current_time = anchor_times.get(anchor_str, 0.0)
            remaining_time = max(0.0, target_time - current_time)
            time_needed = remaining_time * benchmark
            max_time_needed = max(max_time_needed, time_needed)

        print(f"DEBUG: Maximum time needed among incomplete anchors: {max_time_needed} seconds.")
        if max_time_needed <= 0:
            print("DEBUG: No time needed for remaining anchors, using default time limit.")
            return default_time_limit

        # Add 20% safety buffer
        safe_time = max_time_needed * 1.20
        print(f"DEBUG: Safe time with buffer: {safe_time} seconds.")

        # Round up to nearest hour
        hours_needed = math.ceil(safe_time / 3600.0)
        print(f"DEBUG: Hours needed (rounded up): {hours_needed} hours.")

        # Enforce minimum of 1 hour
        hours_needed = max(1, hours_needed)
        print(f"DEBUG: Hours needed after enforcing minimum: {hours_needed} hours.")

        # Format as HH:MM:SS
        calculated_time = f"{hours_needed:02d}:00:00"

        # Parse default_time_limit to compare
        default_seconds = parse_pbs_time(default_time_limit)
        calculated_seconds = hours_needed * 3600

        # Cap at default_time_limit
        print("DEBUG: Calculated time limit:", calculated_time)
        print("DEBUG: Default time limit:", default_time_limit)
        if calculated_seconds > default_seconds:
            return default_time_limit

        return calculated_time

    except Exception as e:
        # If anything goes wrong, fall back to default
        print(f"Warning: Could not calculate optimal time limit: {e}")
        return default_time_limit

def pbs_remote_bd_status_workflow(args):
    """
    Check BD simulation status remotely by examining PBS queue and
    output files directly.

    Returns consistent dict structure:
    {
        "success": bool,
        "error": str or None,
        "manager_status": {...},  # if SLURM job found
        "stage_status": {          # if success
            "bd_finished": bool,
            "milestones_remaining": int,
            "results_files_found": int,
            "details": {...}
        }
    }
    """

    import time
    import json
    import pathlib
    import subprocess
    from dataclasses import dataclass
    import xml.etree.ElementTree as ET
    from typing import Dict, List, Optional

    working_dir = args[0]
    stage_name = "bd"
    n_traj_for_completed = args[1]

    def run(cmd: List[str], check: bool = True) -> tuple:
        """Run command, return (rc, stdout, stderr)"""
        out = subprocess.run(cmd, stdout=subprocess.PIPE, 
                           stderr=subprocess.PIPE, text=True)
        if check and out.returncode != 0:
            raise RuntimeError(f"Command failed: ({out.returncode}): "
                             f"{' '.join(cmd)}\nSTDERR:\n{out.stderr}")
        return out.returncode, out.stdout.strip(), out.stderr.strip()
    
    @dataclass
    class RunState:
        """State tracking for PBS job runs"""
        run_id: str
        workdir: str
        logdir: str
        name: str
        queue: str
        scheduler_options: Optional[str]
        time_limit: str
        account: str
        resource_list: Optional[str]
        cpus_per_task: int
        mem: Optional[str]
        array_spec: Optional[list[int]]
        jobid: str
        submitted_at: float
        n_tasks: int
        cmd_template: str
        attempts: Dict[str, List[int]]
        parent_jobids: List[str]
        worker_init: Optional[str] = None
        state_file: Optional[str] = None
        anchor_times_at_submission: Optional[Dict[str, float]] = None
        last_known_elapsed: Optional[str] = None

        @staticmethod
        def load(path: pathlib.Path) -> "RunState":
            with open(path, "r") as f:
                data = json.load(f)
            return RunState(**data)
        
    def get_state_dir(workdir: pathlib.Path) -> pathlib.Path:
        return workdir / ".pbs_runner"

    def check_bd_completion(
            model_path: pathlib.Path,
            n_traj_target: int
            ) -> Dict:
        """
        Check BD completion by reading model.xml and BD results files directly.
        """
        try:
            # Get the directory of the model_path
            unstarted_dict = {
                "finished": False,
                "progress": "unstarted",
                "n_trajectories": 0,
                "n_trajectories_for_completion": n_traj_target,
                "progress_bar": 0.0,
                "model_xml_found": False,
                "notes": ""
            }
            model_dir = model_path.parent
            if not model_path.exists():
                unstarted_dict["notes"] = f"Model file not found: {model_path}"
                return unstarted_dict
            b_surface_dir = model_dir / "b_surface"
            if not b_surface_dir.exists():
                unstarted_dict["notes"] = f"b_surface directory not found: {b_surface_dir}"
                unstarted_dict["model_xml_found"] = True
                return unstarted_dict
            results_files = list(b_surface_dir.glob("results*.xml"))
            if not results_files:
                unstarted_dict["notes"] = f"No results files found in: {b_surface_dir}"
                unstarted_dict["model_xml_found"] = True
                return unstarted_dict

            # Sort by the number in the filename and get the largest
            def extract_number(filepath):
                """Extract number from results*.xml filename"""
                stem = filepath.stem  # e.g., "results1" or "results"
                if stem == "results":
                    return 0
                try:
                    return int(stem.replace("results", ""))
                except ValueError:
                    return 0
            
            results_files.sort(key=extract_number)
            #latest_results = results_files[-1]
            n_trajectories = 0
            for results_file in reversed(results_files):
                
                # Parse the latest results file
                tree = ET.parse(results_file)
                root = tree.getroot()
            
                # Find the <n_trajectories> tag
                n_traj_elem = root.find(".//n_trajectories")

                if n_traj_elem is None:
                    return {
                        "finished": False,
                        "progress": "started",
                        "n_trajectories": 0,
                        "n_trajectories_for_completion": n_traj_target,
                        "progress_bar": 0.0,
                        "model_xml_found": True,
                        "notes": "n_trajectories tag not found in results file",
                        }
                
                n_trajectories += int(n_traj_elem.text.strip())
            
            # Determine BD stage
            if n_trajectories < n_traj_target:
                progress = "started"
                finished = False
            else:
                progress = "completed"
                finished = True
            
            return {
                "finished": finished,
                "progress": progress,
                "n_trajectories": n_trajectories,
                "n_trajectories_for_completion": n_traj_target,
                "progress_bar": n_trajectories / n_traj_target if n_traj_target > 0 else 0.0,
                "model_xml_found": True,
                "notes": ""
                }
            
        except Exception as e:
            return {
                "notes": f"Unexpected error checking BD completion: {e}",
                "finished": False,
                "progress": "unknown",
                "n_trajectories": 0,
                "n_trajectories_for_completion": n_traj_target,
                "progress_bar": 0.0,
                "model_xml_found": False,
            }
        
    # Initialize result structure
    result = {
        "success": False,
        "error": None,
        "manager_status": None,
        "stage_status": None
    }

    work_dir = pathlib.Path(working_dir)
    state_path = get_state_dir(work_dir) / f"{stage_name}_latest.json"
    
    # Check PBS status if state file exists
    if state_path.exists():
        try:
            st = RunState.load(state_path)
            rc, out, err = run(["bash", "-lc", f"qstat -f -F json {st.jobid}"], check=False)
            
            if rc == 0:
                try:
                    qstat_data = json.loads(out)
                    jobs_dict = qstat_data.get("Jobs", {})
                    
                    if st.jobid in jobs_dict:
                        job_info = jobs_dict[st.jobid]
                        result["manager_status"] = {
                            "tool": "qstat",
                            "timestamp": time.time(),
                            "error": "",
                            "jobs": [job_info]
                        }
                    else:
                        # Job not found in qstat output
                        result["manager_status"] = {
                            "tool": "qstat",
                            "timestamp": time.time(),
                            "error": f"Job {st.jobid} not found in qstat output",
                            "jobs": []
                        }
                except json.JSONDecodeError as e:
                    result["manager_status"] = {
                        "tool": "qstat",
                        "timestamp": time.time(),
                        "error": f"Failed to parse qstat JSON output: {e}",
                        "jobs": []
                    }
            else:
                result["manager_status"] = {
                    "tool": "qstat",
                    "timestamp": time.time(),
                    "error": err,
                    "jobs": []
                }
        except Exception as e:
            result["error"] = f"Failed to check PBS status: {e}"
    
    # Check BD completion status by reading files
    model_path = work_dir / "model.xml"
    stage_status = check_bd_completion(model_path, n_traj_for_completed)
    result["stage_status"] = stage_status
    result["success"] = True
    
    return result


def pbs_remote_hidr_status_workflow(args):
    """
    Check the status of HIDR simulation on PBS by examining PBS queue 
    and checking which anchors have structures assigned in model.xml.
    
    This function is executed remotely via Globus Compute. It checks:
    1. PBS job queue status for running/pending HIDR jobs
    2. model.xml file to determine which anchors have been filled with structures
    
    Parameters
    ----------
    args : list
        List containing:
        - work_dir : str
            Path to the working directory containing model.xml
    
    Returns
    -------
    dict
        Status information containing:
        - success : bool
            Whether the status check completed successfully
        - error : str or None
            Error message if any
        - manager_status : dict or None
            PBS job information if state file exists
        - stage_status : dict
            Dictionary with:
            - filled_anchors : list of int
                Indices of anchors that have structures assigned
            - unfilled_anchors : list of int
                Indices of anchors that don't have structures yet
            - total_anchors : int
                Total number of HIDR anchors
            - progress : str
                One of "unstarted", "started", or "completed"
    """
    import time
    import json
    import pathlib
    import subprocess
    import xml.etree.ElementTree as ET
    from dataclasses import dataclass
    from typing import Dict, List, Optional

    working_dir = args[0]
    stage_name = "hidr"

    def run(cmd: List[str], check: bool = True) -> tuple:
        """Run command, return (rc, stdout, stderr)"""
        out = subprocess.run(cmd, stdout=subprocess.PIPE, 
                           stderr=subprocess.PIPE, text=True)
        if check and out.returncode != 0:
            raise RuntimeError(f"Command failed: ({out.returncode}): "
                             f"{' '.join(cmd)}\nSTDERR:\n{out.stderr}")
        return out.returncode, out.stdout.strip(), out.stderr.strip()
    
    @dataclass
    class RunState:
        """State tracking for PBS job runs"""
        run_id: str
        workdir: str
        logdir: str
        name: str
        queue: str
        scheduler_options: Optional[str]
        time_limit: str
        account: str
        resource_list: Optional[str]
        cpus_per_task: int
        mem: Optional[str]
        array_spec: Optional[list[int]]
        jobid: str
        submitted_at: float
        n_tasks: int
        cmd_template: str
        attempts: Dict[str, List[int]]
        parent_jobids: List[str]
        worker_init: Optional[str] = None
        state_file: Optional[str] = None
        anchor_times_at_submission: Optional[Dict[str, float]] = None
        last_known_elapsed: Optional[str] = None

        @staticmethod
        def load(path: pathlib.Path) -> "RunState":
            with open(path, "r") as f:
                data = json.load(f)
            return RunState(**data)
        
    def get_state_dir(workdir: pathlib.Path) -> pathlib.Path:
        return workdir / ".pbs_runner"
    
    def check_hidr_completion(model_path: pathlib.Path) -> Dict:
        """
        Check HIDR completion by reading model.xml to see which anchors
        have been filled with structures.
        """
        try:
            model_dir = model_path.parent
            if not model_path.exists():
                return {
                    "filled_anchors": [],
                    "unfilled_anchors": [],
                    "total_anchors": 0,
                    "progress": "unstarted",
                    "progress_bar": 0.0,
                    "model_xml_found": False,
                    "notes": f"Model file not found: {model_path}"
                }
            
            tree = ET.parse(str(model_path))
            root = tree.getroot()
            
            # Find all HIDR anchors
            filled_anchors = []
            unfilled_anchors = []
            
            # Find the <anchors> list element
            anchors_list = root.find(".//anchors")
            if anchors_list is None:
                return {
                    "filled_anchors": [],
                    "unfilled_anchors": [],
                    "total_anchors": 0,
                    "progress": "unstarted",
                    "progress_bar": 0.0,
                    "model_xml_found": True,
                    "notes": "No anchors list found in model.xml"
                }
            
            # Iterate through all child elements of <anchors>
            # They will be named like anchors_e0, anchors_e1, etc.
            for anchor_elem in anchors_list:
                # Check if this element has an <index> child
                index_elem = anchor_elem.find("index")
                bulk_elem = anchor_elem.find("bulkstate")
                if index_elem is not None and index_elem.text:
                    anchor_idx = int(index_elem.text.strip())
                    bulk_anchor = True if bulk_elem.text.strip() == "True" else False
                    if bulk_anchor:
                        continue
                    # Check all parameter types for pdb_coordinates_filename
                    has_structure = False
                    for params_type in ["amber_params", "charmm_params", 
                                       "forcefield_params"]:
                        params = anchor_elem.find(params_type)
                        if params is not None:
                            pdb_file = params.find("pdb_coordinates_filename")
                            if pdb_file is not None and pdb_file.text:
                                pdb_text = pdb_file.text.strip()
                                if pdb_text:  # Not empty string
                                    # Check if the file actually exists
                                    anchor_dir_elem = anchor_elem.find("directory")
                                    anchor_build_elem = anchor_elem.find("building_directory")
                                    if anchor_dir_elem is not None and anchor_dir_elem.text:
                                        anchor_dir = model_dir / anchor_dir_elem.text.strip() \
                                            / anchor_build_elem.text.strip()
                                        pdb_path = anchor_dir / pdb_text
                                    else:
                                        pdb_path = model_dir / pdb_text
                                    
                                    if pdb_path.exists():
                                        has_structure = True
                                        break
                    
                    if has_structure:
                        filled_anchors.append(anchor_idx)
                    else:
                        unfilled_anchors.append(anchor_idx)
            
            # Determine HIDR stage
            if len(filled_anchors) == 0:
                progress = "unstarted"
            elif len(unfilled_anchors) == 0:
                progress = "completed"
            else:
                progress = "started"
            
            return {
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
            
        except Exception as e:
            return {
                "filled_anchors": [],
                "unfilled_anchors": [],
                "total_anchors": 0,
                "progress": "unknown",
                "progress_bar": 0.0,
                "model_xml_found": False,
                "notes": f"Unexpected error checking HIDR completion: {e}"
            }
        
    # Initialize result structure
    result = {
        "success": False,
        "error": None,
        "manager_status": None,
        "stage_status": None
    }

    work_dir = pathlib.Path(working_dir)
    state_path = get_state_dir(work_dir) / f"{stage_name}_latest.json"
    
    # Check PBS status if state file exists
    if state_path.exists():
        try:
            st = RunState.load(state_path)
            rc, out, err = run(["bash", "-lc", f"qstat -f -F json {st.jobid}"], check=False)
            
            if rc == 0:
                try:
                    qstat_data = json.loads(out)
                    jobs_dict = qstat_data.get("Jobs", {})
                    
                    if st.jobid in jobs_dict:
                        job_info = jobs_dict[st.jobid]
                        result["manager_status"] = {
                            "tool": "qstat",
                            "timestamp": time.time(),
                            "error": "",
                            "jobs": [job_info]
                        }
                    else:
                        # Job not found in qstat output
                        result["manager_status"] = {
                            "tool": "qstat",
                            "timestamp": time.time(),
                            "error": f"Job {st.jobid} not found in qstat output",
                            "jobs": []
                        }
                except json.JSONDecodeError as e:
                    result["manager_status"] = {
                        "tool": "qstat",
                        "timestamp": time.time(),
                        "error": f"Failed to parse qstat JSON output: {e}",
                        "jobs": []
                    }
            else:
                result["manager_status"] = {
                    "tool": "qstat",
                    "timestamp": time.time(),
                    "error": err,
                    "jobs": []
                }
        except Exception as e:
            result["error"] = f"Failed to check PBS status: {e}"
    
    # Check HIDR completion status by reading files
    model_path = work_dir / "model.xml"
    stage_status = check_hidr_completion(model_path)
    result["stage_status"] = stage_status
    result["success"] = True
    
    return result


def pbs_remote_seekr_status_workflow(args):
    """
    Check the status of SEEKR simulation on PBS by examining PBS queue
    and checking which anchors have MD simulations completed in model.xml.
    
    This function is executed remotely via Globus Compute. It checks:
    1. PBS job queue status for running/pending SEEKR jobs
    2. model.xml and output directories to determine simulation progress
    
    Parameters
    ----------
    args : list
        List containing:
        - work_dir : str
            Path to the working directory containing model.xml
        - anchor_time_for_completion : float
            The simulation time (in ps or ns) required for an anchor to be 
            considered complete
    
    Returns
    -------
    dict
        Status information containing:
        - success : bool
            Whether the status check completed successfully
        - error : str or None
            Error message if any
        - manager_status : dict or None
            PBS job information if state file exists
        - stage_status : dict
            Dictionary with:
            - completed_anchors : list of int
                Indices of anchors with completed simulations
            - incomplete_anchors : list of int
                Indices of anchors needing more simulation
            - total_anchors : int
                Total number of SEEKR anchors
            - progress : str
                One of "unstarted", "started", or "completed"
        - last_known_elapsed : str or None
            Elapsed time (HH:MM:SS) for time optimization
    """
    import time
    import json
    import pathlib
    import subprocess
    import xml.etree.ElementTree as ET
    from dataclasses import dataclass
    from typing import Dict, List, Optional
    import re

    working_dir = args[0]
    stage_name = "seekr"
    anchor_time_for_completion = args[1]

    def run(cmd: List[str], check: bool = True) -> tuple:
        """Run command, return (rc, stdout, stderr)"""
        out = subprocess.run(cmd, stdout=subprocess.PIPE, 
                           stderr=subprocess.PIPE, text=True)
        if check and out.returncode != 0:
            raise RuntimeError(f"Command failed: ({out.returncode}): "
                             f"{' '.join(cmd)}\nSTDERR:\n{out.stderr}")
        return out.returncode, out.stdout.strip(), out.stderr.strip()

    @dataclass
    class RunState:
        """State tracking for PBS job runs"""
        run_id: str
        workdir: str
        logdir: str
        name: str
        queue: str
        scheduler_options: Optional[str]
        time_limit: str
        account: str
        resource_list: Optional[str]
        cpus_per_task: int
        mem: Optional[str]
        array_spec: Optional[list[int]]
        jobid: str
        submitted_at: float
        n_tasks: int
        cmd_template: str
        attempts: Dict[str, List[int]]
        parent_jobids: List[str]
        worker_init: Optional[str] = None
        state_file: Optional[str] = None
        anchor_times_at_submission: Optional[Dict[str, float]] = None
        last_known_elapsed: Optional[str] = None

        @staticmethod
        def load(path: pathlib.Path) -> "RunState":
            with open(path, "r") as f:
                data = json.load(f)
            return RunState(**data)
        
    def get_state_dir(workdir: pathlib.Path) -> pathlib.Path:
        return workdir / ".pbs_runner"
    
    def get_anchor_simulation_time(anchor_dir: pathlib.Path) -> float:
        """
        Get the current simulation time for an anchor by reading the 
        mmvt.restart#.out files.
        
        Returns the simulation time in the same units as stored in the file,
        or 0.0 if no restart files are found or can't be parsed.
        """
        prod_dir = anchor_dir / "prod"
        if not prod_dir.exists():
            return 0.0
        
        # Find all mmvt.restart*.out files
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
    
    def check_seekr_completion(
            model_path: pathlib.Path,
            time_for_completion: float
            ) -> Dict:
        """
        Check SEEKR completion by reading model.xml and checking for
        simulation time in mmvt.restart#.out files.
        """
        try:
            model_dir = model_path.parent
            if not model_path.exists():
                return {
                    "completed_anchors": [],
                    "incomplete_anchors": [],
                    "anchor_times": {},
                    "anchor_time_for_completion": 0.0,
                    "total_anchors": 0,
                    "progress": "unstarted",
                    "progress_bar": 0.0,
                    "model_xml_found": False,
                    "notes": f"Model file not found: {model_path}"
                }
            
            tree = ET.parse(str(model_path))
            root = tree.getroot()
            
            # Find all SEEKR anchors (milestones with MD attribute)
            completed_anchors = []
            incomplete_anchors = []
            
            anchors_list = root.find(".//anchors")
            if anchors_list is None:
                return {
                    "completed_anchors": [],
                    "incomplete_anchors": [],
                    "anchor_times": {},
                    "anchor_time_for_completion": 0.0,
                    "total_anchors": 0,
                    "progress": "unstarted",
                    "progress_bar": 0.0,
                    "model_xml_found": True,
                    "notes": "No anchors list found in model.xml"
                }
            
            # Iterate through all child elements of <anchors>
            # They will be named like anchors_e0, anchors_e1, etc.
            total_sim_time = 0.0
            total_time_for_completion = 0.0
            anchor_times = {}
            for anchor_elem in anchors_list:
                # Check if this element has an <index> child
                index_elem = anchor_elem.find("index")
                bulk_elem = anchor_elem.find("bulkstate")
                if index_elem is not None and index_elem.text:
                    anchor_idx = int(index_elem.text.strip())
                    bulk_anchor = bulk_elem is not None and bulk_elem.text.strip() == "True"
                    if bulk_anchor:
                        continue
                    
                    # Check if this anchor needs MD simulations
                    # The <md> tag is a direct child of the anchor element
                    needs_md = False
                    md_elem = anchor_elem.find("md")
                    if md_elem is not None and md_elem.text:
                        md_attr = md_elem.text.strip()
                        if md_attr.lower() == "true":
                            needs_md = True
                    
                    if not needs_md:
                        # This anchor doesn't need MD, skip it
                        continue
                    
                    # Check simulation time in anchor directory
                    anchor_dir = model_dir / f"anchor_{anchor_idx}"
                    sim_time = get_anchor_simulation_time(anchor_dir)
                    anchor_times[anchor_idx] = sim_time
                    total_sim_time += sim_time
                    total_time_for_completion += time_for_completion
                    
                    if sim_time >= time_for_completion:
                        completed_anchors.append(anchor_idx)
                    elif sim_time > 0.0:
                        # Has some simulation time but not complete
                        incomplete_anchors.append(anchor_idx)
                    else:
                        # No simulation time found - unstarted
                        incomplete_anchors.append(anchor_idx)
            
            # Determine SEEKR stage
            if len(incomplete_anchors) == 0 and len(completed_anchors) > 0:
                progress = "completed"
            elif total_sim_time > 0.0:
                progress = "started"
            else:
                progress = "unstarted"
            
            return {
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
            
        except Exception as e:
            return {
                "completed_anchors": [],
                "incomplete_anchors": [],
                "anchor_times": {},
                "anchor_time_for_completion": 0.0,
                "total_anchors": 0,
                "progress": "unknown",
                "progress_bar": 0.0,
                "model_xml_found": False,
                "notes": f"Unexpected error checking SEEKR completion: {e}"
            }
        
    # Initialize result structure
    result = {
        "success": False,
        "error": None,
        "manager_status": None,
        "stage_status": None
    }

    work_dir = pathlib.Path(working_dir)
    state_path = get_state_dir(work_dir) / f"{stage_name}_latest.json"

    # Track elapsed time for optimization (persists after job finishes)
    last_known_elapsed = None

    # Check PBS status if state file exists
    if state_path.exists():
        try:
            st = RunState.load(state_path)
            rc, out, err = run(["bash", "-lc", f"qstat -f -F json {st.jobid}"], check=False)
            
            if rc == 0:
                try:
                    qstat_data = json.loads(out)
                    jobs_dict = qstat_data.get("Jobs", {})
                    
                    if st.jobid in jobs_dict:
                        job_info = jobs_dict[st.jobid]
                        
                        # Capture elapsed time before it disappears
                        resources_used = job_info.get("resources_used", {})
                        last_known_elapsed = resources_used.get("walltime")
                        
                        # Also save to state file for remote-side tracking
                        if last_known_elapsed:
                            try:
                                state_data = json.loads(state_path.read_text())
                                state_data["last_known_elapsed"] = last_known_elapsed
                                state_path.write_text(json.dumps(state_data, indent=2))
                            except Exception:
                                pass
                        
                        result["manager_status"] = {
                            "tool": "qstat",
                            "timestamp": time.time(),
                            "error": "",
                            "jobs": [job_info]
                        }
                    else:
                        # Job not found in qstat output - try to get elapsed from state file
                        try:
                            state_data = json.loads(state_path.read_text())
                            last_known_elapsed = state_data.get("last_known_elapsed")
                        except Exception:
                            pass
                        
                        result["manager_status"] = {
                            "tool": "qstat",
                            "timestamp": time.time(),
                            "error": f"Job {st.jobid} not found in qstat output",
                            "jobs": []
                        }
                except json.JSONDecodeError as e:
                    # Failed to parse JSON - try to get elapsed from state file
                    try:
                        state_data = json.loads(state_path.read_text())
                        last_known_elapsed = state_data.get("last_known_elapsed")
                    except Exception:
                        pass
                    
                    result["manager_status"] = {
                        "tool": "qstat",
                        "timestamp": time.time(),
                        "error": f"Failed to parse qstat JSON output: {e}",
                        "jobs": []
                    }
            else:
                # No jobs in queue - try to get elapsed from state file
                try:
                    state_data = json.loads(state_path.read_text())
                    last_known_elapsed = state_data.get("last_known_elapsed")
                except Exception:
                    pass
                
                result["manager_status"] = {
                    "tool": "qstat",
                    "timestamp": time.time(),
                    "error": err,
                    "jobs": []
                }
        except Exception as e:
            result["error"] = f"Failed to check PBS status: {e}"

    # Include elapsed time in result so it's transferred back to local machine
    result["last_known_elapsed"] = last_known_elapsed

    # Check SEEKR completion status by reading files
    model_path = work_dir / "model.xml"
    stage_status = check_seekr_completion(model_path, anchor_time_for_completion)
    result["stage_status"] = stage_status
    result["success"] = True
    
    return result


def pbs_remote_run_workflow(args):
    """
    Main PBS workflow execution function.
    Submits jobs to PBS queue and monitors execution.

    Parameters
    ----------
    args : list
        List containing:
        - working_dir : str (args[0])
        - stage_name : str (args[1])
        - logdir : str (args[2])
        - queue : str (args[3])
        - account : str (args[4])
        - resource_list : str (args[5])
        - name : str (args[6])
        - cpus_per_task : int (args[7])
        - mem_per_node : int (args[8])
        - time_limit : str (args[9])
        - scheduler_options : str (args[10])
        - worker_init : str (args[11])
        - command_string : str (args[12])
        - indices : list[int] or None (args[13])
        - model_filename : str (args[14])
        - workflow_type : str (args[15])
        - mps : int (args[16])
        - anchor_times_at_submission : dict or None (args[17])

    Returns
    -------
    dict
        {"success": bool, "error": str or None, "jobid": str, "run_id": str}
    """
    import os
    import json
    import time
    import uuid
    import shlex
    import pathlib
    import subprocess
    from dataclasses import dataclass, asdict
    from typing import Dict, List, Optional

    working_dir = args[0]
    stage_name = args[1]
    logdir = args[2]
    queue = args[3]
    account = args[4]
    resource_list = args[5]
    name = args[6]
    cpus_per_task = args[7]
    mem_per_node = args[8]
    time_limit = args[9]
    scheduler_options = args[10]
    worker_init = args[11]
    command_string = args[12]
    indices = args[13]
    model_filename = args[14]
    workflow_type = args[15]
    mps = args[16]
    anchor_times_at_submission = args[17] if len(args) > 17 else None

    kwargs = {
        "working_dir": working_dir,
        "logdir": logdir,
        "queue": queue,
        "account": account,
        "resource_list": resource_list,
        "name": name,
        "cpus_per_task": cpus_per_task,
        "mem": mem_per_node,
        "time_limit": time_limit,
        "scheduler_options": scheduler_options,
        "worker_init": worker_init,
        "command_string": command_string,
        "indices": indices,
        "model_filename": model_filename,
        "workflow_type": workflow_type,
        "mps": mps,
        "anchor_times_at_submission": anchor_times_at_submission
    }

    @dataclass
    class RunState:
        run_id: str
        workdir: str
        logdir: str
        name: str
        queue: str
        scheduler_options: Optional[str]
        time_limit: str
        account: str
        resource_list: Optional[str]
        cpus_per_task: Optional[int]
        mem: Optional[str]
        array_spec: Optional[list[int]]
        jobid: str
        submitted_at: float
        n_tasks: int
        cmd_template: str
        attempts: Dict[str, List[int]]
        parent_jobids: List[str]
        worker_init: Optional[str] = None
        state_file: Optional[str] = None
        anchor_times_at_submission: Optional[Dict[str, float]] = None
        last_known_elapsed: Optional[str] = None

        def save(self, path: pathlib.Path) -> None:
            ensure_dir(path.parent)
            self.state_file = str(path)
            with open(path, "w") as f:
                json.dump(asdict(self), f, indent=2)

        @staticmethod
        def load(path: pathlib.Path) -> "RunState":
            with open(path, "r") as f:
                data = json.load(f)
            return RunState(**data)

    def run(cmd: List[str], check: bool = True) -> tuple:
        """Run command, return (rc, stdout, stderr)"""
        out = subprocess.run(cmd, stdout=subprocess.PIPE, 
                           stderr=subprocess.PIPE, text=True)
        if check and out.returncode != 0:
            raise RuntimeError(f"Command failed: ({out.returncode}): "
                             f"{' '.join(cmd)}\nSTDERR:\n{out.stderr}")
        return out.returncode, out.stdout.strip(), out.stderr.strip()

    def get_state_dir(workdir: pathlib.Path) -> pathlib.Path:
        return workdir / ".pbs_runner"

    def ensure_dir(p: pathlib.Path) -> None:
        p.mkdir(parents=True, exist_ok=True)

    work_dir = pathlib.Path(working_dir)
    state_path = get_state_dir(work_dir) / f"{stage_name}_latest.json"

    def collapse_indices(idxs: List[int]) -> str:
        """
        Collapse [0,1,2,5,6,9] -> '0-2,5-6,9'
        """
        if not idxs:
            return ""
        idxs = sorted(set(idxs))
        ranges = []
        start = prev = idxs[0]
        for x in idxs[1:]:
            if x == prev + 1:
                prev = x
            else:
                ranges.append(f"{start}-{prev}" if start != prev else str(start))
                start = prev = x
        ranges.append(f"{start}-{prev}" if start != prev else str(start))
        return ",".join(ranges)

    def submit_array(
            name: str,
            queue: str,
            scheduler_options: Optional[str],
            time_limit: str,
            account: Optional[str],
            resource_list: Optional[str],
            array_spec: list[int] | None,
            logdir: pathlib.Path,
            workdir: pathlib.Path,
            cpus_per_task: Optional[int],
            mem: Optional[str],
            cmd_template: str,
            worker_init: str = "",
            mps: int = 1
            ) -> str:
        ensure_dir(logdir)

        # Create new condensed array spec for MPS bundling
        if array_spec is not None and mps > 1:
            n_bundles = (len(array_spec) + mps - 1) // mps  # ceiling division
            condensed_array_spec = list(range(n_bundles))
        else:
            condensed_array_spec = array_spec

        # Build PBS script
        script_lines = ["#!/bin/bash"]
        script_lines.append(f"#PBS -N {name}")
        script_lines.append("#PBS -V")  # Export all environment variables
        if queue:
            script_lines.append(f"#PBS -q {queue}")
        if account:
            script_lines.append(f"#PBS -A {account}")
        script_lines.append(f"#PBS -l walltime={time_limit}")
        script_lines.append(f"#PBS -l nodes=1:ppn={cpus_per_task or 1}")
        if mem:
            script_lines.append(f"#PBS -l mem={mem}")
        if resource_list:
            script_lines.append(f"#PBS -l {resource_list}")

        if condensed_array_spec is None:
            script_lines.append(f"#PBS -o {logdir}/{name}.out")
            script_lines.append(f"#PBS -e {logdir}/{name}.err")
        else:
            array_str = collapse_indices(condensed_array_spec)
            script_lines.append(f"#PBS -J {array_str}")
            # For array jobs, specify output directory only
            # PBS will use default naming: jobname.oJOBID-ARRAYID
            script_lines.append(f"#PBS -o {logdir}/")
            script_lines.append(f"#PBS -e {logdir}/")

        if scheduler_options:
            for line in scheduler_options.split("\n"):
                if line.strip():
                    script_lines.append(line.strip())

        script_lines.append("")
        script_lines.append(f"echo PBS_ARRAY_INDEX: $PBS_ARRAY_INDEX; cd {workdir}")

        # Replace placeholders in cmd_template
        wrap_cmd = cmd_template

        if array_spec is not None and mps > 1:
            # For MPS mode: replace {index0}, {index1}, etc. and {len_array}
            for mps_index in range(mps):
                # Calculate the original index: PBS_ARRAY_INDEX * mps + mps_index
                original_index = f"$(( ${{PBS_ARRAY_INDEX:error}} * {mps} + {mps_index} ))"
                wrap_cmd = wrap_cmd.replace(f"{{index{mps_index}}}", original_index)

            # Calculate len_array: how many indices this PBS task should process
            # Formula: min(mps, total_indices - task_id * mps)
            remaining = f"{len(array_spec)} - ${{PBS_ARRAY_INDEX:error}} * {mps}"
            len_array_formula = f"$(( ({remaining}) < {mps} ? ({remaining}) : {mps} ))"
            wrap_cmd = wrap_cmd.replace("{len_array}", len_array_formula)
        else:
            # Single index mode: just replace {index} with PBS_ARRAY_INDEX
            wrap_cmd = wrap_cmd.replace("{index}", "${PBS_ARRAY_INDEX:error}")

        if worker_init:
            script_lines.append("")
            script_lines.append("# Worker initialization")
            script_lines.append(worker_init)

        script_lines.append("")
        script_lines.append("# Execute command")
        script_lines.append(wrap_cmd)

        script = "\n".join(script_lines)
        
        # Submit job via qsub
        result = subprocess.run(
            ["qsub"],
            input=script,
            capture_output=True,
            text=True,
            timeout=30
        )

        if result.returncode != 0:
            raise RuntimeError(f"qsub failed: {result.stderr}")

        # PBS returns different formats:
        # Non-array job: "12345.panamint4"
        # Array job: "6[].panamint4"
        # For qstat queries, we need just the base numeric job ID
        jobid = result.stdout.strip()
        
        # Remove [] if present (array jobs)
        jobid = jobid.replace('[]', '')
        
        # Remove .hostname suffix to get just numeric ID
        if '.' in jobid:
            jobid = jobid.split('.')[0]
        
        return jobid

    def cmd_submit(args):
        workdir = pathlib.Path(args["working_dir"] or os.getcwd()).resolve()
        logdir = (workdir / (args["logdir"] or "logs")).resolve()
        ensure_dir(workdir)
        ensure_dir(logdir)
        array_spec = args["indices"]
        if array_spec is None:
            n_tasks = 1
        else:
            n_tasks = len(array_spec)
        run_id = uuid.uuid4().hex[:10]
        state_dir = get_state_dir(workdir)
        ensure_dir(state_dir)
        jobid = submit_array(
            name=args["name"],
            queue=args["queue"],
            scheduler_options=args["scheduler_options"],
            time_limit=args["time_limit"],
            account=args["account"],
            resource_list=args["resource_list"],
            array_spec=array_spec,
            logdir=logdir,
            workdir=workdir,
            cpus_per_task=args["cpus_per_task"],
            mem=args["mem"],
            cmd_template=args["command_string"],
            worker_init=args["worker_init"],
            mps=args["mps"]
        )
        st = RunState(
            run_id=run_id,
            workdir=str(workdir),
            logdir=str(logdir),
            name=args["name"],
            queue=args["queue"],
            scheduler_options=args["scheduler_options"],
            time_limit=args["time_limit"],
            account=args["account"],
            resource_list=args["resource_list"],
            cpus_per_task=args["cpus_per_task"],
            mem=args["mem"],
            array_spec=array_spec,
            jobid=jobid,
            submitted_at=time.time(),
            n_tasks=n_tasks,
            cmd_template=args["command_string"],
            attempts={jobid: list(range(n_tasks))},
            parent_jobids=[],
            worker_init=args["worker_init"],
            anchor_times_at_submission=args.get("anchor_times_at_submission")
        )
        state_path = state_dir / f"{stage_name}_{run_id}.json"
        st.save(state_path)
        # symlink latest
        latest = state_dir / f"{stage_name}_latest.json"
        try:
            if latest.exists() or latest.is_symlink():
                latest.unlink()
            latest.symlink_to(state_path.name)
        except Exception:
            (state_dir / f"{stage_name}_latest.json").write_text(
                json.dumps(asdict(st), indent=2))
        return {"success": True, "error": None, "jobid": jobid, "run_id": run_id}

    results = cmd_submit(kwargs)
    return results


def pbs_remote_cancel_workflow(args):
    """
    Cancel PBS job.

    Args:
        args[0]: working_dir
        args[1]: job_id
    """
    import subprocess
    
    working_dir = args[0]
    job_id = args[1]

    try:
        result = subprocess.run(
            ["qdel", job_id],
            capture_output=True,
            text=True,
            timeout=10
        )

        if result.returncode != 0 and "Unknown Job Id" not in result.stderr:
            return {
                "success": False,
                "error": f"qdel failed: {result.stderr}",
                "result": None
            }

        return {
            "success": True,
            "error": None,
            "result": {"canceled_job": job_id}
        }

    except Exception as e:
        return {
            "success": False,
            "error": str(e),
            "result": None
        }


def pbs_remote_force_rerun_workflow(args):
    """
    Force re-run a stage on a PBS system by:
    1. Canceling any running jobs for that stage
    2. Deleting the requisite files to reset the stage
    
    Parameters
    ----------
    args : list
        List containing:
        - working_dir : str
            Path to the working directory containing model.xml
        - stage_name : str
            Name of the stage ("bd", "hidr", or "seekr")
    
    Returns
    -------
    dict
        Result with success status and any errors
    """
    import subprocess
    import os
    import glob
    import pathlib
    import shlex
    import json
    import time
    import uuid
    import xml.etree.ElementTree as ET
    from dataclasses import dataclass
    from typing import Dict, List, Optional
    
    working_dir = args[0]
    stage_name = args[1]
    
    def run(cmd: List[str], check: bool = True) -> tuple:
        """Run command, return (rc, stdout, stderr)"""
        out = subprocess.run(cmd, stdout=subprocess.PIPE, 
                           stderr=subprocess.PIPE, text=True)
        if check and out.returncode != 0:
            raise RuntimeError(f"Command failed: ({out.returncode}): "
                             f"{' '.join(cmd)}\nSTDERR:\n{out.stderr}")
        return out.returncode, out.stdout.strip(), out.stderr.strip()
    
    @dataclass
    class RunState:
        run_id: str
        workdir: str
        logdir: str
        name: str
        queue: str
        scheduler_options: Optional[str]
        time_limit: str
        account: str
        resource_list: Optional[str]
        cpus_per_task: Optional[int]
        mem: Optional[str]
        array_spec: Optional[list[int]]
        jobid: str
        submitted_at: float
        n_tasks: int
        cmd_template: str
        attempts: Dict[str, List[int]]
        parent_jobids: List[str]
        worker_init: Optional[str] = None
        state_file: Optional[str] = None
        anchor_times_at_submission: Optional[Dict[str, float]] = None
        last_known_elapsed: Optional[str] = None

        @staticmethod
        def load(path: pathlib.Path) -> "RunState":
            with open(path, "r") as f:
                data = json.load(f)
            return RunState(**data)
    
    def get_state_dir(workdir: pathlib.Path) -> pathlib.Path:
        return workdir / ".pbs_runner"
    
    def force_bd_rerun_remote(model_dir: pathlib.Path) -> None:
        """Delete BD result files."""
        try:
            model_path = model_dir / "model.xml"
            if not model_path.exists():
                print(f"Model file not found: {model_path}")
                return
            
            tree = ET.parse(str(model_path))
            root = tree.getroot()
            
            # Find b_surface directory from model.xml
            b_surface_elem = root.find(".//b_surface_directory")
            if b_surface_elem is not None and b_surface_elem.text:
                b_surface_dir = model_dir / b_surface_elem.text.strip()
                if b_surface_dir.exists():
                    print(f"Deleting BD results directory: {b_surface_dir}")
                    # Delete all results files
                    for results_file in b_surface_dir.glob("results*.xml"):
                        results_file.unlink()
                        print(f"  Deleted {results_file.name}")
                    # Try to remove directory if empty
                    try:
                        b_surface_dir.rmdir()
                        print(f"  Removed empty directory {b_surface_dir}")
                    except OSError:
                        print(f"  Directory {b_surface_dir} not empty, keeping it")
                else:
                    print(f"BD directory does not exist: {b_surface_dir}")
            else:
                print("No b_surface_directory element found in model.xml")
                    
        except Exception as e:
            print(f"Error in force_bd_rerun_remote: {e}")
    
    def force_hidr_rerun_remote(model_dir: pathlib.Path) -> None:
        """Restore model.xml from pre-HIDR backup."""
        try:
            model_path = model_dir / "model.xml"
            backup_path = model_dir / "model_backup_pre_hidr.xml"
            
            if backup_path.exists():
                print(f"Restoring model.xml from {backup_path}")
                backup_path.replace(model_path)
                print("  Model restored successfully")
            else:
                print(f"No HIDR backup found at {backup_path}, nothing to restore")
                    
        except Exception as e:
            print(f"Error in force_hidr_rerun_remote: {e}")
    
    def force_seekr_rerun_remote(model_dir: pathlib.Path) -> None:
        """Delete production directory contents."""
        try:
            model_path = model_dir / "model.xml"
            if not model_path.exists():
                print(f"Model file not found: {model_path}")
                return
            
            tree = ET.parse(str(model_path))
            root = tree.getroot()
            
            # Find all SEEKR anchors
            anchors_list = root.find(".//anchors")
            if anchors_list is None:
                print("No anchors list found in model.xml")
                return
            
            deleted_count = 0
            for anchor_elem in anchors_list:
                directory_elem = anchor_elem.find("directory")
                if directory_elem is not None and directory_elem.text:
                    anchor_dir = model_dir / directory_elem.text.strip()
                    prod_dir = anchor_dir / "prod"
                    
                    if prod_dir.exists():
                        # Delete all files in prod directory
                        for item in prod_dir.iterdir():
                            if item.is_file():
                                item.unlink()
                                deleted_count += 1
                        print(f"  Deleted {deleted_count} files from {prod_dir}")
            
            if deleted_count > 0:
                print(f"Total files deleted from SEEKR production directories: {deleted_count}")
            else:
                print("No SEEKR production files found to delete")
                                
        except Exception as e:
            print(f"Error in force_seekr_rerun_remote: {e}")
    
    result = {
        "success": False,
        "error": None,
        "canceled_job": None
    }
    
    work_dir = pathlib.Path(working_dir)
    state_path = get_state_dir(work_dir) / f"{stage_name}_latest.json"
    
    # Step 1: Cancel any running job for this stage
    if state_path.exists():
        try:
            st = RunState.load(state_path)
            run(["bash", "-lc", f"qdel {shlex.quote(st.jobid)}"], check=False)
            result["canceled_job"] = st.jobid
            print(f"Canceled job {st.jobid} for {stage_name} stage")
        except Exception as e:
            print(f"Warning: Could not cancel job: {e}")
    
    # Step 2: Delete stage files
    try:
        if stage_name == "bd":
            force_bd_rerun_remote(work_dir)
        elif stage_name == "hidr":
            force_hidr_rerun_remote(work_dir)
        elif stage_name == "seekr":
            force_seekr_rerun_remote(work_dir)
        else:
            result["error"] = f"Unknown stage name: {stage_name}"
            return result
        
        print(f"Deleted files for {stage_name} stage")
        result["success"] = True
        
    except Exception as e:
        result["error"] = f"Failed to delete files: {e}"
        result["success"] = False
    
    return result
