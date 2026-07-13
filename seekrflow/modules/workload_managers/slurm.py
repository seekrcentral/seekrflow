"""
SLURM workload manager functions for remote job submission and status checking.
"""

import math
import pathlib
from typing import Optional


def calculate_optimal_seekr_time_limit(
    job_status_file: str,
    incomplete_anchors: list[int],
    default_time_limit: str
) -> str:
    """
    Calculate optimal SLURM time limit for SEEKR stage based on remaining work
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
                elapsed_str = jobs[0].get("Elapsed")

        if not elapsed_str:
            print("DEBUG: No elapsed time found (no previous job data), using default time limit.")
            return default_time_limit

        def parse_slurm_time(time_str: str) -> float:
            """Convert SLURM time string to seconds"""
            days = 0
            if "-" in time_str:
                day_part, time_part = time_str.split("-")
                days = int(day_part)
                time_str = time_part

            parts = time_str.split(":")
            if len(parts) == 3:
                hours, minutes, seconds = map(int, parts)
                return days * 86400 + hours * 3600 + minutes * 60 + seconds
            return 0.0

        elapsed_seconds = parse_slurm_time(elapsed_str)
        print(f"DEBUG: Previous job elapsed time: {elapsed_seconds} seconds.")

        if elapsed_seconds <= 0:
            print("DEBUG: Elapsed time is zero, using default time limit.")
            return default_time_limit

        # Read anchor_times_at_submission from local state file for baseline
        job_status_path = pathlib.Path(job_status_file)
        root_dir = job_status_path.parent
        slurm_runner_dir = root_dir / ".slurm_runner"
        seekr_latest = slurm_runner_dir / "seekr_latest.json"

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
            # This assumes these are the times accumulated in the last job
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
        # If no time needed (all anchors complete), shouldn't happen but handle it
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
        default_seconds = parse_slurm_time(default_time_limit)
        calculated_seconds = hours_needed * 3600

        # Cap at default_time_limit
        print("DEBUG: Calculated time limit:", calculated_time)
        print("DEBUG: Default time limit:", default_time_limit)
        print("DEBUG: Calculated seconds:", calculated_seconds)
        if calculated_seconds > default_seconds:
            return default_time_limit

        return calculated_time

    except Exception as e:
        # If anything goes wrong, fall back to default
        print(f"Warning: Could not calculate optimal time limit: {e}")
        return default_time_limit


def slurm_remote_status_workflow(args):
    """
    Generic SLURM status workflow.

    Args format:
        [working_dir, stage_name, benchmark_mode=False, anchor="any",
         swarm_id=None, worker_init=""]

    worker_init is a bash snippet that activates the env which has the
    `seekr` package installed (e.g. "source ~/.bashrc; conda activate SEEKR2").
    The seekr progress check is executed via subprocess under this snippet
    rather than imported in-process, because the Globus Compute endpoint
    worker typically runs in its own isolated Python that does not have
    seekr available.
    """
    import json
    import time
    import shlex
    import pathlib
    import subprocess
    from dataclasses import dataclass
    from typing import Dict, List, Optional

    if len(args) < 2:
        return {
            "success": False,
            "error": "Expected args: [working_dir, stage_name, ...]",
            "manager_status": None,
            "stage_status": {
                "finished": False,
                "state": "failed",
                "progress": 0.0,
                "notes": "invalid status workflow arguments",
            },
        }

    working_dir = args[0]
    stage_name = args[1]
    benchmark_mode = bool(args[2]) if len(args) > 2 else False
    anchor = args[3] if len(args) > 3 else "any"
    swarm_id = args[4] if len(args) > 4 else None
    worker_init = args[5] if len(args) > 5 else ""

    @dataclass
    class RunState:
        run_id: str
        workdir: str
        logdir: str
        name: str
        partition: str
        scheduler_options: Optional[str]
        time_limit: str
        account: str
        constraint: Optional[str]
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
        return workdir / ".slurm_runner"

    def find_latest_state_file(workdir: pathlib.Path, stage: str) -> Optional[pathlib.Path]:
        state_dir = get_state_dir(workdir)
        if not state_dir.exists():
            return None
        latest = state_dir / f"{stage}_latest.json"
        if latest.exists():
            return latest
        files = sorted(
            state_dir.glob(f"{stage}_*.json"),
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )
        if files:
            return files[0]
        return None

    def run(cmd: List[str], check: bool = True) -> tuple:
        out = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if check and out.returncode != 0:
            raise RuntimeError(
                f"Command failed: ({out.returncode}): {' '.join(cmd)}\n"
                f"STDERR:\n{out.stderr}"
            )
        return out.returncode, out.stdout.strip(), out.stderr.strip()

    work_dir = pathlib.Path(working_dir)
    manager_status = {
        "tool": "squeue",
        "timestamp": time.time(),
        "error": "",
        "jobs": [],
    }
    last_known_elapsed = None

    state_path = find_latest_state_file(work_dir, stage_name)
    if state_path is not None and state_path.exists():
        try:
            st = RunState.load(state_path)
            cmd = f"squeue -j {st.jobid} -h -o '%i|%P|%j|%u|%T|%M|%l|%D|%C|%R'"
            rc, out, err = run(["bash", "-lc", cmd], check=False)

            if rc == 0:
                rows = []
                for line in out.splitlines():
                    parts = line.split("|")
                    if len(parts) >= 9:
                        rows.append(
                            {
                                "JobID": parts[0],
                                "Partition": parts[1],
                                "Name": parts[2],
                                "User": parts[3],
                                "State": parts[4],
                                "Elapsed": parts[5],
                                "TimeLimit": parts[6],
                                "Nodes": parts[7],
                                "CPUs": parts[8],
                                "Reason": parts[9] if len(parts) > 9 else "",
                            }
                        )

                if rows:
                    last_known_elapsed = rows[0].get("Elapsed")
                    try:
                        state_data = json.loads(state_path.read_text())
                        state_data["last_known_elapsed"] = last_known_elapsed
                        state_path.write_text(json.dumps(state_data, indent=2))
                    except Exception:
                        pass

                manager_status["jobs"] = rows
            else:
                manager_status["error"] = err
                try:
                    state_data = json.loads(state_path.read_text())
                    last_known_elapsed = state_data.get("last_known_elapsed")
                except Exception:
                    pass
        except Exception as e:
            manager_status["error"] = f"Failed to check SLURM status: {e}"

    def _progress_from_stage_progress(stage_progress: dict, partitioned_anchor) -> float:
        progress_map = stage_progress.get("progress", {})
        if partitioned_anchor == "any":
            values = []
            if isinstance(progress_map, dict):
                for _, per_partition in progress_map.items():
                    if isinstance(per_partition, dict):
                        value = per_partition.get("progress")
                        if isinstance(value, (int, float)):
                            values.append(float(value))
            if values:
                return sum(values) / len(values)
            return 0.0

        partitioned_status = {}
        if isinstance(progress_map, dict):
            partitioned_status = progress_map.get(partitioned_anchor, {})
        if isinstance(partitioned_status, dict):
            value = partitioned_status.get("progress", 0.0)
            if isinstance(value, (int, float)):
                return float(value)
        return 0.0

    stage_status = {
        "finished": False,
        "state": "unknown",
        "progress": 0.0,
        "notes": "",
    }

    model_filename = work_dir / "model.json"
    if not model_filename.exists():
        stage_status.update(
            {
                "state": "unstarted",
                "notes": f"Model file not found: {model_filename}",
                "model_xml_found": False,
            }
        )
        return {
            "success": True,
            "error": None,
            "manager_status": manager_status,
            "stage_status": stage_status,
            "last_known_elapsed": last_known_elapsed,
        }

    # Run the seekr progress check via subprocess under worker_init.  The
    # endpoint worker's own Python may not have seekr installed (e.g. when
    # globus-compute-endpoint is in a pipx-isolated venv).  worker_init is
    # the same bash snippet used to activate the env for SLURM jobs.
    snippet_lines = [
        "import json, sys",
        "try:",
        "    import seekr.modules.structures as S",
        "    import seekr.status as st",
        f"    model = S.load_model({str(model_filename)!r})",
        "    stage_names = [s.name for s in model.stages]",
        f"    stage_name = {stage_name!r}",
        "    if stage_name not in stage_names:",
        "        print('__SEEKR_STATUS__' + json.dumps({",
        "            'ok': False, 'unknown_stage': True,",
        "            'stage_names': stage_names}))",
        "        sys.exit(0)",
        "    idx = stage_names.index(stage_name)",
        f"    msg = st.status(model, 'progress', stage_arg=stage_name, anchor_arg={anchor!r}, swarm_id={swarm_id!r}, print_json=True)",
        "    prog = msg.get('progress', {})",
        "    sp = {}",
        "    # seekr keys progress by stage.index (an int), which may differ",
        "    # from this stage's position in model.stages.  Match by name to be",
        "    # robust against either an int or stringified-int key.",
        "    for _k, _v in prog.items():",
        "        if isinstance(_v, dict) and _v.get('stage_name') == stage_name:",
        "            sp = _v",
        "            break",
        "    # Per-anchor time data in picoseconds for the time-limit optimizer.",
        "    try: _ts = float(model.get_timestep())",
        "    except Exception: _ts = 0.0",
        "    _per_anchor = sp.get('progress') if isinstance(sp.get('progress'), dict) else {}",
        "    _anchor_times = {}",
        "    _anchor_targets = []",
        "    _incomplete = []",
        "    for _ak, _av in _per_anchor.items():",
        "        if not isinstance(_av, dict): continue",
        "        _cur = _av.get('current_step')",
        "        _tot = _av.get('total_steps') or _av.get('maximum_steps')",
        "        if _cur is not None and _ts > 0.0:",
        "            _anchor_times[str(_ak)] = float(_cur) * _ts",
        "        if _tot is not None and _ts > 0.0:",
        "            _anchor_targets.append(float(_tot) * _ts)",
        "        if not _av.get('attained', False):",
        "            try: _incomplete.append(int(_ak))",
        "            except Exception: pass",
        "    sp['anchor_times'] = _anchor_times",
        "    sp['anchor_time_for_completion'] = max(_anchor_targets) if _anchor_targets else 0.0",
        "    sp['incomplete_anchors'] = sorted(_incomplete)",
        "    print('__SEEKR_STATUS__' + json.dumps({",
        "        'ok': True, 'stage_progress': sp}))",
        "except Exception as e:",
        "    import traceback",
        "    print('__SEEKR_STATUS__' + json.dumps({",
        "        'ok': False, 'error': str(e), 'traceback': traceback.format_exc()}))",
    ]
    snippet = "\n".join(snippet_lines)
    init = (worker_init or "").strip()
    if init:
        cmd = f"{init}; python -c {shlex.quote(snippet)}"
    else:
        cmd = f"python -c {shlex.quote(snippet)}"

    try:
        proc = subprocess.run(
            ["bash", "-lc", cmd],
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            timeout=120,
        )
        marker = "__SEEKR_STATUS__"
        payload = None
        for line in proc.stdout.splitlines():
            line = line.strip()
            if line.startswith(marker):
                try:
                    payload = json.loads(line[len(marker):])
                except Exception:
                    payload = None
        if payload is None:
            raise RuntimeError(
                f"could not parse seekr status output (rc={proc.returncode}); "
                f"stdout_tail={proc.stdout[-400:]!r}; "
                f"stderr_tail={proc.stderr[-400:]!r}"
            )
        if not payload.get("ok"):
            if payload.get("unknown_stage"):
                stage_status.update(
                    {
                        "state": "failed",
                        "notes": (
                            f"Unknown stage name in model: {stage_name}; "
                            f"available={payload.get('stage_names')}"
                        ),
                    }
                )
                return {
                    "success": False,
                    "error": stage_status["notes"],
                    "manager_status": manager_status,
                    "stage_status": stage_status,
                    "last_known_elapsed": last_known_elapsed,
                }
            raise RuntimeError(
                f"seekr status call failed: {payload.get('error')}\n"
                f"{payload.get('traceback', '')}"
            )
        stage_progress = payload.get("stage_progress") or {}

        # Carry per-anchor data through to stage_status so the time-limit
        # optimizer (which reads .seekrflow_job_status.json) can use it.
        for _k in ("anchor_times", "anchor_time_for_completion", "incomplete_anchors"):
            if _k in stage_progress:
                stage_status[_k] = stage_progress[_k]

        jobs = manager_status.get("jobs", [])
        running_now = len(jobs) > 0

        if benchmark_mode:
            if running_now:
                stage_status["finished"] = False
                stage_status["state"] = "started"
                stage_status["progress"] = 0.0
            elif stage_progress.get("finished", False):
                stage_status["finished"] = True
                stage_status["state"] = "completed"
                stage_status["progress"] = 1.0
            else:
                inferred_progress = _progress_from_stage_progress(stage_progress, anchor)
                if inferred_progress > 0.0:
                    stage_status["finished"] = True
                    stage_status["state"] = "completed"
                    stage_status["progress"] = 1.0
                    stage_status["notes"] = (
                        "benchmark mode inferred completion from recorded progress"
                    )
                else:
                    stage_status["finished"] = False
                    stage_status["state"] = "unstarted"
                    stage_status["progress"] = 0.0
        else:
            if stage_progress.get("finished", False):
                stage_status["finished"] = True
                stage_status["state"] = "completed"
                stage_status["progress"] = 1.0
            else:
                progress = _progress_from_stage_progress(stage_progress, anchor)
                stage_status["finished"] = False
                stage_status["progress"] = progress
                if running_now:
                    stage_status["state"] = "started"
                elif progress > 0.0:
                    stage_status["state"] = "started"
                else:
                    stage_status["state"] = "unstarted"

    except Exception as e:
        stage_status.update(
            {
                "finished": False,
                "state": "failed",
                "progress": 0.0,
                "notes": f"Unexpected error checking stage status: {e}",
            }
        )
        return {
            "success": False,
            "error": str(e),
            "manager_status": manager_status,
            "stage_status": stage_status,
            "last_known_elapsed": last_known_elapsed,
        }

    return {
        "success": True,
        "error": None,
        "manager_status": manager_status,
        "stage_status": stage_status,
        "last_known_elapsed": last_known_elapsed,
    }


def slurm_remote_run_workflow(args):
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
    partition = args[3]
    account = args[4]
    constraint = args[5]
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
    anchor_times_at_submission = args[16] if len(args) > 16 else None

    kwargs = {
        "working_dir": working_dir,
        "logdir": logdir,
        "partition": partition,
        "account": account,
        "constraint": constraint,
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
        "anchor_times_at_submission": anchor_times_at_submission
    }

    @dataclass
    class RunState:
        run_id: str
        workdir: str
        logdir: str
        name: str
        partition: str
        scheduler_options: Optional[str]
        time_limit: str
        account: str
        constraint: Optional[str]
        cpus_per_task: Optional[int]
        mem: Optional[str]
        array_spec: Optional[list[int]]
        jobid: str
        submitted_at: float
        n_tasks: int
        cmd_template: str   # per-member command script
        attempts: Dict[str, List[int]]  # jobid -> indices attempted
        parent_jobids: List[str]  # history
        #slurm_env: Dict[str, str] # extra env to export in script
        worker_init: Optional[str] = None
        state_file: Optional[str] = None
        anchor_times_at_submission: Optional[Dict[str, float]] = None
        last_known_elapsed: Optional[str] = None  # Added to store last known elapsed time

        def save(
                self,
                path: pathlib.Path
                ) -> None:
            ensure_dir(path.parent)
            self.state_file = str(path)
            with open(path, "w") as f:
                json.dump(asdict(self), f, indent=2)
            return

        @staticmethod
        def load(
            path: pathlib.Path
        ) -> "RunState":
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
        return workdir / ".slurm_runner"
    
    def ensure_dir(p: pathlib.Path) -> None:
        p.mkdir(parents=True, exist_ok=True)
    
    work_dir = pathlib.Path(working_dir)
    state_path = get_state_dir(work_dir) / f"{stage_name}_latest.json"

    def collapse_indices(
            idxs: List[int]
            ) -> str:
        """
        Collapse [0,1,2,5,6,9] -> '0-2,5-6,9'
        """
        if not idxs: return ""
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
            partition: str,
            scheduler_options: Optional[str],
            time_limit: str,
            account: Optional[str],
            constraint: Optional[str],
            array_spec: list[int] | None,
            logdir: pathlib.Path,
            workdir: pathlib.Path,
            cpus_per_task: Optional[int],
            mem: Optional[str],
            cmd_template: str,
            worker_init: str = "",
            ) -> str:
        ensure_dir(logdir)

        if array_spec is None:
            args = ["sbatch",
                    "-J", name,
                    "-p", partition,
                    "-t", time_limit,
                    "-o", f"{logdir}/%x_%A.out",
                    "-e", f"{logdir}/%x_%A.err",
                    "-D", f"{workdir}",
                    ]
        else:
            args = ["sbatch",
                    "-J", name,
                    "-p", partition,
                    "-t", time_limit,
                    "--array", collapse_indices(array_spec),
                    "-o", f"{logdir}/%x_%A_%a.out",
                    "-e", f"{logdir}/%x_%A_%a.err",
                    "-D", f"{workdir}",
                    ]
        if scheduler_options: args += [scheduler_options]
        if cpus_per_task: args += ["--cpus-per-task", str(cpus_per_task)]
        if mem: args += ["--mem", str(mem)]
        if account: args += ["-A", account]
        if constraint: args += ["--constraint", constraint]

        wrap_cmd = cmd_template
        full = shlex.quote(wrap_cmd) if not worker_init else shlex.quote(
            f"{worker_init}; {wrap_cmd}")
        args += ["--wrap", full]
        out = run(["bash", "-lc", " ".join(x for x in args)])
        stdout = out[1]
        return stdout.strip().split()[-1]

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
            partition=args["partition"],
            scheduler_options=args["scheduler_options"],
            time_limit=args["time_limit"],
            account=args["account"],
            constraint=args["constraint"],
            array_spec=array_spec,
            logdir=logdir,
            workdir=workdir,
            cpus_per_task=args["cpus_per_task"],
            mem=args["mem"],
            cmd_template=args["command_string"],
            worker_init=args["worker_init"],
        )
        st = RunState(
            run_id=run_id,
            workdir=str(workdir),
            logdir=str(logdir),
            name=args["name"],
            partition=args["partition"],
            scheduler_options=args["scheduler_options"],
            time_limit=args["time_limit"],
            account=args["account"],
            constraint=args["constraint"],
            cpus_per_task=args["cpus_per_task"],
            mem=args["mem"],
            array_spec=array_spec,
            jobid=jobid,
            submitted_at=time.time(),
            n_tasks=n_tasks,
            cmd_template=args["command_string"],
            attempts={jobid: list(range(n_tasks))},  # all indices attempted
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

def slurm_remote_cancel_workflow(args):
    """
    Cancel SLURM job(s) remotely by id and/or scheduler job name.

    Args: [working_dir, job_id, job_name?]
    """
    import shlex
    import subprocess
    from typing import List

    job_id = args[1] if len(args) > 1 else None
    job_name = args[2] if len(args) > 2 else None

    def run(cmd: List[str], check: bool = True) -> tuple:
        """Run command, return (rc, stdout, stderr)"""
        out = subprocess.run(cmd, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, text=True)
        if check and out.returncode != 0:
            raise RuntimeError(f"Command failed: ({out.returncode}): "
                             f"{' '.join(cmd)}\nSTDERR:\n{out.stderr}")
        return out.returncode, out.stdout.strip(), out.stderr.strip()

    if not job_id and not job_name:
        return {"success": False, "error": "job_id or job_name required"}

    results = []
    if job_id:
        results.append(run(
            ["bash", "-lc", f"scancel {shlex.quote(str(job_id))}"],
            check=False))
    if job_name:
        results.append(run(
            ["bash", "-lc", f"scancel --name={shlex.quote(str(job_name))}"],
            check=False))
    return {
        "success": True,
        "error": None,
        "jobid": job_id,
        "job_name": job_name,
        "cancel output:": results,
    }

def slurm_remote_force_rerun_workflow(args):
    """
    Cancel a stage's running SLURM job and reset its runner bookkeeping.

    Cancels any job recorded for the stage, waits until it is fully gone
    (so a resubmission cannot clobber a still-running job), and clears the
    ``.slurm_runner`` state files. Stage output cleanup is NOT done here; it
    is handled by seekr's own ``force_overwrite`` in the run command.

    Parameters
    ----------
    args : list
        List containing:
        - working_dir : str
            Path to the working directory containing model.xml
        - stage_name : str
            Name of the stage to reset

    Returns
    -------
    dict
        Result with success status and any errors
    """
    import os
    import json
    import glob
    import time
    import pathlib
    import subprocess
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
        partition: str
        scheduler_options: Optional[str]
        time_limit: str
        account: str
        constraint: Optional[str]
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
        last_known_elapsed: Optional[str] = None  # Added to store last known elapsed time

        @staticmethod
        def load(path: pathlib.Path) -> "RunState":
            with open(path, "r") as f:
                data = json.load(f)
            return RunState(**data)
    
    def get_state_dir(workdir: pathlib.Path) -> pathlib.Path:
        return workdir / ".slurm_runner"

    def cleanup_stage_state_files(workdir: pathlib.Path, stage: str) -> None:
        state_dir = get_state_dir(workdir)
        if not state_dir.exists():
            return
        for pattern in (f"{stage}_latest.json", f"{stage}_*.json"):
            for filepath in state_dir.glob(pattern):
                try:
                    filepath.unlink()
                except Exception:
                    pass
    
    def wait_until_gone(
            jobid: str,
            max_checks: int = 20,
            poll_seconds: float = 3.0,
            ) -> None:
        """Poll squeue until the job no longer appears (bounded)."""
        for _ in range(max_checks):
            rc, out, err = run(
                ["bash", "-lc", f"squeue -j {jobid} -h -o '%i'"], check=False)
            if rc != 0 or not out.strip():
                return
            time.sleep(poll_seconds)
        print(f"Warning: job {jobid} still visible after cancel; "
              f"proceeding anyway")

    result = {
        "success": False,
        "error": None,
        "canceled_job": None
    }
    
    work_dir = pathlib.Path(working_dir)
    state_path = get_state_dir(work_dir) / f"{stage_name}_latest.json"
    
    # Step 1: Cancel any running job for this stage and confirm it is gone,
    # so a resubmission cannot clobber a still-running job.
    if state_path is not None and state_path.exists():
        try:
            st = RunState.load(state_path)
            print(f"Canceling job {st.jobid} for stage {stage_name}...")
            rc, out, err = run(["bash", "-lc", f"scancel {st.jobid}"], check=False)
            if rc == 0:
                result["canceled_job"] = st.jobid
                print(f"Canceled job {st.jobid}")
            else:
                print(f"Job {st.jobid} may already be finished or canceled")
            wait_until_gone(st.jobid)
        except Exception as e:
            result["error"] = f"Failed to cancel job: {e}"
            print(f"Warning: {result['error']}")
    
    # Step 2: Clear runner state so the stage resubmits fresh. Stage output
    # cleanup is handled by seekr's own force_overwrite in the run command.
    try:
        cleanup_stage_state_files(work_dir, stage_name)
        print(f"Cleared runner state for {stage_name} stage")
        result["success"] = True
    except Exception as e:
        result["error"] = f"Failed to clear runner state: {e}"
        result["success"] = False
    
    return result