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

def pbs_remote_status_workflow(args):
    """
    Generic PBS status workflow.

    Mirrors slurm_remote_status_workflow: a single stage-agnostic status
    function that queries the PBS scheduler for the most recent job and
    runs a seekr-based progress check via a bash subprocess under the
    resource's worker_init.  This avoids depending on whatever Python the
    Globus Compute endpoint daemon happens to use (which may not have the
    seekr package).

    Args format:
        [working_dir, stage_name, benchmark_mode=False, anchor="any",
         swarm_id=None, worker_init=""]

    worker_init is a bash snippet that activates the env which has the
    `seekr` package installed (e.g. "source ~/.bashrc; conda activate SEEKR2").
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

    def find_latest_state_file(
            workdir: pathlib.Path, stage: str) -> Optional[pathlib.Path]:
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

    def _set_nested(d: dict, dotted_key: str, value: str) -> None:
        parts = dotted_key.split(".")
        cur = d
        for p in parts[:-1]:
            nxt = cur.get(p)
            if not isinstance(nxt, dict):
                nxt = {}
                cur[p] = nxt
            cur = nxt
        cur[parts[-1]] = value

    def _parse_qstat_text(text: str) -> dict:
        """
        Parse `qstat -f` plain-text output (universal across PBS versions)
        into the same {job_id: {key: value, ...}} shape that
        `qstat -f -F json` produces under qstat_data["Jobs"].  Nested keys
        like `resources_used.walltime` become nested dicts.
        """
        jobs: dict = {}
        current_id: Optional[str] = None
        current_key: Optional[str] = None
        current_val_parts: List[str] = []

        def _flush():
            if current_id is not None and current_key is not None:
                _set_nested(
                    jobs[current_id],
                    current_key,
                    "".join(current_val_parts).strip(),
                )

        for raw_line in text.splitlines():
            if raw_line.startswith("Job Id:"):
                _flush()
                current_id = raw_line.split(":", 1)[1].strip()
                jobs[current_id] = {}
                current_key = None
                current_val_parts = []
                continue
            if current_id is None:
                continue
            # Continuation lines start with whitespace and have no '='
            # until the parser hits the next "key = value".
            if raw_line.startswith("\t") and current_key is not None:
                current_val_parts.append(raw_line.lstrip())
                continue
            stripped = raw_line.strip()
            if not stripped:
                continue
            if "=" in stripped:
                _flush()
                key, _, val = stripped.partition("=")
                current_key = key.strip()
                current_val_parts = [val.strip()]
            else:
                # Loose continuation without leading tab
                if current_key is not None:
                    current_val_parts.append(stripped)
        _flush()
        return jobs

    def _qstat_jobs(jobid: str) -> tuple[List[dict], str]:
        """
        Query PBS for jobid (and any subjobs of an array job).  Try the
        JSON form first (PBS Pro 19+), fall back to text-form parsing
        which is supported by all PBS variants (Torque, OpenPBS, older
        PBS Pro).  Returns (rows, error_str).
        """
        # -t expands array jobs into one row per subjob; harmless for
        # non-array jobs.
        json_cmd = f"qstat -f -F json -t {shlex.quote(jobid)}"
        rc, out, err = run(["bash", "-lc", json_cmd], check=False)
        if rc == 0 and out:
            try:
                data = json.loads(out)
                jobs_dict = data.get("Jobs", {}) or {}
                rows = []
                for jid, info in jobs_dict.items():
                    row = dict(info) if isinstance(info, dict) else {}
                    row.setdefault("JobID", jid)
                    rows.append(row)
                return rows, ""
            except Exception:
                pass  # fall through to text-form parsing

        text_cmd = f"qstat -f -t {shlex.quote(jobid)}"
        rc, out, err = run(["bash", "-lc", text_cmd], check=False)
        if rc != 0:
            return [], err or f"qstat exited rc={rc}"
        jobs_dict = _parse_qstat_text(out)
        rows = []
        for jid, info in jobs_dict.items():
            info.setdefault("JobID", jid)
            rows.append(info)
        return rows, ""

    def _extract_walltime(row: dict) -> Optional[str]:
        ru = row.get("resources_used")
        if isinstance(ru, dict):
            wt = ru.get("walltime")
            if isinstance(wt, str) and wt.strip():
                return wt.strip()
        return None

    # PBS job_state codes -> SLURM-style uppercase strings expected by
    # seekr_run.determine_stage_manager_status (which reads row["State"]).
    # Reference: PBS Pro / Torque / OpenPBS use single-letter codes.
    _PBS_STATE_MAP = {
        "R": "RUNNING",
        "E": "RUNNING",   # Exiting, still consuming walltime
        "Q": "QUEUED",
        "H": "QUEUED",    # Held
        "W": "QUEUED",    # Waiting
        "S": "QUEUED",    # Suspended
        "T": "QUEUED",    # Transitioning
        "B": "RUNNING",   # Array job has at least one subjob running
        "F": "COMPLETED",
        "C": "COMPLETED", # Torque completed
        "X": "COMPLETED", # Subjob completed (array)
    }

    def _annotate_row_state(row: dict) -> None:
        """
        Mutate row to set a SLURM-style "State" key derived from PBS's
        job_state.  Leaves any existing "State" alone if already set.
        """
        if row.get("State"):
            return
        js = row.get("job_state") or row.get("Job_State") or ""
        if isinstance(js, str) and js.strip():
            code = js.strip()[0].upper()
            row["State"] = _PBS_STATE_MAP.get(code, code or "UNKNOWN")
        else:
            row["State"] = "UNKNOWN"

    work_dir = pathlib.Path(working_dir)
    manager_status = {
        "tool": "qstat",
        "timestamp": time.time(),
        "error": "",
        "jobs": [],
    }
    last_known_elapsed = None

    state_path = find_latest_state_file(work_dir, stage_name)
    if state_path is not None and state_path.exists():
        try:
            st = RunState.load(state_path)
            try:
                rows, qerr = _qstat_jobs(st.jobid)
                if rows:
                    for _r in rows:
                        _annotate_row_state(_r)
                    # Use the first row's walltime as the most recent
                    # benchmarking signal (sufficient for time-limit calc).
                    wt = _extract_walltime(rows[0])
                    if wt is not None:
                        last_known_elapsed = wt
                        try:
                            state_data = json.loads(state_path.read_text())
                            state_data["last_known_elapsed"] = last_known_elapsed
                            state_path.write_text(json.dumps(state_data, indent=2))
                        except Exception:
                            pass
                    manager_status["jobs"] = rows
                else:
                    if qerr:
                        manager_status["error"] = qerr
                    try:
                        state_data = json.loads(state_path.read_text())
                        last_known_elapsed = state_data.get("last_known_elapsed")
                    except Exception:
                        pass
            except Exception as e:
                manager_status["error"] = f"Failed to check PBS status: {e}"
                try:
                    state_data = json.loads(state_path.read_text())
                    last_known_elapsed = state_data.get("last_known_elapsed")
                except Exception:
                    pass
        except Exception as e:
            manager_status["error"] = f"Failed to load PBS state file: {e}"

    def _progress_from_stage_progress(
            stage_progress: dict, partitioned_anchor) -> float:
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
    # the same bash snippet used to activate the env for PBS jobs.
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
                inferred_progress = _progress_from_stage_progress(
                    stage_progress, anchor)
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
