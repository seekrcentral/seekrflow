"""
Direct workload manager functions for direct local/head node job submission, 
status checking, and cancellation.
"""

def direct_run_workload(args):
    """
    Submit a direct workload.
    """
    import os
    import sys
    import time
    import json
    import pathlib
    import subprocess
    import importlib.util
    from typing import List, Optional
    from dataclasses import dataclass, asdict

    manager_payload = args[0]

    root_dir = manager_payload["root_dir"]
    job_name = manager_payload["resource_payload"]["job_name"]
    cpus_per_task = manager_payload["resource_payload"]["cpus"]
    
    #job_specs = manager_payload["job_specs"]
    # Convert the dictionary of job specs to JobSpec objects.
    namespace = {}
    exec(manager_payload["runner_files"]["structures.py"], namespace)
    JobSpec = namespace["JobSpec"]
    job_specs = [
        JobSpec.from_dict(spec) if isinstance(spec, dict) else spec
        for spec in manager_payload["job_specs"]
    ]
    
    job_py_source = manager_payload["runner_files"]["job.py"]
    structures_py_source = manager_payload["runner_files"]["structures.py"]
    
    array_indices: list[int] = []
    stage_indices: set[int] = set()
    anchor_indices: set[int] = set()
    swarm_indices: set[int] = set()
    benchmark_mode = False
    for job_spec in job_specs:
        array_indices.append(job_spec.array_index)
        for stage_spec in job_spec.stage_specs:
            stage_indices.add(stage_spec.stage_index)
            if stage_spec.benchmark:
                benchmark_mode = True
        for run_unit in job_spec.run_units:
            if run_unit.anchor != "any":
                anchor_indices.add(run_unit.anchor)
            if run_unit.swarm_id is not None:
                swarm_indices.add(run_unit.swarm_id)
    stage_indices = list(stage_indices)
    if len(anchor_indices) > 0:
        anchor_indices = list(anchor_indices)
    else:
        anchor_indices = None
    if len(swarm_indices) > 0:
        swarm_indices = list(swarm_indices)
    else:
        swarm_indices = None
    n_tasks = len(array_indices)
    
    if n_tasks <= 0: 
        return {"success": False, "error": "No tasks to submit"}
    root_dir_path = pathlib.Path(root_dir)
    model_filename = root_dir_path / "model.json"
    
    def get_direct_dir(root_dir: pathlib.Path) -> pathlib.Path:
        return root_dir / ".direct_runner"

    def get_stage_dir(root_dir: pathlib.Path) -> pathlib.Path:
        return root_dir / ".stage_states"

    def get_internal_id():
        """
        Since job IDs aren't assigned until the job is submitted, we need a way to
        choose an internal, seekrflow-specific ID to use for Stage and Direct Info
        and State files. This can just be where the direct runner directory
        has its content file names split by underscore and the ID is after the first
        underscore or something.
        """
        direct_runner_dir = get_direct_dir(root_dir_path)
        if not direct_runner_dir.exists():
            return 0
        files = sorted(
            direct_runner_dir.glob("*.json"),
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )
        if files:
            ids = [file.stem.split("_")[-1] for file in files]
            return max(int(id) for id in ids) + 1
        return 0

    @dataclass
    class DirectInfo:
        internal_id: int # The seekrflow-specific ID for this job. Used for info/state file names.
        root_dir: str
        log_dir: str # The directory in root_dir to write log files to.
        cpus_per_task: int # Direct cpus per task
        n_tasks: int # The # runs in the job array. If array_spec is None, must be 1
        array_indices: Optional[list[int]] # The indices of the runs in the job array
        benchmark_mode: bool # Whether to run in benchmark mode
        filename: str
        direct_state_filename: str
        stage_info_filenames: List[str] # Multiple stages possible in a bundled job
        stage_state_filenames: List[str]
        # Below: assigned when the job is submitted
        process_ids: Optional[List[str]] = None
        #process_name: Optional[str] = None
        process_stdout_filenames: Optional[List[str]] = None
        process_stderr_filenames: Optional[List[str]] = None
        submitted_at: Optional[float] = None

        def save(
                self, 
                root_dir_path: pathlib.Path | None = None,
                path: pathlib.Path | None = None
                ) -> None:
            if path is None:
                if self.filename is None:
                    raise ValueError("DirectInfo filename is required.")
                self.filename = pathlib.Path(self.filename).name
                path_basename = pathlib.Path(self.filename)
                assert root_dir_path is not None, "Root directory path is required."
                path = root_dir_path / ".direct_runner" / path_basename
            else:
                path_basename = path.name
                self.filename = str(path_basename)
            # Creates root_dir/.direct_runner; fails if root_dir is missing
            path.parent.mkdir(exist_ok=True)
            path_tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
            with open(path_tmp, "w") as f:
                json.dump(asdict(self), f, indent=2)
            os.rename(path_tmp, path)

        @staticmethod
        def load(path: pathlib.Path) -> "DirectInfo":
            with open(path, "r") as f:
                return DirectInfo(**json.load(f))

        def update(self, **kwargs) -> None:
            for key, value in kwargs.items():
                setattr(self, key, value)
        
    @dataclass
    class DirectState:
        internal_id: int
        filename: str
        direct_info_filename: str
        state: str # 'pending', 'running', 'idle', 'failed'
        notes: Optional[str] = None
        last_timestamp: Optional[float] = None
        last_known_elapsed: Optional[List[str]] = None # One for each of n_tasks or stages?
        last_known_jobs: Optional[List[str]] = None

        def save(
                self, 
                root_dir_path: pathlib.Path | None = None,
                path: pathlib.Path | None = None
                ) -> None:
            if path is None:
                if self.filename is None:
                    raise ValueError("DirectState filename is required.")
                self.filename = pathlib.Path(self.filename).name
                path_basename = pathlib.Path(self.filename)
                assert root_dir_path is not None, "Root directory path is required."
                path = root_dir_path / ".direct_runner" / path_basename
            else:
                path_basename = path.name
                self.filename = str(path_basename)
            # Creates root_dir/.direct_runner; fails if root_dir is missing
            path.parent.mkdir(exist_ok=True)
            path_tmp = path.with_suffix(".tmp")
            with open(path_tmp, "w") as f:
                json.dump(asdict(self), f, indent=2)
            os.rename(path_tmp, path)

        @staticmethod
        def load(path: pathlib.Path) -> "DirectState":
            with open(path, "r") as f:
                return DirectState(**json.load(f))

        def update(self, **kwargs) -> None:
            for key, value in kwargs.items():
                setattr(self, key, value)

    # 1. perform preliminary checks: make sure model.json exists, fill out StageInfo 
    #    for each stage being submitted. Write out a preliminary DirectInfo, 
    #    DirectState, StageInfo, and StateState files.
    # NOTE: keeping asserts on worker - will be propagated back to direct
    if not root_dir_path.exists():
        return {"success": False, "error": "Root directory does not exist"}
    if not root_dir_path.is_dir():
        return {"success": False, "error": "Root directory is not a directory"}
    if not model_filename.exists():
        return {"success": False, "error": "Model file does not exist"}
    direct_dir_path = get_direct_dir(root_dir_path)
    stage_dir_path = get_stage_dir(root_dir_path)
    internal_id = get_internal_id()
    for job_spec in job_specs:
        job_spec.internal_id = internal_id
    LOG_DIR = "logs"
    log_dir_path = root_dir_path / LOG_DIR
    direct_info_basename = f"direct_info_{internal_id}.json"
    direct_state_basename = f"direct_state_{internal_id}.json"
    stage_info_filenames = []
    stage_state_filenames = []
    stage_info_basenames = []
    stage_state_basenames = []
    for stage_index in stage_indices:
        stage_info_filename = stage_dir_path / f"stage_info_{internal_id}_{stage_index}.json"
        stage_info_basenames.append(stage_info_filename.name)
        stage_info_filenames.append(str(stage_info_filename))
        for array_index in array_indices:
            stage_state_filename = stage_dir_path / f"stage_state_{internal_id}_{stage_index}_{array_index}.json"
            stage_state_basenames.append(stage_state_filename.name)
            stage_state_filenames.append(str(stage_state_filename))

    # Write the job.py and structures.py files to the root/.seekr_jobs/ directory.
    seekr_jobs_dir = root_dir_path / ".seekr_jobs"
    seekr_jobs_dir.mkdir(parents=True, exist_ok=True)
    job_py_path = seekr_jobs_dir / "job.py"
    job_py_path.write_text(job_py_source)
    structures_py_path = seekr_jobs_dir / "structures.py"
    structures_py_path.write_text(structures_py_source)

    # Import the job.py and structures.py files.
    spec = importlib.util.spec_from_file_location(
        "seekr_job_structures", structures_py_path)
    structures_module = importlib.util.module_from_spec(spec)
    sys.modules["seekr_job_structures"] = structures_module
    spec.loader.exec_module(structures_module)
    StageInfo = structures_module.StageInfo
    StageState = structures_module.StageState
    job_spec_paths = []
    for job_spec in job_specs:
        job_spec_path = seekr_jobs_dir / f"job_spec_{internal_id}_{job_spec.array_index}.json"
        job_spec.save(job_spec_path)
        job_spec_paths.append(str(job_spec_path))

    # Write out preliminary DirectInfo, DirectState, StageInfo, and StageState files
    direct_info = DirectInfo(
        internal_id=internal_id,
        root_dir=root_dir,
        log_dir=LOG_DIR,
        cpus_per_task=cpus_per_task,
        n_tasks=n_tasks,
        array_indices=array_indices,
        benchmark_mode=benchmark_mode,
        filename=direct_info_basename,
        direct_state_filename=direct_state_basename,
        stage_info_filenames=stage_info_basenames,
        stage_state_filenames=stage_state_basenames,
    )
    direct_info.save(root_dir_path=root_dir_path)
    direct_state = DirectState(
        internal_id=internal_id,
        filename=direct_state_basename,
        direct_info_filename=direct_info_basename,
        state="pending",
    )
    direct_state.save(root_dir_path=root_dir_path)
    
    starting_stage_infos = []
    starting_stage_states = []
    for i, stage_index in enumerate(stage_indices):
        starting_stage_info = StageInfo(
            internal_id=internal_id,
            stage_index=stage_index,
            anchor_indices=anchor_indices,
            swarm_indices=swarm_indices,
            filename=stage_info_basenames[i],
            starting_timestamp=time.time(),
        )
        starting_stage_info.save(path=root_dir_path / ".stage_info_states" / stage_info_basenames[i])
        
        for j, array_index in enumerate(array_indices):
            stage_array_index = i * len(array_indices) + j
            starting_stage_state = StageState(
                internal_id=internal_id,
                state="unknown",
                finished=False,
                stage_anchor_swarm_progress_list=None,
                status_info=None,
                notes=None,
                filename=stage_state_basenames[stage_array_index],
                stage_info_filename=stage_info_basenames[i],
                latest_timestamp=time.time(),
            )
            starting_stage_state.save(root_dir_path=root_dir_path, path=stage_dir_path / stage_state_basenames[stage_array_index])
        
    # 2. submit the job to direct
    log_dir_path.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(cpus_per_task)
    env["OPENMM_CPU_THREADS"] = str(cpus_per_task)

    stdout_basenames = []
    stderr_basenames = []
    process_ids = []
    for array_index in array_indices:
        spec_path = seekr_jobs_dir / f"job_spec_{internal_id}_{array_index}.json"
        stdout_basename = f"{job_name}_{internal_id}_{array_index}.out"
        stderr_basename = f"{job_name}_{internal_id}_{array_index}.err"
        stdout_basenames.append(stdout_basename)
        stderr_basenames.append(stderr_basename)
        stdout_f = open(log_dir_path / stdout_basename, "w")
        stderr_f = open(log_dir_path / stderr_basename, "w")
        try:
            proc = subprocess.Popen(
                [sys.executable, str(job_py_path), str(spec_path)],
                cwd=str(root_dir_path),
                stdout=stdout_f,
                stderr=stderr_f,
                env=env,
                start_new_session=True,
            )
        finally:
            stdout_f.close()
            stderr_f.close()
        process_ids.append(proc.pid)
    
    direct_info.process_ids = process_ids
    direct_info.process_stdout_filenames = stdout_basenames
    direct_info.process_stderr_filenames = stderr_basenames
    direct_info.submitted_at = time.time()
    direct_info.save(root_dir_path=root_dir_path)
    direct_state.state = "running"
    direct_state.last_timestamp = time.time()
    direct_state.last_known_jobs = [str(pid) for pid in process_ids]
    direct_state.save(root_dir_path=root_dir_path)

    return {
        "success": True, 
        "error": None, 
        "internal_id": internal_id,
        "process_ids": process_ids,
        "job_id": internal_id,
        "job_name": job_name,
    }

def direct_status_workload(args):
    """
    Direct status workflow

    Load the DirectInfo and DirectState files and return them.
    """
    import os
    import json
    import time
    import pathlib
    import subprocess
    import importlib.util
    from typing import List, Optional
    from dataclasses import dataclass, asdict

    @dataclass
    class DirectState:
        internal_id: int
        filename: str
        direct_info_filename: str
        state: str # 'pending', 'running', 'idle', 'failed'
        notes: Optional[str] = None
        last_timestamp: Optional[float] = None
        last_known_elapsed: Optional[List[str]] = None # One for each of n_tasks or stages?
        last_known_jobs: Optional[List[str]] = None

        def save(
                self, 
                root_dir_path: pathlib.Path | None = None,
                path: pathlib.Path | None = None
                ) -> None:
            if path is None:
                if self.filename is None:
                    raise ValueError("DirectState filename is required.")
                self.filename = pathlib.Path(self.filename).name
                path_basename = pathlib.Path(self.filename)
                assert root_dir_path is not None, "Root directory path is required."
                path = root_dir_path / ".direct_runner" / path_basename
            else:
                path_basename = path.name
                self.filename = str(path_basename)
            # Creates root_dir/.direct_runner; fails if root_dir is missing
            path.parent.mkdir(exist_ok=True)
            path_tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
            with open(path_tmp, "w") as f:
                json.dump(asdict(self), f, indent=2)
            os.rename(path_tmp, path)

        @staticmethod
        def load(path: pathlib.Path) -> "DirectState":
            with open(path, "r") as f:
                return DirectState(**json.load(f))

        def update(self, **kwargs) -> None:
            for key, value in kwargs.items():
                setattr(self, key, value)

    def parse_json_file_to_dict(path: pathlib.Path) -> dict:
        with open(path, "r") as f:
            return json.load(f)

    def pid_is_running(pid: int) -> bool:
        """True when pid is alive and not a zombie."""
        try:
            os.kill(int(pid), 0)
        except ProcessLookupError:
            return False
        except PermissionError:
            return True
        stat_path = pathlib.Path(f"/proc/{int(pid)}/stat")
        try:
            # /proc/pid/stat: "pid (comm) STATE ..."
            fields = stat_path.read_text().rsplit(")", 1)[-1].split()
            if fields and fields[0] == "Z":
                return False
        except OSError:
            pass
        return True

    def elapsed_since(submitted_at: float | None) -> str | None:
        if submitted_at is None:
            return None
        seconds = max(0, int(time.time() - float(submitted_at)))
        minutes, secs = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        return f"{hours}:{minutes:02d}:{secs:02d}"
    
    incoming_payload = args[0]
    return_payload = {}
    for system_name, payload in incoming_payload["system_payloads"].items():
        root_dir_path = pathlib.Path(payload["root_dir"])
        if not root_dir_path.exists():
            return {"success": False, "error": "Root directory does not exist"}
        if not root_dir_path.is_dir():
            return {"success": False, "error": "Root directory is not a directory"}

        direct_dir_path = root_dir_path / ".direct_runner"
        stage_dir_path = root_dir_path / ".stage_states"
        job_dicts_by_job_id = {}
        for job in payload["jobs"]:
            internal_id = job["internal_id"]
            #process_id = job["job_id"]
            stage_indices = job["stage_indices"]
            array_indices = job["array_indices"] or []
            info_path = direct_dir_path / f"direct_info_{internal_id}.json"
            state_path = direct_dir_path / f"direct_state_{internal_id}.json"
            if not info_path.exists() or not state_path.exists():
                return {
                    "success": False,
                    "error": f"Direct runner files missing for internal_id {internal_id}",
                }
            direct_info = parse_json_file_to_dict(info_path)
            direct_state = DirectState.load(state_path)
            process_ids = direct_info.get("process_ids") or []
            info_array_indices = direct_info.get("array_indices") or []
            pid_by_array_index = {
                index: pid for index, pid in zip(info_array_indices, process_ids)
            }
            elapsed = elapsed_since(direct_info.get("submitted_at"))

            direct_dicts_by_array_index = {}
            live_pids = []
            for array_index in array_indices:
                pid = pid_by_array_index.get(array_index)
                if pid is not None and pid_is_running(pid):
                    live_pids.append(str(pid))
                direct_dicts_by_array_index[array_index] = {
                    "internal_id": internal_id,
                    "state": "running" if pid is not None and pid_is_running(pid) else "idle",
                    "last_timestamp": time.time(),
                    "last_known_elapsed": elapsed,
                    "last_known_jobs": [str(pid)] if pid is not None else [],
                }
            
            array_states = {item["state"] for item in direct_dicts_by_array_index.values()}
            if "running" in array_states:
                direct_state.state = "running"
            else:
                direct_state.state = "idle"
                direct_state.last_timestamp = time.time()
                direct_state.last_known_elapsed = [elapsed] if elapsed is not None else None
                direct_state.last_known_jobs = live_pids
                direct_state.save(root_dir_path=root_dir_path, path=state_path)

            list_fields = (
                "stage_anchor_swarm_progress_list",
                "stage_anchor_swarm_starting_step_list",
                "stage_anchor_swarm_current_step_list",
                "stage_anchor_swarm_total_steps_list",
                "stage_anchor_swarm_time_of_first_progress_list",
                "stage_anchor_swarm_time_of_last_progress_list",
            )
            stage_dicts_by_stage_index = {}
            for stage_index in stage_indices:
                member_dicts = []
                for array_index in array_indices:
                    stage_state_path = stage_dir_path / (
                        f"stage_state_{internal_id}_{stage_index}_{array_index}.json")
                    if stage_state_path.exists():
                        member_dicts.append(parse_json_file_to_dict(stage_state_path))
                    else:
                        member_dicts.append({"state": "unstarted", "finished": False})
                states = [member.get("state") or "unstarted" for member in member_dicts]
                if "error" in states:
                    state = "error"
                elif states and all(member_state == "completed" for member_state in states):
                    state = "completed"
                elif any(member_state in ("started", "completed") for member_state in states):
                    state = "started"
                else:
                    state = "unstarted"
                merged = {
                    "internal_id": internal_id,
                    "state": state,
                    "finished": all(bool(member.get("finished")) for member in member_dicts),
                }
                for field in list_fields:
                    combined = []
                    for member in member_dicts:
                        combined.extend(member.get(field) or [])
                    merged[field] = combined or None
                stage_dicts_by_stage_index[stage_index] = merged
            
            stage_states = {
                entry["state"] for entry in stage_dicts_by_stage_index.values()
            }
            stages_completed = bool(stage_states) and stage_states == {"completed"}
            any_running = any(
                item["state"] == "running"
                for item in direct_dicts_by_array_index.values()
            )
            stderr_by_array_index = {
                index: name
                for index, name in zip(
                    info_array_indices,
                    direct_info.get("process_stderr_filenames") or [],
                )
            }
            log_dir = direct_info.get("log_dir") or "logs"
            if not any_running and not stages_completed:
                for array_index, item in direct_dicts_by_array_index.items():
                    stderr_name = stderr_by_array_index.get(array_index)
                    tail = ""
                    if stderr_name:
                        err_path = root_dir_path / log_dir / stderr_name
                        if err_path.is_file():
                            lines = err_path.read_text(errors="replace").splitlines()
                            tail = "\n".join(lines[-10:])
                    item["state"] = "failed"
                    item["error"] = tail or "process exited before the stage completed"

            job_dicts_by_job_id[job["job_id"]] = {
                "manager_dicts_by_array_index": direct_dicts_by_array_index,
                "stage_dicts_by_stage_index": stage_dicts_by_stage_index,
            }
        # Name it job_dicts_by_internal_id ??
        return_payload[system_name] = job_dicts_by_job_id

    return {
        "success": True,
        "error": None,
        "payload": return_payload,
    }

def direct_cancel_workload(args):
    import os
    import signal
    import json
    import pathlib
    incoming = args[0]
    canceled = []
    errors = []
    for payload in incoming["system_payloads"].values():
        root = pathlib.Path(payload["root_dir"])
        for job in payload["jobs"]:
            internal_id = job.get("internal_id", job.get("job_id"))
            info_path = root / ".direct_runner" / f"direct_info_{internal_id}.json"
            if not info_path.exists():
                errors.append(f"missing {info_path.name}")
                continue
            info = json.loads(info_path.read_text())
            for pid in info.get("process_ids") or []:
                try:
                    os.killpg(int(pid), signal.SIGTERM)
                    canceled.append(int(pid))
                except ProcessLookupError:
                    pass
                except Exception as error:
                    errors.append(f"{pid}: {error}")
    return {
        "success": not errors,
        "error": "; ".join(errors) or None,
        "payload": canceled,
    }