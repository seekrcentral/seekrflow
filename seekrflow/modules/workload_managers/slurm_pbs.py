"""
SLURM/PBS workload manager functions for remote job submission, status checking, 
and cancellation.
"""

# NOTE: notice that this was changed from 'workflow' to 'workload'
# TODO: combine with PBS to make it a single module. (DRY)
def slurm_pbs_run_workload(args):
    """
    Submit a remote workload to SLURM.
    """
    import os
    import sys
    import time
    import json
    import shlex
    import pathlib
    import subprocess
    import importlib.util
    from typing import List, Optional
    from dataclasses import dataclass, asdict

    manager_payload = args[0]

    root_dir = manager_payload["root_dir"]
    scheduler = manager_payload["resource_payload"]["scheduler"]
    partition_queue = manager_payload["resource_payload"]["partition_queue"]
    account = manager_payload["resource_payload"]["account"]
    constraint = manager_payload["resource_payload"].get("constraint")
    job_name = manager_payload["resource_payload"]["job_name"]
    cpus_per_task = manager_payload["resource_payload"]["cpus"]
    memory_mb_per_node = manager_payload["resource_payload"]["memory_mb"]
    time_limit = manager_payload["resource_payload"]["time_limit"]
    scheduler_options = manager_payload["resource_payload"]["scheduler_options"]
    worker_init = manager_payload["resource_payload"]["worker_init"]
    
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

    # Define helper functions for this workload
    def run(cmd: List[str], check: bool = True) -> tuple:
        """Run command, return (rc, stdout, stderr)"""
        out = subprocess.run(cmd, stdout=subprocess.PIPE, 
                           stderr=subprocess.PIPE, text=True)
        if check and out.returncode != 0:
            raise RuntimeError(f"Command failed: ({out.returncode}): "
                             f"{' '.join(cmd)}\nSTDERR:\n{out.stderr}")
        return out.returncode, out.stdout.strip(), out.stderr.strip()

    def get_batch_dir(root_dir: pathlib.Path) -> pathlib.Path:
        return root_dir / ".batch_runner"

    def get_stage_dir(root_dir: pathlib.Path) -> pathlib.Path:
        return root_dir / ".stage_states"

    def get_internal_id():
        """
        Since job IDs aren't assigned until the job is submitted, we need a way to
        choose an internal, seekrflow-specific ID to use for Stage and Batch Info
        and State files. This can just be where the batch runner directory
        has its content file names split by underscore and the ID is after the first
        underscore or something.
        """
        batch_runner_dir = get_batch_dir(root_dir_path)
        if not batch_runner_dir.exists():
            return 0
        files = sorted(
            batch_runner_dir.glob("batch_info_*.json"),
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )
        if files:
            ids = [file.stem.split("_")[-1] for file in files]
            return max(int(id) for id in ids) + 1
        return 0

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

    @dataclass
    class BatchInfo:
        internal_id: int # The seekrflow-specific ID for this job. Used for info/state file names.
        root_dir: str
        log_dir: str # The directory in root_dir to write log files to.
        worker_init: str # what to execute on the worker (head node) before job submit
        partition_queue: str # SLURM/PBS partition/queue
        time_limit: str # SLURM/PBS time limit for this job
        account: str # SLURM/PBS account
        scheduler_options: str # SLURM/PBS scheduler extra options
        constraint: Optional[str] # SLURM/PBS constraint
        cpus_per_task: int # SLURM/PBS cpus per task
        memory_mb: Optional[str] # SLURM/PBS memory per node
        n_tasks: int # The # runs in the job array. If array_spec is None, must be 1
        array_indices: Optional[list[int]] # The indices of the runs in the job array
        benchmark_mode: bool # Whether to run in benchmark mode
        filename: str
        stage_info_filenames: List[str] # Multiple stages possible in a bundled job
        stage_state_filenames: List[str]
        # Below: assigned when the job is submitted
        job_id: Optional[str] = None
        job_name: Optional[str] = None
        job_stdout_filenames: Optional[List[str]] = None
        job_stderr_filenames: Optional[List[str]] = None
        submitted_at: Optional[float] = None

        def save(
                self, 
                root_dir_path: pathlib.Path | None = None,
                path: pathlib.Path | None = None
                ) -> None:
            if path is None:
                if self.filename is None:
                    raise ValueError("BatchInfo filename is required.")
                self.filename = pathlib.Path(self.filename).name
                path_basename = pathlib.Path(self.filename)
                assert root_dir_path is not None, "Root directory path is required."
                path = root_dir_path / ".batch_runner" / path_basename
            else:
                path_basename = path.name
                self.filename = str(path_basename)
            # Creates root_dir/.batch_runner; fails if root_dir is missing
            path.parent.mkdir(exist_ok=True)
            path_tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
            with open(path_tmp, "w") as f:
                json.dump(asdict(self), f, indent=2)
            os.rename(path_tmp, path)

        @staticmethod
        def load(path: pathlib.Path) -> "BatchInfo":
            with open(path, "r") as f:
                return BatchInfo(**json.load(f))

        def update(self, **kwargs) -> None:
            for key, value in kwargs.items():
                setattr(self, key, value)
        
    @dataclass
    class BatchState:
        internal_id: int
        filename: str
        batch_info_filename: str
        state: str # 'pending', 'queued', 'running', 'queued/running', 'idle', 'failed', 'cancelled'
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
                    raise ValueError("BatchState filename is required.")
                self.filename = pathlib.Path(self.filename).name
                path_basename = pathlib.Path(self.filename)
                assert root_dir_path is not None, "Root directory path is required."
                path = root_dir_path / ".batch_runner" / path_basename
            else:
                path_basename = path.name
                self.filename = str(path_basename)
            # Creates root_dir/.batch_runner; fails if root_dir is missing
            path.parent.mkdir(exist_ok=True)
            path_tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
            with open(path_tmp, "w") as f:
                json.dump(asdict(self), f, indent=2)
            os.rename(path_tmp, path)

        @staticmethod
        def load(path: pathlib.Path) -> "BatchState":
            with open(path, "r") as f:
                return BatchState(**json.load(f))

        def update(self, **kwargs) -> None:
            for key, value in kwargs.items():
                setattr(self, key, value)

    # 1. perform preliminary checks: make sure model.json exists, fill out StageInfo 
    #    for each stage being submitted. Write out a preliminary SlurmInfo, 
    #    SlurmState, StageInfo, and StateState files.
    # NOTE: keeping asserts on worker - will be propagated back to local
    if not root_dir_path.exists():
        return {"success": False, "error": "Root directory does not exist"}
    if not root_dir_path.is_dir():
        return {"success": False, "error": "Root directory is not a directory"}
    if not model_filename.exists():
        return {"success": False, "error": "Model file does not exist"}
    batch_dir_path = get_batch_dir(root_dir_path)
    stage_dir_path = get_stage_dir(root_dir_path)
    internal_id = get_internal_id()
    for job_spec in job_specs:
        job_spec.internal_id = internal_id
    LOG_DIR = "logs"
    log_dir_path = root_dir_path / LOG_DIR
    batch_info_basename = f"batch_info_{internal_id}.json"
    batch_state_basenames = []
    for array_index in array_indices:
        batch_state_basenames.append(f"batch_state_{internal_id}_{array_index}.json")

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

    # Write out preliminary SlurmInfo, SlurmState, StageInfo, and StageState files
    batch_info = BatchInfo(
        internal_id=internal_id,
        root_dir=root_dir,
        log_dir=LOG_DIR,
        worker_init=worker_init,
        partition_queue=partition_queue,
        time_limit=time_limit,
        account=account,
        scheduler_options=scheduler_options,
        constraint=constraint,
        cpus_per_task=cpus_per_task,
        memory_mb=memory_mb_per_node,
        n_tasks=n_tasks,
        array_indices=array_indices,
        benchmark_mode=benchmark_mode,
        filename=batch_info_basename,
        stage_info_filenames=stage_info_basenames,
        stage_state_filenames=stage_state_basenames,
    )
    batch_info.save(root_dir_path=root_dir_path)
    batch_states = []
    for array_index, batch_state_basename in zip(array_indices, batch_state_basenames):
        batch_state = BatchState(
            internal_id=internal_id,
            filename=batch_state_basename,
            batch_info_filename=batch_info_basename,
            state="pending",
        )
        batch_state.save(root_dir_path=root_dir_path)
        batch_states.append(batch_state)
    
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
        starting_stage_info.save(path = root_dir_path / ".stage_info_states" / stage_info_basenames[i])
        
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
        
        #stage_info_command = f"python {stage_dir_path}/remote_status_monitor.py "\
        #    f"-r {root_dir} -I {internal_id} -S {stage_index} "\
        #    f"-i {starting_stage_info_json} -s {starting_stage_state_json}"
        #run(["bash", "-lc", stage_info_command])

    # 2. submit the job to SLURM/PBS
    log_dir_path.mkdir(parents=True, exist_ok=True)
    stdout_basenames = []
    stderr_basenames = []
    for array_index in array_indices:
        stdout_basenames.append(f"{job_name}_{internal_id}_{array_index}.out")
        stderr_basenames.append(f"{job_name}_{internal_id}_{array_index}.err")
    array_spec = collapse_indices(array_indices)
    
    if scheduler == "slurm":
        spec_path = seekr_jobs_dir / f"job_spec_{internal_id}_${{SLURM_ARRAY_TASK_ID}}.json"
        wrap_cmd = f"python {job_py_path} {spec_path}"
        if worker_init:
            wrap_cmd = f"{worker_init}; {wrap_cmd}"
        slurm_args = [
            "sbatch",
            "-J", job_name,
            "-p", partition_queue,
            "-t", time_limit,
            "--array", array_spec,
            "-o", f"{log_dir_path}/{job_name}_{internal_id}_%a.out",
            "-e", f"{log_dir_path}/{job_name}_{internal_id}_%a.err",
            "-D", f"{root_dir_path}",
            "--cpus-per-task", str(cpus_per_task),
        ]
        if memory_mb_per_node: 
            slurm_args += ["--mem", str(memory_mb_per_node)]
        if account: 
            slurm_args += ["--account", account]
        if constraint: 
            slurm_args += ["--constraint", constraint]
        if scheduler_options: 
            slurm_args += shlex.split(scheduler_options)

        slurm_args += ["--wrap", "bash -lc " + shlex.quote(f"{wrap_cmd}; rc=$?; exit $rc")]
        out = run(["bash", "-lc", " ".join(shlex.quote(x) for x in slurm_args)])
        job_id = out[1].strip().split()[-1]
    
    elif scheduler == "pbs":
        #spec_path = seekr_jobs_dir / f"job_spec_{internal_id}_${{PBS_ARRAY_INDEX}}.json"
        if "," in array_spec:
            raise RuntimeError(
                f"PBS -J accepts one range, not {array_spec!r}"
            )
        bounds = array_spec if "-" in array_spec else f"{array_spec}-{array_spec}"
        start_s, _, end_s = bounds.partition("-")
        start_i, end_i = int(start_s), int(end_s)
        one_task = start_i == end_i
        script_lines = [
            "#!/bin/bash",
            f"#PBS -N {job_name}",
            "#PBS -V",
            f"#PBS -q {partition_queue}",
            f"#PBS -l walltime={time_limit}",
            f"#PBS -l nodes=1:ppn={cpus_per_task or 1}",
            #f"#PBS -J {pbs_array_spec}",
            f"#PBS -o {log_dir_path}/",
            f"#PBS -e {log_dir_path}/",
        ]
        if not one_task:
            script_lines.append(f"#PBS -J {start_i + 1}-{end_i + 1}")

        if account:
            script_lines.append(f"#PBS -A {account}")
        if memory_mb_per_node:
            script_lines.append(f"#PBS -l mem={memory_mb_per_node}mb")
        if scheduler_options:
            script_lines.extend(
                line.strip() for line in scheduler_options.splitlines() if line.strip())
        script_lines.append(f"cd {root_dir_path}")
        if worker_init:
            script_lines.append(worker_init)
        if one_task:
            spec_name = f"job_spec_{internal_id}_{start_i}.json"
        else:
            spec_name = f"job_spec_{internal_id}_$((PBS_ARRAY_INDEX - 1)).json"
        script_lines.append(f"python {job_py_path} {seekr_jobs_dir}/{spec_name}")
        script_lines.append("rc=$?; exit $rc")
        qsub = subprocess.run(
            ["qsub"],
            input="\n".join(script_lines) + "\n",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if qsub.returncode != 0:
            raise RuntimeError(f"qsub failed: {qsub.stderr}")
        job_id = qsub.stdout.strip().replace("[]", "")
        if "." in job_id:
            job_id = job_id.split(".")[0]
    else:
        raise ValueError(f"Unknown scheduler: {scheduler!r}")


    batch_info.job_id = job_id
    batch_info.job_name = job_name
    batch_info.job_stdout_filenames = stdout_basenames
    batch_info.job_stderr_filenames = stderr_basenames
    batch_info.submitted_at = time.time()
    batch_info.save(root_dir_path=root_dir_path)
    for batch_state in batch_states:
        batch_state.state = "queued"
        batch_state.last_timestamp = time.time()
        batch_state.save(root_dir_path=root_dir_path)
    return {"success": True, "error": None, "internal_id": internal_id, "job_id": job_id,
            "job_name": job_name}

def slurm_pbs_status_workload(args):
    """
    SLURM/PBS status workload.

    Read the SLURM/PBS Stage Info/State files and return them.
    """
    import os
    import json
    import shlex
    import pathlib
    import subprocess
    from dataclasses import dataclass, asdict
    from typing import List, Optional

    @dataclass
    class BatchState:
        internal_id: int
        filename: str
        batch_info_filename: str
        state: str # 'pending', 'queued', 'running', 'queued/running', 'idle', 'failed', 'cancelled'
        notes: Optional[str] = None
        last_timestamp: Optional[float] = None
        last_known_elapsed: Optional[str] = None # One for each of n_tasks or stages?
        last_known_jobs: Optional[List[str]] = None

        def save(
                self, 
                root_dir_path: pathlib.Path | None = None,
                path: pathlib.Path | None = None
                ) -> None:
            if path is None:
                if self.filename is None:
                    raise ValueError("BatchState filename is required.")
                self.filename = pathlib.Path(self.filename).name
                path_basename = pathlib.Path(self.filename)
                assert root_dir_path is not None, "Root directory path is required."
                path = root_dir_path / ".batch_runner" / path_basename
            else:
                path_basename = path.name
                self.filename = str(path_basename)
            # Creates root_dir/.batch_runner; fails if root_dir is missing
            path.parent.mkdir(exist_ok=True)
            path_tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
            with open(path_tmp, "w") as f:
                json.dump(asdict(self), f, indent=2)
            os.rename(path_tmp, path)

        @staticmethod
        def load(path: pathlib.Path) -> "BatchState":
            with open(path, "r") as f:
                return BatchState(**json.load(f))

        def update(self, **kwargs) -> None:
            for key, value in kwargs.items():
                setattr(self, key, value)

    def run(cmd: List[str], check: bool = True) -> tuple:
        """Run command, return (rc, stdout, stderr)"""
        out = subprocess.run(cmd, stdout=subprocess.PIPE, 
                           stderr=subprocess.PIPE, text=True)
        if check and out.returncode != 0:
            raise RuntimeError(f"Command failed: ({out.returncode}): "
                             f"{' '.join(cmd)}\nSTDERR:\n{out.stderr}")
        return out.returncode, out.stdout.strip(), out.stderr.strip()
    # SLURM-specific functions
    def parse_raw_state_slurm(raw_state: str) -> str:
        raw_state = raw_state.strip().lower()
        if raw_state == "running":
            return "running"
        elif (raw_state == "queued") or (raw_state == "pending") or (raw_state == "configuring") \
                or (raw_state == "requeue_fed") or (raw_state == "requeue_hold")\
                or (raw_state == "requeued"):
            return "queued"
        elif (raw_state == "boot_fail") or (raw_state == "failed") or (raw_state == "node_fail") \
                or (raw_state == "out_of_memory"):
            return "failed"
        elif (raw_state == "cancelled") or (raw_state == "preempted") or (raw_state == "stopped"):
            return "cancelled"
        elif (raw_state == "completing") or (raw_state == "completed"):
            return "completed"
        elif (raw_state == "timeout") or (raw_state == "deadline"):
            #return "idle" # what I would have normally done
            return "failed"
        else:
            return "unknown"

    def parse_squeue_rows_slurm(output: str) -> tuple[dict, dict]:
        state_by_task = {}
        elapsed_by_task = {}
        for line in (output or "").splitlines():
            job_id, array_index, raw_state, elapsed = line.split("|")
            state = parse_raw_state_slurm(raw_state)
            if array_index.isdigit():
                index = int(array_index)
            else:
                index = 0
            state_by_task[(job_id, index)] = state
            elapsed_by_task[(job_id, index)] = elapsed
        return state_by_task, elapsed_by_task

    def call_squeue_and_parse_slurm():
        error_notes = ""
        squeue_fmt = "%A|%a|%T|%M"
        try:
            user = os.environ.get("USER") or os.environ.get("LOGNAME") or ""
            user_part = f"-u {shlex.quote(user)} " if user else ""
            cmd = (
                f"squeue {user_part} -r -h -o '{squeue_fmt}'"
            )
            rc, out, err = run(["bash", "-lc", cmd], check=False)
            if rc == 0:
                state_by_task, elapsed_by_task = parse_squeue_rows_slurm(out)
                return {
                    "success": True, "error": None, "state_by_task": state_by_task, 
                    "elapsed_by_task": elapsed_by_task}
            else:
                error_note = err or (f"squeue failed for job {{batch_info.job_name}} (rc={rc})")
                error_notes += f"{error_note}\n"

        except Exception as e:
            error_notes += f"Failed to check SLURM status: {e}\n"

        return {"success": False, "error": error_notes}
    
    # PBS-specific functions
    def parse_pbs_job_id(full_id: str) -> tuple[str, int | None]:
        """'12345[3].server' -> ('12345', 3). A parent id has no index."""
        base = full_id.split(".", 1)[0]
        if "[" not in base:
            return base, None
        parent, _, rest = base.partition("[")
        index = rest.rstrip("]")
        if index.isdigit():
            return parent, int(index)
        return parent, None

    def parse_pbs_state(job_state: str, exit_status: str | None) -> str:
        code = (job_state or "").strip()[:1].upper()
        if code in {"R", "E"}:
            return "running"
        if code in {"Q", "H", "W", "S", "T", "B"}:
            return "queued"
        if code in {"F", "C", "X"}:
            if exit_status not in (None, "", "0"):
                return "failed"
            return "completed"
        return "unknown"

    def parse_qstat_text(text: str) -> dict:
        jobs: dict = {}
        current_id = None
        current_key = None
        current_val: list[str] = []
        def flush():
            if current_id is None or current_key is None:
                return
            parts = current_key.split(".")
            cursor = jobs[current_id]
            for part in parts[:-1]:
                cursor = cursor.setdefault(part, {})
            cursor[parts[-1]] = "".join(current_val).strip()
        for raw_line in text.splitlines():
            if raw_line.startswith("Job Id:"):
                flush()
                current_id = raw_line.split(":", 1)[1].strip()
                jobs[current_id] = {}
                current_key = None
                current_val = []
                continue
            if current_id is None:
                continue
            if raw_line.startswith("\t") and current_key is not None:
                current_val.append(raw_line.strip())
                continue
            stripped = raw_line.strip()
            if not stripped or "=" not in stripped:
                continue
            flush()
            key, _, val = stripped.partition("=")
            current_key = key.strip()
            current_val = [val.strip()]
        flush()
        return jobs

    def call_qstat_and_parse_pbs():
        wanted = []
        for payload in incoming_payload["system_payloads"].values():
            for job in payload.get("jobs", []):
                for array_index in job.get("array_indices") or []:
                    wanted.append((job["job_id"], array_index))
        job_ids = list(dict.fromkeys(str(job_id) for job_id, _ in wanted))
        rows = []
        for job_id in job_ids:
            rc, out, err = run(
                ["bash", "-lc", f"qstat -f -t {shlex.quote(job_id)}"],
                check=False,
            )
            if rc != 0 or not out:
                # Unknown Job Id: leave it out so the shared loop marks it idle.
                continue
            for full_id, info in parse_qstat_text(out).items():
                parent, index = parse_pbs_job_id(full_id)
                raw_state = info.get("job_state") or info.get("Job_State") or ""
                exit_status = info.get("Exit_status") or info.get("exit_status")
                walltime = None
                used = info.get("resources_used")
                if isinstance(used, dict):
                    walltime = used.get("walltime")
                rows.append((parent, index, parse_pbs_state(raw_state, exit_status), walltime))
        state_by_task = {}
        elapsed_by_task = {}
        parent_state = {}
        parent_elapsed = {}
        for parent, index, state, walltime in rows:
            if index is None:
                parent_state[parent] = state
                parent_elapsed[parent] = walltime
                continue
            for job_id, array_index in wanted:
                if str(job_id) == parent and int(array_index) == index-1:
                    state_by_task[(job_id, array_index)] = state
                    elapsed_by_task[(job_id, array_index)] = walltime
        for job_id, array_index in wanted:
            if (job_id, array_index) in state_by_task:
                continue
            parent = parent_state.get(str(job_id))
            if parent in {"queued", "running"}:
                state_by_task[(job_id, array_index)] = parent
                elapsed_by_task[(job_id, array_index)] = parent_elapsed.get(str(job_id))
        return {
            "success": True,
            "error": None,
            "state_by_task": state_by_task,
            "elapsed_by_task": elapsed_by_task,
        }

    def parse_json_file_to_dict(path: pathlib.Path) -> dict:
        with open(path, "r") as f:
            return json.load(f)
    
    incoming_payload = args[0]
    scheduler = incoming_payload["scheduler"]

    if scheduler == "slurm":
        # First, call squeue
        squeue_result = call_squeue_and_parse_slurm()
        if not squeue_result["success"]:
            return {"success": False, "error": squeue_result["error"]}
        state_by_task = squeue_result["state_by_task"]
        elapsed_by_task = squeue_result["elapsed_by_task"]
    elif scheduler == "pbs":
        qstat_result = call_qstat_and_parse_pbs()
        if not qstat_result["success"]:
            return {"success": False, "error": qstat_result["error"]}
        state_by_task = qstat_result["state_by_task"]
        elapsed_by_task = qstat_result["elapsed_by_task"]
    else:
        return {"success": False, "error": "Invalid scheduler"}

    return_payload = {}
    for system_name, payload in incoming_payload["system_payloads"].items():
        root_dir = payload["root_dir"]
        root_dir_path = pathlib.Path(root_dir)
        if not root_dir_path.exists():
            return {"success": False, "error": "Root directory does not exist"}
        if not root_dir_path.is_dir():
            return {"success": False, "error": "Root directory is not a directory"}
        
        job_dicts_by_job_id = {}
        for job in payload["jobs"]:
            internal_id = job["internal_id"]
            job_id = job["job_id"]
            stage_indices = job["stage_indices"]
            array_indices = job["array_indices"]

            batch_dir_path = root_dir_path / ".batch_runner"
            stage_dir_path = root_dir_path / ".stage_states"
            batch_dicts_by_array_index = {}
            for array_index in array_indices:
                #slurm_info_path = slurm_dir_path / f"slurm_info_{internal_id}.json"
                batch_state_path = batch_dir_path / f"batch_state_{internal_id}_{array_index}.json"
                # Modify state based on squeue/qstat output
                batch_state = BatchState.load(batch_state_path)
                if (job_id, array_index) in state_by_task:
                    batch_state.state = state_by_task[(job_id, array_index)]
                    batch_state.last_known_elapsed = elapsed_by_task[(job_id, array_index)]
                    batch_state.last_known_jobs = [job_id]
                else:
                    finished = True
                    for stage_index in stage_indices:
                        stage_state_filename = (
                            stage_dir_path
                            / f"stage_state_{internal_id}_{stage_index}_{array_index}.json"
                        )
                        if not stage_state_filename.exists():
                            finished = False
                            break
                        member = parse_json_file_to_dict(stage_state_filename)
                        if member.get("state") != "completed" or not member.get("finished"):
                            finished = False
                            break
                    batch_state.state = "completed" if finished else "failed"
                    batch_state.last_known_jobs = []

                batch_state.save(root_dir_path=root_dir_path, path=batch_state_path)
                batch_state_dict = asdict(batch_state)
                batch_dicts_by_array_index[array_index] = batch_state_dict
            
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
                    stage_state_filename = stage_dir_path \
                        / f"stage_state_{internal_id}_{stage_index}_{array_index}.json"
                    if not stage_state_filename.exists():
                        member_dicts.append({"state": "unstarted", "finished": False})
                        continue
                    member_dicts.append(parse_json_file_to_dict(stage_state_filename))
                
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
                timestamps = [
                    member.get("latest_timestamp") 
                    for member in member_dicts
                    if member.get("latest_timestamp") is not None
                ]
                merged["latest_timestamp"] = max(timestamps) if timestamps else None
                stage_dicts_by_stage_index[stage_index] = merged

            job_dicts_by_job_id[job_id] = {
                "manager_dicts_by_array_index": batch_dicts_by_array_index,
                "stage_dicts_by_stage_index": stage_dicts_by_stage_index,
            }
    
        return_payload[system_name] = job_dicts_by_job_id
    
    
    result = {
        "success": True,
        "error": None,
        "payload": return_payload,
    }
    return result

def slurm_pbs_cancel_workload(args):
    """
    Cancel SLURM/PBS job(s) remotely by id and/or scheduler job name.

    Args: [root_dir, scheduler, remove_json_files?]
    """
    import shlex
    import pathlib
    import subprocess
    from typing import List

    def run(cmd: List[str], check: bool = True) -> tuple:
        """Run command, return (rc, stdout, stderr)"""
        out = subprocess.run(cmd, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, text=True)
        if check and out.returncode != 0:
            raise RuntimeError(f"Command failed: ({out.returncode}): "
                             f"{' '.join(cmd)}\nSTDERR:\n{out.stderr}")
        return out.returncode, out.stdout.strip(), out.stderr.strip()

    def cleanup_stage_state_files(root_dir_path: pathlib.Path) -> None:
        for dirname in (".batch_runner", ".stage_states"):
            directory = root_dir_path / dirname
            if not directory.exists():
                continue
            for filepath in directory.glob("*.json"):
                filepath.unlink()
    
    incoming_payload = args[0]
    scheduler = incoming_payload["scheduler"]
    remove_json_files = incoming_payload["remove_json_files"]
    
    for system_name, payload in incoming_payload["system_payloads"].items():
        root_dir = payload["root_dir"]
        root_dir_path = pathlib.Path(root_dir)
        if not root_dir_path.exists():
            return {"success": False, "error": "Root directory does not exist"}
        if not root_dir_path.is_dir():
            return {"success": False, "error": "Root directory is not a directory"}
        
        for job in payload["jobs"]:
            job_id = job.get("job_id", None)
            job_name = job.get("job_name", None)

            if not job_id and not job_name:
                return {"success": False, "error": "job_id or job_name required"}

            results = []
            if job_id:
                if scheduler == "slurm":
                    results.append(run(
                        ["bash", "-lc", f"scancel {shlex.quote(str(job_id))}"],
                        check=False))
                elif scheduler == "pbs":
                    results.append(run(
                        ["bash", "-lc", f"qdel {shlex.quote(str(job_id))}"],
                        check=False))
                else:
                    return {"success": False, "error": "Invalid scheduler"}
            if job_name:
                if scheduler == "slurm":
                    results.append(run(
                        ["bash", "-lc", f"scancel --name={shlex.quote(str(job_name))}"],
                        check=False))
                elif scheduler == "pbs":
                    # First look up the job and then cancel
                    rc, out, err = run(
                        ["bash", "-lc", f"qselect -N {shlex.quote(str(job_name))}"],
                        check=False)
                    if rc == 0:
                        for found_id in out.split():
                            results.append(run(
                                ["bash", "-lc", f"qdel {shlex.quote(found_id)}"],
                                check=False))
                else:
                    return {"success": False, "error": "Invalid scheduler"}

    result = {
        "success": True,
        "error": None,
        "payload": results,
    }
    
    for system_name, payload in incoming_payload.items():
        root_dir = payload["root_dir"]
        root_dir_path = pathlib.Path(root_dir)
        if remove_json_files:
            try:
                cleanup_stage_state_files(root_dir_path)
                result["success"] = True
            except Exception as e:
                result["error"] = f"Failed to clear runner state: {e}"
                result["success"] = False

    return result
   