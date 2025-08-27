"""

"""

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
    #import seekr2.modules.common_base as seekr2_base
    #import seekr2.run as seekr2_run
    
    #def bd_finished_wkflow(
    #        model: seekr2_base.Model,
    #        ) -> bool: 
    #    """
    #    Check if the BD stage has finished. This is done by reading the
    #    results XML files and seeing how many BD steps have elapsed.
    #    """
    #    bd_milestone_info_to_run = seekr2_run.choose_next_simulation_browndye2(
    #            model, "any_bd", None, False, None)
    #    if len(bd_milestone_info_to_run) > 0:
    #        return False
    #    else:
    #        return True

    #def hidr_finished_wkflow(
    #        model: seekr2_base.Model
    #        ) -> bool:
    #    """
    #    Check if the Hidr stage has finished. This is done by reading the model files
    #    and seeing if all anchors have starting PDB files.
    #    """
    #    for alpha, anchor in enumerate(model.anchors):
    #        if anchor.bulkstate:
    #            continue
    #        anchor_pdb_filename = seekr2_base.get_anchor_pdb_filename(anchor)
    #        if anchor_pdb_filename == "":
    #            return False
    #        
    #    return True

    #def seekr_anchors_to_run_wkflow(
    #        model: seekr2_base.Model
    #        ) -> bool:
    #    """
    #    Check if SEEKR has finished.
    #    """
    #    md_info_to_run = seekr2_run.choose_next_simulation_openmm(
    #        model, "any_md", None, None, None, None, False, False, None)
    #    anchor_indices_to_run = []
    #    for md_info in md_info_to_run:
    #        anchor_indices_to_run.append(md_info[2])
    #    return anchor_indices_to_run

    working_dir = args[0]
    logdir = args[1]
    partition = args[2]
    account = args[3]
    constraint = args[4]
    name = args[5]
    cpus_per_task = args[6]
    mem_per_node = args[7]
    resubmit_on_timeout = args[8]
    stop_flag = args[9]
    time_limit = args[10]
    scheduler_options = args[11]
    worker_init = args[12]
    #usage_tracking = args[13]
    command_string = args[14]
    indices = args[15]
    model_filename = args[16]
    workflow_type = args[17]

    kwargs = {
        "working_dir": working_dir,
        "logdir": logdir,
        "partition": partition,
        "account": account,
        "constraint": constraint,
        "name": name,
        "cpus_per_task": cpus_per_task,
        "mem": mem_per_node,
        "resubmit_on_timeout": resubmit_on_timeout,
        "stop_flag": stop_flag,
        "time_limit": time_limit,
        "scheduler_options": scheduler_options,
        "worker_init": worker_init,
        #"command_string": command_string,
        "indices": indices,
        "model_filename": model_filename,
        "workflow_type": workflow_type,
    }

    def run(
            cmd: List[str], 
            check: bool = True
        ) -> str:
        """
        Run a shell command, and if check is True, then check the return code.
        If there was an error - output the error to standard error.
        """
        out = subprocess.run(cmd, stdout=subprocess.PIPE, 
                             stderr=subprocess.PIPE, text=True)
        if check and out.returncode != 0:
            raise RuntimeError(f"Command failed: ({out.returncode}): "\
                               +f"{' '.join(cmd)}\nSTDERR:\n{out.stderr}")
        return out.stdout.strip()
    
    def which(
            prog: str,
            ) -> bool:
        """
        See whether a program can be executed without errors.
        """
        return subprocess.run(["bash", "-lc", f"command -v {shlex.quote(prog)} "\
                               "> /dev/null 2>&1"]).returncode == 0

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
    
    def ensure_dir(
            p: pathlib.Path
            ) -> None:
        p.mkdir(parents=True, exist_ok=True)
        return
    
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
        cmd_template: str   # command with {index}
        attempts: Dict[str, List[int]]  # jobid -> indices attempted
        parent_jobids: List[str]  # history
        #slurm_env: Dict[str, str] # extra env to export in script
        worker_init: Optional[str] = None
        state_file: Optional[str] = None

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
        
    #---------------- slurm inspection ----------------
    SLURM_DONE = {"COMPLETED", "CANCELLED", "FAILED", "TIMEOUT", "PREEMPTED", None}
    SLURM_FAILISH = {"FAILED", "TIMEOUT", "NODE_FAIL", "CANCELLED"}
    DEFAULT_INTERVAL = 0.5

    def parse_sacct_array_states(
            jobid: str
            ) -> Dict[int, str]:
        """
        Prefer sacct if available; returns {index: State}
        """
        if not which("sacct"):
            return {}
        out = run(["bash", "-lc", f"sacct -n -P -j {shlex.quote(jobid)} "\
                   f"--format=JobIDRaw,State"], check=False)
        states: Dict[int,str] = {}
        for line in out.splitlines():
            jid_raw, state = line.split("|", 1)
            if "_" in jid_raw: # array element like 12345_3
                part = jid_raw.split("_", 1)[1]
                idx_s = part.split(".")[0]
                if idx_s.isdigit():
                    states[int(idx_s)] = state.strip()
        return states
    
    def scontrol_job_state(
            jobid: str
            ) -> Optional[str]:
        out = run(["bash", "-lc", f"scontrol -o show job {shlex.quote(jobid)}"], 
                  check=False)
        if not out: return None
        # look for 'JobState=RUNNING' etc.
        kv = dict(field.split("=",1) for field in out.split() if "=" in field)
        return kv.get("JobState")
    
    def scontrol_array_index_state(
            jobid: str,
            idx: int
            ) -> Optional[str]:
        out = run(["bash", "-lc", f"scontrol -o show job {shlex.quote(jobid)}_{idx}"], 
                  check=False)
        if not out: return None
        # look for 'JobState=RUNNING' etc.
        kv = dict(field.split("=",1) for field in out.split() if "=" in field)
        return kv.get("JobState")
    
    def squeue_active_indices(
            jobid: str
            ) -> List[int]:
        out = run(["bash", "-lc", f"squeue -h -o '%i %T %K' -j {shlex.quote(jobid)}"], 
                  check=False)
        idxs = []
        pending = []
        for line in out.splitlines():
            # e.g., "12345_7 RUNNING gpu"
            parts = line.split()
            if not parts: continue
            jid = parts[0]
            if "_" in jid:  # array job
                idx_s = jid.split("_", 1)[1]
                if idx_s.isdigit():
                    idxs.append(int(idx_s))
                elif parts[1] == "PENDING":
                    collapsed = idx_s.strip("[]")
                    ranges = collapsed.split(",")
                    for myrange in ranges:
                        results = myrange.split("-")
                        if len(results) == 1:
                            pending.append(int(results[0]))
                        elif len(results) == 2:
                            start, stop = results
                            pending += list(range(int(start), int(stop)+1))
                        else:
                            raise Exception(
                                f"Incorrect PENDING string for parsing: {idx_s}")
                    pending.append(idx_s)
        return sorted(idxs), pending
    
    def get_array_states(
            jobid: str,
            n_tasks: int,
            prefer_sacct: bool = True
            ) -> Dict[int, str]:
        """
        Return per-index state; tries sacct then falls back to scontrol per index.
        """
        if prefer_sacct:
            st = parse_sacct_array_states(jobid)
            if st:
                return st
        # Fallback: query visible active indices from squeue, mark the rest
        #  via scontrol
        active, pending = squeue_active_indices(jobid)
        active = set(active)
        states: Dict[int,str] = {}
        for i in range(n_tasks):
            if i in active:
                states[i] = "ACTIVE" # RUNNING/PENDING; we'll refine
            elif i in pending:
                states[i] = "PENDING"
            else:
                s = scontrol_array_index_state(jobid, i)
                states[i] = s or "UNKNOWN"
        return states

    def submit_array(
            name: str,
            partition: str,
            scheduler_options: Optional[str],
            time_limit: str,
            account: Optional[str],
            constraint: Optional[str],
            array_spec: str,
            logdir: pathlib.Path,
            workdir: pathlib.Path,
            cpus_per_task: Optional[int],
            mem: Optional[str],
            cmd_template: str,
            worker_init: str = ""
            ) -> str:
        ensure_dir(logdir)
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
        
        wrap_cmd = cmd_template.replace("{index}", "${SLURM_ARRAY_TASK_ID:-0}")
        full = shlex.quote(wrap_cmd) if not worker_init else shlex.quote(
            f"{worker_init}; {wrap_cmd}")
        args += ["--wrap", full]
        out = run(["bash", "-lc", " ".join(x for x in args)])
        return out.strip().split()[-1]
    
    def get_state_dir(workdir:pathlib.Path) -> pathlib.Path:
        return workdir/".slurm_runner"

    def cmd_submit(args):
        workdir = pathlib.Path(args["working_dir"] or os.getcwd()).resolve()
        logdir = (workdir / (args["logdir"] or "logs")).resolve()
        ensure_dir(workdir)
        ensure_dir(logdir)
        array_spec = args["indices"]
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
            #slurm_env=dict(kv.split("=",1) for kv in (args["export"] or []))
            worker_init=args["worker_init"]
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
            #slurm_env=dict(kv.split("=",1) for kv in (args.export or []))
            worker_init=args["worker_init"]
        )
        state_path = state_dir / f"{run_id}.json"
        st.save(state_path)
        # symlink latest
        latest = state_dir / "latest.json"
        try:
            if latest.exists() or latest.is_symlink():
                latest.unlink()
            latest.symlink_to(state_path.name)
        except Exception:
            (state_dir / "latest.json").write_text(json.dumps(asdict(st), indent=2))
        
    def cmd_status(args):
        st = RunState.load(pathlib.Path(args["state"]))
        states = get_array_states(st.jobid, st.n_tasks)
        counts = {}
        for s in states.values():
            counts[s] = counts.get(s, 0) + 1
        return states

    def return_index_to_queue(jobid, array_idx):
        out = run(["bash", "-lc", f"scontrol requeue {shlex.quote(jobid)}_{array_idx}"], 
                  check=True)
        return

    def cmd_attach(
            args: Dict
            ) -> bool:
        """
        Attach or reattach to a job: poll until all indices are in a terminal state.
        """
        # NOTE: Keep! There might be a reason to kill jobs from within the workflow someday
        kill_jobs = False
        interval = float(args.get("interval", DEFAULT_INTERVAL))
        st = RunState.load(pathlib.Path(args["state"]))
        while True:
            states = get_array_states(st.jobid, st.n_tasks)
            n_done = sum(1 for s in states.values() if s in SLURM_DONE \
                         or s == "COMPLETED")
            n_fail = sum(1 for s in states.values() if s in SLURM_FAILISH)
            n_act = sum(1 for s in states.values() if s not in SLURM_DONE)
            if n_done == st.n_tasks or (n_act == 0 and n_done + n_fail == st.n_tasks):
                #print(f"[{time.strftime('%H:%M:%S')}] Done: {n_done}, "
                #      f"Failed: {n_fail}, Active: {n_act}", flush=True)
                break
                
            for key, value in states.items():
                if value == "TIMEOUT" and args["resubmit_on_timeout"]:
                    print(f"Job number {st.jobid}_{key} has timed out. "
                          "Resubmitting.")
                    return_index_to_queue(st.jobid, key)

            if os.path.exists(args["stop_flag"]):
                print(f"Detaching from job {st.jobid}.")
                break
                
            time.sleep(interval)
        
        return kill_jobs

    def cmd_cancel(args):
        """
        Cancel a job.
        """
        # NOTE: Keep! There might be a reason to kill jobs from within the workflow someday
        st = RunState.load(pathlib.Path(args["state"]))
        out = run(["bash", "-lc", f"scancel {shlex.quote(st.jobid)}"], check=False)

    """
    # TODO: to be implemented in version 3.0
    #import seekr2.modules.common_base as seekr2_base
    #import seekr2.run as seekr2_run
    #model = seekr2_base.load_model(model_filename)

    @bash_app
    def run_workflow(inputs):
        cmd_str = inputs[0]
        return cmd_str
    """
    work_dir = pathlib.Path(kwargs["working_dir"])
    state_path = get_state_dir(work_dir) / "latest.json"
    some_active = False
    if state_path.exists():
        status_args = {"state":state_path}
        states = cmd_status(status_args)
        
        for key, value in states.items():
            if key in kwargs["indices"]:
                if value in ["ACTIVE", "PENDING"]:
                    kwargs["indices"].remove(key)
                    some_active = True

    if not some_active:
        kwargs["command_string"] = command_string
        cmd_submit(kwargs)

    attach_args = {"state": state_path, "interval": 0.5,
                   "resubmit_on_timeout": kwargs["resubmit_on_timeout"],
                   "stop_flag": kwargs["stop_flag"]}
    kill_jobs = cmd_attach(attach_args)
    if kill_jobs:
        cancel_args = {"state": state_path}
        cmd_cancel(cancel_args)

    return

def slurm_remote_cancel_workflow(args):
    """
    Cancel a SLURM job remotely from within a Globus compute instance.
    """
    import os
    import json
    import shlex
    import pathlib
    import subprocess
    from dataclasses import dataclass, asdict
    from typing import Dict, List, Optional

    working_dir = args[0]

    def run(
            cmd: List[str], 
            check: bool = True
        ) -> str:
        """
        Run a shell command, and if check is True, then check the return code.
        If there was an error - output the error to standard error.
        """
        out = subprocess.run(cmd, stdout=subprocess.PIPE, 
                             stderr=subprocess.PIPE, text=True)
        if check and out.returncode != 0:
            raise RuntimeError(f"Command failed: ({out.returncode}): "\
                               +f"{' '.join(cmd)}\nSTDERR:\n{out.stderr}")
        return out.stdout.strip()
    
    def ensure_dir(
            p: pathlib.Path
            ) -> None:
        p.mkdir(parents=True, exist_ok=True)
        return

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
        cmd_template: str   # command with {index}
        attempts: Dict[str, List[int]]  # jobid -> indices attempted
        parent_jobids: List[str]  # history
        #slurm_env: Dict[str, str] # extra env to export in script
        worker_init: Optional[str] = None
        state_file: Optional[str] = None

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
        
    def get_state_dir(workdir:pathlib.Path) -> pathlib.Path:
        return workdir/".slurm_runner"
    
    def cmd_cancel(args):
        """
        Cancel a job.
        """
        st = RunState.load(pathlib.Path(args["state"]))
        out = run(["bash", "-lc", f"scancel {shlex.quote(st.jobid)}"], 
                  check=True)

    work_dir = pathlib.Path(working_dir)
    state_path = get_state_dir(work_dir) / "latest.json"
    if os.path.exists(state_path):
        cancel_args = {"state": state_path}
        cmd_cancel(cancel_args)
    return

def slurm_remote_status_workflow(args):
    """
    Get the status of a SLURM job remotely from within a Globus compute 
    instance.
    """
    import os
    import time
    import json
    import pathlib
    import subprocess
    from dataclasses import dataclass, asdict
    from typing import Dict, List, Optional

    working_dir = args[0]

    def run(
            cmd: List[str], 
            check: bool = True
        ) -> str:
        """
        Run a shell command, and if check is True, then check the return code.
        If there was an error - output the error to standard error.
        """
        out = subprocess.run(cmd, stdout=subprocess.PIPE, 
                             stderr=subprocess.PIPE, text=True)
        if check and out.returncode != 0:
            raise RuntimeError(f"Command failed: ({out.returncode}): "\
                               +f"{' '.join(cmd)}\nSTDERR:\n{out.stderr}")

        return out.returncode, out.stdout, out.stderr

    def ensure_dir(
            p: pathlib.Path
            ) -> None:
        p.mkdir(parents=True, exist_ok=True)
        return

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
        cmd_template: str   # command with {index}
        attempts: Dict[str, List[int]]  # jobid -> indices attempted
        parent_jobids: List[str]  # history
        #slurm_env: Dict[str, str] # extra env to export in script
        worker_init: Optional[str] = None
        state_file: Optional[str] = None

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
        
    def get_state_dir(workdir:pathlib.Path) -> pathlib.Path:
        return workdir/".slurm_runner"
    
    def cmd_status(args):
        """
        Get the status of a job.
        """
        st = RunState.load(pathlib.Path(args["state"]))
        cmd = f"squeue -j {st.jobid} -h -o '%i|%P|%j|%u|%T|%M|%l|%D|%C|%R'"
        base = ["bash", "-lc", cmd]
        # Try native JSON first (Slurm ≥21.08 supports --json/--yaml)
        limit = None
        # id, part, name, user, state, elapsed, timelimit, nodes, cpus, reason
        fmt = "%i|%P|%j|%u|%T|%M|%l|%D|%C|%R"
        rc, out, err = run(base)
        if rc != 0:
            return {"ok": False, "stderr": err, "rc": rc, "cmd": " ".join(
                base + ["-h","-o",fmt])}
        rows = []
        for line in out.splitlines():
            parts = line.split("|")
            rows.append({
                "JobID": parts[0], "Partition": parts[1], "Name": parts[2], 
                "User": parts[3], "State": parts[4], "Elapsed": parts[5], 
                "TimeLimit": parts[6], "Nodes": parts[7], "CPUs": parts[8], 
                "Reason": parts[9] if len(parts) > 9 else ""
            })
        if limit: rows = rows[:limit]
        return {"tool": "squeue", "ts": time.time(), "rows": rows}
        

    work_dir = pathlib.Path(working_dir)
    state_path = get_state_dir(work_dir) / "latest.json"
    if os.path.exists(state_path):
        status_args = {"state": state_path}
        result = cmd_status(status_args)
    else:
        err = "State file not found."
        result = {"ok": False, "stderr": err}
    return result