"""
modules/seekr_run.py

Using the seekr API, run stages of the seekr calculation.
"""

import os
import sys
import time
import json
import glob
import typing
import select
import signal
import asyncio
import datetime
import tempfile
import multiprocessing
import termios
import tty
from concurrent.futures import ThreadPoolExecutor

# Use "spawn" instead of fork: required for CUDA-safe child processes and
# to avoid inheriting the parent's asyncio/thread state (which produces
# spurious nonzero exit codes during child interpreter shutdown).
multiprocessing.set_start_method("spawn", force=True)

import attrs
from radical.asyncflow import WorkflowEngine, LocalExecutionBackend #, \
                              #ConcurrentExecutionBackend
#from rhapsody.telemetry import define_event
#from rhapsody.telemetry.events import make_event

from rich.console import Console, Group
from rich.live import Live
from rich.table import Table
from rich.text import Text

import seekr.modules.structures as seekr3_structures
import seekr.modules.scales.base as scales_base
import seekr.status as seekr_status

import seekrflow.modules.base as base
import seekrflow.modules.structures as structures
import seekrflow.modules.transfer.base as transfer_base
import seekrflow.modules.workload_managers.local_multiprocessing as workload_local_mp
import seekrflow.modules.workload_managers.remote as workload_remote
import seekrflow.modules.workload_managers.co_schedule_fusion as co_schedule_fusion
import seekrflow.modules.workload_managers.dispatch_lowering as dispatch_lowering
import seekrflow.modules.workload_managers.remote_stage_lifecycle as remote_stage_lifecycle

_RICH_CONSOLE = Console()

KEYSTROKE_CHECK_INTERVAL = 1.0 # seconds
MAIN_LOOP_INTERVAL = 30.0 # 600.0  # seconds
MAX_SUBSEQUENT_NONCOMPLETED_RUNS = 3
IDLE_INCOMPLETE_CHECKS_BEFORE_RESUBMIT = 3
COMPLETED_DRAIN_CHECKS = 3
SHUTDOWN_CANCEL_TIMEOUT = 30.0
ENGINE_SHUTDOWN_TIMEOUT = 15.0

POLLING_INTERVAL = 5.0  # seconds
STATUS_WRITE_INTERVAL = 5.0  # seconds
STATUS_FILE_NAME = ".seekrflow_job_status.json"
STATUS_SCHEMA_VERSION = 1

combine_fused_commands = co_schedule_fusion.combine_fused_commands
populate_fusion_map = co_schedule_fusion.populate_fusion_map
fusion_dependencies_satisfied = co_schedule_fusion.fusion_dependencies_satisfied
skips_remote_submit = co_schedule_fusion.skips_remote_submit
is_fusion_host = co_schedule_fusion.is_fusion_host
fusion_host_name = co_schedule_fusion.fusion_host_name
mark_fused_set_completed = co_schedule_fusion.mark_fused_set_completed
classify_fused_probe_status = co_schedule_fusion.classify_fused_probe_status
remote_scheduler_job_name = remote_stage_lifecycle.remote_scheduler_job_name
classify_remote_probe_status = remote_stage_lifecycle.classify_remote_probe_status
owns_scheduler_job = remote_stage_lifecycle.owns_scheduler_job
remote_cancel_needed = remote_stage_lifecycle.remote_cancel_needed


def build_remote_stage_command_string(
        seekrflow: structures.Seekrflow,
        stage: scales_base.Base_stage,
        resource: structures.Resource_remote_base,
        resolved_execution: structures.Resolved_execution | None,
        destination_path: str,
        *,
        force_overwrite: bool,
        benchmark_mode: bool,
        ) -> tuple[str, list[int] | None]:
    """
    Build the remote command for one stage (single job or array driver).
    """
    output_filename = os.path.join(
        destination_path, f"{stage.name}_run.out")
    force_overwrite_flag = "True" if force_overwrite else "False"
    benchmark_flag = "True" if benchmark_mode else "False"
    dispatch = (
        resolved_execution.dispatch
        if resolved_execution is not None else None)
    dimensions = dispatch.dimensions if dispatch is not None else None
    group_size = dispatch.group_size if dispatch is not None else None
    concurrency = dispatch.concurrency if dispatch is not None else 1

    if dispatch_lowering.needs_unit_enumeration(dimensions):
        unit_counts = workload_remote.fetch_unit_counts_remote(
            seekrflow, stage, resource)
        lowered = dispatch_lowering.lower_dispatch(
            dimensions,
            group_size,
            concurrency,
            unit_counts,
        )
        command_string = dispatch_lowering.build_remote_command_string(
            stage.name,
            lowered,
            force_overwrite=force_overwrite,
            benchmark=benchmark_mode,
            output_filename=output_filename,
        )
        indices = dispatch_lowering.array_member_indices(lowered)
        return command_string, indices

    command_string = (
        "python -c \"import seekr.modules.structures as structures; "
        "import seekr.run as seekr_run; "
        "model = structures.load_model('model.json'); "
        f"seekr_run.run(model, '{stage.name}', 'any', None, "
        f"{force_overwrite_flag}, None, benchmark={benchmark_flag})\""
        f" > {output_filename}"
    )
    return command_string, None


def load_job_status(root_directory: str) -> dict:
    """
    Read the most recent job status snapshot for a pipeline root directory.
    """
    path = os.path.join(root_directory, STATUS_FILE_NAME)
    with open(path, "r") as f:
        return json.load(f)

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


@attrs.define
class StageWorkflow:
    """
    A workflow for a single stage of SEEKR3.
    """
    # Constructor args (no defaults → required, in this order)
    model: seekr3_structures.Seekr_model = attrs.field(repr=False)
    seekrflow: structures.Seekrflow = attrs.field(repr=False)
    stage: scales_base.Base_stage = attrs.field(repr=False)
    workflow_engine: WorkflowEngine = attrs.field(repr=False)
    resource_name: str = attrs.field(default="local")
    resource: structures.Resource_remote_base | None = attrs.field(
        default=None, repr=False)
    resolved_execution: structures.Resolved_execution | None = attrs.field(
        default=None, repr=False)

    # Derived / mutable state
    dependency_indices: list[int] = attrs.field(factory=list)
    dependency_tasks: list = attrs.field(factory=list)
    task: typing.Any = attrs.field(default=None, repr=False)
    process: multiprocessing.Process | None = attrs.field(
        default=None, repr=False)
    state: str = attrs.field(
        default="unknown", validator=attrs.validators.in_(
            {"unknown", "unstarted", "started", "completed", "failed"}))
    semaphore: str = attrs.field(
        default="go", validator=attrs.validators.in_({"go", "wait", "stop"}))
    progress: float = attrs.field(
        default=0.0, validator=attrs.validators.ge(0.0))
    transfer_from: str | None = attrs.field(default=None)
    manager_status: str = attrs.field(
        default="idle", validator=attrs.validators.in_(
            {"idle", "running", "queued", "running/queued", "unknown"}))
    job_ids: set[str] = attrs.field(factory=set)
    subsequent_noncompleted_runs: int = attrs.field(default=0)
    last_progress: float = attrs.field(default=0.0)
    running_start_time: float | None = attrs.field(default=None)
    quick_failure_warned: bool = attrs.field(default=False)
    force_overwrite: bool = attrs.field(default=False)
    benchmark_mode: bool = attrs.field(default=False)
    detached_requested: bool = attrs.field(default=False)
    transfer_status: str = attrs.field(default="idle")
    transfer_error: str | None = attrs.field(default=None)
    last_raw_status: dict | None = attrs.field(default=None, repr=False)
    co_schedule_with: str | None = attrs.field(default=None)
    fusion_host: str | None = attrs.field(default=None)
    fused_before: list[str] = attrs.field(factory=list)
    fused_after: list[str] = attrs.field(factory=list)
    peer_workflows: dict[str, "StageWorkflow"] = attrs.field(
        factory=dict, repr=False)

    def __attrs_post_init__(self) -> None:
        """
        Normalize resource fields so StageWorkflow is the single source of
        truth for resource resolution.
        """
        if self.resource is not None:
            if self.resource_name == "local":
                self.resource_name = self.resource.name
            elif self.resource.name != self.resource_name:
                raise ValueError(
                    f"Stage {self.stage.name!r} resource mismatch: "
                    f"resource_name={self.resource_name!r}, "
                    f"resource.name={self.resource.name!r}"
                )
            return

        if self.resource_name == "local":
            return

        self.resource = self.seekrflow.run_settings.get_resource_by_name(
            self.resource_name
        )

    def assign_dependency(
            self, 
            dependency_index: int, 
            dependency_task: typing.Any, 
            dependency_resource_name: str
            ) -> None:
        """ Assign a dependency to this stage workflow. """
        self.dependency_indices.append(dependency_index)
        self.dependency_tasks.append(dependency_task)
        if dependency_resource_name != self.resource_name:
            self.transfer_from = dependency_resource_name
        return

    async def create_monitor_only_task(self) -> None:
        """Register only the monitor loop (no transfer or submission)."""
        @self.workflow_engine.function_task
        async def monitor_only(*args):
            await self._monitor_stage_loop()

        self.task = monitor_only()
        await asyncio.sleep(0)

    async def probe_remote_launch(
            self,
            ) -> tuple[str, dict | None]:
        """
        One-shot remote probe for fresh launch from ``unstarted`` / ``unknown``.

        Returns ``(action, status)`` where action is ``completed``,
        ``reattach``, or ``submit``.
        """
        if self.resource_name == "local" or self.fusion_host is not None:
            return "submit", None
        if is_fusion_host(self):
            member_statuses: dict[str, dict | None] = {}
            for stage_name in co_schedule_fusion.fused_set_members(self):
                try:
                    member_statuses[stage_name] = workload_remote.status_remote(
                        self.seekrflow,
                        stage_name,
                        self.model,
                        benchmark_mode=self.benchmark_mode,
                    )
                except Exception as e:
                    print(
                        f"[remote-probe] fused set member {stage_name}: "
                        f"probe failed ({e}); will submit fresh")
                    member_statuses[stage_name] = None
            action = classify_fused_probe_status(
                member_statuses, self.stage.name)
            return action, member_statuses.get(self.stage.name)
        try:
            status = workload_remote.status_remote(
                self.seekrflow,
                self.stage.name,
                self.model,
                benchmark_mode=self.benchmark_mode,
            )
        except Exception as e:
            print(
                f"[remote-probe] stage {self.stage.name}: probe failed "
                f"({e}); submitting fresh")
            return "submit", None
        return classify_remote_probe_status(status), status

    async def create_tasks(self):
        """
        Register tasks that will be python functions and shell command strings.
        """
        if skips_remote_submit(self):
            await self.create_monitor_only_task()
            return

        # Transfer files if dependent stage resource is different
        @self.workflow_engine.function_task
        async def transfer_files(*args):
            if self.transfer_from is None:
                # For the first remote stage (no dependency), seed remote workdir
                # from local so model/config files exist before submission.
                dep_index = getattr(self.stage, "input_stage_index", 0)
                if self.resource_name != "local" and dep_index <= 0:
                    src_resource = None
                    dst_resource = self.resource
                else:
                    self.transfer_status = "skipped"
                    self.transfer_error = None
                    return
            else:
                src_resource = self.seekrflow.run_settings.get_resource_by_name(
                    self.transfer_from)
                dst_resource = self.resource

            local_directory = self.model.directory

            if src_resource is None and dst_resource is None:
                self.transfer_status = "skipped"
                self.transfer_error = None
                return

            try:
                self.transfer_status = "running"
                self.transfer_error = None

                # Conservative resource hop: remote source -> local -> remote destination.
                if src_resource is not None:
                    transfer_base.transfer_files_to_from_remote_resource(
                        self.seekrflow.name,
                        src_resource,
                        local_directory,
                        backwards=True,
                    )
                if dst_resource is not None:
                    transfer_base.transfer_files_to_from_remote_resource(
                        self.seekrflow.name,
                        dst_resource,
                        local_directory,
                        backwards=False,
                    )

                self.transfer_status = "completed"
            except Exception as e:
                self.transfer_status = "failed"
                self.transfer_error = str(e)
                self.semaphore = "wait"
                self.state = "failed"
                print(
                    f"[transfer] stage {self.stage.name} transfer failed; "
                    f"setting semaphore=wait. error={e}"
                )
                raise
            
        dep_index = getattr(self.stage, "input_stage_index", 0)
        should_transfer = bool(self.transfer_from) or (
            self.resource_name != "local" and dep_index <= 0
        )
        if should_transfer:
            self.dependency_tasks.append(transfer_files(*self.dependency_tasks))

        # Run stage
        # TODO: this becomes a combo of executing task submission, then also
        #   monitoring the task and updating progress, state, etc. based 
        #   on telemetry.
        @self.workflow_engine.function_task
        async def run_stage(*args):
            if self.resource_name == "local":
                # Capture the force-overwrite request before clearing it; it is
                # a one-shot flag but must still be propagated to run_locally().
                force_overwrite_now = self.force_overwrite
                if self.force_overwrite:
                    workload_local_mp.kill_existing_local_stage_processes(
                        self.stage.name, self.model.directory)
                    # Force-rerun is a one-shot request.
                    self.force_overwrite = False
                existing_state = workload_local_mp\
                    .check_for_existing_local_processes(
                        self.model.directory, self.stage.name)
                if existing_state and not force_overwrite_now:
                    # Note: We can't truly "reattach" to a multiprocessing.Process object,
                    # but we can track the PID and monitor/kill it via the state file
                    print(f"  Reattached to {self.stage.name} process (PID: {existing_state.pid})")
                    # self.process remains None - we'll use the PID from state file
                    # TODO: does more need to be done here?
                else:
                    self.process = multiprocessing.Process(
                        target=workload_local_mp.run_locally,
                        args=(self.model.directory, self.stage.name,),
                        kwargs={
                            "force_overwrite": force_overwrite_now,
                            "benchmark_mode": self.benchmark_mode,
                        },
                    )
                    self.process.start()
            else:
                if self.resource is None:
                    raise Exception(
                        f"Remote resource config missing for {self.resource_name!r}")
                try:
                    destination_path = os.path.join(
                        self.resource.remote_working_directory, self.seekrflow.name)
                    destination_model_filename = os.path.join(
                        destination_path, "model.json")

                    # Capture and apply each stage's one-shot force_overwrite.
                    # For remote stages this both (a) cancels any running job
                    # and clears runner state via cancel_and_reset_remote_stage and
                    # (b) threads the flag into the seekr run command so seekr
                    # wipes the stage output, mirroring the local path.
                    stages_in_job = (
                        list(self.fused_before)
                        + [self.stage.name]
                        + list(self.fused_after))
                    force_map: dict[str, bool] = {}
                    for stage_name in stages_in_job:
                        target_sw = (
                            self if stage_name == self.stage.name
                            else self.peer_workflows[stage_name])
                        force_now = target_sw.force_overwrite
                        force_map[stage_name] = force_now
                        if force_now:
                            workload_remote.cancel_and_reset_remote_stage(
                                self.seekrflow, stage_name)
                            target_sw.force_overwrite = False

                    stage_commands: list[str] = []
                    indices = None
                    for fused_name in self.fused_before:
                        fused_sw = self.peer_workflows[fused_name]
                        fused_cmd, _ = build_remote_stage_command_string(
                            self.seekrflow,
                            fused_sw.stage,
                            self.resource,
                            fused_sw.resolved_execution,
                            destination_path,
                            force_overwrite=force_map[fused_name],
                            benchmark_mode=fused_sw.benchmark_mode,
                        )
                        stage_commands.append(fused_cmd)

                    host_cmd, host_indices = build_remote_stage_command_string(
                        self.seekrflow,
                        self.stage,
                        self.resource,
                        self.resolved_execution,
                        destination_path,
                        force_overwrite=force_map[self.stage.name],
                        benchmark_mode=self.benchmark_mode,
                    )
                    stage_commands.append(host_cmd)
                    indices = host_indices

                    for fused_name in self.fused_after:
                        fused_sw = self.peer_workflows[fused_name]
                        fused_cmd, _ = build_remote_stage_command_string(
                            self.seekrflow,
                            fused_sw.stage,
                            self.resource,
                            fused_sw.resolved_execution,
                            destination_path,
                            force_overwrite=force_map[fused_name],
                            benchmark_mode=fused_sw.benchmark_mode,
                        )
                        stage_commands.append(fused_cmd)

                    if len(stage_commands) == 1:
                        command_string = host_cmd
                    else:
                        command_string = combine_fused_commands(stage_commands)

                    # Determine effective walltime:
                    #   - benchmark mode: short fixed cap
                    #   - SEEKR stage with prior progress: dynamic optimizer
                    #   - otherwise: resource default (via override=None)
                    time_limit_override = None
                    anchor_times_for_submit = None
                    if self.benchmark_mode:
                        time_limit_override = base.BENCHMARK_REMOTE_TIME_LIMIT
                        print(
                            f"[seekr-time] stage {self.stage.name}: benchmark mode "
                            f"-> requesting {time_limit_override}"
                        )
                    elif self.stage.name == "seekr":
                        job_status_filename = os.path.join(
                            self.model.directory, STATUS_FILE_NAME)
                        incomplete_anchors: list[int] = []
                        try:
                            with open(job_status_filename, "r") as _f:
                                _status_data = json.load(_f)
                            _stage_status = (
                                _status_data.get(self.stage.name, {})
                                .get("stage_status", {})
                            )
                            incomplete_anchors = list(
                                _stage_status.get("incomplete_anchors") or []
                            )
                            anchor_times_for_submit = (
                                _stage_status.get("anchor_times") or None
                            )
                        except Exception:
                            pass
                        if incomplete_anchors:
                            optimized = workload_remote.calculate_optimal_time_limit(
                                self.seekrflow,
                                self.stage.name,
                                self.resource,
                                incomplete_anchors,
                                job_status_filename,
                            )
                            if optimized and optimized != self.resource.time_limit:
                                time_limit_override = optimized
                                print(
                                    f"[seekr-time] stage {self.stage.name}: "
                                    f"requesting {optimized} (cap: "
                                    f"{self.resource.time_limit})"
                                )

                    run_result = workload_remote.submit_remote_run_workflow(
                        self.seekrflow,
                        self.stage.name,
                        destination_path,
                        self.resource,
                        command_string,
                        destination_model_filename,
                        workflow_type=self.stage.name,
                        indices=indices,
                        anchor_times=anchor_times_for_submit,
                        time_limit_override=time_limit_override,
                    )
                    self.last_raw_status = run_result

                    if not isinstance(run_result, dict):
                        raise RuntimeError(
                            f"Remote submit returned invalid payload: {run_result!r}"
                        )
                    if run_result.get("success") is False:
                        raise RuntimeError(
                            f"Remote submit failed for stage {self.stage.name}: "
                            f"{run_result.get('error')}"
                        )
                    if "success" not in run_result:
                        raise RuntimeError(
                            f"Remote submit returned payload without success flag: {run_result!r}"
                        )

                    self.state = "started"
                    self.manager_status = "running"
                    job_id = None
                    if run_result is not None:
                        job_id = run_result.get("job_id", run_result.get("jobid"))
                    if job_id is not None:
                        self.job_ids.add(str(job_id))
                    else:
                        raise RuntimeError(
                            f"Remote submit did not return a job id: {run_result!r}"
                        )
                except Exception as e:
                    self.last_raw_status = {"success": False, "error": str(e)}
                    self.state = "failed"
                    self.manager_status = "idle"
                    if self.semaphore != "wait":
                        self.semaphore = "wait"
                    print(
                        f"[remote-submit] stage {self.stage.name} submit failed; "
                        f"setting semaphore=wait. error={e}"
                    )
                    raise
            return

        self.task = run_stage(*self.dependency_tasks)
        await asyncio.sleep(0) # let run_stage's async_wrapper register

        @self.workflow_engine.function_task
        async def monitor_stage(task):
            await self._monitor_stage_loop()

        self.task = monitor_stage(self.task)
        await asyncio.sleep(0) # let monitor_stage's async_wrapper register
        return

    async def _monitor_stage_loop(self) -> None:
        idle_incomplete_checks = 0
        completed_drain_checks = 0
        remote_scheduler_active = False
        while True:
            if self.detached_requested:
                break
            if self.semaphore == "stop":
                break
            fused_host_name = None
            if self.resource_name == "local":
                process_alive = (
                    self.process is not None and self.process.is_alive())
                try:
                    status = workload_local_mp.status_local(
                        self.model.directory, self.stage.index,
                        self.stage.name, self.process)
                except Exception as e:
                    if process_alive:
                        print(f"[local-status] stage {self.stage.name}: "
                              f"status check raised (process still alive, "
                              f"will retry): {e}")
                        await asyncio.sleep(POLLING_INTERVAL)
                        continue
                    print(f"[local-status] stage {self.stage.name}: status "
                          f"check raised and process not alive: {e}")
                    status = None
                if status is not None:
                    self.last_raw_status = status
                try:
                    workload_local_mp.check_and_raise_if_process_failed(
                        self.stage.name, self.process,
                        self.model.directory)
                except Exception as e:
                    self.state = "failed"
                    self.manager_status = "idle"
                    print(f"[local-run] stage {self.stage.name} failed:\n{e}")
                if status is not None and self.state != "failed":
                    self.state = status["stage_status"]["state"]
                    self.progress = status["stage_status"]["progress"]
                    manager_status = status["manager_status"]
                    self.manager_status = determine_stage_manager_status(
                        manager_status)
                    if manager_status and manager_status.get("jobs"):
                        for job in manager_status["jobs"]:
                            self.job_ids.add(job["JobID"])
            else:
                fused_host_name = fusion_host_name(self)
                host_sw = None
                if self.fusion_host is not None:
                    host_sw = self.peer_workflows.get(self.fusion_host)
                fused_host_name = fusion_host_name(self)

                status = None
                try:
                    status = workload_remote.status_remote(
                        self.seekrflow,
                        self.stage.name,
                        self.model,
                        benchmark_mode=self.benchmark_mode,
                    )
                except Exception as e:
                    print(
                        f"[remote-status] stage {self.stage.name}: "
                        f"status call raised an exception: {e}"
                    )
                    await asyncio.sleep(POLLING_INTERVAL)
                    continue
                if status is not None:
                    self.last_raw_status = status
                    manager_status = status.get("manager_status")
                    self.manager_status = determine_stage_manager_status(
                        manager_status)
                    if manager_status and manager_status.get("jobs"):
                        for job in manager_status["jobs"]:
                            job_id = job.get("JobID")
                            if job_id:
                                self.job_ids.add(str(job_id))
                    host_active = False
                    if host_sw is not None:
                        host_active = host_sw.manager_status in {
                            "running", "queued", "running/queued",
                        }
                        if host_sw.job_ids:
                            self.job_ids.update(host_sw.job_ids)
                    status_ok = status.get("success", True)
                    stage_status = status.get("stage_status", {})
                    scheduler_active = self.manager_status in {
                        "running", "queued", "running/queued",
                    }
                    if host_sw is not None:
                        scheduler_active = scheduler_active or host_active
                    remote_scheduler_active = scheduler_active
                    if not status_ok:
                        remote_err = (
                            status.get("error")
                            or stage_status.get("notes", ""))
                        notes_msg = stage_status.get("notes", "")
                        display_err = notes_msg if notes_msg else remote_err
                        if scheduler_active:
                            print(
                                f"[remote-status] stage {self.stage.name}: "
                                f"progress check failed (job still in "
                                f"scheduler): {display_err}"
                            )
                        else:
                            print(
                                f"[remote-status] stage {self.stage.name}: "
                                f"progress check failed and no active jobs: "
                                f"{display_err}"
                            )
                            self.state = "failed"
                    else:
                        self.state = stage_status.get("state", self.state)
                        self.progress = stage_status.get(
                            "progress", self.progress)
                        # Only submit-capable remote stages may resubmit.
                        if self.fusion_host is None:
                            if fused_host_name is not None:
                                set_incomplete = not (
                                    co_schedule_fusion.fused_set_completed(
                                        fused_host_name,
                                        self.peer_workflows))
                                progress_for_resubmit = (
                                    co_schedule_fusion.fused_set_progress(
                                        fused_host_name,
                                        self.peer_workflows))
                            else:
                                set_incomplete = self.state != "completed"
                                progress_for_resubmit = self.progress
                            if set_incomplete and not scheduler_active:
                                idle_incomplete_checks += 1
                                if (idle_incomplete_checks
                                        >= IDLE_INCOMPLETE_CHECKS_BEFORE_RESUBMIT):
                                    if (progress_for_resubmit
                                            > self.last_progress):
                                        self.last_progress = (
                                            progress_for_resubmit)
                                        self.subsequent_noncompleted_runs = 0
                                        print(
                                            f"[remote-resubmit] stage "
                                            f"{self.stage.name}: job gone "
                                            f"with progress advance "
                                            f"({progress_for_resubmit:.1%}); "
                                            f"resubmitting")
                                        self.state = "unstarted"
                                        break
                                    self.subsequent_noncompleted_runs += 1
                                    if (self.subsequent_noncompleted_runs
                                            > MAX_SUBSEQUENT_NONCOMPLETED_RUNS):
                                        print(
                                            f"[remote-resubmit] stage "
                                            f"{self.stage.name}: no progress "
                                            f"after "
                                            f"{self.subsequent_noncompleted_runs} "
                                            f"idle resubmits; giving up")
                                        self.state = "failed"
                                        self.semaphore = "wait"
                                        break
                                    print(
                                        f"[remote-resubmit] stage "
                                        f"{self.stage.name}: job gone without "
                                        f"progress "
                                        f"({self.subsequent_noncompleted_runs}/"
                                        f"{MAX_SUBSEQUENT_NONCOMPLETED_RUNS}); "
                                        f"resubmitting")
                                    self.state = "unstarted"
                                    break
                            else:
                                idle_incomplete_checks = 0

            if self.resource_name == "local":
                terminal_completed = self.state == "completed"
            elif fused_host_name is not None:
                terminal_completed = co_schedule_fusion.fused_set_completed(
                    fused_host_name, self.peer_workflows)
            else:
                terminal_completed = self.state == "completed"

            if terminal_completed:
                if (self.resource_name != "local"
                        and remote_scheduler_active):
                    completed_drain_checks += 1
                    if completed_drain_checks >= COMPLETED_DRAIN_CHECKS:
                        print(
                            f"[remote-drain] stage {self.stage.name}: "
                            f"scheduler still active after "
                            f"{completed_drain_checks} completion checks; "
                            f"forcing manager=idle")
                        self.manager_status = "idle"
                        break
                    await asyncio.sleep(POLLING_INTERVAL)
                    continue
                if self.resource_name != "local":
                    self.manager_status = "idle"
                break
            elif self.state == "failed":
                if self.semaphore != "wait":
                    self.semaphore = "wait"
                    print(f"[semaphore] stage {self.stage.name} failed; "
                          f"setting semaphore=wait")
                break
            else:
                await asyncio.sleep(POLLING_INTERVAL)
        self.task = None
        return

    def status_snapshot(self) -> dict:
        """
        Return a JSON-serializable snapshot of this stage's current state
        for the pipeline-level status file.
        """
        process_info = None
        if self.process is not None:
            process_info = {
                "pid": self.process.pid,
                "alive": self.process.is_alive(),
                "exitcode": self.process.exitcode,
            }
        return {
            "index": self.stage.index,
            "name": self.stage.name,
            "resource_name": self.resource_name,
            "state": self.state,
            "semaphore": self.semaphore,
            "progress": self.progress,
            "transfer_status": self.transfer_status,
            "transfer_error": self.transfer_error,
            "last_progress": self.last_progress,
            "manager_status": self.manager_status,
            "job_ids": sorted(self.job_ids),
            "force_overwrite": self.force_overwrite,
            "transfer_from": self.transfer_from,
            "benchmark_mode": self.benchmark_mode,
            "subsequent_noncompleted_runs":
                self.subsequent_noncompleted_runs,
            "running_start_time": self.running_start_time,
            "quick_failure_warned": self.quick_failure_warned,
            "process": process_info,
            "raw_status": self.last_raw_status,
        }

    def summary_line(self) -> str:
        """
        Short one-line human summary for the 'i' keystroke command.
        Format: name | state | progress | manager_status | resource | pid
        """
        pid = self.process.pid if self.process is not None else "-"
        try:
            progress_str = f"{self.progress:.2%}"
        except (TypeError, ValueError):
            progress_str = str(self.progress)
        return (f"{self.stage.name:<10}  state={self.state:<10}  "
                f"progress={progress_str:>7}  manager={self.manager_status:<14}  "
                f"resource={self.resource_name:<10}  pid={pid}")

    def kill(self) -> None:
        """
        Stop this stage's monitor loop and terminate any running job.
        Idempotent; safe to call on stages that never started.
        """
        self.semaphore = "stop"
        if self.resource_name != "local":
            if self.resource is None:
                return
            if not owns_scheduler_job(self.fusion_host):
                self.manager_status = "idle"
                return
            if (is_fusion_host(self)
                    and co_schedule_fusion.fused_set_completed(
                        self.stage.name, self.peer_workflows)):
                self.manager_status = "idle"
                return
            try:
                status = workload_remote.status_remote(
                    self.seekrflow,
                    self.stage.name,
                    self.model,
                    benchmark_mode=self.benchmark_mode,
                )
            except Exception as e:
                print(f"Warning: failed to query remote jobs for cancel: {e}")
                status = None
            if not remote_cancel_needed(
                    status,
                    local_state=self.state,
                    tracked_job_ids=self.job_ids):
                print(
                    f"[shutdown] skipping cancel for {self.stage.name}: "
                    f"no active scheduler jobs")
                self.manager_status = "idle"
                return
            manager_status = (status or {}).get("manager_status") or {}
            for job in manager_status.get("jobs", []):
                job_id = job.get("JobID")
                if job_id:
                    self.job_ids.add(str(job_id))
            job_name = remote_scheduler_job_name(
                self.seekrflow.name, self.stage.name)
            print(
                f"[shutdown] cancelling remote jobs for {self.stage.name}: "
                f"ids={sorted(self.job_ids)}, name={job_name!r}")
            for job_id in sorted(self.job_ids):
                try:
                    workload_remote.submit_remote_cancel_workflow(
                        self.seekrflow,
                        self.resource_name,
                        job_id=job_id,
                    )
                except Exception as e:
                    print(f"Warning: failed to cancel remote job {job_id}: {e}")
            try:
                workload_remote.submit_remote_cancel_workflow(
                    self.seekrflow,
                    self.resource_name,
                    job_name=job_name,
                )
            except Exception as e:
                print(
                    f"Warning: failed to cancel remote jobs by name "
                    f"{job_name!r}: {e}")
            self.manager_status = "idle"
            return
        process = self.process
        if process is None or not process.is_alive():
            return
        pid = process.pid
        print(f"Terminating {self.stage.name} process group (PID: {pid})...")
        try:
            os.killpg(pid, signal.SIGTERM)
        except ProcessLookupError:
            return
        process.join(timeout=5)
        if process.is_alive():
            print(f"Force killing {self.stage.name} process group (PID: {pid})...")
            try:
                os.killpg(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.join()
        return


def _semaphore_text(semaphore: str) -> Text:
    colors = {"go": "green", "wait": "yellow", "stop": "red"}
    return Text(semaphore, style=colors.get(semaphore, "white"))


def _state_text(state: str) -> Text:
    colors = {
        "completed": "green",
        "failed": "bold red",
        "started": "yellow",
        "transferring": "cyan",
        "unstarted": "dim",
        "unknown": "dim",
    }
    return Text(state, style=colors.get(state, "white"))


def _manager_text(manager: str) -> Text:
    colors = {
        "running": "green",
        "queued": "cyan",
        "running/queued": "blue",
        "transfer": "cyan",
        "idle": "dim",
        "unknown": "dim",
    }
    return Text(manager, style=colors.get(manager, "white"))


def _build_stage_table(stage_workflows: list) -> Table:
    table = Table(show_header=True, header_style="bold cyan")
    table.add_column("Stage", style="bold")
    table.add_column("State")
    table.add_column("Progress", justify="right")
    table.add_column("Semaphore")
    table.add_column("Manager")
    table.add_column("Resource")
    table.add_column("PID", justify="right")
    for sw in stage_workflows:
        pid_str = str(sw.process.pid) if sw.process is not None else "-"
        try:
            progress_str = f"{sw.progress:.1%}"
        except (TypeError, ValueError):
            progress_str = str(sw.progress)

        display_state = sw.state
        display_manager = sw.manager_status
        if sw.transfer_status == "running":
            display_state = "transferring"
            display_manager = "transfer"

        table.add_row(
            sw.stage.name,
            _state_text(display_state),
            progress_str,
            _semaphore_text(sw.semaphore),
            _manager_text(display_manager),
            sw.resource_name,
            pid_str,
        )
    return table


@attrs.define
class SeekrPipeline:
    """
    The entire set of workflows for a seekr run.
    """
    # Constructor args (no defaults → required, in this order)
    model: seekr3_structures.Seekr_model = attrs.field(repr=False)
    seekrflow: structures.Seekrflow = attrs.field(repr=False)
    workflow_engine: WorkflowEngine = attrs.field(repr=False)
    
    # Derived / mutable state
    stage_names: list[str] = attrs.field(factory=list)
    stage_workflows: list[StageWorkflow] \
        = attrs.field(factory=list, repr=False)
    stage_tasks: list = attrs.field(factory=list, repr=False)
    task_id_to_stage: dict = attrs.field(factory=dict, repr=False)
    telemetry: typing.Any = attrs.field(default=None, repr=False)
    force_rerun_stages: set[str] = attrs.field(factory=set)
    semaphore_overrides: dict[str, str] = attrs.field(factory=dict)
    benchmark_stage: str | None = attrs.field(default=None)
    keystrokes_enabled: bool = attrs.field(default=True)
    _input_buffer: str = attrs.field(default="", repr=False)
    _live_display: typing.Any = attrs.field(default=None, repr=False)
    _detached_requested: bool = attrs.field(default=False, repr=False)
    _shutting_down: bool = attrs.field(default=False, repr=False)
    _shutdown_task: asyncio.Task | None = attrs.field(
        default=None, repr=False)
    _stop_event: asyncio.Event | None = attrs.field(
        default=None, repr=False)

    def __attrs_post_init__(
            self, 
            ) -> None:
        for stage in self.model.stages:
            self.add_stage(stage)
        by_name = {sw.stage.name: sw for sw in self.stage_workflows}
        for sw in self.stage_workflows:
            sw.peer_workflows = by_name
        populate_fusion_map(self.stage_workflows, self.benchmark_stage)
        for sw in self.stage_workflows:
            if sw.stage.name in self.semaphore_overrides:
                sw.semaphore = self.semaphore_overrides[sw.stage.name]

    def add_stage(
            self, 
            stage: scales_base.Base_stage
            ) -> None:
        """
        Add a stage to the pipeline.
        """
        self.stage_names.append(stage.name)
        try:
            resolved = self.seekrflow.run_settings.resolve_stage_execution(
                stage.name, self.seekrflow.workflow.procedure)
        except ValueError:
            resolved = None
        resource_name = "local" if resolved is None else resolved.resource_name
        stage_workflow = StageWorkflow(
            self.model,
            self.seekrflow,
            stage,
            self.workflow_engine,
            resource_name,
        )
        if resolved is not None:
            stage_workflow.resolved_execution = resolved
            if resolved.resource is not None:
                stage_workflow.resource = resolved.resource
        stage_workflow.force_overwrite = stage.name in self.force_rerun_stages
        stage_workflow.benchmark_mode = (stage.name == self.benchmark_stage)
        self.stage_workflows.append(stage_workflow)

    """
    async def telemetry_setup(self):
        ""
        Set up telemetry for the pipeline.
        ""
        # NOTE: make sure this telemetry watches the job 
        #  'monitors' not the job 'launcher'.
        self.telemetry = await self.flow.start_telemetry(resource_poll_interval=5.0)

        def on_event(event):
            et = event.event_type
            stage_workflow = self.task_id_to_stage.get(event.task_id)
            if et == "asyncflow.TaskResolved":
                self.task_id_time_resolved[event.task_id] = event.event_time
            elif et == "TaskStarted":
                dep_wait = event.event_time - self.task_id_time_resolved.get(event.task_id, event.event_time)
                print(f"[telemetry] TaskStarted   {event.task_id}  dep_wait={dep_wait*1000:.0f} ms")
                if stage_workflow:
                    stage_workflow.state = "started"
            elif et == "TaskCompleted":
                print(f"[telemetry] TaskCompleted {event.task_id}  duration={event.duration_seconds*1000:.0f} ms")
                if stage_workflow:
                    stage_workflow.state = "completed"
            elif et == "TaskFailed":
                print(f"[telemetry] TaskFailed    {event.task_id}  error={event.error_type}")
                if stage_workflow:
                    stage_workflow.state = "failed"
            elif et == "TaskCanceled":
                print(f"[telemetry] TaskCanceled  {event.task_id}")
                if stage_workflow:
                    stage_workflow.state = "canceled"
            
        self.telemetry.subscribe(on_event)
    """

    def _on_signal(self) -> None:
        """ Signal-safe callback: schedule async shutdown. """
        if self._shutting_down:
            print("[shutdown] forcing immediate exit")
            os._exit(1)
        if self._shutdown_task is not None and not self._shutdown_task.done():
            print("[shutdown] forcing immediate exit")
            os._exit(1)
        print("\n\n=== Received interrupt signal ===")
        self._shutdown_task = asyncio.create_task(self.shutdown())

    def _cancel_all_stages(self) -> None:
        """Synchronous remote/local cancellation for all stages."""
        for stage_workflow in self.stage_workflows:
            stage_workflow.kill()

    async def shutdown(self) -> None:
        """ Kill all stages and shut down the workflow engine. Idempotent. """
        if self._shutting_down:
            return
        self._shutting_down = True
        timed_out = False
        loop = asyncio.get_running_loop()
        try:
            await asyncio.wait_for(
                loop.run_in_executor(None, self._cancel_all_stages),
                timeout=SHUTDOWN_CANCEL_TIMEOUT,
            )
        except asyncio.TimeoutError:
            timed_out = True
            print(
                f"[shutdown] remote job cancellation timed out after "
                f"{SHUTDOWN_CANCEL_TIMEOUT}s (endpoint unreachable?); "
                f"abandoning remote cancels and exiting")
        # Wake up background loops (status writer, keystroke watcher) so they
        # can exit promptly instead of waiting for the next tick.
        if self._stop_event is not None:
            self._stop_event.set()
        # Capture terminal state before tearing down the engine.
        self.write_status_snapshot()
        try:
            await asyncio.wait_for(
                self.workflow_engine.shutdown(),
                timeout=ENGINE_SHUTDOWN_TIMEOUT,
            )
        except asyncio.TimeoutError:
            timed_out = True
            print(
                f"[shutdown] workflow engine shutdown timed out after "
                f"{ENGINE_SHUTDOWN_TIMEOUT}s; forcing exit")
        if timed_out:
            os._exit(1)
        return

    async def detach_shutdown(self) -> None:
        """
        Soft shutdown for detach: stop local monitoring/UI loops and close the
        workflow engine, but do NOT kill running stage processes.
        """
        if self._shutting_down:
            return
        self._shutting_down = True
        if self._stop_event is not None:
            self._stop_event.set()
        self.write_status_snapshot()
        await self.workflow_engine.shutdown()
        return

    def status_snapshot(self) -> dict:
        """
        Build the full pipeline status snapshot (pipeline-level metadata
        plus per-stage snapshots).
        """
        root = self.model.directory
        try:
            local_state_dir = workload_local_mp.get_local_state_dir(root)
        except Exception:
            local_state_dir = None
        now = time.time()
        return {
            "schema_version": STATUS_SCHEMA_VERSION,
            "timestamp": now,
            "timestamp_iso": datetime.datetime.fromtimestamp(
                now, tz=datetime.timezone.utc).isoformat(),
            "pid": os.getpid(),
            "work_directory": self.seekrflow.work_directory,
            "root_directory": root,
            "local_state_directory": local_state_dir,
            "status_file": os.path.join(root, STATUS_FILE_NAME),
            "shutting_down": self._shutting_down,
            "detached_requested": self._detached_requested,
            "force_rerun_stages": sorted(self.force_rerun_stages),
            "benchmark_stage": self.benchmark_stage,
            "stages": {
                sw.stage.name: sw.status_snapshot()
                for sw in self.stage_workflows
            },
        }

    def write_status_snapshot(self) -> None:
        """
        Atomically write the pipeline status snapshot to disk.
        Best-effort: never raises; logs and continues on failure.
        """
        root = self.model.directory
        target_path = os.path.join(root, STATUS_FILE_NAME)
        tmp_path = None
        try:
            snapshot = self.status_snapshot()
            os.makedirs(root, exist_ok=True)
            fd, tmp_path = tempfile.mkstemp(
                dir=root, prefix=".seekrflow_job_status.", suffix=".tmp")
            with os.fdopen(fd, "w") as f:
                json.dump(snapshot, f, indent=4, default=str)
            os.replace(tmp_path, target_path)
            tmp_path = None
        except Exception as e:
            print(f"[warning] failed to write job status file: {e}")
        finally:
            if tmp_path is not None and os.path.exists(tmp_path):
                try:
                    os.unlink(tmp_path)
                except OSError:
                    pass

    async def _status_writer_loop(
            self, stop_event: asyncio.Event) -> None:
        """
        Periodically write the status snapshot until stop_event is set.
        """
        while not stop_event.is_set():
            self.write_status_snapshot()
            try:
                await asyncio.wait_for(
                    stop_event.wait(), STATUS_WRITE_INTERVAL)
            except asyncio.TimeoutError:
                pass
        return

    # ------------------------------------------------------------------
    # Keystroke command handling
    # ------------------------------------------------------------------

    def _instruction_str(self) -> str:
        return (
            "Keystroke commands (one per line, then Enter):\n"
            "    s [STAGE]  - (s)top semaphore for STAGE or all\n"
            "    g [STAGE]  - (g)o semaphore for STAGE or all\n"
            "    w [STAGE]  - (w)ait semaphore for STAGE or all\n"
            "    t [STAGE]  - (t)ransfer remote files back for STAGE or all (not yet implemented)\n"
            "    q          - (q)uit / detach (leave running jobs active)\n"
            "    h | ?      - (h)elp: reprint this message"
        )

    def _resolve_stage_arg(
            self, stage_arg: str | None) -> list[StageWorkflow] | None:
        """
        Resolve a keystroke STAGE argument to a list of stage workflows.
        None or "all" → every stage. Unknown name → None (caller prints error).
        """
        if stage_arg is None or stage_arg == "all":
            return list(self.stage_workflows)
        for sw in self.stage_workflows:
            if sw.stage.name == stage_arg:
                return [sw]
        return None

    def request_transfer(self, stage_arg: str | None) -> None:
        """ 't' command: transfer files from remote resource(s) to local. """
        targets = self._resolve_stage_arg(stage_arg)
        if targets is None:
            print(f"[keystroke] unknown stage: {stage_arg!r}")
            return

        root_directory = self.model.directory
        transferred_resources: set[str] = set()

        for sw in targets:
            if sw.resource_name == "local" or sw.resource is None:
                sw.transfer_status = "skipped"
                sw.transfer_error = None
                continue

            resource_name = sw.resource.name
            if resource_name in transferred_resources:
                sw.transfer_status = "completed"
                sw.transfer_error = None
                continue

            try:
                sw.transfer_status = "running"
                sw.transfer_error = None
                transfer_base.transfer_files_to_from_remote_resource(
                    self.seekrflow.name,
                    sw.resource,
                    root_directory,
                    backwards=True,
                )
                transferred_resources.add(resource_name)
                sw.transfer_status = "completed"
                print(
                    f"[keystroke] transferred from remote resource "
                    f"{resource_name} (stage={sw.stage.name})"
                )
            except Exception as e:
                sw.transfer_status = "failed"
                sw.transfer_error = str(e)
                sw.semaphore = "wait"
                print(
                    f"[keystroke] transfer failed for stage {sw.stage.name}; "
                    f"setting semaphore=wait. error={e}"
                )

        if len(transferred_resources) == 0:
            print("[keystroke] no non-local resources selected for transfer.")

    def set_semaphore(
            self, stage_arg: str | None, value: str) -> None:
        """ Set semaphore value for one stage or all stages. """
        targets = self._resolve_stage_arg(stage_arg)
        if targets is None:
            print(f"[keystroke] unknown stage: {stage_arg!r}")
            return
        for sw in targets:
            sw.semaphore = value
            print(f"[keystroke] set {sw.stage.name} semaphore={value}")
            if value == "stop":
                sw.kill()

    def detach(self) -> None:
        """ 'q' command: soft-detach without killing running jobs. """
        if self._detached_requested:
            print("[detach] detach already requested.")
            return
        self._detached_requested = True
        for stage_workflow in self.stage_workflows:
            stage_workflow.detached_requested = True
        self._disown_local_processes_for_detach()
        print("[detach] detaching now; running jobs will continue in background.")
        print("[detach] run flow.py again to reattach monitoring.")
        if self._stop_event is not None:
            self._stop_event.set()

    def _disown_local_processes_for_detach(self) -> None:
        """
        Remove active local multiprocessing children from Python's 
        internal child registry so interpreter shutdown doesn't 
        block on implicit joins.

        This is required for true detach semantics: leave stage processes
        running after the parent flow.py process exits.
        """
        try:
            child_registry = multiprocessing.process._children
        except Exception:
            child_registry = None

        detached = []
        for stage_workflow in self.stage_workflows:
            process = stage_workflow.process
            if process is None:
                continue
            if process.is_alive():
                detached.append((stage_workflow.stage.name, process.pid))
            if child_registry is not None:
                try:
                    child_registry.discard(process)
                except Exception:
                    pass
            # After detach, rely on state files for reattachment and status.
            stage_workflow.process = None

        if detached:
            details = ", ".join(
                f"{stage_name}(pid={pid})" for stage_name, pid in detached)
            print(f"[detach] disowned local process handles: {details}")

    def _make_live_renderable(self) -> Group:
        hint = Text(
            "  Commands: s  g  w  t  q  h  [+ STAGE]   (h = help)",
            style="dim italic",
        )
        return Group(
            _build_stage_table(self.stage_workflows),
            hint,
            Text(f"> {self._input_buffer}", style="bold"),
            Text(""),
        )

    def _refresh_live_display(self) -> None:
        live = self._live_display
        if live is None:
            return
        try:
            live.update(self._make_live_renderable())
            live.refresh()
        except Exception:
            # Display refresh is best-effort; command handling should continue.
            pass

    def _handle_command(self, line: str) -> None:
        """ Parse a single keystroke command line and dispatch. """
        line = line.strip()
        if not line:
            return
        parts = line.split()
        cmd = parts[0].lower()
        stage_arg = parts[1] if len(parts) > 1 else None
        if cmd == "t":
            self.request_transfer(stage_arg)
        elif cmd == "s":
            self.set_semaphore(stage_arg, "stop")
        elif cmd == "g":
            self.set_semaphore(stage_arg, "go")
        elif cmd == "w":
            self.set_semaphore(stage_arg, "wait")
        elif cmd == "q":
            self.detach()
        elif cmd in ("h", "?"):
            print(self._instruction_str())
        else:
            print(f"[keystroke] unknown command: {cmd!r}. Type 'h' for help.")

    async def _keystroke_loop(self, stop_event: asyncio.Event) -> None:
        """
        Watch stdin for character-buffered keystroke commands until stop_event
        is set or stdin reaches EOF. No-op if stdin isn't a TTY or the
        platform doesn't support add_reader on stdin.
        """
        if not self.keystrokes_enabled:
            return
        try:
            if not sys.stdin.isatty():
                return
        except Exception:
            return
        loop = asyncio.get_running_loop()
        fd = sys.stdin.fileno()
        old_tty_attrs = None

        # Use cbreak mode so we can track partially typed input and keep it
        # visible in the live footer across table redraws.
        try:
            old_tty_attrs = termios.tcgetattr(fd)
            tty.setcbreak(fd)
        except Exception as e:
            print(f"[keystroke] failed to configure terminal input mode: {e}")
            return

        def on_stdin() -> None:
            try:
                data = os.read(fd, 1024)
            except Exception as e:
                print(f"[keystroke] read error: {e}")
                try:
                    loop.remove_reader(fd)
                except Exception:
                    pass
                return
            if not data:  # EOF
                try:
                    loop.remove_reader(fd)
                except Exception:
                    pass
                return
            try:
                text = data.decode(errors="ignore")
                dirty = False
                for char in text:
                    if char in ("\n", "\r"):
                        line = self._input_buffer
                        self._input_buffer = ""
                        dirty = True
                        if line.strip():
                            self._handle_command(line)
                    elif char in ("\x7f", "\b"):
                        self._input_buffer = self._input_buffer[:-1]
                        dirty = True
                    elif char == "\x03":
                        # Ctrl-C: mirror signal behavior and clear pending input.
                        self._input_buffer = ""
                        dirty = True
                        self._on_signal()
                    elif char.isprintable():
                        self._input_buffer += char
                        dirty = True
                if dirty:
                    self._refresh_live_display()
            except Exception as e:
                print(f"[keystroke] command error: {e}")

        try:
            loop.add_reader(fd, on_stdin)
        except (NotImplementedError, OSError) as e:
            print(f"[keystroke] stdin watcher unavailable: {e}")
            return
        print(self._instruction_str(), flush=True)
        try:
            await stop_event.wait()
        finally:
            try:
                loop.remove_reader(fd)
            except Exception:
                pass
            if old_tty_attrs is not None:
                try:
                    termios.tcsetattr(fd, termios.TCSADRAIN, old_tty_attrs)
                except Exception:
                    pass
        return

    async def _live_display_loop(self, stop_event: asyncio.Event) -> None:
        """
        Continuously refresh the stage-status table in-place until
        stop_event is set. No-op when stdout is not a TTY.

        Uses auto_refresh=False so redraws only happen on our schedule
        (every POLLING_INTERVAL seconds).  A hint line and a blank line are
        appended below the table so the terminal cursor parks there after
        each redraw, giving the user a visible input zone that is not
        overwritten between refreshes.
        """
        try:
            if not sys.stdout.isatty():
                return
        except Exception:
            return

        with Live(
            self._make_live_renderable(),
            console=_RICH_CONSOLE,
            auto_refresh=False,
            redirect_stdout=True,
            redirect_stderr=True,
        ) as live:
            self._live_display = live
            while not stop_event.is_set():
                self._refresh_live_display()
                try:
                    await asyncio.wait_for(
                        stop_event.wait(), POLLING_INTERVAL)
                except asyncio.TimeoutError:
                    pass
            self._refresh_live_display()
        self._live_display = None
        return

    async def run_workflows(self):
        """
        Run the workflow for the entire pipeline.
        """
        loop = asyncio.get_running_loop()
        for sig in (signal.SIGINT, signal.SIGTERM):
            loop.add_signal_handler(sig, self._on_signal)
        try:
            await self._run_workflows_body()
        finally:
            for sig in (signal.SIGINT, signal.SIGTERM):
                loop.remove_signal_handler(sig)

    async def _run_workflows_body(self):
        stage_by_name = {sw.stage.name: sw for sw in self.stage_workflows}
        launched_tasks: list[typing.Any] = []

        def stage_dependencies_completed(stage_workflow: StageWorkflow) -> bool:
            dep_index = getattr(stage_workflow.stage, "input_stage_index", 0)
            if dep_index <= 0:
                return True
            dep_stage = self.model.stages[dep_index - 1]
            dep_sw = stage_by_name.get(dep_stage.name)
            if dep_sw is None:
                return False
            return dep_sw.state == "completed"

        def stage_running(stage_workflow: StageWorkflow) -> bool:
            return stage_workflow.process is not None \
                and stage_workflow.process.is_alive()

        def stage_task_active(stage_workflow: StageWorkflow) -> bool:
            # A scheduled stage has an active monitor task (local or remote).
            return stage_workflow.task is not None

        def stage_is_terminal(stage_workflow: StageWorkflow) -> bool:
            if stage_workflow.state == "completed":
                return True
            if stage_workflow.semaphore == "stop" and not stage_running(stage_workflow):
                return True
            return False

        def stage_allowed_to_run(stage_workflow: StageWorkflow) -> bool:
            if self.benchmark_stage is None:
                return True
            return stage_workflow.stage.name == self.benchmark_stage

        for stage_workflow in self.stage_workflows:
            stage_workflow.state = "unstarted"
            if not stage_allowed_to_run(stage_workflow):
                stage_workflow.semaphore = "stop"
                print(f"[benchmark] skipping non-benchmark stage: "
                      f"{stage_workflow.stage.name}")

        print("[reattach] scanning for existing local processes from prior sessions...")
        for stage_workflow in self.stage_workflows:
            if stage_workflow.resource_name != "local":
                continue
            existing_state = workload_local_mp.check_for_existing_local_processes(
                self.model.directory, stage_workflow.stage.name)
            if existing_state is not None:
                print(
                    f"[reattach] found stage {stage_workflow.stage.name} "
                    f"running as PID {existing_state.pid}; monitoring will resume "
                    f"when this stage is scheduled."
                )

        # Start background JSON file polling and keystroke watcher
        stop_event = asyncio.Event()
        self._stop_event = stop_event
        status_task = asyncio.create_task(
            self._status_writer_loop(stop_event))
        keystroke_task = asyncio.create_task(
            self._keystroke_loop(stop_event))
        live_display_task = asyncio.create_task(
            self._live_display_loop(stop_event))
        #monitor_task = asyncio.create_task(poll_job_status(self.telemetry, stop_event))

        try:
            while True:
                if self._detached_requested or stop_event.is_set():
                    break
                all_terminal = True
                for stage_workflow in self.stage_workflows:
                    if self._detached_requested or stop_event.is_set():
                        all_terminal = True
                        break
                    if not stage_allowed_to_run(stage_workflow):
                        continue

                    if stage_workflow.semaphore == "stop":
                        stage_workflow.kill()
                        continue

                    if stage_workflow.semaphore != "go":
                        all_terminal = False
                        continue

                    if stage_workflow.fusion_host is not None:
                        deps_ok = fusion_dependencies_satisfied(
                            stage_workflow, stage_by_name)
                    else:
                        deps_ok = stage_dependencies_completed(stage_workflow)
                    if not deps_ok:
                        all_terminal = False
                        continue

                    if stage_workflow.state == "completed":
                        continue

                    if stage_running(stage_workflow):
                        all_terminal = False
                        continue

                    if stage_task_active(stage_workflow):
                        all_terminal = False
                        continue

                    if stage_workflow.state in {"failed", "unstarted", "unknown"}:
                        # (Re)launch stage when permitted by semaphore and deps.
                        prior_state = stage_workflow.state
                        stage_workflow.dependency_indices = []
                        stage_workflow.dependency_tasks = []
                        stage_workflow.task = None
                        stage_workflow.process = None
                        stage_workflow.state = "unstarted"
                        stage_workflow.progress = 0.0
                        stage_workflow.manager_status = "idle"
                        # Wire cross-resource transfer intent from the upstream
                        # stage (input_stage_index is 1-indexed; 0 means none).
                        upstream_idx = getattr(
                            stage_workflow.stage, "input_stage_index", 0)
                        stage_workflow.transfer_from = None
                        if upstream_idx and upstream_idx > 0:
                            upstream_sw = next(
                                (sw for sw in self.stage_workflows
                                 if sw.stage.index == upstream_idx),
                                None,
                            )
                            if (upstream_sw is not None
                                    and upstream_sw.resource_name
                                    != stage_workflow.resource_name):
                                stage_workflow.transfer_from = \
                                    upstream_sw.resource_name
                        if (prior_state in {"unstarted", "unknown"}
                                and stage_workflow.resource_name != "local"
                                and stage_workflow.fusion_host is None):
                            probe_action, probe_status = (
                                await stage_workflow.probe_remote_launch())
                            if probe_action == "completed":
                                if is_fusion_host(stage_workflow):
                                    mark_fused_set_completed(
                                        stage_workflow, stage_by_name)
                                else:
                                    stage_workflow.state = "completed"
                                    stage_workflow.progress = 1.0
                                continue
                            if probe_action == "reattach":
                                manager_status = (
                                    (probe_status or {})
                                    .get("manager_status") or {})
                                stage_workflow.state = "started"
                                stage_workflow.manager_status = (
                                    determine_stage_manager_status(
                                        manager_status))
                                for job in manager_status.get("jobs", []):
                                    job_id = job.get("JobID")
                                    if job_id:
                                        stage_workflow.job_ids.add(str(job_id))
                                await stage_workflow.create_monitor_only_task()
                                if stage_workflow.task is not None:
                                    launched_tasks.append(stage_workflow.task)
                                    self.stage_tasks.append(
                                        stage_workflow.task)
                                all_terminal = False
                                continue
                        await stage_workflow.create_tasks()
                        if stage_workflow.task is not None:
                            launched_tasks.append(stage_workflow.task)
                            self.stage_tasks.append(stage_workflow.task)
                        all_terminal = False
                        continue

                    if not stage_is_terminal(stage_workflow):
                        all_terminal = False

                if all_terminal:
                    break
                await asyncio.sleep(POLLING_INTERVAL)
        finally:
            # Stop the background monitor, keystroke watcher, and live display
            stop_event.set()
            await status_task
            await keystroke_task
            await live_display_task

        if (not self._detached_requested) and len(launched_tasks) > 0:
            await asyncio.gather(*launched_tasks, return_exceptions=True)
        #monitor_task.cancel()
        #try:
        #    await monitor_task
        #except asyncio.CancelledError:
        #    pass

        # Transfer files back

        if self._detached_requested:
            await self.detach_shutdown()
        elif self._shutdown_task is not None:
            await self._shutdown_task
        else:
            await self.shutdown()
        # NOTE: could summarize telemetry here.

async def launch_seekr_pipeline(
        model: seekr3_structures.Seekr_model,
        seekrflow: structures.Seekrflow,
        force_rerun_stages: set[str] | None = None,
    semaphore_overrides: dict[str, str] | None = None,
        benchmark_stage: str | None = None,
        keystrokes_enabled: bool = True,
        ) -> None:
    backend = await LocalExecutionBackend(ThreadPoolExecutor())
    workflow_engine = await WorkflowEngine.create(backend=backend)
    pipeline = SeekrPipeline(
        model, seekrflow, workflow_engine,
        force_rerun_stages=force_rerun_stages or set(),
        semaphore_overrides=semaphore_overrides or {},
        benchmark_stage=benchmark_stage,
        keystrokes_enabled=keystrokes_enabled)
    await pipeline.run_workflows()
    return

def run_model(
        seekrflow: structures.Seekrflow,
        transfer_before: str | None = None,
        transfer_from_remote_only: str | None = None,
        force_rerun: list[str] | None = None,
        benchmark_stage: str | None = None,
        placement_resource_overrides: dict[str, str] | None = None,
        semaphore_dict: dict[str, str] | None = None,
        keystrokes_enabled: bool = True,
        ) -> None:
    """
    Run the SEEKR calculation using remote and local resources.
    """
    curdir = os.getcwd()
    if seekrflow.work_directory is not None:
        seekrflow.work_directory = os.path.abspath(seekrflow.work_directory)
        os.chdir(seekrflow.work_directory)
    if seekrflow.root_directory is None:
        root_directory = os.path.join(seekrflow.work_directory, structures.ROOT)
    else:
        root_directory = seekrflow.root_directory
    model_filename = os.path.join(root_directory, "model.json")
    model = seekr3_structures.load_model(model_filename)
    structures.validate_run_settings(seekrflow, model)
    if benchmark_stage is not None:
        all_stage_names_for_bench = [s.name for s in model.stages]
        if benchmark_stage not in all_stage_names_for_bench:
            raise ValueError(
                f"Unknown benchmark_stage {benchmark_stage!r}. "
                f"Available stages: {all_stage_names_for_bench}")
        benchmark_index = all_stage_names_for_bench.index(benchmark_stage)
        # Require all dependency ancestors of the stage to be completed.
        ancestor_indices: list[int] = []
        seen_indices: set[int] = {benchmark_index}
        current_index = benchmark_index
        while True:
            parent_one_based = getattr(
                model.stages[current_index], "input_stage_index", 0)
            if parent_one_based <= 0:
                break
            parent_index = parent_one_based - 1
            if parent_index in seen_indices:
                raise ValueError(
                    f"Cycle detected in stage dependencies while "
                    f"resolving benchmark stage {benchmark_stage!r}.")
            seen_indices.add(parent_index)
            ancestor_indices.append(parent_index)
            current_index = parent_index
        if len(ancestor_indices) > 0:
            status_msg = seekr_status.status(
                model, "progress", print_json=True)
            progress_by_stage = status_msg.get("progress", {})
            unfinished_ancestors: list[str] = []
            for ancestor_index in ancestor_indices:
                ancestor_name = model.stages[ancestor_index].name
                ancestor_entry = progress_by_stage.get(ancestor_index, {})
                if not ancestor_entry.get("finished", False):
                    unfinished_ancestors.append(ancestor_name)
            if len(unfinished_ancestors) > 0:
                raise ValueError(
                    f"Cannot benchmark stage {benchmark_stage!r}: "
                    f"unfinished dependency stage(s): "
                    f"{unfinished_ancestors}. Complete dependencies "
                    f"before running benchmark mode.")
        print(f"Benchmark mode enabled for stage {benchmark_stage!r}: "
              f"strict mode is active; only this stage will be run.")
        # Set remote submission times to small values 
        # TODO: move to beginning of seekr stage?
        #for resource in seekrflow.run_settings.resources:
        #    if resource.type == "slurm_remote":
        #        resource.time_limit = "00:30:00"
        
    # Track consecutive empty job checks to avoid premature resubmission
    MAX_EMPTY_CHECKS_BEFORE_RESUBMIT = 10
    #empty_jobs_check_count = {stage: MAX_EMPTY_CHECKS_BEFORE_RESUBMIT 
    #  + 1 for stage in stage_names}
    
    # Track quick failures (jobs that start running but fail within 2 cycles)
    # 2 check cycles (60 seconds)
    MIN_RUNNING_TIME_BEFORE_IDLE = 2 * MAIN_LOOP_INTERVAL  

    # Override resource names if provided (dotted target paths).
    if placement_resource_overrides is not None:
        for dotted_target, resource_name in placement_resource_overrides.items():
            if resource_name is not None:
                seekrflow.run_settings.set_placement_resource(
                    dotted_target.split("."), resource_name)
    
    # Signal handling is installed by SeekrPipeline.run_workflows().

    # Resolve force_rerun: None=nothing, []=all stages, list=those stages.
    all_stage_names = [s.name for s in model.stages]

    def _resolve_transfer_stage_targets(stage_arg: str) -> list[str]:
        if stage_arg == "all":
            return list(all_stage_names)
        if stage_arg not in all_stage_names:
            raise ValueError(
                f"Unknown stage in transfer request: {stage_arg!r}. "
                f"Available stages: {all_stage_names}")
        return [stage_arg]

    def _transfer_unique_stage_resources(
            stage_names: list[str], backwards: bool) -> list[str]:
        resources_by_name: dict[str, structures.Resource_remote_base] = {}
        for stage_name in stage_names:
            try:
                resource = seekrflow.run_settings.get_stage_resource(
                    stage_name, seekrflow.workflow.procedure)
            except ValueError:
                resource = None
            if resource is None:
                continue
            resources_by_name[resource.name] = resource

        for resource in resources_by_name.values():
            transfer_base.transfer_files_to_from_remote_resource(
                seekrflow.name,
                resource,
                root_directory,
                backwards=backwards,
            )
        return list(resources_by_name.keys())
    if force_rerun is None:
        force_rerun_stages: set[str] = set()
    elif len(force_rerun) == 0:
        force_rerun_stages = set(all_stage_names)
    else:
        unknown = [s for s in force_rerun if s not in all_stage_names]
        if unknown:
            raise ValueError(
                f"Unknown stage(s) in force_rerun: {unknown}. "
                f"Available stages: {all_stage_names}")
        force_rerun_stages = set(force_rerun)
    if force_rerun_stages:
        print(f"Force-rerun requested for stages: "
              f"{sorted(force_rerun_stages)}")

    if semaphore_dict is None:
        semaphore_dict = {}
    unknown_semaphore_stages = [
        stage_name for stage_name in semaphore_dict
        if stage_name not in all_stage_names]
    if unknown_semaphore_stages:
        raise ValueError(
            f"Unknown stage(s) in semaphore settings: "
            f"{unknown_semaphore_stages}. "
            f"Available stages: {all_stage_names}")
    invalid_semaphore_values = [
        (stage_name, value)
        for stage_name, value in semaphore_dict.items()
        if value not in {"go", "wait", "stop"}]
    if invalid_semaphore_values:
        raise ValueError(
            f"Invalid semaphore value(s): {invalid_semaphore_values}. "
            "Allowed values are: go, wait, stop.")

    perform_final_transfer = True

    if transfer_before is not None:
        target_stages = _resolve_transfer_stage_targets(transfer_before)
        transferred_resources = _transfer_unique_stage_resources(
            target_stages, backwards=False)
        if len(transferred_resources) > 0:
            print(
                f"[transfer_before] transferred to resource(s): "
                f"{transferred_resources}")
        else:
            print("[transfer_before] no non-local stage resources to transfer.")

    if transfer_from_remote_only is not None:
        target_stages = _resolve_transfer_stage_targets(transfer_from_remote_only)
        transferred_resources = _transfer_unique_stage_resources(
            target_stages, backwards=True)
        if len(transferred_resources) > 0:
            print(
                f"[transfer_from_remote_only] transferred from resource(s): "
                f"{transferred_resources}")
        else:
            print(
                "[transfer_from_remote_only] no non-local stage resources "
                "to transfer.")
        os.chdir(curdir)
        return
    
    # Detect jobs that start running but fail quickly (within 2 check cycles)
            
    # handle transfer_from_host/remote_only
    
    # Optimize SEEKR time limit based on remaining work and benchmark data
    
    # Only cleanup processes if NOT detaching (perform_final_transfer is our indicator)
    # When detaching (perform_final_transfer=False from 'd' command), leave processes running
    
    # Perform final transfer

    asyncio.run(launch_seekr_pipeline(
        model, seekrflow, force_rerun_stages, semaphore_dict,
        benchmark_stage,
        keystrokes_enabled=keystrokes_enabled))

    if perform_final_transfer:
        detached_requested = False
        completed_stage_resources: dict[str, structures.Resource_remote_base] = {}
        try:
            snapshot = load_job_status(root_directory)
            detached_requested = bool(snapshot.get("detached_requested", False))
            stages_snapshot = snapshot.get("stages", {})
            for _, stage_info in stages_snapshot.items():
                if stage_info.get("state") != "completed":
                    continue
                resource_name = stage_info.get("resource_name", "local")
                if resource_name == "local":
                    continue
                resource = seekrflow.run_settings.get_resource_by_name(resource_name)
                if resource is not None:
                    completed_stage_resources[resource_name] = resource
        except Exception as e:
            print(f"[warning] final transfer discovery failed: {e}")

        if detached_requested:
            print("[transfer] skipping final backward transfer due to detach.")
        else:
            for resource in completed_stage_resources.values():
                try:
                    transfer_base.transfer_files_to_from_remote_resource(
                        seekrflow.name,
                        resource,
                        root_directory,
                        backwards=True,
                    )
                except Exception as e:
                    print(
                        f"[warning] final backward transfer failed for "
                        f"resource {resource.name}: {e}")

    os.chdir(curdir)
    