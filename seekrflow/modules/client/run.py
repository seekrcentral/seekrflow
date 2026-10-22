"""
modules/client/run.py

Handle client run processes.
"""

import os
import math
import time
import json
import fcntl
import signal
import typing
import asyncio
import datetime
from typing import Dict, List, Optional
from concurrent.futures import ThreadPoolExecutor

import attrs
import seekr.modules.structures as seekr_structures
from radical.asyncflow import WorkflowEngine, LocalExecutionBackend
import seekr.modules.structures as seekr_structures
import seekr.modules.scales.base as scales_base

import seekrflow.modules.base as base
import seekrflow.modules.structures as structures
import seekrflow.modules.transfer.base as transfer_base
import seekrflow.modules.client.structures as client_structures
import seekrflow.modules.client.validation as client_validation
import seekrflow.modules.workload_managers.remote_local as workload_remote_local
import seekrflow.modules.job.structures as job_structures

MAIN_LOOP_INTERVAL = 5.0 # seconds
STATUS_WRITE_INTERVAL = 5.0 # seconds
CONTROL_FILE_READ_INTERVAL = 5.0 # seconds
STATUS_SCHEMA_VERSION = 2
SHUTDOWN_CANCEL_TIMEOUT = 30.0
ENGINE_SHUTDOWN_TIMEOUT = 15.0
REMOTE_STATUS_WRITE_INTERVAL = 15.0 # seconds

async def _run_blocking(fn: typing.Callable[..., typing.Any], *args, **kwargs):
    """
    Run a blocking callable off the asyncio event loop (thread pool).
    Used for Globus / remote status and submit calls so client stays 
    responsive while endpoints are slow.
    """
    loop = asyncio.get_running_loop()
    return await loop.run_in_executor(
        None, lambda: fn(*args, **kwargs))

def _format_hms(seconds: float) -> str:
    """
    Convert seconds to a human-readable string in the format of HH:MM:SS.
    """
    whole = max(0, math.ceil(seconds))
    hours, rem = divmod(whole, 3600)
    minutes, secs = divmod(rem, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"

def _elapsed_seconds(values: list) -> float | None:
    """
    Longest scheduler elapsed time, or None when none is known.
    """
    best = None
    for value in values:
        if not isinstance(value, str) or not value:
            continue
        days, _, clock = value.partition("-")
        if not clock:
            days, clock = "0", days
        parts = [int(part) for part in clock.split(":")]
        while len(parts) < 3:
            parts.insert(0, 0)
        seconds = ((int(days) * 24 + parts[0]) * 60 + parts[1]) * 60 + parts[2]
        if best is None or seconds > best:
            best = seconds
    return best

def run_units_from_info(
        dimensions: list[dict],
        number_of_anchors: int,
        number_of_swarms: int,
    ) -> list[job_structures.RunUnit]:
    """
    Generate run units grid from dimensions and numbers of anchors and swarms.
    """
    if len(dimensions) == 0:
        return [job_structures.RunUnit(anchor="any", swarm_id=None)]

    if "swarm" in dimensions:
        if number_of_swarms <= 1:
            raise ValueError(
                f"Dispatch dimension 'swarm' requires num_swarms > 1, "
                f"got {number_of_swarms}.")
    # Assign the grid of run units
    if dimensions == ["swarm"]:
        return [
            job_structures.RunUnit(anchor="any", swarm_id=swarm)
            for swarm in range(number_of_swarms)
        ]

    if "anchor" in dimensions:
        if (number_of_anchors is None) or (number_of_anchors <= 0):
            raise ValueError(
                "Dispatch dimension 'anchor' requires a positive num_anchors."
                f"got {number_of_anchors}.")
    
    if dimensions == ["anchor"]:
        return [
            job_structures.RunUnit(anchor=anchor, swarm_id=None)
            for anchor in range(number_of_anchors)
        ]
    
    if dimensions == ["anchor", "swarm"]:
        units: list[job_structures.RunUnit] = []
        for anchor in range(number_of_anchors):
            for swarm in range(number_of_swarms):
                units.append(job_structures.RunUnit(anchor=anchor, swarm_id=swarm))
        return units

    raise ValueError(f"Unsupported dispatch dimensions: {dimensions!r}.")
    
# TODO: move to a client_stage_workflow.py module?
@attrs.define
class StageWorkflow:
    """
    A workflow for a single stage of SEEKR3.
    """
    # Constructor args (no defaults → required, in this order)
    model: seekr_structures.Seekr_model = attrs.field(repr=False)
    seekrflow: structures.Seekrflow = attrs.field(repr=False)
    stage_list: List[scales_base.Base_stage] = attrs.field(repr=False)
    workflow_engine: WorkflowEngine = attrs.field(repr=False)
    resource: structures.Resource_base | None = attrs.field(repr=False)
    resolved_execution: structures.Resolved_execution | None = attrs.field(
        default=None, repr=False)
    force_overwrite: Dict[str, bool] = attrs.field(factory=dict)
    benchmark_mode: Dict[str, bool] = attrs.field(factory=dict)
    
    # Derived / mutable state
    dependency_tasks: list = attrs.field(factory=list)
    task: typing.Any = attrs.field(default=None)
    stage_state: Dict[str, str] = attrs.field(factory=dict)
    progress: Dict[str, dict[str, float]] = attrs.field(factory=dict)
    start_step: Dict[str, dict[str, int]] = attrs.field(factory=dict)
    current_step: Dict[str, dict[str, int]] = attrs.field(factory=dict)
    total_steps: Dict[str, dict[str, int]] = attrs.field(factory=dict)
    time_of_first_progress: Dict[str, dict[str, float]] = attrs.field(factory=dict)
    time_of_last_progress: Dict[str, dict[str, float]] = attrs.field(factory=dict)
    # Accumulated steps and seconds over multiple jobs for a stage
    accounted_steps: Dict[str, dict[str, float]] = attrs.field(factory=dict)
    accounted_seconds: Dict[str, dict[str, float]] = attrs.field(factory=dict)
    accounted_internal_id: Optional[int] = attrs.field(default=None)

    manager_status: str = attrs.field(
        default="idle", validator=attrs.validators.in_(
            {"pending", "queued", "running", "queued/running", "idle", "failed", 
            "cancelled", "pushing", "pulling"}))
    semaphore: str = attrs.field(
        default="go", validator=attrs.validators.in_(
            {"go", "wait", "stop"}))
    internal_id: Optional[int] = attrs.field(default=None)
    job_id: Optional[str] = attrs.field(default=None)
    job_name: Optional[str] = attrs.field(default=None)
    array_indices: Optional[List[int]] = attrs.field(default=None)

    transfer_requested: bool = attrs.field(default=False)
    transfer_status: str = attrs.field(default="idle")
    transfer_direction: str | None = attrs.field(default=None)
    transfer_error: str | None = attrs.field(default=None)
    transfer_relaunch_count: int = attrs.field(default=0)
    transfer_from: str | None = attrs.field(default=None)
    last_error: str | None = attrs.field(default=None)
    telemetry_poll_interval: float | None = attrs.field(default=None)
    detached_requested: bool = attrs.field(default=False)
    control_file_semaphore: str | None = attrs.field(default=None, repr=False)

    submitted_at: float | None = attrs.field(default=None)
    progress_at_submit: dict[str, dict] = attrs.field(factory=dict)
    monitor_generation: int = attrs.field(default=0)
    failed_job_id: str | None = attrs.field(default=None)
    
    async def probe_launch(self) -> str:
        """
        Probe for any running jobs for our stages and return the action and status.
        """
        if self.internal_id is None:
            return "submit"
        root_dir = workload_remote_local.resolve_model_directory(
            self.seekrflow, self.resource)
        manager_payload = {
            "system_payloads": {
                self.seekrflow.name: {
                    "root_dir": root_dir,
                    "jobs": [
                        {
                            "internal_id": self.internal_id,
                            "job_id": self.job_id,
                            "job_name": self.job_name,
                            "stage_indices": [stage.index for stage in self.stage_list],
                            "array_indices": self.array_indices,
                        },
                    ],
                },
            }
        }

        try:
            result = await _run_blocking(
                workload_remote_local.status,
                self.resource,
                manager_payload,
                silent=True,
            )
        except Exception as e:
            if isinstance(e, ImportError):
                self.last_error = f"probe failed: {e}"
                self.control_file_semaphore = "wait"
                self.semaphore = "wait"
                self.manager_status = "failed"
                for stage in self.stage_list:
                    self.stage_state[stage.name] = "error"
                print(
                    f"[remote-probe] stage {self.stage_list[0].name}: "
                    f"{self.last_error}; setting semaphore=wait"
                )
                return "wait"
            print(
                f"[remote-probe] stage {self.stage_list[0].name}: "
                f"probe failed ({e}); not submitting until state can be read"
            )
            if self.job_id is not None:
                return "reattach"
            return "defer"

        job_dict = (
            (result.get("payload") or {})
            .get(self.seekrflow.name, {})
            .get(self.job_id)
        )
        if not job_dict:
            return "submit"
        self.apply_stage_dicts(job_dict.get("stage_dicts_by_stage_index") or {})
        states = {
            entry.get("state")
            for entry in (job_dict.get("manager_dicts_by_array_index", {}).values())
        }
        if states & {"running", "queued", "pending"}:
            return "reattach"
        stage_dicts = job_dict.get("stage_dicts_by_stage_index") or {}
        if stage_dicts and all(
                entry.get("finished") or entry.get("state") == "completed" 
                for entry in stage_dicts.values()):
            return "completed"
        return "submit"
        
    def _outbound_transfer_resources(self):
        """
        Source/destination resources for the outbound copy, or (None, None).
        """
        if self.transfer_from is None:
            dep_index = getattr(self.stage_list[0], "input_stage_index", 0)
            if self.resource.name != "local" and dep_index <= 0:
                return None, self.resource
            return None, None
        src_resource = self.systemrun.seekrflow.run_settings.get_resource_by_name(
            self.transfer_from)
        return src_resource, self.resource

    def _fold_latest_job(self) -> None:
        """Add this internal_id's steps and seconds once."""
        if self.internal_id is None or self.internal_id == self.accounted_internal_id:
            return
        for stage in self.stage_list:
            name = stage.name
            starts = self.start_step.get(name) or {}
            currents = self.current_step.get(name) or {}
            firsts = self.time_of_first_progress.get(name) or {}
            lasts = self.time_of_last_progress.get(name) or {}
            step_dest = self.accounted_steps.setdefault(name, {})
            time_dest = self.accounted_seconds.setdefault(name, {})
            for key, current in currents.items():
                start = starts.get(key)
                t0 = firsts.get(key)
                t1 = lasts.get(key)
                if start is None or t0 is None or t1 is None:
                    continue
                steps = current - start
                elapsed = t1 - t0
                if steps > 0 and elapsed > 0:
                    step_dest[key] = step_dest.get(key, 0) + steps
                    time_dest[key] = time_dest.get(key, 0) + elapsed
        self.accounted_internal_id = self.internal_id
    
    def _hint_rate_per_second(
            self, 
            stage, 
            hint: float
            ) -> float | None:
        """
        Account for the estimated_performance is in ns/day for MD and 
        trajectories/day for BD.
        """
        kind = getattr(stage, "scale_type", None)
        per_day = hint
        if kind == "molecular_dynamics":
            try:
                timestep_ps = self.model.get_timestep_by_type(kind)
            except (ValueError, AssertionError):
                return None
            if not timestep_ps:
                return None
            per_day *= 1000.0 / timestep_ps  # ns/day → steps/day
        elif kind != "brownian_dynamics":
            return None
        return per_day / 86400.0

    def _compute_submit_time_limit_override(
            self,
            units_grid: list,
            group_size: int | None,
            concurrency: int | None,
            ) -> str | None:
        if self.resolved_execution is None:
            return None
        if any(self.benchmark_mode.values()):
            # TODO: the only problem here is that there might be a benchmark on a
            # later stage in the list, which will need a longer time limit.
            return base.BENCHMARK_REMOTE_TIME_LIMIT
        policy = self.resolved_execution.time_policy
        if not isinstance(policy, structures.Time_policy_adaptive):
            return None
        self._fold_latest_job()
        unit_seconds: list[float] = []
        for unit in units_grid:
            key = f"anchor_{unit.anchor}_swarm_{unit.swarm_id}"
            seconds = 0.0
            for stage in self.stage_list:
                name = stage.name
                total = (self.total_steps.get(name) or {}).get(key)
                current = (self.current_step.get(name) or {}).get(key, 0)
                done = (self.accounted_steps.get(name) or {}).get(key)
                elapsed = (self.accounted_seconds.get(name) or {}).get(key)
                if total is None:
                    return None
                remaining = max(0, total - current)
                if remaining == 0:
                    continue
                if done and elapsed:
                    rate = done / elapsed
                else:
                    if policy.estimated_performance is None:
                        return None
                    rate = self._hint_rate_per_second(stage, policy.estimated_performance)
                    if rate is None:
                        return None
                seconds += remaining / rate
            unit_seconds.append(seconds)
        if not unit_seconds:
            return None
        width = concurrency or 1
        size = group_size or len(unit_seconds)
        member_times = []
        for start in range(0, len(unit_seconds), size):
            chunk = unit_seconds[start:start + size]
            member = 0.0
            for wave in range(0, len(chunk), width):
                member += max(chunk[wave:wave + width])
            member_times.append(member)
        seconds = max(member_times) * (1.0 + policy.safety_factor)
        cap_text = policy.max_time_limit or getattr(self.resource, "max_time_limit", None)
        if cap_text is not None:
            seconds = min(seconds, client_validation.time_limit_to_seconds(cap_text))
        if policy.min_time_limit is not None:
            seconds = max(seconds, client_validation.time_limit_to_seconds(policy.min_time_limit))
        return _format_hms(seconds)

    def move_files(
            self,
            backwards: bool = False
            ) -> None:
        """
        Copy this system's files to/from the remote workdir.
        """
        if backwards:
            steps = [(self.resource, True)] if self.resource is not None else []
        else:
            src, dst = self._outbound_transfer_resources()
            steps = []
            if src is not None:
                steps.append((src, True))
            if dst is not None:
                steps.append((dst, False))
        if not steps:
            self.transfer_status = "skipped"
            self.transfer_error = None
            return
        self.transfer_status = "running"
        self.transfer_direction = "in" if backwards else "out"
        prior_status = self.manager_status
        self.manager_status = "pulling" if backwards else "pushing"
        try:
            for resource, pull in steps:
                transfer_base.transfer_files_to_from_remote_resource(
                    self.seekrflow.name, resource, self.model.directory,
                    backwards=pull, seekrflow_name=self.seekrflow.name)
        except Exception as e:
            self.transfer_error = str(e)
            self.control_file_semaphore = "wait"
            self.semaphore = "wait"
            self.manager_status = "failed"
            for stage in self.stage_list:
                self.stage_state[stage.name] = "error"
            self.last_error = str(e)
            return

        self.transfer_status = "completed"
        self.transfer_direction = None
        self.manager_status = prior_status # if backwards else "idle"
        self.transfer_relaunch_count = 0

    def apply_stage_dicts(
            self, 
            stage_dicts: dict
            ) -> None:
        """
        Copy one status payload onto this workflow.
        """
        series = (
            ("stage_anchor_swarm_progress_list", self.progress),
            ("stage_anchor_swarm_starting_step_list", self.start_step),
            ("stage_anchor_swarm_current_step_list", self.current_step),
            ("stage_anchor_swarm_total_steps_list", self.total_steps),
            ("stage_anchor_swarm_time_of_first_progress_list", self.time_of_first_progress),
            ("stage_anchor_swarm_time_of_last_progress_list", self.time_of_last_progress),
        )
        for stage in self.stage_list:
            stage_dict = stage_dicts.get(stage.index)
            if stage_dict is None:
                stage_dict = stage_dicts.get(str(stage.index))
            if not stage_dict:
                continue
            self.stage_state[stage.name] = stage_dict.get("state", "unknown")
            for field, destination in series:
                per_swarm = {}
                for row in stage_dict.get(field) or []:
                    key_stage_index, key_anchor, key_swarm_id, value = row
                    if key_stage_index == stage.index:
                        per_swarm[f"anchor_{key_anchor}_swarm_{key_swarm_id}"] = value
                destination[stage.name] = per_swarm
        
    async def create_tasks(self) -> None:
        """
        Register tasks for running this stage_workflow.
        """
        self.manager_status = "queued"
        # Transfer files if dependent stage resource is different
        @self.workflow_engine.function_task
        async def transfer_files(*args):
            if self.transfer_from is None:
                # For the first remote stage (no dependency), seed remote workdir
                # from local so model/config files exist before submission.
                dep_index = getattr(self.stage_list[0], "input_stage_index", 0)
                if self.resource.name != "local" and dep_index <= 0:
                    pass
                else:
                    self.transfer_status = "skipped"
                    self.transfer_error = None
                    return
            await _run_blocking(self.move_files, backwards=False)

        dep_index = getattr(self.stage_list[0], "input_stage_index", 0)
        should_transfer = bool(self.transfer_from) or (
            self.resource.name != "local" and dep_index <= 0
        )
        if should_transfer:
            self.dependency_tasks.append(transfer_files(*self.dependency_tasks))

        # Run stage
        @self.workflow_engine.function_task
        async def run_stage(*args):
            if self.semaphore in ("wait", "stop"):
                return
            stage_names = [stage.name for stage in self.stage_list]
            if self.resource is None:
                raise Exception(
                    f"Remote resource config missing for {self.resource.name!r}")
            try:
                remote_or_cloud = workload_remote_local.resource_kind(
                    self.resource)
                #if remote_or_cloud == "remote":
                destination_path = (
                    workload_remote_local.resolve_model_directory(
                        self.seekrflow, self.resource))
                destination_model_filename = os.path.join(
                    destination_path, "model.json")

                # If any stage is force-overwrite, cancel the whole job
                for stage_name, force_now in self.force_overwrite.items():
                    if force_now and self.job_id is not None:
                        manager_payload = {
                            "remove_json_files": False,
                            "system_payloads": {
                                self.seekrflow.name: {
                                    "root_dir": destination_path,
                                    "jobs": [
                                        {
                                            "job_id": self.job_id,
                                            "job_name": self.job_name,
                                        },
                                    ],
                                },
                            },
                        }
                        workload_remote_local.submit_cancel_workload(
                            self.resource,
                            manager_payload,
                            silent=True,
                        )
                        self.job_id = None
                        break

                # prepare stage specs, run units, and job spec
                stage_specs = []
                for stage in self.stage_list:
                    stage_spec = job_structures.StageSpec(
                        stage_index=stage.index,
                        stage_name=stage.name,
                        force_overwrite=self.force_overwrite[stage.name],
                        benchmark=self.benchmark_mode[stage.name],
                    )
                    stage_specs.append(stage_spec)
                
                force_any = any(self.force_overwrite.values())
                for stage_name in self.force_overwrite:
                    self.force_overwrite[stage_name] = False

                dispatch = (
                    self.resolved_execution.dispatch \
                        if self.resolved_execution is not None else None)
                dimensions = dispatch.dimensions if dispatch is not None else None
                group_size = dispatch.group_size if dispatch is not None else None
                concurrency = dispatch.concurrency if dispatch is not None else None
                
                if bool(dimensions):
                    # Unit enumeration is needed.
                    unit_counts = workload_remote_local.fetch_unit_counts(
                        self.seekrflow,
                        self.stage_list[0],
                        self.resource,
                        silent=True,
                    )
                    number_of_anchors = unit_counts.num_anchors
                    number_of_swarms = unit_counts.num_swarms
                else:
                    number_of_anchors = 0
                    number_of_swarms = 1
                
                units_grid = run_units_from_info(
                    dimensions or [], number_of_anchors, number_of_swarms)
                size = group_size or len(units_grid)
                array_size = 0
                job_specs = []
                self.array_indices = []
                for array_index, start in enumerate(range(0, len(units_grid), size)):
                    units_chunk = units_grid[start:start + size]
                    # one JobSpec: same stage_specs, run_units=chunk,
                    # concurrency=min(concurrency, len(chunk)), array_index=array_index
                    job_spec = job_structures.JobSpec(
                        internal_id=0, # Assign on the remote
                        remote_root_dir=destination_path,
                        status_write_interval=REMOTE_STATUS_WRITE_INTERVAL,
                        telemetry_poll_interval=self.telemetry_poll_interval,
                        stage_specs=stage_specs,
                        run_units=units_chunk,
                        concurrency=min(concurrency, len(units_chunk)),
                        array_index=array_index,
                    )
                    self.array_indices.append(array_index)
                    job_specs.append(job_spec)
                    array_size += 1

                # Retrieve existing jobs, if any
                if not force_any:
                    # If our stages have been run (or running), and no force overwrite,
                    pre_action = await self.probe_launch()
                    if pre_action == "completed":
                        for stage in self.stage_list:
                            self.stage_state[stage.name] = "completed"
                        return
                    if pre_action == "reattach":
                        self.manager_status = "running"
                        return
                    if pre_action in ("defer", "wait"):
                        return

                time_limit_override = self._compute_submit_time_limit_override(
                    units_grid, group_size, concurrency)

                # Finally, submit the job
                run_result = await _run_blocking(
                    workload_remote_local.submit_job,
                    self.seekrflow,
                    self.resource,
                    stage_names,
                    job_specs,
                    self.resolved_execution,
                    time_limit_override,
                )
                if not isinstance(run_result, dict):
                    raise RuntimeError(
                        f"Remote submit returned invalid payload: {run_result!r}"
                    )
                if run_result.get("success") is False:
                    raise RuntimeError(
                        f"Remote submit failed for stage {self.stage_list[0].name}: "
                        f"{run_result.get('error')}"
                    )
                if "success" not in run_result:
                    raise RuntimeError(
                        f"Remote submit returned payload without success flag: {run_result!r}"
                    )
                internal_id = run_result["internal_id"]
                job_id = run_result["job_id"]
                job_name = run_result["job_name"]
                for job_spec in job_specs:
                    job_spec.internal_id = internal_id
                self.internal_id = internal_id
                self.job_id = job_id
                self.submitted_at = time.time()
                self.progress_at_submit = {
                    name: dict(values) for name, values in self.progress.items()
                }
                self.job_name = job_name
                self.manager_status = "queued"
                    
            except Exception as e:
                for stage_name in stage_names:
                    self.stage_state[stage_name] = "error"
                self.last_error = f"remote submit failed: {e}"
                self.control_file_semaphore = "wait"
                self.semaphore = "wait"
                self.manager_status = "failed"
                print(
                    f"[remote-submit] stage {stage_names[0]} submit failed; "
                    f"setting semaphore=wait. error={e}"
                )
                raise
            return

        self.task = run_stage(*self.dependency_tasks)
        await asyncio.sleep(0)

        @self.workflow_engine.function_task
        async def monitor_stage(task):
            await self._monitor_stage_loop()

        self.task = monitor_stage(self.task)
        await asyncio.sleep(0) # let monitor_stage's async_wrapper register
        return
    
    async def _monitor_stage_loop(self) -> None:
        """
        Monitor for transfer requests.
        """
        generation = self.monitor_generation
        while True:
            if (
                self.detached_requested
                or self.semaphore in ("stop", "wait")
                or self.monitor_generation != generation
                or self.is_terminal()
            ):
                break
            if self.transfer_requested:
                self.transfer_requested = False
                await _run_blocking(self.move_files, backwards=True)
            await asyncio.sleep(STATUS_WRITE_INTERVAL)

    def status_snapshot(self) -> dict:
        """
        Return the JSON-serializable status snapshot of this stage's current
        state.
        """
        process_info = None
        stages = {}
        for stage in self.stage_list:
            stages[stage.name] = {
                "name": stage.name,
                "index": stage.index,
                "scale": getattr(stage, "scale_type", None),
                "stage_state": self.stage_state.get(stage.name, "unknown"),
                "force_overwrite": self.force_overwrite[stage.name],
                "benchmark_mode": self.benchmark_mode[stage.name],
                "semaphore": self.semaphore,
                "progress": self.progress.get(stage.name, {}),
            }
        return {
            "resource_name": self.resource.name if self.resource is not None else "local",
            "manager_status": self.manager_status,
            "transfer_status": self.transfer_status,
            "transfer_direction": self.transfer_direction,
            "transfer_error": self.transfer_error,
            "transfer_relaunch_count": self.transfer_relaunch_count,
            "last_error": self.last_error,
            "internal_id": self.internal_id,
            "job_id": self.job_id,
            "job_name": self.job_name,
            "array_indices": self.array_indices,
            "stages": stages,
        }

    def kill(self) -> None:
        """
        Stop this stage's monitor loop and terminate any running job.
        Idempotent; safe to call on stages that never started.
        """
        self.semaphore = "stop"
        #self.manager_status = "idle"
        if self.resource is None or self.job_id is None:
            return
        root_dir = workload_remote_local.resolve_model_directory(
            self.seekrflow, self.resource)
        manager_payload = {
            "remove_json_files": False,
            "system_payloads": {
                self.seekrflow.name: {
                    "root_dir": root_dir,
                    "jobs": [
                        {
                            "job_id": self.job_id,
                            "job_name": self.job_name,
                        },
                    ],
                },
            },
        }
        try:
            workload_remote_local.submit_cancel_workload(
                self.resource, manager_payload, silent=True)
        except Exception as error:
            print(f"[cancel] stage {self.stage_list[0].name}: {error}")
            return
        self.job_id = None

    def has_running_task_or_process(self) -> bool:
        """
        Check if this stage has any running tasks or processes.
        """
        return self.manager_status in {"running", "queued", "queued/running"}

    def is_terminal(self) -> bool:
        """
        Check if this stage is terminal.
        """
        if self.stage_list and all(
                self.stage_state.get(stage.name) == "completed"
                for stage in self.stage_list):
            return True
        if self.semaphore == "stop" and not self.has_running_task_or_process():
            return True
        return False

def _link_to_previous(prev_resolved, cur_resolved) -> bool:
    """
    Check if the current stage is linked to the previous stage.
    """
    if cur_resolved.co_schedule_with == "predecessor":
        return True
    if prev_resolved.co_schedule_with == "successor":
        return True
    return False

@attrs.define
class SeekrPipeline:
    """
    The entire set of workflows for a seekr run.
    """
    # Constructor args (no defaults → required, in this order)
    systemrun: client_structures.SystemRun = attrs.field(repr=False)
    workflow_engine: WorkflowEngine = attrs.field(repr=False)
    
    # Derived / mutable state
    stage_workflows: list[StageWorkflow] \
        = attrs.field(factory=list, repr=False)
    stage_tasks: list = attrs.field(factory=list, repr=False)
    detached_requested: bool = attrs.field(default=False, repr=False)
    shutting_down: bool = attrs.field(default=False, repr=False)
    _shutdown_task: asyncio.Task | None = attrs.field(
        default=None, repr=False)
    _stop_event: asyncio.Event | None = attrs.field(
        default=None, repr=False)
    limit_semaphores_by_resource: dict = attrs.field(factory=dict, repr=False)
    
    def __attrs_post_init__(
            self, 
            ) -> None:
        # Assign stage lists by resolved.co_schedule_with fields
        new_stage_lists:List[List[scales_base.Base_stage]] = []
        new_resolved_lists:List[List[structures.Resolved_execution]] = []
        current: List[scales_base.Base_stage] = []
        current_resolved: List[structures.Resolved_execution] = []
        for stage in self.systemrun.model.stages:
            resolved = client_validation.resolve_stage_execution(
                self.systemrun.seekrflow.run_settings, stage.name, 
                self.systemrun.seekrflow.workflow.procedure)
            if current and _link_to_previous(prev_resolved, resolved):
                current.append(stage)
                current_resolved.append(resolved)
            else:
                if current:
                    new_stage_lists.append(current)
                    new_resolved_lists.append(current_resolved)
                current = [stage]
                current_resolved = [resolved]
            prev_resolved = resolved
        if current:
            new_stage_lists.append(current)
            new_resolved_lists.append(current_resolved)

        for stage_list, resolved_list \
                in zip(
                    new_stage_lists, new_resolved_lists):
            self.add_stage_list(stage_list, resolved_list)
            
    
    def add_stage_list(
            self, 
            stage_list: List[scales_base.Base_stage],
            resolved_list: List[structures.Resolved_execution]
            ) -> None:
        """
        Add a seekr stage to the pipeline.
        """
        group_resource = None
        for resolved in resolved_list:
            resource = resolved.resource
            if group_resource is None:
                group_resource = resource
            elif group_resource.name != resource.name:
                raise ValueError(
                    f"All stages in a group must have the same resource: "\
                    f"{group_resource.name} != {resource.name}")
        stage_workflow = StageWorkflow(
            self.systemrun.model, self.systemrun.seekrflow, stage_list, self.workflow_engine, 
            resource)
        stage_workflow.telemetry_poll_interval = self.systemrun.telemetry_poll_interval
        stage_workflow.resolved_execution = resolved
        stage_workflow_semaphore = None
        for stage in stage_list:
            stage_workflow.stage_state[stage.name] = "unstarted"
            stage_workflow.force_overwrite[stage.name] \
                = stage.name in self.systemrun.force_rerun_stages
            stage_workflow.benchmark_mode[stage.name] \
                = (stage.name == self.systemrun.benchmark_stage)
            stage_semaphore = self.systemrun.semaphore_dict.get(stage.name, "go")
            if stage_workflow_semaphore is None:
                stage_workflow_semaphore = stage_semaphore
            elif stage_workflow_semaphore != stage_semaphore:
                    raise ValueError(
                        f"All stages in a group must have the same semaphore: "\
                        f"{stage_workflow_semaphore} != {stage_semaphore}")
        if stage_workflow_semaphore is None:
            stage_workflow_semaphore = "go"
        stage_workflow.semaphore = stage_workflow_semaphore
        self.stage_workflows.append(stage_workflow)
    
    def _cancel_all_stages(self) -> None:
        """
        Cancel all stages in the pipeline.
        """
        for stage_workflow in self.stage_workflows:
            stage_workflow.kill()
        return

    async def shutdown(self) -> None:
        """
        Kill all stages and shut down the workflow engine.
        """
        if self.shutting_down:
            return
        self.shutting_down = True
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
        if self._stop_event is not None:
            self._stop_event.set()
        # Capture terminal state before tearing down the engine.
        #self.write_status_snapshot()
        #try:
        #    await asyncio.wait_for(
        #        self.workflow_engine.shutdown(),
        #        timeout=ENGINE_SHUTDOWN_TIMEOUT,
        #    )
        #except asyncio.TimeoutError:
        #    timed_out = True
        #    print(
        #        f"[shutdown] workflow engine shutdown timed out after "
        #        f"{ENGINE_SHUTDOWN_TIMEOUT}s; forcing exit")
        if timed_out:
            os._exit(1)
        return

    async def detach_shutdown(self) -> None:
        """
        Soft shutdown for detach: stop local monitoring/UI loops and close
        the workflow engine, but do not kill running stage processes.
        """
        if self.shutting_down:
            return
        self.shutting_down = True
        if self._stop_event is not None:
            self._stop_event.set()
        #self.write_status_snapshot()
        #await self.workflow_engine.shutdown()
        return

    def status_snapshot(self) -> dict:
        """
        Get the current status snapshot of the pipeline.
        """
        root = self.systemrun.model.directory
        now = time.time()
        stage_workflows = {
            sw_index: sw.status_snapshot() \
                for sw_index, sw in enumerate(self.stage_workflows)
        }
        return {
            "schema_version": STATUS_SCHEMA_VERSION,
            "timestamp": now,
            "timestamp_iso": datetime.datetime.fromtimestamp(
                now, tz=datetime.timezone.utc).isoformat(),
            "pid": os.getpid(),
            "work_directory": self.systemrun.seekrflow.work_directory,
            "root_directory": root,
            "shutting_down": self.shutting_down,
            "detached_requested": self.detached_requested,
            "stage_workflows": stage_workflows,
        }

    async def run_workflows(self) -> None:
        """
        Run the workflows for the entire pipeline.
        This involves:
        - Check for detach or stop conditions
        - For each stage workflow
            - Check for remaining dependency tasks
            - Launch the create_tasks

        """
        #stage_by_name = {sw.stages.name: sw for sw in self.stage_workflows}
        stage_workflow_by_name = {}
        for stage_workflow in self.stage_workflows:
            for stage in stage_workflow.stage_list:
                stage_workflow_by_name[stage.name] = stage_workflow

        def stages_dependencies_satisfied(stage_workflow: StageWorkflow) -> bool:
            """
            Check if all stages in this stage_workflow have their dependencies satisfied.
            """
            own_names = {stage.name for stage in stage_workflow.stage_list}
            for stage in stage_workflow.stage_list:
                dep_index = stage.input_stage_index
                if dep_index <= 0:
                    continue
                dep_stage = self.systemrun.model.stages[dep_index - 1]
                if dep_stage.name in own_names:
                    continue
                sw = stage_workflow_by_name[dep_stage.name]
                if sw.stage_state[dep_stage.name] != "completed":
                    return False
            return True

        stop_event = asyncio.Event()
        self._stop_event = stop_event
        held_slots: list[StageWorkflow] = []
        def release_slot(stage_workflow: StageWorkflow) -> None:
            """
            Release a slot for the stage workflow.
            """
            if stage_workflow not in held_slots:
                return
            held_slots.remove(stage_workflow)
            sem = self.limit_semaphores_by_resource.get(stage_workflow.resource.name)
            if sem is not None:
                sem.release()

        try:
            while True:
                if self.detached_requested or stop_event.is_set():
                    break
                all_terminal = True
                for stage_workflow in self.stage_workflows:
                    # If we need to detach or kill the client
                    if self.detached_requested or stop_event.is_set():
                        all_terminal = True
                        break
                    
                    if stage_workflow in held_slots and (
                            stage_workflow.is_terminal()
                            or (stage_workflow.semaphore == "wait"
                                and not stage_workflow.has_running_task_or_process())):
                        release_slot(stage_workflow)

                    # If all stage semaphores in this stage_workflow are stopped
                    if stage_workflow.semaphore == "stop":
                        await _run_blocking(stage_workflow.kill)
                        continue

                    # If stage dependencies are not satisfied
                    if not stages_dependencies_satisfied(stage_workflow):
                        all_terminal = False
                        continue

                    # If stage_workflow.stage_state is "completed" for all stages
                    if all(stage_workflow.stage_state[stage.name] == "completed" \
                            for stage in stage_workflow.stage_list):
                        continue

                    # If a stage in the workflow process or task is already running
                    if stage_workflow.has_running_task_or_process():
                        all_terminal = False
                        continue

                    if stage_workflow.task is not None and not stage_workflow.task.done():
                        all_terminal = False
                        continue

                    prior_state = stage_workflow.manager_status
                    if prior_state in {
                           "idle", "failed", "unstarted", "unknown", "queued", "started"} \
                            and stage_workflow.semaphore == "go":
                        # Nothing is currently running:
                        # (Re)launch stage when permitted by semaphore and deps.
                        stage_workflow.dependency_tasks = []
                        if prior_state != "failed":
                            stage_workflow.manager_status = "idle"

                        # Transfer from remote if necessary
                        upstream_idx = getattr(
                            stage_workflow.stage_list[0], "input_stage_index", 0)
                        stage_workflow.transfer_from = None
                        if upstream_idx and upstream_idx > 0:
                            upstream_sw = next(
                                (sw for sw in self.stage_workflows
                                 if sw.stage_list[0].index == upstream_idx),
                                None,
                            )
                            if (upstream_sw is not None
                                    and upstream_sw.resource.name
                                    != stage_workflow.resource.name):
                                stage_workflow.transfer_from = \
                                    upstream_sw.resource.name

                        #await stage_workflow.create_tasks()
                        sem = self.limit_semaphores_by_resource.get(
                            stage_workflow.resource.name)
                        if sem is not None and stage_workflow not in held_slots:
                            if sem.locked():
                                all_terminal = False
                                continue
                            await sem.acquire()
                            held_slots.append(stage_workflow)
                        try:
                            await stage_workflow.create_tasks()
                        except Exception:
                            release_slot(stage_workflow)
                            raise

                        if stage_workflow.task is not None:
                            self.stage_tasks.append(stage_workflow.task)
                        all_terminal = False
                        continue

                    if not stage_workflow.is_terminal():
                        all_terminal = False

                if all_terminal:
                    break

                #if all(
                #        sw.semaphore == "wait"
                #        and not sw.has_running_task_or_process()
                #        for sw in self.stage_workflows):
                #    return

                await asyncio.sleep(MAIN_LOOP_INTERVAL)

        finally:
            stop_event.set()
            for stage_workflow in list(held_slots):
                release_slot(stage_workflow)

        if (not self.detached_requested) and len(self.stage_tasks) > 0:
            await asyncio.gather(*self.stage_tasks, return_exceptions=True)

        if self.detached_requested:
            await self.detach_shutdown()
        elif self._shutdown_task is not None:
            await self._shutdown_task
        else:
            await self.shutdown()

async def assign_backend_workflow_engine_to_session(
        session: client_structures.RunSession,
    ) -> None:
    """
    Assign a backend and workflow engine to the session.
    """
    backend = await LocalExecutionBackend(ThreadPoolExecutor())
    workflow_engine = await WorkflowEngine.create(backend=backend)
    session.backend = backend
    session.workflow_engine = workflow_engine
    return

async def monitor_session_pipelines(
        pipelines: list[SeekrPipeline]
    ) -> None:
    """
    Monitor the session pipelines.
    """
    stop = False
    dead_resources: set[str] = set()
    while not stop:
        # Make status payload for all pipelines in the session
        ran_nothing = True
        status_payloads_by_resource = {}
        resource_dict = {}
        # Initialize the dicts first
        for pipeline in pipelines:
            system_name = pipeline.systemrun.name
            if pipeline.shutting_down:
                continue
            ran_nothing = False
            for sw in pipeline.stage_workflows:
                if sw.internal_id is None:
                    continue
                
                resource_name = sw.resource.name
                root_directory = workload_remote_local.resolve_model_directory(
                        sw.seekrflow, sw.resource)

                if resource_name not in resource_dict:
                    resource_dict[resource_name] = sw.resource
                if resource_name not in status_payloads_by_resource:
                    status_payloads_by_resource[resource_name] = {"system_payloads": {}}
                if system_name not in status_payloads_by_resource[resource_name]["system_payloads"]:
                    status_payloads_by_resource[resource_name]["system_payloads"][system_name] = {}
                    status_payloads_by_resource[resource_name]["system_payloads"][system_name]["root_dir"] \
                        = root_directory
                    status_payloads_by_resource[resource_name]["system_payloads"][system_name]["jobs"] = []

                job_dict = {}
                job_dict["internal_id"] = sw.internal_id
                job_dict["job_id"] = sw.job_id
                job_dict["job_name"] = sw.job_name
                job_dict["stage_indices"] = [stage.index for stage in sw.stage_list]
                job_dict["array_indices"] = sw.array_indices
                status_payloads_by_resource[resource_name]["system_payloads"]\
                    [system_name]["jobs"].append(job_dict)

        for resource_name, payload in status_payloads_by_resource.items():
            if resource_name in dead_resources:
                continue
            resource = resource_dict[resource_name]
            status = None
            try:
                status = await _run_blocking(
                    workload_remote_local.status,
                    resource,
                    manager_payload=payload,
                )
            except Exception as e:
                if isinstance(e, ImportError):
                    dead_resources.add(resource_name)
                    for pipeline in pipelines:
                        for sw in pipeline.stage_workflows:
                            if sw.resource is None or sw.resource.name != resource_name:
                                continue
                            if sw.semaphore != "go":
                                continue
                            sw.last_error = f"status failed: {e}"
                            sw.control_file_semaphore = "wait"
                            sw.semaphore = "wait"
                            sw.manager_status = "failed"
                            for stage in sw.stage_list:
                                sw.stage_state[stage.name] = "error"
                    print(
                        f"[remote-status] {resource_name}: {e}; setting semaphore=wait"
                    )
                    continue
                print(
                    f"[remote-status] {resource_name}: "
                    f"status call raised an exception: {e}"
                )
                await asyncio.sleep(STATUS_WRITE_INTERVAL)
                continue
            if status is not None:
                if not isinstance(status, dict) or status.get("success") is False or "payload" not in status:
                    print(f"[remote-status] {resource_name}: {status.get('error') if isinstance(status, dict) else status}")
                    await asyncio.sleep(STATUS_WRITE_INTERVAL)
                    continue
                return_payload = status["payload"]
                for pipeline in pipelines:
                    system_name = pipeline.systemrun.name
                    system_payload = return_payload.get(system_name, None)
                    if system_payload is None:
                        continue
                    for sw in pipeline.stage_workflows:
                        job_dict = system_payload.get(sw.job_id, None)
                        if job_dict is None:
                            continue
                        manager_dicts_by_array_index = job_dict["manager_dicts_by_array_index"]
                        manager_state_set = set()
                        for array_index, manager_dict in manager_dicts_by_array_index.items():
                            manager_state_set.add(manager_dict["state"])
                            
                        if len(manager_state_set) == 0:
                            sw.manager_status = "idle"
                        else:
                            if "running" in manager_state_set and (("queued" in manager_state_set) \
                                    or ("pending" in manager_state_set)):
                                sw.manager_status = "queued/running"
                            elif "running" in manager_state_set:
                                sw.manager_status = "running"
                            elif ("queued" in manager_state_set) or ("pending" in manager_state_set):
                                sw.manager_status = "queued"
                            elif "failed" in manager_state_set:
                                sw.manager_status = "failed"
                                newly_failed = sw.failed_job_id != sw.job_id
                                if newly_failed:
                                    sw.monitor_generation += 1
                                    sw.failed_job_id = sw.job_id
                                stage_dicts = job_dict.get("stage_dicts_by_stage_index") or {}

                                # Check if we made any progress since previous submission
                                made_progress = False
                                for stage in sw.stage_list:
                                    entry = stage_dicts.get(stage.index) or stage_dicts.get(str(stage.index)) or {}
                                    baseline = (sw.progress_at_submit or {}).get(stage.name) or {}
                                    for key_stage_index, key_anchor, key_swarm_id, value in (
                                            entry.get("stage_anchor_swarm_progress_list") or []):
                                        if key_stage_index != stage.index:
                                            continue
                                        key = f"anchor_{key_anchor}_swarm_{key_swarm_id}"
                                        if value > baseline.get(key, 0.0):
                                            made_progress = True
                                            break
                                    if made_progress:
                                        break
                                elapsed_values = [
                                    item.get("last_known_elapsed")
                                    for item in manager_dicts_by_array_index.values()
                                ]
                                run_seconds = _elapsed_seconds(elapsed_values)
                                # We don't want to count the queue time.
                                if run_seconds is None:
                                    run_seconds = 0
                                #if run_seconds is None and sw.submitted_at is not None:
                                #    run_seconds = time.time() - sw.submitted_at
                                resubmit = (
                                    sw.semaphore == "go"
                                    and made_progress
                                    and run_seconds is not None
                                    and run_seconds > 120
                                )
                                if newly_failed and not resubmit and sw.semaphore == "go":
                                    errors = [
                                        item.get("error") or item.get("notes")
                                        for item in manager_dicts_by_array_index.values()
                                        if item.get("error") or item.get("notes")
                                    ]
                                    sw.last_error = (
                                        errors[0] if errors
                                        else "job exited before completion"
                                    )
                                    #for stage in sw.stage_list:
                                    #    sw.stage_state[stage.name] = "error"
                                    sw.control_file_semaphore = "wait"
                                    sw.semaphore = "wait"
                                    print(
                                        f"[remote-job] stage {sw.stage_list[0].name} "
                                        f"failed; setting semaphore=wait. "
                                        f"error={sw.last_error}"
                                    )
                            else:
                                sw.manager_status = "idle"

                        stage_dicts_by_stage_index = job_dict["stage_dicts_by_stage_index"]
                        sw.apply_stage_dicts(stage_dicts_by_stage_index)
                        if sw.semaphore == "wait" and sw.manager_status == "failed":
                            for stage in sw.stage_list:
                                sw.stage_state[stage.name] = "error"

        if ran_nothing:
            stop = True
            return
        await asyncio.sleep(STATUS_WRITE_INTERVAL)
        
    return

async def status_writer_loop(
        output_file: str,
        pipelines: list[SeekrPipeline],
        ) -> None:
    """
    Periodically write the status snapshot until stop_event is set.
    """
    stop = False
    while not stop:
        session_statuses = {}
        all_shutting_down = True
        for pipeline in pipelines:
            pipeline_status = pipeline.status_snapshot()
            session_statuses[pipeline.systemrun.name] = pipeline_status
            if not pipeline_status.get("shutting_down", False):
                all_shutting_down = False
        if all_shutting_down:
            stop = True
            
        output_file_tmp = output_file + f".{os.getpid()}.tmp"
        with open(output_file_tmp, "w") as f:
            json.dump(session_statuses, f)
        os.rename(output_file_tmp, output_file)
        
        await asyncio.sleep(STATUS_WRITE_INTERVAL)

async def control_file_reader_loop(
        control_file: str,
        pipelines: list[SeekrPipeline],
    ) -> None:
    """
    Read the control file and update the pipelines.
    """
    stop = False
    lock_path = control_file + ".lock"
    while not stop:
        all_shutting_down = True
        made_changes = False
        with open(lock_path, "a+") as lockf:
            # Need to apply a lock to prevent the control file from being read 
            # and written to concurrently by the GUI program and the client.
            fcntl.flock(lockf, fcntl.LOCK_EX)
            try:
                with open(control_file, "r") as f:
                    control_dict = json.load(f)
                pending = []
                for pipeline in pipelines:
                    if pipeline.shutting_down:
                        continue
                    all_shutting_down = False
                    system_name = pipeline.systemrun.name
                    system_control = control_dict[system_name]
                    detach = system_control["detach"]
                    if detach:
                        pipeline.detached_requested = True
                    for sw in pipeline.stage_workflows:
                        stage_workflow_control = system_control["stage_workflows"][sw.stage_list[0].name]
                        transfer = stage_workflow_control["transfer"]
                        if transfer and not sw.transfer_requested and sw.transfer_status != "running":
                            sw.transfer_requested = True
                        if sw.control_file_semaphore is not None:
                            system_control["stage_workflows"][sw.stage_list[0].name]["semaphore"] \
                                = sw.control_file_semaphore
                            sw.semaphore = sw.control_file_semaphore
                            pending.append((sw, sw.control_file_semaphore))
                            made_changes = True
                        else:
                            sw.semaphore = stage_workflow_control["semaphore"]
                            if sw.semaphore == "stop":
                                transfer_base.cancel_running_transfers()
                    
                    control_dict[system_name] = system_control
                
                if made_changes:
                    tmp_file = control_file + f".{os.getpid()}.tmp"
                    with open(tmp_file, "w") as f:
                        json.dump(control_dict, f, indent=2)
                    os.rename(tmp_file, control_file)
                    for sw, value in pending:
                        if sw.control_file_semaphore == value:
                            sw.control_file_semaphore = None
            finally:
                fcntl.flock(lockf, fcntl.LOCK_UN)

        if all_shutting_down:
            stop = True
        await asyncio.sleep(CONTROL_FILE_READ_INTERVAL)
    return

def initialize_control_file(
        control_file: str,
        pipelines: list[SeekrPipeline],
    ) -> None:
    """
    Initialize the control file to have all "go" semaphores and no detach 
    or transfer upon client start.
    """
    control_dict = {}
    for pipeline in pipelines:
        system_name = pipeline.systemrun.name
        control_dict[system_name] = {"detach": False, "stage_workflows": {}}
        for sw in pipeline.stage_workflows:
            stage_workflow_control = {
                "semaphore": sw.semaphore,
                "transfer": False,
            }
            control_dict[system_name]["stage_workflows"][sw.stage_list[0].name] \
                = stage_workflow_control
    lock_path = control_file + ".lock"
    with open(lock_path, "a+") as lockf:
        fcntl.flock(lockf, fcntl.LOCK_EX)
        try:
            tmp_file = control_file + f".{os.getpid()}.tmp"
            with open(tmp_file, "w") as f:
                json.dump(control_dict, f, indent=2)
            os.rename(tmp_file, control_file)
        finally:
            fcntl.flock(lockf, fcntl.LOCK_UN)
    return

def restore_prior_jobs(
        workflows: list, 
        output_file: str
        ) -> None:
    """
    Restore job identity from the last status snapshot if a previous session was 
    interrupted.
    """
    try:
        with open(output_file) as f:
            session_statuses = json.load(f)
    except (OSError, json.JSONDecodeError):
        return
    by_stages = {}
    for system_name, system_status in session_statuses.items():
        for sw_status in (system_status.get("stage_workflows") or {}).values():
            names = frozenset((sw_status.get("stages") or {}))
            if names and sw_status.get("internal_id") is not None:
                by_stages[(system_name, names)] = sw_status
    for sw in workflows:
        names = frozenset(stage.name for stage in sw.stage_list)
        prior = by_stages.get((sw.seekrflow.name, names))
        if prior is None:
            continue
        sw.internal_id = prior["internal_id"]
        sw.job_id = prior.get("job_id")
        sw.job_name = prior.get("job_name")
        sw.array_indices = prior.get("array_indices") or []

async def run_pipeline_pool(
        pipelines: list[SeekrPipeline],
        #limit: int | None,
        ) -> None:
    """
    Run pipelines, with at most `limit` active. None runs them all.
    """
    #if limit is None or limit >= len(pipelines):
    await asyncio.gather(*(pipeline.run_workflows() for pipeline in pipelines))
    return
    '''pending = list(pipelines)
    async def worker() -> None:
        while pending:
            pipeline = pending.pop(0)
            if pipeline.shutting_down:
                continue
            await pipeline.run_workflows()
            if not pipeline.shutting_down:
                pending.append(pipeline)
                await asyncio.sleep(MAIN_LOOP_INTERVAL)
    await asyncio.gather(*(worker() for _ in range(limit)))
    '''

async def drop_dead_restored_jobs(
        pipelines: list[SeekrPipeline],
        ) -> None:
    """
    Forget a restored job that has already exited. Keep one that is still
    running, queued, or finished.
    """
    for pipeline in pipelines:
        for sw in pipeline.stage_workflows:
            if sw.internal_id is None:
                continue
            action = await sw.probe_launch()
            if action == "completed":
                for stage in sw.stage_list:
                    sw.stage_state[stage.name] = "completed"
            elif action == "submit":
                sw.internal_id = None
                sw.job_id = None
                sw.job_name = None
                sw.array_indices = []

async def launch_session_pipelines(
        session: client_structures.RunSession,
        ) -> None:
    """
    Launch a session pipeline.
    """
    await assign_backend_workflow_engine_to_session(session)
    pipelines = [
        SeekrPipeline(systemrun=systemrun, workflow_engine=session.workflow_engine)
        for systemrun in session.systemrun_objects
    ]
    initialize_control_file(session.control_file, pipelines)
    restore_prior_jobs(
        [sw for pipeline in pipelines for sw in pipeline.stage_workflows],
        session.output_file,
    )
    await drop_dead_restored_jobs(pipelines)
    # Propagate interrupt signals to the pipelines.
    def _on_interrupt() -> None:
        if any(pipeline.shutting_down for pipeline in pipelines):
            print("[shutdown] forcing immediate exit")
            os._exit(1)
        print("\n\n=== Received interrupt signal ===")
        transfer_base.cancel_running_transfers()
        for pipeline in pipelines:
            pipeline.shutting_down = True
            for stage_workflow in pipeline.stage_workflows:
                stage_workflow.semaphore = "stop"
    loop = asyncio.get_running_loop()
    for sig in (signal.SIGINT, signal.SIGTERM):
        loop.add_signal_handler(sig, _on_interrupt)
    # Run the pipelines.
    session.limit_semaphores_by_resource = {}
    for systemrun in session.systemrun_objects:
        for resource in systemrun.seekrflow.run_settings.resources:
            if resource.max_concurrent_systems is None:
                continue
            if resource.name not in session.limit_semaphores_by_resource:
                asyncio_semaphore = asyncio.Semaphore(resource.max_concurrent_systems)
                session.limit_semaphores_by_resource[resource.name] = asyncio_semaphore
                for pipeline in pipelines:
                    pipeline.limit_semaphores_by_resource[resource.name] = asyncio_semaphore
            else:
                asyncio_semaphore = session.limit_semaphores_by_resource[resource.name]
    try:
        await asyncio.gather(
            control_file_reader_loop(session.control_file, pipelines),
            status_writer_loop(session.output_file, pipelines),
            monitor_session_pipelines(pipelines),
            run_pipeline_pool(pipelines) #, limit),
        )
    finally:
        await session.workflow_engine.shutdown()
    return
