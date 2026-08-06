"""
modules/workload_managers/remote.py

An intermediate module between the seekrflow runner and the remote workload
managers (e.g., SLURM). This module handles common functionality needed for
remote workload management.
"""

import os
import typing

import seekrflow.modules.structures as structures
import seekrflow.modules.workload_managers.slurm as workload_slurm
import seekrflow.modules.workload_managers.pbs as workload_pbs
import seekrflow.modules.workload_managers.aws as workload_aws
import seekrflow.modules.workload_managers.dispatch_lowering as dispatch_lowering
import seekrflow.modules.remote_interfaces.globus_compute_sdk as remote_globus
import seekrflow.modules.remote_interfaces.ssh as remote_ssh
import seekrflow.modules.remote_interfaces.local_shell as remote_local_shell


def resource_kind(resource: structures.Resource_base | None) -> str:
    """Return ``local``, ``remote``, or ``cloud`` for a resource object."""
    if resource is None:
        return "local"
    if resource.type in ("slurm_remote", "pbs_remote"):
        return "remote"
    if resource.type == "aws_cloud":
        return "cloud"
    return "local"

def _get_stage_resource_name(
        seekrflow: structures.Seekrflow,
        stage_name: str,
        ) -> str:
    """
    Resolve the configured resource name for a stage via placement policies.
    """
    return seekrflow.run_settings.resolve_stage_execution(
        stage_name, seekrflow.workflow.procedure).resource_name


def _normalize_stage_status(
        status_result: dict,
        ) -> dict:
    """
    Normalize stage status payload to align with local status contract:
    finished(bool), state(str), progress(float), notes(str).
    """
    stage_status = dict(status_result.get("stage_status") or {})
    if len(stage_status) == 0:
        stage_status = {
            "finished": False,
            "state": "unknown",
            "progress": 0.0,
            "notes": "",
        }

    finished = bool(stage_status.get("finished", False))
    notes = str(stage_status.get("notes", ""))

    progress_raw = stage_status.get("progress", 0.0)
    progress: float
    if isinstance(progress_raw, (int, float)):
        progress = float(progress_raw)
    elif isinstance(progress_raw, str):
        lowered = progress_raw.strip().lower()
        if lowered in {"completed", "complete", "finished"}:
            progress = 1.0
        elif lowered in {"started", "running", "queued"}:
            progress = 0.0
        else:
            progress = 0.0
    else:
        progress = 0.0

    state = stage_status.get("state")
    if state is None:
        if finished:
            state = "completed"
        else:
            manager_status = status_result.get("manager_status") or {}
            jobs = manager_status.get("jobs") or []
            if len(jobs) > 0:
                state = "started"
            else:
                state = "unstarted"

    stage_status["finished"] = finished
    stage_status["state"] = state
    stage_status["progress"] = progress
    stage_status["notes"] = notes
    return stage_status


def calculate_optimal_time_limit(
        seekrflow: structures.Seekrflow,
        stage_name: str,
        resource: structures.Resource_remote_base,
        incomplete_anchors: list[int],
        job_status_filename: str
) -> str:
    """
    Deprecated legacy entry point.

    Dynamic walltime now lives in ``walltime.estimate_submit_time_limit`` and
    is applied at submit via ``Resolved_execution.time_policy``. This wrapper
    returns the Resource ``time_limit`` unchanged.
    """
    del seekrflow, stage_name, incomplete_anchors, job_status_filename
    return resource.time_limit


def submit_remote_workflow(
        seekrflow: structures.Seekrflow,
        resource_name: str,
        workflow: typing.Any,
        extra_args: list | None = None,
        silent: bool = False
        ) -> dict:
    """
    Check if the SEEKR calculation has finished remotely. Submit a Globus Compute
    or SSH workflow to obtain this.
    """
    resource: structures.Resource_remote_base \
        = seekrflow.run_settings.get_resource_by_name(resource_name)
    destination_path = os.path.join(
        resource.remote_working_directory, seekrflow.name)
    if resource.type == "slurm_remote":
        args = [destination_path]
        if extra_args is not None:
            args.extend(extra_args)
    elif resource.type == "pbs_remote":
        args = [destination_path]
        if extra_args is not None:
            args.extend(extra_args)
    else:
        raise NotImplementedError(
            f"Resource type {resource.type} is not implemented.")
    
    if resource.remote_interface.type == "globus_compute_sdk":
        endpoint = resource.remote_interface.endpoint_id
        result = remote_globus\
            .submit_remote_workflow_with_globus_compute(
                resource.name, workflow, endpoint, args, silent)
    elif resource.remote_interface.type == "ssh":
        hostname = resource.remote_interface.hostname
        username = resource.remote_interface.username
        password = resource.remote_interface.password
        port = resource.remote_interface.port
        private_key_filename = resource.remote_interface.private_key_filename
        private_key_passphrase = resource.remote_interface.private_key_passphrase
        result = remote_ssh\
            .submit_remote_workflow_with_ssh(
                resource.name, workflow, args,  
                hostname, username, password, port,
                private_key_filename, private_key_passphrase)
    elif resource.remote_interface.type == "local_shell":
        result = remote_local_shell\
            .submit_remote_workflow_with_local_shell(
                resource.name, workflow, args,
                python_executable=resource.remote_interface.python_executable,
                silent=silent)
    else:
        raise NotImplementedError(
            f"Remote interface type not implemented: "\
            f"{resource.remote_interface.type}")
    
    return result

def status_remote(
        seekrflow: structures.Seekrflow,
        stage_name: str,
        model: typing.Any,
        silent: bool = False,
        benchmark_mode: bool = False,
        host_stage_name: str | None = None,
        announce_failure: bool = True,
        ) -> dict:
    """
    Generalized remote/cloud stage status query for any stage whose producing
    procedure is configured in run_settings.placements.

    For fused AWS members, pass ``host_stage_name`` so Batch liveness is read
    from the host job while Seekr progress still comes from this stage's
    S3 status object.

    ``announce_failure=False`` suppresses console Batch failure dumps (used by
    quiet pre-submit probes that may still see a stale FAILED job id).
    """
    resource_name = _get_stage_resource_name(seekrflow, stage_name)
    if resource_name == "local":
        raise ValueError(
            f"Stage {stage_name!r} is configured for local execution; "
            f"status_remote is for remote stages only."
        )

    resource = seekrflow.run_settings.get_resource_by_name(resource_name)
    kind = resource_kind(resource)

    if kind == "cloud":
        model_directory = getattr(model, "directory", None) or "."
        result = workload_aws.status_aws(
            seekrflow,
            stage_name,
            resource,
            model_directory,
            host_stage_name=host_stage_name,
            benchmark_mode=benchmark_mode,
            announce_failure=announce_failure,
        )
    elif kind == "remote":
        # Note: positional args slot — the SLURM status workflow expects
        # args[0]=working_dir, [1]=stage_name, [2]=benchmark_mode, [3]=anchor,
        # [4]=swarm_id, [5]=worker_init.  worker_init is shelled into a bash -lc
        # invocation inside the workflow so it can activate the env that has seekr.
        extra_args: list[typing.Any] = [
            stage_name,
            benchmark_mode,
            "any",
            None,
            getattr(resource, "worker_init", "") or "",
        ]
        if resource.type == "slurm_remote":
            status_workflow = workload_slurm.slurm_remote_status_workflow
        else:
            status_workflow = workload_pbs.pbs_remote_status_workflow
        result = submit_remote_workflow(
            seekrflow,
            resource_name,
            status_workflow,
            extra_args=extra_args,
            silent=silent,
        )
    else:
        raise NotImplementedError(
            f"Resource type {resource.type} not implemented"
        )

    # Normalize stage_status regardless of success so callers always have
    # manager_status / job list available (e.g. for display while seekr module
    # is unavailable in the status-check environment).
    result["stage_status"] = _normalize_stage_status(result)
    return result

def remote_run_python_snippet_workflow(args):
    """Run a python snippet on the remote worker; return raw stdout/stderr.
    Args: [working_dir, snippet, worker_init]. Stdlib only (no seekrflow)."""
    import shlex, pathlib, subprocess
    working_dir, snippet = args[0], args[1]
    worker_init = args[2] if len(args) > 2 else ""
    init = (worker_init or "").strip()
    cmd = (f"{init}; python -c {shlex.quote(snippet)}" if init
           else f"python -c {shlex.quote(snippet)}")
    try:
        proc = subprocess.run(["bash", "-lc", cmd],
                              cwd=str(pathlib.Path(working_dir)),
                              capture_output=True, text=True, timeout=120)
        return {"success": True, "error": None, "stdout": proc.stdout,
                "stderr": proc.stderr, "returncode": proc.returncode}
    except Exception as e:
        return {"success": False, "error": str(e), "stdout": "", "stderr": ""}


def fetch_unit_counts_remote(
        seekrflow: structures.Seekrflow,
        launching_stage: typing.Any,
        resource: structures.Resource_remote_base,
        silent: bool = False,
        ) -> dispatch_lowering.StageUnitCounts:
    """
    Query seekr ``info`` on a remote worker for launch-time unit enumeration.
    """
    if resource.type not in ("slurm_remote", "pbs_remote"):
        raise NotImplementedError(
            f"Resource type {resource.type} is not implemented.")

    input_stage_index = getattr(launching_stage, "input_stage_index", 0)
    launching_stage_index = getattr(
        launching_stage, "index", launching_stage.name)
    snippet = dispatch_lowering.build_info_fetch_snippet(
        "model.json",
        launching_stage.name,
        input_stage_index,
        launching_stage_index,
    )
    extra_args = [
        snippet,
        getattr(resource, "worker_init", "") or "",
    ]

    result = submit_remote_workflow(
        seekrflow,
        resource.name,
        remote_run_python_snippet_workflow,
        extra_args=extra_args,
        silent=silent,
    )
    if not result.get("success"):
        raise RuntimeError(
            f"Remote unit-count fetch failed: {result.get('error')}")
    return dispatch_lowering.parse_info_fetch_output(
        result.get("stdout", ""))


def submit_remote_run_workflow(
        seekrflow: structures.Seekrflow,
        stage_name: str,
        destination_path: str,
        resource: structures.Resource_base,
        command_string: str,
        model_filename: str,
        workflow_type: str,
        indices: list | None = None,
        silent: bool = False,
        anchor_times: dict | None = None,
        time_limit_override: str | None = None,
        stage_specs: list[dict] | None = None,
        model_directory: str | None = None,
        resolved_execution: structures.Resolved_execution | None = None,
    ) -> dict:
    """
    Run the stage on a remote or cloud resource.

    Compute sizing prefers ``resolved_execution`` (Placement overrides +
    Resource defaults). ``time_limit_override`` still wins for walltime
    (benchmark / adaptive optimizer).

    For ``aws_cloud``, ``stage_specs`` and ``model_directory`` are required and
    ``command_string`` is ignored (the container command is built in aws.py).
    """
    kind = resource_kind(resource)
    if kind == "cloud":
        if not stage_specs:
            raise ValueError(
                "stage_specs is required for aws_cloud submit")
        if model_directory is None:
            raise ValueError(
                "model_directory is required for aws_cloud submit")
        timeout_seconds = None
        if time_limit_override is not None:
            timeout_seconds = structures.time_limit_to_seconds(
                time_limit_override)
        elif (resolved_execution is not None
                and resolved_execution.time_limit is not None):
            timeout_seconds = structures.time_limit_to_seconds(
                resolved_execution.time_limit)
        cpus = (
            resolved_execution.cpus
            if resolved_execution is not None else None)
        memory_mb = (
            resolved_execution.memory_mb
            if resolved_execution is not None else None)
        array_size = None
        if indices is not None and len(indices) > 1:
            array_size = len(indices)
        return workload_aws.submit_aws_job(
            seekrflow,
            stage_name,
            resource,
            stage_specs,
            model_directory,
            timeout_seconds=timeout_seconds,
            cpus=cpus,
            memory_mb=memory_mb,
            array_size=array_size,
        )

    if resolved_execution is not None:
        cpus = (
            resolved_execution.cpus
            if resolved_execution.cpus is not None
            else resource.cpus_per_task)
        memory_mb = (
            resolved_execution.memory_mb
            if resolved_execution.memory_mb is not None
            else resource.memory_per_node)
        time_limit = (
            resolved_execution.time_limit
            if resolved_execution.time_limit is not None
            else resource.time_limit)
    else:
        cpus = resource.cpus_per_task
        memory_mb = resource.memory_per_node
        time_limit = resource.time_limit

    log_directory = os.path.join(destination_path, "logs")
    name = seekrflow.name + "_" + stage_name
    effective_time_limit = time_limit_override or time_limit
    if resource.type == "slurm_remote":
        args = [
            stage_name,
            log_directory,
            resource.partition,
            resource.account,
            resource.constraint,
            name,
            cpus,
            memory_mb,
            effective_time_limit,
            resource.scheduler_options,
            resource.worker_init,
            command_string,
            indices,
            model_filename,
            workflow_type,
            anchor_times
        ]
        run_workflow = workload_slurm.slurm_remote_run_workflow
    elif resource.type == "pbs_remote":
        args = [
            stage_name,
            log_directory,
            resource.queue,              # PBS equivalent of partition
            resource.account,
            resource.resource_list,      # PBS equivalent of constraint
            name,
            cpus,
            memory_mb,
            effective_time_limit,
            resource.scheduler_options,
            resource.worker_init,
            command_string,
            indices,
            model_filename,
            workflow_type,
            anchor_times
        ]
        run_workflow = workload_pbs.pbs_remote_run_workflow
    else:
        raise NotImplementedError(
            f"Resource type {resource.type} is not implemented.")

    result = submit_remote_workflow(
        seekrflow, resource.name, run_workflow, extra_args=args, silent=silent)
    return result


def submit_remote_cancel_workflow(
        seekrflow: structures.Seekrflow,
        resource_name: str,
        job_id: str | None = None,
        job_name: str | None = None,
        silent: bool = False
    ) -> None:
    """
    Cancel a remote or cloud stage job.

    Cancel by ``job_id``, by scheduler ``job_name``, or both when provided.
    """
    if job_id is None and job_name is None:
        raise ValueError("job_id or job_name is required for remote cancel")
    resource = seekrflow.run_settings.get_resource_by_name(resource_name)
    kind = resource_kind(resource)
    if kind == "cloud":
        workload_aws.cancel_aws_job(
            resource, job_id=job_id, job_name=job_name)
        return

    extra_args: list[str] = []
    if job_id is not None:
        extra_args.append(job_id)
    if job_name is not None:
        if job_id is None:
            extra_args = ["", job_name]
        else:
            extra_args.append(job_name)

    if resource.type == "slurm_remote":
        cancel_workflow = workload_slurm.slurm_remote_cancel_workflow
    elif resource.type == "pbs_remote":
        cancel_workflow = workload_pbs.pbs_remote_cancel_workflow
    else:
        raise NotImplementedError(f"Resource type {resource.type} not implemented")

    submit_remote_workflow(seekrflow, resource_name,
                           cancel_workflow, extra_args=extra_args, silent=silent)
    return


def cancel_and_reset_remote_stage(
        seekrflow: structures.Seekrflow,
        stage_name: str,
        silent: bool = False,
        model_directory: str | None = None,
        monitor_stage_names: list[str] | None = None,
        ) -> None:
    """
    Cancel a stage's remote/cloud job and reset its scheduler/runner bookkeeping.

    Cancels any running job for the stage, waits until it is fully gone, and
    clears the workload-manager runner state so the stage resubmits fresh.
    Stage output cleanup is handled separately by seekr's own
    ``force_overwrite`` in the run command, so this does not delete trajectory
    data. For ``aws_cloud``, S3 monitor JSON (``*_status.json`` /
    ``*_dispatch.json``) is removed so a force-rerun is not treated as already
    completed; pass ``monitor_stage_names`` to clear fused members too.
    """
    resource_name = _get_stage_resource_name(seekrflow, stage_name)
    
    if resource_name == "local":
        raise ValueError(f"Stage {stage_name} is configured to run locally, not remotely")

    print(f"Canceling and resetting remote stage: {stage_name} on resource {resource_name}")

    resource = seekrflow.run_settings.get_resource_by_name(resource_name)
    kind = resource_kind(resource)

    if kind == "cloud":
        if model_directory is None:
            raise ValueError(
                "model_directory is required to reset an aws_cloud stage")
        result = workload_aws.cancel_and_reset_aws_stage(
            seekrflow,
            stage_name,
            resource,
            model_directory,
            monitor_stage_names=monitor_stage_names,
        )
        if result.get("canceled_job"):
            print(f"  Canceled job: {result['canceled_job']}")
        cleared = result.get("cleared_s3_keys") or []
        if cleared:
            print(f"  Cleared {len(cleared)} S3 monitor object(s) for force-rerun")
        print(f"  Cleared runner state for {stage_name} stage")
        return

    if resource.type == "slurm_remote":
        cancel_reset_workflow = workload_slurm.slurm_remote_force_rerun_workflow
    elif resource.type == "pbs_remote":
        cancel_reset_workflow = workload_pbs.pbs_remote_force_rerun_workflow
    else:
        raise NotImplementedError(f"Resource type {resource.type} not implemented")

    result = submit_remote_workflow(
        seekrflow, resource_name, cancel_reset_workflow,
        extra_args=[stage_name], silent=silent)

    if result.get("canceled_job"):
        print(f"  Canceled job: {result['canceled_job']}")

    print(f"  Cleared runner state for {stage_name} stage")

    return
