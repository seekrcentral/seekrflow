"""
modules/workload_managers/remote.py

An intermediate module between the seekrflow runner and the remote workload
managers (e.g., SLURM). This module handles common functionality needed for
remote workload management.
"""

import os
import typing

import seekr2.modules.common_base as seekr2_base

import seekrflow.modules.structures as structures
import seekrflow.modules.workload_managers.slurm as workload_slurm
import seekrflow.modules.workload_managers.pbs as workload_pbs
import seekrflow.modules.remote_interfaces.globus_compute_sdk as remote_globus
import seekrflow.modules.remote_interfaces.ssh as remote_ssh


def _get_stage_resource_name(
        seekrflow: structures.Seekrflow,
        stage_name: str,
        ) -> str:
    """
    Resolve the configured resource name for a stage using the canonical
    run_settings.stage_resource_name mapping.
    """
    if stage_name not in seekrflow.run_settings.stage_resource_name:
        raise ValueError(
            f"Unknown stage name {stage_name!r} in run_settings.stage_resource_name"
        )
    return seekrflow.run_settings.stage_resource_name[stage_name]


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
    Calculate optimal time limit for a stage based on workload manager type.

    This function dispatches to the appropriate workload manager's optimization
    function based on the resource type. Currently supports SLURM with extensibility
    for future workload managers (PBS, LSF, etc.).

    Args:
        seekrflow: Seekrflow object containing configuration
        stage_name: Name of the stage (e.g., "seekr", "hidr", "bd")
        resource: Remote resource configuration object
        incomplete_anchors: List of anchor indices that need more simulation
        job_status_filename: Path to .seekrflow_job_status.json file

    Returns:
        Optimized time limit string in HH:MM:SS format, or original time limit
        if optimization is not available/applicable for this workload manager
    """
    if resource.type == "slurm_remote":
        # SLURM-specific time limit optimization
        if stage_name == "seekr":
            return workload_slurm.calculate_optimal_seekr_time_limit(
                job_status_filename,
                incomplete_anchors,
                resource.time_limit
            )
        else:
            # Currently only SEEKR stage is optimized
            return resource.time_limit
    elif resource.type == "pbs_remote":
        # PBS-specific time limit optimization
        if stage_name == "seekr":
            return workload_pbs.calculate_optimal_seekr_time_limit(
                job_status_filename,
                incomplete_anchors,
                resource.time_limit
            )
        else:
            # Currently only SEEKR stage is optimized
            return resource.time_limit
    else:
        # Unknown or unsupported workload manager type
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
    else:
        raise NotImplementedError(
            f"Remote interface type not implemented: "\
            f"{resource.remote_interface.type}")
    
    return result

def status_remote(
        seekrflow: structures.Seekrflow,
        stage_name: str,
        model: seekr2_base.Model,
        silent: bool = False,
        benchmark_mode: bool = False,
        ) -> dict:
    """
    Generalized remote stage status query for any stage configured in
    run_settings.stage_resource_name.
    """
    resource_name = _get_stage_resource_name(seekrflow, stage_name)
    if resource_name == "local":
        raise ValueError(
            f"Stage {stage_name!r} is configured for local execution; "
            f"status_remote is for remote stages only."
        )

    resource = seekrflow.run_settings.get_resource_by_name(resource_name)

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
    elif resource.type == "pbs_remote":
        status_workflow = workload_pbs.pbs_remote_status_workflow
    else:
        raise NotImplementedError(
            f"Resource type {resource.type} not implemented"
        )

    result = submit_remote_workflow(
        seekrflow,
        resource_name,
        status_workflow,
        extra_args=extra_args,
        silent=silent,
    )
    # Normalize stage_status regardless of success so callers always have
    # manager_status / job list available (e.g. for display while seekr module
    # is unavailable in the status-check environment).
    result["stage_status"] = _normalize_stage_status(result)
    return result

def submit_remote_run_workflow(
        seekrflow: structures.Seekrflow,
        stage_name: str,
        destination_path: str,
        resource: structures.Resource_remote_base,
        command_string: str,
        model_filename: str,
        workflow_type: str,
        indices: list | None = None,
        silent: bool = False,
        anchor_times: dict | None = None,
        time_limit_override: str | None = None,
    ) -> list:
    """
    Run the stage remotely. Submit a Globus Compute or SSH workflow to 
    do this.

    Parameters
    ----------
    time_limit_override : str or None
        If provided, this walltime (HH:MM:SS) is submitted to the workload
        manager in place of ``resource.time_limit``.  Used for benchmark-mode
        runs and for dynamic time-limit optimization on subsequent SEEKR
        submissions.  ``resource.time_limit`` is never mutated.
    """
    log_directory = os.path.join(destination_path, "logs")
    name = seekrflow.name + "_" + stage_name
    effective_time_limit = time_limit_override or resource.time_limit
    if resource.type == "slurm_remote":
        args = [
            stage_name,
            log_directory,
            resource.partition,
            resource.account,
            resource.constraint,
            name,
            resource.cpus_per_task,
            resource.memory_per_node,
            effective_time_limit,
            resource.scheduler_options,
            resource.worker_init,
            command_string,
            indices,
            model_filename,
            workflow_type,
            getattr(resource, "mps", 1),
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
            resource.cpus_per_task,
            resource.memory_per_node,
            effective_time_limit,
            resource.scheduler_options,
            resource.worker_init,
            command_string,
            indices,
            model_filename,
            workflow_type,
            getattr(resource, "mps", 1),
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
        job_id: str,
        silent: bool = False
    ) -> None:
    """
    Cancel the remote stage. Submit a Globus Compute or SSH workflow to
    signal the cancellation.
    """
    args = [job_id]
    resource = seekrflow.run_settings.get_resource_by_name(resource_name)

    if resource.type == "slurm_remote":
        cancel_workflow = workload_slurm.slurm_remote_cancel_workflow
    elif resource.type == "pbs_remote":
        cancel_workflow = workload_pbs.pbs_remote_cancel_workflow
    else:
        raise NotImplementedError(f"Resource type {resource.type} not implemented")

    submit_remote_workflow(seekrflow, resource_name,
                           cancel_workflow, extra_args=args, silent=silent)
    return

def force_rerun_stage_remote(
        seekrflow: structures.Seekrflow,
        stage_name: str,
        silent: bool = False
        ) -> None:
    """
    Force re-run a stage on a remote resource by canceling jobs and deleting files.
    
    Parameters
    ----------
    seekrflow : structures.Seekrflow
        The seekrflow configuration
    stage_name : str
        Name of the stage ("bd", "hidr", or "seekr")
    silent : bool
        Whether to suppress output
    """
    resource_name = _get_stage_resource_name(seekrflow, stage_name)
    
    if resource_name == "local":
        raise ValueError(f"Stage {stage_name} is configured to run locally, not remotely")

    print(f"Forcing re-run of remote stage: {stage_name} on resource {resource_name}")

    resource = seekrflow.run_settings.get_resource_by_name(resource_name)

    if resource.type == "slurm_remote":
        force_rerun_workflow = workload_slurm.slurm_remote_force_rerun_workflow
    elif resource.type == "pbs_remote":
        force_rerun_workflow = workload_pbs.pbs_remote_force_rerun_workflow
    else:
        raise NotImplementedError(f"Resource type {resource.type} not implemented")

    result = submit_remote_workflow(
        seekrflow, resource_name, force_rerun_workflow,
        extra_args=[stage_name], silent=silent)

    if result.get("canceled_job"):
        print(f"  Canceled job: {result['canceled_job']}")

    print(f"  Deleted files for {stage_name} stage")

    return
