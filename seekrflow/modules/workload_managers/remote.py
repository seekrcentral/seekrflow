"""
modules/workload_managers/remote.py

An intermediate module between the seekrflow runner and the remote workload
managers (e.g., SLURM). This module handles common functionality needed for
remote workload management.
"""

import os
import typing

import seekr2.modules.common_base as seekr2_base

import seekrflow.modules.base as base
import seekrflow.modules.structures as structures
import seekrflow.modules.workload_managers.slurm as workload_slurm
import seekrflow.modules.remote_interfaces.globus_compute_sdk as remote_globus
import seekrflow.modules.remote_interfaces.ssh as remote_ssh

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
    else:
        raise NotImplementedError(
            f"Resource type {resource.type} is not implemented.")
    
    if resource.remote_interface.type == "globus_compute_sdk":
        endpoint = resource.remote_interface.endpoint_id
        result = remote_globus\
            .submit_remote_workflow_with_globus_compute(
                resource.name, workflow, endpoint, args, silent)
        if not result["success"]:
            print("result:", result)
        assert result["success"], "Workflow failed on remote resource."
    elif resource.remote_interface.type == "ssh":
        hostname = resource.remote_interface.hostname
        username = resource.remote_interface.username
        password = resource.remote_interface.password
        port = resource.remote_interface.port
        #private_key_filename = resource.remote_interface.private_key_filename
        #private_key_passphrase = resource.remote_interface.private_key_passphrase
        result = remote_ssh\
            .submit_remote_workflow_with_ssh(
                resource.name, workflow, args,  
                hostname, username, password, port)
    else:
        raise NotImplementedError(
            f"Remote interface type not implemented: "\
            f"{resource.remote_interface.type}")
    
    return result

def bd_status_remote(
        seekrflow: structures.Seekrflow,
        model: seekr2_base.Model,
        silent: bool = False
        ) -> str:
    """
    Check if the BD stage has finished remotely. Submit a Globus Compute
    or SSH workflow to obtain this.
    """
    assert model.k_on_info is not None
    n_traj_for_completed = model.k_on_info.b_surface_num_trajectories
    bd_status_workflow = workload_slurm.slurm_remote_bd_status_workflow
    result = submit_remote_workflow(seekrflow, seekrflow.run_settings.bd_stage_resource_name, 
                           bd_status_workflow, extra_args=[n_traj_for_completed], silent=silent)
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
        mps: int = 1,
        silent: bool = False
    ) -> list:
    """
    Run the stage remotely. Submit a Globus Compute or SSH workflow to 
    do this.
    """
    log_directory = os.path.join(destination_path, "logs")
    name = seekrflow.name + "_" + stage_name
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
            resource.time_limit,
            resource.scheduler_options,
            resource.worker_init,
            command_string,
            indices,
            model_filename,
            workflow_type,
            mps
        ]
        run_workflow = workload_slurm.slurm_remote_run_workflow
        
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
    cancel_workflow = workload_slurm.slurm_remote_cancel_workflow
    submit_remote_workflow(seekrflow, resource_name, 
                           cancel_workflow, extra_args=args, silent=silent)
    return

def hidr_status_remote(
        seekrflow: structures.Seekrflow,
        model: seekr2_base.Model,
        silent: bool = False
        ) -> str:
    """
    Check if the HIDR stage has finished remotely. Submit a Globus Compute
    or SSH workflow to obtain this.
    """
    hidr_status_workflow = workload_slurm.slurm_remote_hidr_status_workflow
    result = submit_remote_workflow(seekrflow, seekrflow.run_settings.hidr_stage_resource_name,
                                    hidr_status_workflow, silent=silent)
    return result

def seekr_status_remote(
        seekrflow: structures.Seekrflow,
        model: seekr2_base.Model,
        silent: bool = False,
        benchmark_mode: bool = False
        ) -> str:
    """
    Check if the HIDR stage has finished remotely. Submit a Globus Compute
    or SSH workflow to obtain this.
    """
    if benchmark_mode:
        time_for_completed = base.BENCHMARK_MIN_SIMULATION_LENGTH * model.get_timestep()
    else:
        time_for_completed = model.calculation_settings.num_production_steps \
            * model.get_timestep()
    seekr_status_workflow = workload_slurm.slurm_remote_seekr_status_workflow
    result = submit_remote_workflow(
        seekrflow, seekrflow.run_settings.seekr_stage_resource_name, 
        seekr_status_workflow, extra_args=[time_for_completed], 
        silent=silent)
    return result

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
    # Determine which resource this stage runs on
    if stage_name == "bd":
        resource_name = seekrflow.run_settings.bd_stage_resource_name
    elif stage_name == "hidr":
        resource_name = seekrflow.run_settings.hidr_stage_resource_name
    elif stage_name == "seekr":
        resource_name = seekrflow.run_settings.seekr_stage_resource_name
    else:
        raise ValueError(f"Unknown stage name: {stage_name}")
    
    if resource_name == "local":
        raise ValueError(f"Stage {stage_name} is configured to run locally, not remotely")
    
    print(f"Forcing re-run of remote stage: {stage_name} on resource {resource_name}")
    
    force_rerun_workflow = workload_slurm.slurm_remote_force_rerun_workflow
    result = submit_remote_workflow(
        seekrflow, resource_name, force_rerun_workflow, 
        extra_args=[stage_name], silent=silent)
    
    if result.get("canceled_job"):
        print(f"  Canceled job: {result['canceled_job']}")
    
    print(f"  Deleted files for {stage_name} stage")
    
    return
