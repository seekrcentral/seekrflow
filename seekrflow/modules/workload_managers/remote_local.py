"""
modules/workload_managers/remote_local.py

An intermediate module between the seekrflow runner and the 
remote/local workload managers (e.g., SLURM/PBS/AWS/local). 
This module handles common functionality needed for
remote/local workload management.
"""

from math import e
import os
import json
import typing
import pathlib
from dataclasses import dataclass

import seekrflow.modules.structures as structures
import seekrflow.modules.client.validation as client_validation
import seekrflow.modules.workload_managers.direct as workload_direct
import seekrflow.modules.workload_managers.slurm_pbs as workload_slurm_pbs
import seekrflow.modules.workload_managers.aws as workload_aws
import seekrflow.modules.remote_interfaces.globus_compute_sdk as remote_globus
import seekrflow.modules.remote_interfaces.ssh as remote_ssh
import seekrflow.modules.remote_interfaces.local_shell as remote_local_shell

@dataclass(frozen=True)
class StageUnitCounts:
    """Counts used to size a stage launch."""
    scope: str
    num_anchors: int | None
    num_swarms: int

def resource_kind(resource: structures.Resource_base | None) -> str:
    """Return ``local``, ``remote``, or ``cloud`` for a resource object."""
    if resource is None:
        return "local"
    if resource.type in ("slurm_remote", "pbs_remote"):
        return "remote"
    if resource.type in ("aws_cloud", ):
        return "cloud"
    return "local"

def resolve_model_directory(
        seekrflow: structures.Seekrflow,
        resource: structures.Resource_base,
        ) -> str:
    """
    Return the model root used by remote SLURM/PBS workflows or local_shell.
    """
    if hasattr(resource, "remote_interface"):
        if resource.remote_interface.type == "local_shell":
            return str(seekrflow.get_root_directory())
    if hasattr(resource, "remote_working_directory"):
        return os.path.join(resource.remote_working_directory, seekrflow.name)
    else:
        return seekrflow.name

def submit_workload(
        resource: structures.Resource_remote_base,
        workload: typing.Any,
        payload: dict,
        silent: bool = False,
        ) -> dict:
    """
    Submit a workflow to a remote resource.
    """
    # TODO: figure out what to put into 'args'
    if resource.remote_interface.type == "globus_compute_sdk":
        endpoint = resource.remote_interface.endpoint_id
        result = remote_globus\
            .submit_remote_workload_with_globus_compute(
                resource.name, workload, endpoint, payload, silent)
    elif resource.remote_interface.type == "ssh":
        hostname = resource.remote_interface.hostname
        username = resource.remote_interface.username
        password = resource.remote_interface.password
        port = resource.remote_interface.port
        private_key_filename = resource.remote_interface.private_key_filename
        private_key_passphrase = resource.remote_interface.private_key_passphrase
        result = remote_ssh\
            .submit_remote_workload_with_ssh(
                resource.name, workload, payload,  
                hostname, username, password, port,
                private_key_filename, private_key_passphrase)
    # TODO: add ORBIT
    elif resource.remote_interface.type == "local_shell":
        result = remote_local_shell\
            .submit_remote_workflow_with_local_shell(
                resource.name, workload, payload,
                python_executable=resource.remote_interface.python_executable,
                silent=silent)
    else:
        raise NotImplementedError(
            f"Remote interface type not implemented: "\
            f"{resource.remote_interface.type}")
    
    return result

# NOTE: used in fetch_unit_counts
def run_python_snippet_workload(payload: dict):
    """
    Run a python snippet on the remote worker - return raw stdout/stderr
    """
    import shlex, pathlib, subprocess
    root_dir = payload["root_dir"]
    snippet = payload["snippet"]
    worker_init = payload.get("worker_init", "")
    init = (worker_init or "").strip()
    cmd = (f"{init}; python -c {shlex.quote(snippet)}" if init
           else f"python -c {shlex.quote(snippet)}")
    try:
        proc = subprocess.run(["bash", "-lc", cmd],
                              cwd=str(pathlib.Path(root_dir)),
                              capture_output=True, text=True, timeout=90)
        return {"success": True, "error": None, "stdout": proc.stdout,
                "stderr": proc.stderr, "returncode": proc.returncode}
    except subprocess.TimeoutExpired as e:
        return {
            "success": False, 
            "error": str(e), 
            "stdout": e.stdout or "", 
            "stderr": e.stderr or "",
            "returncode": None,
            }
    except Exception as e:
        return {
            "success": False, 
            "error": str(e), 
            "stdout": "", 
            "stderr": "",
            "returncode": None,
            }

def parse_info_fetch_output(stdout: str) -> StageUnitCounts:
    """
    Parse the output of a seekr info fetch and return the unit counts.
    """
    marker = "__SEEKR_INFO__"
    payload = None
    for line in stdout.splitlines():
        line = line.strip()
        if line.startswith(marker):
            payload = json.loads(line[len(marker):])
    if payload is None or not payload.get("ok"):
        raise RuntimeError(
            f"seekr info fetch failed: {payload!r}")
    launching_info = payload.get("launching_info")
    if not isinstance(launching_info, dict):
        raise RuntimeError(
            f"seekr info fetch missing launching_info: {payload!r}")
    num_swarms = int(launching_info.get("num_swarms", 0))
    assert num_swarms > 0, "No swarms found in seekr info fetch"
    num_anchors = launching_info.get("num_anchors", 0)
    if num_anchors == "null":
        num_anchors = None
    scope = str(launching_info.get("scope"))
    if scope.upper() == "N/A":
        scope = "unpartitioned"
    return StageUnitCounts(
        scope=scope,
        num_swarms=num_swarms,
        num_anchors=num_anchors,
    )

# NOTE: this function is used to obtain numbers of swarms, and maybe anchors,
# after a previous generating stage has finished.
def fetch_unit_counts(
        seekrflow: structures.Seekrflow,
        launching_stage: typing.Any,
        resource: structures.Resource_remote_base,
        silent: bool = False,
        ) -> StageUnitCounts:
    """
    Run seekr status 'info' on a remote worker for launch-time unit enumeration.
    """
    #if resource.type not in ("slurm_remote", "pbs_remote"):
    #    raise NotImplementedError(
    #        f"Resource type {resource.type} is not implemented.")

    launching_stage_index = getattr(
        launching_stage, "index", launching_stage.name)
    model_filename = "model.json"
    lines = [
        "import json, sys",
        "try:",
        "    import seekr.modules.structures as S",
        "    import seekr.status as st",
        f"    model = S.load_model({model_filename!r})",
        f"    launching_name = {launching_stage.name!r}",
        f"    launching_index = {launching_stage_index!r}",
        "    def _entry(msg, stage_name):",
        "        info = msg.get('info') or {}",
        "        for v in info.values():",
        "            if isinstance(v, dict) and v.get('stage_name') == stage_name:",
        "                return v",
        "        raise KeyError(stage_name)",
        "    launch_msg = st.status(model, instruction='info',",
        "        stage_arg=launching_index, print_json=True)",
        "    launching_info = _entry(launch_msg, launching_name)",
        "    payload = {'ok': True, 'launching_info': launching_info}",
        "    print('__SEEKR_INFO__' + json.dumps(payload))",
        "except Exception as e:",
        "    import traceback",
        "    print('__SEEKR_INFO__' + json.dumps({",
        "        'ok': False, 'error': str(e), 'traceback': traceback.format_exc()}))",
    ]
    snippet = "\n".join(lines)
    payload = {
        "root_dir": resolve_model_directory(seekrflow, resource),
        "snippet": snippet,
        "worker_init": getattr(resource, "worker_init", "") or "",
    }
    
    if hasattr(resource, "remote_interface"):
        result = submit_workload(
            resource,
            run_python_snippet_workload,
            payload=payload,
            silent=silent,
        )
    else:
        payload["root_dir"] = str(seekrflow.get_root_directory())
        result = run_python_snippet_workload(payload)
        
    if not result.get("success"):
        raise RuntimeError(
            f"Remote unit-count fetch failed: {result.get('error')}"
            f"stderr={result.get('stderr', '')!r} stdout={result.get('stdout', '')!r}")
    stdout = result.get("stdout", "")
    stderr = result.get("stderr") or ""
    if result.get("returncode") not in (0, None) or "__SEEKR_INFO__" not in stdout:
        raise RuntimeError(
            "Remote unit-count fetch failed: "
            f"returncode={result.get('returncode')} "
            f"stderr={stderr[-2000:]!r} stdout={stdout[-500:]!r}"
        )
    return parse_info_fetch_output(stdout)

def submit_job(
        seekrflow: structures.Seekrflow,
        resource: typing.Any,
        stage_names: list[str],
        job_specs: list[dict],
        resolved_execution: structures.Resolved_execution | None = None,
        time_limit_override: str | None = None,
        silent: bool = False,
        ) -> dict:
    """
    Run the stage on a local, remote, or cloud resource.
    """
    destination_path = resolve_model_directory(seekrflow, resource)
    job_name = seekrflow.name + "_" + stage_names[0]
    kind = resource_kind(resource)
    if resolved_execution is not None:
        if resolved_execution.cpus is not None:
            cpus = resolved_execution.cpus
        elif hasattr(resource, "cpus_per_task"):
            cpus = resource.cpus_per_task
        elif hasattr(resource, "n_vcpus"):
            cpus = resource.n_vcpus
        else:
            cpus = 1
        
    else:
        if hasattr(resource, "cpus_per_task"):
            cpus = resource.cpus_per_task
        elif hasattr(resource, "n_vcpus"):
            cpus = resource.n_vcpus
        else:
            cpus = 1
        

    if resolved_execution is not None and resolved_execution.time_limit is not None:
        time_limit = resolved_execution.time_limit
    elif hasattr(resource, "max_time_limit"):
        time_limit = resource.max_time_limit
    else:
        time_limit = resource.job_timeout_seconds

    if resolved_execution is not None and resolved_execution.memory_mb is not None:
        memory_mb = resolved_execution.memory_mb
    elif hasattr(resource, "memory_per_node"):
        memory_mb = resource.memory_per_node
    elif hasattr(resource, "memory_mb"):
        memory_mb = resource.memory_mb
    else:
        memory_mb = None
    # Pass the contents of job.py and structures.py to the remote worker.
    _JOB_DIR = pathlib.Path(__file__).resolve().parent.parent / "job"
    job_py_source = (_JOB_DIR / "job.py").read_text()
    structures_py_source = (_JOB_DIR / "structures.py").read_text()

    manager_payload = {
        "resource_payload": None,
        "job_specs": [spec.to_dict() for spec in job_specs],
        "runner_files": {
            "job.py": job_py_source,
            "structures.py": structures_py_source,
        }
    }
    effective_time_limit = time_limit_override or time_limit
        
    if resource.type == "direct":
        run_workload = workload_direct.direct_run_workload
        manager_payload["root_dir"] = destination_path
        resource_payload = {
            "job_name": job_name,
            "cpus": cpus,
            "time_limit": effective_time_limit,
            "worker_init": resource.worker_init,
            "n_gpus": resource.n_gpus,
        }
    elif resource.type == "slurm_remote":
        manager_payload["root_dir"] = destination_path
        resource_payload = {
            "scheduler": "slurm",
            "partition_queue": resource.partition,
            "account": resource.account,
            "constraint": resource.constraint,
            "job_name": job_name,
            "cpus": cpus,
            "memory_mb": memory_mb,
            "time_limit": effective_time_limit,
            "scheduler_options": resource.scheduler_options,
            "worker_init": resource.worker_init,
        }
        run_workload = workload_slurm_pbs.slurm_pbs_run_workload
    elif resource.type == "pbs_remote":
        manager_payload["root_dir"] = destination_path
        resource_payload = {
            "scheduler": "pbs",
            "partition_queue": resource.queue,
            "account": resource.account,
            "job_name": job_name,
            "cpus": cpus,
            "memory_mb": memory_mb,
            "time_limit": effective_time_limit,
            "scheduler_options": resource.scheduler_options,
            "worker_init": resource.worker_init,
        }
        run_workload = workload_slurm_pbs.slurm_pbs_run_workload
    elif resource.type == "aws_cloud":
        manager_payload["seekrflow_name"] = seekrflow.name
        resource_payload = {
            "name": resource.name,
            "scheduler": "aws",
            "region": resource.region,
            "transfer_settings_s3_uri": resource.transfer_settings.get_uri(),
            "seekr_image_uri": resource.get_seekr_image_uri(),
            "transfer_settings_input_tarball_name": resource.transfer_settings\
                .get_input_tarball_name(),
            "transfer_settings_bucket": resource.transfer_settings.bucket,
            "job_name": job_name,
            "cpus": cpus,
            "memory_mb": memory_mb,
            "n_gpus": resource.n_gpus,
            "time_limit": (
                effective_time_limit if isinstance(effective_time_limit, int)
                else client_validation.time_limit_to_seconds(effective_time_limit)
            ),
            "job_queue_name": resource.job_queue_name,
            #"scheduler_options": resource.scheduler_options,
        }
        manager_payload["resource_payload"] = resource_payload
        result = workload_aws.submit_aws_job(manager_payload)
        return result
    else:
        raise NotImplementedError(
            f"Resource type {resource.type} is not implemented.")
    
    manager_payload["resource_payload"] = resource_payload
    result = submit_workload(
        resource, run_workload, manager_payload, silent)
    return result

def status(
        resource: structures.Resource_base | None,
        manager_payload: dict,
        silent: bool = False,
        ) -> dict:
    """
    Generalized remote/cloud stage status query for any stage whose producing
    procedure is configured in run_settings.placements.
    """
    
    if resource.type == "direct":
        status_workflow = workload_direct.direct_status_workload
    elif resource.type in ("slurm_remote", "pbs_remote"):
        status_workflow = workload_slurm_pbs.slurm_pbs_status_workload
        if resource.type == "slurm_remote":
            manager_payload["scheduler"] = "slurm"
        elif resource.type == "pbs_remote":
            manager_payload["scheduler"] = "pbs"
    elif resource.type == "aws_cloud":
        manager_payload["scheduler"] = "aws"
        manager_payload["resource_payload"] = {
            "region": resource.region,
            "transfer_settings_s3_uri": resource.transfer_settings.get_uri(),
            "transfer_settings_bucket": resource.transfer_settings.bucket,
        }
        result = workload_aws.status_aws(manager_payload)
        return result
    else:
        raise NotImplementedError(f"Resource type {resource.type} not implemented")
    result = submit_workload(resource, status_workflow, payload=manager_payload, 
                             silent=silent)
    return result

def submit_cancel_workload(
        resource: structures.Resource_base,
        manager_payload: dict,
        silent: bool = False
    ) -> None:
    """
    Cancel a remote or cloud stage job.
    """
    if resource.type == "direct":
        cancel_workload = workload_direct.direct_cancel_workload
        manager_payload["scheduler"] = "direct"
    elif resource.type in ("slurm_remote", "pbs_remote"):
        cancel_workload = workload_slurm_pbs.slurm_pbs_cancel_workload
        if resource.type == "slurm_remote":
            manager_payload["scheduler"] = "slurm"
        elif resource.type == "pbs_remote":
            manager_payload["scheduler"] = "pbs"
    elif resource.type == "aws_cloud":
        resource_payload = {
            "region": resource.region,
            "transfer_settings_s3_uri": resource.transfer_settings.get_uri(),
            "transfer_settings_bucket": resource.transfer_settings.bucket,
            "job_queue_name": resource.job_queue_name,
        }
        manager_payload["resource_payload"] = resource_payload
        result = workload_aws.cancel_aws_job(manager_payload)
        return result
        
    submit_workload(
        resource, cancel_workload, payload=manager_payload, silent=silent)
    return