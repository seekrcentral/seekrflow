"""
AWS Batch workload manager (client-side boto3).

Submit / status / cancel for "aws_cloud" resources. Seekr progress is
published by the container to S3; manager liveness comes from
"batch.describe_jobs".
"""
from __future__ import annotations

import os
import re
import json
import time
import typing
from urllib.parse import urlparse

from seekrflow.modules.transfer.aws_s3 import LOCAL_RUNNER_DIRNAMES

if typing.TYPE_CHECKING:
    from seekrflow.modules import structures
    import boto3

AWS_RUNNER_DIRNAME = ".aws_runner"
STATUS_POLL_SECONDS = 30
CONTAINER_ROOT_DIR = "/work" # Needs to be 'work' so as not to conflict with Linux /root
# Default CloudWatch log group for Batch jobs using the awslogs driver.
CLOUDWATCH_LOG_GROUP = "/aws/batch/job"
FAILURE_LOG_TAIL_LINES = 100

def sanitize_job_definition_name(name: str) -> str:
    """
    AWS job definition names: letters, numbers, underscores, hyphens.
    """
    cleaned = re.sub(r"[^A-Za-z0-9_-]+", "-", name.strip()) or "seekrflow"
    return cleaned[:128]

def job_definition_name_for_resource(resource_payload: dict) -> str:
    """
    Return the AWS job definition name for a resource.
    """
    return sanitize_job_definition_name(f"seekrflow-{resource_payload['name']}")

def s3_status_key(
        resource: "structures.Resource_cloud_aws",
        seekrflow_name: str,
        stage_name: str,
        ) -> str:
    prefix = resource.transfer_settings.prefix.strip("/")
    parts = [p for p in (prefix, seekrflow_name, AWS_RUNNER_DIRNAME) if p]
    return "/".join(parts + [f"{stage_name}_status.json"])

def s3_dispatch_key(
        resource: "structures.Resource_cloud_aws",
        seekrflow_name: str,
        stage_name: str,
        ) -> str:
    """Seekr ``info`` + ``progress`` snapshot for cloud dispatch sizing."""
    prefix = resource.transfer_settings.prefix.strip("/")
    parts = [p for p in (prefix, seekrflow_name, AWS_RUNNER_DIRNAME) if p]
    return "/".join(parts + [f"{stage_name}_dispatch.json"])


def s3_dispatch_uri(
        resource: "structures.Resource_cloud_aws",
        seekrflow_name: str,
        stage_name: str,
        ) -> str:
    return (
        f"s3://{resource.transfer_settings.bucket}/"
        f"{s3_dispatch_key(resource, seekrflow_name, stage_name)}"
    )


def s3_run_script_key(
        resource_payload: dict,
        seekrflow_name: str,
        job_name: str,
        ) -> str:
    """
    S3 object key for the Batch run script (kept out of containerOverrides).
    """
    prefix = urlparse(resource_payload['transfer_settings_s3_uri']).path.strip("/")
    parts = [p for p in (prefix, seekrflow_name, AWS_RUNNER_DIRNAME) if p]
    return "/".join(parts + [f"{job_name}_run.sh"])


def s3_run_script_uri(
        resource_payload: dict,
        seekrflow_name: str,
        job_name: str,
        ) -> str:
    return (
        f"s3://{resource_payload['transfer_settings_bucket']}/"
        f"{s3_run_script_key(resource_payload, seekrflow_name, job_name)}"
    )

def local_runner_dir(model_directory: str) -> str:
    path = os.path.join(model_directory, AWS_RUNNER_DIRNAME)
    os.makedirs(path, exist_ok=True)
    return path

def local_failure_log_path(model_directory: str, stage_name: str) -> str:
    return os.path.join(
        local_runner_dir(model_directory), f"{stage_name}_failure.log")


def local_failure_reported_path(model_directory: str, stage_name: str) -> str:
    return os.path.join(
        local_runner_dir(model_directory), f"{stage_name}_failure_reported.json")


def clear_stage_failure_artifacts(
        model_directory: str, 
        stage_name: str) -> None:
    """
    Remove local failure log / reported marker before a fresh submit.
    """
    for path in (
            local_failure_log_path(model_directory, stage_name),
            local_failure_reported_path(model_directory, stage_name),
    ):
        if os.path.isfile(path):
            try:
                os.remove(path)
            except OSError:
                pass

def _batch_client(resource_payload: dict):
    import boto3
    return boto3.client("batch", region_name=resource_payload["region"])


def _s3_client(resource_payload: dict):
    import boto3
    return boto3.client("s3", region_name=resource_payload["region"])

def map_batch_state(status: str, reason: str = "") -> str:
        code = (status or "").upper()
        if code in ("RUNNING", "STARTING"):
            return "running"
        if code in ("SUBMITTED", "PENDING", "RUNNABLE"):
            return "queued"
        if code == "FAILED":
            if "cancel" in (reason or "").lower():
                return "cancelled"
            return "failed"
        # SUCCEEDED, or the id is already gone: same as a finished job leaving squeue.
        return "idle"

def build_seekr_stage_python_command(
        stage_name: str,
        *,
        force_overwrite: bool = False,
        benchmark_mode: bool = False,
        ) -> str:
    """Seekr invocation inside the container (cwd=/work)."""
    force_flag = "True" if force_overwrite else "False"
    bench_flag = "True" if benchmark_mode else "False"
    return (
        "python -c \"import seekr.modules.structures as structures; "
        "import seekr.run as seekr_run; "
        "model = structures.load_model('model.json'); "
        f"seekr_run.run(model, '{stage_name}', 'any', None, "
        f"{force_flag}, None, benchmark={bench_flag})\""
        f" > {stage_name}_run.out"
    )

def build_aws_container_script(manager_payload: dict) -> str:
    """
    Container script: sync the model tree, run job.py, sync results back.
    """
    resource_payload = manager_payload["resource_payload"]
    s3_uri = (
        f"{resource_payload['transfer_settings_s3_uri'].rstrip('/')}/"
        f"{manager_payload['seekrflow_name']}"
    )
    tarball = resource_payload["transfer_settings_input_tarball_name"]
    region = resource_payload["region"]
    cpus = resource_payload.get("cpus") or 1
    # NOTE: I don't think there needs to be a worker_init for AWS.
    #worker_init = resource_payload.get("worker_init") or ""
    spec0 = manager_payload["job_specs"][0]
    internal_id = spec0["internal_id"]
    interval = spec0["status_write_interval"]
    root = CONTAINER_ROOT_DIR
    jobs = f"{root}/.seekr_jobs"
    # Left for the shell. Batch sets this to 0..N-1.
    index = "${AWS_BATCH_JOB_ARRAY_INDEX}" \
        if len(manager_payload["job_specs"]) > 1 else "0"
    spec = f"{jobs}/job_spec_{internal_id}_{index}.json"
    wrap = f"python {jobs}/job.py {spec}"
    #if worker_init:
    #    wrap = f"{worker_init}\n{wrap}"
    runner_excludes = " ".join(
        f"--exclude '{name}/*'" for name in LOCAL_RUNNER_DIRNAMES
    )
    return f"""
set -euo pipefail
pip install --quiet --no-cache-dir awscli
mkdir -p {root}
aws s3 cp {s3_uri}/{tarball} /tmp/{tarball} --region {region} --only-show-errors
tar xzf /tmp/{tarball} -C {root} || {{
  _tar_rc=$?
  if [ "$_tar_rc" -gt 1 ]; then exit "$_tar_rc"; fi
}}
aws s3 sync {s3_uri} {root} --region {region} --only-show-errors \\
  --exclude {tarball} {runner_excludes}
rm -f {root}/{tarball}
cd {root}
export OMP_NUM_THREADS={cpus}
export OPENMM_CPU_THREADS={cpus}
(
  while true; do
    sleep {interval}
    aws s3 sync {root}/.stage_states {s3_uri}/.stage_states --region {region} --only-show-errors || true
  done
) &
sync_pid=$!
set +e
{wrap}
rc=$?
set -e
kill "$sync_pid" 2>/dev/null || true
wait "$sync_pid" 2>/dev/null || true
aws s3 sync {root} {s3_uri} --region {region} --only-show-errors \
  --exclude {tarball} {runner_excludes}
exit $rc
""".strip()


def upload_aws_run_script(
        resource_payload: dict,
        seekrflow_name: str,
        job_name: str,
        script: str,
        ) -> str:
    """Upload the run script to S3; return its ``s3://`` URI."""
    s3 = _s3_client(resource_payload)
    key = s3_run_script_key(resource_payload, seekrflow_name, job_name)
    s3.put_object(
        Bucket=resource_payload["transfer_settings_bucket"],
        Key=key,
        Body=script.encode("utf-8"),
        ContentType="text/x-shellscript",
    )
    return s3_run_script_uri(resource_payload, seekrflow_name, job_name)


def build_aws_container_command(
        resource_payload: dict,
        seekrflow_name: str,
        job_name: str,
        ) -> list[str]:
    """
    Short Batch ``containerOverrides.command``: fetch run script from S3 and exec.

    Must stay well under AWS's 8192-byte containerOverrides limit.
    """
    script_uri = s3_run_script_uri(resource_payload, seekrflow_name, job_name)
    region = resource_payload["region"]
    stub = f"""
set -euo pipefail
pip install --quiet --no-cache-dir awscli
aws s3 cp {script_uri} /tmp/seekrflow_run.sh --region {region} --only-show-errors
exec bash /tmp/seekrflow_run.sh
""".strip()
    encoded_len = len(json.dumps(["bash", "-lc", stub]))
    if encoded_len > 8192:
        raise ValueError(
            f"AWS containerOverrides stub is {encoded_len} bytes "
            f"(limit 8192); shorten build_aws_container_command.")
    return ["bash", "-lc", stub]

def _matching_job_definition_arn(batch, name, image, requirements):
    """Active revision with this image and resource requirements, if any."""
    wanted = {(item["type"], str(item["value"])) for item in requirements}
    token = None
    while True:
        kwargs = {"jobDefinitionName": name, "status": "ACTIVE"}
        if token:
            kwargs["nextToken"] = token
        resp = batch.describe_job_definitions(**kwargs)
        for defn in resp.get("jobDefinitions") or []:
            props = defn.get("containerProperties") or {}
            got = {
                (item["type"], str(item["value"]))
                for item in props.get("resourceRequirements") or []
            }
            if props.get("image") == image and got == wanted:
                return defn["jobDefinitionArn"]
        token = resp.get("nextToken")
        if not token:
            return None

def register_job_definition(resource_payload: dict) -> str:
    """
    Register a new job-definition revision; return its ARN.
    """
    batch = _batch_client(resource_payload)
    name = job_definition_name_for_resource(resource_payload)
    n_vcpus = resource_payload["cpus"]
    mem = resource_payload["memory_mb"]
    requirements = [
        {"type": "VCPU", "value": str(n_vcpus)},
        {"type": "MEMORY", "value": str(mem)},
    ]
    if resource_payload["n_gpus"] > 0:
        requirements.insert(0, {"type": "GPU", "value": str(resource_payload["n_gpus"])})
    image = resource_payload["seekr_image_uri"]
    arn = _matching_job_definition_arn(batch, name, image, requirements)
    if arn is not None:
        return arn
    try:
        resp = batch.register_job_definition(
            jobDefinitionName=name,
            type="container",
            platformCapabilities=["EC2"],
            containerProperties={
                "image": resource_payload["seekr_image_uri"],
                "resourceRequirements": requirements,
                "command": ["bash", "-lc", "echo overridden-at-submit"],
                "logConfiguration": {"logDriver": "awslogs"},
            },
        )
    except Exception as exc:
        if "versions do not match" not in str(exc):
            raise
        arn = _matching_job_definition_arn(batch, name, image, requirements)
        if arn is None:
            raise
        return arn
    return resp["jobDefinitionArn"]

def next_internal_id(
        s3: typing.Any, 
        bucket: str, 
        s3_prefix: str
        ) -> int:
    """
    Next seekrflow id from ".batch_runner/batch_info_{id}.json" on S3.
    """
    prefix = f"{s3_prefix.strip('/')}/.batch_runner/batch_info_".lstrip("/")
    ids: list[int] = []
    for page in s3.get_paginator("list_objects_v2").paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents") or []:
            stem = obj["Key"].rsplit("/", 1)[-1]
            if not (stem.startswith("batch_info_") and stem.endswith(".json")):
                continue
            token = stem[len("batch_info_"):-len(".json")]
            if token.isdigit():
                ids.append(int(token))
    return max(ids) + 1 if ids else 0

def submit_aws_job(
        manager_payload: dict,
        ) -> dict:
    """
    Register job def, submit Batch job for the host (possibly fused) stage set.

    Returns "{success, job_id, job_name, error}".
    """
    resource_payload = manager_payload["resource_payload"]
    seekrflow_name = manager_payload["seekrflow_name"]
    job_name = resource_payload["job_name"]
    #job_specs = manager_payload["job_specs"]
    # Convert the dictionary of job specs to JobSpec objects.
    namespace = {}
    exec(manager_payload["runner_files"]["structures.py"], namespace)
    JobSpec = namespace["JobSpec"]
    job_specs = [
        JobSpec.from_dict(spec) if isinstance(spec, dict) else spec
        for spec in manager_payload["job_specs"]
    ]
    
    try:
        # Drop stale failure diagnostics from a prior Batch job for this stage
        # so startup probes / leftover FAILED ids do not keep alarming.

        # TODO: restore?
        #clear_stage_failure_artifacts(model_directory, host_stage_name)

        job_def_arn = register_job_definition(resource_payload)
        
        # Create the state files for the job.
        s3 = _s3_client(resource_payload)
        bucket = resource_payload["transfer_settings_bucket"]
        s3_prefix = (
            f"{urlparse(resource_payload['transfer_settings_s3_uri']).path.strip('/')}/"
            f"{seekrflow_name}"
        ).strip("/")
        internal_id = next_internal_id(s3, bucket, s3_prefix)
        # Set the job spec to have the correct internal id - needed for job.py run
        for job_spec in job_specs:
            job_spec.internal_id = internal_id
            job_spec.remote_root_dir = CONTAINER_ROOT_DIR

        info_key = f"{s3_prefix}/.batch_runner/batch_info_{internal_id}.json".lstrip("/")
        s3.put_object(
            Bucket=bucket,
            Key=info_key,
            Body=json.dumps({"internal_id": internal_id}).encode("utf-8"),
            ContentType="application/json",
        )

        jobs_prefix = f"{s3_prefix}/.seekr_jobs".strip("/")
        runner_files = manager_payload["runner_files"]
        for file_name in ("job.py", "structures.py"):
            s3.put_object(
                Bucket=bucket,
                Key=f"{jobs_prefix}/{file_name}",
                Body=runner_files[file_name].encode("utf-8"),
                ContentType="text/x-python",
            )
        for job_spec in job_specs:
            s3.put_object(
                Bucket=bucket,
                Key=f"{jobs_prefix}/job_spec_{internal_id}_{job_spec.array_index}.json",
                Body=json.dumps(job_spec.to_dict()).encode("utf-8"),
                ContentType="application/json",
            )
        
        manager_payload["job_specs"] = [spec.to_dict() for spec in job_specs]
        script = build_aws_container_script(manager_payload)
        script_uri = upload_aws_run_script(
            resource_payload, seekrflow_name, job_name, script)
        print(f"[aws-batch] uploaded run script to {script_uri}")
        command = build_aws_container_command(
            resource_payload, seekrflow_name, job_name)
        timeout = resource_payload["time_limit"]
        n_array = len(job_specs)
        if n_array < 1:
            raise ValueError("AWS job array size must be at least 1")
        submit_kwargs: dict = {
            "jobName": job_name,
            "jobQueue": resource_payload["job_queue_name"],
            "jobDefinition": job_def_arn,
            "containerOverrides": {"command": command},
            "timeout": {"attemptDurationSeconds": int(timeout)},
        }
        if n_array > 1:
            submit_kwargs["arrayProperties"] = {"size": n_array}
            print(
                f"[aws-batch] submitting array job size={n_array} "
                f"for job {job_name}"
            )
        batch = _batch_client(resource_payload)
        resp = batch.submit_job(**submit_kwargs)
        job_id = resp["jobId"]
        return {
            "success": True,
            "error": None,
            "internal_id": internal_id,
            "job_id": job_id,
            "job_name": job_name,
        }
    except Exception as e:
        return {
            "success": False,
            "error": str(e),
        }

def format_elapsed(
        job: dict) -> str | None:
    started = job.get("startedAt")
    if not isinstance(started, (int, float)) or started <= 0:
        return None
    stopped = job.get("stoppedAt")
    end = stopped if isinstance(stopped, (int, float)) and stopped >= started \
        else time.time() * 1000
    seconds = int((end - started) / 1000)
    hours, rem = divmod(seconds, 3600)
    minutes, secs = divmod(rem, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"

def describe(
        batch: "boto3.client",
        job_ids: list[str]) -> dict[str, dict]:
    found = {}
    for start in range(0, len(job_ids), 100):
        chunk = job_ids[start:start + 100]
        for job in batch.describe_jobs(jobs=chunk).get("jobs") or []:
            found[job["jobId"]] = job
    return found

def load_stage_state(
        s3: "boto3.client",
        bucket: str,
        s3_prefix: str,
        internal_id: int, 
        stage_index: int,
        array_index: int
        ) -> dict:
    key = (
        f"{s3_prefix}/.stage_states/"
        f"stage_state_{internal_id}_{stage_index}_{array_index}.json"
    ).lstrip("/")
    try:
        body = s3.get_object(Bucket=bucket, Key=key)["Body"].read()
        return json.loads(body)
    except Exception as error:
        code = getattr(getattr(error, "response", {}), "get", lambda *_: {})("Error", {})
        if isinstance(code, dict) and code.get("Code") in ("404", "NoSuchKey"):
            return {
                "internal_id": internal_id,
                "state": "unstarted",
                "finished": False,
                "stage_anchor_swarm_progress_list": [],
            }
        raise

def _failure_notes(resource_payload: dict, child: dict | None, state: str) -> str | None:
    """Batch statusReason plus the CloudWatch tail when the job failed."""
    if child is None:
        return None
    container = child.get("container") or {}
    reason = child.get("statusReason") or ""
    exit_code = container.get("exitCode")
    notes = f"{reason} (exit {exit_code})".strip() if exit_code is not None else reason
    log_stream = container.get("logStreamName")
    if state != "failed" or not log_stream:
        return notes or None
    import boto3
    logs = boto3.client("logs", region_name=resource_payload["region"])
    events = logs.get_log_events(
        logGroupName=CLOUDWATCH_LOG_GROUP,
        logStreamName=log_stream,
        limit=FAILURE_LOG_TAIL_LINES,
        startFromHead=False,
    ).get("events") or []
    tail = "\n".join(event["message"] for event in events).strip()
    if tail:
        notes = f"{notes}\n{tail}" if notes else tail
    return notes or None

def status_aws(manager_payload: dict) -> dict:
    """
    Batch describe_jobs for manager state, S3 stage-state JSON for seekr state.
    Returns the standard result shape.
    """
    resource_payload = manager_payload["resource_payload"]
    bucket = resource_payload["transfer_settings_bucket"]
    batch = _batch_client(resource_payload)
    s3 = _s3_client(resource_payload)
    return_payload = {}
    try:
        for system_name, payload in manager_payload["system_payloads"].items():
            job_dicts_by_job_id = {}
            s3_prefix = (
                f"{urlparse(resource_payload['transfer_settings_s3_uri']).path.strip('/')}/"
                f"{system_name}"
            ).strip("/")
            for job in payload["jobs"]:
                internal_id = job["internal_id"]
                job_id = job.get("job_id")
                array_indices = list(job.get("array_indices") or [])
                stage_indices = list(job.get("stage_indices") or [])
                if not job_id:
                    continue
                parent_id = str(job_id)
                ids = [parent_id]
                if len(array_indices) > 1:
                    ids.extend(f"{parent_id}:{index}" for index in array_indices)
                described = describe(batch, ids)
                parent = described.get(parent_id, {})
                parent_state = map_batch_state(
                    parent.get("status", ""), parent.get("statusReason", ""))
                batch_dicts_by_array_index = {}
                for array_index in array_indices:
                    if len(array_indices) > 1:
                        child = described.get(f"{parent_id}:{array_index}")
                    else:
                        child = parent or None
                    if child is None:
                        state = "queued" if parent_state in ("queued", "running") else "idle"
                        elapsed = None
                        known = []
                    else:
                        state = map_batch_state(
                            child.get("status", ""), child.get("statusReason", ""))
                        elapsed = format_elapsed(child)
                        known = [parent_id] if state in ("running", "queued") else []
                    batch_dicts_by_array_index[array_index] = {
                        "internal_id": internal_id,
                        "filename": f"batch_state_{internal_id}_{array_index}.json",
                        "batch_info_filename": f"batch_info_{internal_id}.json",
                        "state": state,
                        "notes": _failure_notes(resource_payload, child, state),
                        "last_timestamp": time.time(),
                        "last_known_elapsed": elapsed,
                        "last_known_jobs": known,
                    }
                    
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
                    member_dicts = [
                        load_stage_state(
                            s3, bucket, s3_prefix, internal_id, stage_index, array_index)
                        for array_index in array_indices
                    ]
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

                job_dicts_by_job_id[parent_id] = {
                    "manager_dicts_by_array_index": batch_dicts_by_array_index,
                    "stage_dicts_by_stage_index": stage_dicts_by_stage_index,
                }
            return_payload[system_name] = job_dicts_by_job_id
    except Exception as error:
        return {"success": False, "error": str(error), "payload": None}
    return {"success": True, "error": None, "payload": return_payload}

def cancel_aws_job(manager_payload: dict) -> dict:
    """
    Terminate Batch jobs by id and/or name.

    Payload matches status_aws, plus "remove_json_files".
    Returns "{success, error, payload}" where payload is the canceled ids.
    """
    resource_payload = manager_payload["resource_payload"]
    remove_json_files = bool(manager_payload.get("remove_json_files", False))
    batch = _batch_client(resource_payload)
    s3 = _s3_client(resource_payload)
    bucket = resource_payload["transfer_settings_bucket"]
    queue = resource_payload["job_queue_name"]
    canceled: list[str] = []
    errors: list[str] = []

    def terminate(job_id: str) -> None:
        if job_id in canceled:
            return
        try:
            batch.terminate_job(jobId=job_id, reason="seekrflow cancel")
            canceled.append(job_id)
        except Exception as error:
            errors.append(f"{job_id}: {error}")

    def ids_for_name(job_name: str) -> list[str]:
        found: list[str] = []
        for status in ("RUNNING", "RUNNABLE", "STARTING", "SUBMITTED", "PENDING"):
            resp = batch.list_jobs(
                jobQueue=queue,
                jobStatus=status,
                filters=[{"name": "JOB_NAME", "values": [job_name]}],
            )
            for summary in resp.get("jobSummaryList") or []:
                job_id = summary.get("jobId")
                if job_id and job_id not in found:
                    found.append(job_id)
        return found

    def delete_prefix(prefix: str) -> None:
        paginator = s3.get_paginator("list_objects_v2")
        for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
            keys = [{"Key": obj["Key"]} for obj in page.get("Contents") or []]
            for start in range(0, len(keys), 1000):
                s3.delete_objects(
                    Bucket=bucket, Delete={"Objects": keys[start:start + 1000]})

    try:
        for name, payload in manager_payload["system_payloads"].items():
            for job in payload["jobs"]:
                job_id = job.get("job_id")
                job_name = job.get("job_name")
                if not job_id and not job_name:
                    return {
                        "success": False,
                        "error": "job_id or job_name required",
                        "payload": None,
                    }
                if job_id:
                    terminate(str(job_id))
                elif job_name:
                    for found_id in ids_for_name(str(job_name)):
                        terminate(found_id)
            
            if remove_json_files:
                s3_prefix = (
                    f"{urlparse(resource_payload['transfer_settings_s3_uri']).path.strip('/')}/"
                    f"{name}"
                ).strip("/")
                for dirname in (".batch_runner", ".stage_states"):
                    delete_prefix(f"{s3_prefix}/{dirname}/".lstrip("/"))
    except Exception as error:
        return {"success": False, "error": str(error), "payload": None}

    return {
        "success": not errors,
        "error": "; ".join(errors) if errors else None,
        "payload": canceled,
    }
