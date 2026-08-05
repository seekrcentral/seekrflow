"""
AWS Batch workload manager (client-side boto3).

Submit / status / cancel for ``aws_cloud`` resources. Seekr progress is
published by the container to S3; manager liveness comes from
``batch.describe_jobs``.
"""
from __future__ import annotations

import json
import os
import re
import typing

from seekrflow.modules.transfer.aws_s3 import LOCAL_RUNNER_DIRNAMES

if typing.TYPE_CHECKING:
    from seekrflow.modules import structures

AWS_RUNNER_DIRNAME = ".aws_runner"
STATUS_POLL_SECONDS = 30
CONTAINER_WORK_DIR = "/work"
# Default CloudWatch log group for Batch jobs using the awslogs driver.
CLOUDWATCH_LOG_GROUP = "/aws/batch/job"
FAILURE_LOG_TAIL_LINES = 100


def sanitize_job_definition_name(name: str) -> str:
    """AWS job definition names: letters, numbers, underscores, hyphens."""
    cleaned = re.sub(r"[^A-Za-z0-9_-]+", "-", name.strip()) or "seekrflow"
    return cleaned[:128]


def job_definition_name_for_resource(resource: "structures.Resource_cloud_aws") -> str:
    return sanitize_job_definition_name(f"seekrflow-{resource.name}")


def s3_model_uri(resource: "structures.Resource_cloud_aws", seekrflow_name: str) -> str:
    """S3 URI for this run's model tree (tarball + synced outputs)."""
    base = resource.transfer_settings.get_uri().rstrip("/")
    return f"{base}/{seekrflow_name}"


def s3_status_key(
        resource: "structures.Resource_cloud_aws",
        seekrflow_name: str,
        stage_name: str,
        ) -> str:
    prefix = resource.transfer_settings.prefix.strip("/")
    parts = [p for p in (prefix, seekrflow_name, AWS_RUNNER_DIRNAME) if p]
    return "/".join(parts + [f"{stage_name}_status.json"])


def s3_status_uri(
        resource: "structures.Resource_cloud_aws",
        seekrflow_name: str,
        stage_name: str,
        ) -> str:
    return (
        f"s3://{resource.transfer_settings.bucket}/"
        f"{s3_status_key(resource, seekrflow_name, stage_name)}"
    )


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
        resource: "structures.Resource_cloud_aws",
        seekrflow_name: str,
        host_stage_name: str,
        ) -> str:
    """S3 object key for the Batch run script (kept out of containerOverrides)."""
    prefix = resource.transfer_settings.prefix.strip("/")
    parts = [p for p in (prefix, seekrflow_name, AWS_RUNNER_DIRNAME) if p]
    return "/".join(parts + [f"{host_stage_name}_run.sh"])


def s3_run_script_uri(
        resource: "structures.Resource_cloud_aws",
        seekrflow_name: str,
        host_stage_name: str,
        ) -> str:
    return (
        f"s3://{resource.transfer_settings.bucket}/"
        f"{s3_run_script_key(resource, seekrflow_name, host_stage_name)}"
    )


def local_runner_dir(model_directory: str) -> str:
    path = os.path.join(model_directory, AWS_RUNNER_DIRNAME)
    os.makedirs(path, exist_ok=True)
    return path


def local_latest_state_path(model_directory: str, stage_name: str) -> str:
    return os.path.join(
        local_runner_dir(model_directory), f"{stage_name}_latest.json")


def local_failure_log_path(model_directory: str, stage_name: str) -> str:
    return os.path.join(
        local_runner_dir(model_directory), f"{stage_name}_failure.log")


def local_failure_reported_path(model_directory: str, stage_name: str) -> str:
    return os.path.join(
        local_runner_dir(model_directory), f"{stage_name}_failure_reported.json")


def clear_stage_failure_artifacts(model_directory: str, stage_name: str) -> None:
    """Remove local failure log / reported marker before a fresh submit."""
    for path in (
            local_failure_log_path(model_directory, stage_name),
            local_failure_reported_path(model_directory, stage_name),
    ):
        if os.path.isfile(path):
            try:
                os.remove(path)
            except OSError:
                pass


def _failure_already_reported(
        model_directory: str,
        stage_name: str,
        job_id: str,
        ) -> bool:
    path = local_failure_reported_path(model_directory, stage_name)
    if not os.path.isfile(path):
        return False
    try:
        with open(path, "r") as f:
            data = json.load(f)
        return str(data.get("job_id") or "") == str(job_id)
    except Exception:
        return False


def _mark_failure_reported(
        model_directory: str,
        stage_name: str,
        job_id: str,
        ) -> None:
    path = local_failure_reported_path(model_directory, stage_name)
    try:
        with open(path, "w") as f:
            json.dump({"job_id": str(job_id)}, f)
    except OSError:
        pass


def save_job_state(
        model_directory: str,
        stage_name: str,
        job_id: str,
        job_name: str,
        *,
        array_size: int | None = None,
        ) -> None:
    path = local_latest_state_path(model_directory, stage_name)
    payload: dict = {
        "job_id": job_id,
        "job_name": job_name,
        "stage_name": stage_name,
    }
    if array_size is not None and array_size > 1:
        payload["array_size"] = int(array_size)
    # Preserve walltime baselines written before submit created this file.
    pending = os.path.join(
        model_directory, ".seekrflow_walltime", f"{stage_name}.json")
    for source_path in (path, pending):
        if not os.path.isfile(source_path):
            continue
        try:
            with open(source_path, "r") as f:
                prior = json.load(f)
            if isinstance(prior, dict):
                for key in (
                        "work_at_submission",
                        "last_known_elapsed",
                        "last_known_elapsed_seconds",
                ):
                    if key in prior and key not in payload:
                        payload[key] = prior[key]
        except Exception:
            pass
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)


def load_job_state(
        model_directory: str,
        stage_name: str,
        ) -> dict | None:
    path = local_latest_state_path(model_directory, stage_name)
    if not os.path.isfile(path):
        return None
    try:
        with open(path, "r") as f:
            return json.load(f)
    except Exception:
        return None


def clear_job_state(model_directory: str, stage_name: str) -> None:
    path = local_latest_state_path(model_directory, stage_name)
    if os.path.isfile(path):
        try:
            os.remove(path)
        except OSError:
            pass


def _batch_client(resource: "structures.Resource_cloud_aws"):
    import boto3
    return boto3.client("batch", region_name=resource.region)


def _s3_client(resource: "structures.Resource_cloud_aws"):
    import boto3
    return boto3.client("s3", region_name=resource.region)


def _logs_client(resource: "structures.Resource_cloud_aws"):
    import boto3
    return boto3.client("logs", region_name=resource.region)


def _map_batch_status(batch_status: str) -> str:
    """Map AWS Batch status to seekrflow manager job State strings."""
    status = (batch_status or "").upper()
    if status in ("RUNNING", "STARTING"):
        return "RUNNING"
    if status in ("SUBMITTED", "PENDING", "RUNNABLE"):
        return "PENDING"
    if status == "SUCCEEDED":
        return "COMPLETED"
    if status in ("FAILED",):
        return "FAILED"
    return status or "UNKNOWN"


def _container_from_batch_job(job: dict) -> dict:
    """Prefer the latest attempt container that has a log stream or reason."""
    container = dict(job.get("container") or {})
    if container.get("logStreamName") or container.get("reason"):
        return container
    for attempt in reversed(job.get("attempts") or []):
        attempt_container = dict(attempt.get("container") or {})
        if (attempt_container.get("logStreamName")
                or attempt_container.get("reason")
                or attempt_container.get("exitCode") is not None):
            return attempt_container
    return container


def format_batch_failure_summary(job: dict) -> str:
    """Short human-readable Batch failure summary (no CloudWatch body)."""
    job_id = job.get("jobId") or "?"
    status_reason = (job.get("statusReason") or "").strip()
    container = _container_from_batch_job(job)
    exit_code = container.get("exitCode")
    container_reason = (container.get("reason") or "").strip()
    parts = [f"Batch job FAILED (jobId={job_id})"]
    if status_reason:
        parts.append(f"statusReason={status_reason}")
    if exit_code is not None:
        parts.append(f"exitCode={exit_code}")
    if container_reason:
        parts.append(f"containerReason={container_reason}")
    return "; ".join(parts)


def fetch_cloudwatch_log_tail(
        resource: "structures.Resource_cloud_aws",
        log_stream_name: str,
        *,
        limit: int = FAILURE_LOG_TAIL_LINES,
        log_group: str = CLOUDWATCH_LOG_GROUP,
        ) -> str:
    """
    Return the last ``limit`` CloudWatch log lines for a Batch container stream.

    Empty string if the stream is missing or logs cannot be read (permissions,
    container never started, etc.).
    """
    if not log_stream_name:
        return ""
    try:
        logs = _logs_client(resource)
        # startFromHead=False walks from the end of the stream.
        resp = logs.get_log_events(
            logGroupName=log_group,
            logStreamName=log_stream_name,
            limit=max(1, int(limit)),
            startFromHead=False,
        )
        messages = [
            (event.get("message") or "").rstrip()
            for event in (resp.get("events") or [])
        ]
        return "\n".join(m for m in messages if m)
    except Exception as e:
        return f"(CloudWatch log fetch failed: {e})"


def collect_batch_failure_diagnostics(
        resource: "structures.Resource_cloud_aws",
        job: dict,
        *,
        model_directory: str,
        stage_name: str,
        announce: bool = True,
        ) -> tuple[str, str]:
    """
    Build failure summary + optional CloudWatch tail; write ``*_failure.log``.

    Returns ``(short_error, notes_with_log_excerpt)``.
    When ``announce`` is True, print diagnostics at most once per Batch job id.
    """
    summary = format_batch_failure_summary(job)
    container = _container_from_batch_job(job)
    stream = container.get("logStreamName") or ""
    log_tail = ""
    if stream:
        log_tail = fetch_cloudwatch_log_tail(resource, stream)

    failure_path = local_failure_log_path(model_directory, stage_name)
    body_lines = [
        summary,
        f"logStreamName={stream or '(none)'}",
        f"logGroup={CLOUDWATCH_LOG_GROUP}",
        "",
        "=== CloudWatch log tail ===",
        log_tail or (
            "(No CloudWatch log stream — container may never have started.)"
            if not stream else "(empty)"
        ),
        "",
    ]
    try:
        with open(failure_path, "w") as f:
            f.write("\n".join(body_lines))
    except OSError:
        failure_path = "(could not write local failure log)"

    # Only surface a CloudWatch excerpt when we actually have stream content.
    has_useful_logs = bool(stream and log_tail and not log_tail.startswith(
        "(CloudWatch log fetch failed"))
    excerpt = (
        "\n".join(log_tail.splitlines()[-40:]) if has_useful_logs else "")

    notes_parts = [summary, f"Failure details written to {failure_path}"]
    if not stream:
        notes_parts.append(
            "No CloudWatch log stream (container may never have started).")
    elif excerpt:
        notes_parts.append("CloudWatch log tail:\n" + excerpt)
    notes = "\n".join(notes_parts)

    # Print at most once per Batch job id. Startup probes often re-discover a
    # previously FAILED job (stale *_latest.json) before a fresh submit; those
    # calls pass announce=False. Monitor calls announce=True.
    job_id = str(job.get("jobId") or "")
    if (
            announce
            and job_id
            and not _failure_already_reported(
                model_directory, stage_name, job_id)
    ):
        print(f"[aws-batch] stage {stage_name}: {summary}")
        print(f"[aws-batch] failure log: {failure_path}")
        if excerpt:
            print(f"[aws-batch] CloudWatch log tail:\n{excerpt}")
        _mark_failure_reported(model_directory, stage_name, job_id)
    return summary, notes


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


def build_status_publisher_snippet(stage_names: list[str]) -> str:
    """
    Python run inside the container to write per-stage status + dispatch JSON.

    Progress handling mirrors local/SLURM: read seekr's ``progress`` field
    (stage-level float, or mean of per-anchor/swarm ``progress`` values).

    Writes:
      ``/tmp/<stage>_status.json`` — monitor payload (success/stage_status)
      ``/tmp/<stage>_dispatch.json`` — seekr ``info`` entry + raw progress
        for client-side unit enumeration without a full S3 tree pull
    """
    stages_repr = repr(list(stage_names))
    return f"""
import json, os
import seekr.modules.structures as structures
import seekr.status as seekr_status

def _progress_from_stage_progress(stage_progress):
    # Same contract as local_multiprocessing / slurm: seekr always exposes
    # a ``progress`` value at stage level and per partition (anchor/swarm).
    progress_map = stage_progress.get("progress", {{}})
    if isinstance(progress_map, (int, float)):
        return float(progress_map)
    values = []
    if isinstance(progress_map, dict):
        for _, per_partition in progress_map.items():
            if isinstance(per_partition, dict):
                value = per_partition.get("progress")
                if isinstance(value, (int, float)):
                    values.append(float(value))
            elif isinstance(per_partition, (int, float)):
                values.append(float(per_partition))
    if values:
        return sum(values) / len(values)
    return 0.0

def _match_stage_entry(root, stage_name):
    if not isinstance(root, dict):
        return {{}}
    for _k, _v in root.items():
        if isinstance(_v, dict) and (
                _v.get("stage_name") == stage_name
                or _v.get("name") == stage_name
                or _k == stage_name):
            return _v
    return {{}}

model = structures.load_model("model.json")
_started_path = "/tmp/seekrflow_stages_started"
_started = set()
if os.path.isfile(_started_path):
    with open(_started_path) as _sf:
        _started = {{line.strip() for line in _sf if line.strip()}}

for stage_name in {stages_repr}:
    info_entry = None
    stage_progress = {{}}
    # Fused jobs sync prior outputs into /work before later stages run.
    # Do not publish seekr progress for a stage until its command has started,
    # or stale on-disk completion from a previous run looks "done" early.
    if stage_name not in _started:
        payload = {{
            "success": True,
            "error": None,
            "stage_status": {{
                "finished": False,
                "state": "unstarted",
                "progress": 0.0,
                "notes": "waiting for prior stage(s) in this Batch job",
            }},
        }}
        with open(f"/tmp/{{stage_name}}_status.json", "w") as f:
            json.dump(payload, f)
        with open(f"/tmp/{{stage_name}}_dispatch.json", "w") as f:
            json.dump({{"version": 1, "info": None, "progress": None}}, f)
        continue

    # ``info`` is optional for monitor health: older seekr images lack it.
    # Keep progress publishing working even when info is unsupported.
    try:
        info_msg = seekr_status.status(
            model, instruction="info", stage_arg=stage_name, print_json=True)
        if isinstance(info_msg, str):
            info_msg = json.loads(info_msg)
        info_entry = _match_stage_entry(info_msg.get("info", {{}}), stage_name)
        if not info_entry and isinstance(info_msg.get("info"), dict):
            _vals = [v for v in info_msg["info"].values() if isinstance(v, dict)]
            if len(_vals) == 1:
                info_entry = _vals[0]
    except Exception as _info_err:
        info_entry = None
        _info_note = f"info unavailable: {{_info_err}}"
    else:
        _info_note = ""

    try:
        msg = seekr_status.status(
            model, instruction="progress", stage_arg=stage_name, print_json=True)
        if isinstance(msg, str):
            msg = json.loads(msg)
        stage_progress = _match_stage_entry(msg.get("progress", {{}}), stage_name)
        finished = bool(stage_progress.get("finished", False))
        if finished:
            progress = 1.0
            state = "completed"
        else:
            progress = _progress_from_stage_progress(stage_progress)
            state = "started" if progress > 0.0 else "unstarted"
        payload = {{
            "success": True,
            "error": None,
            "stage_status": {{
                "finished": finished,
                "state": state,
                "progress": progress,
                "notes": _info_note,
            }},
        }}
    except Exception as e:
        payload = {{
            "success": False,
            "error": str(e),
            "stage_status": {{
                "finished": False,
                "state": "unknown",
                "progress": 0.0,
                "notes": str(e),
            }},
        }}
    with open(f"/tmp/{{stage_name}}_status.json", "w") as f:
        json.dump(payload, f)
    with open(f"/tmp/{{stage_name}}_dispatch.json", "w") as f:
        json.dump({{
            "version": 1,
            "info": info_entry,
            "progress": stage_progress if stage_progress else None,
        }}, f)
""".strip()


def build_aws_container_script(
        resource: "structures.Resource_cloud_aws",
        seekrflow_name: str,
        stage_specs: list[dict],
        ) -> str:
    """
    Full container run script (download inputs, status loop, seekr, upload).

    Uploaded to S3 at submit time; the Batch ``containerOverrides.command`` only
    fetches and executes this file (AWS caps overrides at 8192 bytes).

    Each ``stage_specs`` entry needs ``name``, ``force_overwrite``,
    ``benchmark_mode``, and optionally ``seekr_command`` (pre-built bash
    fragment from dispatch lowering). When ``seekr_command`` is absent, a
    whole-stage ``seekr_run`` with ``'any'`` is used.
    """
    if not stage_specs:
        raise ValueError("stage_specs must be non-empty")
    s3_uri = s3_model_uri(resource, seekrflow_name)
    tarball = resource.transfer_settings.get_input_tarball_name()
    region = resource.region
    stage_names = [s["name"] for s in stage_specs]
    status_snippet = build_status_publisher_snippet(stage_names)
    seekr_cmds = []
    for spec in stage_specs:
        custom = spec.get("seekr_command")
        if custom:
            seekr_cmds.append(str(custom).strip())
        else:
            seekr_cmds.append(
                build_seekr_stage_python_command(
                    spec["name"],
                    force_overwrite=bool(spec.get("force_overwrite", False)),
                    benchmark_mode=bool(spec.get("benchmark_mode", False)),
                )
            )
    # Mark each stage started immediately before its command so the status
    # publisher ignores stale on-disk completion for later fused members.
    wrapped_cmds = []
    for spec, cmd in zip(stage_specs, seekr_cmds):
        name = spec["name"]
        wrapped_cmds.append(
            f'echo {name!r} >> /tmp/seekrflow_stages_started && ( {cmd} )'
        )
    seekr_chain = " && ".join(wrapped_cmds)
    status_uploads = "\n".join(
        (
            f"    aws s3 cp /tmp/{name}_status.json "
            f"{s3_status_uri(resource, seekrflow_name, name)} "
            f"--region {region} --only-show-errors || true\n"
            f"    aws s3 cp /tmp/{name}_dispatch.json "
            f"{s3_dispatch_uri(resource, seekrflow_name, name)} "
            f"--region {region} --only-show-errors || true"
        )
        for name in stage_names
    )
    # Keep client/local runner state out of the shared model tree on S3.
    runner_excludes = " ".join(
        f"--exclude '{name}/*'" for name in LOCAL_RUNNER_DIRNAMES
    )

    return f"""
set -euo pipefail
echo '=== [1/6] Installing awscli ==='
pip install --quiet --no-cache-dir awscli
echo '=== [2/6] Checking seekr status instructions in this image ==='
python - <<'SEEKRCHECK_EOF'
import seekr, seekr.status as st
print("seekr module:", getattr(seekr, "__file__", seekr))
print("AVAILABLE_INSTRUCTIONS:", st.AVAILABLE_INSTRUCTIONS)
if "info" not in st.AVAILABLE_INSTRUCTIONS:
    print(
        "WARNING: this image's seekr lacks the 'info' instruction. "
        "Progress monitoring still works; cloud unit enumeration "
        "(dispatch artifacts) will be incomplete until the image is rebuilt "
        "with a newer seekr."
    )
else:
    print("OK: 'info' instruction is available.")
SEEKRCHECK_EOF
echo '=== [3/6] Seeding /work from tarball, then syncing S3 tree ==='
mkdir -p {CONTAINER_WORK_DIR}
# Tarball preserves empty dirs that S3 sync cannot represent. Extract first,
# then overlay the live cloud prefix so later Batch jobs see prior stage outputs
# without a client-side re-tar / re-transfer.
aws s3 cp {s3_uri}/{tarball} /tmp/{tarball} --region {region} --only-show-errors
# GNU tar exits 1 on recoverable warnings (e.g. dangling symlinks).
tar xzf /tmp/{tarball} -C {CONTAINER_WORK_DIR} || {{
  _tar_rc=$?
  if [ "$_tar_rc" -gt 1 ]; then exit "$_tar_rc"; fi
}}
aws s3 sync {s3_uri} {CONTAINER_WORK_DIR} --region {region} --only-show-errors \
  --exclude {tarball} {runner_excludes}
rm -f {CONTAINER_WORK_DIR}/{tarball}
cd {CONTAINER_WORK_DIR}
cat > /tmp/publish_seekr_status.py << 'SEEKRSTATUS_EOF'
{status_snippet}
SEEKRSTATUS_EOF
echo '=== [4/6] Starting status publisher ==='
: > /tmp/seekrflow_stages_started
(
  while true; do
    python /tmp/publish_seekr_status.py || true
{status_uploads}
    sleep {STATUS_POLL_SECONDS}
  done
) &
STATUS_PID=$!
cleanup() {{ kill "$STATUS_PID" 2>/dev/null || true; }}
trap cleanup EXIT
echo '=== [5/6] Running SEEKR stage(s): {" ".join(stage_names)} ==='
echo "AWS_BATCH_JOB_ARRAY_INDEX=${{AWS_BATCH_JOB_ARRAY_INDEX:-}}"
{seekr_chain}
echo '=== [6/6] Final status publish + upload results to {s3_uri} ==='
python /tmp/publish_seekr_status.py || true
{status_uploads}
aws s3 sync {CONTAINER_WORK_DIR} {s3_uri} --region {region} --only-show-errors \
  --exclude {tarball} {runner_excludes}
echo '=== DONE ==='
""".strip()


def upload_aws_run_script(
        resource: "structures.Resource_cloud_aws",
        seekrflow_name: str,
        host_stage_name: str,
        script: str,
        ) -> str:
    """Upload the run script to S3; return its ``s3://`` URI."""
    s3 = _s3_client(resource)
    key = s3_run_script_key(resource, seekrflow_name, host_stage_name)
    s3.put_object(
        Bucket=resource.transfer_settings.bucket,
        Key=key,
        Body=script.encode("utf-8"),
        ContentType="text/x-shellscript",
    )
    return s3_run_script_uri(resource, seekrflow_name, host_stage_name)


def build_aws_container_command(
        resource: "structures.Resource_cloud_aws",
        seekrflow_name: str,
        host_stage_name: str,
        ) -> list[str]:
    """
    Short Batch ``containerOverrides.command``: fetch run script from S3 and exec.

    Must stay well under AWS's 8192-byte containerOverrides limit.
    """
    script_uri = s3_run_script_uri(resource, seekrflow_name, host_stage_name)
    region = resource.region
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


def register_job_definition(
        resource: "structures.Resource_cloud_aws",
        *,
        cpus: int | None = None,
        memory_mb: int | None = None,
        ) -> str:
    """Register a new job-definition revision; return its ARN."""
    batch = _batch_client(resource)
    name = job_definition_name_for_resource(resource)
    n_vcpus = cpus if cpus is not None else resource.n_vcpus
    mem = memory_mb if memory_mb is not None else resource.memory_mb
    requirements = [
        {"type": "VCPU", "value": str(n_vcpus)},
        {"type": "MEMORY", "value": str(mem)},
    ]
    if resource.n_gpus > 0:
        requirements.insert(0, {"type": "GPU", "value": str(resource.n_gpus)})
    resp = batch.register_job_definition(
        jobDefinitionName=name,
        type="container",
        platformCapabilities=["EC2"],
        containerProperties={
            "image": resource.get_seekr_image_uri(),
            "resourceRequirements": requirements,
            "command": ["bash", "-lc", "echo overridden-at-submit"],
            "logConfiguration": {"logDriver": "awslogs"},
        },
    )
    return resp["jobDefinitionArn"]


def submit_aws_job(
        seekrflow: "structures.Seekrflow",
        host_stage_name: str,
        resource: "structures.Resource_cloud_aws",
        stage_specs: list[dict],
        model_directory: str,
        job_name: str | None = None,
        timeout_seconds: int | None = None,
        cpus: int | None = None,
        memory_mb: int | None = None,
        array_size: int | None = None,
        ) -> dict:
    """
    Register job def, submit Batch job for the host (possibly fused) stage set.

    ``cpus`` / ``memory_mb`` / ``timeout_seconds`` override Resource defaults
    when provided (typically from ``Resolved_execution``).
    ``array_size`` > 1 submits a Batch array job (one child per member).

    Returns ``{success, job_id, job_name, error}``.
    """
    try:
        # Drop stale failure diagnostics from a prior Batch job for this stage
        # so startup probes / leftover FAILED ids do not keep alarming.
        clear_stage_failure_artifacts(model_directory, host_stage_name)
        job_def_arn = register_job_definition(
            resource, cpus=cpus, memory_mb=memory_mb)
        name = job_name or f"{seekrflow.name}_{host_stage_name}"
        script = build_aws_container_script(
            resource, seekrflow.name, stage_specs)
        script_uri = upload_aws_run_script(
            resource, seekrflow.name, host_stage_name, script)
        print(f"[aws-batch] uploaded run script to {script_uri}")
        command = build_aws_container_command(
            resource, seekrflow.name, host_stage_name)
        timeout = (
            timeout_seconds
            if timeout_seconds is not None
            else resource.job_timeout_seconds)
        n_array = int(array_size) if array_size is not None else 1
        if n_array < 1:
            n_array = 1
        submit_kwargs: dict = {
            "jobName": name,
            "jobQueue": resource.job_queue_name,
            "jobDefinition": job_def_arn,
            "containerOverrides": {"command": command},
            "timeout": {"attemptDurationSeconds": int(timeout)},
        }
        if n_array > 1:
            submit_kwargs["arrayProperties"] = {"size": n_array}
            print(
                f"[aws-batch] submitting array job size={n_array} "
                f"for stage {host_stage_name}"
            )
        batch = _batch_client(resource)
        resp = batch.submit_job(**submit_kwargs)
        job_id = resp["jobId"]
        save_job_state(
            model_directory,
            host_stage_name,
            job_id,
            name,
            array_size=n_array if n_array > 1 else None,
        )
        return {
            "success": True,
            "job_id": job_id,
            "job_name": name,
            "array_size": n_array if n_array > 1 else None,
            "error": None,
        }
    except Exception as e:
        return {
            "success": False,
            "job_id": None,
            "job_name": None,
            "array_size": None,
            "error": str(e),
        }


def _download_stage_status_from_s3(
        resource: "structures.Resource_cloud_aws",
        seekrflow_name: str,
        stage_name: str,
        ) -> dict | None:
    try:
        s3 = _s3_client(resource)
        key = s3_status_key(resource, seekrflow_name, stage_name)
        resp = s3.get_object(
            Bucket=resource.transfer_settings.bucket, Key=key)
        body = resp["Body"].read().decode("utf-8")
        return json.loads(body)
    except Exception:
        return None


def download_stage_dispatch_from_s3(
        resource: "structures.Resource_cloud_aws",
        seekrflow_name: str,
        stage_name: str,
        ) -> dict | None:
    """Return ``{version, info, progress}`` for ``stage_name``, or None."""
    try:
        s3 = _s3_client(resource)
        key = s3_dispatch_key(resource, seekrflow_name, stage_name)
        resp = s3.get_object(
            Bucket=resource.transfer_settings.bucket, Key=key)
        body = resp["Body"].read().decode("utf-8")
        payload = json.loads(body)
        return payload if isinstance(payload, dict) else None
    except Exception:
        return None


def _stage_name_for_index(model: typing.Any, stage_index: int) -> str:
    for stage in getattr(model, "stages", []) or []:
        if getattr(stage, "index", None) == stage_index:
            return stage.name
    stages = list(getattr(model, "stages", []) or [])
    if 0 <= stage_index < len(stages):
        return stages[stage_index].name
    raise ValueError(
        f"No stage with index {stage_index!r} in model "
        f"({[getattr(s, 'name', s) for s in stages]!r}).")


def fetch_unit_counts_cloud(
        resource: "structures.Resource_cloud_aws",
        seekrflow_name: str,
        model: typing.Any,
        launching_stage: typing.Any,
        ) -> "dispatch_lowering.StageUnitCounts":
    """
    Unit counts for cloud dispatch from the input stage's S3 dispatch artifact.

    When the input is the prepared initial stage (index 0), fall back to local
    seekr ``info`` — that tree is created on the client before the first upload.
    """
    from seekrflow.modules.workload_managers import dispatch_lowering

    input_index = getattr(launching_stage, "input_stage_index", 0)
    if input_index == 0:
        return dispatch_lowering.fetch_unit_counts_local(model, launching_stage)

    input_name = _stage_name_for_index(model, input_index)
    key = s3_dispatch_key(resource, seekrflow_name, input_name)
    uri = s3_dispatch_uri(resource, seekrflow_name, input_name)
    payload = download_stage_dispatch_from_s3(
        resource, seekrflow_name, input_name)
    info = (payload or {}).get("info") if payload else None
    if not isinstance(info, dict) or not info:
        status_payload = _download_stage_status_from_s3(
            resource, seekrflow_name, input_name)
        status_state = None
        if status_payload and isinstance(status_payload.get("stage_status"), dict):
            status_state = status_payload["stage_status"].get("state")
        detail = (
            f"object missing at {uri}"
            if payload is None
            else f"object exists at {uri} but info is empty/null"
        )
        hint = ""
        if status_state == "completed":
            hint = (
                f" Stage {input_name!r} has a completed status JSON but no "
                f"usable dispatch artifact — it likely finished under an older "
                f"seekrflow run script. Re-run that stage (e.g. -f {input_name} "
                f"or its fusion host) with current seekrflow so the container "
                f"publishes {input_name}_dispatch.json."
            )
        else:
            hint = (
                f" Re-run input stage {input_name!r} on AWS with current "
                f"seekrflow so the container publishes {input_name}_dispatch.json."
            )
        raise RuntimeError(
            f"Missing S3 dispatch info for input stage {input_name!r} "
            f"({detail}; key={key}).{hint}")
    return dispatch_lowering.resolve_unit_counts(
        info, None, input_is_initial=False)


def fetch_stage_progress_cloud(
        resource: "structures.Resource_cloud_aws",
        seekrflow_name: str,
        stage_name: str,
        ) -> tuple[bool, dict]:
    """
    ``(finished, partition_progress_map)`` from the stage's S3 dispatch artifact.

    Missing artifact ⇒ not finished, empty map (treat all units incomplete).
    """
    payload = download_stage_dispatch_from_s3(
        resource, seekrflow_name, stage_name)
    stage_progress = (payload or {}).get("progress")
    if not isinstance(stage_progress, dict):
        return False, {}
    finished = bool(stage_progress.get("finished", False))
    partition_map = stage_progress.get("progress", {})
    if not isinstance(partition_map, dict):
        partition_map = {}
    return finished, partition_map


def publish_stage_completed_status(
        resource: "structures.Resource_cloud_aws",
        seekrflow_name: str,
        stage_name: str,
        *,
        notes: str = "all dispatch units already attained; skipped Batch submit",
        ) -> None:
    """
    Write a completed stage status JSON to S3 (client-side).

    Used when incomplete-unit filtering leaves nothing to run, so probes and
    the monitor treat the stage as done without a Batch job.
    """
    payload = {
        "success": True,
        "error": None,
        "stage_status": {
            "finished": True,
            "state": "completed",
            "progress": 1.0,
            "notes": notes,
        },
    }
    s3 = _s3_client(resource)
    key = s3_status_key(resource, seekrflow_name, stage_name)
    s3.put_object(
        Bucket=resource.transfer_settings.bucket,
        Key=key,
        Body=json.dumps(payload).encode("utf-8"),
        ContentType="application/json",
    )


def _find_job_id(
        resource: "structures.Resource_cloud_aws",
        model_directory: str,
        stage_name: str,
        seekrflow_name: str,
        ) -> str | None:
    state = load_job_state(model_directory, stage_name)
    if state and state.get("job_id"):
        return str(state["job_id"])
    # Fall back: list recent jobs by name
    job_name = f"{seekrflow_name}_{stage_name}"
    try:
        batch = _batch_client(resource)
        for status in ("RUNNING", "RUNNABLE", "STARTING", "SUBMITTED", "PENDING"):
            resp = batch.list_jobs(
                jobQueue=resource.job_queue_name,
                jobStatus=status,
                filters=[{"name": "JOB_NAME", "values": [job_name]}],
            )
            jobs = resp.get("jobSummaryList") or []
            if jobs:
                return jobs[0]["jobId"]
    except Exception:
        pass
    return None


def _array_child_job_ids(parent_job_id: str, array_size: int) -> list[str]:
    return [f"{parent_job_id}:{i}" for i in range(int(array_size))]


def _collect_active_batch_jobs(
        batch,
        parent_job: dict,
        *,
        array_size: int | None,
        ) -> tuple[str | None, list[dict], dict | None]:
    """
    Return ``(batch_status, manager_jobs, failure_job)``.

    For array jobs, ``manager_jobs`` lists active children; ``failure_job`` is
    the first FAILED child (else the parent) for diagnostics.
    """
    parent_id = parent_job.get("jobId")
    parent_status = parent_job.get("status")
    mapped_parent = _map_batch_status(parent_status or "")
    manager_jobs: list[dict] = []
    failure_job: dict | None = None

    n_array = int(array_size) if array_size and array_size > 1 else 0
    if n_array > 1 and parent_id:
        child_ids = _array_child_job_ids(str(parent_id), n_array)
        # describe_jobs accepts up to 100 ids per call
        children: list[dict] = []
        for i in range(0, len(child_ids), 100):
            chunk = child_ids[i:i + 100]
            children.extend(
                batch.describe_jobs(jobs=chunk).get("jobs") or [])
        for child in children:
            mapped = _map_batch_status(child.get("status") or "")
            if mapped in ("RUNNING", "PENDING"):
                manager_jobs.append({
                    "JobID": child.get("jobId"),
                    "State": mapped,
                    "JobName": child.get("jobName"),
                    "StatusReason": child.get("statusReason", ""),
                })
            if child.get("status") == "FAILED" and failure_job is None:
                failure_job = child
        # Parent array job itself is often not a useful "running" row; if any
        # child is active, treat overall status as RUNNING for coarse checks.
        if manager_jobs:
            batch_status = "RUNNING"
        elif any(c.get("status") == "FAILED" for c in children):
            batch_status = "FAILED"
        elif children and all(c.get("status") == "SUCCEEDED" for c in children):
            batch_status = "SUCCEEDED"
        else:
            batch_status = parent_status
        if failure_job is None and parent_status == "FAILED":
            failure_job = parent_job
        return batch_status, manager_jobs, failure_job

    if mapped_parent in ("RUNNING", "PENDING"):
        manager_jobs.append({
            "JobID": parent_id,
            "State": mapped_parent,
            "JobName": parent_job.get("jobName"),
            "StatusReason": parent_job.get("statusReason", ""),
        })
    if parent_status == "FAILED":
        failure_job = parent_job
    return parent_status, manager_jobs, failure_job


def status_aws(
        seekrflow: "structures.Seekrflow",
        stage_name: str,
        resource: "structures.Resource_cloud_aws",
        model_directory: str,
        *,
        host_stage_name: str | None = None,
        benchmark_mode: bool = False,
        announce_failure: bool = True,
        ) -> dict:
    """
    Combine Batch ``describe_jobs`` (manager) with S3 status JSON (stage).

    For fused members, pass ``host_stage_name`` so manager liveness is taken
    from the host Batch job while Seekr progress is still this stage's S3 file.

    On Batch ``FAILED``, attaches statusReason / container exit info and a
    CloudWatch log tail (also written under ``.aws_runner/{stage}_failure.log``).
    Set ``announce_failure=False`` for quiet pre-submit probes of stale jobs.
    """
    del benchmark_mode  # reserved; container publishes completion itself
    manager_stage = host_stage_name or stage_name
    manager_status = {
        "tool": "aws_batch",
        "timestamp": None,
        "error": "",
        "jobs": [],
    }
    state = load_job_state(model_directory, manager_stage)
    array_size = None
    if state:
        try:
            array_size = int(state["array_size"]) if state.get("array_size") else None
        except (TypeError, ValueError):
            array_size = None
    job_id = _find_job_id(
        resource, model_directory, manager_stage, seekrflow.name)
    batch_status = None
    batch_job: dict | None = None
    failure_job: dict | None = None
    if job_id:
        try:
            batch = _batch_client(resource)
            jobs = batch.describe_jobs(jobs=[job_id]).get("jobs") or []
            if jobs:
                batch_job = jobs[0]
                batch_status, active_jobs, failure_job = (
                    _collect_active_batch_jobs(
                        batch, batch_job, array_size=array_size))
                manager_status["jobs"].extend(active_jobs)
                save_job_state(
                    model_directory,
                    manager_stage,
                    job_id,
                    batch_job.get("jobName")
                    or f"{seekrflow.name}_{manager_stage}",
                    array_size=array_size,
                )
                started_at = batch_job.get("startedAt")
                stopped_at = batch_job.get("stoppedAt")
                if (isinstance(started_at, (int, float))
                        and isinstance(stopped_at, (int, float))
                        and stopped_at >= started_at):
                    elapsed_s = float(stopped_at - started_at) / 1000.0
                    state_path = local_latest_state_path(
                        model_directory, manager_stage)
                    try:
                        with open(state_path, "r") as handle:
                            state_data = json.load(handle)
                        if not isinstance(state_data, dict):
                            state_data = {}
                        state_data["last_known_elapsed_seconds"] = elapsed_s
                        with open(state_path, "w") as handle:
                            json.dump(state_data, indent=2, fp=handle)
                    except Exception:
                        pass
        except Exception as e:
            manager_status["error"] = str(e)

    s3_payload = _download_stage_status_from_s3(
        resource, seekrflow.name, stage_name)
    if s3_payload and isinstance(s3_payload.get("stage_status"), dict):
        stage_status = dict(s3_payload["stage_status"])
        success = bool(s3_payload.get("success", True))
        error = s3_payload.get("error")
    else:
        # No status object yet: infer coarse state from Batch
        success = True
        error = None
        if batch_status in ("RUNNING", "STARTING"):
            stage_status = {
                "finished": False,
                "state": "started",
                "progress": 0.0,
                "notes": "waiting for first S3 status publish",
            }
        elif batch_status == "SUCCEEDED":
            stage_status = {
                "finished": True,
                "state": "completed",
                "progress": 1.0,
                "notes": "Batch SUCCEEDED before S3 status appeared",
            }
        elif batch_status == "FAILED":
            stage_status = {
                "finished": False,
                "state": "failed",
                "progress": 0.0,
                "notes": "Batch job FAILED",
            }
            success = False
            error = "Batch job FAILED"
        elif manager_status["jobs"]:
            stage_status = {
                "finished": False,
                "state": "started",
                "progress": 0.0,
                "notes": "",
            }
        else:
            stage_status = {
                "finished": False,
                "state": "unstarted",
                "progress": 0.0,
                "notes": "",
            }

    # Batch FAILED wins over a stale/incomplete S3 progress snapshot: early
    # container crashes often never publish status JSON.
    diag_job = failure_job or (
        batch_job if batch_status == "FAILED" else None)
    if diag_job is not None and (
            batch_status == "FAILED" or failure_job is not None):
        summary, notes = collect_batch_failure_diagnostics(
            resource,
            diag_job,
            model_directory=model_directory,
            stage_name=manager_stage,
            announce=announce_failure,
        )
        stage_status = {
            "finished": False,
            "state": "failed",
            "progress": float(stage_status.get("progress") or 0.0),
            "notes": notes,
        }
        success = False
        error = summary

    return {
        "success": success,
        "error": error,
        "manager_status": manager_status,
        "stage_status": stage_status,
    }


def cancel_aws_job(
        resource: "structures.Resource_cloud_aws",
        *,
        job_id: str | None = None,
        job_name: str | None = None,
        reason: str = "seekrflow cancel",
        ) -> dict:
    """Terminate Batch job(s) by id and/or name."""
    batch = _batch_client(resource)
    canceled: list[str] = []
    errors: list[str] = []
    ids: list[str] = []
    if job_id:
        ids.append(str(job_id))
    if job_name and not job_id:
        try:
            for status in (
                    "RUNNING", "RUNNABLE", "STARTING", "SUBMITTED", "PENDING"):
                resp = batch.list_jobs(
                    jobQueue=resource.job_queue_name,
                    jobStatus=status,
                    filters=[{"name": "JOB_NAME", "values": [job_name]}],
                )
                for summary in resp.get("jobSummaryList") or []:
                    jid = summary.get("jobId")
                    if jid and jid not in ids:
                        ids.append(jid)
        except Exception as e:
            errors.append(str(e))
    for jid in ids:
        try:
            batch.terminate_job(jobId=jid, reason=reason)
            canceled.append(jid)
        except Exception as e:
            errors.append(f"{jid}: {e}")
    return {
        "success": len(errors) == 0,
        "canceled": canceled,
        "error": "; ".join(errors) if errors else None,
    }


def clear_stage_s3_monitor_artifacts(
        resource: "structures.Resource_cloud_aws",
        seekrflow_name: str,
        stage_name: str,
        ) -> list[str]:
    """
    Delete S3 ``*_status.json`` / ``*_dispatch.json`` for ``stage_name``.

    Model trajectory data is left alone (seekr ``force_overwrite`` handles that).
    Stale monitor JSON with ``state=completed`` would otherwise make force-rerun
    look finished before the new Batch job publishes fresh status/dispatch.
    """
    s3 = _s3_client(resource)
    bucket = resource.transfer_settings.bucket
    deleted: list[str] = []
    for key in (
            s3_status_key(resource, seekrflow_name, stage_name),
            s3_dispatch_key(resource, seekrflow_name, stage_name),
            ):
        try:
            s3.delete_object(Bucket=bucket, Key=key)
            deleted.append(key)
        except Exception:
            pass
    return deleted


def cancel_and_reset_aws_stage(
        seekrflow: "structures.Seekrflow",
        stage_name: str,
        resource: "structures.Resource_cloud_aws",
        model_directory: str,
        *,
        monitor_stage_names: list[str] | None = None,
        ) -> dict:
    """
    Terminate any Batch job for the stage, clear local runner state, and
    remove S3 monitor artifacts for this stage (and optional fused members).
    """
    state = load_job_state(model_directory, stage_name)
    job_id = state.get("job_id") if state else None
    job_name = (
        (state.get("job_name") if state else None)
        or f"{seekrflow.name}_{stage_name}"
    )
    result = cancel_aws_job(
        resource, job_id=job_id, job_name=job_name,
        reason="seekrflow force overwrite / reset")
    clear_job_state(model_directory, stage_name)
    stages_to_clear = list(monitor_stage_names) if monitor_stage_names else [stage_name]
    # Always include the host stage name.
    if stage_name not in stages_to_clear:
        stages_to_clear.insert(0, stage_name)
    cleared_keys: list[str] = []
    for name in stages_to_clear:
        cleared_keys.extend(
            clear_stage_s3_monitor_artifacts(resource, seekrflow.name, name))
    result["canceled_job"] = (
        result["canceled"][0] if result.get("canceled") else None)
    result["cleared_s3_keys"] = cleared_keys
    return result


# Back-compat stub name used by early remote.py wiring
def aws_cloud_status(*args, **kwargs):
    raise NotImplementedError(
        "Use status_aws(seekrflow, stage_name, resource, model_directory)"
    )
