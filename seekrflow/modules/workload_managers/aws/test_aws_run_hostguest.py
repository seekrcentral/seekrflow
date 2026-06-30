"""
test_aws_run_hostguest.py

This script submits an example host-guest SEEKR calculation
to AWS Batch.

What is needed before running this script:
1. AWS account through Cloudbank (Or just have one yourself).
2. Have created an S3 bucket. Assumed to be named f"seekrflow-data-{ACCOUNT_ID}-{REGION}-an"
3. The docker image named "seekr-engines-bd" with tag "latest" pushed to the ECR
4. A batch compute environment named "seekrflow-gpu-compute"
5. A batch queue named "seekrflow-gpu-queue"
6. Locally: the AWS CLI configured (`aws configure`) and boto3 installed
   (`pip install boto3`).

End-to-end flow:
  upload local model dir -> S3
    -> submit Batch job (job downloads from S3, runs SEEKR, uploads to S3)
      -> wait for SUCCEEDED/FAILED, print CloudWatch logs
        -> download results from S3 for verification.
"""

import os
import sys
import time
import shutil
import tempfile
import subprocess

REGION = "us-west-2"
ACCOUNT_ID = "628031635908"
SEEKR_IMAGE_URI = f"{ACCOUNT_ID}.dkr.ecr.{REGION}.amazonaws.com/seekr-engines-bd:latest"
COMPUTE_ENV_NAME = "seekrflow-gpu-compute"
JOB_QUEUE_NAME = "seekrflow-gpu-queue"
JOB_DEF_GENERIC = "seekrflow-job-generic"

# ---------------------------------------------------------------------------
# Things you may need to fill in / adjust per run.
# ---------------------------------------------------------------------------

# Local directory holding the prepared SEEKR calculation. It MUST contain a
# model.json (plus its anchor_*/ subdirectories). This is the only value you
# are guaranteed to need to change.
# TODO: replace with the real path once the prepared calculation exists.
#LOCAL_MODEL_DIR = "/Users/lvotapka/tmp/seekr3_test_toy_metad"
LOCAL_MODEL_DIR = "/home/lvotapka/tmp/seekr3_test_toy_metad"

# S3 bucket (created in Phase 1) and a prefix (folder) for this test run.
# Inputs are uploaded here; results are synced back to the same prefix
# (overwrite semantics, matching how a real resume would behave).
S3_BUCKET = f"seekrflow-data-{ACCOUNT_ID}-{REGION}-an"
S3_PREFIX = "hostguest-test"
S3_URI = f"s3://{S3_BUCKET}/{S3_PREFIX}"

# Local directory the results are downloaded into for verification.
#LOCAL_RESULTS_DIR = "aws_results"
LOCAL_RESULTS_DIR = LOCAL_MODEL_DIR

# Name of the tarball used to ship the input tree (preserves empty dirs that
# 'aws s3 sync' would otherwise drop).
INPUT_TARBALL_NAME = "model_input.tar.gz"

# SEEKR stage number to run (start with 1).
if len(sys.argv) > 1:
    try:
        SEEKR_STAGE = int(sys.argv[1])
    except ValueError:
        sys.exit(f"ERROR: SEEKR stage must be an integer, got {sys.argv[1]!r}.")
else:
    SEEKR_STAGE = 1

# Container resources. NOTE: g4dn.xlarge / g5.xlarge have 16 GiB RAM, so keep
# MEMORY safely below that (ECS reserves some) or the job will only fit the
# pricey p3.2xlarge -- or get stuck in RUNNABLE.
JOB_GPUS = 1
JOB_VCPUS = 4
JOB_MEMORY_MB = 14000

# Hard wall-clock cap for the Batch job (seconds); the job is killed if exceeded.
JOB_TIMEOUT_SECONDS = 2 * 60 * 60  # 2 hours

# Seconds between local status polls.
POLL_INTERVAL_SECONDS = 20

# Default CloudWatch log group used by the awslogs driver.
LOG_GROUP_NAME = "/aws/batch/job"


def _run(cmd):
    """Run a local command, echoing it; raise CalledProcessError on failure."""
    print(f"$ {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def check_prerequisites():
    """Fail early, with clear messages, before spending any cloud time."""
    if LOCAL_MODEL_DIR == "/path/to/prepared/seekr/calculation" \
            or not os.path.isdir(LOCAL_MODEL_DIR):
        sys.exit(f"ERROR: LOCAL_MODEL_DIR does not exist: {LOCAL_MODEL_DIR!r}\n"
                 "Edit this script and set it to your prepared SEEKR directory "
                 "(the folder containing model.json).")
    if not os.path.isfile(os.path.join(LOCAL_MODEL_DIR, "model.json")):
        sys.exit(f"ERROR: no model.json found in {LOCAL_MODEL_DIR!r}.")
    if shutil.which("aws") is None:
        sys.exit("ERROR: the 'aws' CLI is not on your PATH. Install it and run "
                 "'aws configure' (or 'aws login') first.")


def upload_inputs():
    """Tar the local model directory and upload it as a single S3 object.

    NOTE: 'aws s3 sync' silently drops empty directories because S3 has no
    real directories -- only objects (keys). The prepared SEEKR calculation
    relies on empty dirs (e.g. .../stage_1_equilibration/ where run.log is
    written), so we package the tree as a tarball, which preserves empty
    directories, and extract it inside the container.
    """
    print(f"\n=== Uploading {LOCAL_MODEL_DIR} -> {S3_URI}/{INPUT_TARBALL_NAME} ===")
    tarball_path = os.path.join(tempfile.gettempdir(), INPUT_TARBALL_NAME)
    # -C so the archive contains the directory's contents at the top level
    # (model.json, anchor_*/, ...) rather than the LOCAL_MODEL_DIR name itself.
    _run(["tar", "czf", tarball_path, "-C", LOCAL_MODEL_DIR, "."])
    _run(["aws", "s3", "cp", tarball_path, f"{S3_URI}/{INPUT_TARBALL_NAME}",
          "--region", REGION, "--only-show-errors"])
    # The local tarball has served its purpose once uploaded; don't leave it
    # cluttering the temp directory.
    os.remove(tarball_path)


def register_job_definition(batch):
    """Register a new revision of the generic job definition; return its ARN.

    register_job_definition is idempotent in the sense that it simply creates a
    new revision each call, so we never have to check whether one already exists.
    The command here is a placeholder -- it is overridden at submit time.
    """
    print(f"\n=== Registering job definition {JOB_DEF_GENERIC} ===")
    resp = batch.register_job_definition(
        jobDefinitionName=JOB_DEF_GENERIC,
        type="container",
        platformCapabilities=["EC2"],
        containerProperties={
            "image": SEEKR_IMAGE_URI,
            "resourceRequirements": [
                {"type": "GPU", "value": str(JOB_GPUS)},
                {"type": "VCPU", "value": str(JOB_VCPUS)},
                {"type": "MEMORY", "value": str(JOB_MEMORY_MB)},
            ],
            "command": ["bash", "-lc", "echo overridden-at-submit"],
            "logConfiguration": {"logDriver": "awslogs"},
        },
    )
    arn = resp["jobDefinitionArn"]
    print(f"Registered: {arn}")
    return arn


def build_container_command():
    """Bash run inside the container: pull from S3, run SEEKR, push to S3.

    awscli is pip-installed at runtime so the shared HPC image stays lean (an
    HPC user never needs it). The container reaches S3 using the ecsInstanceRole
    credentials from the instance metadata -- no job role required.

    'set -euo pipefail' guarantees the job exits non-zero (and is marked FAILED)
    if any single step fails, so a permissions or data problem is easy to spot
    in the CloudWatch logs printed by this script.
    """
    script = "\n".join([
        "set -euo pipefail",
        "echo '=== [1/4] Installing awscli ==='",
        "pip install --quiet --no-cache-dir awscli",
        f"echo '=== [2/4] Downloading inputs from {S3_URI}/{INPUT_TARBALL_NAME} ==='",
        "mkdir -p /work",
        f"aws s3 cp {S3_URI}/{INPUT_TARBALL_NAME} /tmp/{INPUT_TARBALL_NAME} "
        f"--region {REGION} --only-show-errors",
        # Extracting the tarball recreates empty directories that 'aws s3 sync'
        # would otherwise have dropped.
        f"tar xzf /tmp/{INPUT_TARBALL_NAME} -C /work",
        f"echo '=== [3/4] Running SEEKR stage {SEEKR_STAGE} ==='",
        "cd /work",
        f"python -m seekr.run /work/model.json {SEEKR_STAGE}",
        f"echo '=== [4/4] Uploading results to {S3_URI} ==='",
        f"aws s3 sync /work {S3_URI} --region {REGION} --only-show-errors",
        "echo '=== DONE ==='",
    ])
    return ["bash", "-lc", script]


def submit_job(batch, job_def_arn):
    """Submit the job with a per-run command override and a wall-clock cap."""
    print(f"\n=== Submitting job to queue {JOB_QUEUE_NAME} ===")
    resp = batch.submit_job(
        jobName=f"seekr-hostguest-stage{SEEKR_STAGE}",
        jobQueue=JOB_QUEUE_NAME,
        jobDefinition=job_def_arn,
        containerOverrides={"command": build_container_command()},
        timeout={"attemptDurationSeconds": JOB_TIMEOUT_SECONDS},
    )
    job_id = resp["jobId"]
    print(f"Submitted jobId={job_id}")
    return job_id


def wait_for_job(batch, job_id):
    """Poll until the job reaches a terminal state; return the final job dict."""
    print("\n=== Waiting for job (Ctrl-C stops watching, job keeps running) ===")
    last_status = None
    while True:
        job = batch.describe_jobs(jobs=[job_id])["jobs"][0]
        status = job["status"]
        if status != last_status:
            reason = job.get("statusReason", "")
            print(f"  status={status}" + (f" ({reason})" if reason else ""))
            last_status = status
        if status in ("SUCCEEDED", "FAILED"):
            return job
        time.sleep(POLL_INTERVAL_SECONDS)


def print_logs(logs, job):
    """Print the container's CloudWatch log stream (the actual SEEKR output)."""
    stream = (job.get("container") or {}).get("logStreamName")
    if not stream:
        for attempt in reversed(job.get("attempts", [])):
            stream = (attempt.get("container") or {}).get("logStreamName")
            if stream:
                break
    if not stream:
        print("No log stream found (the container may never have started). "
              "Check the AWS Batch console for this job.")
        return
    print(f"\n=== CloudWatch logs ({stream}) ===")
    token = None
    while True:
        kwargs = {"logGroupName": LOG_GROUP_NAME, "logStreamName": stream,
                  "startFromHead": True}
        if token:
            kwargs["nextToken"] = token
        resp = logs.get_log_events(**kwargs)
        for event in resp["events"]:
            print(event["message"])
        if not resp["events"] or resp["nextForwardToken"] == token:
            break
        token = resp["nextForwardToken"]


def download_results():
    """Sync the (now updated) S3 prefix back down for local verification."""
    print(f"\n=== Downloading results {S3_URI} -> {LOCAL_RESULTS_DIR} ===")
    os.makedirs(LOCAL_RESULTS_DIR, exist_ok=True)
    _run(["aws", "s3", "sync", S3_URI, LOCAL_RESULTS_DIR,
          "--region", REGION, "--only-show-errors"])
    # The input tarball lives in the same S3 prefix, so the sync pulls it back
    # down too. It is not a result -- drop it from the local results directory.
    downloaded_tarball = os.path.join(LOCAL_RESULTS_DIR, INPUT_TARBALL_NAME)
    if os.path.isfile(downloaded_tarball):
        os.remove(downloaded_tarball)


def main():
    check_prerequisites()
    try:
        import boto3
    except ImportError:
        sys.exit("ERROR: boto3 is required locally. Run 'pip install boto3'.")

    batch = boto3.client("batch", region_name=REGION)
    logs = boto3.client("logs", region_name=REGION)

    upload_inputs()
    job_def_arn = register_job_definition(batch)
    job_id = submit_job(batch, job_def_arn)
    job = wait_for_job(batch, job_id)
    print_logs(logs, job)

    if job["status"] == "SUCCEEDED":
        download_results()
        print(f"\nSUCCESS: SEEKR stage {SEEKR_STAGE} ran end-to-end on AWS Batch. "
              f"Results downloaded to {LOCAL_RESULTS_DIR!r}.")
    else:
        sys.exit(f"\nFAILED: job {job_id} ended with status {job['status']}. "
                 "See the CloudWatch logs above for the cause.")


if __name__ == "__main__":
    main()
