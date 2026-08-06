#!/usr/bin/env python3
"""
Check whether the AWS SEEKR container image provides seekr.status 'info'.

Two modes:

  docker  (default, fastest) — run the check in a local/ECR image via docker
  batch   — submit a one-shot AWS Batch job using the same image URI that
            seekrflow registers, then print CloudWatch logs

Examples (from the repo root):

  # Fastest — pull/run the ECR image locally with docker:
  python seekrflow/modules/workload_managers/aws/check_seekr_image_instructions.py \\
      --image 628031635908.dkr.ecr.us-west-2.amazonaws.com/seekr-engines-bd:latest

  # Or force a one-shot Batch job (same image URI seekrflow uses):
  python seekrflow/modules/workload_managers/aws/check_seekr_image_instructions.py \\
      --mode batch \\
      --region us-west-2 \\
      --account-id 628031635908 \\
      --job-queue seekrflow-gpu-queue

Exit codes: 0 if 'info' is present, 1 otherwise (or on tool failure).
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time

CHECK_PY = r"""
import seekr, seekr.status as st
print("seekr module:", getattr(seekr, "__file__", seekr))
print("AVAILABLE_INSTRUCTIONS:", st.AVAILABLE_INSTRUCTIONS)
ok = "info" in st.AVAILABLE_INSTRUCTIONS
print("HAS_INFO:", ok)
raise SystemExit(0 if ok else 1)
"""


def _run(cmd: list[str], *, check: bool = True) -> subprocess.CompletedProcess:
    print("$", " ".join(cmd), flush=True)
    return subprocess.run(cmd, check=check)


def check_via_docker(image: str) -> int:
    # Ensure we can talk to ECR if this is an ECR URI
    if ".dkr.ecr." in image and ".amazonaws.com/" in image:
        # 123.dkr.ecr.us-west-2.amazonaws.com/repo:tag
        host = image.split("/")[0]
        region = host.split(".")[3]
        login = subprocess.run(
            ["aws", "ecr", "get-login-password", "--region", region],
            check=True,
            capture_output=True,
            text=True,
        )
        subprocess.run(
            [
                "docker", "login", "--username", "AWS",
                "--password-stdin", host,
            ],
            input=login.stdout,
            text=True,
            check=True,
        )
        _run(["docker", "pull", image])

    result = subprocess.run(
        ["docker", "run", "--rm", "--entrypoint", "python", image, "-c", CHECK_PY],
        check=False,
    )
    return int(result.returncode)


def check_via_batch(
        *,
        image: str,
        region: str,
        job_queue: str,
        job_def_name: str,
        n_gpus: int,
        n_vcpus: int,
        memory_mb: int,
        timeout_seconds: int,
        ) -> int:
    import boto3

    batch = boto3.client("batch", region_name=region)
    logs = boto3.client("logs", region_name=region)

    requirements = [
        {"type": "VCPU", "value": str(n_vcpus)},
        {"type": "MEMORY", "value": str(memory_mb)},
    ]
    if n_gpus > 0:
        requirements.insert(0, {"type": "GPU", "value": str(n_gpus)})

    print(f"Registering job definition {job_def_name!r} with image {image}", flush=True)
    reg = batch.register_job_definition(
        jobDefinitionName=job_def_name,
        type="container",
        platformCapabilities=["EC2"],
        containerProperties={
            "image": image,
            "resourceRequirements": requirements,
            "command": ["python", "-c", CHECK_PY],
            "logConfiguration": {"logDriver": "awslogs"},
        },
    )
    job_def_arn = reg["jobDefinitionArn"]
    print(f"jobDefinitionArn={job_def_arn}", flush=True)

    sub = batch.submit_job(
        jobName="seekrflow-check-seekr-info",
        jobQueue=job_queue,
        jobDefinition=job_def_arn,
        timeout={"attemptDurationSeconds": timeout_seconds},
    )
    job_id = sub["jobId"]
    print(f"Submitted jobId={job_id}", flush=True)

    terminal = {"SUCCEEDED", "FAILED"}
    job = None
    while True:
        desc = batch.describe_jobs(jobs=[job_id])["jobs"][0]
        status = desc.get("status")
        print(f"  status={status}", flush=True)
        if status in terminal:
            job = desc
            break
        time.sleep(15)

    # Print CloudWatch log tail if present
    container = (job or {}).get("container") or {}
    log_stream = container.get("logStreamName")
    if log_stream:
        print(f"\n=== CloudWatch log stream: {log_stream} ===", flush=True)
        try:
            events = logs.get_log_events(
                logGroupName="/aws/batch/job",
                logStreamName=log_stream,
                startFromHead=True,
            ).get("events") or []
            for ev in events:
                print(ev.get("message", ""), end="" if str(ev.get("message", "")).endswith("\n") else "\n")
        except Exception as e:
            print(f"(could not fetch logs: {e})", flush=True)

    if (job or {}).get("status") == "SUCCEEDED":
        # Job success means CHECK_PY exited 0 ⇒ info present
        print("\nOK: Batch job succeeded — 'info' is available in this image.", flush=True)
        return 0
    print("\nFAIL: Batch job did not succeed (image likely lacks 'info', or job error).", flush=True)
    reason = container.get("reason") or (job or {}).get("statusReason")
    if reason:
        print(f"reason: {reason}", flush=True)
    return 1


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--mode", choices=("docker", "batch"), default="docker",
        help="How to run the check (default: docker)",
    )
    p.add_argument(
        "--image",
        default=None,
        help="Full image URI (default: <account>.dkr.ecr.<region>.amazonaws.com/seekr-engines-bd:latest)",
    )
    p.add_argument("--region", default="us-west-2")
    p.add_argument("--account-id", default=None, help="Required unless --image is given")
    p.add_argument("--job-queue", default="seekrflow-gpu-queue")
    p.add_argument("--job-def-name", default="seekrflow-check-seekr-info")
    p.add_argument("--gpus", type=int, default=1)
    p.add_argument("--vcpus", type=int, default=4)
    p.add_argument("--memory-mb", type=int, default=14000)
    p.add_argument("--timeout-seconds", type=int, default=600)
    args = p.parse_args(argv)

    image = args.image
    if not image:
        if not args.account_id:
            p.error("--account-id is required when --image is omitted")
        image = (
            f"{args.account_id}.dkr.ecr.{args.region}.amazonaws.com/"
            f"seekr-engines-bd:latest"
        )

    print(f"Checking image: {image}", flush=True)
    if args.mode == "docker":
        return check_via_docker(image)
    return check_via_batch(
        image=image,
        region=args.region,
        job_queue=args.job_queue,
        job_def_name=args.job_def_name,
        n_gpus=args.gpus,
        n_vcpus=args.vcpus,
        memory_mb=args.memory_mb,
        timeout_seconds=args.timeout_seconds,
    )


if __name__ == "__main__":
    sys.exit(main())
