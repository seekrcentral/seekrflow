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
"""

REGION = "us-west-2"
ACCOUNT_ID = "628031635908"
SEEKR_IMAGE_URI = f"{ACCOUNT_ID}.dkr.ecr.{REGION}.amazonaws.com/seekr-engines-bd:latest"
COMPUTE_ENV_NAME = "seekrflow-gpu-compute"
JOB_QUEUE_NAME = "seekrflow-gpu-queue"
JOB_DEF_GENERIC = "seekrflow-job-generic"
