"""
modules/transfer/aws_s3.py

Handle file transfers via AWS S3 - the required method for Amazon web services.
"""

import os
import subprocess
import tempfile

# Workload-manager state that must stay on the client (or be recreated).
# Absolute symlinks under these dirs break container ``tar`` extract.
LOCAL_RUNNER_DIRNAMES = (
    ".multiprocessing",
    ".aws_runner",
    ".slurm_runner",
    ".pbs_runner",
)


def _run(cmd: list[str]) -> None:
    print(f"$ {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        parts = [f"Command {cmd!r} failed with exit status {result.returncode}"]
        if result.stderr and result.stderr.strip():
            parts.append(f"stderr:\n{result.stderr.rstrip()}")
        if result.stdout and result.stdout.strip():
            parts.append(f"stdout:\n{result.stdout.rstrip()}")
        raise RuntimeError("\n".join(parts))


def build_model_tarball_command(
        tarball_path: str,
        local_path: str,
        ) -> list[str]:
    """``tar czf`` command that omits local runner-state directories."""
    cmd = ["tar", "czf", tarball_path]
    for name in LOCAL_RUNNER_DIRNAMES:
        cmd.append(f"--exclude=./{name}")
    cmd.extend(["-C", local_path, "."])
    return cmd


def transfer_files_with_aws_s3(
        local_path: str,
        s3_uri: str,
        input_tarball_name: str,
        region: str,
        backwards: bool,
        ) -> None:
    """Tar the local model directory and upload it as a single S3 object.

    NOTE: 'aws s3 sync' silently drops empty directories because S3 has no
    real directories -- only objects (keys). The prepared SEEKR calculation
    relies on empty dirs (e.g. .../stage_1_equilibration/ where run.log is
    written), so we package the tree as a tarball, which preserves empty
    directories, and extract it inside the container.
    """
    s3_uri = s3_uri.rstrip("/")
    if backwards:
        os.makedirs(local_path, exist_ok=True)
        _run(["aws", "s3", "sync", s3_uri, local_path,
              "--region", region, "--only-show-errors"])
        downloaded_tarball = os.path.join(local_path, input_tarball_name)
        if os.path.isfile(downloaded_tarball):
            os.remove(downloaded_tarball)
    else:
        tarball_path = os.path.join(tempfile.gettempdir(), input_tarball_name)
        _run(build_model_tarball_command(tarball_path, local_path))
        _run(["aws", "s3", "cp", tarball_path, f"{s3_uri}/{input_tarball_name}",
              "--region", region, "--only-show-errors"])
        os.remove(tarball_path)
        print("Transfer complete!")
