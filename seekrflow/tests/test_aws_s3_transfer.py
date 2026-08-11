"""Tests for AWS S3 transfer helpers."""
from __future__ import annotations

import os
import subprocess
import tarfile

from seekrflow.modules.transfer import aws_s3


def test_build_model_tarball_excludes_runner_dirs(tmp_path):
    root = tmp_path / "root"
    root.mkdir()
    (root / "model.json").write_text("{}")
    (root / "stage_a").mkdir()
    (root / ".multiprocessing").mkdir()
    latest = root / ".multiprocessing" / "bd_latest.json"
    latest.symlink_to("/nonexistent/host/path/bd_state.json")
    (root / ".aws_runner").mkdir()
    (root / ".aws_runner" / "x_status.json").write_text("{}")

    tarball = tmp_path / "model_input.tar.gz"
    cmd = aws_s3.build_model_tarball_command(str(tarball), str(root))
    assert "--exclude=./.multiprocessing" in cmd
    subprocess.run(cmd, check=True)

    names = set(tarfile.open(tarball, "r:gz").getnames())
    assert "./model.json" in names or "model.json" in names
    assert not any(".multiprocessing" in n for n in names)
    assert not any(".aws_runner" in n for n in names)


def test_local_runner_dirnames_covers_known_managers():
    assert ".multiprocessing" in aws_s3.LOCAL_RUNNER_DIRNAMES
    assert ".aws_runner" in aws_s3.LOCAL_RUNNER_DIRNAMES
    assert ".slurm_runner" in aws_s3.LOCAL_RUNNER_DIRNAMES
    assert ".pbs_runner" in aws_s3.LOCAL_RUNNER_DIRNAMES
