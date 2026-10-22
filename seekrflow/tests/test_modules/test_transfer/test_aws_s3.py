"""
Tests for seekrflow/modules/transfer/aws_s3.py
"""

import subprocess
import tarfile

from seekrflow.modules.transfer import aws_s3

def test_build_model_tarball_excludes_runner_dirs(tmp_path):
    root = tmp_path / "root"
    root.mkdir()
    (root / "model.json").write_text("{}")
    (root / "stage_a").mkdir()
    for name in aws_s3.LOCAL_RUNNER_DIRNAMES:
        runner = root / name
        runner.mkdir()
        (runner / "bd_latest.json").symlink_to(
            "/nonexistent/host/path/bd_state.json")

    tarball = tmp_path / "model_input.tar.gz"
    cmd = aws_s3.build_model_tarball_command(str(tarball), str(root))
    for name in aws_s3.LOCAL_RUNNER_DIRNAMES:
        assert f"--exclude=./{name}" in cmd
    subprocess.run(cmd, check=True)

    names = set(tarfile.open(tarball, "r:gz").getnames())
    assert "./model.json" in names or "model.json" in names
    assert not any(
        dirname in n
        for n in names
        for dirname in aws_s3.LOCAL_RUNNER_DIRNAMES
    )