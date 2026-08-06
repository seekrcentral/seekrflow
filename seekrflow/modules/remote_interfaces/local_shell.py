"""
modules/remote_interfaces/local_shell.py

Provide workflow submission on the local machine (e.g. an HPC login node)
via subprocess, mirroring the SSH python -c transport without Fabric.
"""

from __future__ import annotations

import subprocess
import typing

import seekrflow.modules.remote_interfaces.python_c_runner as python_c_runner


def submit_remote_workflow_with_local_shell(
        name: str,
        workflow: typing.Any,
        args: list | tuple,
        python_executable: str = "python3",
        silent: bool = False,
        ) -> dict:
    """
    Run ``workflow`` locally with the same serialization/parse path as SSH.

    ``name`` is accepted for API parity with other remote interfaces; it is
    unused except for optional diagnostics.
    """
    del silent  # reserved for future quiet logging parity with Globus
    cmd = python_c_runner.build_python_c_command(
        workflow, args, python_executable=python_executable)
    completed = subprocess.run(
        cmd,
        shell=True,
        capture_output=True,
        text=True,
    )
    stdout = completed.stdout or ""
    stderr = completed.stderr or ""
    if completed.returncode != 0:
        stderr_tail = "\n".join(stderr.splitlines()[-10:])
        stdout_tail = "\n".join(stdout.splitlines()[-10:])
        raise RuntimeError(
            f"Local shell workflow {name!r} exited with code "
            f"{completed.returncode}. stdout_tail={stdout_tail!r}. "
            f"stderr_tail={stderr_tail!r}"
        )
    return python_c_runner.parse_workflow_stdout(
        stdout, stderr, transport_label="local_shell")
