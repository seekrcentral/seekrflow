"""
Tests for seekrflow/modules/remote_interfaces/local_shell.py
"""

import sys

import pytest

from seekrflow.modules.remote_interfaces import local_shell


def echo_local_shell_payload(args):
    return {"ok": True, "n": args[0]["n"]}


def failing_local_shell_workflow(args):
    raise RuntimeError("boom")


def test_submit_remote_workflow_with_local_shell_round_trip():
    result = local_shell.submit_remote_workflow_with_local_shell(
        "test",
        echo_local_shell_payload,
        {"n": 42},
        python_executable=sys.executable,
    )
    assert result == {"ok": True, "n": 42}


def test_submit_remote_workflow_with_local_shell_nonzero_exit():
    with pytest.raises(RuntimeError, match="exited with code"):
        local_shell.submit_remote_workflow_with_local_shell(
            "fail",
            failing_local_shell_workflow,
            {},
            python_executable=sys.executable,
        )


def test_submit_remote_workflow_with_local_shell_uses_python_executable():
    with pytest.raises(RuntimeError):
        local_shell.submit_remote_workflow_with_local_shell(
            "missing",
            echo_local_shell_payload,
            {"n": 1},
            python_executable="/nonexistent/seekrflow-python",
        )
