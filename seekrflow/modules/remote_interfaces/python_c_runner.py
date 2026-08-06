"""
modules/remote_interfaces/python_c_runner.py

Shared helpers for running a workflow via ``python -c`` (SSH and local_shell).
"""

from __future__ import annotations

import ast
import inspect
import shlex
import typing


def build_python_c_command(
        workflow: typing.Any,
        args: list | tuple,
        python_executable: str = "python3",
        ) -> str:
    """
    Build a shell command that defines ``workflow`` and prints its return value.

    Args are embedded with ``repr`` so strings, dicts, None, and booleans
    round-trip safely through ``python -c``.
    """
    args_list = [repr(arg) for arg in args]
    if len(args_list) == 1:
        arg_str = "(" + args_list[0] + ",)"
    else:
        arg_str = "(" + ",".join(args_list) + ")"

    workflow_source_code = inspect.getsource(workflow)
    # The wanted function is on the first line; remote workflows must be coded
    # that way.
    func_name = workflow_source_code.split("\n")[0].split("def ")[1].split("(")[0]
    call_raw_str = f"{workflow_source_code}print({func_name}({arg_str}))"
    call_str = shlex.quote(call_raw_str)
    return f"{python_executable} -u -c {call_str}"


def parse_workflow_stdout(
        stdout: str,
        stderr: str = "",
        transport_label: str = "remote",
        ) -> dict:
    """
    Parse the structured return value from workflow stdout.

    Workflows may print diagnostics before the return payload; the last
    non-empty line is treated as the structured return value.
    """
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    try:
        if lines:
            val = ast.literal_eval(lines[-1])
        else:
            val = ast.literal_eval(stdout)
    except Exception as e:
        stderr_tail = "\n".join(stderr.splitlines()[-10:])
        stdout_tail = "\n".join(stdout.splitlines()[-10:])
        raise RuntimeError(
            f"Failed to parse remote workflow result via {transport_label}. "
            f"parse_error={e}. stdout_tail={stdout_tail!r}. "
            f"stderr_tail={stderr_tail!r}"
        ) from e
    if not isinstance(val, dict):
        raise RuntimeError(
            f"Remote workflow via {transport_label} returned non-dict value: "
            f"{type(val).__name__}"
        )
    return val
