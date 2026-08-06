"""
modules/remote_interfaces/ssh.py

Provide workflow submission over SSH via Fabric.
"""

import typing
import fabric
from io import StringIO

import seekrflow.modules.remote_interfaces.python_c_runner as python_c_runner


def submit_remote_workflow_with_ssh(
        name: str,
        workflow: typing.Any,
        args: tuple,
        hostname: str,
        username: str = "",
        password: str = "",
        port: int = 22,
        private_key_filename: str | None = None,
        private_key_passphrase: str | None = None,
        ) -> dict:
    connect_kwargs = {}
    if password:
        connect_kwargs["password"] = password
    if private_key_filename:
        connect_kwargs["key_filename"] = private_key_filename
    if private_key_passphrase:
        connect_kwargs["passphrase"] = private_key_passphrase
    c = fabric.Connection(
        host=hostname,
        user=username,
        port=port,
        connect_kwargs=connect_kwargs,
    )
    cmd = python_c_runner.build_python_c_command(workflow, args)
    # Buffers where Fabric’s I/O threads will write output in real-time
    out_buf, err_buf = StringIO(), StringIO()

    # Run in async mode; don't allocate a PTY so stdout/stderr stay distinct
    promise = c.run(
        cmd,
        asynchronous=True,     # returns immediately with a Promise
        pty=False,             # keep separate stdout/stderr streams
        hide=True,             # don't echo to our local terminal
        out_stream=out_buf,    # capture live stdout here
        err_stream=err_buf,    # capture live stderr here
    )
    # Track how much we’ve already printed from the buffers
    seen_out = seen_err = 0

    try:
        # Process finished — grab the final Result (exit code, full buffers, etc.)
        promise.join()
        so = out_buf.getvalue()
        if len(so) > seen_out:
            seen_out = len(so)

        se = err_buf.getvalue()
        if len(se) > seen_err:
            seen_err = len(se)

        val = python_c_runner.parse_workflow_stdout(
            so, se, transport_label="SSH")
        c.close()
        return val

    finally:
        # Best-effort cleanup if the loop exits unexpectedly
        if not promise.runner.process_is_finished:
            promise.runner.stop()
        c.close()
