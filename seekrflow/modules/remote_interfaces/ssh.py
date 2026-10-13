"""
modules/remote_interfaces/ssh.py

Provide workflow submission over SSH via Fabric.
"""

import typing
import fabric
from io import StringIO

import seekrflow.modules.remote_interfaces.base as remote_interface_base

# Timeout for a SSH promise future. On expiry, raise TimeoutError
# so callers can retry on the next poll without treating it as a hard failure.
SSH_RESULT_TIMEOUT_S = 120.0

def submit_remote_workload_with_ssh(
        resource_name: str,
        workload: typing.Any,
        manager_payload: dict,
        hostname: str,
        username: str = "",
        password: str = "",
        port: int = 22,
        private_key_filename: str | None = None,
        private_key_passphrase: str | None = None,
        ) -> dict:
    args = [manager_payload]

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
    cmd = remote_interface_base.build_python_c_command(workload, args)
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
        timeout=SSH_RESULT_TIMEOUT_S,
    )
    # Track how much we’ve already printed from the buffers
    seen_out = seen_err = 0

    try:
        # Process finished — grab the final Result (exit code, full buffers, etc.)
        try:
            promise.join()
        except Exception:
            so = out_buf.getvalue()
            se = err_buf.getvalue()
            raise RuntimeError(
                "SSH command failed.\n"
                f"stdout:\n{so[-4000:]}\n"
                f"stderr:\n{se[-4000:]}"
            ) from None
            
        so = out_buf.getvalue()
        if len(so) > seen_out:
            seen_out = len(so)

        se = err_buf.getvalue()
        if len(se) > seen_err:
            seen_err = len(se)

        val = remote_interface_base.parse_workflow_stdout(
            so, se, transport_label="SSH")
        c.close()
        return val

    finally:
        # Best-effort cleanup if the loop exits unexpectedly
        if not promise.runner.process_is_finished:
            promise.runner.stop()
        c.close()
