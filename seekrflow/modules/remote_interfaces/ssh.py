"""
modules/remote_interfaces/globus_compute_sdk.py

Provide workflow submission with globus compute SDK.
"""

import ast
import shlex
import typing
import fabric
import inspect
from io import StringIO

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
        ) -> bool:
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
    args_list = []
    for arg in args:
        # Use Python literal representation so strings with quotes, dicts,
        # None, and booleans round-trip safely through python -c.
        args_list.append(repr(arg))
    if len(args_list) == 1:
        arg_str = "(" + args_list[0] + ",)"
    else:
        arg_str = "(" + ",".join(args_list) + ")"
    
    workflow_source_code = inspect.getsource(workflow)
    # NOTE: can we always assume that the wanted function is in the first line?
    #  Answer: yes - but we have to be intentional about coding the other workflows
    #  that way.
    func_name = workflow_source_code.split("\n")[0].split("def ")[1].split("(")[0]
    call_raw_str = f"{workflow_source_code}print({func_name}({arg_str}))"
    call_str = shlex.quote(call_raw_str)
    cmd = f"python3 -u -c {call_str}"
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
        result = promise.join()
        so = out_buf.getvalue()
        if len(so) > seen_out:
            #print(so[seen_out:], end="")
            seen_out = len(so)

        se = err_buf.getvalue()
        if len(se) > seen_err:
            #print(se[seen_err:], end="", file=sys.stderr)
            seen_err = len(se)
        
        val = {}
        # Workflows may print diagnostics before the return payload.
        # Parse the last non-empty line as the structured return value.
        lines = [line.strip() for line in so.splitlines() if line.strip()]
        try:
            if lines:
                val = ast.literal_eval(lines[-1])
            else:
                val = ast.literal_eval(so)
        except Exception as e:
            stderr_tail = "\n".join(se.splitlines()[-10:])
            stdout_tail = "\n".join(so.splitlines()[-10:])
            raise RuntimeError(
                "Failed to parse remote workflow result via SSH. "
                f"parse_error={e}. stdout_tail={stdout_tail!r}. "
                f"stderr_tail={stderr_tail!r}"
            )
        c.close()
        return val

    finally:
        # Best-effort cleanup if the loop exits unexpectedly
        if not promise.runner.process_is_finished:
            promise.runner.stop()
        c.close()
