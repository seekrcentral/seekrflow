"""
modules/remote_interfaces/globus_compute_sdk.py

Provide workflow submission with globus compute SDK.
"""

import sys
import ast
import shlex
import select
import typing
import fabric
import inspect
from io import StringIO
import seekrflow.modules.prettify as prettify

def submit_remote_workflow_with_ssh(
        name: str,
        workflow: typing.Any,
        args: tuple,
        hostname: str,
        username: str = "",
        password: str = "",
        port: int = 22,
        ) -> bool:
    c = fabric.Connection(host=hostname, user=username, port=port,
                          connect_kwargs={"password": password})
    args_list = []
    for arg in args:
        if isinstance(arg, str):
            args_list.append(f"\"{arg}\"")
        elif isinstance(arg, list):
            args_list.append(f"[{','.join(map(str, arg))}]")
        else:
            args_list.append(f"{arg}")

    arg_str = "(" + ",".join(args_list) + ")"
    workflow_source_code = inspect.getsource(workflow)
    # NOTE: can we always assume that the wanted function is in the first line?
    #  Answer: yes - but we have to be intentional about coding the other workflows
    #  that way.
    func_name = run_workflow_source_code.split("\n")[0].split("def ")[1].split("(")[0]
    call_raw_str = f"{workflow_source_code}{func_name}({arg_str})"
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
        so = out_buf.getvalue()
        if len(so) > seen_out:
            print(so[seen_out:], end="")
            seen_out = len(so)

        se = err_buf.getvalue()
        if len(se) > seen_err:
            print(se[seen_err:], end="", file=sys.stderr)
            seen_err = len(se)
        result = promise.join()
        c.close()
        return 

    finally:
        # Best-effort cleanup if the loop exits unexpectedly
        if not promise.runner.process_is_finished:
            promise.runner.stop()
        c.close()
