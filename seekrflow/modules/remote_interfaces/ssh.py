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
        run_workflow: typing.Any,
        cancel_workflow: typing.Any,
        status_workflow: typing.Any,
        args: tuple,
        hostname: str,
        username: str = "",
        password: str = "",
        port: int = 22,
        ) -> bool:
    
    INTERVAL = 0.5
    stop_flag = args[9]

    c = fabric.Connection(host=hostname, user=username, port=port,
                          connect_kwargs={"password": password})
    touch_cmd = "mkdir -p \"$(dirname '{path}')\" && touch '{path}'"
    args_list = []
    for arg in args:
        if isinstance(arg, str):
            args_list.append(f"\"{arg}\"")
        elif isinstance(arg, list):
            args_list.append(f"[{','.join(map(str, arg))}]")
        else:
            args_list.append(f"{arg}")

    arg_str = "(" + ",".join(args_list) + ")"
    run_workflow_source_code = inspect.getsource(run_workflow)
    # NOTE: can we always assume that the wanted function is in the first line?
    #  Answer: yes - but we have to be intentional about coding the other workflows
    #  that way.
    run_func_name = run_workflow_source_code.split("\n")[0].split("def ")[1].split("(")[0]
    status_workflow_source_code = inspect.getsource(status_workflow)
    status_func_name = status_workflow_source_code.split("\n")[0].split("def ")[1].split("(")[0]
    cancel_workflow_source_code = inspect.getsource(cancel_workflow)
    cancel_func_name = cancel_workflow_source_code.split("\n")[0].split("def ")[1].split("(")[0]

    run_call_raw_str = f"{run_workflow_source_code}{run_func_name}({arg_str})"
    run_call_str = shlex.quote(run_call_raw_str)
    run_cmd = f"python3 -u -c {run_call_str}"

    status_call_raw_str = f"{status_workflow_source_code}print({status_func_name}({arg_str}))"
    status_call_str = shlex.quote(status_call_raw_str)
    status_cmd = f"python3 -u -c {status_call_str}"

    cancel_call_raw_str = f"{cancel_workflow_source_code}{cancel_func_name}({arg_str})"
    cancel_call_str = shlex.quote(cancel_call_raw_str)
    cancel_cmd = f"python3 -u -c {cancel_call_str}"
    
    # Buffers where Fabric’s I/O threads will write output in real-time
    out_buf_run, err_buf_run = StringIO(), StringIO()
    out_buf_status, err_buf_status = StringIO(), StringIO()
    
    # Run in async mode; don't allocate a PTY so stdout/stderr stay distinct
    promise_run = c.run(
        run_cmd,
        asynchronous=True,     # returns immediately with a Promise
        pty=False,             # keep separate stdout/stderr streams
        hide=True,             # don't echo to our local terminal
        out_stream=out_buf_run,    # capture live stdout here
        err_stream=err_buf_run,    # capture live stderr here
    )

    # Track how much we’ve already printed from the buffers
    seen_out = seen_err = 0

    print("Jobs queued...")
    print(f"If seekrflow is killed or loses connection, SLURM job(s) "
            "should continue to run until: completion, the time limit is "
            "exceeded, or an error is encountered. In such a case, re-running "
            "seekrflow 'run' will re-establish contact with SLURM and reattach.")
    killing_job = False
    file_transfer_back = True
    print(f"Options: print (i)nformation, (d)etach, or (k)ill with file "
            "transfer back or (K)ill without file transfer back: ", end="", 
            flush=True)
    try:
        while not promise_run.runner.process_is_finished:  # poll the remote process
            rlist, _, _ = select.select([sys.stdin], [], [], INTERVAL)
            if rlist:
                user_input = sys.stdin.readline().strip()
                if user_input == "i":
                    
                    result_status = c.run(
                        status_cmd,
                        asynchronous=False,
                        pty=False,
                        out_stream=out_buf_status,
                        err_stream=err_buf_status,
                    )
                    stdout_str = out_buf_status.getvalue()
                    #print("result_status:", result_status.stdout)
                    try:
                        val = ast.literal_eval(stdout_str)
                        if isinstance(val, dict):
                            prettify.prettify_job_info(val, color=True)
                    except Exception as e:
                        print(f"Error: {e}")

                    print(f"Options: print (i)nformation, (d)etach, or (k)ill with file "
                          f"transfer back or (K)ill without file transfer back: ", end="",
                          flush=True)
                elif user_input == "d":
                    print(f"Detaching job...")
                    result_detach = c.run(
                        touch_cmd.format(path=stop_flag),
                        asynchronous=False,     # returns immediately with a Promise
                    )
                    file_transfer_back = False

                elif user_input == "k":
                    print(f"Killing job with file transfer...")
                    killing_job = True

                elif user_input == "K":
                    print(f"Killing job without file transfer...")
                    killing_job = True
                    file_transfer_back = False

                if killing_job:
                    result_status = c.run(
                        cancel_cmd,
                        asynchronous=False,     # returns immediately with a Promise
                    )
                    
        # Process finished — grab the final Result (exit code, full buffers, etc.)
        so = out_buf_run.getvalue()
        if len(so) > seen_out:
            print(so[seen_out:], end="")
            seen_out = len(so)

        se = err_buf_run.getvalue()
        if len(se) > seen_err:
            print(se[seen_err:], end="", file=sys.stderr)
            seen_err = len(se)
        promise_run.join()
        c.close()
        return file_transfer_back

    finally:
        # Best-effort cleanup if the loop exits unexpectedly
        if not promise_run.runner.process_is_finished:
            promise_run.runner.stop()
        c.close()