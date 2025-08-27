"""
modules/remote_interfaces/globus_compute_sdk.py

Provide workflow submission with globus compute SDK.
"""

import sys
import select
import typing
import seekrflow.modules.prettify as prettify

def submit_remote_workflow_with_globus_compute(
        run_workflow: typing.Any,
        cancel_workflow: typing.Any,
        status_workflow: typing.Any,
        endpoint: str,
        args: tuple,
) -> bool:
    from globus_compute_sdk import Client, Executor, ShellFunction
    from globus_compute_sdk.serialize import ComputeSerializer, CombinedCode
    INTERVAL = 0.5

    c = Client()
    status = c.get_endpoint_status(endpoint)
    s = (status.get("status") or "unknown").lower()
    d = status.get("details") or {}
    idle = d.get("idle_workers")
    total = d.get("total_workers")
    pending = d.get("pending_tasks") or d.get("outstanding_tasks")

    if s == "online":
        if total is not None and idle == 0:
            print(
                "WARNING: Endpoint is ONLINE but currently has 0 idle workers "
                f"(total_workers={total}, pending_tasks={pending or 0}). "
                "Your task will queue until workers become available."
            )
    elif s in {"offline", "disconnected", "stopped"}:
        print(
            f"WARNING: Endpoint is {s.upper()}. Start it on the host: "
            "`globus-compute-endpoint start <ENDPOINT_NAME>`."
        )
    elif s in {"initializing", "starting"}:
        print("Endpoint is STARTING; tasks will queue until workers connect.")

    else:
        print("WARNING: Endpoint is in an UNKNOWN state. Jobs will not run.")

    stop_flag = args[9]
    touch = ShellFunction("mkdir -p \"$(dirname '{path}')\" && touch '{path}'")

    with Executor(endpoint) as gcx:
        gcx.serializer = ComputeSerializer(strategy_code=CombinedCode())
        function_id = gcx.register_function(run_workflow, description="run")
        future = gcx.submit_to_registered_function(function_id=function_id, 
                                                   args=(args,))
        function_id_status = gcx.register_function(status_workflow, 
                                                   description="status")
        function_id_scancel = gcx.register_function(cancel_workflow, 
                                                    description="cancelling")
        print("Jobs queued...")
        print(f"If seekrflow is killed or loses connection, SLURM job(s) "
              "should continue to run until: completion, the time limit is "
              "exceeded, or an error is encountered. In such a case, re-running "
              "seekrflow 'run' will re-establish contact with SLURM and reattach.")
        
        keep_running = True
        killing_job = False
        file_transfer_back = True
        print(f"Options: print (i)nformation, (d)etach, or (k)ill with file "
              "transfer back or (K)ill without file transfer back: ", end="", 
              flush=True)
        while keep_running:
            rlist, _, _ = select.select([sys.stdin], [], [], INTERVAL)
            if rlist:
                user_input = sys.stdin.readline().strip()
                if user_input == "i":
                    
                    future_status = gcx.submit_to_registered_function(
                        function_id=function_id_status, args=(args,))
                    result = future_status.result()
                    if not result:
                        print("No job information available.")
                    elif "rows" not in result:
                        print("No job information available.")
                        if "std_err" in result:
                            print(f"STDERR: {result['std_err']}")
                    else:
                        prettify.prettify_job_info(result, color=True)
                    print(f"Options: print (i)nformation, (d)etach, or (k)ill with file "
                          f"transfer back or (K)ill without file transfer back: ", end="",
                          flush=True)
                elif user_input == "d":
                    print(f"Detaching job...")
                    gcx.submit(touch, path=stop_flag)
                    keep_running = False
                    file_transfer_back = False

                elif user_input == "k":
                    print(f"Killing job with file transfer...")
                    killing_job = True

                elif user_input == "K":
                    print(f"Killing job without file transfer...")
                    killing_job = True
                    file_transfer_back = False

                if killing_job:
                    
                    future_scancel = gcx.submit_to_registered_function(
                        function_id=function_id_scancel, args=(args,))
                    future_scancel.result()
                    keep_running = False

                if future.done():
                    print("Task is DONE.")
                    keep_running = False

    future.result()
    return file_transfer_back