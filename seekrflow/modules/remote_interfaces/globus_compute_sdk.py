"""
modules/remote_interfaces/globus_compute_sdk.py

Provide workflow submission with globus compute SDK.
"""

import typing

def submit_remote_workflow_with_globus_compute(
        name: str,
        workflow: typing.Any,
        endpoint: str,
        args: tuple,
        silent: bool = False
) -> dict:
    from globus_compute_sdk import Client, Executor
    from globus_compute_sdk.serialize import ComputeSerializer, CombinedCode
    c = Client()
    status = c.get_endpoint_status(endpoint)
    s = (status.get("status") or "unknown").lower()
    d = status.get("details") or {}
    idle = d.get("idle_workers")
    total = d.get("total_workers")
    pending = d.get("pending_tasks") or d.get("outstanding_tasks")
    
    if s == "online":
        if total is not None and idle == 0:
            if not silent:
                print(
                    f"WARNING: Globus endpoint '{name}' is ONLINE but currently has 0 idle workers "
                    f"(total_workers={total}, pending_tasks={pending or 0}). "
                    "Your task will queue until workers become available."
                )
    elif s in {"offline", "disconnected", "stopped"}:
        if not silent:
            print(
                f"WARNING: Globus endpoint '{name}' is {s.upper()}. Start it on the host: "
                "`globus-compute-endpoint start <ENDPOINT_NAME>`."
            )
    elif s in {"initializing", "starting"}:
        if not silent:
            print(f"Globus endpoint '{name}' is STARTING; tasks will queue until workers connect.")

    else:
        if not silent:
            print(f"WARNING: Globus endpoint '{name}' is in an UNKNOWN state. Jobs will not run.")
    with Executor(endpoint) as gcx:
        gcx.serializer = ComputeSerializer(strategy_code=CombinedCode())
        function_id = gcx.register_function(workflow, description="run")
        future = gcx.submit_to_registered_function(function_id=function_id, 
                                                   args=(args,))
        
    return future.result()
    