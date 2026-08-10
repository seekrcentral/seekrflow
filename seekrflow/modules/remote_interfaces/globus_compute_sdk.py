"""
modules/remote_interfaces/globus_compute_sdk.py

Provide workflow submission with globus compute SDK.
"""

from __future__ import annotations

import time
import typing

# Soft wait for a submitted Globus Compute future. On expiry, raise TimeoutError
# so callers can retry on the next poll without treating it as a hard failure.
GLOBUS_RESULT_TIMEOUT_S = 1200.0

# Dedupe endpoint health warnings: same (endpoint, warn_key) at most once per
# this many seconds; a different warn_key prints immediately.
_ENDPOINT_WARN_INTERVAL_S = 300.0
_last_endpoint_warn: dict[str, tuple[str, float]] = {}


def _warn_endpoint(name: str, warn_key: str, message: str) -> None:
    now = time.monotonic()
    previous = _last_endpoint_warn.get(name)
    if (
            previous is not None
            and previous[0] == warn_key
            and (now - previous[1]) < _ENDPOINT_WARN_INTERVAL_S
    ):
        return
    _last_endpoint_warn[name] = (warn_key, now)
    print(message)


def submit_remote_workflow_with_globus_compute(
        name: str,
        workflow: typing.Any,
        endpoint: str,
        args: tuple,
        silent: bool = False,
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

    if not silent:
        if s == "online":
            if total is not None and idle == 0:
                _warn_endpoint(
                    name,
                    "online_idle0",
                    f"WARNING: Globus endpoint '{name}' is ONLINE but currently "
                    f"has 0 idle workers (total_workers={total}, "
                    f"pending_tasks={pending or 0}). Your task will queue until "
                    f"workers become available.",
                )
        elif s in {"offline", "disconnected", "stopped"}:
            _warn_endpoint(
                name,
                f"state_{s}",
                f"WARNING: Globus endpoint '{name}' is {s.upper()}. Start it "
                f"on the host: `globus-compute-endpoint start <ENDPOINT_NAME>`.",
            )
        elif s in {"initializing", "starting"}:
            _warn_endpoint(
                name,
                f"state_{s}",
                f"Globus endpoint '{name}' is STARTING; tasks will queue until "
                f"workers connect.",
            )
        else:
            _warn_endpoint(
                name,
                f"state_{s}",
                f"WARNING: Globus endpoint '{name}' is in an UNKNOWN state. "
                f"Jobs will not run.",
            )

    with Executor(endpoint) as gcx:
        gcx.serializer = ComputeSerializer(strategy_code=CombinedCode())
        function_id = gcx.register_function(workflow, description="run")
        future = gcx.submit_to_registered_function(
            function_id=function_id, args=(args,))

    try:
        return future.result(timeout=GLOBUS_RESULT_TIMEOUT_S)
    except TimeoutError as e:
        raise TimeoutError(
            f"Globus Compute task on endpoint {name!r} exceeded "
            f"{GLOBUS_RESULT_TIMEOUT_S:.0f}s waiting for a result; "
            f"will retry on the next poll"
        ) from e
