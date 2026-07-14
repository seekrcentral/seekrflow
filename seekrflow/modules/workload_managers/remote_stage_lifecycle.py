"""
Pure helpers for remote stage launch probing and scheduler naming.

Kept separate from seekr_run so unit tests do not require radical.asyncflow.
"""
from __future__ import annotations


def remote_scheduler_job_name(seekrflow_name: str, stage_name: str) -> str:
    """Workload-manager job name for a stage launch."""
    return f"{seekrflow_name}_{stage_name}"


def classify_remote_probe_status(status: dict | None) -> str:
    """
    Classify a remote status payload for fresh-launch probing.

    Returns one of ``completed``, ``reattach``, or ``submit``.
    """
    if status is None:
        return "submit"
    stage_status = status.get("stage_status") or {}
    if stage_status.get("state") == "completed":
        return "completed"
    jobs = (status.get("manager_status") or {}).get("jobs") or []
    if jobs:
        return "reattach"
    return "submit"


def owns_scheduler_job(fusion_host: str | None) -> bool:
    """True when this stage may own a scheduler submission (not a fused member)."""
    return fusion_host is None


def remote_cancel_needed(
        status: dict | None,
        *,
        local_state: str,
        tracked_job_ids: set[str] | frozenset[str],
        ) -> bool:
    """
    Whether shutdown should attempt scheduler cancellation for this stage.

    Cancel only when the fresh status shows jobs in the queue, or when the
    status probe failed and we still believe a job may be active locally.
    Completed/failed/idle stages with an empty queue are skipped even if
    ``tracked_job_ids`` contains stale ids.
    """
    jobs = ((status or {}).get("manager_status") or {}).get("jobs") or []
    if jobs:
        return True
    if status is None:
        return local_state == "started" or bool(tracked_job_ids)
    return False
