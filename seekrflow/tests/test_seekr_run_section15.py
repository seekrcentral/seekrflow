"""
Section 15 tests: resubmit remote jobs that time out / die before completion.
"""
from __future__ import annotations

import asyncio
import sys
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

if "fabric" not in sys.modules:
    sys.modules["fabric"] = MagicMock()


def _ensure_radical_mock() -> None:
    if "radical.asyncflow" not in sys.modules:
        radical = MagicMock()
        asyncflow = MagicMock()
        radical.asyncflow = asyncflow
        sys.modules["radical"] = radical
        sys.modules["radical.asyncflow"] = asyncflow


def _idle_incomplete_status(
        *,
        progress: float = 0.5,
        state: str = "started",
        ) -> dict:
    return {
        "success": True,
        "manager_status": {"jobs": []},
        "stage_status": {"state": state, "progress": progress},
    }


def _make_remote_stage_workflow(
        seekr_run,
        *,
        name: str = "bd",
        fusion_host: str | None = None,
        fused_after: list[str] | None = None,
        last_progress: float = 0.0,
        subsequent_noncompleted_runs: int = 0,
        ):
    from seekrflow.modules import structures

    sw = seekr_run.StageWorkflow(
        model=SimpleNamespace(directory="/work"),
        seekrflow=SimpleNamespace(name="run1"),
        stage=SimpleNamespace(name=name, index=1),
        workflow_engine=object(),
        resource_name="cluster",
        resource=structures.Resource_remote_slurm(name="cluster"),
    )
    sw.fusion_host = fusion_host
    sw.fused_after = list(fused_after or [])
    sw.last_progress = last_progress
    sw.subsequent_noncompleted_runs = subsequent_noncompleted_runs
    sw.state = "started"
    sw.semaphore = "go"
    sw.peer_workflows = {}
    return sw


async def _run_monitor_with_statuses(
        sw,
        monkeypatch,
        statuses: list[dict],
        ) -> int:
    from seekrflow.modules import seekr_run
    from seekrflow.modules.workload_managers import remote as workload_remote

    calls = {"n": 0}

    def fake_status(*args, **kwargs):
        idx = calls["n"]
        calls["n"] += 1
        if idx >= len(statuses):
            sw.detached_requested = True
            return statuses[-1]
        return statuses[idx]

    monkeypatch.setattr(workload_remote, "status_remote", fake_status)
    monkeypatch.setattr(seekr_run, "POLLING_INTERVAL", 0)
    await sw._monitor_stage_loop()
    return calls["n"]


@pytest.fixture
def seekr_run_module():
    _ensure_radical_mock()
    from seekrflow.modules import seekr_run
    return seekr_run


class TestRemoteTimeoutResubmit:
    def test_job_gone_with_progress_resubmits(self, monkeypatch, seekr_run_module):
        sw = _make_remote_stage_workflow(
            seekr_run_module, last_progress=0.0)
        statuses = [
            _idle_incomplete_status(progress=0.5),
            _idle_incomplete_status(progress=0.5),
            _idle_incomplete_status(progress=0.5),
        ]
        n_calls = asyncio.run(
            _run_monitor_with_statuses(sw, monkeypatch, statuses))
        assert n_calls == 3
        assert sw.state == "unstarted"
        assert sw.semaphore == "go"
        assert sw.subsequent_noncompleted_runs == 0
        assert sw.last_progress == 0.5

    def test_job_gone_without_progress_gives_up(self, monkeypatch, seekr_run_module):
        sw = _make_remote_stage_workflow(
            seekr_run_module,
            last_progress=0.2,
            subsequent_noncompleted_runs=3,
        )
        statuses = [
            _idle_incomplete_status(progress=0.2),
            _idle_incomplete_status(progress=0.2),
            _idle_incomplete_status(progress=0.2),
        ]
        asyncio.run(_run_monitor_with_statuses(sw, monkeypatch, statuses))
        assert sw.state == "failed"
        assert sw.semaphore == "wait"
        assert sw.subsequent_noncompleted_runs == 4

    def test_single_idle_check_does_not_resubmit(
            self, monkeypatch, seekr_run_module):
        sw = _make_remote_stage_workflow(seekr_run_module)
        statuses = [
            _idle_incomplete_status(progress=0.1),
            {
                "success": True,
                "manager_status": {
                    "jobs": [{"JobID": "99", "State": "RUNNING"}],
                },
                "stage_status": {"state": "started", "progress": 0.1},
            },
            _idle_incomplete_status(progress=0.1),
            _idle_incomplete_status(progress=0.1),
            _idle_incomplete_status(progress=0.1),
        ]
        n_calls = asyncio.run(
            _run_monitor_with_statuses(sw, monkeypatch, statuses))
        assert n_calls == 5
        assert sw.state == "unstarted"
        assert sw.semaphore == "go"

    def test_fused_stage_does_not_self_resubmit(
            self, monkeypatch, seekr_run_module):
        host = _make_remote_stage_workflow(seekr_run_module, name="host")
        member = _make_remote_stage_workflow(
            seekr_run_module, fusion_host="host", name="member")
        host.fused_after = ["member"]
        by_name = {"host": host, "member": member}
        host.peer_workflows = by_name
        member.peer_workflows = by_name
        statuses = [
            _idle_incomplete_status(progress=0.1),
            _idle_incomplete_status(progress=0.1),
            _idle_incomplete_status(progress=0.1),
            _idle_incomplete_status(progress=0.1),
        ]
        n_calls = asyncio.run(
            _run_monitor_with_statuses(member, monkeypatch, statuses))
        assert n_calls >= 4
        assert member.state == "started"
        assert member.semaphore == "go"
        assert member.subsequent_noncompleted_runs == 0

    def test_fused_host_waits_for_set_completion(
            self, monkeypatch, seekr_run_module):
        host = _make_remote_stage_workflow(
            seekr_run_module, name="host", fused_after=["member"])
        member = _make_remote_stage_workflow(
            seekr_run_module, fusion_host="host", name="member")
        by_name = {"host": host, "member": member}
        host.peer_workflows = by_name
        member.peer_workflows = by_name
        member.state = "started"
        statuses = [
            {
                "success": True,
                "manager_status": {"jobs": []},
                "stage_status": {"state": "completed", "progress": 1.0},
            },
            {
                "success": True,
                "manager_status": {"jobs": []},
                "stage_status": {"state": "completed", "progress": 1.0},
            },
            {
                "success": True,
                "manager_status": {"jobs": []},
                "stage_status": {"state": "completed", "progress": 1.0},
            },
        ]
        n_calls = asyncio.run(
            _run_monitor_with_statuses(host, monkeypatch, statuses))
        assert n_calls >= 3
        assert member.state == "started"

    def test_fused_host_exits_when_set_complete(
            self, monkeypatch, seekr_run_module):
        host = _make_remote_stage_workflow(
            seekr_run_module, name="host", fused_after=["member"])
        member = _make_remote_stage_workflow(
            seekr_run_module, fusion_host="host", name="member")
        by_name = {"host": host, "member": member}
        host.peer_workflows = by_name
        member.peer_workflows = by_name
        member.state = "completed"
        statuses = [
            {
                "success": True,
                "manager_status": {"jobs": []},
                "stage_status": {"state": "completed", "progress": 1.0},
            },
        ]
        n_calls = asyncio.run(
            _run_monitor_with_statuses(host, monkeypatch, statuses))
        assert n_calls == 1
        assert host.task is None

    def test_fused_host_resubmits_when_member_incomplete(
            self, monkeypatch, seekr_run_module):
        host = _make_remote_stage_workflow(
            seekr_run_module, name="host", fused_after=["member"])
        member = _make_remote_stage_workflow(
            seekr_run_module, fusion_host="host", name="member")
        by_name = {"host": host, "member": member}
        host.peer_workflows = by_name
        member.peer_workflows = by_name
        host.state = "completed"
        member.state = "started"
        member.progress = 0.3
        statuses = [
            _idle_incomplete_status(progress=1.0, state="completed"),
            _idle_incomplete_status(progress=1.0, state="completed"),
            _idle_incomplete_status(progress=1.0, state="completed"),
        ]
        asyncio.run(_run_monitor_with_statuses(host, monkeypatch, statuses))
        assert host.state == "unstarted"
