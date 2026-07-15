"""
Section 16/17 tests: manager drain after completion + bounded shutdown.
"""
from __future__ import annotations

import asyncio
import sys
import time
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


def _completed_with_jobs(*, progress: float = 1.0) -> dict:
    return {
        "success": True,
        "manager_status": {
            "jobs": [{"JobID": "99", "State": "RUNNING"}],
        },
        "stage_status": {"state": "completed", "progress": progress},
    }


def _completed_idle(*, progress: float = 1.0) -> dict:
    return {
        "success": True,
        "manager_status": {"jobs": []},
        "stage_status": {"state": "completed", "progress": progress},
    }


def _make_remote_stage_workflow(seekr_run):
    from seekrflow.modules import structures

    sw = seekr_run.StageWorkflow(
        model=SimpleNamespace(directory="/work"),
        seekrflow=SimpleNamespace(name="run1"),
        stage=SimpleNamespace(name="bd", index=1),
        workflow_engine=object(),
        resource_name="cluster",
        resource=structures.Resource_remote_slurm(name="cluster"),
    )
    sw.state = "started"
    sw.semaphore = "go"
    sw.peer_workflows = {}
    return sw


async def _run_monitor_with_statuses(sw, monkeypatch, statuses: list[dict]) -> int:
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


class TestManagerDrainAfterCompletion:
    def test_drains_scheduler_before_idle(self, monkeypatch, seekr_run_module):
        sw = _make_remote_stage_workflow(seekr_run_module)
        statuses = [
            _completed_with_jobs(),
            _completed_idle(),
        ]
        asyncio.run(_run_monitor_with_statuses(sw, monkeypatch, statuses))
        assert sw.state == "completed"
        assert sw.manager_status == "idle"

    def test_forces_idle_after_drain_cap(self, monkeypatch, seekr_run_module):
        sw = _make_remote_stage_workflow(seekr_run_module)
        statuses = [
            _completed_with_jobs(),
            _completed_with_jobs(),
            _completed_with_jobs(),
        ]
        asyncio.run(_run_monitor_with_statuses(sw, monkeypatch, statuses))
        assert sw.state == "completed"
        assert sw.manager_status == "idle"


class TestBoundedShutdown:
    def test_shutdown_times_out_and_force_exits(
            self, monkeypatch, seekr_run_module):
        pipeline = seekr_run_module.SeekrPipeline.__new__(
            seekr_run_module.SeekrPipeline)
        pipeline._shutting_down = False
        pipeline._stop_event = None
        pipeline.stage_workflows = []
        monkeypatch.setattr(
            seekr_run_module.SeekrPipeline,
            "write_status_snapshot",
            lambda self: None,
        )

        exits: list[int] = []

        class FakeEngine:
            async def shutdown(self):
                return None

        pipeline.workflow_engine = FakeEngine()
        monkeypatch.setattr(
            seekr_run_module.SeekrPipeline,
            "_cancel_all_stages",
            lambda self: time.sleep(0.2),
        )
        monkeypatch.setattr(seekr_run_module, "SHUTDOWN_CANCEL_TIMEOUT", 0.05)

        def fake_exit(code):
            exits.append(code)
            raise SystemExit(code)

        monkeypatch.setattr(seekr_run_module.os, "_exit", fake_exit)

        with pytest.raises(SystemExit):
            asyncio.run(pipeline.shutdown())
        assert exits == [1]

    def test_second_ctrl_c_force_exits(self, monkeypatch, seekr_run_module):
        pipeline = seekr_run_module.SeekrPipeline.__new__(
            seekr_run_module.SeekrPipeline)
        pipeline._shutting_down = True
        pipeline._shutdown_task = None
        exits: list[int] = []

        def fake_exit(code):
            exits.append(code)
            raise SystemExit(code)

        monkeypatch.setattr(seekr_run_module.os, "_exit", fake_exit)

        with pytest.raises(SystemExit):
            pipeline._on_signal()
        assert exits == [1]
