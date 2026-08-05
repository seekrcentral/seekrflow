"""
Section 14 tests: Ctrl-C remote cancel + pre-submission remote probe.
"""
from __future__ import annotations

import asyncio
import sys
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

if "fabric" not in sys.modules:
    sys.modules["fabric"] = MagicMock()

from seekrflow.modules import structures
from seekrflow.modules.workload_managers import pbs as workload_pbs
from seekrflow.modules.workload_managers import remote as workload_remote
from seekrflow.modules.workload_managers import remote_stage_lifecycle
from seekrflow.modules.workload_managers import slurm as workload_slurm


def _ensure_radical_mock() -> None:
    if "radical.asyncflow" not in sys.modules:
        radical = MagicMock()
        asyncflow = MagicMock()
        radical.asyncflow = asyncflow
        sys.modules["radical"] = radical
        sys.modules["radical.asyncflow"] = asyncflow


class TestClassifyRemoteProbeStatus:
    def test_completed(self):
        status = {
            "stage_status": {"state": "completed", "progress": 1.0},
            "manager_status": {"jobs": []},
        }
        assert remote_stage_lifecycle.classify_remote_probe_status(
            status) == "completed"

    def test_reattach_when_jobs_present(self):
        status = {
            "stage_status": {"state": "started", "progress": 0.5},
            "manager_status": {"jobs": [{"JobID": "42", "State": "RUNNING"}]},
        }
        assert remote_stage_lifecycle.classify_remote_probe_status(
            status) == "reattach"

    def test_submit_when_idle(self):
        status = {
            "stage_status": {"state": "unstarted", "progress": 0.0},
            "manager_status": {"jobs": []},
        }
        assert remote_stage_lifecycle.classify_remote_probe_status(
            status) == "submit"

    def test_submit_on_none(self):
        assert remote_stage_lifecycle.classify_remote_probe_status(
            None) == "submit"


class TestForceOverwriteSkipsLaunchProbe:
    def test_true_when_force(self):
        assert remote_stage_lifecycle.force_overwrite_skips_launch_probe(True)

    def test_false_when_not_force(self):
        assert not remote_stage_lifecycle.force_overwrite_skips_launch_probe(
            False)


class TestRemoteSchedulerJobName:
    def test_name_format(self):
        assert remote_stage_lifecycle.remote_scheduler_job_name(
            "tryp_ben", "bd") == "tryp_ben_bd"


class TestRemoteCancelNeeded:
    def test_jobs_present(self):
        status = {
            "manager_status": {"jobs": [{"JobID": "1"}]},
            "stage_status": {"state": "started"},
        }
        assert remote_stage_lifecycle.remote_cancel_needed(
            status, local_state="started", tracked_job_ids=set())

    def test_completed_no_jobs(self):
        status = {
            "manager_status": {"jobs": []},
            "stage_status": {"state": "completed", "finished": True},
        }
        assert not remote_stage_lifecycle.remote_cancel_needed(
            status, local_state="completed", tracked_job_ids={"999"})

    def test_failed_no_jobs(self):
        status = {
            "manager_status": {"jobs": []},
            "stage_status": {"state": "failed"},
        }
        assert not remote_stage_lifecycle.remote_cancel_needed(
            status, local_state="failed", tracked_job_ids=set())

    def test_unstarted_no_jobs(self):
        status = {
            "manager_status": {"jobs": []},
            "stage_status": {"state": "unstarted"},
        }
        assert not remote_stage_lifecycle.remote_cancel_needed(
            status, local_state="unstarted", tracked_job_ids=set())

    def test_status_none_started_conservative(self):
        assert remote_stage_lifecycle.remote_cancel_needed(
            None, local_state="started", tracked_job_ids=set())

    def test_status_none_completed_no_ids(self):
        assert not remote_stage_lifecycle.remote_cancel_needed(
            None, local_state="completed", tracked_job_ids=set())

    def test_stale_ids_but_queue_empty(self):
        status = {
            "manager_status": {"jobs": []},
            "stage_status": {"state": "started"},
        }
        assert not remote_stage_lifecycle.remote_cancel_needed(
            status, local_state="started", tracked_job_ids={"12345"})


class TestOwnsSchedulerJob:
    def test_standalone_or_host(self):
        assert remote_stage_lifecycle.owns_scheduler_job(None)

    def test_fused_member(self):
        assert not remote_stage_lifecycle.owns_scheduler_job("host_stage")


class TestRemoteCancelWorkflows:
    def test_slurm_cancel_by_id_and_name(self, monkeypatch):
        commands: list[str] = []

        def fake_run(cmd, stdout=None, stderr=None, text=None):
            commands.append(" ".join(cmd))
            result = MagicMock()
            result.returncode = 0
            result.stdout = ""
            result.stderr = ""
            return result

        import subprocess
        monkeypatch.setattr(subprocess, "run", fake_run)

        result = workload_slurm.slurm_remote_cancel_workflow(
            ["/work", "12345", "myrun_bd"])
        assert result["success"] is True
        joined = " ".join(commands)
        assert "scancel 12345" in joined
        assert "scancel --name=myrun_bd" in joined

    def test_pbs_cancel_by_name(self, monkeypatch):
        calls: list[list[str]] = []

        def fake_run(cmd, capture_output=None, text=None, timeout=None):
            calls.append(list(cmd))
            result = MagicMock()
            result.returncode = 0
            if cmd[:2] == ["qselect", "-N"]:
                result.stdout = "100 101\n"
            else:
                result.stdout = ""
            result.stderr = ""
            return result

        import subprocess
        monkeypatch.setattr(subprocess, "run", fake_run)

        result = workload_pbs.pbs_remote_cancel_workflow(
            ["/work", "", "myrun_bd"])
        assert result["success"] is True
        assert ["qselect", "-N", "myrun_bd"] in calls
        assert ["qdel", "100"] in calls
        assert ["qdel", "101"] in calls


class TestSubmitRemoteCancelWorkflow:
    def test_cancel_by_name_passes_empty_id_slot(self, monkeypatch):
        captured: dict = {}

        def fake_submit(seekrflow, resource_name, workflow, extra_args=None, silent=False):
            captured["extra_args"] = extra_args
            captured["workflow"] = workflow

        monkeypatch.setattr(
            workload_remote, "submit_remote_workflow", fake_submit)

        seekrflow = SimpleNamespace(
            name="run1",
            run_settings=SimpleNamespace(
                get_resource_by_name=lambda name: structures.Resource_remote_slurm(
                    name="cluster",
                    remote_working_directory="/remote",
                    remote_interface=structures.Remote_interface_ssh(
                        hostname="h", username="u"),
                ),
            ),
        )
        workload_remote.submit_remote_cancel_workflow(
            seekrflow, "cluster", job_name="run1_bd")
        assert captured["extra_args"] == ["", "run1_bd"]
        assert captured["workflow"] is workload_slurm.slurm_remote_cancel_workflow


class TestStageWorkflowKill:
    def test_kill_cancels_by_ids_and_name(self, monkeypatch):
        _ensure_radical_mock()
        from seekrflow.modules import seekr_run

        cancel_calls: list[tuple[str | None, str | None]] = []

        def fake_status(*args, **kwargs):
            return {
                "manager_status": {
                    "jobs": [{"JobID": "999", "State": "RUNNING"}],
                },
            }

        def fake_cancel(
                seekrflow, resource_name, job_id=None, job_name=None, silent=False):
            cancel_calls.append((job_id, job_name))

        monkeypatch.setattr(workload_remote, "status_remote", fake_status)
        monkeypatch.setattr(
            workload_remote, "submit_remote_cancel_workflow", fake_cancel)

        stage = SimpleNamespace(name="bd", index=1)
        seekrflow = SimpleNamespace(name="tryp_ben")
        sw = seekr_run.StageWorkflow(
            model=SimpleNamespace(directory="/work"),
            seekrflow=seekrflow,
            stage=stage,
            workflow_engine=object(),
            resource_name="cluster",
            resource=structures.Resource_remote_slurm(name="cluster"),
        )
        sw.job_ids = {"12345"}

        sw.kill()

        assert ("12345", None) in cancel_calls
        assert ("999", None) in cancel_calls
        assert (None, "tryp_ben_bd") in cancel_calls

    def test_kill_skips_when_completed(self, monkeypatch):
        _ensure_radical_mock()
        from seekrflow.modules import seekr_run

        cancel_calls: list[tuple[str | None, str | None]] = []

        def fake_status(*args, **kwargs):
            return {
                "stage_status": {"state": "completed", "progress": 1.0},
                "manager_status": {"jobs": []},
            }

        def fake_cancel(
                seekrflow, resource_name, job_id=None, job_name=None, silent=False):
            cancel_calls.append((job_id, job_name))

        monkeypatch.setattr(workload_remote, "status_remote", fake_status)
        monkeypatch.setattr(
            workload_remote, "submit_remote_cancel_workflow", fake_cancel)

        sw = seekr_run.StageWorkflow(
            model=SimpleNamespace(directory="/work"),
            seekrflow=SimpleNamespace(name="tryp_ben"),
            stage=SimpleNamespace(name="bd", index=1),
            workflow_engine=object(),
            resource_name="cluster",
            resource=structures.Resource_remote_slurm(name="cluster"),
        )
        sw.job_ids = {"12345"}
        sw.state = "completed"

        sw.kill()

        assert cancel_calls == []

    def test_fused_member_kill_skips_cancel(self, monkeypatch):
        _ensure_radical_mock()
        from seekrflow.modules import seekr_run

        cancel_calls: list[tuple[str | None, str | None]] = []
        status_calls = {"n": 0}

        def fake_status(*args, **kwargs):
            status_calls["n"] += 1
            return {
                "stage_status": {"state": "started", "progress": 0.5},
                "manager_status": {
                    "jobs": [{"JobID": "999", "State": "RUNNING"}],
                },
            }

        def fake_cancel(
                seekrflow, resource_name, job_id=None, job_name=None, silent=False):
            cancel_calls.append((job_id, job_name))

        monkeypatch.setattr(workload_remote, "status_remote", fake_status)
        monkeypatch.setattr(
            workload_remote, "submit_remote_cancel_workflow", fake_cancel)

        sw = seekr_run.StageWorkflow(
            model=SimpleNamespace(directory="/work"),
            seekrflow=SimpleNamespace(name="tryp_ben"),
            stage=SimpleNamespace(name="logistic", index=2),
            workflow_engine=object(),
            resource_name="cluster",
            resource=structures.Resource_remote_slurm(name="cluster"),
        )
        sw.fusion_host = "metadynamics_seed"
        sw.state = "started"

        sw.kill()

        assert cancel_calls == []
        assert status_calls["n"] == 0


def test_detach_shutdown_does_not_kill(monkeypatch):
    _ensure_radical_mock()
    from seekrflow.modules import seekr_run

    killed: list[str] = []

    class FakeStageWorkflow:
        stage = SimpleNamespace(name="bd")

        def kill(self):
            killed.append(self.stage.name)

    class FakeEngine:
        async def shutdown(self):
            return None

    pipeline = seekr_run.SeekrPipeline.__new__(seekr_run.SeekrPipeline)
    pipeline._shutting_down = False
    pipeline._stop_event = None
    pipeline.stage_workflows = [FakeStageWorkflow()]
    pipeline.workflow_engine = FakeEngine()
    monkeypatch.setattr(
        seekr_run.SeekrPipeline, "write_status_snapshot", lambda self: None)

    asyncio.run(pipeline.detach_shutdown())
    assert killed == []


def test_probe_remote_launch_completed(monkeypatch):
    _ensure_radical_mock()
    from seekrflow.modules import seekr_run

    def fake_status(*args, **kwargs):
        return {
            "stage_status": {"state": "completed", "progress": 1.0},
            "manager_status": {"jobs": []},
        }

    monkeypatch.setattr(workload_remote, "status_remote", fake_status)

    sw = seekr_run.StageWorkflow(
        model=SimpleNamespace(directory="/work"),
        seekrflow=SimpleNamespace(name="run1"),
        stage=SimpleNamespace(name="bd", index=1),
        workflow_engine=object(),
        resource_name="cluster",
        resource=structures.Resource_remote_slurm(name="cluster"),
    )
    action, status = asyncio.run(sw.probe_remote_launch())
    assert action == "completed"
    assert status is not None


def test_probe_remote_launch_failure_falls_through(monkeypatch):
    _ensure_radical_mock()
    from seekrflow.modules import seekr_run

    def fake_status(*args, **kwargs):
        raise RuntimeError("ssh down")

    monkeypatch.setattr(workload_remote, "status_remote", fake_status)

    sw = seekr_run.StageWorkflow(
        model=SimpleNamespace(directory="/work"),
        seekrflow=SimpleNamespace(name="run1"),
        stage=SimpleNamespace(name="bd", index=1),
        workflow_engine=object(),
        resource_name="cluster",
        resource=structures.Resource_remote_slurm(name="cluster"),
    )
    action, status = asyncio.run(sw.probe_remote_launch())
    assert action == "submit"
    assert status is None
