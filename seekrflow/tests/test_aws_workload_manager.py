"""
Unit tests for AWS Batch workload manager helpers (mocked boto3).
"""
from __future__ import annotations

import json
import os
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from seekrflow.modules import structures
from seekrflow.modules.workload_managers import aws as workload_aws


def _aws_resource(**kwargs):
    defaults = dict(
        name="aws_gpu",
        region="us-west-2",
        account_id="123456789012",
        compute_env_name="ce",
        job_queue_name="queue",
        job_timeout_seconds=7200,
        n_gpus=1,
        n_vcpus=4,
        memory_mb=14000,
        transfer_settings=structures.Transfer_settings_aws_s3(
            bucket="my-bucket", prefix="runs"),
    )
    defaults.update(kwargs)
    return structures.Resource_cloud_aws(**defaults)


class TestAwsNamingAndUris:
    def test_job_definition_name(self):
        resource = _aws_resource(name="My GPU!")
        assert workload_aws.job_definition_name_for_resource(resource) == (
            "seekrflow-My-GPU-")

    def test_s3_model_uri_includes_seekrflow_name(self):
        resource = _aws_resource()
        assert workload_aws.s3_model_uri(resource, "hg") == (
            "s3://my-bucket/runs/hg")

    def test_s3_status_key(self):
        resource = _aws_resource()
        assert workload_aws.s3_status_key(resource, "hg", "equil") == (
            "runs/hg/.aws_runner/equil_status.json")


class TestBuildContainerCommand:
    def test_full_script_includes_fused_stages_and_status_loop(self):
        resource = _aws_resource()
        specs = [
            {"name": "seed", "force_overwrite": False, "benchmark_mode": False},
            {
                "name": "logistic",
                "force_overwrite": True,
                "benchmark_mode": False,
            },
        ]
        script = workload_aws.build_aws_container_script(
            resource, "hg", specs)
        assert "seed_status.json" in script
        assert "logistic_status.json" in script
        assert "seekr_run.run(model, 'seed'" in script
        assert "seekr_run.run(model, 'logistic'" in script
        assert "True, None, benchmark=False)" in script  # logistic force
        assert "STATUS_PID" in script
        assert "_progress_from_stage_progress" in script
        assert "Checking seekr status instructions" in script
        assert "AVAILABLE_INSTRUCTIONS" in script
        assert "seekrflow_stages_started" in script
        assert "waiting for prior stage(s) in this Batch job" in script
        assert "current_step" not in script
        assert "maximum_steps" not in script
        assert "total_steps" not in script
        assert "aws s3 sync /work s3://my-bucket/runs/hg" in script
        assert "aws s3 sync s3://my-bucket/runs/hg /work" in script
        assert f"--exclude {resource.transfer_settings.get_input_tarball_name()}" in script \
            or "--exclude model_input.tar.gz" in script
        assert "--exclude '.multiprocessing/*'" in script
        assert "--exclude '.aws_runner/*'" in script
        assert "_tar_rc" in script
        assert "dispatch.json" in script
        assert "_dispatch.json" in script
        assert '"version": 1' in script or "'version': 1" in script
        # Each fused stage is marked started before its seekr command.
        assert "echo 'seed' >> /tmp/seekrflow_stages_started" in script
        assert "echo 'logistic' >> /tmp/seekrflow_stages_started" in script

    def test_container_overrides_stub_is_short(self):
        resource = _aws_resource()
        cmd = workload_aws.build_aws_container_command(
            resource, "hg", "equil_stage_1")
        assert cmd[0] == "bash"
        stub = cmd[2]
        assert "seekrflow_run.sh" in stub
        assert "equil_stage_1_run.sh" in stub
        assert "seekr_run.run" not in stub
        assert len(json.dumps(cmd)) < 8192

    def test_script_embeds_custom_seekr_command_and_array_env(self):
        resource = _aws_resource()
        specs = [{
            "name": "MMVT",
            "force_overwrite": False,
            "benchmark_mode": False,
            "seekr_command": (
                "python -c \"print('array-driver')\" > MMVT_run.out"
            ),
        }]
        script = workload_aws.build_aws_container_script(
            resource, "hg", specs)
        assert "AWS_BATCH_JOB_ARRAY_INDEX" in script
        assert "array-driver" in script
        assert "seekr_run.run(model, 'MMVT'" not in script


class TestSubmitStatusCancel:
    def test_submit_registers_and_submits(self, monkeypatch, tmp_path):
        resource = _aws_resource()
        batch = MagicMock()
        batch.register_job_definition.return_value = {
            "jobDefinitionArn": "arn:jobdef"}
        batch.submit_job.return_value = {"jobId": "job-1"}
        s3 = MagicMock()
        monkeypatch.setattr(workload_aws, "_batch_client", lambda r: batch)
        monkeypatch.setattr(workload_aws, "_s3_client", lambda r: s3)

        seekrflow = SimpleNamespace(name="hg")
        result = workload_aws.submit_aws_job(
            seekrflow,
            "seed",
            resource,
            [{"name": "seed", "force_overwrite": False,
              "benchmark_mode": False}],
            str(tmp_path),
            timeout_seconds=3600,
            cpus=8,
            memory_mb=16000,
        )
        assert result["success"] is True
        assert result["job_id"] == "job-1"
        s3.put_object.assert_called_once()
        put_kw = s3.put_object.call_args.kwargs
        assert put_kw["Key"].endswith("seed_run.sh")
        assert b"seekr_run.run" in put_kw["Body"]
        batch.register_job_definition.assert_called_once()
        reqs = batch.register_job_definition.call_args.kwargs[
            "containerProperties"]["resourceRequirements"]
        by_type = {r["type"]: r["value"] for r in reqs}
        assert by_type["VCPU"] == "8"
        assert by_type["MEMORY"] == "16000"
        command = batch.submit_job.call_args.kwargs[
            "containerOverrides"]["command"]
        assert "seed_run.sh" in command[2]
        assert "seekr_run.run" not in command[2]
        assert "arrayProperties" not in batch.submit_job.call_args.kwargs
        assert batch.submit_job.call_args.kwargs["timeout"] == {
            "attemptDurationSeconds": 3600}
        state = workload_aws.load_job_state(str(tmp_path), "seed")
        assert state["job_id"] == "job-1"

    def test_submit_array_job(self, monkeypatch, tmp_path):
        resource = _aws_resource()
        batch = MagicMock()
        batch.register_job_definition.return_value = {
            "jobDefinitionArn": "arn:jobdef"}
        batch.submit_job.return_value = {"jobId": "parent-1"}
        s3 = MagicMock()
        monkeypatch.setattr(workload_aws, "_batch_client", lambda r: batch)
        monkeypatch.setattr(workload_aws, "_s3_client", lambda r: s3)
        seekrflow = SimpleNamespace(name="hg")
        result = workload_aws.submit_aws_job(
            seekrflow,
            "MMVT",
            resource,
            [{"name": "MMVT", "force_overwrite": False,
              "benchmark_mode": False,
              "seekr_command": "python -c 'pass'"}],
            str(tmp_path),
            array_size=4,
        )
        assert result["success"] is True
        assert result["array_size"] == 4
        assert batch.submit_job.call_args.kwargs["arrayProperties"] == {
            "size": 4}
        state = workload_aws.load_job_state(str(tmp_path), "MMVT")
        assert state["array_size"] == 4

    def test_status_merges_batch_and_s3(self, monkeypatch, tmp_path):
        resource = _aws_resource()
        workload_aws.save_job_state(str(tmp_path), "seed", "job-1", "hg_seed")
        batch = MagicMock()
        batch.describe_jobs.return_value = {
            "jobs": [{
                "jobId": "job-1",
                "jobName": "hg_seed",
                "status": "RUNNING",
                "statusReason": "",
            }],
        }
        monkeypatch.setattr(workload_aws, "_batch_client", lambda r: batch)
        monkeypatch.setattr(
            workload_aws,
            "_download_stage_status_from_s3",
            lambda *a, **k: {
                "success": True,
                "stage_status": {
                    "finished": False,
                    "state": "started",
                    "progress": 0.4,
                    "notes": "",
                },
            },
        )
        seekrflow = SimpleNamespace(name="hg")
        status = workload_aws.status_aws(
            seekrflow, "seed", resource, str(tmp_path))
        assert status["success"] is True
        assert status["stage_status"]["progress"] == 0.4
        assert status["manager_status"]["jobs"][0]["JobID"] == "job-1"
        assert status["manager_status"]["jobs"][0]["State"] == "RUNNING"

    def test_status_failed_includes_batch_reason_and_logs(
            self, monkeypatch, tmp_path):
        resource = _aws_resource()
        workload_aws.save_job_state(str(tmp_path), "seed", "job-9", "hg_seed")
        batch = MagicMock()
        batch.describe_jobs.return_value = {
            "jobs": [{
                "jobId": "job-9",
                "jobName": "hg_seed",
                "status": "FAILED",
                "statusReason": "Essential container in task exited",
                "container": {
                    "exitCode": 1,
                    "reason": "Essential container in task exited",
                    "logStreamName": "seekrflow/default/abc",
                },
            }],
        }
        monkeypatch.setattr(workload_aws, "_batch_client", lambda r: batch)
        monkeypatch.setattr(
            workload_aws,
            "_download_stage_status_from_s3",
            lambda *a, **k: None,
        )
        monkeypatch.setattr(
            workload_aws,
            "fetch_cloudwatch_log_tail",
            lambda *a, **k: "=== [2/5] Downloading\nAccessDenied",
        )
        seekrflow = SimpleNamespace(name="hg")
        status = workload_aws.status_aws(
            seekrflow, "seed", resource, str(tmp_path))
        assert status["success"] is False
        assert "job-9" in status["error"]
        assert "Essential container" in status["error"]
        assert "AccessDenied" in status["stage_status"]["notes"]
        failure_log = workload_aws.local_failure_log_path(str(tmp_path), "seed")
        assert os.path.isfile(failure_log)
        with open(failure_log) as f:
            body = f.read()
        assert "AccessDenied" in body
        assert "job-9" in body

    def test_status_array_lists_active_children(self, monkeypatch, tmp_path):
        resource = _aws_resource()
        workload_aws.save_job_state(
            str(tmp_path), "MMVT", "parent-1", "hg_MMVT", array_size=2)
        batch = MagicMock()

        def describe_jobs(jobs):
            if jobs == ["parent-1"]:
                return {"jobs": [{
                    "jobId": "parent-1",
                    "jobName": "hg_MMVT",
                    "status": "RUNNING",
                }]}
            out = []
            for jid in jobs:
                idx = jid.split(":")[-1]
                out.append({
                    "jobId": jid,
                    "jobName": "hg_MMVT",
                    "status": "RUNNING" if idx == "0" else "SUCCEEDED",
                })
            return {"jobs": out}

        batch.describe_jobs.side_effect = describe_jobs
        monkeypatch.setattr(workload_aws, "_batch_client", lambda r: batch)
        monkeypatch.setattr(
            workload_aws,
            "_download_stage_status_from_s3",
            lambda *a, **k: {
                "success": True,
                "stage_status": {
                    "finished": False,
                    "state": "started",
                    "progress": 0.5,
                    "notes": "",
                },
            },
        )
        seekrflow = SimpleNamespace(name="hg")
        status = workload_aws.status_aws(
            seekrflow, "MMVT", resource, str(tmp_path))
        assert status["success"] is True
        active_ids = [j["JobID"] for j in status["manager_status"]["jobs"]]
        assert active_ids == ["parent-1:0"]

    def test_cancel_terminates_by_id(self, monkeypatch):
        resource = _aws_resource()
        batch = MagicMock()
        monkeypatch.setattr(workload_aws, "_batch_client", lambda r: batch)
        result = workload_aws.cancel_aws_job(resource, job_id="job-9")
        assert result["success"] is True
        batch.terminate_job.assert_called_once_with(
            jobId="job-9", reason="seekrflow cancel")

    def test_cancel_and_reset_clears_s3_monitor_artifacts(
            self, monkeypatch, tmp_path):
        resource = _aws_resource()
        batch = MagicMock()
        s3 = MagicMock()
        monkeypatch.setattr(workload_aws, "_batch_client", lambda r: batch)
        monkeypatch.setattr(workload_aws, "_s3_client", lambda r: s3)
        seekrflow = SimpleNamespace(name="hg")
        workload_aws.save_job_state(
            str(tmp_path), "seed", "job-1", "hg_seed")
        result = workload_aws.cancel_and_reset_aws_stage(
            seekrflow,
            "seed",
            resource,
            str(tmp_path),
            monitor_stage_names=["seed", "assign"],
        )
        assert result["canceled_job"] == "job-1"
        batch.terminate_job.assert_called()
        deleted_keys = [
            c.kwargs["Key"] for c in s3.delete_object.call_args_list]
        assert any(k.endswith("seed_status.json") for k in deleted_keys)
        assert any(k.endswith("seed_dispatch.json") for k in deleted_keys)
        assert any(k.endswith("assign_status.json") for k in deleted_keys)
        assert any(k.endswith("assign_dispatch.json") for k in deleted_keys)
        assert workload_aws.load_job_state(str(tmp_path), "seed") is None

    def test_publish_stage_completed_status(self, monkeypatch):
        resource = _aws_resource()
        s3 = MagicMock()
        monkeypatch.setattr(workload_aws, "_s3_client", lambda r: s3)
        workload_aws.publish_stage_completed_status(
            resource, "hg", "MMVT", notes="skip test")
        s3.put_object.assert_called_once()
        kwargs = s3.put_object.call_args.kwargs
        assert kwargs["Bucket"] == "my-bucket"
        assert kwargs["Key"].endswith("MMVT_status.json")
        body = json.loads(kwargs["Body"].decode("utf-8"))
        assert body["success"] is True
        assert body["stage_status"]["state"] == "completed"
        assert body["stage_status"]["progress"] == 1.0
        assert body["stage_status"]["notes"] == "skip test"

    def test_fetch_unit_counts_cloud_from_dispatch(self, monkeypatch):
        resource = _aws_resource()

        def _no_local(*a, **k):
            raise RuntimeError("force S3 fallback")

        monkeypatch.setattr(
            "seekrflow.modules.workload_managers.dispatch_lowering"
            ".fetch_unit_counts_local",
            _no_local,
        )
        monkeypatch.setattr(
            workload_aws,
            "download_stage_dispatch_from_s3",
            lambda *a, **k: {
                "version": 1,
                "info": {
                    "2": {
                        "stage_name": "metadynamics_assign_structures",
                        "scope": "partitioned",
                        "num_anchors": 8,
                        "num_swarms": 1,
                    },
                    "3": {
                        "stage_name": "MMVT",
                        "scope": "partitioned",
                        "num_anchors": 8,
                        "num_swarms": 1,
                    },
                },
                "progress": {"finished": True, "progress": {}},
            },
        )
        monkeypatch.setattr(
            workload_aws,
            "_download_stage_status_from_s3",
            lambda *a, **k: None,
        )
        model = SimpleNamespace(stages=[
            SimpleNamespace(name="metadynamics_assign_structures", index=2),
            SimpleNamespace(name="MMVT", index=3),
        ])
        launching = SimpleNamespace(
            name="MMVT", index=3, input_stage_index=2)
        counts = workload_aws.fetch_unit_counts_cloud(
            resource, "hg", model, launching)
        assert counts.num_anchors == 8
        assert counts.num_swarms == 1
        assert counts.scope == "partitioned"

    def test_fetch_unit_counts_cloud_ramd_uses_launching_swarms(self, monkeypatch):
        resource = _aws_resource()

        def _no_local(*a, **k):
            raise RuntimeError("force S3 fallback")

        monkeypatch.setattr(
            "seekrflow.modules.workload_managers.dispatch_lowering"
            ".fetch_unit_counts_local",
            _no_local,
        )
        monkeypatch.setattr(
            workload_aws,
            "download_stage_dispatch_from_s3",
            lambda *a, **k: {
                "version": 1,
                "info": {
                    "2": {
                        "stage_name": "ramd_procedure_assign_swarm",
                        "scope": "unpartitioned",
                        "num_anchors": None,
                        "num_swarms": 1,
                    },
                    "3": {
                        "stage_name": "ramd_procedure_ramd",
                        "scope": "unpartitioned",
                        "num_anchors": None,
                        "num_swarms": 100,
                    },
                },
            },
        )
        model = SimpleNamespace(stages=[
            SimpleNamespace(name="ramd_procedure_assign_swarm", index=2),
            SimpleNamespace(name="ramd_procedure_ramd", index=3),
        ])
        launching = SimpleNamespace(
            name="ramd_procedure_ramd", index=3, input_stage_index=2)
        counts = workload_aws.fetch_unit_counts_cloud(
            resource, "hg", model, launching)
        assert counts.num_swarms == 100
        assert counts.scope == "unpartitioned"

    def test_fetch_unit_counts_cloud_missing_raises(self, monkeypatch):
        resource = _aws_resource()

        def _no_local(*a, **k):
            raise RuntimeError("force S3 fallback")

        monkeypatch.setattr(
            "seekrflow.modules.workload_managers.dispatch_lowering"
            ".fetch_unit_counts_local",
            _no_local,
        )
        monkeypatch.setattr(
            workload_aws,
            "download_stage_dispatch_from_s3",
            lambda *a, **k: None,
        )
        monkeypatch.setattr(
            workload_aws,
            "_download_stage_status_from_s3",
            lambda *a, **k: None,
        )
        model = SimpleNamespace(stages=[
            SimpleNamespace(name="assign", index=1),
        ])
        launching = SimpleNamespace(name="MMVT", index=2, input_stage_index=1)
        with pytest.raises(RuntimeError, match="Missing S3 dispatch info"):
            workload_aws.fetch_unit_counts_cloud(
                resource, "hg", model, launching)

    def test_fetch_stage_progress_cloud(self, monkeypatch):
        resource = _aws_resource()
        monkeypatch.setattr(
            workload_aws,
            "download_stage_dispatch_from_s3",
            lambda *a, **k: {
                "version": 1,
                "info": {},
                "progress": {
                    "finished": False,
                    "progress": {
                        0: {"attained": True, "progress": 1.0},
                        1: {"attained": False, "progress": 0.2},
                    },
                },
            },
        )
        finished, pmap = workload_aws.fetch_stage_progress_cloud(
            resource, "hg", "MMVT")
        assert finished is False
        assert pmap[0]["attained"] is True

    def test_no_job_def_generic_on_resource(self):
        attr_names = [
            a.name for a in structures.Resource_cloud_aws.__attrs_attrs__]
        assert "job_def_generic" not in attr_names
        resource = _aws_resource()
        assert not hasattr(resource, "job_def_generic")
