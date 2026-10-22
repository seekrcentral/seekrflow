"""
Tests for seekrflow/modules/workload_managers/aws.py
"""

from __future__ import annotations

import json
import pathlib
import sys
import types
from unittest.mock import MagicMock

import pytest

import seekrflow.modules.job.structures as job_structures
from seekrflow.modules.workload_managers import aws as workload_aws

_STRUCTURES_SOURCE = pathlib.Path(job_structures.__file__).read_text()
_IMAGE = "123456789012.dkr.ecr.us-west-2.amazonaws.com/seekr-engines-bd:latest"


def _resource_payload(**overrides) -> dict:
    payload = {
        "name": "aws_gpu",
        "region": "us-west-2",
        "transfer_settings_s3_uri": "s3://my-bucket/runs",
        "seekr_image_uri": _IMAGE,
        "transfer_settings_input_tarball_name": "model_input.tar.gz",
        "transfer_settings_bucket": "my-bucket",
        "job_name": "hg_seed",
        "cpus": 8,
        "memory_mb": 16000,
        "n_gpus": 1,
        "time_limit": 3600,
        "job_queue_name": "queue",
    }
    payload.update(overrides)
    return payload


def _job_spec(array_index: int, internal_id: int = 0, interval: float = 15.0) -> dict:
    spec = job_structures.JobSpec(
        internal_id=internal_id,
        remote_root_dir="/work",
        status_write_interval=interval,
        telemetry_poll_interval=None,
        stage_specs=[job_structures.StageSpec(
            stage_index=0,
            stage_name="seed",
            force_overwrite=False,
            benchmark=False,
        )],
        run_units=[job_structures.RunUnit(anchor=0, swarm_id=None)],
        concurrency=1,
        array_index=array_index,
    )
    return spec.to_dict()


def _manager_payload(n_specs: int = 1, **resource_overrides) -> dict:
    return {
        "seekrflow_name": "hg",
        "resource_payload": _resource_payload(**resource_overrides),
        "job_specs": [_job_spec(index) for index in range(n_specs)],
        "runner_files": {
            "job.py": "print('job')\n",
            "structures.py": _STRUCTURES_SOURCE,
        },
    }


def _status_payload(
        job_id: str = "job-1",
        array_indices: list[int] | None = None,
        stage_indices: list[int] | None = None,
        ) -> dict:
    return {
        "resource_payload": _resource_payload(),
        "system_payloads": {
            "hg": {
                "jobs": [{
                    "internal_id": 0,
                    "job_id": job_id,
                    "job_name": "hg_seed",
                    "stage_indices": [1] if stage_indices is None else stage_indices,
                    "array_indices": [0] if array_indices is None else array_indices,
                }],
            },
        },
    }


def _cancel_payload(jobs: list[dict], remove_json_files: bool = False) -> dict:
    return {
        "remove_json_files": remove_json_files,
        "resource_payload": _resource_payload(),
        "system_payloads": {"hg": {"jobs": jobs}},
    }


class _NoSuchKey(Exception):
    def __init__(self) -> None:
        super().__init__("missing")
        self.response = {"Error": {"Code": "NoSuchKey"}}


def _s3_body(payload: dict) -> dict:
    body = MagicMock()
    body.read.return_value = json.dumps(payload).encode("utf-8")
    return {"Body": body}


def _clients(
        monkeypatch,
        pages: list | None = None,
        job_definitions: dict | None = None,
        ):
    s3 = MagicMock()
    paginator = MagicMock()
    paginator.paginate.return_value = pages if pages is not None else [{"Contents": []}]
    s3.get_paginator.return_value = paginator
    batch = MagicMock()
    batch.describe_job_definitions.return_value = (
        job_definitions if job_definitions is not None else {"jobDefinitions": []})
    batch.register_job_definition.return_value = {"jobDefinitionArn": "arn:jobdef"}
    batch.submit_job.return_value = {"jobId": "job-1"}
    monkeypatch.setattr(workload_aws, "_s3_client", lambda resource: s3)
    monkeypatch.setattr(workload_aws, "_batch_client", lambda resource: batch)
    return s3, batch


def test_job_definition_name_sanitizes_resource_name():
    assert workload_aws.job_definition_name_for_resource(
        {"name": "My GPU!"}) == "seekrflow-My-GPU-"


def test_container_script_syncs_model_and_runs_job_py():
    single = workload_aws.build_aws_container_script(_manager_payload(1))
    assert "s3://my-bucket/runs/hg" in single
    assert "model_input.tar.gz" in single
    assert "_tar_rc" in single
    assert "--exclude model_input.tar.gz" in single
    assert "--exclude '.seekr_runner/*'" in single
    assert "--exclude '.batch_runner/*'" in single
    assert "python /work/.seekr_jobs/job.py /work/.seekr_jobs/job_spec_0_0.json" in single
    assert "aws s3 sync /work/.stage_states s3://my-bucket/runs/hg/.stage_states" in single
    assert "sleep 15.0" in single
    assert "AWS_BATCH_JOB_ARRAY_INDEX" not in single
    assert "seekr_run.run" not in single
    assert "aws s3 sync s3://my-bucket/runs/hg /work" in single
    assert "aws s3 sync /work s3://my-bucket/runs/hg" in single

    array_script = workload_aws.build_aws_container_script(_manager_payload(2))
    assert "job_spec_0_${AWS_BATCH_JOB_ARRAY_INDEX}.json" in array_script


def test_container_command_fetches_script_and_rejects_oversize(monkeypatch):
    command = workload_aws.build_aws_container_command(
        _resource_payload(), "hg", "hg_seed")
    assert command[0] == "bash"
    assert command[1] == "-lc"
    stub = command[2]
    assert "aws s3 cp s3://my-bucket/runs/hg/.aws_runner/hg_seed_run.sh /tmp/seekrflow_run.sh" in stub
    assert "exec bash /tmp/seekrflow_run.sh" in stub
    assert "seekr_run.run" not in stub

    monkeypatch.setattr(
        workload_aws, "s3_run_script_uri",
        lambda *args, **kwargs: "s3://bucket/" + ("x" * 9000))
    with pytest.raises(ValueError, match="8192"):
        workload_aws.build_aws_container_command(
            _resource_payload(), "hg", "hg_seed")


def test_next_internal_id_increments_past_numeric_info_files():
    s3 = MagicMock()
    paginator = MagicMock()
    paginator.paginate.side_effect = [
        [{"Contents": [
            {"Key": "runs/hg/.batch_runner/batch_info_2.json"},
            {"Key": "runs/hg/.batch_runner/batch_info_7.json"},
            {"Key": "runs/hg/.batch_runner/batch_info_notes.json"},
        ]}],
        [{}],
    ]
    s3.get_paginator.return_value = paginator
    assert workload_aws.next_internal_id(s3, "my-bucket", "runs/hg") == 8
    assert workload_aws.next_internal_id(s3, "my-bucket", "runs/hg") == 0


def test_map_batch_state():
    assert workload_aws.map_batch_state("RUNNING") == "running"
    assert workload_aws.map_batch_state("STARTING") == "running"
    assert workload_aws.map_batch_state("SUBMITTED") == "queued"
    assert workload_aws.map_batch_state("PENDING") == "queued"
    assert workload_aws.map_batch_state("RUNNABLE") == "queued"
    assert workload_aws.map_batch_state("FAILED", "user cancel") == "cancelled"
    assert workload_aws.map_batch_state("FAILED", "Essential container") == "failed"
    assert workload_aws.map_batch_state("SUCCEEDED") == "idle"


def _requirements() -> list[dict]:
    return [
        {"type": "GPU", "value": "1"},
        {"type": "VCPU", "value": "8"},
        {"type": "MEMORY", "value": "16000"},
    ]


def test_submit_registers_single_job_and_uploads_runner_files(monkeypatch):
    s3, batch = _clients(monkeypatch)
    result = workload_aws.submit_aws_job(_manager_payload(1))
    assert result["success"] is True
    assert result["job_id"] == "job-1"
    keys = [call.kwargs["Key"] for call in s3.put_object.call_args_list]
    assert "runs/hg/.batch_runner/batch_info_0.json" in keys
    assert "runs/hg/.seekr_jobs/job.py" in keys
    assert "runs/hg/.seekr_jobs/structures.py" in keys
    assert "runs/hg/.seekr_jobs/job_spec_0_0.json" in keys
    assert "runs/hg/.aws_runner/hg_seed_run.sh" in keys
    batch.register_job_definition.assert_called_once()
    requirements = batch.register_job_definition.call_args.kwargs[
        "containerProperties"]["resourceRequirements"]
    by_type = {item["type"]: item["value"] for item in requirements}
    assert by_type["VCPU"] == "8"
    assert by_type["MEMORY"] == "16000"
    submit_kwargs = batch.submit_job.call_args.kwargs
    assert "arrayProperties" not in submit_kwargs
    assert submit_kwargs["timeout"] == {"attemptDurationSeconds": 3600}
    assert "seekr_run.run" not in submit_kwargs["containerOverrides"]["command"][2]


def test_submit_reuses_matching_job_definition(monkeypatch):
    s3, batch = _clients(monkeypatch, job_definitions={
        "jobDefinitions": [{
            "jobDefinitionArn": "arn:existing",
            "containerProperties": {
                "image": _IMAGE,
                "resourceRequirements": _requirements(),
            },
        }],
    })
    result = workload_aws.submit_aws_job(_manager_payload(1))
    assert result["success"] is True
    batch.register_job_definition.assert_not_called()
    assert batch.submit_job.call_args.kwargs["jobDefinition"] == "arn:existing"
    del s3


def test_submit_array_sets_size_from_job_specs(monkeypatch):
    s3, batch = _clients(monkeypatch)
    result = workload_aws.submit_aws_job(_manager_payload(2))
    assert result["success"] is True
    assert batch.submit_job.call_args.kwargs["arrayProperties"] == {"size": 2}
    del s3


def test_status_merges_running_parent_with_stage_state(monkeypatch):
    s3, batch = _clients(monkeypatch)
    batch.describe_jobs.return_value = {"jobs": [{
        "jobId": "job-1",
        "status": "RUNNING",
        "statusReason": "",
    }]}
    s3.get_object.return_value = _s3_body({"state": "started", "finished": False})
    result = workload_aws.status_aws(_status_payload())
    assert result["success"] is True
    batch.describe_jobs.assert_called_with(jobs=["job-1"])
    job = result["payload"]["hg"]["job-1"]
    assert job["manager_dicts_by_array_index"][0]["state"] == "running"
    assert job["stage_dicts_by_stage_index"][1]["state"] == "started"


def test_status_missing_stage_state_is_unstarted(monkeypatch):
    s3, batch = _clients(monkeypatch)
    batch.describe_jobs.return_value = {"jobs": [{
        "jobId": "job-1",
        "status": "RUNNING",
    }]}
    s3.get_object.side_effect = _NoSuchKey()
    result = workload_aws.status_aws(_status_payload())
    stage = result["payload"]["hg"]["job-1"]["stage_dicts_by_stage_index"][1]
    assert stage["state"] == "unstarted"
    assert stage["finished"] is False


def test_status_failed_notes_include_cloudwatch_tail(monkeypatch):
    s3, batch = _clients(monkeypatch)
    batch.describe_jobs.return_value = {"jobs": [{
        "jobId": "job-9",
        "status": "FAILED",
        "statusReason": "Essential container in task exited",
        "container": {
            "exitCode": 1,
            "logStreamName": "seekrflow/default/abc",
        },
    }]}
    s3.get_object.side_effect = _NoSuchKey()
    logs = MagicMock()
    logs.get_log_events.return_value = {"events": [{"message": "AccessDenied"}]}
    fake_boto3 = types.ModuleType("boto3")
    fake_boto3.client = lambda *args, **kwargs: logs
    monkeypatch.setitem(sys.modules, "boto3", fake_boto3)
    result = workload_aws.status_aws(_status_payload(job_id="job-9"))
    notes = result["payload"]["hg"]["job-9"]["manager_dicts_by_array_index"][0]["notes"]
    assert "Essential container in task exited" in notes
    assert "AccessDenied" in notes


def test_status_array_reports_each_child(monkeypatch):
    s3, batch = _clients(monkeypatch)
    described = {
        "parent-1": {"jobId": "parent-1", "status": "RUNNING"},
        "parent-1:0": {"jobId": "parent-1:0", "status": "RUNNING"},
        "parent-1:1": {"jobId": "parent-1:1", "status": "SUCCEEDED"},
    }

    def describe_jobs(jobs):
        return {"jobs": [described[job_id] for job_id in jobs]}

    batch.describe_jobs.side_effect = describe_jobs
    s3.get_object.side_effect = _NoSuchKey()
    result = workload_aws.status_aws(_status_payload(
        job_id="parent-1", array_indices=[0, 1]))
    batch.describe_jobs.assert_called_with(
        jobs=["parent-1", "parent-1:0", "parent-1:1"])
    managers = result["payload"]["hg"]["parent-1"]["manager_dicts_by_array_index"]
    assert managers[0]["state"] == "running"
    assert managers[1]["state"] == "idle"


def test_cancel_terminates_by_id_and_rejects_unnamed_job(monkeypatch):
    s3, batch = _clients(monkeypatch)
    result = workload_aws.cancel_aws_job(_cancel_payload([{"job_id": "job-9"}]))
    assert result["success"] is True
    assert result["payload"] == ["job-9"]
    batch.terminate_job.assert_called_once_with(
        jobId="job-9", reason="seekrflow cancel")

    missing = workload_aws.cancel_aws_job(_cancel_payload([{}]))
    assert missing["success"] is False
    assert missing["error"] == "job_id or job_name required"
    del s3


def test_cancel_remove_json_files_deletes_runner_prefixes(monkeypatch):
    s3, batch = _clients(monkeypatch)
    prefixes = []

    def paginate(Bucket, Prefix):
        prefixes.append(Prefix)
        return [{"Contents": [{"Key": Prefix + "kept.json"}]}]

    s3.get_paginator.return_value.paginate.side_effect = paginate
    result = workload_aws.cancel_aws_job(_cancel_payload(
        [{"job_id": "job-9"}], remove_json_files=True))
    assert result["success"] is True
    assert "runs/hg/.batch_runner/" in prefixes
    assert "runs/hg/.stage_states/" in prefixes
    assert s3.delete_objects.called
    s3.delete_object.assert_not_called()
    del batch
