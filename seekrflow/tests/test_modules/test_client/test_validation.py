"""
Tests for seekrflow/modules/client/validation.py
"""

import pytest

import seekrflow.modules.client.validation as validation
import seekrflow.modules.structures as structures
import seekrflow.modules.workflows.stage_procedures as stage_procedures
from seekrflow.tests.test_modules.test_workflows.test_stage_procedures import (
    host_guest_procedure,
)


class _FakeStage:
    def __init__(self, name, input_stage_index=0):
        self.name = name
        self.input_stage_index = input_stage_index


class _FakeModel:
    def __init__(self, stages):
        self.stages = stages


def _resolve(run_settings, stage_name, procedure=None):
    return validation.resolve_stage_execution(
        run_settings, stage_name, procedure or host_guest_procedure())


def test_role_default_when_no_placement():
    resolved = _resolve(
        structures.Run_settings(), "metadynamics_assign_structures")
    assert isinstance(resolved.resource, structures.Resource_direct)
    assert resolved.resource.name == "local"
    assert resolved.co_schedule_with == "predecessor"


def test_broad_placement_applies_to_subroles():
    run_settings = structures.Run_settings(
        resources=[structures.Resource_remote_slurm(name="delta")],
        placements=[
            structures.Placement(target=["metadynamics"], resource="delta"),
        ])
    assert _resolve(run_settings, "metadynamics_seed").resource.name == "delta"
    assert _resolve(
        run_settings, "metadynamics_assign_structures").resource.name == "delta"


def test_specific_placement_overrides_broader_target():
    run_settings = structures.Run_settings(
        resources=[structures.Resource_remote_slurm(name="delta")],
        placements=[
            structures.Placement(
                target=["metadynamics", "logistic"], resource="local"),
            structures.Placement(target=["metadynamics"], resource="delta"),
        ])
    assert _resolve(run_settings, "metadynamics_seed").resource.name == "delta"
    assert _resolve(
        run_settings, "metadynamics_assign_structures").resource.name == "local"


def test_duplicate_placement_targets_raise():
    run_settings = structures.Run_settings(placements=[
        structures.Placement(target=["bd"], resource="local"),
        structures.Placement(target=["bd"], resource="delta"),
    ])
    with pytest.raises(ValueError, match="Duplicate placement target"):
        _resolve(run_settings, "bd")


def test_dispatch_merge_preserves_role_dimensions():
    ramd = stage_procedures.RAMD_stage_procedure(name="ramd_proc")
    procedure = stage_procedures.Composite_stage_procedure(
        name="root", procedures=[ramd])
    run_settings = structures.Run_settings(
        resources=[structures.Resource_remote_slurm(name="cluster")],
        placements=[
            structures.Placement(
                target=["ramd_proc"],
                resource="cluster",
                dispatch=stage_procedures.Dispatch(group_size=2)),
        ])
    resolved = _resolve(run_settings, "ramd_proc_ramd", procedure)
    assert resolved.dispatch.dimensions == ["swarm"]
    assert resolved.dispatch.group_size == 2


def test_local_collapses_group_size_and_keeps_dimensions():
    run_settings = structures.Run_settings(placements=[
        structures.Placement(
            target=["MMVT"],
            resource="local",
            dispatch=stage_procedures.Dispatch(group_size=4)),
    ])
    resolved = _resolve(run_settings, "MMVT")
    assert resolved.dispatch.group_size is None
    assert resolved.dispatch.dimensions == ["anchor"]


def test_concurrency_above_mps_raises():
    run_settings = structures.Run_settings(
        resources=[structures.Resource_remote_slurm(name="cluster", mps=2)],
        placements=[
            structures.Placement(
                target=["MMVT"],
                resource="cluster",
                dispatch=stage_procedures.Dispatch(concurrency=8)),
        ])
    with pytest.raises(ValueError, match="dispatch.concurrency=8"):
        _resolve(run_settings, "MMVT")


def test_time_limit_parsers():
    assert validation.time_limit_to_seconds("02:30:00") == 9000
    assert validation.seconds_to_time_limit(9000) == "02:30:00"
    assert validation.time_limit_to_seconds("1-00:00:00") == 86400
    assert validation.seconds_to_time_limit(86400) == "24:00:00"


def test_slurm_defaults_and_placement_cpu_memory_override():
    resource = structures.Resource_remote_slurm(
        name="cluster",
        cpus_per_task=16,
        memory_per_node=32000,
        mps=2,
    )
    run_settings = structures.Run_settings(
        resources=[resource],
        placements=[
            structures.Placement(target=["bd"], resource="cluster"),
            structures.Placement(
                target=["MMVT"],
                resource="cluster",
                cpus=8,
                memory_mb=16000,
            ),
        ])
    bd = _resolve(run_settings, "bd")
    assert bd.cpus == 16
    assert bd.memory_mb == 32000
    assert bd.mps == 2
    assert bd.time_limit is None
    assert isinstance(bd.time_policy, structures.Time_policy_adaptive)
    mmvt = _resolve(run_settings, "MMVT")
    assert mmvt.cpus == 8
    assert mmvt.memory_mb == 16000
    assert mmvt.mps == 2


def test_aws_mmvt_keeps_array_group_size():
    resource = structures.Resource_cloud_aws(
        name="aws_gpu",
        region="us-west-2",
        account_id="123",
        compute_env_name="ce",
        job_queue_name="q",
        n_vcpus=8,
        memory_mb=1000,
        mps=2,
    )
    run_settings = structures.Run_settings(
        resources=[resource],
        placements=[
            structures.Placement(target=["MMVT"], resource="aws_gpu"),
        ])
    resolved = _resolve(run_settings, "MMVT")
    assert resolved.cpus == 8
    assert resolved.memory_mb == 1000
    assert resolved.mps == 2
    assert resolved.dispatch.dimensions == ["anchor"]
    assert resolved.dispatch.group_size == 1


def _seekrflow_with(procedure, placements, resources=None):
    seekrflow = structures.Seekrflow()
    seekrflow.workflow.procedure = procedure
    seekrflow.run_settings.resources = list(resources or [])
    seekrflow.run_settings.placements = placements
    return seekrflow


def test_validate_unknown_target_and_resource():
    procedure = host_guest_procedure()
    with pytest.raises(ValueError, match="Unknown placement target"):
        validation.validate_run_settings(_seekrflow_with(
            procedure,
            [structures.Placement(target=["nonexistent"], resource="local")]))
    with pytest.raises(ValueError, match="not found"):
        validation.validate_run_settings(_seekrflow_with(
            procedure,
            [structures.Placement(target=["bd"], resource="missing")]))


def test_validate_co_schedule_resource_and_dimensions():
    procedure = host_guest_procedure()
    model = _FakeModel([
        _FakeStage("metadynamics_seed", input_stage_index=0),
        _FakeStage("metadynamics_assign_structures", input_stage_index=1),
    ])
    delta = structures.Resource_remote_slurm(name="delta")
    with pytest.raises(ValueError, match="same resource"):
        validation.validate_run_settings(
            _seekrflow_with(
                procedure,
                [
                    structures.Placement(
                        target=["metadynamics", "sampling"], resource="delta"),
                    structures.Placement(
                        target=["metadynamics", "logistic"], resource="local"),
                ],
                resources=[delta],
            ),
            model,
        )
    with pytest.raises(ValueError, match="dispatch.dimensions"):
        validation.validate_run_settings(
            _seekrflow_with(
                procedure,
                [
                    structures.Placement(
                        target=["metadynamics", "sampling"],
                        resource="delta",
                        dispatch=stage_procedures.Dispatch(
                            dimensions=["anchor"])),
                    structures.Placement(
                        target=["metadynamics", "logistic"], resource="delta"),
                ],
                resources=[delta],
            ),
            model,
        )


def test_validate_fixed_time_limit_against_resource_cap():
    procedure = host_guest_procedure()
    resource = structures.Resource_remote_slurm(
        name="cluster", max_time_limit="01:00:00")
    over = _seekrflow_with(
        procedure,
        [structures.Placement(
            target=["bd"],
            resource="cluster",
            time_policy=structures.Time_policy_fixed(time_limit="02:00:00"),
        )],
        resources=[resource],
    )
    with pytest.raises(ValueError, match="exceeds resource"):
        validation.validate_run_settings(over)
    within = _seekrflow_with(
        procedure,
        [structures.Placement(
            target=["bd"],
            resource="cluster",
            time_policy=structures.Time_policy_fixed(time_limit="00:30:00"),
        )],
        resources=[resource],
    )
    validation.validate_run_settings(within)
