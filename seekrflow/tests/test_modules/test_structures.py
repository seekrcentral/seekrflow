"""
Tests for seekrflow/modules/structures.py
"""

import json

import pytest

import seekrflow.modules.workflows.stage_procedures as stage_procedures
from seekrflow.modules import structures


def test_local_shell_requires_transfer_none():
    with pytest.raises(ValueError, match="must be 'none'"):
        structures.Resource_remote_slurm(
            name="login_slurm",
            remote_interface=structures.Remote_interface_local_shell(),
            transfer_settings=structures.Transfer_settings_globus(),
        )
    with pytest.raises(ValueError, match="must be 'none'"):
        structures.Resource_remote_pbs(
            name="login_pbs",
            remote_interface=structures.Remote_interface_local_shell(),
            transfer_settings=structures.Transfer_settings_rsync(),
        )
    with pytest.raises(ValueError, match="must be 'none'"):
        structures.Resource_direct(
            name="local",
            remote_interface=structures.Remote_interface_local_shell(),
            transfer_settings=structures.Transfer_settings_globus(),
        )


def test_local_shell_resource_round_trip(tmp_path):
    executable = "/opt/seekr/bin/python"
    resource = structures.Resource_remote_slurm(
        name="login_slurm",
        remote_working_directory="/scratch/user/seekrflow_work",
        remote_interface=structures.Remote_interface_local_shell(
            python_executable=executable),
        transfer_settings=structures.Transfer_settings_none(),
    )
    flow = structures.Seekrflow(name="round_trip")
    flow.run_settings.resources = [resource]
    path = tmp_path / "seekrflow.json"
    path.write_text(json.dumps(flow.to_dict()))
    loaded = structures.load_seekrflow(str(path))
    loaded_resource = next(
        item for item in loaded.run_settings.resources
        if item.name == "login_slurm")
    assert isinstance(
        loaded_resource.remote_interface,
        structures.Remote_interface_local_shell)
    assert loaded_resource.remote_interface.python_executable == executable
    assert isinstance(
        loaded_resource.transfer_settings,
        structures.Transfer_settings_none)


def test_set_placement_resource_upserts_target():
    run_settings = structures.Run_settings(
        resources=[structures.Resource_remote_slurm(name="delta")])
    run_settings.set_placement_resource(
        ["metadynamics", "logistic"], "local")
    assert len(run_settings.placements) == 1
    assert run_settings.placements[0].target == ["metadynamics", "logistic"]
    assert run_settings.placements[0].resource == "local"
    run_settings.set_placement_resource(
        ["metadynamics", "logistic"], "delta")
    assert len(run_settings.placements) == 1
    assert run_settings.placements[0].resource == "delta"
    with pytest.raises(ValueError, match="not found"):
        run_settings.set_placement_resource(["bd"], "missing")


def test_placement_dispatch_and_time_policy_round_trip(tmp_path):
    flow = structures.Seekrflow(name="placement_round_trip")
    flow.workflow.procedure = stage_procedures.MMVT_stage_procedure(name="MMVT")
    flow.run_settings.placements = [
        structures.Placement(
            target=["MMVT"],
            resource="delta",
            dispatch=stage_procedures.Dispatch(
                dimensions=["anchor"], group_size=2, concurrency=1),
            time_policy=structures.Time_policy_fixed(time_limit="01:00:00"),
        ),
    ]
    path = tmp_path / "seekrflow.json"
    path.write_text(json.dumps(flow.to_dict()))
    loaded = structures.load_seekrflow(str(path))
    placement = next(
        item for item in loaded.run_settings.placements
        if item.target == ["MMVT"])
    assert placement.dispatch.group_size == 2
    assert isinstance(placement.time_policy, structures.Time_policy_fixed)
    assert placement.time_policy.time_limit == "01:00:00"
