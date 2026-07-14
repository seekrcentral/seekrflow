"""
Tests for placement/dispatch policy resolution (Phases 1-7).

Run with: pytest --noconftest seekrflow/tests/test_placement_dispatch.py
"""

import tempfile
import os

import pytest

import seekrflow.modules.structures as structures
import seekrflow.modules.workflows.stage_procedures as stage_procedures


def _host_guest_procedure() -> stage_procedures.Composite_stage_procedure:
    equil = stage_procedures.Explicit_stage_procedure(
        name="equilibration",
        items=[stage_procedures.MD_stage_item(name="equilibration")])
    seeding = stage_procedures.Seeding_stage_procedure(
        name="metadynamics", cv_names=["distance"])
    mmvt = stage_procedures.MMVT_stage_procedure(name="MMVT")
    bd = stage_procedures.BD_stage_procedure(name="bd")
    return stage_procedures.Composite_stage_procedure(
        name="full_procedure",
        procedures=[equil, seeding, mmvt, bd])


def _nested_composite_procedure() -> stage_procedures.Composite_stage_procedure:
    grandchild = stage_procedures.MMVT_stage_procedure(name="mmvt_inner")
    child = stage_procedures.Composite_stage_procedure(
        name="child", procedures=[grandchild])
    return stage_procedures.Composite_stage_procedure(
        name="root", procedures=[child])


class TestBuildStageAddressMap:
    def test_host_guest_paths(self):
        proc = _host_guest_procedure()
        address_map = stage_procedures.build_stage_address_map(proc)
        assert address_map["equilibration"] == (
            ["equilibration", "equilibration"], "equilibration")
        assert address_map["metadynamics_seed"] == (
            ["metadynamics", "sampling"], "sampling")
        assert address_map["metadynamics_assign_structures"] == (
            ["metadynamics", "logistic"], "logistic")
        assert address_map["MMVT"] == (["MMVT"], None)
        assert address_map["bd"] == (["bd"], None)

    def test_nested_composite(self):
        proc = _nested_composite_procedure()
        address_map = stage_procedures.build_stage_address_map(proc)
        assert address_map["mmvt_inner"] == (["child", "mmvt_inner"], None)

    def test_duplicate_stage_names_raise(self):
        dup_a = stage_procedures.BD_stage_procedure(name="bd")
        dup_b = stage_procedures.MMVT_stage_procedure(name="bd")
        proc = stage_procedures.Composite_stage_procedure(
            name="root", procedures=[dup_a, dup_b])
        with pytest.raises(ValueError, match="Duplicate stage name"):
            stage_procedures.build_stage_address_map(proc)


class TestBuildStagePolicyMap:
    def test_mmvt_anchor_dimensions(self):
        proc = stage_procedures.MMVT_stage_procedure(name="MMVT")
        policy = stage_procedures.build_stage_policy_map(proc)
        assert policy["MMVT"].dispatch.dimensions == ["anchor"]
        assert policy["MMVT"].dispatch.group_size == 1

    def test_bd_generic_fallback(self):
        proc = stage_procedures.BD_stage_procedure(name="bd")
        policy = stage_procedures.build_stage_policy_map(proc)
        emitted = policy["bd"]
        assert emitted.dispatch.dimensions is None
        assert emitted.dispatch.group_size is None
        assert emitted.dispatch.concurrency == 1
        assert emitted.co_schedule_with is None

    def test_seeding_logistic_co_schedule(self):
        proc = stage_procedures.Seeding_stage_procedure(
            name="metadynamics", cv_names=["distance"])
        policy = stage_procedures.build_stage_policy_map(proc)
        assert policy["metadynamics_seed"].co_schedule_with is None
        assert policy["metadynamics_assign_structures"].co_schedule_with == (
            "predecessor")

    def test_ramd_swarm_dimensions(self):
        proc = stage_procedures.RAMD_stage_procedure(name="ramd_proc")
        policy = stage_procedures.build_stage_policy_map(proc)
        assert policy["ramd_proc_ramd"].dispatch.dimensions == ["swarm"]
        assert policy["ramd_proc_ramd"].dispatch.group_size == 1

    def test_explicit_and_equilibration_generic(self):
        equil = stage_procedures.Equilibration_globular_stage_procedure(
            name="equilibration")
        explicit = stage_procedures.Explicit_stage_procedure(
            name="explicit",
            items=[stage_procedures.MD_stage_item(name="sampling")])
        for proc in (equil, explicit):
            policy = stage_procedures.build_stage_policy_map(proc)
            for emitted in policy.values():
                assert emitted.dispatch.dimensions is None
                assert emitted.dispatch.group_size is None
                assert emitted.co_schedule_with is None


class TestResolverPrecedence:
    def _run_settings(self, **kwargs) -> structures.Run_settings:
        return structures.Run_settings(**kwargs)

    def test_role_default_when_no_placement(self):
        proc = _host_guest_procedure()
        rs = self._run_settings()
        resolved = rs.resolve_stage_execution("metadynamics_assign_structures", proc)
        assert resolved.resource_name == "local"
        assert resolved.co_schedule_with == "predecessor"

    def test_broad_placement_applies_to_subroles(self):
        proc = _host_guest_procedure()
        rs = self._run_settings(
            resources=[structures.Resource_remote_slurm(name="delta")],
            placements=[
                structures.Placement(target=["metadynamics"], resource="delta"),
            ])
        seed = rs.resolve_stage_execution("metadynamics_seed", proc)
        assign = rs.resolve_stage_execution("metadynamics_assign_structures", proc)
        assert seed.resource_name == "delta"
        assert assign.resource_name == "delta"

    def test_specific_placement_overrides_subrole(self):
        proc = _host_guest_procedure()
        rs = self._run_settings(
            resources=[structures.Resource_remote_slurm(name="delta")],
            placements=[
                structures.Placement(target=["metadynamics"], resource="delta"),
                structures.Placement(
                    target=["metadynamics", "logistic"], resource="local"),
            ])
        seed = rs.resolve_stage_execution("metadynamics_seed", proc)
        assign = rs.resolve_stage_execution("metadynamics_assign_structures", proc)
        assert seed.resource_name == "delta"
        assert assign.resource_name == "local"

    def test_duplicate_targets_raise(self):
        proc = _host_guest_procedure()
        rs = self._run_settings(placements=[
            structures.Placement(target=["bd"], resource="local"),
            structures.Placement(target=["bd"], resource="delta"),
        ])
        with pytest.raises(ValueError, match="Duplicate placement target"):
            rs.resolve_stage_execution("bd", proc)

    def test_dispatch_merge_preserves_role_dimensions(self):
        ramd = stage_procedures.RAMD_stage_procedure(name="ramd_proc")
        proc = stage_procedures.Composite_stage_procedure(
            name="root", procedures=[ramd])
        rs = self._run_settings(
            resources=[structures.Resource_remote_slurm(name="cluster")],
            placements=[
                structures.Placement(
                    target=["ramd_proc"],
                    resource="cluster",
                    dispatch=stage_procedures.Dispatch(group_size=2)),
            ])
        resolved = rs.resolve_stage_execution("ramd_proc_ramd", proc)
        assert resolved.dispatch.dimensions == ["swarm"]
        assert resolved.dispatch.group_size == 2


class TestCapabilityClamp:
    def test_local_collapses_group_size(self):
        proc = _host_guest_procedure()
        rs = structures.Run_settings(placements=[
            structures.Placement(
                target=["MMVT"],
                resource="local",
                dispatch=stage_procedures.Dispatch(group_size=4)),
        ])
        resolved = rs.resolve_stage_execution("MMVT", proc)
        assert resolved.dispatch.group_size is None

    def test_concurrency_exceeds_mps_raises(self):
        proc = _host_guest_procedure()
        resource = structures.Resource_remote_slurm(name="cluster", mps=2)
        rs = structures.Run_settings(
            resources=[resource],
            placements=[
                structures.Placement(
                    target=["MMVT"],
                    resource="cluster",
                    dispatch=stage_procedures.Dispatch(concurrency=8)),
            ])
        with pytest.raises(ValueError, match="dispatch.concurrency=8"):
            rs.resolve_stage_execution("MMVT", proc)


class _FakeStage:
    def __init__(self, name: str, input_stage_index: int = 0):
        self.name = name
        self.input_stage_index = input_stage_index


class _FakeModel:
    def __init__(self, stages: list[_FakeStage]):
        self.stages = stages


class TestValidateRunSettings:
    def test_unknown_target_raises(self):
        proc = _host_guest_procedure()
        seekrflow = structures.Seekrflow()
        seekrflow.workflow.procedure = proc
        seekrflow.run_settings.placements = [
            structures.Placement(target=["nonexistent"], resource="local"),
        ]
        with pytest.raises(ValueError, match="Unknown placement target"):
            structures.validate_run_settings(seekrflow)

    def test_unknown_resource_raises(self):
        proc = _host_guest_procedure()
        seekrflow = structures.Seekrflow()
        seekrflow.workflow.procedure = proc
        seekrflow.run_settings.placements = [
            structures.Placement(target=["bd"], resource="missing"),
        ]
        with pytest.raises(ValueError, match="not found"):
            structures.validate_run_settings(seekrflow)

    def test_co_schedule_requires_same_resource(self):
        proc = _host_guest_procedure()
        seekrflow = structures.Seekrflow()
        seekrflow.workflow.procedure = proc
        seekrflow.run_settings.resources = [
            structures.Resource_remote_slurm(name="delta")]
        seekrflow.run_settings.placements = [
            structures.Placement(
                target=["metadynamics", "sampling"], resource="delta"),
            structures.Placement(
                target=["metadynamics", "logistic"], resource="local"),
        ]
        model = _FakeModel([
            _FakeStage("metadynamics_seed", input_stage_index=0),
            _FakeStage("metadynamics_assign_structures", input_stage_index=1),
        ])
        with pytest.raises(ValueError, match="same resource"):
            structures.validate_run_settings(seekrflow, model)


class TestRoundTrip:
    def test_save_and_load_placements(self):
        proc = _host_guest_procedure()
        seekrflow = structures.Seekrflow()
        seekrflow.workflow.procedure = proc
        seekrflow.run_settings.placements = [
            structures.Placement(target=[], resource="local"),
            structures.Placement(
                target=["metadynamics", "sampling"],
                resource="delta",
                dispatch=stage_procedures.Dispatch(
                    dimensions=["anchor"], group_size=2, concurrency=1)),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "seekrflow.json")
            seekrflow.save(path)
            loaded = structures.load_seekrflow(path)
        assert loaded.run_settings.placements[1].target == [
            "metadynamics", "sampling"]
        assert loaded.run_settings.placements[1].dispatch.group_size == 2


class TestSetPlacementResource:
    def test_cli_style_override(self):
        rs = structures.Run_settings(
            resources=[structures.Resource_remote_slurm(name="delta")])
        rs.set_placement_resource(
            ["metadynamics", "logistic"], "local")
        assert rs.placements[0].target == ["metadynamics", "logistic"]
        assert rs.placements[0].resource == "local"
        rs.set_placement_resource(
            ["metadynamics", "logistic"], "delta")
        assert len(rs.placements) == 1
        assert rs.placements[0].resource == "delta"
