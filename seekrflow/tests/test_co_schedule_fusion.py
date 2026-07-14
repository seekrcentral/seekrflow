"""
Tests for co_schedule_with job-script fusion (Phase 8d).
"""
from __future__ import annotations

import pytest

from seekrflow.modules import structures
from seekrflow.modules.workflows import stage_procedures
from seekrflow.modules.workload_managers import co_schedule_fusion


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


class _MockStage:
    def __init__(
            self,
            name: str,
            index: int,
            input_stage_index: int = 0,
            ):
        self.name = name
        self.index = index
        self.input_stage_index = input_stage_index


class _MockStageWorkflow:
    def __init__(
            self,
            stage: _MockStage,
            *,
            resource_name: str = "cluster",
            co_schedule_with: str | None = None,
            benchmark_mode: bool = False,
            resolved_execution: structures.Resolved_execution | None = None,
            ):
        self.stage = stage
        self.resource_name = resource_name
        self.co_schedule_with = co_schedule_with
        self.benchmark_mode = benchmark_mode
        self.resolved_execution = resolved_execution
        self.resource = (
            structures.Resource_remote_slurm(name="cluster")
            if resource_name != "local" else None)
        self.fusion_host = None
        self.fused_before: list[str] = []
        self.fused_after: list[str] = []
        self.peer_workflows: dict = {}
        self.state = "unstarted"
        self.manager_status = "idle"


class TestCombineFusedCommands:
    def test_predecessor_order(self):
        host = "python -c host > /work/host_run.out"
        fused = "python -c fused > /work/fused_run.out"
        combined = co_schedule_fusion.combine_fused_commands([host, fused])
        assert combined == f"({host}) && ({fused})"

    def test_successor_order(self):
        fused = "python -c fused > /work/fused_run.out"
        host = "python -c host > /work/host_run.out"
        combined = co_schedule_fusion.combine_fused_commands([fused, host])
        assert combined == f"({fused}) && ({host})"

    def test_each_command_keeps_redirect(self):
        cmds = [
            "python -c a > /work/a_run.out",
            "python -c b > /work/b_run.out",
        ]
        combined = co_schedule_fusion.combine_fused_commands(cmds)
        assert "> /work/a_run.out" in combined
        assert "> /work/b_run.out" in combined


class TestPopulateFusionMap:
    def _resolved(
            self,
            resource_name: str = "cluster",
            co_schedule_with: str | None = None,
            ) -> structures.Resolved_execution:
        resource = (
            structures.Resource_remote_slurm(name=resource_name)
            if resource_name != "local" else None)
        return structures.Resolved_execution(
            stage_name="unused",
            resource_name=resource_name,
            resource=resource,
            dispatch=stage_procedures.Dispatch(),
            co_schedule_with=co_schedule_with,
            cpus_per_task=None,
            memory_per_node=None,
            time_limit=None,
            mps=None,
        )

    def test_remote_logistic_fuses_into_predecessor(self):
        seed = _MockStageWorkflow(
            _MockStage("metadynamics_seed", index=1, input_stage_index=0),
            resolved_execution=self._resolved(),
        )
        logistic = _MockStageWorkflow(
            _MockStage(
                "metadynamics_assign_structures",
                index=2,
                input_stage_index=1,
            ),
            resolved_execution=self._resolved(
                co_schedule_with="predecessor"),
        )
        workflows = [seed, logistic]
        co_schedule_fusion.populate_fusion_map(workflows)
        assert seed.fused_after == ["metadynamics_assign_structures"]
        assert logistic.fusion_host == "metadynamics_seed"
        assert co_schedule_fusion.skips_remote_submit(logistic)
        assert not co_schedule_fusion.skips_remote_submit(seed)

    def test_local_co_schedule_is_no_op(self):
        seed = _MockStageWorkflow(
            _MockStage("metadynamics_seed", index=1, input_stage_index=0),
            resource_name="local",
            resolved_execution=self._resolved(resource_name="local"),
        )
        logistic = _MockStageWorkflow(
            _MockStage(
                "metadynamics_assign_structures",
                index=2,
                input_stage_index=1,
            ),
            resource_name="local",
            resolved_execution=self._resolved(
                resource_name="local",
                co_schedule_with="predecessor"),
        )
        workflows = [seed, logistic]
        co_schedule_fusion.populate_fusion_map(workflows)
        assert seed.fused_after == []
        assert logistic.fusion_host is None
        assert not co_schedule_fusion.skips_remote_submit(logistic)

    def test_successor_fusion(self):
        fused = _MockStageWorkflow(
            _MockStage("prep_stage", index=1, input_stage_index=0),
            resolved_execution=self._resolved(
                co_schedule_with="successor"),
        )
        host = _MockStageWorkflow(
            _MockStage("host_stage", index=2, input_stage_index=1),
            resolved_execution=self._resolved(),
        )
        workflows = [fused, host]
        co_schedule_fusion.populate_fusion_map(workflows)
        assert host.fused_before == ["prep_stage"]
        assert fused.fusion_host == "host_stage"

    def test_chained_predecessor_fusion(self):
        host = _MockStageWorkflow(
            _MockStage("host", index=1, input_stage_index=0),
            resolved_execution=self._resolved(),
        )
        middle = _MockStageWorkflow(
            _MockStage("middle", index=2, input_stage_index=1),
            resolved_execution=self._resolved(
                co_schedule_with="predecessor"),
        )
        tail = _MockStageWorkflow(
            _MockStage("tail", index=3, input_stage_index=2),
            resolved_execution=self._resolved(
                co_schedule_with="predecessor"),
        )
        workflows = [host, middle, tail]
        co_schedule_fusion.populate_fusion_map(workflows)
        assert host.fused_after == ["middle", "tail"]
        assert middle.fusion_host == "host"
        assert tail.fusion_host == "host"
        assert co_schedule_fusion.fusion_host_name(host) == "host"
        assert co_schedule_fusion.is_fusion_host(host)
        assert co_schedule_fusion.is_fusion_member(tail)


class TestFusedSetHelpers:
    def _by_name(self, *workflows):
        for sw in workflows:
            sw.peer_workflows = {w.stage.name: w for w in workflows}
        return {sw.stage.name: sw for sw in workflows}

    def test_fused_set_completed_all_done(self):
        host = _MockStageWorkflow(_MockStage("host", index=1))
        fused = _MockStageWorkflow(_MockStage("fused", index=2))
        fused.fusion_host = "host"
        host.fused_after = ["fused"]
        host.state = "completed"
        fused.state = "completed"
        by_name = self._by_name(host, fused)
        assert co_schedule_fusion.fused_set_completed("host", by_name)

    def test_fused_set_completed_partial(self):
        host = _MockStageWorkflow(_MockStage("host", index=1))
        fused = _MockStageWorkflow(_MockStage("fused", index=2))
        fused.fusion_host = "host"
        host.fused_after = ["fused"]
        host.state = "completed"
        fused.state = "started"
        by_name = self._by_name(host, fused)
        assert not co_schedule_fusion.fused_set_completed("host", by_name)

    def test_fused_set_progress_excludes_logistic(self):
        host = _MockStageWorkflow(_MockStage("host", index=1))
        host.stage.scale_type = "molecular_dynamics"
        host.progress = 1.0
        logistic = _MockStageWorkflow(_MockStage("logistic", index=2))
        logistic.co_schedule_with = "predecessor"
        logistic.progress = 0.0
        host.fused_after = ["logistic"]
        by_name = self._by_name(host, logistic)
        assert co_schedule_fusion.fused_set_progress("host", by_name) == 1.0

    def test_classify_fused_probe_all_complete(self):
        statuses = {
            "host": {
                "stage_status": {"state": "completed"},
                "manager_status": {"jobs": []},
            },
            "fused": {
                "stage_status": {"state": "completed"},
                "manager_status": {"jobs": []},
            },
        }
        assert co_schedule_fusion.classify_fused_probe_status(
            statuses, "host") == "completed"

    def test_classify_fused_probe_reattach(self):
        statuses = {
            "host": {
                "stage_status": {"state": "completed"},
                "manager_status": {"jobs": [{"JobID": "1"}]},
            },
            "fused": {
                "stage_status": {"state": "started"},
                "manager_status": {"jobs": []},
            },
        }
        assert co_schedule_fusion.classify_fused_probe_status(
            statuses, "host") == "reattach"

    def test_classify_fused_probe_submit(self):
        statuses = {
            "host": {
                "stage_status": {"state": "completed"},
                "manager_status": {"jobs": []},
            },
            "fused": {
                "stage_status": {"state": "started"},
                "manager_status": {"jobs": []},
            },
        }
        assert co_schedule_fusion.classify_fused_probe_status(
            statuses, "host") == "submit"


class TestFusionDependenciesSatisfied:
    def test_waits_for_host_started(self):
        host = _MockStageWorkflow(
            _MockStage("host", index=1, input_stage_index=0))
        fused = _MockStageWorkflow(
            _MockStage("fused", index=2, input_stage_index=1))
        fused.fusion_host = "host"
        by_name = {"host": host, "fused": fused}
        assert not co_schedule_fusion.fusion_dependencies_satisfied(fused, by_name)
        host.state = "started"
        assert co_schedule_fusion.fusion_dependencies_satisfied(fused, by_name)


class TestValidateRunSettingsFusionGuards:
    def test_array_dimension_on_fused_stage_raises(self):
        proc = _host_guest_procedure()
        seekrflow = structures.Seekrflow()
        seekrflow.workflow.procedure = proc
        seekrflow.run_settings.resources = [
            structures.Resource_remote_slurm(name="cluster")]
        seekrflow.run_settings.placements = [
            structures.Placement(target=["metadynamics"], resource="cluster"),
        ]
        model = type("M", (), {"stages": [
            type("S", (), {
                "name": "metadynamics_seed",
                "input_stage_index": 0,
            })(),
            type("S", (), {
                "name": "metadynamics_assign_structures",
                "input_stage_index": 1,
            })(),
        ]})()
        seekrflow.run_settings.placements.append(
            structures.Placement(
                target=["metadynamics", "logistic"],
                resource="cluster",
                dispatch=stage_procedures.Dispatch(dimensions=["swarm"]),
            ))
        with pytest.raises(ValueError, match="cannot be co-scheduled"):
            structures.validate_run_settings(seekrflow, model)

    def test_predecessor_without_predecessor_raises(self):
        proc = _host_guest_procedure()
        seekrflow = structures.Seekrflow()
        seekrflow.workflow.procedure = proc
        seekrflow.run_settings.resources = [
            structures.Resource_remote_slurm(name="cluster")]
        seekrflow.run_settings.placements = [
            structures.Placement(target=["metadynamics"], resource="cluster"),
        ]
        model = type("M", (), {"stages": [
            type("S", (), {
                "name": "metadynamics_assign_structures",
                "input_stage_index": 0,
            })(),
        ]})()
        with pytest.raises(ValueError, match="no predecessor"):
            structures.validate_run_settings(seekrflow, model)

    def test_valid_co_schedule_passes(self):
        proc = _host_guest_procedure()
        seekrflow = structures.Seekrflow()
        seekrflow.workflow.procedure = proc
        seekrflow.run_settings.resources = [
            structures.Resource_remote_slurm(name="cluster")]
        seekrflow.run_settings.placements = [
            structures.Placement(target=["metadynamics"], resource="cluster"),
        ]
        model = type("M", (), {"stages": [
            type("S", (), {
                "name": "metadynamics_seed",
                "input_stage_index": 0,
            })(),
            type("S", (), {
                "name": "metadynamics_assign_structures",
                "input_stage_index": 1,
            })(),
        ]})()
        structures.validate_run_settings(seekrflow, model)
