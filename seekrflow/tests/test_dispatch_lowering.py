"""Tests for Phase 8 dispatch lowering (pure logic)."""
import math

import pytest

from seekrflow.modules.workload_managers import dispatch_lowering as dl


def _counts(
        scope: str = "partitioned",
        num_anchors: int | None = 8,
        num_swarms: int = 1,
        ) -> dl.StageUnitCounts:
    return dl.StageUnitCounts(
        scope=scope, num_anchors=num_anchors, num_swarms=num_swarms)


class TestBuildRunUnits:
    def test_anchor_dimension(self):
        units = dl.build_run_units(["anchor"], _counts(num_anchors=8))
        assert len(units) == 8
        assert units[0].anchor == 0 and units[0].swarm_id is None

    def test_swarm_dimension_none_anchors(self):
        units = dl.build_run_units(
            ["swarm"], _counts(scope="unpartitioned", num_anchors=None, num_swarms=100))
        assert len(units) == 100
        assert all(u.anchor == "any" for u in units)

    def test_anchor_swarm_cross_product(self):
        units = dl.build_run_units(
            ["anchor", "swarm"], _counts(num_anchors=4, num_swarms=3))
        assert len(units) == 12
        assert units[0] == dl.RunUnit(anchor=0, swarm_id=0)
        assert units[3] == dl.RunUnit(anchor=1, swarm_id=0)
        assert units[-1] == dl.RunUnit(anchor=3, swarm_id=2)

    def test_empty_dimensions_single_unit(self):
        units = dl.build_run_units(None, _counts())
        assert units == [dl.RunUnit(anchor="any", swarm_id=None)]


class TestValidation:
    def test_anchor_on_unpartitioned_raises(self):
        with pytest.raises(ValueError, match="partitioned"):
            dl.build_run_units(
                ["anchor"],
                _counts(scope="unpartitioned", num_anchors=None))

    def test_swarm_requires_multiple(self):
        with pytest.raises(ValueError, match="num_swarms"):
            dl.build_run_units(["swarm"], _counts(num_swarms=1))

    def test_zero_swarms_raises(self):
        with pytest.raises(ValueError, match="num_swarms=0"):
            dl.resolve_unit_counts(
                {"stage_name": "seed", "num_swarms": 0},
                None,
                input_is_initial=False)


class TestArraySpecMath:
    def test_member_count(self):
        assert dl.member_count(8, 2) == 4
        assert dl.member_count(8, 3) == 3
        assert dl.member_count(8, None) == 1

    def test_member_slice(self):
        assert list(range(8)[dl.member_unit_slice(0, 8, 3)]) == [0, 1, 2]
        assert list(range(8)[dl.member_unit_slice(2, 8, 3)]) == [6, 7]

    def test_lower_dispatch_members(self):
        lowered = dl.lower_dispatch(
            ["anchor"], group_size=3, concurrency=1, counts=_counts(num_anchors=8))
        assert lowered.n_members == math.ceil(8 / 3)
        assert dl.array_member_indices(lowered) == [0, 1, 2]


class TestConcurrencyWaves:
    def test_three_waves_of_two(self):
        waves = dl.wave_slices(6, 2)
        assert len(waves) == 3
        assert list(range(6))[waves[0]] == [0, 1]
        assert list(range(6))[waves[1]] == [2, 3]
        assert list(range(6))[waves[2]] == [4, 5]

    def test_single_wave_when_equal(self):
        waves = dl.wave_slices(4, 4)
        assert len(waves) == 1
        assert list(range(4))[waves[0]] == [0, 1, 2, 3]


class TestInfoEnumeration:
    def test_match_by_stage_name_not_key(self):
        message = {
            "info": {
                "4": {
                    "stage_name": "MMVT",
                    "scope": "partitioned",
                    "num_anchors": 8,
                    "num_swarms": 1,
                }
            }
        }
        entry = dl.extract_info_entry(message, stage_name="MMVT")
        assert entry["num_anchors"] == 8

    def test_input_stage_query(self):
        class _Stage:
            name = "MMVT"
            index = 4
            input_stage_index = 3

        calls = []

        def fake_status(model, instruction, stage_arg, print_json=True):
            calls.append(stage_arg)
            if stage_arg == 3:
                return {"info": {"3": {
                    "stage_name": "seed",
                    "scope": "partitioned",
                    "num_anchors": 8,
                    "num_swarms": 2,
                }}}
            raise AssertionError(stage_arg)

        counts = dl.fetch_unit_counts_local(None, _Stage(), status_fn=fake_status)
        assert calls == [3]
        assert counts.num_swarms == 2
        assert counts.num_anchors == 8

    def test_initial_stage_uses_launching_scope(self):
        class _Stage:
            name = "MMVT"
            index = 2
            input_stage_index = 0

        def fake_status(model, instruction, stage_arg, print_json=True):
            if stage_arg == 0:
                return {"info": {"0": {
                    "stage_name": "initial",
                    "scope": "N/A",
                    "num_anchors": 12,
                    "num_swarms": 3,
                }}}
            if stage_arg == 2:
                return {"info": {"2": {
                    "stage_name": "MMVT",
                    "scope": "partitioned",
                    "num_anchors": 8,
                    "num_swarms": 1,
                }}}
            raise AssertionError(stage_arg)

        counts = dl.fetch_unit_counts_local(None, _Stage(), status_fn=fake_status)
        assert counts.num_swarms == 3
        assert counts.scope == "partitioned"
        assert counts.num_anchors == 8

    def test_initial_unpartitioned_launch_rejects_anchor(self):
        class _Stage:
            name = "bd"
            index = 1
            input_stage_index = 0

        def fake_status(model, instruction, stage_arg, print_json=True):
            if stage_arg == 0:
                return {"info": {"0": {
                    "stage_name": "initial",
                    "scope": "N/A",
                    "num_anchors": 12,
                    "num_swarms": 2,
                }}}
            return {"info": {"1": {
                "stage_name": "bd",
                "scope": "unpartitioned",
                "num_anchors": None,
                "num_swarms": 1,
            }}}

        counts = dl.fetch_unit_counts_local(None, _Stage(), status_fn=fake_status)
        with pytest.raises(ValueError, match="partitioned"):
            dl.build_run_units(["anchor"], counts)


class TestCommandString:
    def test_single_unit_default(self):
        lowered = dl.lower_dispatch(None, None, 1, _counts())
        cmd = dl.build_remote_command_string(
            "MMVT", lowered, force_overwrite=False, benchmark=False,
            output_filename="out.log")
        assert "'any'" in cmd
        assert "swarm_id" not in cmd or ", None," in cmd

    def test_multi_unit_uses_driver(self):
        lowered = dl.lower_dispatch(
            ["anchor"], group_size=2, concurrency=1, counts=_counts(num_anchors=4))
        cmd = dl.build_remote_command_string(
            "ramd", lowered, force_overwrite=False, benchmark=False,
            output_filename="out.log")
        assert "SLURM_ARRAY_TASK_ID" in cmd
        assert "group_size = 2" in cmd
