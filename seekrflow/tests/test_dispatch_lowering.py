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


class TestFilterIncompleteUnits:
    def test_drops_attained_anchors(self):
        units = dl.build_run_units(["anchor"], _counts(num_anchors=4))
        progress = {
            0: {"progress": 1.0, "attained": True},
            1: {"progress": 0.5, "attained": False},
            "2": {"progress": 0.0, "attained": False},
            # 3 missing → incomplete
        }
        kept = dl.filter_incomplete_units(units, progress)
        assert [u.anchor for u in kept] == [1, 2, 3]

    def test_stage_finished_returns_empty(self):
        units = dl.build_run_units(["anchor"], _counts(num_anchors=3))
        assert dl.filter_incomplete_units(
            units, {}, stage_finished=True) == []

    def test_force_overwrite_keeps_all(self):
        units = dl.build_run_units(["anchor"], _counts(num_anchors=3))
        progress = {0: {"attained": True}, 1: {"attained": True}, 2: {"attained": True}}
        kept = dl.filter_incomplete_units(
            units, progress, force_overwrite=True)
        assert len(kept) == 3

    def test_progress_one_without_attained(self):
        units = [dl.RunUnit(anchor=0, swarm_id=None)]
        assert dl.filter_incomplete_units(
            units, {0: {"progress": 1.0}}) == []

    def test_swarm_keys(self):
        units = dl.build_run_units(
            ["swarm"],
            _counts(scope="unpartitioned", num_anchors=None, num_swarms=3))
        progress = {
            0: {"attained": True},
            1: {"progress": 0.2},
        }
        kept = dl.filter_incomplete_units(units, progress)
        assert [u.swarm_id for u in kept] == [1, 2]

    def test_empty_progress_keeps_all(self):
        units = dl.build_run_units(["anchor"], _counts(num_anchors=2))
        assert len(dl.filter_incomplete_units(units, None)) == 2
        assert len(dl.filter_incomplete_units(units, {})) == 2


class TestLowerFromFilteredUnits:
    def test_dense_array_from_subset(self):
        units = [
            dl.RunUnit(anchor=1, swarm_id=None),
            dl.RunUnit(anchor=3, swarm_id=None),
            dl.RunUnit(anchor=7, swarm_id=None),
        ]
        lowered = dl.lower_dispatch_from_units(units, group_size=1, concurrency=1)
        assert lowered.n_members == 3
        assert lowered.use_array
        assert dl.array_member_indices(lowered) == [0, 1, 2]
        assert [u.anchor for u in lowered.units] == [1, 3, 7]

    def test_empty_units(self):
        lowered = dl.lower_dispatch_from_units([], group_size=2, concurrency=1)
        assert lowered.n_members == 0
        assert lowered.units == ()
        assert not lowered.use_array


class TestFetchStageProgressLocal:
    def test_parses_finished_and_map(self):
        class _Stage:
            name = "MMVT"
            index = 2

        def fake_status(model, instruction, stage_arg=None, print_json=True):
            assert instruction == "progress"
            return {
                "progress": {
                    2: {
                        "finished": False,
                        "progress": {
                            0: {"attained": True, "progress": 1.0},
                            1: {"attained": False, "progress": 0.1},
                        },
                    }
                }
            }

        finished, pmap = dl.fetch_stage_progress_local(
            None, _Stage(), status_fn=fake_status)
        assert finished is False
        assert 0 in pmap and pmap[0]["attained"] is True


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
        assert "out_${AWS_BATCH_JOB_ARRAY_INDEX" in cmd
        assert cmd.rstrip().endswith(".log") or ".log" in cmd.split(">")[-1]

    def test_filtered_units_embedded_in_driver(self):
        units = [
            dl.RunUnit(anchor=2, swarm_id=None),
            dl.RunUnit(anchor=5, swarm_id=None),
        ]
        lowered = dl.lower_dispatch_from_units(units, group_size=1, concurrency=1)
        cmd = dl.build_remote_command_string(
            "MMVT", lowered, force_overwrite=False, benchmark=False,
            output_filename="out.log")
        # json.dumps embeds null for swarm_id=None
        assert "[[2, null], [5, null]]" in cmd
        assert "group_size = 1" in cmd
