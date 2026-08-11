"""Tests for adaptive walltime helpers and Placement time_policy resolve."""

import pytest

import seekrflow.modules.structures as structures
import seekrflow.modules.workload_managers.walltime as walltime
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


def _slurm_resource(name: str = "cluster", time_limit: str = "24:00:00"):
    return structures.Resource_remote_slurm(
        name=name, time_limit=time_limit, mps=4)


def _aws_resource(name: str = "aws_gpu", timeout: int = 86400):
    return structures.Resource_cloud_aws(
        name=name,
        job_timeout_seconds=timeout,
        transfer_settings=structures.Transfer_settings_aws_s3(
            bucket="test-bucket",
            prefix="test",
        ),
    )


class TestTimePolicyResolve:
    def test_unset_policy_is_adaptive_with_none_time_limit(self):
        rs = structures.Run_settings(
            resources=[_slurm_resource()],
            placements=[
                structures.Placement(target=["MMVT"], resource="cluster"),
            ],
        )
        resolved = rs.resolve_stage_execution(
            "MMVT", _host_guest_procedure())
        assert isinstance(
            resolved.time_policy, structures.Time_policy_adaptive)
        assert resolved.time_limit is None

    def test_fixed_uses_resource_cap_when_unset(self):
        rs = structures.Run_settings(
            resources=[_slurm_resource(time_limit="12:00:00")],
            placements=[
                structures.Placement(
                    target=["MMVT"],
                    resource="cluster",
                    time_policy=structures.Time_policy_fixed(),
                ),
            ],
        )
        resolved = rs.resolve_stage_execution(
            "MMVT", _host_guest_procedure())
        assert isinstance(resolved.time_policy, structures.Time_policy_fixed)
        assert resolved.time_limit == "12:00:00"

    def test_fixed_override_and_cap_validation(self):
        rs = structures.Run_settings(
            resources=[_slurm_resource(time_limit="12:00:00")],
            placements=[
                structures.Placement(
                    target=["MMVT"],
                    resource="cluster",
                    time_policy=structures.Time_policy_fixed(
                        time_limit="48:00:00"),
                ),
            ],
        )
        with pytest.raises(ValueError, match="exceeds resource"):
            rs.resolve_stage_execution("MMVT", _host_guest_procedure())

    def test_adaptive_placement_time_limit_becomes_max(self):
        rs = structures.Run_settings(
            resources=[_slurm_resource(time_limit="24:00:00")],
            placements=[
                structures.Placement(
                    target=["MMVT"],
                    resource="cluster",
                    time_limit="06:00:00",
                    time_policy=structures.Time_policy_adaptive(
                        estimated_performance=50.0),
                ),
            ],
        )
        resolved = rs.resolve_stage_execution(
            "MMVT", _host_guest_procedure())
        assert resolved.time_limit is None
        assert isinstance(
            resolved.time_policy, structures.Time_policy_adaptive)
        assert resolved.time_policy.max_time_limit == "06:00:00"
        assert resolved.time_policy.estimated_performance == 50.0


class TestWalltimeMath:
    def test_timestep_prefers_get_timestep_by_scale(self):
        class _Model:
            def get_timestep_by_scale(self, scale_type):
                assert scale_type == "molecular_dynamics"
                return 0.002

        assert walltime.molecular_dynamics_timestep_ps(_Model()) == 0.002

    def test_steps_to_ns_uses_picoseconds(self):
        # 500_000 steps * 0.002 ps = 1000 ps = 1 ns
        assert walltime.steps_to_ns(500_000, 0.002) == pytest.approx(1.0)
        progress = {
            "0": {
                "completion_criteria": "number_of_steps",
                "attained": False,
                "swarms": {
                    "0": {
                        "completion_criteria": "number_of_steps",
                        "attained": False,
                        "current_step": 250000,
                        "total_steps": 500000,
                    }
                },
            }
        }
        # 250000 steps * 0.002 ps = 500 ps = 0.5 ns
        kind, amount = walltime.extract_remaining_work(
            progress, timestep_ps=0.002)
        assert kind == "ns"
        assert amount == pytest.approx(0.5)

    def test_extract_remaining_bd(self):
        progress = {
            "null": {
                "completion_criteria": "number_of_trajectories",
                "attained": False,
                "trajectories_needed": 10000,
                "trajectories_completed": 2500,
            }
        }
        kind, amount = walltime.extract_remaining_work(progress)
        assert kind == "trajectories"
        assert amount == 7500.0

    def test_reporter_progress_not_countable(self):
        progress = {
            "null": {
                "completion_criteria": "reporter_progress",
                "attained": False,
                "swarms": {
                    "0": {
                        "completion_criteria": "reporter_progress",
                        "attained": False,
                        "current_step": 100,
                        "maximum_steps": 1_000_000_000,
                    }
                },
            }
        }
        assert walltime.extract_remaining_work(progress) is None

    def test_estimate_with_hint_pads_and_rounds_hour(self):
        resource = _slurm_resource(time_limit="24:00:00")
        resolved = structures.Resolved_execution(
            stage_name="MMVT",
            resource_name="cluster",
            resource=resource,
            dispatch=stage_procedures.Dispatch(),
            co_schedule_with=None,
            cpus=1,
            memory_mb=4000,
            time_limit=None,
            mps=1,
            time_policy=structures.Time_policy_adaptive(
                estimated_performance=12.0,  # ns/day
                safety_factor=0.2,
            ),
        )
        progress = {
            "0": {
                "completion_criteria": "number_of_steps",
                "attained": False,
                "swarms": {
                    "0": {
                        "completion_criteria": "number_of_steps",
                        "attained": False,
                        "current_step": 0,
                        "total_steps": 500_000,  # 1 ns at 0.002 ps
                    }
                },
            }
        }
        # 1 ns / 12 ns/day = 2h; +20% = 2.4h → round up to 3h
        seconds = walltime.estimate_stage_walltime_seconds(
            resolved=resolved,
            partition_progress=progress,
            timestep_ps=0.002,
        )
        assert seconds == 3 * 3600

    def test_aws_min_floor_ten_minutes(self):
        resource = _aws_resource(timeout=86400)
        resolved = structures.Resolved_execution(
            stage_name="equilibration",
            resource_name="aws_gpu",
            resource=resource,
            dispatch=stage_procedures.Dispatch(),
            co_schedule_with=None,
            cpus=4,
            memory_mb=16000,
            time_limit=None,
            mps=1,
            time_policy=structures.Time_policy_adaptive(
                estimated_performance=1000.0,
            ),
        )
        progress = {
            "null": {
                "completion_criteria": "number_of_steps",
                "attained": False,
                "swarms": {
                    "0": {
                        "completion_criteria": "number_of_steps",
                        "attained": False,
                        "current_step": 0,
                        "total_steps": 1000,  # tiny
                    }
                },
            }
        }
        seconds = walltime.estimate_stage_walltime_seconds(
            resolved=resolved,
            partition_progress=progress,
            timestep_ps=0.002,
        )
        assert seconds == walltime.MIN_AWS_SECONDS

    def test_measured_rate_preferred_over_hint(self):
        resource = _slurm_resource()
        resolved = structures.Resolved_execution(
            stage_name="MMVT",
            resource_name="cluster",
            resource=resource,
            dispatch=stage_procedures.Dispatch(),
            co_schedule_with=None,
            cpus=1,
            memory_mb=4000,
            time_limit=None,
            mps=1,
            time_policy=structures.Time_policy_adaptive(
                estimated_performance=1.0,  # very slow hint
            ),
        )
        progress = {
            "0": {
                "completion_criteria": "number_of_steps",
                "attained": False,
                "swarms": {
                    "0": {
                        "completion_criteria": "number_of_steps",
                        "attained": False,
                        "current_step": 250_000,
                        "total_steps": 500_000,
                    }
                },
            }
        }
        # Prior job: 250k steps in 1 hour at 0.002 ps = 0.5 ns/hour = 12 ns/day
        baseline = {"0:0": 0.0}
        seconds = walltime.estimate_stage_walltime_seconds(
            resolved=resolved,
            partition_progress=progress,
            elapsed_seconds=3600.0,
            work_baseline=baseline,
            timestep_ps=0.002,
        )
        # remaining 0.5 ns / 12 ns/day = 1h; +20% = 1.2h → 2h
        assert seconds == 2 * 3600

    def test_fused_sum(self):
        resource = _slurm_resource(time_limit="48:00:00")
        policy = structures.Time_policy_adaptive(estimated_performance=12.0)
        host = structures.Resolved_execution(
            stage_name="a",
            resource_name="cluster",
            resource=resource,
            dispatch=stage_procedures.Dispatch(),
            co_schedule_with=None,
            cpus=1,
            memory_mb=4000,
            time_limit=None,
            mps=1,
            time_policy=policy,
        )
        progress = {
            "null": {
                "completion_criteria": "number_of_steps",
                "attained": False,
                "swarms": {
                    "0": {
                        "completion_criteria": "number_of_steps",
                        "attained": False,
                        "current_step": 0,
                        "total_steps": 500_000,  # 1 ns
                    }
                },
            }
        }
        # Each stage: 1 ns / 12 ns/day = 2h; +20% = 2.4h raw.
        # Fused sum 4.8h → round up once to 5h (min/round applied once).
        result = walltime.estimate_submit_time_limit(
            host_resolved=host,
            stage_inputs=[
                {
                    "resolved": host,
                    "partition_progress": progress,
                    "timestep_ps": 0.002,
                    "scale_type": "molecular_dynamics",
                },
                {
                    "resolved": host,
                    "partition_progress": progress,
                    "timestep_ps": 0.002,
                    "scale_type": "molecular_dynamics",
                },
            ],
        )
        assert result == "05:00:00"

    def test_fused_logistic_ignored_seed_forces_full_max(self):
        resource = _slurm_resource(time_limit="24:00:00")
        policy = structures.Time_policy_adaptive(estimated_performance=12.0)
        host = structures.Resolved_execution(
            stage_name="metadynamics_seed",
            resource_name="cluster",
            resource=resource,
            dispatch=stage_procedures.Dispatch(),
            co_schedule_with=None,
            cpus=1,
            memory_mb=4000,
            time_limit=None,
            mps=1,
            time_policy=policy,
        )
        seed_progress = {
            "null": {
                "completion_criteria": "reporter_progress",
                "attained": False,
                "swarms": {
                    "0": {
                        "completion_criteria": "reporter_progress",
                        "attained": False,
                        "current_step": 100,
                        "maximum_steps": 2_000_000_000,
                    }
                },
            }
        }
        logistic_progress = {
            "null": {
                "completion_criteria": "reporter_progress",
                "attained": False,
            }
        }
        result = walltime.estimate_submit_time_limit(
            host_resolved=host,
            stage_inputs=[
                {
                    "resolved": host,
                    "partition_progress": seed_progress,
                    "timestep_ps": 0.002,
                    "scale_type": "molecular_dynamics",
                    "completion_criteria_type": "reporter_progress",
                },
                {
                    "resolved": host,
                    "partition_progress": logistic_progress,
                    "scale_type": "logistic",
                    "completion_criteria_type": None,
                },
            ],
        )
        assert result == "24:00:00"

    def test_fused_countable_plus_logistic_ignores_logistic(self):
        resource = _slurm_resource(time_limit="48:00:00")
        policy = structures.Time_policy_adaptive(estimated_performance=12.0)
        host = structures.Resolved_execution(
            stage_name="equilibration",
            resource_name="cluster",
            resource=resource,
            dispatch=stage_procedures.Dispatch(),
            co_schedule_with=None,
            cpus=1,
            memory_mb=4000,
            time_limit=None,
            mps=1,
            time_policy=policy,
        )
        equil_progress = {
            "null": {
                "completion_criteria": "number_of_steps",
                "attained": False,
                "swarms": {
                    "0": {
                        "completion_criteria": "number_of_steps",
                        "attained": False,
                        "current_step": 0,
                        "total_steps": 500_000,  # 1 ns → ~3h
                    }
                },
            }
        }
        result = walltime.estimate_submit_time_limit(
            host_resolved=host,
            stage_inputs=[
                {
                    "resolved": host,
                    "partition_progress": equil_progress,
                    "timestep_ps": 0.002,
                    "scale_type": "molecular_dynamics",
                },
                {
                    "resolved": host,
                    "partition_progress": None,
                    "scale_type": "logistic",
                },
            ],
        )
        assert result == "03:00:00"

    def test_fallback_total_steps_when_swarms_missing(self):
        resource = _aws_resource(timeout=86400)
        policy = structures.Time_policy_adaptive(estimated_performance=1000.0)
        host = structures.Resolved_execution(
            stage_name="equil_stage_2",
            resource_name="aws_gpu",
            resource=resource,
            dispatch=stage_procedures.Dispatch(),
            co_schedule_with=None,
            cpus=4,
            memory_mb=14000,
            time_limit=None,
            mps=1,
            time_policy=policy,
        )
        # Progress knows criteria but has no swarm step counts yet.
        progress = {
            "null": {
                "completion_criteria": "number_of_steps",
                "attained": False,
                "progress": 0.0,
                "num_swarms": 0,
                "error": "no starting structures (swarms) found yet",
            }
        }
        # 100000 steps * 0.004 ps = 0.4 ns / 1000 ns/day → tiny → AWS 10 min floor
        result = walltime.estimate_submit_time_limit(
            host_resolved=host,
            stage_inputs=[{
                "resolved": host,
                "partition_progress": progress,
                "timestep_ps": 0.004,
                "scale_type": "molecular_dynamics",
                "fallback_total_steps": 100000.0,
                "completion_criteria_type": "number_of_steps",
            }],
        )
        assert result == "00:10:00"

    def test_fused_aws_min_applied_once_not_per_stage(self):
        resource = _aws_resource(timeout=86400)
        policy = structures.Time_policy_adaptive(estimated_performance=1000.0)
        host = structures.Resolved_execution(
            stage_name="equil_stage_1",
            resource_name="aws_gpu",
            resource=resource,
            dispatch=stage_procedures.Dispatch(),
            co_schedule_with=None,
            cpus=4,
            memory_mb=14000,
            time_limit=None,
            mps=1,
            time_policy=policy,
        )
        # Seven tiny equil-like stages: each ~seconds raw; must not become 70 min.
        progress = {
            "null": {
                "completion_criteria": "number_of_steps",
                "attained": False,
                "swarms": {
                    "0": {
                        "completion_criteria": "number_of_steps",
                        "attained": False,
                        "current_step": 0,
                        "total_steps": 20000,  # 0.08 ns at 0.004 ps
                    }
                },
            }
        }
        inputs = [{
            "resolved": host,
            "partition_progress": progress,
            "timestep_ps": 0.004,
            "scale_type": "molecular_dynamics",
            "completion_criteria_type": "number_of_steps",
        } for _ in range(7)]
        result = walltime.estimate_submit_time_limit(
            host_resolved=host,
            stage_inputs=inputs,
        )
        assert result == "00:10:00"

    def test_no_rate_with_remaining_requests_full_max_not_min_floor(self):
        """
        Remaining work without a performance rate must not become a silent
        0s estimate that clamps to the SLURM 1h minimum.
        """
        resource = _slurm_resource(time_limit="48:00:00")
        host = structures.Resolved_execution(
            stage_name="metadynamics_seed",
            resource_name="cluster",
            resource=resource,
            dispatch=stage_procedures.Dispatch(),
            co_schedule_with=None,
            cpus=1,
            memory_mb=4000,
            time_limit=None,
            mps=1,
            time_policy=structures.Time_policy_adaptive(
                estimated_performance=None,
                safety_factor=0.2,
            ),
        )
        progress = {
            "null": {
                "completion_criteria": "number_of_steps",
                "attained": False,
                "current_step": 0,
                "total_steps": 18_000_000,  # 36 ns at 0.002 ps
            }
        }
        assert walltime.estimate_stage_walltime_seconds(
            resolved=host,
            partition_progress=progress,
            timestep_ps=0.002,
            finalize=False,
        ) is None
        result = walltime.estimate_submit_time_limit(
            host_resolved=host,
            stage_inputs=[{
                "resolved": host,
                "partition_progress": progress,
                "timestep_ps": 0.002,
                "scale_type": "molecular_dynamics",
                "completion_criteria_type": "number_of_steps",
            }],
        )
        assert result == "48:00:00"
