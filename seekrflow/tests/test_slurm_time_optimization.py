"""
Test the SLURM time limit optimization functionality.
"""

import json
import os
import tempfile
import pytest

import seekrflow.modules.workload_managers.slurm as workload_slurm


@pytest.fixture
def mock_status_dir():
    """Create temporary directory for mock status files"""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


def create_job_status(
    tmpdir,
    elapsed="00:15:00",
    anchor_times=None,
    target_time=100.0,
    incomplete_anchors=None,
    has_manager_status=True
):
    """Helper to create mock .seekrflow_job_status.json file"""
    job_status_file = os.path.join(tmpdir, ".seekrflow_job_status.json")

    anchor_times = anchor_times or {"0": 0.0, "1": 0.0}
    incomplete_anchors = incomplete_anchors or [0, 1]

    manager_status = {
        "jobs": [{"Elapsed": elapsed, "State": "COMPLETED"}]
    } if has_manager_status else None

    job_status = {
        "seekr": {
            "manager_status": manager_status,
            "stage_status": {
                "anchor_times": anchor_times,
                "anchor_time_for_completion": target_time,
                "incomplete_anchors": incomplete_anchors
            }
        }
    }

    with open(job_status_file, "w") as f:
        json.dump(job_status, f)

    return job_status_file


def create_state_file(tmpdir, anchor_times_at_submission):
    """Helper to create mock .slurm_runner/seekr_latest.json state file"""
    slurm_runner_dir = os.path.join(tmpdir, ".slurm_runner")
    os.makedirs(slurm_runner_dir, exist_ok=True)

    state_file = os.path.join(slurm_runner_dir, "seekr_latest.json")
    state_data = {
        "jobid": "12345",
        "anchor_times_at_submission": anchor_times_at_submission
    }

    with open(state_file, "w") as f:
        json.dump(state_data, f)

    return state_file


def test_normal_case_without_baseline(mock_status_dir):
    """Test with previous run data, no baseline (uses current anchor times as progress)"""
    job_status_file = create_job_status(
        mock_status_dir,
        elapsed="00:15:00",  # 900 seconds
        anchor_times={"0": 100.0, "1": 100.0, "2": 50.0, "3": 25.0},
        incomplete_anchors=[2, 3]
    )

    # Progress: 50 + 25 = 75 ps in 900s → 12 s/ps
    # Anchor 3 needs 75 ps → 900s → with 20% buffer → 1080s → rounds to 1 hour
    result = workload_slurm.calculate_optimal_seekr_time_limit(
        job_status_file, [2, 3], "48:00:00"
    )
    assert result == "01:00:00"


def test_with_baseline_state_file(mock_status_dir):
    """Test with anchor_times_at_submission baseline for accurate progress calculation"""
    job_status_file = create_job_status(
        mock_status_dir,
        elapsed="00:30:00",  # 1800 seconds
        anchor_times={"0": 100.0, "1": 100.0, "2": 60.0, "3": 40.0},
        incomplete_anchors=[2, 3]
    )

    create_state_file(
        mock_status_dir,
        anchor_times_at_submission={"0": 100.0, "1": 100.0, "2": 40.0, "3": 10.0}
    )

    # Progress: (60-40) + (40-10) = 50 ps in 1800s → 36 s/ps
    # Anchor 3 needs 60 ps → 2160s → with buffer → 2592s → rounds to 1 hour
    result = workload_slurm.calculate_optimal_seekr_time_limit(
        job_status_file, [2, 3], "48:00:00"
    )
    assert result == "01:00:00"


def test_cap_at_max_time_limit(mock_status_dir):
    """Test that calculated time never exceeds default_time_limit"""
    job_status_file = create_job_status(
        mock_status_dir,
        elapsed="02:00:00",
        anchor_times={"0": 5.0, "1": 0.0},
        incomplete_anchors=[0, 1]
    )

    result = workload_slurm.calculate_optimal_seekr_time_limit(
        job_status_file, [0, 1], "03:00:00"
    )
    assert result == "03:00:00"


def test_minimum_one_hour(mock_status_dir):
    """Test that minimum time limit is enforced even for nearly complete jobs"""
    job_status_file = create_job_status(
        mock_status_dir,
        elapsed="00:10:00",
        anchor_times={"0": 90.0, "1": 95.0},
        incomplete_anchors=[0, 1]
    )

    result = workload_slurm.calculate_optimal_seekr_time_limit(
        job_status_file, [0, 1], "48:00:00"
    )
    assert result == "01:00:00"


def test_multiday_elapsed_time(mock_status_dir):
    """Test parsing of SLURM's multi-day elapsed time format (D-HH:MM:SS)"""
    job_status_file = create_job_status(
        mock_status_dir,
        elapsed="2-04:30:00",  # 2 days, 4h, 30m = 187800 seconds
        anchor_times={"0": 50.0, "1": 50.0},
        incomplete_anchors=[0, 1]
    )

    # Progress: 100 ps in 187800s → 1878 s/ps
    # Need 50 ps → 93900s → with buffer → 112680s → 31.3h → rounds to 32 hours
    result = workload_slurm.calculate_optimal_seekr_time_limit(
        job_status_file, [0, 1], "48:00:00"
    )
    assert result == "32:00:00"


@pytest.mark.parametrize("condition,expected", [
    ("no_manager_status", "24:00:00"),  # First run, no benchmark data
    ("no_progress", "24:00:00"),        # Previous run made no progress
    ("file_not_exist", "24:00:00"),     # File doesn't exist
])
def test_fallback_cases(mock_status_dir, condition, expected):
    """Test various fallback scenarios that should return default time limit"""
    if condition == "no_manager_status":
        job_status_file = create_job_status(
            mock_status_dir, has_manager_status=False
        )
    elif condition == "no_progress":
        job_status_file = create_job_status(
            mock_status_dir,
            elapsed="00:30:00",
            anchor_times={"0": 0.0, "1": 0.0, "2": 0.0},
            incomplete_anchors=[0, 1, 2]
        )
    else:  # file_not_exist
        job_status_file = "/nonexistent/path/.seekrflow_job_status.json"

    result = workload_slurm.calculate_optimal_seekr_time_limit(
        job_status_file, [0, 1], "24:00:00"
    )
    assert result == expected
