"""Legacy SLURM calculator is a deprecated no-op; see test_walltime.py."""

import seekrflow.modules.workload_managers.slurm as workload_slurm


def test_deprecated_calculator_returns_default(tmp_path):
    status = tmp_path / ".seekrflow_job_status.json"
    status.write_text("{}")
    assert workload_slurm.calculate_optimal_seekr_time_limit(
        str(status), [0, 1], "12:00:00"
    ) == "12:00:00"
