"""
Tests for seekrflow/modules/workload_managers/remote_local.py
"""

import os

from seekrflow.modules import structures
from seekrflow.modules.workload_managers import remote_local


def test_resolve_model_directory_local_shell_uses_work_root(tmp_path):
    work = tmp_path / "work"
    work.mkdir()
    seekrflow = structures.Seekrflow(
        name="host_guest_example",
        work_directory=str(work),
    )
    resource = structures.Resource_direct(
        name="local",
        remote_working_directory="/should/be/ignored",
    )
    resolved = remote_local.resolve_model_directory(seekrflow, resource)
    assert resolved == str(work / "root")
    assert "/should/be/ignored" not in resolved
    assert "host_guest_example" not in os.path.basename(resolved)


def test_resolve_model_directory_ssh_uses_remote_working_directory():
    seekrflow = structures.Seekrflow(name="my_run")
    remote_working_directory = "/scratch/user/seekrflow_work"
    resource = structures.Resource_remote_slurm(
        name="remote_slurm",
        remote_working_directory=remote_working_directory,
        remote_interface=structures.Remote_interface_ssh(
            hostname="login.example.edu"),
        transfer_settings=structures.Transfer_settings_globus(),
    )
    resolved = remote_local.resolve_model_directory(seekrflow, resource)
    assert resolved == os.path.join(
        remote_working_directory,
        "my_run")
