"""Tests for local_shell remote interface and none transfer."""
from __future__ import annotations

import json
import os
import sys

import cattrs
import pytest
from cattrs.strategies import include_subclasses

from seekrflow.modules import structures
from seekrflow.modules.remote_interfaces import local_shell
from seekrflow.modules.transfer import base as transfer_base
from seekrflow.modules.workload_managers import remote as workload_remote


def trivial_local_shell_workflow(args):
    """Must be top-level so inspect.getsource works under python -c."""
    return {"ok": True, "n": args[0]}


def failing_local_shell_workflow(args):
    raise RuntimeError("boom")


def test_submit_remote_workflow_with_local_shell_round_trip():
    result = local_shell.submit_remote_workflow_with_local_shell(
        "test",
        trivial_local_shell_workflow,
        [42],
        python_executable=sys.executable,
    )
    assert result == {"ok": True, "n": 42}


def test_submit_remote_workflow_with_local_shell_nonzero_exit():
    with pytest.raises(RuntimeError, match="exited with code"):
        local_shell.submit_remote_workflow_with_local_shell(
            "fail",
            failing_local_shell_workflow,
            [],
            python_executable=sys.executable,
        )


def test_resource_slurm_local_shell_and_none_json_round_trip():
    original = structures.Resource_remote_slurm(
        name="login_slurm",
        remote_working_directory="/scratch/user/seekrflow_work",
        remote_interface=structures.Remote_interface_local_shell(
            python_executable="python3",
        ),
        transfer_settings=structures.Transfer_settings_none(),
    )

    converter = cattrs.Converter()
    include_subclasses(structures.Transfer_settings_base, converter)
    include_subclasses(structures.Resource_base, converter)

    resource_dict = converter.unstructure(original)
    loaded = converter.structure(
        json.loads(json.dumps(resource_dict)),
        structures.Resource_remote_slurm,
    )

    assert loaded.name == original.name
    assert loaded.remote_working_directory == original.remote_working_directory
    assert isinstance(loaded.remote_interface, structures.Remote_interface_local_shell)
    assert loaded.remote_interface.type == "local_shell"
    assert loaded.remote_interface.python_executable == "python3"
    assert isinstance(loaded.transfer_settings, structures.Transfer_settings_none)
    assert loaded.transfer_settings.type == "none"


def test_transfer_none_is_noop(tmp_path):
    resource = structures.Resource_remote_slurm(
        name="login_slurm",
        remote_working_directory=str(tmp_path),
        remote_interface=structures.Remote_interface_local_shell(),
        transfer_settings=structures.Transfer_settings_none(),
    )
    # Should return without raising even if directories do not exist remotely.
    transfer_base.transfer_files_to_from_remote_resource(
        "root", resource, str(tmp_path), backwards=False)
    transfer_base.transfer_files_to_from_remote_resource(
        "root", resource, str(tmp_path), backwards=True)


def test_resolve_remote_model_directory_local_shell_uses_root(tmp_path):
    work = tmp_path / "work"
    work.mkdir()
    seekrflow = structures.Seekrflow(
        name="host_guest_example",
        work_directory=str(work),
    )
    resource = structures.Resource_remote_slurm(
        name="login_slurm",
        remote_working_directory="/should/be/ignored",
        remote_interface=structures.Remote_interface_local_shell(),
        transfer_settings=structures.Transfer_settings_none(),
    )
    resolved = workload_remote.resolve_remote_model_directory(
        seekrflow, resource)
    assert resolved == str(seekrflow.get_root_directory())
    assert resolved == str(work / "root")
    assert "/should/be/ignored" not in resolved
    assert "host_guest_example" not in os.path.basename(resolved)


def test_resolve_remote_model_directory_ssh_uses_remote_wd():
    seekrflow = structures.Seekrflow(name="my_run")
    resource = structures.Resource_remote_slurm(
        name="remote_slurm",
        remote_working_directory="/scratch/user/seekrflow_work",
        remote_interface=structures.Remote_interface_ssh(hostname="login.example.edu"),
        transfer_settings=structures.Transfer_settings_globus(),
    )
    resolved = workload_remote.resolve_remote_model_directory(
        seekrflow, resource)
    assert resolved == os.path.join(
        "/scratch/user/seekrflow_work", "my_run")


def test_local_shell_requires_transfer_none_slurm():
    with pytest.raises(ValueError, match="must be 'none'"):
        structures.Resource_remote_slurm(
            name="login_slurm",
            remote_interface=structures.Remote_interface_local_shell(),
            transfer_settings=structures.Transfer_settings_globus(),
        )


def test_local_shell_requires_transfer_none_pbs():
    with pytest.raises(ValueError, match="must be 'none'"):
        structures.Resource_remote_pbs(
            name="login_pbs",
            remote_interface=structures.Remote_interface_local_shell(),
            transfer_settings=structures.Transfer_settings_rsync(),
        )
