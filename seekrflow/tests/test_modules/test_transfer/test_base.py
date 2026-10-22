"""
Tests for seekrflow/modules/transfer/base.py
"""

from seekrflow.modules import structures
from seekrflow.modules.transfer import base as transfer_base


def test_transfer_none_leaves_local_files(tmp_path):
    local = tmp_path / "local"
    local.mkdir()
    marker = local / "marker.txt"
    marker.write_text("keep")
    resource = structures.Resource_direct(
        name="local",
        remote_working_directory=str(tmp_path / "remote"),
    )
    transfer_base.transfer_files_to_from_remote_resource(
        "root", resource, str(local), backwards=False)
    transfer_base.transfer_files_to_from_remote_resource(
        "root", resource, str(local), backwards=True)
    assert marker.read_text() == "keep"
