"""
Tests for PBS/Torque workload manager support.

This module tests PBS resource configuration, serialization, and basic workflow functions.
"""

import cattrs
import pytest
import json
from seekrflow.modules.structures import (
    Resource_remote_pbs,
    Run_settings
)
from cattrs.strategies import include_subclasses


def test_resource_pbs_creation():
    """Test that PBS resource can be created with default values"""
    resource = Resource_remote_pbs(
        name="test_pbs",
        queue="batch",
        account="test_account",
    )

    assert resource.type == "pbs_remote"
    assert resource.name == "test_pbs"
    assert resource.queue == "batch"
    assert resource.account == "test_account"
    assert resource.cpus_per_task == 1  # default
    assert resource.memory_per_node == 4  # default
    assert resource.time_limit == "00:30:00"  # default


def test_resource_pbs_with_custom_values():
    """Test PBS resource with custom values"""
    resource = Resource_remote_pbs(
        name="hpc_cluster",
        queue="gpu",
        account="MY_ALLOCATION",
        resource_list="qos=high",
        cpus_per_task=8,
        memory_per_node=64000,
        time_limit="04:00:00",
        scheduler_options="-m abe -M user@email.com",
        worker_init="module load python; conda activate SEEKR2",
    )

    assert resource.queue == "gpu"
    assert resource.resource_list == "qos=high"
    assert resource.cpus_per_task == 8
    assert resource.memory_per_node == 64000
    assert resource.time_limit == "04:00:00"
    assert resource.scheduler_options == "-m abe -M user@email.com"
    assert resource.worker_init == "module load python; conda activate SEEKR2"


def test_resource_pbs_serialization():
    """Test that PBS resource can be serialized to dict"""
    resource = Resource_remote_pbs(
        name="test_pbs",
        queue="batch",
        account="test_account",
        cpus_per_task=4,
        memory_per_node=8000,
        time_limit="01:00:00",
    )

    # Create converter with subclass registration
    from seekrflow.modules.structures import Resource_base
    import cattrs
    test_converter = cattrs.Converter()
    include_subclasses(Resource_base, test_converter)

    resource_dict = test_converter.unstructure(resource)

    assert resource_dict["type"] == "pbs_remote"
    assert resource_dict["name"] == "test_pbs"
    assert resource_dict["queue"] == "batch"
    assert resource_dict["account"] == "test_account"
    assert resource_dict["cpus_per_task"] == 4
    assert resource_dict["memory_per_node"] == 8000
    assert resource_dict["time_limit"] == "01:00:00"


def test_resource_pbs_deserialization():
    """Test that PBS resource can be deserialized from dict"""
    resource_dict = {
        "type": "pbs_remote",
        "name": "test_pbs",
        "remote_seekr2_directory": "$HOME/seekr2/seekr2/",
        "remote_seekrtools_directory": "$HOME/seekrtools/seekrtools/",
        "remote_working_directory": "/scratch/user/seekrflow",
        "queue": "batch",
        "account": "test_account",
        "resource_list": None,
        "nodes_per_block": 1,
        "cpus_per_task": 4,
        "memory_per_node": 8000,
        "time_limit": "01:00:00",
        "scheduler_options": "",
        "worker_init": "",
        "remote_interface": {
            "type": "globus_compute_sdk",
            "endpoint_id": "",
        },
        "transfer_settings": {
            "type": "globus",
            "source_endpoint_id": "",
            "destination_endpoint_id": "",
        },
    }

    # Create converter with subclass registration
    from seekrflow.modules.structures import Resource_base
    import cattrs
    test_converter = cattrs.Converter()
    include_subclasses(Resource_base, test_converter)

    resource = test_converter.structure(resource_dict, Resource_remote_pbs)

    assert isinstance(resource, Resource_remote_pbs)
    assert resource.type == "pbs_remote"
    assert resource.name == "test_pbs"
    assert resource.queue == "batch"
    assert resource.account == "test_account"
    assert resource.cpus_per_task == 4
    assert resource.memory_per_node == 8000
    assert resource.time_limit == "01:00:00"


def test_resource_pbs_json_round_trip():
    """Test full JSON serialization/deserialization cycle"""
    original = Resource_remote_pbs(
        name="cluster1",
        queue="gpu",
        account="allocation123",
        resource_list="qos=high",
        cpus_per_task=16,
        memory_per_node=128000,
        time_limit="08:00:00",
    )

    # Create converter with subclass registration
    from seekrflow.modules.structures import Resource_base
    import cattrs
    test_converter = cattrs.Converter()
    include_subclasses(Resource_base, test_converter)

    # Serialize to dict then JSON
    resource_dict = test_converter.unstructure(original)
    json_str = json.dumps(resource_dict)

    # Deserialize back
    loaded_dict = json.loads(json_str)
    loaded = test_converter.structure(loaded_dict, Resource_remote_pbs)

    # Verify all key fields match
    assert loaded.name == original.name
    assert loaded.queue == original.queue
    assert loaded.account == original.account
    assert loaded.resource_list == original.resource_list
    assert loaded.cpus_per_task == original.cpus_per_task
    assert loaded.memory_per_node == original.memory_per_node
    assert loaded.time_limit == original.time_limit


def test_run_settings_with_pbs_resource():
    """Test that Run_settings can include PBS resources"""
    pbs_resource = Resource_remote_pbs(
        name="pbs_cluster",
        queue="batch",
        account="test",
    )

    run_settings = Run_settings(
        resources=[pbs_resource],
        bd_stage_resource_name="local",
        hidr_stage_resource_name="pbs_cluster",
        seekr_stage_resource_name="pbs_cluster",
    )

    assert len(run_settings.resources) == 1
    assert isinstance(run_settings.resources[0], Resource_remote_pbs)
    assert run_settings.resources[0].name == "pbs_cluster"


def test_run_settings_mixed_resources():
    """Test Run_settings with both SLURM and PBS resources"""
    from seekrflow.modules.structures import Resource_remote_slurm

    slurm_resource = Resource_remote_slurm(
        name="slurm_cluster",
        partition="gpu",
        account="test1",
    )

    pbs_resource = Resource_remote_pbs(
        name="pbs_cluster",
        queue="batch",
        account="test2",
    )

    run_settings = Run_settings(
        resources=[slurm_resource, pbs_resource],
        bd_stage_resource_name="local",
        hidr_stage_resource_name="slurm_cluster",
        seekr_stage_resource_name="pbs_cluster",
    )

    assert len(run_settings.resources) == 2
    assert isinstance(run_settings.resources[0], Resource_remote_slurm)
    assert isinstance(run_settings.resources[1], Resource_remote_pbs)


def test_get_resource_by_name_pbs():
    """Test getting PBS resource by name from Run_settings"""
    resource = Resource_remote_pbs(
        name="my_pbs",
        queue="batch",
        account="test",
    )

    run_settings = Run_settings(
        resources=[resource],
    )

    retrieved = run_settings.get_resource_by_name("my_pbs")
    assert retrieved is not None
    assert isinstance(retrieved, Resource_remote_pbs)
    assert retrieved.name == "my_pbs"


@pytest.mark.needs_pbs
def test_pbs_workflow_functions_available():
    """Test that PBS workflow functions can be imported"""
    from seekrflow.modules.workload_managers import pbs

    # Check that all required functions exist
    assert hasattr(pbs, "calculate_optimal_seekr_time_limit")
    assert hasattr(pbs, "pbs_remote_run_workflow")
    assert hasattr(pbs, "pbs_remote_bd_status_workflow")
    assert hasattr(pbs, "pbs_remote_hidr_status_workflow")
    assert hasattr(pbs, "pbs_remote_seekr_status_workflow")
    assert hasattr(pbs, "pbs_remote_cancel_workflow")
    assert hasattr(pbs, "pbs_remote_force_rerun_workflow")

@pytest.mark.needs_pbs
def test_pbs_time_calculation_no_file():
    """Test PBS time calculation falls back to default when file missing"""
    from seekrflow.modules.workload_managers.pbs import calculate_optimal_seekr_time_limit

    result = calculate_optimal_seekr_time_limit(
        job_status_file="/nonexistent/path.json",
        incomplete_anchors=[0, 1, 2],
        default_time_limit="02:00:00"
    )

    # Should return default when file doesn't exist
    assert result == "02:00:00"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
