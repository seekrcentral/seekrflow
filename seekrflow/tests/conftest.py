"""
conftest.py

configurations for SEEKR2 tests
"""

import os
import pytest
import copy

import seekrflow.tests.create_seekrflow as create_seekrflow

TEST_DIRECTORY = os.path.dirname(__file__)

@pytest.fixture(scope="session")
def tryp_ben_seekrflow_amber_unparameterized_persistent(tmpdir_factory):
    """
    Create a model object that is persistent across the tests in this file.
    """
    workdir = tmpdir_factory.mktemp("tryp_ben_seekrflow_amber_unparam")
    tryp_ben_seekrflow_amber_unparam_persistent_obj \
        = create_seekrflow.create_unparameterized_tryp_ben_seekrflow()
    tryp_ben_seekrflow_amber_unparam_persistent_obj.work_directory = str(workdir)
    return tryp_ben_seekrflow_amber_unparam_persistent_obj

@pytest.fixture()
def tryp_ben_seekrflow_amber_unparameterized(tryp_ben_seekrflow_amber_unparameterized_persistent):
    """
    Create a model object that is not persistent. But this at least
    doesn't require us to generate an entirely new model.
    """
    tryp_ben_seekrflow_amber_unparameterized_obj = copy.deepcopy(
        tryp_ben_seekrflow_amber_unparameterized_persistent)
    return tryp_ben_seekrflow_amber_unparameterized_obj

@pytest.fixture(scope="session")
def tryp_ben_seekrflow_system_xml_parameterized_persistent(tmpdir_factory):
    """
    Create a model object that is persistent across the tests in this file.
    """
    workdir = tmpdir_factory.mktemp("tryp_ben_seekrflow_system_xml_param")
    tryp_ben_seekrflow_system_xml_param_persistent_obj \
        = create_seekrflow.create_parameterized_tryp_ben_openmm_xml_seekrflow()
    tryp_ben_seekrflow_system_xml_param_persistent_obj.work_directory = str(workdir)
    return tryp_ben_seekrflow_system_xml_param_persistent_obj

@pytest.fixture()
def tryp_ben_seekrflow_system_xml_parameterized(tryp_ben_seekrflow_system_xml_parameterized_persistent):
    """
    Create a model object that is not persistent. But this at least
    doesn't require us to generate an entirely new model.
    """
    tryp_ben_seekrflow_system_xml_param_obj = copy.deepcopy(
        tryp_ben_seekrflow_system_xml_parameterized_persistent)
    return tryp_ben_seekrflow_system_xml_param_obj

# TODO test espaloma parameterization
@pytest.fixture(scope="session")
def tryp_ben_seekrflow_espaloma_unparameterized_persistent(tmpdir_factory):
    """
    Create a model object that is persistent across the tests in this file.
    """
    workdir = tmpdir_factory.mktemp("tryp_ben_seekrflow_espaloma_unparam")
    espaloma_ptr_path = os.getenv("ESPALOMA_PTR_PATH")
    assert espaloma_ptr_path is not None, \
        "ESPALOMA_PTR_PATH environment variable must be set to run this test."
    tryp_ben_seekrflow_espaloma_unparam_persisent_obj \
        = create_seekrflow.create_unparameterized_tryp_ben_seekrflow(ff=espaloma_ptr_path)
    tryp_ben_seekrflow_espaloma_unparam_persisent_obj.work_directory = str(workdir)
    return tryp_ben_seekrflow_espaloma_unparam_persisent_obj

@pytest.fixture()
def tryp_ben_seekrflow_espaloma_unparameterized(tryp_ben_seekrflow_espaloma_unparameterized_persistent):
    """
    Create a model object that is not persistent. But this at least
    doesn't require us to generate an entirely new model.
    """
    tryp_ben_seekrflow_espaloma_unparam_obj = copy.deepcopy(
        tryp_ben_seekrflow_espaloma_unparameterized_persistent)
    return tryp_ben_seekrflow_espaloma_unparam_obj

@pytest.fixture(scope="session")
def host_guest_seekrflow_amber_parameterized_persistent(tmpdir_factory):
    """
    Create a model object that is persistent across the tests in this file.
    """
    workdir = tmpdir_factory.mktemp("host_guest_seekrflow_amber_param")
    host_guest_seekrflow_amber_parameterized_persistent_obj \
        = create_seekrflow.create_parameterized_host_guest_amber_seekrflow()
    host_guest_seekrflow_amber_parameterized_persistent_obj.work_directory = str(workdir)
    host_guest_seekrflow_amber_parameterized_persistent_obj.ligand_resname = "APN"
    host_guest_seekrflow_amber_parameterized_persistent_obj.receptor_selection = "resname MGO"
    return host_guest_seekrflow_amber_parameterized_persistent_obj

@pytest.fixture()
def host_guest_seekrflow_amber_parameterized(host_guest_seekrflow_amber_parameterized_persistent):
    """
    Create a model object that is not persistent. But this at least
    doesn't require us to generate an entirely new model.
    """
    host_guest_seekrflow_amber_param_obj = copy.deepcopy(
        host_guest_seekrflow_amber_parameterized_persistent)
    return host_guest_seekrflow_amber_param_obj