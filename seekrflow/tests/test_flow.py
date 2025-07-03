"""
Test the prepare functions of the flow.py script.
"""

import pytest

import seekrflow.modules.base as base
import seekrflow.flow as flow

def test_prepare_system_xml_tryp_ben(tryp_ben_seekrflow_system_xml_parametrized):
    myflow = tryp_ben_seekrflow_system_xml_parametrized
    instruction = "prepare"
    seekrflow_glob, seekrflow_base = flow.flow(myflow, instruction)
    return

def test_prepare_amber_host_guest(host_guest_seekrflow_amber_parametrized):
    myflow = host_guest_seekrflow_amber_parametrized
    instruction = "prepare"
    seekrflow_glob, seekrflow_base = flow.flow(myflow, instruction)
    return
