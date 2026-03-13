"""
Test the prepare functions of the flow.py script.
"""

import seekrflow.flow as flow

def test_prepare_system_xml_tryp_ben(tryp_ben_seekrflow_system_xml_parameterized):
    myflow = tryp_ben_seekrflow_system_xml_parameterized
    instruction = "prepare"
    src_pdb_filename = myflow.workflow.solvated_system_for_md.solvated_pdb
    flow.flow(myflow, instruction, src_pdb_filename)
    return

def test_prepare_amber_host_guest(host_guest_seekrflow_amber_parameterized):
    myflow = host_guest_seekrflow_amber_parameterized
    instruction = "prepare"
    src_pdb_filename = myflow.workflow.solvated_system_for_md.solvated_pdb
    flow.flow(myflow, instruction, src_pdb_filename)
    return
