# batch_job.py

import os
import seekrflow.modules.base as seekrflow_base
import seekrflow.modules.structures as seekrflow_structures
import seekrflow.flow as seekrflow_flow

# A list of system information. This could, in theory, be read
#  from a data file, a Pandas object, or a spreadsheet.

# ["name", "parm7 filename", "pdb filename", "ligand pqr filename"]
job_list = [
["ligand_1", "BCD_1-butanol.parm7", "BCD_1-butanol.pdb", 
    "BCD_1-butanol_ligand.pqr"],
["ligand_2", "BCD_1-naphthylethanol.parm7", "BCD_1-naphthylethanol.pdb", 
    "BCD_1-naphthylethanol_ligand.pqr"],
["ligand_3", "BCD_1-propanol.parm7", "BCD_1-propanol.pdb", 
    "BCD_1-propanol_ligand.pqr"]]

input_json = "seekrflow_1_butanol.json"

for name, parm7, pdb, pqr in job_list:
    seekrflow = seekrflow_structures.load_seekrflow(input_json)
    seekrflow.name = name
    seekrflow.workflow.solvated_system_for_md.solvated_pdb = os.path.join(
        "../params_and_structures", pdb)
    seekrflow.workflow.solvated_system_for_md.parameters_topology.prmtop_filename\
          = os.path.join(
            "../params_and_structures", parm7)
    seekrflow.work_directory = f"work_{name}"
    seekrflow.workflow.ligand_pqr_filename_for_bd = os.path.join(
        "../params_and_structures", pqr)
    seekrflow.workflow.parameterizer_information.receptor_ligand_pdb_filename\
        = os.path.join(
            "../params_and_structures", pdb)    
    seekrflow.make_work_directory()
    curdir = os.getcwd()
    os.chdir(seekrflow.work_directory)
    src_pdb_filename = seekrflow_flow.assign_seekrflow_defaults(
        seekrflow, "")
    os.chdir(curdir)
    seekrflow.workflow.ligand_indices = seekrflow_base.get_ligand_indices(
        os.path.join("params_and_structures", pdb), 
        seekrflow.workflow.parameterizer_information.ligand_resname)
    seekrflow_flow.flow(seekrflow, "prepare", src_pdb_filename)
    seekrflow_flow.flow(seekrflow, "run")
    
