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

input_json = "seekrflow.json"

for name, parm7, pdb, pqr in job_list:
    seekrflow = seekrflow_structures.load_seekrflow(input_json)
    seekrflow.name = name
    seekrflow.receptor_ligand_pdb = os.path.join(
        "../params_and_structures", pdb)
    seekrflow.work_directory = f"work_{name}"
    seekrflow.bd_settings.ligand_pqr_filename = os.path.join(
        "../params_and_structures", pqr)
    seekrflow.md_parameters_topology.prmtop_filename = os.path.join(
        "../params_and_structures", parm7)
    seekrflow.starting_pdb_filename = os.path.join(
        "../params_and_structures", pdb)
        
    seekrflow.make_work_directory()
    seekrflow.ligand_indices = seekrflow_base.get_ligand_indices(
        os.path.join("params_and_structures", pdb), seekrflow.ligand_resname)
    seekrflow_flow.flow(seekrflow, "prepare")
    seekrflow_flow.flow(seekrflow, "run")
    
