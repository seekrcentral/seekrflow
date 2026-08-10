Using the seekrflow API to Run a Batch Job
==========================================

Sometimes, one may have a set of SEEKR calculations to perform, and one would 
like to run them all in a batch job. One could, of course, apply the 
principles from the previous tutorials to each system one at a time. But if 
one wants to script seekrflow in some way to, for instance, run a set of 
host-guest systems in a batch job, one can use the seekrflow API - which this 
tutorial will demonstrate.

Prerequisites
-------------

- seekrflow installed and working (see :doc:`getting_started`) along with all required dependencies, including:
  - seekr
  - Browndye2

Step 1: Make the Batch Job Python Script
----------------------------------------

One could use a Python script such as the following to run a batch job of SEEKR calculations
using the seekrflow API:

.. code-block:: python

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
    
A batch method such as this, using the API, could enable a high degree of automation within
one's SEEKR calculations. 
