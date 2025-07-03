Tutorials
=========

This section provides step-by-step tutorials to help you get started with seekrflow. 

.. contents::
   :local:
   :depth: 2

Tutorial 1: Basic Protein-Ligand Setup
--------------------------------------

This tutorial walks you through setting up a basic protein-ligand binding simulation using 
the trypsin-benzamidine example system.

Prerequisites
~~~~~~~~~~~~~

- seekrflow installed and working (see :doc:`getting_started`) along with all required, recommended, and optional dependencies, including:
  - SEEKR2
  - SeekrTools
  - Browndye2
  - OpenEye Toolkits
  - PDBFixer
  - PDB2PQR
  - openmmforcefields
- Basic understanding of molecular dynamics simulations and SEEKR (https://github.com/seekrcentral/seekr2)

Step 1.1: Prepare Your System
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, let's examine the example system provided with seekrflow:

.. code-block:: bash

    # Navigate to the example directory
    cd seekrflow/examples/trypsin_benzamidine/
    
    # List the contents
    ls -la

You should see several files including:

- ``protein_ligand.pdb`` - An initial structure containing a protein receptor and a bound small molecule ligand
- ``seekrflow.json`` - The main input settings file in JSON format

Step 1.2: Understanding the Configuration and Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One may open the ``seekrflow.json`` file to see the configuration settings. Inside,
there are several key settings that one would define based on their preferences, such
as the *name* of the project, the type of *workflow* being employed, as well as physical
parameters or simulation settings, like *temperature*, *nonbonded cutoff*, etc.

In this tutorial, we will be using seekrflow to parametrize a system for which only the
structure is known. Therefore, the *parametrizer* section of the configuration file
has been filled out. We also have some settings for outside programs like *pdb2pqr*
and *PDBFixer* which are used for configuration.

Also, take a look at the ``protein_ligand.pdb`` file. This file contains the initial structure
of the protein and ligand. One may view it in your favorite standard molecular viewer such as
VMD, PyMOL, Chimera, Maestro, or NGLView.

Step 1.3: Run the Parametrization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Enter the example trypsin-benzamidine directory and execute the parametrization workflow

.. code-block:: bash

    cd seekrflow/seekrflow/examples/trypsin_benzamidine
    mamba activate SEEKRFLOW_PARAM
    python ~/seekrflow/seekrflow/parametrize.py protein_ligand.pdb --input_json seekrflow.json

Within a few minutes, a directory named "work" will be created, which will contain a directory
"parametrize", along with some other files. Within this directory, the forcefield parameters
specified within the configuration JSON file will be used to generate a system for MD simulations,
as well as creating PQR files for the ligand and protein for a BD simulation. Using these defaults,
the protein will be parametrized with AMBER FF14SB, the ligand with AMBER GAFF2.11, and the 
system will be solvated with waters and ions from TIP3P.

Once the parametrization is complete, one ought to check the system carefully to ensure that
parameters were assigned correctly. The parametrization script will automatically perform
minimization and a short equilibration to check for system stability. However, one should also
load their solvated structure into a molecular viewer to ensure that the charges are assigned
reasonably, and that the solvent doesn't contain any artifacts, for instance. One should ideally
also run additional equilibrations to monitor the system stability, and to ensure that
the system is ready for production runs.

.. caution::

    This parametrization feature in seekrflow is convenient, but relies on experimental tools
    such as OpenFF, and could potentially introduce incorrect parameters (hence why it is
    placed in a separate script from the main workflow). Ideally, one should
    carefully produce their own parametrized system by hand, using tools such as AMBER's LEAP,
    CharmmGUI, or OpenFF's tools step-by-step. However, if you're feeling adventurous, the
    seekrflow automated parametrization feature is here for your convenience.

An additional file has also been created inside work/ - a new copy of the ``seekrflow.json`` file
with the parametrization settings filled in. This file is used to run the SEEKR calculation, and 
should be used from now on - not the original ``seekrflow.json`` file outside of the work/ directory.

Step 1.4: Run the SEEKR stages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once the parametrization is complete, you can proceed to run the rest of the workflow
and prepare and run the SEEKR calculation:

.. code-block:: bash

    mamba activate SEEKR2
    python ~/seekrflow/seekrflow/flow.py work/seekrflow.json prepare
    python ~/seekrflow/seekrflow/flow.py work/seekrflow.json run
    python ~/seekr2/seekr2/analyze.py work/root/model.xml

Normally, within SEEKR, one would need to define the ligand atom indices, as well as the
atom indices defining the binding site. The seekrflow "protein_ligand" workflow automates 
this process to save us time, by using a residue name of the ligand (assumed to be unique).
Many other setup steps are also automated based on a best-practices workflow for proteins
and small molecules. Other types of workflows, such as "protein_protein", or perhaps involving
membranes or nucleic acids, will be developed in the future with their own best practices.

The previous commands will take a while (probably about an hour) to run, and when it is complete, 
you should see a (very coarse) approximation of the k-off and k-on of trypsin-benzamidine binding.

.. warning::

    These settings were designed to allow one to quickly and easily run the seekrflow
    workflow, but they would require substantial modifications to be used for
    an accurate calculation on this system, or any other biomedically-relevant system.
    For instance, the step count is way too low at 1,000,000 steps (2 ns) per anchor. In the 
    original SEEKR publication involving this system, we ran each anchor for 250,000,000 steps 
    (500 ns). For one's own SEEKR and seekflow calculations, one must carefully and thoughtfully
    choose all settings and validate the correctness and optimality of all settings, force 
    field parameters, and starting structures.

Tutorial 2: Parametrizing the System with Espaloma
--------------------------------------------------

This tutorial shows how one would parametrize a system using the Espaloma force field.

Prerequisites
~~~~~~~~~~~~~

- seekrflow installed and working (see :doc:`getting_started`) along with all required, recommended, and optional dependencies, including:
  - SEEKR2
  - SeekrTools
  - Browndye2
  - OpenEye Toolkits
  - PDBFixer
  - PDB2PQR
  - openmmforcefields
  - espaloma
- The espaloma force field ".pt" file downloaded and available somewhere on your system. Download from: https://github.com/choderalab/espaloma/releases/download/0.3.2/espaloma-0.3.2.pt.

Step 2.1: Espaloma Overview
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Espaloma is a force field of a similar functional for as AMBER or CHARMM, yet whose 
valence parameters have been trained on quantum mechanical calculations, in many cases,
providing a more accurate description of molecular interactions. Espaloma uses a graph
convolutional neural network to predict bond, angle, and dihedral parameters. Charges
are chosen based on either the AM1-BCC method, or a neural network trained to reproduce
AM1-BCC charges. This approach will entirely replace parameters for both the ligand as well
as the protein, although TIP3P will continue to be used for the solvent.

Step 2.2: Parametrizing with Espaloma and Seekrflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The only change will be to the arguments to the parametrize.py script, (although these changes
could also be made at the level of the ```seekrflow.json``` configuration file). We must
point to the location of the espaloma force field file, and also, specify a new working directory.

. code-block:: bash

    cd seekrflow/seekrflow/examples/trypsin_benzamidine
    mamba activate SEEKRFLOW_PARAM
    python ~/seekrflow/seekrflow/parametrize.py protein_ligand.pdb --input_json seekrflow.json --external_ff_file /path/to/espaloma-0.3.2.pt --work_directory work_espaloma

For more information about espaloma, see the Github repository at https://github.com/choderalab/espaloma.

One may then run the rest of the workflow as before, using the new configuration file
``work_espaloma/seekrflow.json``:

Tutorial 3: Running on a Remote HPC System
------------------------------------------

This tutorial walks you through setting up and running our trypsin-benzamidine example
on a remote HPC system using Globus and Parsl

Prerequisites
~~~~~~~~~~~~~

- seekrflow installed and working (see :doc:`getting_started`) along with all required, recommended, and optional dependencies, including:
  - SEEKR2
  - SeekrTools
  - Browndye2
  - OpenEye Toolkits
  - PDBFixer
  - PDB2PQR
  - openmmforcefields
  - Globus Endpoints
  - Globus Compute SDK
- Access to a remote HPC system, where you can submit jobs and manage resources.

Step 3.1: Justification
~~~~~~~~~~~~~~~~~~~~~~~
Full SEEKR calculations almost always require a power GPU cluster or supercomputer, yet transferring
files to and from a remote system, as well as managing job submissions with SLURM/PBS scripts can be 
slow and cumbersome. This tutorial shows how to use seekrflow to streamline this process, although
care must be taken to ensure that the remote system is configured correctly, and that the
Globus endpoints are set up properly. So make sure that all dependences are installed on both the
local and remote machines, and that the remote machine is set up as defined in :doc:`getting_started`.

Step 3.2: Prepare the Configuration JSON Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A couple of example configuration files are provided in the directory, named ```seekrflow_delta.json``` 
and ``seekrflow_anvil.json``. If one opens these files, one will see a large section filled out
titled *run_settings*. This section contains *Parsl* and *Globus* settings, which are used to
transfer files to and from the remote system, as well as to submit jobs to the remote system. One 
will need to configure these settings to their own system - they will be different for everyone, and
it's impossible for me to anticipate the changes you will make, so you will need to be proactive and
resourceful in order to get this tutorial completed. Yet, if you can complete this tutorial, you
should be all set for running your own systems on HPC. Let us consider some of the settings that
one will likely need to modify in order to complete this tutorial.

- "type": This should be set to either "slurm_remote" or "pbs_remote", depending on the job scheduler 
  used by your HPC system. NOTE: at this time only "slurm_remote" is supported.

- "name": Choose any name for your resource, which will be referenced by the "_resource_name" fields
  lower in the configuration file.

- "remote_seekr2_directory": As the name suggests, enter the path to SEEKR2 on your remote system.

- "remote_seekrtools_directory": As the name suggests, enter the path to SEEKRTools on your remote system.

- "remote_working_directory": This is the directory on the remote system where the SEEKR workflow will be 
  copied into and run. Typically, HPC resources have a "scratch", "work", or "projects" directory where 
  intensive read/write operations can be performed. Make sure that you have write permissions to this 
  directory.

- "max_workers_per_node": This will probably always be set to 1. However, there might be some conceivable
  situations where one might want more than one Parsl worker per node. Consult the Parsl documentation
  to explore other possible settings for this parameter.

- "partition": This is the partition on the remote system where the jobs will be submitted. This is often 
  something like "gpu" or "compute". Check with your HPC documentation to find the correct 
  partition name.

- "account": This is the account name that you were assigned for job submissions on the remote system. 
  You should check with any online portal or HPC documentation to find the correct account name.

- "nodes_per_block": This is the number of nodes that will be requested per Parsl "block", and will
  probably usually be kept at 1. Consult the Parsl documentation to explore other possible settings 
  for this parameter.

- "cores_per_node": This is the number of cores that will be requested per node. This should be set to
  the number of cores that you would like to use for each job. Note that seekrflow is designed to 
  request shared resources, so this should not exceed the proportional number of cores that you
  would like to use for a shared job (using a single GPU, for instance). Consult the Parsl documentation to explore other possible settings 
  for this parameter, as well as your HPC documentation to find the correct number of cores
  to request for your jobs.

- "memory_per_node": This is the amount of memory that will be requested per node. This should be set to
  the amount of memory that you would like to use for each job. Note that seekrflow is designed to 
  request shared resources, so this should not exceed the proportional amount of memory that you
  would like to use for a shared job (using a single GPU, for instance). Consult the Parsl documentation to explore other possible settings 
  for this parameter, as well as your HPC documentation to find the correct amount of memory
  to request for your jobs.

- "time_limit": This is the maximum amount of time that the job will be allowed to run on the remote system.
  This should be set to a reasonable amount of time for your jobs, and should be set according to your HPC
  documentation. Note that this is not the same as the total simulation time, but rather the maximum time
  that the job will be allowed to run before it is killed. Example: "time_limit": "24:00:00" would be 24 
  hours.

- "scheduler_options": These are settings that will be passed to the job scheduler when submitting jobs. 
  This can include things like job names, output files, error files, etc. Most importantly, this line
  will probably be used to assign GPU settings. Consult your HPC documentation to find the correct 
  settings for your system.

- "worker_init": These settings define which commands will be run upon the creation of a new Parsl "worker".
  This might include the loading of important modules, setting environment variables, or
  activating a conda/mamba environment. This will depend on which HPC resource is being used. One
  should consult the HPC documentation, and probably experiment with debug/test job submissions in order
  to find the correct settings for their system.

- "globus_compute_endpoint_id": This is the Globus Compute endpoint ID that will be used to submit jobs 
  to the remote system. This should be set to the endpoint ID that you created for your remote system. 
  You can find this ID by running ```globus compute endpoint list``` in the terminal of the remote HCP
  resource.

- "transfer_settings":
  
  - "type": This should be set to "globus" to use Globus for file transfers. At present, only
    Globus transfers are supported.
  - "local_collection_id": The Globus collection (endpoint) ID for the local system. One can find
     it by getting globus_connect_personal running on one's own machine, and then using the Globus
     web portal "Collections" page to find its UUID.
  - "remote_collection_id": The Globus collection (endpoint) ID for the remote system. One can find
     it by searching for the Globus collection UUID of the HPC resource in the Globus web portal 
     "Collections" page.

- "bd_stage_resource_name", "hidr_stage_resource_name", and "seekr_stage_resource_name":
  These are the names of the resources that will be used to stage the Browndye, Hidr, and SEEKR 
  calculations, respectively. These should match the names defined in the "resources" section of the 
  configuration file. Note that multiple resources can be defined and used, including just a local
  computer.

- "allow_parsl_usage_tracking": If set to True, Parsl will collect usage statistics and send them 
  to the Parsl team. This is optional, but helps improve the library.

Step 3.3: Run Seekrflow on a Remote HPC System
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Once these settings are configured, one can run the seekrflow workflow on the remote HPC system
in the same way as before:

.. code-block:: bash

    mamba activate SEEKR2
    python ~/seekrflow/seekrflow/flow.py work/seekrflow.json prepare
    python ~/seekrflow/seekrflow/flow.py work/seekrflow.json run
    python ~/seekr2/seekr2/analyze.py work/root/model.xml

The job will probably take quite a long time to run, depending on the resources available
on the remote system, as well as the backlog in the remote job queue. However, the BD 
simulations should still be run remotely and synchronously with the rest of the jobs.
In this configuration, HIDR will be run on the remote system first, and then the SEEKR 
anchor calculations will be run synchronously with each other. All file transfers should
be automatically handled to and from the remote resource.


Tutorial 4: Host-Guest System: Existing Force Field Parameters  
--------------------------------------------------------------

This tutorial shows how to run a calculation when the force field parameters already 
exist for a molecular system - in this case, a host-guest system.

This host-guest system (where the host is beta-cyclodextrin (BCD) and the guest is one of a
collection of small molecules like 1-butanol), has been parametrized much more optimally than with a
generic small-molecular force field like AMBER GAFF. We want to use these existing
parameters, not anything we would make with parametrize.py in seekrflow.

Prerequisites
~~~~~~~~~~~~~

- seekrflow installed and working (see :doc:`getting_started`) along with all required dependencies, including:
  - SEEKR2
  - SeekrTools
  - Browndye2

Step 4.1: Prepare and Run The System
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, find the host-guest example directory:

.. code-block:: bash

    mamba activate SEEKR2
    cd ~/seekrflow/seekrflow/examples/host_guest/

In this directory, you will find a file named ``seekrflow_1_butanol.json``. This file contains the
configuration settings for the host-guest system, including the receptor/host (BCD) and the 
guest/ligand (butanol). There is also a directory "params_and_structures", which contains many files,
including a PDB starting structure for the BCD/1-butanol system, as well as the force field parameters
within a AMBER-formatted ``.parm7`` file, as well as PQR files to use for BD. These are all
defined within the ``seekrflow_1_butanol.json`` file, and one can open this file to see how they 
are featured.

Next, run the parametrization workflow:

.. code-block:: bash

    python ~/seekrflow/seekrflow/flow.py seekrflow_1_butanol.json prepare
    python ~/seekrflow/seekrflow/flow.py work_butanol/seekrflow.json run
    python ~/seekr2/seekr2/analyze.py work_butanol/root/model.xml

This tutorial gave an example for how to run SEEKR with a system that has already been parametrized.


Tutorial 5: Host-Guest System: Set up a Batch Job with the seekrflow API  
------------------------------------------------------------------------

Sometimes, one may have a set of SEEKR calculations to perform, and one would like to
run them all in a batch job. One could, of course, apply the principles from the previous
tutorials to each system one at a time. But if one wants to script seekrflow in some way to,
for instance, run a set of host-guest systems in a batch job, one can use the seekrflow API - 
which this tutorial will demonstrate.

Prerequisites
~~~~~~~~~~~~~

- seekrflow installed and working (see :doc:`getting_started`) along with all required dependencies, including:
  - SEEKR2
  - SeekrTools
  - Browndye2

Step 5.1: Make the Batch Job Python Script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
    
A batch method such as this, using the API, could enable a high degree of automation within
one's SEEKR calculations. 

Next Steps
----------

After completing these tutorials, you should be able to:

- Set up basic protein-ligand binding simulations
- Customize simulation parameters for your systems
- Manage complex workflows that may use combinations of local or remote resources
- Apply seekrflow to your own systems of interest

For more advanced topics, see:

- :doc:`user_guide` - Comprehensive usage documentation
- :doc:`api` - Complete API reference
- :doc:`developer_guide` - Contributing to seekrflow

