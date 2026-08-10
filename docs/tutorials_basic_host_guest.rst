Host-Guest System: Basic Protein-Ligand Setup
=============================================

This tutorial demonstrates how to run a calculation when the force field 
parameters already exist for a molecular system - in this case, a host-guest 
system.

This host-guest system (where the host is beta-cyclodextrin (BCD) and the 
guest is one of a collection of small molecules like 1-butanol), has been 
parameterized much more optimally than with a generic small-molecular force 
field like AMBER GAFF. We want to use these existing parameters.

Prerequisites
-------------

- seekrflow installed and working (see :doc:`getting_started`) along with all required dependencies, including:
  - seekr
  - Browndye2

Step 1: Prepare and Run The System
----------------------------------

First, find the host-guest example directory:

.. code-block:: bash

    mamba activate SEEKR
    cd ~/seekrflow/seekrflow/examples/host_guest/

In this directory, you will find a file named 
``seekrflow_1_butanol_local.json``. This file contains the
configuration settings for the host-guest system, including the receptor/host 
(BCD) and the guest/ligand (butanol). There is also a directory 
"params_and_structures", which contains many files, including an equilibrated 
PDB starting structure for the BCD/1-butanol system, as well as the force 
field parameters within a AMBER-formatted ``.parm7`` file, as well as PQR 
files to use for BD. These are all defined within the 
``seekrflow_1_butanol_local.json`` file, and one can open this file to see 
how they are featured.

Next, run the parameterization workflow:

.. code-block:: bash

    python ~/seekrflow/seekrflow/flow.py prepare-i seekrflow_1_butanol_local.json
    python ~/seekrflow/seekrflow/flow.py run -i seekrflow_1_butanol_local.json
    python ~/seekr/seekr/analyze.py work_local/root/model.xml

This tutorial gave an example for how to run SEEKR with a system that has already been parameterized.
