Trypsin-Benzamidine System: Basic Protein-Ligand Setup
======================================================

This tutorial walks you through setting up a basic protein-ligand binding 
simulation using the trypsin-benzamidine example system.

Prerequisites
-------------

- seekrflow installed and working (see :doc:`getting_started`) along with all required, recommended, and optional dependencies, including:

  - seekr
  
  - Browndye2
  
  - OpenEye Toolkits
  
  - PDBFixer
  
  - PDB2PQR
  
  - openmmforcefields
  
- Basic understanding of molecular dynamics simulations and seekr (https://github.com/seekrcentral/seekr)

Step 1: Prepare Your System
---------------------------

First, let's examine the example system provided with seekrflow:

.. code-block:: bash

    # Navigate to the example directory
    cd seekrflow/examples/trypsin_benzamidine/
    
    # List the contents
    ls -la

You should see several files including:

- ``protein_ligand.pdb`` - An initial structure containing a protein receptor and a bound small molecule ligand
- ``seekrflow.json`` - The main input settings file in JSON format

Step 2: Understanding the Configuration and Structure
-----------------------------------------------------

One may open the ``seekrflow.json`` file to see the configuration settings. 
Inside, there are several key settings that one would define based on their 
preferences, such as the *name* of the project, the *workflow* being 
employed, as well as physical parameters or simulation settings, like 
*temperature*, *nonbonded cutoff*, etc.

In this tutorial, we will be using seekrflow to parameterize a system for 
which only the structure is known. Therefore, the *parameterizer* section 
of the configuration file has been filled out. We also have some settings 
for outside programs like *pdb2pqr* and *PDBFixer* which are used for 
configuration.

Also, take a look at the ``protein_ligand.pdb`` file. This file contains 
the initial structure of the protein and ligand. One may view it in one's 
favorite standard molecular viewer such as VMD, PyMOL, Chimera, Maestro, 
or NGLView.

Step 3: Run the Parameterization
--------------------------------

Enter the example trypsin-benzamidine directory and execute the 
parameterization workflow

.. code-block:: bash

    cd seekrflow/seekrflow/examples/trypsin_benzamidine
    mamba activate SEEKRFLOW_PARAM
    python ~/seekrflow/seekrflow/flow.py parameterize -i seekrflow.json

Within a few minutes, a directory named "work" will be created, which will 
contain a directory "parameterize", along with some other files. Within this 
directory, the forcefield parameters, specified within the configuration JSON 
file, will be used to generate a system for MD simulations, as well as 
creating PQR files for the ligand and protein for a BD simulation. Using these 
defaults, the protein will be parameterized with AMBER FF14SB, the ligand with 
AMBER GAFF2.11, and the system will be solvated with waters and ions from TIP3P,
and the Joung and Cheatham recommended salt models and divalent counterion
parameters.

.. note::
    The force fields supported by seekrflow are whatever are provided 
    by openmmforcefields (https://github.com/openmm/openmmforcefields).
    This includes all major AMBER force fields, CHARMM non-polarizable
    protein, nucleic acid, and pre-parameterized small-molecule force fields.
    The espaloma force field may also be used, and is covered in a different
    tutorial.

.. note::
    Small molecule force fields include GAFF1.81 to GAFF2.2.20, as well as 
    OpenFF force fields up to openff-2.x.y "Sage".

Once the parameterization is complete, one ought to check the system carefully 
to ensure that parameters were assigned correctly. The parameterization 
script will automatically perform minimization and a short equilibration to 
check for system stability. However, one should also load their solvated 
structure into a molecular viewer to ensure that the charges are assigned
reasonably, and that the solvent doesn't contain any artifacts, for instance. 
One should ideally also run additional equilibrations to monitor the system 
stability, and to ensure that the system is ready for production runs.

.. caution::

    This parameterization feature in seekrflow is convenient, but relies 
    on experimental tools such as OpenFF, and could potentially introduce 
    incorrect parameters. Ideally, one should carefully produce their own 
    parameterized system by hand, using tools such as AMBER's LEAP, 
    CharmmGUI, or OpenFF's tools step-by-step. However, if you're feeling 
    adventurous, the seekrflow automated parameterization feature is here 
    for your convenience. The parameterization feature in seekrflow is also 
    not likely to be able to automatically handle covalent metal parameters 
    correctly at present.

The files are now created in work/parameterize/.

Step 4: Run the SEEKR stages
----------------------------

Once the parameterization is complete, you can proceed to run the rest of 
the workflow and prepare and run the seekr calculation:

.. code-block:: bash

    mamba activate SEEKR
    python ~/seekrflow/seekrflow/flow.py prepare -i seekrflow.json
    python ~/seekrflow/seekrflow/flow.py run -i seekrflow.json

.. warning::

    These settings were designed to allow one to quickly and easily run the 
    seekrflow workflow, but they would require substantial modifications 
    to be used for an accurate calculation on this system, or any other 
    biomedically-relevant system. For instance, the step count is way too 
    low at 1,000,000 steps (2 ns) per anchor. In the original seekr 
    publication involving this system, we ran each anchor for 250,000,000 
    steps (500 ns). For one's own seekr and seekflow calculations, one must 
    carefully and thoughtfully choose all settings and validate the 
    correctness and optimality of all settings, force  field parameters, 
    and starting structures.

For those who are curious, the experimentally-measured k-off and k-on of 
trypsin-benzamidine binding are (6 ± 3) * 10^+02 and (2.4 ± 0.2) * 10^+7 
M^-1 s^-1, respectively.