Parameterizing the Trypsin-Benzamidine System with Espaloma
===========================================================

This tutorial shows how one would parameterize a system using the 
Espaloma force field.

Prerequisites
-------------
    
- seekrflow installed and working (see :doc:`getting_started`) along with all required, recommended, and optional dependencies, including:
  - seekr
  - Browndye2
  - OpenEye Toolkits
  - PDBFixer
  - PDB2PQR
  - openmmforcefields
  - espaloma
- The espaloma force field ".pt" file downloaded and available somewhere on your system. Download from: https://github.com/choderalab/espaloma/releases/download/0.3.2/espaloma-0.3.2.pt.

Step 1: Espaloma Overview
-------------------------
Espaloma is a force field of a similar functional for as AMBER or CHARMM, 
yet whose valence parameters have been trained on quantum mechanical 
calculations, in many cases, providing a more accurate description of 
molecular interactions. Espaloma uses a graph-convolutional neural network 
to predict bond, angle, and dihedral parameters. Charges are chosen based 
on either the AM1-BCC method, or a neural network trained to reproduce
AM1-BCC charges. This approach will entirely replace parameters for both 
the ligand as well as the protein, although TIP3P will continue to be used 
for the solvent.

Step 2: Parameterizing with Espaloma and Seekrflow
--------------------------------------------------
The only change will be to the arguments to the parameterize.py script, 
(although these changes could also be made at the level of the 
```seekrflow.json``` configuration file). We must point to the location of 
the espaloma force field file.

.. code-block:: bash

    cd seekrflow/seekrflow/examples/trypsin_benzamidine

Open the ``seekrflow_espaloma.json`` file and modify the 
"small_molecule_forcefield" field to point to the location of the espaloma 
force field file you have downloaded.

.. code-block:: bash

    mamba activate SEEKRFLOW_PARAM
    python ~/seekrflow/seekrflow/flow.py parameterize -i seekrflow_espaloma.json

.. note::

    Even though the espaloma force field file is located in the 
    "small_molecule_forcefield" field of the seekrflow.json configuration file,
    it will be used for both the protein and the ligand. The Amber ff14SB 
    entries will only be used to construct a preliminary construct, but will
    be overwritten by the espaloma force field parameters.

For more information about espaloma, see the Github repository at 
https://github.com/choderalab/espaloma.

One may then run the rest of the workflow as before:

.. code-block:: bash

    mamba activate SEEKR
    python ~/seekrflow/seekrflow/flow.py prepare -i seekrflow_espaloma.json
    python ~/seekrflow/seekrflow/flow.py run -i seekrflow_espaloma.json
    python ~/seekr/seekr/analyze.py work_espaloma/root/model.xml