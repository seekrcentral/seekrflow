Getting Started
===============

Welcome to seekrflow! This guide will walk you through installing and setting up seekrflow for molecular 
dynamics simulation workflows.

What is seekrflow?
------------------

seekrflow is a workflow management system designed for molecular dynamics simulations using the SEEKR2 
package. It provides:

- Automated parameter generation for molecular systems
- Streamlined workflow execution for binding/unbinding simulations
- Integration with high-performance computing resources
- Support for complex molecular systems including protein-ligand interactions

Installation
------------

The easiest, quickest way to install the seekrflow is to use Mamba. If you don't already have 
Mamba installed, Download the Miniforge install script and run.

.. code-block:: bash
    curl -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh
    bash Miniforge3-$(uname)-$(uname -m).sh


Fill out the prompts as they appear.

Once this has been done, set up a new environment:

.. code-block:: bash

    mamba create -n SEEKR2 python=3.11 --yes

If you plan to use seekrflow for parametrization, you will probably need a second environment to 
avoid conflicts.

.. code-block:: bash

    mamba create -n SEEKRFLOW_PARAM python=3.11 --yes

Dependencies
~~~~~~~~~~~~

Many of the dependencies for seekrflow will be installed alongside seekrflow, but some must be 
installed separately, and are installed before seekrflow

Git (required)
++++++++++++++

Make sure git is installed to clone repositories. If git isn't already installed on your 
computer, run:

.. code-block:: bash

    mamba install git --yes

PDBFixer (recommended; needed for parametrization)
+++++++++++++++++++++++++++++++++++++++++++++++++++

PDBfixer is recommended to install in case one wants to run the parametrization using seekrflow.

.. code-block:: bash

    mamba activate SEEKRFLOW_PARAM
    mamba install pdbfixer --yes

pdb2pqr (recommended; needed for parametrization)
+++++++++++++++++++++++++++++++++++++++++++++++++++

pdb2pqr is used in the parametrization workflow in order to choose protonation states, and
to optionally produce PQR files for BD simulations

.. code-block:: bash

    pip install pdb2pqr

OpenEye Toolkits (recommended; possibly needed for parametrization)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The OpenEye Toolkits are used for quantum chemistry-based force field parameterization. 
The OpenEye toolkits require a valid OpenEye academic license, free for academic users 
but must be obtained directly from https://www.eyesopen.com/academic-licensing. 
Seekrflow uses it to generate SDF files from PDB files of small molecules. If you do not 
plan to use seekrflow for parametrization, or will already have SDF files for your small 
molecules, then you do not need to install OpenEye. But if you wish to install it, then 
follow these steps.

.. code-block:: bash

    mamba install openeye::openeye-toolkits --yes

After obtaining an OpenEye academic license, save the provided oe_license.txt file in a secure 
location on your computer system. For example, you may place it in:

.. code-block:: bash

    /home/USERNAME/licenses/oe_license.txt

To ensure that OpenEye toolkits can find the license file at runtime, export the license path 
by adding the following line to your ~/.bashrc.

.. code-block:: bash

    export OE_LICENSE="/home/USERNAME/licenses/oe_license.txt"

Then apply the change within the current terminal session.

.. code-block:: bash

    source ~/.bashrc

OpenMM Forcefields (recommended; needed for parametrization)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

If you wish to parametrize your molecular system with a common forcefield, such as AMBER FF14SB, GAFF2, 
CHARMM 36, OpenFF's SMIRNOFF, or espaloma, you will need to install the OpenMM forcefields package.

.. code-block:: bash

    mamba install openmmforcefields --yes

Espaloma Machine-learned Forcefield (optional)
+++++++++++++++++++++++++++++++++++++++++++++++

If you want to parametrize your molecular system with the machine-learned forcefield espaloma, 
you will need to install it.

.. code-block:: bash

    mamba install espaloma=0.3.2 --yes

You will also need to download the correct espaloma .pt file, and save it somewhere on your 
computer system.

.. code-block:: bash

    curl -O https://github.com/choderalab/espaloma/releases/download/0.3.2/espaloma-0.3.2.pt

Browndye2 (recommended)
+++++++++++++++++++++++

SEEKR2 can use Browndye2 if Brownian dynamics (BD) simulations will be run (necessary for 
k-on calculations). Please see (https://browndye.ucsd.edu/) for Browndye2 installation 
instructions. Some of these steps require sudo privileges (administrator access). 
If you do not have sudo access, contact your system administrator.

SEEKR2 plugin, SEEKR2, and SeekrTools (required)
+++++++++++++++++++++++++++++++++++++++++++++++++

This step installs OpenMM plugin for SEEKR2 package.

.. code-block:: bash

    mamba activate SEEKR2
    mamba install seekr2_openmm_plugin openmm=8.1 --yes

Run the following command to check if the SEEKR2 OpenMM plugin is correctly installed. If no error message appears, the installation was successful.

.. code-block:: bash

    python -c "import seekr2plugin"

Install SEEKR2.

.. code-block:: bash

    cd ~
    git clone https://github.com/seekrcentral/seekr2.git
    cd seekr2
    python -m pip install .

Optionally, one may run the tests for SEEKR2.

.. code-block:: bash

    pytest

Next, clone and install SeekrTools.

.. code-block:: bash

    cd ~
    git clone https://github.com/seekrcentral/seekrtools.git
    cd seekrtools
    python -m pip install .

Optionally, one may run the tests for SeekrTools.

.. code-block:: bash

    pytest

Install seekrflow
~~~~~~~~~~~~~~~~~

Finally, with the dependencies out of the way, we can install seekrflow.

.. code-block:: bash

    git clone https://github.com/seekrcentral/seekrflow.git
    cd seekrflow
    python -m pip install .

Testing seekrflow (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To test seekrflow, run the following command in the seekrflow/ directory:

.. code-block:: bash

    pytest

Installing Dependencies for Remote Execution on HPC Systems
-----------------------------------------------------------

Install Globus Compute SDK
~~~~~~~~~~~~~~~~~~~~~~~~~~

To use remote execution on high-performance computing (HPC) systems, you will need to install 
the Globus Compute SDK. This allows you to run SEEKR2 jobs on remote resources.

.. important::
    The python versions must be the same (or very similar at least) across the local and remote machines.

On the local machine:

.. code-block:: bash

    pip install globus-compute-sdk --yes

On the remote HPC system, you will need to follow the steps above to create the Mamba
environment (on the head or login node; the environment is assumed to be named "SEEKR2").
Into this environment, install SEEKR, SeekrTools, and seekrflow. Then, make sure it is
activated. (There is probably no need to make a parametrization environment on the remote machine,
as it is presumed that parametrization, if any, will be performed on the local machine.)

.. code-block:: bash

    mamba activate SEEKR2

Install pipx for environmental isolation.

.. code-block:: bash

    python3 -m pip install --user pipx

Once this is completed, install globus-compute-endpoint on the remote machine:

.. code-block:: bash

    python3 -m pipx install globus-compute-endpoint

Configure the endpoint:

.. code-block:: bash

    globus-compute-endpoint configure my_seekr_endpoint

At this point, one should modify the file at ~/.globus-compute/my_seekr_endpoint/config.yaml 
in order to properly make full use of the HSP resource's full capabilities.

Here is an example Globus Compute Endpoint configuration file that I used for the NCSA Delta supercomputer:

.. code-block::

    display_name: null
    engine:
    max_workers_per_node: 1
    provider:
        worker_init: "source $HOME/.bashrc; conda activate SEEKR2; export OPENMM_CUDA_COMPILER=`which nvcc`"
        init_blocks: 1
        max_blocks: 1
        min_blocks: 0
        type: LocalProvider
    type: GlobusComputeEngine

.. note::
    This configuration applies to the head/login node *before* SLURM/PBS job submission. The "parsl"
    Python library will be used to actually submit the jobs to the HPC system using SLURM or PBS.

Start the endpoint.

.. code-block:: bash

    globus-compute-endpoint start my_seekr_endpoint

You will need to authenticate with Globus on a browser to start the endpoint, and enter the 
Authorization code given in the browser into the terminal.

One can see the endpoint, as well as its endpoint ID, and those of any other endpoints, 
by listing them:

.. code-block:: bash

    globus-compute-endpoint list

Copy the endpoint ID from the 'start' or 'list' commands above, and save for future reference.

One can also stop the endpoint if/when desired:

.. code-block:: bash

    globus-compute-endpoint stop my_seekr_endpoint

Parsl must be installed on the remote system (or install seekrflow on the remote system - 
which will automatically install parsl):

.. code-block:: bash

    pip install parsl

Before submitting remote jobs, always check the endpoints to make sure they are started
on the remote machine, and start them if they are not.

.. code-block:: bash

    globus-compute-endpoint list
    globus-compute-endpoint start my_seekr_endpoint

More information about Globus endpoints and the Globus SDK can be found here: 
https://globus-compute.readthedocs.io/en/latest/endpoints/endpoints.html

Install Globus Connect Personal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install the latest version of Globus Connect Personal on your local machine. This allows you to
transfer files to and from the remote HPC system using Globus.

Follow the instructions here: https://www.globus.org/globus-connect-personal

--or-- For a typical Linux/Unix installation:
.. code-block:: bash

    curl -O https://downloads.globus.org/globus-connect-personal/linux/stable/globusconnectpersonal-latest.tgz
    tar xzf globusconnectpersonal-latest.tgz
    cd globusconnectpersonal-x.y.z
    ./globusconnectpersonal

And follow the on-screen instructions to set up your Globus Connect Personal instance.

Find the Collection IDs Using the Globus Web Portal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now, one must find the Collection IDs for the Globus Connect Personal instance on the local machine,
and the Globus Compute Endpoint on the remote HPC system. This is done by logging into the
Globus Web Portal (https://app.globus.org) and navigating to the "Collections" tab.

Find your local Globus Connect Personal instance under "Administered By You", click it's name,
and search the page for "UUID". This is the Collection ID for your local Globus Connect Personal instance.
Save this UUID for your seekrflow configuration input file.

Find your remote Globus Collection ID by backtracking to the "Collections" tab, and searching for your
remote Globus Compute Endpoint. Click on the endpoint name, and search the page for "UUID". 
This is the Collection ID for your remote Globus Compute Endpoint. Save this UUID for your seekrflow 
configuration input file.

Quick Start Example System - Trypsin/Benzamidine
------------------------------------------------

A seekrflow calculation will need a input settings JSON file to run, as well as a starting PDB 
file containing a bound complex of the receptor and ligand molecules. An example system can be 
found in seekrflow/seekrflow/examples/trypsin_benzamidine. Here, the bound complex of the 
receptor and ligand is found in the file "protein_ligand.pdb", and an example input JSON 
file can be found in "seekrflow.json". A full workflow run can be done with the following 
steps.

.. code-block:: bash

    mamba activate SEEKRFLOW_PARAM
    python ~/seekrflow/seekrflow/parametrize.py protein_ligand.pdb --input_json seekrflow.json --ligand_resname BEN
    mamba deactivate
    mamba activate SEEKR2
    python ~/seekrflow/seekrflow/flow.py work/seekrflow.json prepare
    python ~/seekrflow/seekrflow/flow.py work/seekrflow.json run
    python ~/seekr2/seekr2/analyze.py work/root/model.xml

Important Options and Hints
---------------------------

* In general, seekrflow, SEEKR2, and SeekrTools programs can be run with the '-h' argument to see 
all available options. Please see https://seekr2.readthedocs.io/en/latest for a detailed 
description of programs and options.

For a complete tutorial, see the :doc:`tutorials` section.

Troubleshooting
---------------

Getting Help
~~~~~~~~~~~~

If you encounter issues:

1. Check the :doc:`user_guide` for detailed usage instructions
2. Review the :doc:`api` for complete API documentation
3. Look at the example in ``seekrflow/examples/trypsin_benzamidine/``
4. See https://seekr2.readthedocs.io/en/latest for SEEKR2-specific help
5. Submit issues to the project repository

Next Steps
----------

Now that you have seekrflow installed, you can:

- Follow the :doc:`tutorials` for step-by-step examples
- Read the :doc:`user_guide` for detailed usage information
- Explore the :doc:`api` reference for complete documentation
- Check out the :doc:`developer_guide` if you want to contribute
