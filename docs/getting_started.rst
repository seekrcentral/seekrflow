Getting Started
===============

Welcome to seekrflow! This guide will walk you through installing and setting up seekrflow for molecular 
dynamics simulation workflows.

What is seekrflow?
------------------

seekrflow is a workflow management system designed for molecular dynamics simulations using the seekr 
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

If you plan to use seekrflow for parameterization, you will probably need a second environment to 
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

PDBFixer (recommended; needed for parameterization)
+++++++++++++++++++++++++++++++++++++++++++++++++++

PDBfixer is recommended to install in case one wants to run the parameterization using seekrflow.

.. code-block:: bash

    mamba activate SEEKRFLOW_PARAM
    mamba install pdbfixer --yes

pdb2pqr (recommended; needed for parameterization)
++++++++++++++++++++++++++++++++++++++++++++++++++

pdb2pqr is used in the parameterization workflow in order to choose protonation states, and
to optionally produce PQR files for BD simulations

.. code-block:: bash

    pip install pdb2pqr

OpenEye Toolkits (recommended; possibly needed for parameterization)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The OpenEye Toolkits are used for quantum chemistry-based force field parameterization. 
The OpenEye toolkits require a valid OpenEye academic license, free for academic users 
but must be obtained directly from https://www.eyesopen.com/academic-licensing. 
Seekrflow uses it to generate SDF files from PDB files of small molecules. If you do not 
plan to use seekrflow for parameterization, or will already have SDF files for your small 
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

OpenMM Forcefields (recommended; needed for parameterization)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

If you wish to parameterize your molecular system with a common forcefield, such as AMBER FF14SB, GAFF2, 
CHARMM 36, OpenFF's SMIRNOFF, or espaloma, you will need to install the OpenMM forcefields package.

.. code-block:: bash

    mamba install openmmforcefields --yes

Espaloma Machine-learned Forcefield (optional)
++++++++++++++++++++++++++++++++++++++++++++++

If you want to parameterize your molecular system with the machine-learned forcefield espaloma, 
you will need to install it.

.. code-block:: bash

    mamba install espaloma=0.3.2 --yes

You will also need to download the correct espaloma .pt file, and save it somewhere on your 
computer system.

.. code-block:: bash

    curl -O https://github.com/choderalab/espaloma/releases/download/0.3.2/espaloma-0.3.2.pt

Browndye2 (recommended)
+++++++++++++++++++++++

Seekr can use Browndye2 if Brownian dynamics (BD) simulations will be run (necessary for 
k-on calculations). Please see (https://browndye.ucsd.edu/) for Browndye2 installation 
instructions. Some of these steps require sudo privileges (administrator access). 
If you do not have sudo access, contact your system administrator.

SEEKR2 plugin, SEEKR2, and SeekrTools (required)
++++++++++++++++++++++++++++++++++++++++++++++++

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To test seekrflow, run the following command in the seekrflow/ directory:

.. code-block:: bash

    pytest

Installing Dependencies for Remote Execution on HPC Systems
-----------------------------------------------------------

If one want to use seekrflow to run SEEKR on a remote system, one must decide whether
to use SSH or Globus Compute SDK in order to control the jobs. They both function
very similarly, but here are the pros and cons of each.

.. list-table:: SSH vs Globus Compute SDK Comparison
   :widths: 20 20 20
   :header-rows: 1

   * - Remote System Approach
     - Pros
     - Cons
   * - SSH
     - Simple, supported on most Linux systems
     - Requires password or key-based authentication, passwords in plain text, cannot bypass 2-factor authentication
   * - Globus Compute SDK
     - Only an endpoint id is needed, avoiding password and 2-factor authentication.
     - More complicated setup of software and settings

Install Globus Compute SDK (optional - if not using SSH)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To use remote execution on high-performance computing (HPC) systems with Globus Compute SDK, 
you will need to install the package. This allows you to run seekr 
jobs on remote resources while avoiding annoying things like passwords and 2-factor 
authentication.

.. important::
    The python versions must be the same (or very similar at least) across the local and remote machines.

On the local machine:

.. code-block:: bash

    pip install globus-compute-sdk

On the remote HPC system, you will need to follow the steps above to create the Mamba
environment (on the head or login node; the environment is assumed to be named "SEEKR2").
Into this environment, install SEEKR, SeekrTools, and seekrflow. Then, make sure it is
activated. (There is probably no need to make a parameterization environment on the remote machine,
as it is presumed that parameterization, if any, will be performed on the local machine.)

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
in order to properly make full use of the HPC resource's full capabilities.

Here is an example Globus Compute Endpoint configuration file that I used for the NCSA Delta supercomputer:

.. code-block::

    display_name: null
    engine:
    max_workers_per_node: 2
    provider:
        worker_init: "source /path/to/your/.bashrc; conda activate SEEKR2; export OPENMM_CUDA_COMPILER=`which nvcc`"
        init_blocks: 1
        max_blocks: 1
        min_blocks: 0
        type: LocalProvider
    type: GlobusComputeEngine

.. note::
    This configuration applies to the head/login node *before* SLURM/PBS job submission. 
    Notice that the **max_workers_per_node** is set to **2** for proper handling of the
    Globus compute SDK client in seekrflow.

.. note::
    It may be a better idea to define the OPENMM_CUDA_COMPILER environment variable in your
    shell initialization file (e.g., ~/.bashrc) rather than in the worker_init command

Start the endpoint.

.. code-block:: bash

    globus-compute-endpoint start my_seekr_endpoint

You will need to authenticate with Globus on a browser to start the endpoint, and enter the 
Authorization code given in the browser into the terminal. This authentication session 
should last for a number of days or weeks, depending on your Globus account settings.

One can see the endpoint, as well as its endpoint ID, and those of any other endpoints, 
by listing them:

.. code-block:: bash

    globus-compute-endpoint list

Copy the endpoint ID from the 'start' or 'list' commands above, and save for future reference.

One can also stop the endpoint if/when desired:

.. code-block:: bash

    globus-compute-endpoint stop my_seekr_endpoint

Before submitting remote jobs, always check the endpoints to make sure they are started
on the remote machine, and start them if they are not.

.. code-block:: bash

    globus-compute-endpoint list
    globus-compute-endpoint start my_seekr_endpoint

More information about Globus endpoints and the Globus SDK can be found here: 
https://globus-compute.readthedocs.io/en/latest/endpoints/endpoints.html

Install Globus Connect Personal (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

File transfers may be done with Globus or with rsync. Some resources prefer the use
of Globus for large file transfers, instead of rsync.

If you want to use Globus to transfer files to and from your remote machine, install 
the latest version of Globus Connect Personal on your local machine. This allows you to
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
    python ~/seekrflow/seekrflow/parameterize.py --input_json seekrflow.json --ligand_resname BEN
    mamba deactivate
    mamba activate SEEKR2
    python ~/seekrflow/seekrflow/flow.py prepare --input_json seekrflow.json
    python ~/seekrflow/seekrflow/flow.py run --input_json seekrflow.json
    python ~/seekr2/seekr2/analyze.py work/root/model.xml

Set up your .seekrflow_resources.json (Optional)
------------------------------------------------

One can save a lot of time preparing to submit remote jobs if one sets up their `.seekrflow_resources.json`
file. Instead of always including one's run_settings.resources list within the seekrflow.json files,
one may define all one's resources in the `.seekrflow_resources.json` file. The file may exist either
in one's home directory (recommended), within one's seekrflow work/ directory, or within one's 
seekrflow work/root/ directory. The file is structured as follows:

.. code-block:: bash

    {
    "resources": [
    {
    "type": "slurm_remote",
    "name": "delta",
    "remote_seekr2_directory": "/PATH/TO/seekr2/seekr2/",
    "remote_seekrtools_directory": "/PATH/TO/seekrtools/seekrtools/",
    "remote_working_directory": "/PATH/TO/seekrflow_playground/",
    "max_workers_per_node": 1,
    "partition": "gpuA100x4,gpuA40x4",
    "account": "YOUR ACCOUNT NAME",
    "constraint": "scratch",
    "nodes_per_block": 1,
    "cpus_per_task": 32,
    "memory_per_node": 5000,
    "time_limit": "00:30:00",
    "scheduler_options": "--gpus-per-node=1 --gpu-bind=closest",
    "worker_init": "source $HOME/.bashrc; conda activate SEEKR2; export OPENMM_CUDA_COMPILER=`which nvcc`",
    "remote_interface": {
        "type": "globus_compute_sdk",
        "endpoint_id": "YOUR GLOBUS COMPUTE ENDPOINT"
    },
    "transfer_settings": {
        "type": "globus",
        "local_collection_id": "YOUR LOCAL COLLECTION ID",
        "remote_collection_id": "7e936164-de58-4e3d-85da-21aa23c07169"
    },
    { 
    "type": "slurm_remote",
    "name": "anvil",
    "remote_seekr2_directory": "/PATH/TO/seekr2/seekr2/",
    "remote_seekrtools_directory": "/PATH/TO/seekrtools/seekrtools/",
    "remote_working_directory": "/PATH/TO/seekrflow_playground",
    "max_workers_per_node": 1,
    ...
    }
    ]
    }
    
As one utilizes resources to run seekrflow, one may populate this file with various resources
that one may use to run seekr remotely.

Important Options and Hints
---------------------------

* In general, seekrflow, seekr, and SeekrTools programs can be run with the '-h' argument to see 
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
