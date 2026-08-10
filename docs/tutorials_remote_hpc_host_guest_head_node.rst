Running Host-Guest on a Remote HPC System on the Head Node
==========================================================

This tutorial walks you through setting up and running our host-guest 
example system on a remote HPC system on the head node.

Prerequisites
-------------

- seekrflow installed and working on a HPC head node(see :doc:`getting_started`) 
  along with all required dependencies, including:
  - seekr
  - Browndye2

Step 1: Prepare and Run The System
----------------------------------

First, find the host-guest example directory:

.. code-block:: bash

    mamba activate SEEKR
    cd ~/seekrflow/seekrflow/examples/host_guest/

In this directory, you will find a file named 
``seekrflow_1_butanol_head_node_slurm.json``. This file contains the
configuration settings for running the host-guest system on the head node 
using SLURM and one can open this file to see how they are featured.

Next, run the parameterization workflow:

.. code-block:: bash

    python ~/seekrflow/seekrflow/flow.py prepare -i seekrflow_1_butanol_head_node_slurm.json
    python ~/seekrflow/seekrflow/flow.py run -i seekrflow_1_butanol_head_node_slurm.json
    python ~/seekr/seekr/analyze.py work_head_node_slurm/root/model.xml

This tutorial gave an example for how to run SEEKR with a system that has already been parameterized on the head node.
