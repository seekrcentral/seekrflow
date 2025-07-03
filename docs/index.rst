.. seekrflow documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to seekrflow's documentation!
=========================================================

Seekrflow is a workflow management system for molecular simulations using seekr. It automates the complex process of setting up, running, and analyzing molecular kinetics calculations.

**Key Features:**

- Automated workflow from system parametrization to running simulations
- Seekr stages may be run locally or on remote high-performance computing resources
- Automatically parametrize a molecular system with a standard force field (if you dare!)
- API may be used for batch seekr job creation and submission

.. grid:: 1 1 2 2

    .. grid-item-card:: Getting Started
      :margin: 0 3 0 0
      
      Install seekrflow and get up and running quickly with a basic example.

      .. button-link:: ./getting_started.html
         :color: primary
         :outline:
         :expand:

         To the Getting Started Guide

      

    .. grid-item-card:: Tutorials
      :margin: 0 3 0 0
      
      Step-by-step tutorials covering common workflows and advanced features.

      .. button-link:: ./tutorials.html
         :color: primary
         :outline:
         :expand:

         To the Tutorials

      

    .. grid-item-card::  User Guide
      :margin: 0 3 0 0
      
      Comprehensive guide covering all aspects of using seekrflow.

      .. button-link:: ./user_guide.html
         :color: primary
         :outline:
         :expand:

         To the User Guide
      
      

    .. grid-item-card:: API Reference
      :margin: 0 3 0 0
      
      Complete technical reference for all modules, classes, and functions.

      .. button-link:: ./api.html
         :color: primary
         :outline:
         :expand:

         To the API Reference

      

    .. grid-item-card::  Developer Guide
      :margin: 0 3 0 0
      
      Learn how to contribute to seekrflow and extend its functionality.

      .. button-link:: ./developer_guide.html
         :color: primary
         :outline:
         :expand:

         To the Developer Guide

Quick Start Example
===================

Once seekrflow, SEEKR2, and SeekrTools are installed, a full example workflow 
requires only a PDB files containing a protein and a bound ligand, a configuration
JSON file, and the resname of the ligand. The example can be run with the following commands:

.. code-block:: bash

   cd ~/seekrflow/examples/trypsin_benzamidine/
   python ~/seekrflow/seekrflow/parametrize.py protein_ligand.pdb --input_json seekrflow.json --ligand_resname BEN
   python ~/seekrflow/seekrflow/flow.py seekrflow.json prepare
   python ~/seekrflow/seekrflow/flow.py seekrflow.json run
   python ~/seekr2/seekr2/analyze.py work/root/model.xml

For detailed installation instructions, see the :doc:`getting_started` guide.

.. toctree::
   :maxdepth: 2
   :hidden:
   :titlesonly:

   getting_started
   tutorials
   user_guide
   api
   developer_guide

