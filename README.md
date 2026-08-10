seekrflow
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/seekrcentral/seekrflow/workflows/CI/badge.svg)](https://github.com/seekrcentral/seerkflow/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/seekrcentral/seekrflow/branch/master/graph/badge.svg)](https://codecov.io/gh/seekrcentral/seekrflow/branch/master)


[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) 
[![python](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/)

## Overview
Automate seekr handling with seekrflow to decrease time and manual intervention to complete seekr calculations.

The seekr (simulation-enabled estimation of kinetic rates) program facilitates the preparation, running, and
analysis of molecular kinetics calculations using a multiscale molecular dynamics (MD) / Brownian dynamics (BD) 
simulation framework along with a combination of enhanced sampling methods such as milestoning, steered molecular dynamics (SMD) and metadynamics (MetaD). 
Other types of enhanced sampling methods are possible as well through plugins. Nevertheless, running a seekr calculation can be difficult, time-consuming, and error-prone. Parameterization of a system for MD requires time and expertise. Preparing the input files for a seekr calculation is similarly tedious and time-consuming. Running a SEEKR calculation often requires transferring files to a high-performance-computing (HPC) resource, and then writing and submitting several SLURM or PBS scripts per system in order to efficiently conduct the calculations. In addition, all the previously mentioned steps can become exponentially difficult and time-consuming if several systems or batches of seekr jobs must be run.

The seekrflow program aims to remedy the difficulty of using seekr by establishing commonly-used workflows that facilitate the preparation of a seekr calculation. Parameterization (careful!), preparation, and running of seekr calculations are all automated in seekrflow.

This README is only a quickstart guide to get seekrflow up and running as soon as possible on a local computer. To see more detailed **instructions** and **tutorials**, including how to get seekrflow up and running on a remote HPC machine, please see https://seekrflow.readthedocs.io/en/latest or the docs/ subfolder.

## Install

The easiest, quickest way to install seekrflow is to use Mamba. If you don't already have 
Mamba installed, Download the Miniforge install script and run.

```sh
curl -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh
bash Miniforge3-$(uname)-$(uname -m).sh
```

Once this has been done, set up a new environment:

```sh
mamba create -n SEEKR python=3.12 --yes
```

If you plan to use seekrflow for parameterization, you will probably need a second environment to avoid conflicts.
```sh
mamba create -n SEEKRFLOW_PARAM python=3.12 --yes
```

### Bare Minimum Installation
One may install seekrflow with the bare minimum dependencies, which are required for the most basic functionality

This step installs seekrflow from mamba.
```sh
mamba activate SEEKR
mamba install seekrflow --yes
```

Or one can install it from github.
```sh
git clone https://github.com/seekrcentral/seekrflow.git
cd seekrflow
python -m pip install .
```

But then one must install seekr separate (instructions below). The installation of extra dependencies (below) are recommended for certain capabilities.

### Extra Dependencies
Many of the dependencies for seekrflow will be installed alongside seekrflow, but some must be installed separately, and are installed before seekrflow

#### PDBFixer (recommended; needed for parameterization)
PDBfixer is recommended to install in case one wants to run the parameterization using seekrflow.

```sh
mamba activate SEEKRFLOW_PARAM
mamba install pdbfixer --yes
```

#### OpenEye Toolkits (recommended; possibly needed for parameterization)
The OpenEye Toolkits are used for quantum chemistry-based force field parameterization. The OpenEye toolkits require a valid OpenEye academic license, free for academic users but must be obtained directly from https://www.eyesopen.com/academic-licensing. Seekrflow uses it to generate SDF files from PDB files of small molecules. If you do not plan to use seekrflow for parameterization, if will already have SDF files for your small molecules, or if you plan to use RDKit to generate SDF files, then you do not need to install OpenEye.

```sh
mamba install openeye::openeye-toolkits --yes
```

After obtaining an OpenEye academic license, save the provided oe_license.txt file in a secure location on your computer system. For example, you may place it in:
```sh
/home/USERNAME/licenses/oe_license.txt
```

To ensure that OpenEye toolkits can find the license file at runtime, export the license path by adding the following line to your ~/.bashrc. 

```sh
export OE_LICENSE="/home/USERNAME/licenses/oe_license.txt"
```

Then apply the change within the current terminal session.

```sh
source ~/.bashrc
```

#### OpenMM Forcefields(recommended; needed for parameterization)

```sh
mamba install openmmforcefields --yes
```

#### Espaloma Machine-learned Forcefield (optional)
If you want to parameterize your molecular system with the machine-learned forcefield espaloma, you will need to install it.
```sh
mamba install espaloma=0.3.2 --yes
```

You will also need to download the correct espaloma .pt file, and save it somewhere on your computer system.
```sh
curl -O https://github.com/choderalab/espaloma/releases/download/0.3.2/espaloma-0.3.2.pt
```

#### Browndye (recommended)
SEEKR2 can use Browndye if Brownian dynamics (BD) simulations will be run (necessary for k-on calculations). Please see (https://browndye.ucsd.edu/) for Browndye installation instructions. Some of these steps require sudo privileges (administrator access). If you do not have sudo access, contact your system administrator.

#### seekr (required)
This step installs OpenMM plugin for seekr package if you installed seekrflow from github. (A seekrflow installation from Mamba will automatically install seekr.)
```sh
mamba activate SEEKR
mamba install seekr --yes
```

or install seekr from git.
```sh
cd ~
git clone https://github.com/seekrcentral/seekr.git
cd seekr
python -m pip install .
```

Optionally, one may run the tests for seekr.
```sh
pytest
```

### Testing seekrflow (Optional)
To test seekrflow, run the following command in the seekrflow/ directory:

```
pytest
```

## Example System - Trypsin/Benzamidine

A seekrflow calculation will need a input JSON file to run, as well as a starting PDB file containing a bound complex of the receptor and ligand molecules. An example system can be found in seekrflow/seekrflow/examples/trypsin_benzamidine. Here, the bound complex of the receptor and ligand is found in the file "protein_ligand.pdb", and an example input JSON file can be found in "seekrflow.json". A full workflow run can be done with the following steps. (This example runs all calculations locally, and assumes that you have a GPU available.)

```sh
mamba activate SEEKRFLOW_PARAM
python ~/seekrflow/seekrflow/flow.py parameterize -i seekrflow.json
mamba deactivate
mamba activate SEEKR
python ~/seekrflow/seekrflow/flow.py prepare -i seekrflow.json
python ~/seekrflow/seekrflow/flow.py run -i seekrflow.json
python ~/seekr/seekr/analyze.py work/root/model.json
```

### Important Options and Hints

* In general, seekr's and seekrflow's programs can be run with the '-h' argument to see all available options. Please see https://seekrflow.readthedocs.io/en/latest for a detailed description of programs and options.

## Authors and Contributors

The following people have contributed directly to the coding and validation
efforts of seekrflow (listed an alphabetical order of last name). 
Thanks also to everyone who has helped or will help improve this project by 
providing feedback, bug reports, or other comments.

* Rommie Amaro (principal investigator)
* Anand Ojha (developer)
* Lane Votapka (developer)


### Citing seekrflow and Dependencies

For BibTex files of many of the following citations, please visit: https://seekr2.readthedocs.io/en/latest/citations.html

If you use seekrflow, please cite the following paper:

* Anupam A. Ojha, Lane W. Votapka, Shiksha Dutta, Anson F. Noland, Sonya M. Hanson, Rommie E. Amaro, seekrflow: Towards end-to-end automated simulation pipeline with machine-learned force fields for accelerated drug-target kinetic and thermodynamic predictions, bioRxiv 2025.08.13.669965; doi: https://doi.org/10.1101/2025.08.13.669965

One should also cite SEEKR2 and its dependencies:

* Votapka, L. W.; Stokely, A. M.; Ojha, A. A.; Amaro, R. E. SEEKR2: Versatile Multiscale Milestoning Utilizing the OpenMM Molecular Dynamics Engine. J. Chem. Inf. Mod. 2022 62 (13), 3253-3262. DOI: 10.1021/acs.jcim.2c00501

* Van Der Walt, S., Colbert, S.C. & Varoquaux, G., 2011. The NumPy array: a structure for efficient numerical computation. Computing in Science & Engineering, 13(2), pp.22–30.

* Jones, E. et al., 2001. SciPy: Open source scientific tools for Python.

* Hunter, J.D., 2007. Matplotlib: A 2D graphics environment. Computing in Science & Engineering, 9(3), pp.90–95.

* T.D. Swinburne and D.J. Wales, Defining, Calculating, and Converging Observables of a Kinetic Transition Network, J. Chemical Theory and Computation (2020), https://doi.org/10.1021/acs.jctc.9b01211

You may also optionally cite the following papers related to SEEKR2:

* Ojha, A. A., Votapka L. W., Amaro, R. E. QMrebind: incorporating quantum mechanical force field reparameterization at the ligand binding site for improved drug-target kinetics through milestoning simulations. Chemical Science 14 (45), 13159-13175

* Ojha A. A., Srivastava A., Votapka L. W., and Amaro R. E. Selectivity and Ranking of Tight-Binding JAK-STAT Inhibitors Using Markovian Milestoning with Voronoi Tessellations. Journal of Chemical Information and Modeling 2023 63 (8), 2469-2482. DOI: 10.1021/acs.jcim.2c01589

* Votapka, L. W.; Jagger, B. R.; Heyneman, A. L.; Amaro, R. E. SEEKR: Simulation Enabled Estimation of Kinetic Rates, A Computational Tool to Estimate Molecular Kinetics and Its Application to Trypsin–Benzamidine Binding. J. Phys. Chem. B 2017, 121 (15), 3597–3606. https://doi.org/10.1021/acs.jpcb.6b09388. 

* Jagger, B. R.; Ojha, A. A.; Amaro, R. E. Predicting Ligand Binding Kinetics Using a Markovian Milestoning with Voronoi Tessellations Multiscale Approach. J. Chem. Theory Comput. 2020. https://doi.org/10.1021/acs.jctc.0c00495. 

* Jagger, B. R.; Lee, C. T.; Amaro, R. E. Quantitative Ranking of Ligand Binding Kinetics with a Multiscale Milestoning Simulation Approach. J. Phys. Chem. Lett. 2018, 9 (17), 4941–4948. https://doi.org/10.1021/acs.jpclett.8b02047. 

* Votapka LW, Amaro RE (2015) Multiscale Estimation of Binding Kinetics Using Brownian Dynamics, Molecular Dynamics and Milestoning. PLOS Computational Biology 11(10): e1004381. https://doi.org/10.1371/journal.pcbi.1004381

### Copyright

Copyright (c) 2025, Lane Votapka


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.
