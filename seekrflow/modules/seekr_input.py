"""
modules/seekr_input.py

Using the seekr API, create an input object and generate the model.
"""

import os
import pathlib
from shutil import copyfile

import parmed
import seekr2.modules.common_base as seekr2_base
import seekr2.modules.common_prepare as seekr2_common_prepare
import seekr2.modules.common_cv as seekr2_common_cv
import seekr2.modules.check as seekr2_check
import seekr2.prepare as seekr2_prepare

import seekrflow.modules.structures as structures
import seekrflow.modules.cvs as cvs

APBS_GRID_SPACING = 0.5

# NOTE: These are valid for SEEKR2 - will need to change for SEEKR3
def assign_amber_params(input_anchor, prmtop_filename, pdb_filename):
    input_anchor.starting_amber_params = seekr2_base.Amber_params()
    input_anchor.starting_amber_params.prmtop_filename = prmtop_filename
    input_anchor.starting_amber_params.pdb_coordinates_filename = pdb_filename
    return

def assign_forcefield_params(input_anchor, built_in_ff_list, custom_ff_list, 
                             pdb_filename):
    input_anchor.starting_forcefield_params = seekr2_base.Forcefield_params()
    input_anchor.starting_forcefield_params.built_in_forcefield_filenames \
        = built_in_ff_list
    input_anchor.starting_forcefield_params.custom_forcefield_filenames \
        = custom_ff_list
    input_anchor.starting_forcefield_params.pdb_coordinates_filename \
        = pdb_filename
    return

def assign_system_params(input_anchor, system_filename, pdb_filename):
    input_anchor.starting_forcefield_params = seekr2_base.Forcefield_params()
    input_anchor.starting_forcefield_params.system_filename \
        = system_filename
    input_anchor.starting_forcefield_params.pdb_coordinates_filename \
        = pdb_filename
    return

def assign_charmm_params(input_anchor, psf_filename, charmm_ff_filenames, 
                         pdb_filename):
    input_anchor.starting_charmm_params = seekr2_base.Charmm_params()
    input_anchor.starting_charmm_params.psf_filename = psf_filename
    input_anchor.starting_charmm_params.charmm_ff_files = charmm_ff_filenames
    input_anchor.starting_charmm_params.pdb_coordinates_filename = pdb_filename
    return

def create_mmvt_receptor_ligand_com_com_model_input_seekr2(
        seekrflow: structures.Seekrflow,
        alpha_carbon_ligand_threshold: float = 0.6
        ) -> seekr2_common_prepare.Model_input:
    """
    Create the input object for the MMVT model in SEEKR2.
    """
    # Copy starting PDB file over to the prepare directory
    root_directory = seekrflow.get_root_directory()
    model_input = seekr2_common_prepare.Model_input()
    model_input.calculation_type = "mmvt"
    model_input.calculation_settings = seekr2_common_prepare.MMVT_input_settings()
    model_input.calculation_settings.md_output_interval = seekrflow.seekr_settings.md_output_interval
    model_input.calculation_settings.md_steps_per_anchor = seekrflow.seekr_settings.md_steps_per_anchor
    model_input.temperature = seekrflow.temperature
    if seekrflow.pressure is None:
        model_input.pressure = 1.0
        model_input.ensemble = "nvt"
    else:
        model_input.pressure = seekrflow.pressure
        model_input.ensemble = "npt"
    model_input.root_directory = root_directory
    model_input.md_program = "openmm"
    model_input.constraints = "HBonds"
    model_input.rigidWater = True
    if seekrflow.hmass == 1.008:
        model_input.hydrogenMass = None
    else:
        model_input.hydrogenMass = seekrflow.hmass
    model_input.timestep = seekrflow.stepsize
    model_input.nonbonded_cutoff = seekrflow.nonbonded_cutoff

    assert seekrflow.starting_pdb_filename != "", \
        "Starting PDB filename is not set in the input JSON file."
    starting_pdb_full_path = os.path.join(
        seekrflow.work_directory, seekrflow.starting_pdb_filename)
    cv_input1 = seekr2_common_cv.Spherical_cv_input()
    receptor_atoms, ligand_atoms = cvs.get_receptor_ligand_com_com_selections(
        seekrflow, starting_pdb_full_path, alpha_carbon_ligand_threshold)
    cv_input1.group1 = receptor_atoms
    cv_input1.group2 = ligand_atoms
    cv_input1.input_anchors = []
    radius_list = seekrflow.seekr_settings.anchor_radius_list
    assert len(radius_list) > 0, \
        "Anchor radius list is empty in the input JSON file."
    assert seekrflow.md_parameters_topology is not None, \
        "md_parameters_topology is not set in the input JSON file."
    for i, radius in enumerate(radius_list):
        input_anchor = seekr2_common_cv.Spherical_cv_anchor()
        input_anchor.radius = radius
        if seekrflow.md_parameters_topology.type == "Amber":
            amber_prmtop_filename = seekrflow.md_parameters_topology.prmtop_filename
            assign_amber_params(input_anchor, amber_prmtop_filename, "")
        elif seekrflow.md_parameters_topology.type == "Gromacs":
            raise NotImplementedError(
                "Gromacs topology type is not implemented yet.")
        elif seekrflow.md_parameters_topology.type == "Charmm":
            charmm_psf_filename = seekrflow.md_parameters_topology.psf_filename
            charmm_ff_filenames = seekrflow.md_parameters_topology.param_filename_list
            assign_charmm_params(input_anchor, charmm_psf_filename, 
                                 charmm_ff_filenames, "")
        elif seekrflow.md_parameters_topology.type == "OpenMM_forcefield":
            built_in_ff_list = seekrflow.md_parameters_topology.built_in_forcefield_filenames
            custom_ff_list = seekrflow.md_parameters_topology.custom_forcefield_filenames
            assign_forcefield_params(input_anchor, built_in_ff_list, 
                                     custom_ff_list, "")
        elif seekrflow.md_parameters_topology.type == "OpenMM_system":
            system_filename = seekrflow.md_parameters_topology.system_filename
            assign_system_params(input_anchor, system_filename, "")
        else:
            raise Exception("Type not supported for seekrflow.md_parameters_topology: "\
                            f"{seekrflow.md_parameters_topology.type}")

        if i == 0:
            input_anchor.bound_state = True
        else:
            input_anchor.bound_state = False
            
        if i == len(radius_list)-1:
            input_anchor.bulk_anchor = True
        else:
            input_anchor.bulk_anchor = False
    
        cv_input1.input_anchors.append(input_anchor)
    
    model_input.cv_inputs = [cv_input1]
    if seekrflow.bd_settings is not None:
        # TODO: fill out from parametrize
        bd_rec_indices, bd_lig_indices = cvs.get_bd_receptor_ligand_selections(
            seekrflow, starting_pdb_full_path, receptor_atoms, ligand_atoms)
        model_input.browndye_settings_input \
            = seekr2_common_prepare.Browndye_settings_input()
        model_input.browndye_settings_input.binary_directory = seekrflow.bd_settings.binary_directory
        model_input.browndye_settings_input.receptor_pqr_filename \
            = seekrflow.bd_settings.receptor_pqr_filename
        model_input.browndye_settings_input.ligand_pqr_filename \
            = seekrflow.bd_settings.ligand_pqr_filename
        model_input.browndye_settings_input.apbs_grid_spacing = APBS_GRID_SPACING
        model_input.browndye_settings_input.receptor_indices = bd_rec_indices
        model_input.browndye_settings_input.ligand_indices = bd_lig_indices
        
        ion1 = seekr2_base.Ion()
        ion1.radius = 1.2
        ion1.charge = -1.0
        ion1.conc = seekrflow.ionic_strength
        ion2 = seekr2_base.Ion()
        ion2.radius = 0.9
        ion2.charge = 1.0
        ion2.conc = seekrflow.ionic_strength
        model_input.browndye_settings_input.ions = [ion1, ion2]
        model_input.browndye_settings_input.num_b_surface_trajectories = seekrflow.bd_settings.num_trajectories
        model_input.browndye_settings_input.n_threads = 1
    else:
        model_input.browndye_settings_input = None

    return model_input


def prepare_model(
        seekrflow: structures.Seekrflow,
        ) -> None:
    """
    Construct and prepare the model for seekr simulation.
    """
    seekrflow.work_directory = os.path.abspath(seekrflow.work_directory)
    curdir = os.getcwd()
    os.chdir(seekrflow.work_directory)
    if seekrflow.workflow_type ==  "protein_ligand":
        # Create the input object for the MMVT model
        model_input = create_mmvt_receptor_ligand_com_com_model_input_seekr2(
            seekrflow, alpha_carbon_ligand_threshold=0.6)

    model, xml_path = seekr2_prepare.prepare(model_input, force_overwrite=True)
    if model.anchor_rootdir == ".":
        model_dir = os.path.dirname(xml_path)
        model.anchor_rootdir = os.path.abspath(model_dir)
    seekr2_check.check_pre_simulation_all(model)
    starting_pdb_basename = os.path.basename(seekrflow.starting_pdb_filename)
    src_pdb_filename = os.path.join(seekrflow.work_directory, seekrflow.starting_pdb_filename)
    dest_pdb_relative_filename = os.path.join(structures.ROOT, starting_pdb_basename)
    dest_pdb_filename = os.path.join(seekrflow.work_directory, dest_pdb_relative_filename)
    if not os.path.exists(dest_pdb_filename):
        copyfile(src_pdb_filename, dest_pdb_filename)
    seekrflow.starting_pdb_filename = dest_pdb_relative_filename
    os.chdir(curdir)