"""
modules/seekr_input.py

Using the seekr API, create an input object and generate the model.
"""

import os
from shutil import copyfile

import seekr2.modules.check as seekr2_check
import seekr2.prepare as seekr2_prepare

import seekrflow.modules.structures as structures

def prepare_model(
        seekrflow: structures.Seekrflow,
        force_overwrite: bool = False
        ) -> None:
    """
    Construct and prepare the model for seekr simulation.
    """
    seekrflow.work_directory = os.path.abspath(seekrflow.work_directory)
    curdir = os.getcwd()
    os.chdir(seekrflow.work_directory)
    root_directory = seekrflow.get_root_directory()
    param_directory = seekrflow.get_parameterize_directory()
    if seekrflow.workflow.type ==  "protein_ligand_seekr2":
        # Create the input object for the MMVT model
        if seekrflow.workflow.mmvt_settings.cv_type == "com_com_distance":
            model_input = seekrflow.workflow\
                .create_mmvt_receptor_ligand_com_com_model_input_seekr2(
                    root_directory=root_directory,
                    param_directory=param_directory,
                    physical_attributes=seekrflow.physical_attributes, 
                    alpha_carbon_ligand_threshold=0.6)
        else:
            raise NotImplementedError(
                f"Collective variable type "
                f"'{seekrflow.workflow.mmvt_settings.cv_type}' "
                f"is not implemented in protein_ligand_seekr2 workflow.")
    elif seekrflow.workflow.type == "protein_ligand_membrane_seekr2":
        if seekrflow.workflow.mmvt_settings.cv_type == "com_com_distance":
            model_input = seekrflow.workflow\
                .create_mmvt_receptor_ligand_com_com_model_input_seekr2(
                    root_directory=root_directory,
                    param_directory=param_directory,
                    physical_attributes=seekrflow.physical_attributes, 
                    alpha_carbon_ligand_threshold=0.6)
        else:
            raise NotImplementedError(
                f"Collective variable type "
                f"'{seekrflow.workflow.mmvt_settings.cv_type}' "
                f"is not implemented in protein_ligand_membrane_seekr2 workflow.")
    elif seekrflow.workflow.type == "protein_protein_seekr2":
        raise NotImplementedError(
            "Protein-protein SEEKR2 workflow is not implemented yet.")
    else:
        raise NotImplementedError(
            f"Workflow type {seekrflow.workflow.type} is not implemented.")

    model, xml_path = seekr2_prepare.prepare(model_input, force_overwrite=force_overwrite)
    if model.anchor_rootdir == ".":
        model_dir = os.path.dirname(xml_path)
        model.anchor_rootdir = os.path.abspath(model_dir)
    seekr2_check.check_pre_simulation_all(model)
    if seekrflow.workflow.solvated_system_for_md is None:
        starting_pdb_basename = "complex-equil.pdb"
        src_pdb_filename = os.path.join(structures.PARAMETERIZE, starting_pdb_basename)
    elif seekrflow.workflow.solvated_system_for_md.solvated_pdb == "":
        starting_pdb_basename = "complex-equil.pdb"
        src_pdb_filename = os.path.join(structures.PARAMETERIZE, starting_pdb_basename)
    else:
        starting_pdb_basename = os.path.basename(
            seekrflow.workflow.solvated_system_for_md.solvated_pdb)
        src_pdb_filename = seekrflow.workflow.solvated_system_for_md.solvated_pdb
    dest_pdb_relative_filename = os.path.join(structures.ROOT, starting_pdb_basename)
    dest_pdb_filename = os.path.join(seekrflow.work_directory, dest_pdb_relative_filename)
    if not os.path.exists(dest_pdb_filename):
        copyfile(src_pdb_filename, dest_pdb_filename)
    os.chdir(curdir)
