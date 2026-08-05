"""
modules/seekr_input.py

Using the seekr API, create an input object and generate the model.
"""

import os
from shutil import copyfile

import seekrflow.modules.structures as structures
import seekrflow.modules.workflows.prepare as prepare

def prepare_model(
        seekrflow: structures.Seekrflow,
        #src_pdb_filename: str | None = None,
        force_overwrite: bool = False,
        skip_checks: bool = False,
        ) -> None:
    """
    Construct and prepare the model for seekr simulation.
    """
    #assert src_pdb_filename is not None, \
    #    "Source PDB filename must be provided to prepare the model."
    seekrflow.work_directory = os.path.abspath(
        os.path.expanduser(seekrflow.work_directory))
    curdir = os.getcwd()
    os.chdir(seekrflow.work_directory)
    # If MD system was left blank for parameterization, load canonical outputs.
    import seekrflow.parameterize as parameterize_module
    parameterize_module.fill_md_system_from_parameterize(seekrflow)
    root_directory = str(seekrflow.get_root_directory())

    model, model_json_path = prepare.prepare_workflow(
        seekrflow.workflow,
        root_directory,
        #src_pdb_filename,
        seekrflow.physical_attributes,
        force_overwrite=force_overwrite,
        skip_checks=skip_checks)
    
    #starting_pdb_basename = os.path.basename(src_pdb_filename)
    #dest_pdb_relative_filename = os.path.join(
    #    structures.ROOT, starting_pdb_basename)
    #dest_pdb_filename = os.path.join(
    #    seekrflow.work_directory, dest_pdb_relative_filename)
    #if not os.path.exists(dest_pdb_filename):
    #    copyfile(src_pdb_filename, dest_pdb_filename)
    os.chdir(curdir)
    return

    # TODO: work on force-overwrite functionality...
