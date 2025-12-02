"""
parameterize.py

Given a starting PDB complex with a ligand and a receptor, this script
will parameterize the components of the system for input to a seekr
calculation.
"""

import os
import typing
import argparse
#import warnings
from shutil import copyfile

import mdtraj
from openmmforcefields.generators import SystemGenerator

import seekrflow.modules.base as base
import seekrflow.modules.structures as structures
import seekrflow.modules.workflows.structures as workflow_structures
import seekrflow.modules.parameters_topology_structures as parameters_topology_structures

PARAM_SEEKRFLOW_GLOB = "seekrflow_pre_param_*.json"
PARAM_SEEKRFLOW_BASE = "seekrflow_pre_param_{}.json"

# Lots of warnings are thrown from Espaloma and OpenFF, so we suppress them.
#warnings.filterwarnings("ignore")

def make_ligand_sdf_file(
        seekrflow: structures.Seekrflow,
        param_directory: str = structures.PARAMETERIZE
        ) -> None:
    """
    Make an SDF file for the ligand from the PDB file.
    """
    ligand_pdb_filename = os.path.join(
        param_directory, seekrflow.workflow.get_parametrizer_ligand_pdb_filename())
    ligand_sdf_filename = os.path.join(
        param_directory, seekrflow.workflow.get_parametrizer_default_sdf_filename())
    seekrflow.workflow.set_parametrizer_sdf_filename(ligand_sdf_filename)
    assert os.path.exists(ligand_pdb_filename), \
        f"Ligand PDB file {ligand_pdb_filename} does not exist."
    # Create input and output molecule streams
    from openeye import oechem
    ifs = oechem.oemolistream()
    ofs = oechem.oemolostream()
    # Open input PDB file and output SDF file
    if not ifs.open(str(ligand_pdb_filename)):
        oechem.OEThrow.Fatal(f"Unable to open the input PDB file: {ligand_pdb_filename}")
    if not ofs.open(ligand_sdf_filename):
        oechem.OEThrow.Fatal(f"Unable to create the output SDF file: {ligand_sdf_filename}")
    # Convert PDB to SDF using a generator for reading molecules
    for mol in ifs.GetOEGraphMols():
        oechem.OEWriteMolecule(ofs, mol)
    # Close the streams
    ifs.close()
    ofs.close()
    return

def parameterize(
    seekrflow: structures.Seekrflow,
    ) -> typing.Tuple[typing.Optional[str], typing.Optional[str]]:
    """
    Given a Seekrflow object including a PDB file receptor-ligand complex 
    as a starting structure, attempt to parameterize the file based on input 
    parameters.
    """
    assert seekrflow.workflow.type \
        in ["protein_ligand_seekr2", "protein_ligand_membrane_seekr2",
            "protein_protein_seekr2"], \
        "Only the following workflows are supported at this time:"\
        " - protein-ligand"\
        " - protein-ligand-membrane"\
        " - protein-protein"\
        " supported at this time for parameterize.py"
    work_param_dir = seekrflow.get_parameterize_directory()
    param_dir = structures.PARAMETERIZE
    curdir = os.getcwd()
    pdb_with_system = seekrflow.workflow.get_parametrizer_pdb_filename()
    work_copy_pdb_with_system = os.path.join(
        seekrflow.work_directory, os.path.basename(pdb_with_system))
    copyfile(pdb_with_system, work_copy_pdb_with_system)
    # Create a no-hydrogens structure version of work_copy_pdb_with_system
    work_copy = mdtraj.load(work_copy_pdb_with_system)
    non_h_indices = work_copy.topology.select("not element H")
    traj_noh = work_copy.atom_slice(non_h_indices)
    base, ext = os.path.splitext(work_copy_pdb_with_system)
    work_copy_pdb_with_system_noh = f"{base}_noh{ext}"
    print(f"Saving noh file at: {work_copy_pdb_with_system_noh}")
    traj_noh.save(work_copy_pdb_with_system_noh)
    
    seekrflow.workflow.set_parametrizer_pdb_filename(os.path.basename(work_copy_pdb_with_system))
    ligand_sdf_filename = seekrflow.workflow.get_parametrizer_sdf_filename()
    if ligand_sdf_filename != "":
        work_copy_ligand_sdf = os.path.join(
            work_param_dir, seekrflow.workflow.get_parametrizer_default_sdf_filename())
        copyfile(ligand_sdf_filename, work_copy_ligand_sdf)
    os.chdir(seekrflow.work_directory)
    seekrflow.workflow.split_molecules(param_dir)
    if seekrflow.workflow.has_small_molecule_ligand():
        if ligand_sdf_filename == "":
            make_ligand_sdf_file(seekrflow)

    if seekrflow.parameterizer is None:
        # TODO: logger?
        print("seekrflow.parameterizer is None, Skipping full parameterization.")
        if (seekrflow.workflow.solvated_system_for_md is not None) \
                and (seekrflow.workflow.solvated_system_for_md.solvated_pdb != "") \
                and (seekrflow.workflow.solvated_system_for_md.parameters_topology is not None):
            parmed_complex = seekrflow.workflow.solvated_system_for_md.parameters_topology\
                .make_parmed(seekrflow.workflow.solvated_system_for_md.solvated_pdb, 
                             seekrflow.work_directory)
            if seekrflow.workflow.bd_settings is not None:
                seekrflow.workflow.write_component_pqr_files(parmed_complex, structures.PARAMETERIZE)
        else:
            print("No solvated system for MD is set; cannot write PQR files. "
                  "There is nothing for parameterize.py to do.")
        os.chdir(curdir)
        return None, None
    else:
        assert seekrflow.workflow.solvated_system_for_md is None,\
            "Attempting to parametrize a system that has already be parametrized:"\
            "the solvated_system_for_md should be set to None for parameterize.py."
        serialized_xml, output_pdb_filename = seekrflow.workflow.create_complex(
            seekrflow.parameterizer,
            seekrflow.physical_attributes,
            seekrflow.workflow.md_settings,
            param_dir
        )
        starting_pdb_basename = os.path.basename(output_pdb_filename)
        seekrflow.workflow.solvated_system_for_md = workflow_structures.Solvated_system_for_md()
        seekrflow.workflow.solvated_system_for_md.solvated_pdb \
            = os.path.join(param_dir, starting_pdb_basename)
        system_basename = os.path.basename(serialized_xml)
        seekrflow.workflow.solvated_system_for_md.parameters_topology \
            = parameters_topology_structures.Openmm_system(
                system_filename=os.path.join(structures.PARAMETERIZE, system_basename),)
        parmed_complex = seekrflow.workflow.solvated_system_for_md.parameters_topology.make_parmed(
            output_pdb_filename)
        if seekrflow.workflow.bd_settings is not None:
            seekrflow.workflow.write_component_pqr_files(parmed_complex, structures.PARAMETERIZE)
        os.chdir(curdir)
        return serialized_xml, output_pdb_filename

def main() -> None:
    argparser = argparse.ArgumentParser(
        description="Given a PDB file of a system, attempt to parameterize " \
        "the system based on some input parameters.")
    argparser.add_argument(
        "-i", "--input_json", dest="input_json",
        metavar="INPUT_JSON", type=str, default="",
        help="Path to the input JSON file containing the parameters for "\
        "the parameterization. Default is empty, which will use the default "\
        "seekrflow parameters. If this is provided, then any other "\
        "command-line arguments will override the values in the JSON file.")
    argparser.add_argument(
        "-p", "--pdb_with_system", dest="pdb_with_system", 
        metavar="PDB_WITH_SYSTEM", type=str, default="",
        help="Path to the input PDB file that contains the parameterizable molecules.")
    argparser.add_argument(
        "-s", "--ligand_sdf_file", dest="ligand_sdf_file",
        metavar="LIGAND_SDF_FILE", type=str, default="",
        help="Optional path to the input SDF file containing the ligand molecule ")
    argparser.add_argument(
        "-L", "--ligand_resname", dest="ligand_resname", 
        metavar="LIGAND_RESNAME", type=str, default="",
        help="The residue name of the ligand molecule for automatic index "\
        "selection.")
    argparser.add_argument(
        "-l", "--ligand_indices", dest="ligand_indices", 
        metavar="LIGAND_INDICES", type=str, default="",
        help="A comma-separated list of integers defining site within the "\
        "ref_pdb structure.")
    argparser.add_argument(
        "-w", "--work_directory", dest="work_directory",
        metavar="WORK_DIRECTORY", type=str, default=structures.WORK,
        help="Path to the work directory for the parameterization.")
    argparser.add_argument(
        "-e", "--external_ff_file", dest="external_ff_file",
        metavar="EXTERNAL_FF_FILE", type=str, default=None,
        help="Any external force field file to use for the parameterization. "\
            "For example, 'espaloma-0.3.2.pt' can be entered here.")
    args = argparser.parse_args()
    args = vars(args)
    input_json = args["input_json"]
    pdb_with_system = args["pdb_with_system"]
    ligand_sdf_file = args["ligand_sdf_file"]
    ligand_indices = args["ligand_indices"]
    ligand_resname = args["ligand_resname"]
    work_dir = args["work_directory"]
    external_ff_file = args["external_ff_file"]

    # If the input JSON file is provided, then extract the necessary values
    # from it, but override with any command-line arguments provided.
    if input_json != "":
        if not os.path.exists(input_json):
            raise FileNotFoundError(f"Input JSON file {input_json} does not exist.")
        seekrflow = structures.load_seekrflow(input_json)
    else:
        seekrflow = structures.Seekrflow()

    if pdb_with_system == "":
        pdb_with_system = seekrflow.workflow.get_parametrizer_pdb_filename()
    else:
        seekrflow.workflow.set_parametrizer_pdb_filename(pdb_with_system)

    assert os.path.exists(seekrflow.workflow.get_parametrizer_pdb_filename()), \
        f"Input PDB file {seekrflow.workflow.get_parametrizer_pdb_filename()} does not exist."

    if ligand_sdf_file != "":
        assert os.path.exists(ligand_sdf_file), \
            "Nonexistent ligand SDF file: {}".format(ligand_sdf_file)
        seekrflow.workflow.parameterizer_information.ligand_sdf_file\
            = ligand_sdf_file
    if work_dir == "":
        work_dir = seekrflow.work_directory
    seekrflow.make_work_directory(work_dir)
    if external_ff_file is not None:
        if not os.path.exists(external_ff_file):
            raise FileNotFoundError(f"External force field file {external_ff_file} "
                                    "does not exist.")
        seekrflow.parameterizer.forcefield = external_ff_file
    seekrflow.handle_ligand_indices(ligand_indices, ligand_resname, pdb_with_system)
 
    system_filename, positions_filename = parameterize(seekrflow)
    seekrflow.work_directory = str(work_dir)
    # TODO: not going to save a new seekrflow file - attempt to use default values.
    # structures.save_new_seekrflow(seekrflow, PARAM_SEEKRFLOW_GLOB, PARAM_SEEKRFLOW_BASE,
    #                               directory=seekrflow.work_directory)
    return

if __name__ == "__main__":
    main()
