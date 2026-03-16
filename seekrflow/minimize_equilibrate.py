"""
Minimize and equilibrate systems that have come from parameterize.py, or 
another source. This performs a 'gentle' procedure, where restraints are
put on the ligand and protein that are gradually released, and heating is
performed in stages.

This is a temporary script, used while seekrflow is focused on SEEKR2. Once
seekr 3 is available, this will occur within a stage.
"""

import os
import argparse

import seekrflow.modules.structures as structures

EQUIL_SEEKRFLOW_GLOB = "seekrflow_pre_equil_*.json"
EQUIL_SEEKRFLOW_BASE = "seekrflow_pre_equil_{}.json"

def minimize_equilibrate(
        seekrflow,
        reference_pdb_filename: str = "",
        ) -> str:
    startdir = os.path.abspath(os.getcwd())
    os.chdir(seekrflow.work_directory)
    assert seekrflow.workflow.type \
        in ["protein_ligand_seekr2", "protein_ligand_membrane_seekr2",
            "protein_protein_seekr2"], \
        "Only the following workflows are supported at this time:"\
        " - protein-ligand"\
        " - protein-ligand-membrane"\
        " - protein-protein"\
        " supported at this time for parameterize.py"
    equilibrated_pdb_filename = seekrflow.workflow.minimize_equilibrate(
        seekrflow.physical_attributes, unsolvated_pdb_filename=reference_pdb_filename)
    os.chdir(startdir)
    return equilibrated_pdb_filename

def main() -> None:
    argparser = argparse.ArgumentParser(
        description="Attempt to parameterize the system based on some input parameters.")
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
        help="Path to the input PDB file that contains the system to be "\
        "minimized/equilibrated.")
    argparser.add_argument(
        "-P", "--pdb_reference_system", dest="pdb_reference_system", 
        metavar="PDB_REFERENCE_SYSTEM", type=str, default="",
        help="Path to a PDB file that will be used as a reference structure."\
        "That is, the atomic positions within the reference PDB will be used "\
        "to assign restraints to 'known' heavy atoms to gradually release during "\
        "minimization and equilibration. If left empty, a nonsolvated structure "\
        "will be used.")
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
    
    args = argparser.parse_args()
    args = vars(args)
    input_json = args["input_json"]
    pdb_with_system = args["pdb_with_system"]
    ligand_indices = args["ligand_indices"]
    ligand_resname = args["ligand_resname"]
    work_dir = args["work_directory"]

    if input_json != "":
        if not os.path.exists(input_json):
            raise FileNotFoundError(f"Input JSON file {input_json} does not exist.")
        seekrflow = structures.load_seekrflow(input_json)
    else:
        seekrflow = structures.Seekrflow()

    if pdb_with_system == "":
        pdb_with_system = seekrflow.workflow.get_parameterizer_pdb_filename()
    else:
        seekrflow.workflow.set_parameterizer_pdb_filename(pdb_with_system)

    assert os.path.exists(seekrflow.workflow.get_parameterizer_pdb_filename()), \
        f"Input PDB file {seekrflow.workflow.get_parameterizer_pdb_filename()} does not exist."

    if work_dir == "":
        work_dir = seekrflow.work_directory
    seekrflow.make_work_directory(work_dir)
    seekrflow.handle_ligand_indices(ligand_indices, ligand_resname, pdb_with_system)
    equilibrated_pdb_filename = minimize_equilibrate(seekrflow)
    seekrflow.workflow.solvated_system_for_md.solvated_pdb = equilibrated_pdb_filename
    structures.save_new_seekrflow(seekrflow, EQUIL_SEEKRFLOW_GLOB, EQUIL_SEEKRFLOW_BASE,
                                  directory=seekrflow.work_directory)
    
if __name__ == "__main__":
    main()
