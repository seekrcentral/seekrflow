"""
modules/base.py

Base routines for seekrflow.
"""

import typing

import mdtraj

def initialize_ref_indices(
        ref_indices: str
    ) -> typing.List[int]:
    """
    Convert a string with comma-separated integers into a true list of integers.
    """
    assert ref_indices != ""
    ref_integers = ref_indices.split(",")
    for ref_integer in ref_integers:
        # Variable not used - just catches integer conversion errors
        ref_integer_int = int(ref_integer)
        assert ref_integer_int >= 0, \
            "Reference indices must be non-negative integers. Found: {}".format(ref_integer)
    
    return list(map(int, ref_integers))

def get_ligand_indices(
        protein_ligand_pdb: str, 
        ligand_resname: str, 
        include_H: bool = False
        ) -> typing.List[int]:
    """
    Given a PDB file and a ligand residue name, return the indices of the ligand
    atoms in the PDB file.
    """
    traj = mdtraj.load(protein_ligand_pdb)
    if include_H:
        selection_string = "resname {}"
    else:
        selection_string = "resname {} and not type H"
    ligand_indices = list(traj.topology.select(selection_string.format(ligand_resname)))
    assert len(ligand_indices), "Selection not found: resname {}".format(ligand_resname)
    ligand_indices_ints = [int(index) for index in ligand_indices]
    return ligand_indices_ints