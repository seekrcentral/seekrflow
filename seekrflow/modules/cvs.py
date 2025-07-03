"""
modules/cvs.py

Handle the creation of CVs for SEEKR calculations.
"""

import typing

import numpy as np
import mdtraj
import parmed

import seekrflow.modules.base as base
import seekrflow.modules.structures as structures

def alpha_carbon_selection_within_cutoff(
        traj: mdtraj.Trajectory, 
        ligand_indices: list, 
        distance_cutoff: float, 
        ligand_resname: str = "",
        receptor_selection: str = "protein"
        ) -> list:
    """
    Select alpha carbon indices of protein residues that are within a certain
    distance cutoff from any ligand atom.
    """
    contact_prot_indices = []
    contact_resids = set()
    if ligand_resname == "":
        protein_indices = traj.topology.select(receptor_selection)
    else:
        protein_indices = traj.topology.select(f"{receptor_selection} and not resname {ligand_resname}")

    for protein_index in protein_indices:
        if protein_index in ligand_indices:
            continue
        
        for ligand_index in ligand_indices:
            prot_atom_pos = traj.xyz[0, protein_index]
            lig_atom_pos = traj.xyz[0, ligand_index]
            dist = np.linalg.norm(prot_atom_pos - lig_atom_pos)
            if dist <= distance_cutoff:
                contact_prot_indices.append(protein_index)
    
    for contact_prot_index in contact_prot_indices:
        contact_atom = traj.topology.atom(contact_prot_index)
        resid = contact_atom.residue.index
        resname = contact_atom.residue.name
        contact_resids.add(resid)
    
    alpha_carbon_indices = []
    for resid in contact_resids:
        ca_selection_str = "name CA and resid {}".format(resid)
        ca_selection_list = traj.topology.select(ca_selection_str)
        if len(ca_selection_list) == 0:
            # then just continue this loop, no CA found for this residue
            continue
        ca_index = ca_selection_list[0]
        alpha_carbon_indices.append(ca_index)

    if len(alpha_carbon_indices) == 0 and receptor_selection != "protein":
        # No alpha carbons found, maybe it's not a protein? Try a different selection
        for resid in contact_resids:
            ca_selection_str = "not element H and resid {}".format(resid)
            ca_selection_list = list(traj.topology.select(ca_selection_str))
            if len(ca_selection_list) > 0:
                alpha_carbon_indices += ca_selection_list
    
    alpha_carbon_indices.sort()
    return alpha_carbon_indices

def get_receptor_ligand_com_com_selections(
        seekrflow: structures.Seekrflow,
        complex_pdb_filename: str,
        alpha_carbon_ligand_threshold: float = 0.6,
        ) -> typing.Tuple[list, list]:
    """
    Get the selections for the receptor and ligand in a com-com calculation.
    
    Parameters
    ----------
    complex_pdb_filename : str
        The filename of the complex PDB file.
    
    Returns
    -------
    tuple
        A tuple containing the receptor selection and ligand selection.
    """
    # This is a placeholder implementation. The actual implementation would
    # depend on the specific structure of the PDB file and how the receptor
    # and ligand are defined.
    
    if len(seekrflow.ligand_indices) == 0:
        assert seekrflow.ligand_resname != "", \
            "ligand_resname must be set in the input JSON file if ligand_indices is empty."
        seekrflow.ligand_indices = base.get_ligand_indices(complex_pdb_filename, seekrflow.ligand_resname)
    traj = mdtraj.load(complex_pdb_filename)
    if seekrflow.receptor_selection == "":
        receptor_selection = "protein"
    else:
        receptor_selection = seekrflow.receptor_selection
    receptor_selection = alpha_carbon_selection_within_cutoff(
        traj, 
        seekrflow.ligand_indices, 
        alpha_carbon_ligand_threshold, 
        ligand_resname=seekrflow.ligand_resname,
        receptor_selection=receptor_selection
    )
    return receptor_selection, seekrflow.ligand_indices

def get_bd_receptor_ligand_selections(
        seekrflow: structures.Seekrflow,
        complex_pdb_filename: str,
        md_receptor_selection: list,
        md_ligand_selection: list,
        ) -> typing.Tuple[list, list]:
    traj = mdtraj.load(complex_pdb_filename)
    receptor_pqr_filename = seekrflow.bd_settings.receptor_pqr_filename
    ligand_pqr_filename = seekrflow.bd_settings.ligand_pqr_filename
    receptor_pqr_parmed = parmed.load_file(receptor_pqr_filename)
    ligand_pqr_parmed = parmed.load_file(ligand_pqr_filename)
    receptor_selection_atom_name_list = []
    receptor_selection_resid_list = []
    for md_receptor_index in md_receptor_selection:
        rec_atom_name = traj.topology.atom(md_receptor_index).name
        rec_atom_resid = traj.topology.atom(md_receptor_index).residue.index
        receptor_selection_atom_name_list.append(rec_atom_name)
        receptor_selection_resid_list.append(str(rec_atom_resid))

    bd_receptor_indices = []

    for name, resid in zip(
            receptor_selection_atom_name_list, receptor_selection_resid_list):
        for i, atom in enumerate(receptor_pqr_parmed.atoms):
            if atom.name == name and str(atom.residue.number-1) == resid:
                bd_receptor_indices.append(i)
        
    ligand_selection_atom_name_list = []
    for md_ligand_index in md_ligand_selection:
        lig_atom_name = traj.topology.atom(md_ligand_index).name
        ligand_selection_atom_name_list.append(lig_atom_name)

    bd_ligand_indices = []
    for name in ligand_selection_atom_name_list:
        for i, atom in enumerate(ligand_pqr_parmed.atoms):
            if atom.name == name:
                bd_ligand_indices.append(i)

    assert len(bd_receptor_indices) == len(md_receptor_selection), \
        "BD receptor indices do not match MD receptor selection. " \
        "This may indicate duplicate structure in the complex or in the "\
        "PQR files. Be sure to check all structures."
    assert len(bd_ligand_indices) == len(md_ligand_selection), \
        "BD ligand indices do not match MD ligand selection. " \
        "This may indicate duplicate structure in the complex or in the "\
        "PQR files. Be sure to check all structures."
    return bd_receptor_indices, bd_ligand_indices