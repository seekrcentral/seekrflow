"""
modules/workflows/utilities.py

Useful miscellaneous utilities for workflow construction and execution.
"""

import typing

import numpy as np
import mdtraj

from Bio import Align
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment

d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

OPEN_GAP_SCORE = -10.0
EXTEND_GAP_SCORE = -0.5
END_OPEN_GAP_SCORE = 0.0
END_EXTEND_GAP_SCORE = 0.0

def align_sequence_from_mdtraj(
        traj1: typing.Any, 
        traj2: typing.Any, 
        atom_selection: str = "protein and not type H"):
    """
    Align two MDTraj trajectories based on their sequences, and provide
    the atom indices of the selected atom within the aligned residues.
    """
    atom_counter = 0
    residue_counter = 0
    traj1_resids = []
    traj1_resnames = []
    traj1_indices = traj1.topology.select(atom_selection)
    traj1_atom_info_list = []
    traj1_resids_in_selection = []
    for residue in traj1.topology.residues:
        #if not residue.is_protein:
        #    continue
        traj1_resids.append(residue_counter)
        traj1_resnames.append(residue.name)
        for atom in residue.atoms:
            if atom_counter in traj1_indices:
                traj1_atom_info = [atom_counter, atom.name, residue_counter]
                traj1_atom_info_list.append(traj1_atom_info)
                traj1_resids_in_selection.append(residue_counter)
            
            atom_counter += 1
        
        residue_counter += 1

    traj1_seq_list = []
    for resname in traj1_resnames:
        if resname in d3to1:
            traj1_seq_list.append(d3to1[resname])
        else:
            traj1_seq_list.append("?")
    traj1_sequence = "".join(traj1_seq_list)
    
    atom_counter = 0
    residue_counter = 0
    traj2_resids = []
    traj2_resnames = []
    traj2_indices = traj2.topology.select(atom_selection)
    traj2_atom_info_list = []
    traj2_resids_in_selection = []
    for residue in traj2.topology.residues:
        #if not residue.is_protein:
        #    continue
        traj2_resids.append(residue_counter)
        traj2_resnames.append(residue.name)
        for atom in residue.atoms:
            if atom_counter in traj2_indices:
                traj2_atom_info = [atom_counter, atom.name, residue_counter]
                traj2_atom_info_list.append(traj2_atom_info)
                traj2_resids_in_selection.append(residue_counter)
            
            atom_counter += 1
        
        residue_counter += 1

    traj2_seq_list = []
    for resname in traj2_resnames:
        if resname in d3to1:
            traj2_seq_list.append(d3to1[resname])
        else:
            #print("not found:", resname)
            traj2_seq_list.append("?")
    traj2_sequence = "".join(traj2_seq_list)

    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = OPEN_GAP_SCORE
    aligner.extend_gap_score = EXTEND_GAP_SCORE
    aligner.open_end_gap_score = END_OPEN_GAP_SCORE
    aligner.extend_end_gap_score = END_EXTEND_GAP_SCORE
    alignments = aligner.align(traj1_sequence, traj2_sequence)
    #print("alignments[0][0,:]:", alignments[0][0,:])
    #print("alignments[0][1,:]:", alignments[0][1,:])
    assert len(alignments[0][0,:]) == len(alignments[0][1,:])
    traj1_resid = 0
    traj2_resid = 0
    traj1_selected_indices = []
    traj2_selected_indices = []
    for i, (charA, charB) in enumerate(zip(alignments[0][0,:], alignments[0][1,:])):
        if (charA != "-") and (charB != "-"):
            traj1_atom_name_to_index = {}
            for atom in traj1.topology.residue(traj1_resid).atoms:
                if atom.index in traj1_indices:
                    traj1_atom_name_to_index[atom.name] = atom.index

            
            traj2_atom_name_to_index = {}
            for atom in traj2.topology.residue(traj2_resid).atoms:
                if atom.index in traj2_indices:
                    traj2_atom_name_to_index[atom.name] = atom.index

            
            # Now find common atom names in both residues
            common_atom_names = set(traj1_atom_name_to_index.keys()).intersection(
                set(traj2_atom_name_to_index.keys()))
            for atom_name in common_atom_names:
                traj1_selected_indices.append(traj1_atom_name_to_index[atom_name])
                traj2_selected_indices.append(traj2_atom_name_to_index[atom_name])
        if charA != "-":
            traj1_resid += 1
            if traj1_resid >= len(traj1_seq_list):
                break
            
        if charB != "-":
            traj2_resid += 1
            if traj2_resid >= len(traj2_seq_list):
                break

        
    assert len(traj1_selected_indices) == len(traj2_selected_indices), \
        "Number of selected aligned atom indices should be the same."
    
    return traj1_selected_indices, traj2_selected_indices

def obtain_md_ref_structure_map(
        ref_mdtraj: typing.Any,
        md_mdtraj: typing.Any,
        align_selection_str: str = "protein and not type H"
        ) -> tuple[typing.List[int], typing.List[int]]:
    """
    Since the ref structure is assumed to have been used to construct
    the md structure using, say, LEAP, then loop through the ref structure
    atom-by-atom and see which atom corresponds is in the md structure.
    These will be used to restrain "known" atoms - the ones that were in the 
    experimental structure.
    """
    SAME_ATOM_MAX_CUTOFF_DIST = 0.05  # half of an angstrom
    # alpha carbons of a particular residue - try changing residues if
    # error occurs
    #align_by_selection = "name CA and resname MET" 
    
    # These indices need to be chosen by sequence alignment of protein
    align_indices_solvated, align_indices_unsolvated = align_sequence_from_mdtraj(
        md_mdtraj, ref_mdtraj, align_selection_str)
    assert len(align_indices_solvated) == len(align_indices_unsolvated), \
        "Number of alignment atoms should be the same in solvated and unsolvated "\
        "structures."
    assert len(align_indices_solvated) > 0, \
        "No alignment atoms found between solvated and unsolvated structures."
    
    # Align the solvated structure using the unsolvated as a reference
    ref_mdtraj.superpose(
        md_mdtraj, atom_indices=align_indices_unsolvated, 
        ref_atom_indices=align_indices_solvated)

    # Compute distance matrix
    distances = cdist(ref_mdtraj.xyz[0,:,:], md_mdtraj.xyz[0,:,:])  # shape (N, M)
    
    # Check how many atoms are within cutoff
    within_cutoff = distances <= SAME_ATOM_MAX_CUTOFF_DIST
    num_within_cutoff = np.sum(within_cutoff)
    
    if num_within_cutoff == 0:
        # No atoms within cutoff - try relaxing the cutoff
        print(f"WARNING: No atoms within {SAME_ATOM_MAX_CUTOFF_DIST} cutoff.")
        print(f"Minimum distance found: {np.min(distances):.6f}")
        print("Consider increasing SAME_ATOM_MAX_CUTOFF_DIST or checking alignment.")
        
        # Try with a larger cutoff
        RELAXED_CUTOFF = 0.1  # 1 angstrom
        within_cutoff_relaxed = distances <= RELAXED_CUTOFF
        num_within_relaxed = np.sum(within_cutoff_relaxed)
        print(f"With relaxed cutoff of {RELAXED_CUTOFF}: {num_within_relaxed} pairs found")
        
        if num_within_relaxed == 0:
            raise ValueError("No matching atoms found even with relaxed cutoff. Check structures.")
    
    # Instead of using Hungarian algorithm on the full matrix, 
    # find all valid pairs within cutoff and then solve a feasible assignment
    within_cutoff_pairs = np.where(distances <= SAME_ATOM_MAX_CUTOFF_DIST)
    unsolved_rows = within_cutoff_pairs[0]
    unsolved_cols = within_cutoff_pairs[1]
    
    if len(unsolved_rows) == 0:
        raise ValueError("No atoms within cutoff distance found.")
    
    # Create a mapping of unique rows and columns that have matches
    unique_rows = np.unique(unsolved_rows)
    unique_cols = np.unique(unsolved_cols)
    
    # Create a reduced distance matrix with only rows/cols that have matches
    row_mapper = {old_idx: new_idx for new_idx, old_idx in enumerate(unique_rows)}
    col_mapper = {old_idx: new_idx for new_idx, old_idx in enumerate(unique_cols)}
    
    # Build reduced matrix
    n_rows = len(unique_rows)
    n_cols = len(unique_cols)
    distances_reduced = np.full((n_rows, n_cols), np.inf)
    
    for i, j in zip(unsolved_rows, unsolved_cols):
        new_i = row_mapper[i]
        new_j = col_mapper[j]
        distances_reduced[new_i, new_j] = distances[i, j]
    
    # Now solve the assignment on the reduced matrix
    row_ind_reduced, col_ind_reduced = linear_sum_assignment(distances_reduced)
    
    # Map back to original indices
    row_ind = unique_rows[row_ind_reduced]
    col_ind = unique_cols[col_ind_reduced]
    
    # Filter to only valid pairs (within original cutoff)
    valid_pairs = distances[row_ind, col_ind] <= SAME_ATOM_MAX_CUTOFF_DIST
    indices_unsolv = row_ind[valid_pairs]
    indices_solv = col_ind[valid_pairs]
    
    assert len(indices_unsolv) == len(indices_solv), \
        "Length of unsolvated and solvated index lists should be the same."
    assert len(indices_unsolv) > 0, \
        "No matching atoms found between solvated and unsolvated structures."
    # Could check for same atom names, but I cannot assume that during
    #  parametrization and solvation that atom names will be preserved.
    if len(indices_unsolv) < 100:
        print("WARNING: Less than 100 matching atoms found between solvated "
              "and unsolvated structures.")
    
    # Create the mapping list
    map_unsolv_to_solv_list = list(zip(indices_unsolv, indices_solv))
    unzipped_mapped_lists = list(zip(*map_unsolv_to_solv_list))
    ref_atom_indices = unzipped_mapped_lists[0]
    md_atom_indices = unzipped_mapped_lists[1]
    # Convert from int64 to int
    ref_atom_indices = [int(idx) for idx in ref_atom_indices]
    md_atom_indices = [int(idx) for idx in md_atom_indices]
    assert len(ref_atom_indices) == len(md_atom_indices), \
        "Number of mapped atoms between solvated and unsolvated structures should be the same"
    return ref_atom_indices, md_atom_indices

# NOTE: no longer used?
def select_ref_indices_from_md_indices(
        md_atom_indices: typing.List[int],
        ref_atom_indices: typing.List[int],
        md_selection_indices: typing.List[int]
        ) -> tuple[typing.List[int], typing.List[int]]:
    """
    Given a list of md atom indices and corresponding ref atom indices, and a selection of md atom indices,
    return the corresponding lists of md and ref atom indices that correspond to the selection.
    """
    selected_md_indices = []
    selected_ref_indices = []
    for md_idx, ref_idx in zip(md_atom_indices, ref_atom_indices):
        if md_idx in md_selection_indices:
            selected_md_indices.append(md_idx)
            selected_ref_indices.append(ref_idx)
    
    assert len(selected_md_indices) == len(selected_ref_indices), \
        "Number of selected MD and reference indices should be the same"
    return selected_md_indices, selected_ref_indices


def select_both_from_mdtraj_and_indices(
        md_traj: mdtraj.Trajectory, 
        md_indices: typing.List[int], 
        ref_traj: mdtraj.Trajectory, 
        ref_indices: typing.List[int], 
        selection: str) -> tuple[typing.List[int], typing.List[int]]:
    """
    Given an MDTraj trajectory and a list of atom indices,
    apply an MDTraj selection string to filter those indices.
    
    Parameters
    ----------
    traj : mdtraj.Trajectory
        The MDTraj trajectory object
    indices : list of int
        List of atom indices to filter
    selection : str
        MDTraj selection string
        
    Returns
    -------
    list of int
        Filtered list of atom indices matching the selection
    """
    filtered_md_indices = []
    filtered_ref_indices = []
    md_selected_indices = md_traj.topology.select(selection)
    ref_selected_indices = ref_traj.topology.select(selection)
    #print(f"MDTraj md_selection '{selection}' returned {len(md_selected_indices)} atoms.")
    #print(f"MDTraj ref_selection '{selection}' returned {len(ref_selected_indices)} atoms.")
    for md_idx, ref_idx in zip(md_indices, ref_indices):
        if md_idx in md_selected_indices and ref_idx in ref_selected_indices:
            filtered_md_indices.append(md_idx)
            filtered_ref_indices.append(ref_idx)
        
    return filtered_md_indices, filtered_ref_indices

#####################################################################################################
    #print("nonprotein_solvated_indices:", nonprotein_solvated_indices)
    #print("nonprotein_unsolvated_indices:", nonprotein_unsolvated_indices)
    # Loop through the solvated and unsolvated atoms in pair using zip, and
    # assert that elements are the same, at least
#    for solv_idx, unsolv_idx in zip(nonprotein_solvated_indices, nonprotein_unsolvated_indices):
#        solv_atom = md_mdtraj.topology.atom(solv_idx)
#        unsolv_atom = ref_mdtraj.topology.atom(unsolv_idx)
#        if solv_atom.element != unsolv_atom.element:
#            print("ERROR: Non-matching elements found between solvated and unsolvated "
#                  "structures for non-protein atom:")
#            print(f"  Solvated atom index: {solv_idx}, name: {solv_atom.name}, "
#                  f"element: {solv_atom.element}")
#            print(f"  Unsolvated atom index: {unsolv_idx}, name: {unsolv_atom.name}, "
#                  f"element: {unsolv_atom.element}")