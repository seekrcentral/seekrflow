"""
OpenMM Membrane/Globular System Equilibration

This module provides gentle minimization and equilibration for globular protein
or protein-membrane systems following best practices. It handles both soluble 
proteins and membrane-embedded systems with specific restraints.

Key Features:
- stage-by-stage progressive restraint relaxation (AMBER/CHARMM-GUI standard)
- Backbone and sidechain specific restraints on 'known' portions of protein
- Z-only phosphorus restraints (maintains bilayer thickness)
- Lipid headgroup and dihedral restraints
- Semi-isotropic NPT for membranes
- Automatic lipid atom detection from topology
- Bilayer property monitoring (thickness, area per lipid, tail extension)

Globular Equilibration Protocol:
    Stage 1 (E1): Heat 100K → target T with strong 'known' protein restraints
    Stage 2 (E2): NPT at target T with full restraints until energy stabilizes
    Stage 3 (E3): Heat 100K → target T with reduced restraints
    Stage 4 (E4): NPT with backbone - only restraints
    Stage 5 (E5): Unrestrained simulation

Membrane Equilibration Protocol:
    Stage 1 (E1): Heat 100K → target T with strong protein + lipid restraints
    Stage 2 (E2): NPT at target T with moderate restraints
    Stage 3 (E3): NPT with reduced restraints
    Stage 4 (E4): NPT with light restraints
    Stage 5 (E5): NPT with backbone-only restraints
    Stage 6 (E6): NPT with no positional restraints (dihedral only)
"""

import os
from sys import stdout

from Bio import Align
import openmm
import openmm.app as openmm_app
import openmm.unit as unit
import mdtraj
import numpy as np
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

FORCE_CONST_NAME_BACKBONE = "backbone_positional_force_const"
FORCE_CONST_NAME_SIDECHAIN = "sidechain_positional_force_const"
FORCE_CONST_NAME_NONPROTEIN = "nonprotein_positional_force_const"
FORCE_CONST_NAME_PHOS = "phosphorus_z_force_const"
FORCE_CONST_NAME_DIHE_PROPER = "lipid_dihedral_proper_force_const"
FORCE_CONST_NAME_DIHE_IMPROPER = "lipid_dihedral_improper_force_const"

DEFAULT_RANDOM_SEED = 111

# Equilibration constants
# Stage 1 (E1): Strong restraints during initial heating
DEFAULT_STAGE_1_K_BACKBONE = 100.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_1_K_SIDECHAIN = 100.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_1_K_NONPROTEIN = 100.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_1_K_HEAD = 5.0*unit.kilocalories_per_mole/unit.angstroms**2  # Z-restraint
DEFAULT_STAGE_1_K_DIHE = 250.0*unit.kilocalories_per_mole
DEFAULT_STAGE_1_TIMESTEP = 0.001 * unit.picoseconds
DEFAULT_STAGE_1_START_TEMP = 100.0*unit.kelvin
DEFAULT_STAGE_1_STEPS_PER_RAMP = 10000  # 10 ps with 1 fs time step

# Stage 2 (E2): Strong restraints at target temperature
DEFAULT_STAGE_2_K_BACKBONE = 100.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_2_K_SIDECHAIN = 100.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_2_K_NONPROTEIN = 100.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_2_K_HEAD = 5.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_2_K_DIHE = 100.0*unit.kilocalories_per_mole
DEFAULT_STAGE_2_TIMESTEP = 0.001 * unit.picoseconds
DEFAULT_STAGE_2_TOTAL_STEPS = 100000  # 100 ps with 1 fs time step

# Stage 3 (E3): Moderate restraints during second heating
DEFAULT_STAGE_3_K_BACKBONE = 10.0*unit.kilocalories_per_mole/unit.angstroms**2 
DEFAULT_STAGE_3_K_SIDECHAIN = 10.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_3_K_NONPROTEIN = 10.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_3_K_HEAD = 2.5*unit.kilocalories_per_mole/unit.angstroms**2  # Z-restraint
DEFAULT_STAGE_3_K_DIHE = 100.0*unit.kilocalories_per_mole 
DEFAULT_STAGE_3_TIMESTEP = 0.001 * unit.picoseconds
DEFAULT_STAGE_3_START_TEMP = 100.0*unit.kelvin
DEFAULT_STAGE_3_STEPS_PER_RAMP = 10000  # 10 ps with 1 fs time step

# Stage 4 (E4): Moderate restraints, longer equilibration
DEFAULT_STAGE_4_K_BACKBONE = 10.0*unit.kilocalories_per_mole/unit.angstroms**2 
DEFAULT_STAGE_4_K_SIDECHAIN = 10.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_4_K_NONPROTEIN = 10.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_4_K_HEAD = 2.5*unit.kilocalories_per_mole/unit.angstroms**2  # Z-restraint
DEFAULT_STAGE_4_K_DIHE = 100.0*unit.kilocalories_per_mole 
DEFAULT_STAGE_4_TIMESTEP = 0.001 * unit.picoseconds
DEFAULT_STAGE_4_TOTAL_STEPS = 100000  # 100 ps with 1 fs time step

# Stage 5 (E5): Moderate restraints on protein, lower restraints on membrane
DEFAULT_STAGE_5_K_BACKBONE = 10.0*unit.kilocalories_per_mole/unit.angstroms**2 
DEFAULT_STAGE_5_K_SIDECHAIN = 10.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_5_K_NONPROTEIN = 10.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_5_K_HEAD = 1.0*unit.kilocalories_per_mole/unit.angstroms**2  # Z-restraint
DEFAULT_STAGE_5_K_DIHE = 50.0*unit.kilocalories_per_mole 
DEFAULT_STAGE_5_TIMESTEP = 0.001 * unit.picoseconds
DEFAULT_STAGE_5_TOTAL_STEPS = 100000  # 100 ps with 1 fs time step

# Stage 6 (E6): Lower restraints on protein, further lower restraints on membrane
DEFAULT_STAGE_6_K_BACKBONE = 1.0*unit.kilocalories_per_mole/unit.angstroms**2 
DEFAULT_STAGE_6_K_SIDECHAIN = 1.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_6_K_NONPROTEIN = 1.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_6_K_HEAD = 0.5*unit.kilocalories_per_mole/unit.angstroms**2  # Z-restraint
DEFAULT_STAGE_6_K_DIHE = 50.0*unit.kilocalories_per_mole 
DEFAULT_STAGE_6_TIMESTEP = 0.002 * unit.picoseconds
DEFAULT_STAGE_6_TOTAL_STEPS = 200000  # 400 ps with 2 fs time step

# Stage 7 (E7): Lower restraints on protein backbone only, further lower restraints on membrane
DEFAULT_STAGE_7_K_BACKBONE = 1.0*unit.kilocalories_per_mole/unit.angstroms**2 
DEFAULT_STAGE_7_K_SIDECHAIN = 0.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_7_K_NONPROTEIN = 0.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_7_K_HEAD = 0.1*unit.kilocalories_per_mole/unit.angstroms**2  # Z-restraint
DEFAULT_STAGE_7_K_DIHE = 25.0*unit.kilocalories_per_mole 
DEFAULT_STAGE_7_TIMESTEP = 0.002 * unit.picoseconds
DEFAULT_STAGE_7_TOTAL_STEPS = 200000  # 400 ps with 2 fs time step

# Stage 8 (E8): No restraints (production-like)
DEFAULT_STAGE_8_K_BACKBONE = 0.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_8_K_SIDECHAIN = 0.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_8_K_NONPROTEIN = 0.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_8_K_PHOS = 0.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_8_K_HEAD = 0.0*unit.kilocalories_per_mole/unit.angstroms**2
DEFAULT_STAGE_8_K_DIHE = 0.0*unit.kilocalories_per_mole
DEFAULT_STAGE_8_TIMESTEP = 0.002 * unit.picoseconds
DEFAULT_STAGE_8_TOTAL_STEPS = 500000  # 1 ns with 2 fs time step

# Trajectory and data output filenames
#EQ_TRAJ_FILENAME = "eq_traj.pdb"
DEFAULT_EQ_TRAJ_FILENAME = "eq_traj.dcd"
DEFAULT_EQ_TRAJ_INTERVAL = 10000  # in steps
DEFAULT_EQ_STATE_FILENAME = "eq_state.dat"
DEFAULT_EQ_STATE_INTERVAL = 1000  # in steps

# Final structure output
OUTPUT_PDB_FILE = "equilibrated.pdb"
OUTPUT_PDB_FILE_IMAGED = "equilibrated_imaged.pdb"

"""
import os
import glob

import openmm
import openmm.app as openmm_app
import openmm.unit as unit
import md_min_equil_gentle

# For membrane protein systems:
# We will need the same topology
psf = openmm_app.CharmmPsfFile("charmm-gui-6540578083/step5_assembly.psf")

# NOTE: box vectors will be changed later by the state file.
psf.setBox(115.039972*unit.angstroms, 115.039972*unit.angstroms, 154.53262*unit.angstroms)

unsolvated_pdb_filename = "9ARZ.pdb"
solvated_pdb_filename = "charmm-gui-6540578083/step5_assembly.pdb"
pdb = openmm_app.PDBFile(solvated_pdb_filename)

all_charmm_files = glob.glob("charmm-gui-6540578083/toppar/*")
#print("all_charmm_files:")
#print([os.path.abspath(file) for file in all_charmm_files])
#exit()

params = openmm_app.CharmmParameterSet(*all_charmm_files)
nonbonded_cutoff = 1.2*unit.nanometer
hydrogenMass = 3.0*unit.daltons
topology = psf.topology
system = psf.createSystem(
    params,
    nonbondedMethod=openmm_app.PME,
    nonbondedCutoff=nonbonded_cutoff,
    constraints=openmm_app.HBonds,
    hydrogenMass=hydrogenMass
)

temperature = 300.00*unit.kelvin
target_pressure = 1.0*unit.bar

md_min_equil_gentle.main(
    unsolvated_pdb_filename,
    solvated_pdb_filename,
    system,
    topology,
    temperature,
    target_pressure,
    is_membrane_system=True
)

"""

# Custom reporter class for detailed energy components
class DetailedEnergyReporter(object):
    """Custom reporter to output detailed energy components"""
    
    def __init__(self, file, reportInterval, phosphorus_indices=None):
        self._file = file
        self._reportInterval = reportInterval
        self._phosphorus_indices = phosphorus_indices
        self._hasInitialized = False
        
    def __del__(self):
        if hasattr(self, '_file') and self._file is not None:
            self._file.close()
    
    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, True, False, None)
    
    def report(self, simulation, state):            
        # Get energy by force group
        energies = {}
        for i, name in enumerate(['Bond', 'Angle', 'Dihedral', 'Nonbonded']):
            try:
                energy = simulation.context.getState(getEnergy=True, groups={i+1}).getPotentialEnergy()
                energies[name] = energy.value_in_unit(unit.kilocalories_per_mole)
            except:
                energies[name] = 0.0
        
        # Get total energy and other properties
        total_state = simulation.context.getState(getEnergy=True, getPositions=True)
        total_energy = total_state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
        kinetic_energy = total_state.getKineticEnergy().value_in_unit(unit.kilocalories_per_mole)
        temp = simulation.context.getIntegrator().getTemperature().value_in_unit(unit.kelvin)
        
        # Get volume and XY area if periodic
        try:
            box_vectors = total_state.getPeriodicBoxVectors()
            if box_vectors is not None:
                volume = (box_vectors[0][0] * box_vectors[1][1] * box_vectors[2][2]).value_in_unit(unit.nanometer**3)
                
                # Calculate XY area using cross product to handle sheared boxes
                vec_x = np.array([box_vectors[0][0].value_in_unit(unit.nanometer),
                                  box_vectors[0][1].value_in_unit(unit.nanometer),
                                  box_vectors[0][2].value_in_unit(unit.nanometer)])
                vec_y = np.array([box_vectors[1][0].value_in_unit(unit.nanometer),
                                  box_vectors[1][1].value_in_unit(unit.nanometer),
                                  box_vectors[1][2].value_in_unit(unit.nanometer)])
                cross_product = np.cross(vec_x, vec_y)
                xy_area = np.linalg.norm(cross_product)
            else:
                volume = 0.0
                xy_area = 0.0
        except:
            volume = 0.0
            xy_area = 0.0
        
        # Calculate bilayer thickness if phosphorus atoms are provided
        bilayer_thickness = 0.0
        if self._phosphorus_indices is not None and len(self._phosphorus_indices) > 0:
            try:
                positions = total_state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                p_z_coords = positions[self._phosphorus_indices, 2]
                
                if len(p_z_coords) >= 2:
                    # Compute overall P COM z-coordinate
                    p_com_z = np.mean(p_z_coords)
                    
                    # Split into upper and lower leaflets
                    upper_leaflet_z = p_z_coords[p_z_coords > p_com_z]
                    lower_leaflet_z = p_z_coords[p_z_coords < p_com_z]
                    
                    # Calculate thickness if both leaflets have atoms
                    if len(upper_leaflet_z) > 0 and len(lower_leaflet_z) > 0:
                        upper_com_z = np.mean(upper_leaflet_z)
                        lower_com_z = np.mean(lower_leaflet_z)
                        bilayer_thickness = abs(upper_com_z - lower_com_z)
            except Exception as e:
                bilayer_thickness = 0.0
            
        # Write to file
        if simulation.currentStep == 0:
            # Write header
            self._file.write("#Step,Potential_Energy_kcal_mol,Kinetic_Energy_kcal_mol,Total_Energy_kcal_mol,Temperature_K,Box_Volume_nm3,Box_XY_Area_nm2,Bilayer_Thickness_nm,Bond_Energy_kcal_mol,Angle_Energy_kcal_mol,Dihedral_Energy_kcal_mol,Nonbonded_Energy_kcal_mol\n")
        
        self._file.write(f'{simulation.currentStep},{total_energy:.4f},{kinetic_energy:.4f},{total_energy+kinetic_energy:.4f},{temp:.2f},{volume:.4f},{xy_area:.4f},{bilayer_thickness:.4f},{energies["Bond"]:.4f},{energies["Angle"]:.4f},{energies["Dihedral"]:.4f},{energies["Nonbonded"]:.4f}\n')
        self._file.flush()

def align_sequence_from_mdtraj(traj1, traj2, atom_selection="protein and not type H"):
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
        if not residue.is_protein:
            continue
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
            #print("not found:", resname)
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
        if not residue.is_protein:
            continue
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
    aligner.end_open_gap_score = END_OPEN_GAP_SCORE
    aligner.end_extend_gap_score = END_EXTEND_GAP_SCORE
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

def find_known_atom_indices(solvated_struct, unsolvated_struct):
    """
    Since the unsolvated structure is assumed to have been used to construct
    the solvated one using, say, LEAP, then loop through the unsolvated structure
    atom-by-atom and see which atom in corresponds is in the solvated structure.
    These will be used to restrain "known" atoms - the ones that were in the 
    experimental structure.
    """
    SAME_ATOM_MAX_CUTOFF_DIST = 0.05  # half of an angstrom
    # alpha carbons of a particular residue - try changing residues if
    # error occurs
    #align_by_selection = "name CA and resname MET" 
    
    # These indices need to be chosen by sequence alignment of protein
    align_indices_solvated, align_indices_unsolvated = align_sequence_from_mdtraj(
        solvated_struct, unsolvated_struct)
    
    assert len(align_indices_solvated) == len(align_indices_unsolvated), \
        "Number of alignment atoms should be the same in solvated and unsolvated "\
        "structures."
    
    # Align the solvated structure using the unsolvated as a reference
    unsolvated_struct.superpose(
        solvated_struct, atom_indices=align_indices_unsolvated, 
        ref_atom_indices=align_indices_solvated)

    #solvated_struct.save("/home/lvotapka/tmp/align1.pdb")
    #unsolvated_struct.save("/home/lvotapka/tmp/align2.pdb")
    
    # But these only include protein atoms
    #map_unsolv_to_solv_list = list(zip(align_indices_unsolvated, align_indices_solvated))
    #return map_unsolv_to_solv_list

    # NOTE: No longer Hungarian algorithm

    # Compute distance matrix
    distances = cdist(unsolvated_struct.xyz[0,:,:], solvated_struct.xyz[0,:,:])  # shape (N, M)
    
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
    
    return map_unsolv_to_solv_list

def is_lipid_residue(residue):
    """
    Determine if a residue is a lipid based on multiple criteria.
    
    Uses three filtering approaches:
    1. Segment ID check (CHARMM-GUI typically uses "MEMB" for membrane)
    2. Exclusion of protein and water residues
    3. Common lipid residue name patterns
    
    Parameters
    ----------
    residue : mdtraj.Topology.Residue
        The residue to check
        
    Returns
    -------
    bool
        True if residue is a lipid, False otherwise
    """
    # Common lipid residue names in CHARMM force field
    LIPID_NAMES = {
        'POPC', 'POPE', 'POPS', 'POPG',  # PC, PE, PS, PG lipids
        'PSPC', 'SDPC', 'PAPC', 'SSM',
        'DOPC', 'DOPE', 'DOPS', 'DOPG',
        'DPPC', 'DPPE', 'DPPS', 'DPPG',
        'DLPC', 'DLPE', 'DLPS', 'DLPG',
        'DMPC', 'DMPE', 'DMPS', 'DMPG',
        'CHOL', 'CHL1',  # Cholesterol
        'PALM', 'STEA', 'OLEO', 'LINO',  # Fatty acids
        'SAPI', 'SAPI24', 'SAPI25',  # Saposin
        'POPA', 'DOPA', 'DPPA', 'DLPA', 'DMPA',  # PA lipids
    }
    
    # Check if residue is protein or water (exclude these)
    if residue.is_protein or residue.is_water:
        return False
    
    # Check segment ID (CHARMM-GUI uses "MEMB" for membrane)
    if hasattr(residue, 'segment_id') and residue.segment_id == "MEMB":
        return True
    
    # Check against known lipid names
    if residue.name in LIPID_NAMES:
        return True
    
    return False

def identify_lipid_atoms_with_mdtraj(solvated_mdtraj):
    """
    Identify lipid atom groups (phosphorus, headgroup heavy atoms, cis double bonds)
    using MDTraj topology analysis and CHARMM naming conventions.
    
    Only processes atoms from lipid residues (excludes protein, water, ligands).
    
    Uses atom naming patterns to identify:
    - Headgroup atoms: P and atoms named [CNOS]1[0-9] (e.g., O11, C12, N11)
    - Tail carbons: C2X and C3X patterns
    - Double bonds: sp2 carbons (2 carbon neighbors + 1 hydrogen)
    - Glycerol atoms for improper dihedral.
    
    Returns empty dict if no phosphorus atoms detected.
    
    Parameters
    ----------
    solvated_mdtraj : mdtraj.Trajectory
        MDTraj trajectory object of solvated system
        
    Returns
    -------
    dict
        Dictionary with keys:
        - 'phosphorus': list of P atom indices
        - 'headgroup_heavy': list of headgroup heavy atom indices
        - 'cis_double_bonds': list of 4-atom tuples for cis double bonds
        - 'glycerol_impropers': list of 4-atom tuples for glycerol impropers
    """
    import re
    
    result = {
        'phosphorus': [],
        'headgroup_heavy': [],
        'cis_double_bonds': [],
        'glycerol_impropers': []
    }
    
    print("Detecting lipid atoms using CHARMM naming conventions...")
    
    # Step 1: Find phosphorus atoms
    try:
        p_indices = solvated_mdtraj.topology.select("element P")
        if len(p_indices) == 0:
            print("WARNING: No phosphorus atoms found. Skipping lipid restraints.")
            return result
        
        print(f"Found {len(p_indices)} phosphorus atoms in system.")
        result['phosphorus'] = list(p_indices)
        
    except Exception as e:
        print(f"WARNING: Phosphorus detection failed: {e}")
        print("Proceeding without lipid restraints.")
        return result
    
    # Step 2: Identify headgroup atoms by naming pattern
    # Pattern: [CNOS]1[0-9] - element symbol, then "1", then a digit
    # Only process atoms from lipid residues
    headgroup_pattern = re.compile(r'^[CNOS]1\d$')
    headgroup_indices = set(p_indices)
    
    for atom in solvated_mdtraj.topology.atoms:
        if is_lipid_residue(atom.residue) and headgroup_pattern.match(atom.name):
            headgroup_indices.add(atom.index)
    
    result['headgroup_heavy'] = sorted(list(headgroup_indices))
    print(f"Identified {len(result['headgroup_heavy'])} headgroup heavy atoms (P + [CNOS]1X pattern).")
    
    # Step 3: Build neighbor map for topology analysis
    print("Building atom connectivity map...")
    neighbor_map = {i: {'C': [], 'H': [], 'other': []} for i in range(solvated_mdtraj.n_atoms)}
    
    for bond in solvated_mdtraj.topology.bonds:
        atom1, atom2 = bond[0], bond[1]
        idx1, idx2 = atom1.index, atom2.index
        elem1, elem2 = atom1.element.symbol, atom2.element.symbol

        # Categorize neighbors by element
        if elem2 == 'C':
            neighbor_map[idx1]['C'].append(idx2)
        elif elem2 == 'H':
            neighbor_map[idx1]['H'].append(idx2)
        else:
            neighbor_map[idx1]['other'].append(idx2)
        
        if elem1 == 'C':
            neighbor_map[idx2]['C'].append(idx1)
        elif elem1 == 'H':
            neighbor_map[idx2]['H'].append(idx1)
        else:
            neighbor_map[idx2]['other'].append(idx1)
    
    # Step 4: Identify tail carbons (C2X or C3X pattern)
    # Only process atoms from lipid residues
    tail_pattern = re.compile(r'^C[23]\d+$')
    tail_carbons = []
    
    for atom in solvated_mdtraj.topology.atoms:
        if is_lipid_residue(atom.residue) and tail_pattern.match(atom.name) and atom.element.symbol == 'C':
            tail_carbons.append(atom.index)
    
    print(f"Found {len(tail_carbons)} tail carbon atoms (C2X/C3X pattern).")
    
    # Step 5: Identify sp2 carbons (double bond carbons)
    # sp2 carbons have exactly 2 carbon neighbors and 1 hydrogen neighbor
    sp2_carbons = []
    
    for carbon_idx in tail_carbons:
        atom = solvated_mdtraj.topology.atom(carbon_idx)
        neighbors = neighbor_map[carbon_idx]
        num_C = len(neighbors['C'])
        num_H = len(neighbors['H'])
        
        # sp2 hybridization: 2 carbons + 1 hydrogen (3 total bonds, ignoring implicit H)
        if num_C == 2 and num_H == 1:
            sp2_carbons.append(carbon_idx)
    
    print(f"Identified {len(sp2_carbons)} sp2 carbons (potential double bond carbons).")
    
    # Step 6: Find adjacent pairs of sp2 carbons (these form C=C double bonds)
    double_bond_pairs = []
    
    for sp2_idx in sp2_carbons:
        carbon_neighbors = neighbor_map[sp2_idx]['C']
        for neighbor_idx in carbon_neighbors:
            if neighbor_idx in sp2_carbons:
                # Found a C=C bond - store as sorted tuple to avoid duplicates
                pair = tuple(sorted([sp2_idx, neighbor_idx]))
                if pair not in double_bond_pairs:
                    double_bond_pairs.append(pair)
    
    print(f"Found {len(double_bond_pairs)} C=C double bonds.")
    
    # Step 7: Build 4-atom dihedral tuples (C-C=C-C) for each double bond
    for c_i, c_j in double_bond_pairs:
        # c_i and c_j are the double bond carbons
        # Find carbon bonded to c_i (but not c_j) = atom h
        # Find carbon bonded to c_j (but not c_i) = atom k
        
        c_i_neighbors = [n for n in neighbor_map[c_i]['C'] if n != c_j]
        c_j_neighbors = [n for n in neighbor_map[c_j]['C'] if n != c_i]
        
        if len(c_i_neighbors) > 0 and len(c_j_neighbors) > 0:
            # Take the first carbon neighbor on each side
            atom_h = c_i_neighbors[0]
            atom_k = c_j_neighbors[0]
            
            # Dihedral is (h, i, j, k)
            dihedral = (atom_h, c_i, c_j, atom_k)
            result['cis_double_bonds'].append(dihedral)
    
    print(f"Created {len(result['cis_double_bonds'])} dihedral restraints for cis double bonds.")
    
    # Step 8: Identify glycerol backbone stereocenter improper dihedrals
    # For non-cholesterol lipids, maintain chirality at C2 by restraining
    # the improper dihedral between C1, C3, C2, and O21
    print("\nDetecting glycerol backbone stereocenters...")
    
    # Cholesterol residue names to exclude
    CHOLESTEROL_NAMES = {'CHOL', 'CHL1'}
    
    # Group atoms by residue for efficient lookup
    residue_atoms = {}
    for atom in solvated_mdtraj.topology.atoms:
        if is_lipid_residue(atom.residue) and atom.residue.name not in CHOLESTEROL_NAMES:
            res_key = (atom.residue.chain.index, atom.residue.index)
            if res_key not in residue_atoms:
                residue_atoms[res_key] = {'atoms': [], 'residue': atom.residue}
            residue_atoms[res_key]['atoms'].append(atom)
    
    print(f"Processing {len(residue_atoms)} non-cholesterol lipid residues...")
    
    # Process each lipid residue
    impropers_found = 0
    for res_key, res_data in residue_atoms.items():
        atoms = res_data['atoms']
        residue = res_data['residue']
        
        # Find the glycerol backbone atoms by name
        c1_atom = None
        c2_atom = None
        c3_atom = None
        o21_atom = None
        
        for atom in atoms:
            if atom.name == 'C1':
                c1_atom = atom
            elif atom.name == 'C2':
                c2_atom = atom
            elif atom.name == 'C3':
                c3_atom = atom
            elif atom.name == 'O21':
                o21_atom = atom
        
        # Check if all required atoms are present
        if c1_atom is None or c2_atom is None or c3_atom is None or o21_atom is None:
            continue  # Skip this lipid if missing any required atoms
        
        # Verify that C1, C3, and O21 are all bonded to C2
        c2_neighbors = neighbor_map[c2_atom.index]['C'] + neighbor_map[c2_atom.index]['other']
        
        c1_bonded = c1_atom.index in c2_neighbors
        c3_bonded = c3_atom.index in c2_neighbors
        o21_bonded = o21_atom.index in c2_neighbors
        
        if not (c1_bonded and c3_bonded and o21_bonded):
            # Skip if not all atoms are bonded to C2
            continue
        
        # Create the improper dihedral tuple: (C1, C3, C2, O21)
        # This ordering maintains the tetrahedral geometry at C2:
        # - C2 is the central chiral atom
        # - C1, C3, O21 are the three substituents
        # - The improper angle should be ~120° for proper stereochemistry
        improper = (c1_atom.index, c3_atom.index, c2_atom.index, o21_atom.index)
        result['glycerol_impropers'].append(improper)
        impropers_found += 1
    
    print(f"Created {impropers_found} improper dihedral restraints for glycerol stereocenters.")
    print(f"  (Maintains chirality at C2 with C1-C3-C2-O21 improper dihedral)")
    
    return result

def create_lipid_proper_dihedral_restraints(system, cis_double_bond_atoms, force_constant):
    """
    Create dihedral restraints to maintain cis stereochemistry in lipid double bonds.
    
    Parameters
    ----------
    system : openmm.System
        The OpenMM system to add restraints to
    cis_double_bond_atoms : list of tuples
        List of 4-atom tuples (i, j, k, l) defining dihedral angles
    force_constant : openmm.unit.Quantity
        Force constant for dihedral restraints (energy units)
        
    Returns
    -------
    openmm.CustomTorsionForce or None
        The created force object, or None if no dihedrals provided
    """
    if not cis_double_bond_atoms or len(cis_double_bond_atoms) == 0:
        print("No cis double bond dihedrals to restrain.")
        return None
    
    print(f"Creating dihedral restraints for {len(cis_double_bond_atoms)} cis double bonds.")
    
    # Create custom torsion force with harmonic restraint around 0 degrees (cis)
    # Energy = k * theta^2 where theta is the dihedral angle
    force = openmm.CustomTorsionForce(f"{FORCE_CONST_NAME_DIHE_PROPER} * theta^2")
    force.addGlobalParameter(FORCE_CONST_NAME_DIHE_PROPER, force_constant)
    
    # Add each dihedral
    for atom_tuple in cis_double_bond_atoms:
        if len(atom_tuple) == 4:
            force.addTorsion(int(atom_tuple[0]), int(atom_tuple[1]), 
                           int(atom_tuple[2]), int(atom_tuple[3]), [])
    
    # Assign to force group 5 for energy tracking
    force.setForceGroup(5)
    
    # Add to system
    system.addForce(force)
    
    print(f"Added dihedral restraints with force constant: {force_constant}")
    
    return force

def create_lipid_improper_dihedral_restraints(system, glycerol_improper_atoms, force_constant):
    """
    Create improper dihedral restraints to maintain glycerol backbone stereochemistry.
    
    Restrains the improper dihedral at the C2 stereocenter to maintain proper
    R-configuration chirality. The improper dihedral is defined as (C1, C3, C2, O21)
    and is restrained to approximately 120 degrees.
    
    Parameters
    ----------
    system : openmm.System
        The OpenMM system to add restraints to
    glycerol_improper_atoms : list of tuples
        List of 4-atom tuples (C1, C3, C2, O21) defining improper dihedrals
    force_constant : openmm.unit.Quantity
        Force constant for improper dihedral restraints (energy units)
        
    Returns
    -------
    openmm.CustomTorsionForce or None
        The created force object, or None if no impropers provided
    """
    if not glycerol_improper_atoms or len(glycerol_improper_atoms) == 0:
        print("No glycerol improper dihedrals to restrain.")
        return None
    
    print(f"Creating improper dihedral restraints for {len(glycerol_improper_atoms)} glycerol stereocenters.")
    
    # Create custom torsion force with harmonic restraint around -120 degrees
    # Energy = k * (theta - theta0)^2 where theta0 = -120° = -2.0944 radians
    # Note: OpenMM torsion angles are in radians
    force = openmm.CustomTorsionForce(f"{FORCE_CONST_NAME_DIHE_IMPROPER} * (theta - theta0)^2")
    force.addGlobalParameter(FORCE_CONST_NAME_DIHE_IMPROPER, force_constant)
    force.addGlobalParameter("theta0", -2.0944)  # -120 degrees in radians
    
    # Add each improper dihedral
    for atom_tuple in glycerol_improper_atoms:
        if len(atom_tuple) == 4:
            # Order is (C1, C3, C2, O21) where C2 is the chiral center
            force.addTorsion(int(atom_tuple[0]), int(atom_tuple[1]), 
                           int(atom_tuple[2]), int(atom_tuple[3]), [])
    
    # Assign to force group 6 for energy tracking
    force.setForceGroup(6)
    
    # Add to system
    system.addForce(force)
    
    print(f"Added improper dihedral restraints with force constant: {force_constant}")
    print(f"  Target angle: -120° (maintains R-stereochemistry at C2)")
    
    return force

def create_lipid_phosphorus_restraints(system, phosphorus_atoms, reference_positions, force_constant):
    """
    Create positional restraints on lipid phosphorus atoms in Z-direction only.
    
    This allows lipids to diffuse laterally (XY plane) while maintaining bilayer thickness.
    
    Parameters
    ----------
    system : openmm.System
        The OpenMM system to add restraints to
    phosphorus_atoms : list
        List of phosphorus atom indices
    reference_positions : np.ndarray or list
        Reference positions for restraints (shape: [n_atoms, 3])
    force_constant : openmm.unit.Quantity
        Force constant for restraints (energy/length^2)
        
    Returns
    -------
    openmm.CustomExternalForce or None
        The created force object, or None if no atoms provided
    """
    if not phosphorus_atoms or len(phosphorus_atoms) == 0:
        print("No phosphorus atoms to restrain.")
        return None
    
    print(f"Creating Z-only positional restraints for {len(phosphorus_atoms)} phosphorus atoms.")
    
    # Create custom external force with harmonic restraint only in Z direction
    # This allows lateral diffusion while maintaining bilayer structure
    expr = f"0.5*{FORCE_CONST_NAME_PHOS}*periodicdistance(0, 0, z, 0, 0, z0)^2"
    force = openmm.CustomExternalForce(expr)
    force.addGlobalParameter(FORCE_CONST_NAME_PHOS, force_constant)
    force.addPerParticleParameter("z0")
    
    # Add each phosphorus atom with its reference Z position
    for p_idx in phosphorus_atoms:
        position = reference_positions[p_idx]
        # Only store z0, since we only restrain in Z
        force.addParticle(int(p_idx), [position[2]])
    
    # Assign to force group 7 for energy tracking
    force.setForceGroup(7)
    
    # Add to system
    system.addForce(force)
    
    print(f"Added Z-only phosphorus restraints with force constant: {force_constant}")
    
    return force

def select_both_from_mdtraj_and_indices(solv_traj, solv_indices, unsolv_traj, 
                                        unsolv_indices, selection):
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
    filtered_solv_indices = []
    filtered_unsolv_indices = []
    solv_selected_indices = solv_traj.topology.select(selection)
    unsolv_selected_indices = unsolv_traj.topology.select(selection)
    print(f"MDTraj solvated_selection '{selection}' returned {len(solv_selected_indices)} atoms.")
    print(f"MDTraj unsolvated_selection '{selection}' returned {len(unsolv_selected_indices)} atoms.")
    for solv_idx, unsolv_idx in zip(solv_indices, unsolv_indices):
        if solv_idx in solv_selected_indices and unsolv_idx in unsolv_selected_indices:
            filtered_solv_indices.append(solv_idx)
            filtered_unsolv_indices.append(unsolv_idx)
        
    return filtered_solv_indices, filtered_unsolv_indices

def create_positional_restraints(
        system,
        unsolvated_mdtraj,
        solvated_atom_indices,
        unsolvated_atom_indices,
        restraint_force_constant,
        force_const_name
        ):
    print("Imposing positional restraints to defined atoms.")
    reference_positions = unsolvated_mdtraj.xyz[0,:,:]
    expr = f"0.5*{force_const_name}*periodicdistance(x, y, z, x0, y0, z0)^2"
    force = openmm.CustomExternalForce(expr)
    force.addGlobalParameter(force_const_name, restraint_force_constant)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    for i, (solv_index, unsolv_index) in enumerate(zip(solvated_atom_indices, unsolvated_atom_indices)):
        position = reference_positions[unsolv_index]
        force.addParticle(solv_index, [position[0], position[1], position[2]])

    system.addForce(force)
    return
    
def create_new_simulation_objects(
        topology,
        system, 
        start_temp,
        temperature_freq=1.0,
        constant_pressure=False,
        target_pressure=1.0*unit.bar,
        pressure_frequency=25,
        time_step=0.002*unit.picoseconds,
        cuda_index="0",
        use_membrane_barostat=False,
        surface_tension=0*unit.bar*unit.nanometer,
        random_seed=DEFAULT_RANDOM_SEED
        ):
    # Set up force groups for detailed energy reporting
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        force_name = force.__class__.__name__
        
        if 'Bond' in force_name:
            force.setForceGroup(1)  # Bond
        elif 'Angle' in force_name:
            force.setForceGroup(2)  # Angle  
        elif 'Torsion' in force_name or 'Dihedral' in force_name:
            force.setForceGroup(3)  # Dihedral
        elif 'Nonbonded' in force_name or 'Electrostatic' in force_name or 'LJ' in force_name:
            force.setForceGroup(4)  # Nonbonded
        # CustomExternalForce (restraints) will stay in default group

    if constant_pressure:
        if use_membrane_barostat:
            print("Using membrane barostat (anisotropic pressure control).")
            barostat = openmm.MonteCarloMembraneBarostat(
                target_pressure, 
                surface_tension, 
                start_temp,
                openmm.MonteCarloMembraneBarostat.XYIsotropic,
                openmm.MonteCarloMembraneBarostat.ZFree,
                pressure_frequency
            )
        else:
            barostat = openmm.MonteCarloBarostat(target_pressure, start_temp, pressure_frequency)
        barostat.setRandomNumberSeed(random_seed)
        system.addForce(barostat)
    integrator = openmm.LangevinIntegrator(start_temp, temperature_freq/unit.picosecond, time_step)
    integrator.setRandomNumberSeed(random_seed)
    platform = openmm.Platform.getPlatformByName('CUDA')
    properties = {"CudaDeviceIndex": cuda_index, "CudaPrecision": "mixed"}
    simulation = openmm_app.Simulation(topology, system, integrator, platform, properties)
    return simulation

def main(
        system: openmm.System,
        solvated_pdb_filename: str,
        unsolvated_pdb_filename: str,
        temperature: unit.Quantity,
        target_pressure: unit.Quantity,
        is_membrane_system: bool = False,
        surface_tension: unit.Quantity = 0*unit.bar*unit.nanometer,
        cuda_index: str = "0",
        working_directory: str = "minimization_equilibration",
        stage_1_steps_per_ramp: int = DEFAULT_STAGE_1_STEPS_PER_RAMP,
        stage_2_total_steps: int = DEFAULT_STAGE_2_TOTAL_STEPS,
        stage_3_steps_per_ramp: int = DEFAULT_STAGE_3_STEPS_PER_RAMP,
        stage_4_total_steps: int = DEFAULT_STAGE_4_TOTAL_STEPS,
        stage_5_total_steps: int = DEFAULT_STAGE_5_TOTAL_STEPS,
        stage_6_total_steps: int = DEFAULT_STAGE_6_TOTAL_STEPS,
        stage_7_total_steps: int = DEFAULT_STAGE_7_TOTAL_STEPS,
        stage_8_total_steps: int = DEFAULT_STAGE_8_TOTAL_STEPS,
        stage_1_timestep: int = DEFAULT_STAGE_1_TIMESTEP,
        stage_2_timestep: int = DEFAULT_STAGE_2_TIMESTEP,
        stage_3_timestep: int = DEFAULT_STAGE_3_TIMESTEP,
        stage_4_timestep: int = DEFAULT_STAGE_4_TIMESTEP,
        stage_5_timestep: int = DEFAULT_STAGE_5_TIMESTEP,
        stage_6_timestep: int = DEFAULT_STAGE_6_TIMESTEP,
        stage_7_timestep: int = DEFAULT_STAGE_7_TIMESTEP,
        stage_8_timestep: int = DEFAULT_STAGE_8_TIMESTEP,
        eq_traj_interval: int = DEFAULT_EQ_TRAJ_INTERVAL,
        eq_traj_filename: str = DEFAULT_EQ_TRAJ_FILENAME,
        eq_state_interval: int = DEFAULT_EQ_STATE_INTERVAL,
        eq_state_filename: str = DEFAULT_EQ_STATE_FILENAME,
        random_seed: int = DEFAULT_RANDOM_SEED
        # TODO: support the ability to change from force const defaults? ...
        ) -> str:

    solvated_pdb = openmm_app.PDBFile(solvated_pdb_filename)
    topology = solvated_pdb.topology
    solvated_mdtraj = mdtraj.load(solvated_pdb_filename)
    solvated_mdtraj.topology = mdtraj.Topology.from_openmm(topology)
    unsolvated_mdtraj = mdtraj.load(unsolvated_pdb_filename)
    startdir = os.path.abspath(os.getcwd())
    os.chdir(working_directory)
    # Detect lipid atoms if this is a membrane system
    lipid_dihedral_k = DEFAULT_STAGE_1_K_DIHE
    lipid_p_restraint_k = DEFAULT_STAGE_1_K_HEAD
    lipid_atoms = {'phosphorus': [], 'headgroup_heavy': [], 'cis_double_bonds': [],
                   'glycerol_impropers': []}
    if is_membrane_system:
        print("\n" + "="*60)
        print("MEMBRANE SYSTEM DETECTED - Identifying lipid atoms")
        print("="*60)
        lipid_atoms = identify_lipid_atoms_with_mdtraj(solvated_mdtraj)
        print(f"Lipid detection summary:")
        print(f"  - Phosphorus atoms: {len(lipid_atoms['phosphorus'])}")
        print(f"  - Headgroup heavy atoms: {len(lipid_atoms['headgroup_heavy'])}")
        print(f"  - Cis double bonds: {len(lipid_atoms['cis_double_bonds'])}")
        print(f"  - Glycerol impropers: {len(lipid_atoms['glycerol_impropers'])}")
        print("="*60 + "\n")
        
        # Create dihedral restraints for cis double bonds if any were detected
        if len(lipid_atoms['cis_double_bonds']) > 0:
            create_lipid_proper_dihedral_restraints(
                system, 
                lipid_atoms['cis_double_bonds'], 
                lipid_dihedral_k
            )
        
        # Create improper dihedral restraints for glycerol stereocenters
        if len(lipid_atoms['glycerol_impropers']) > 0:
            create_lipid_improper_dihedral_restraints(
                system,
                lipid_atoms['glycerol_impropers'],
                lipid_dihedral_k  # Using same force constant as proper dihedrals for now
            )
        # Create positional restraints for phosphorus atoms
        if len(lipid_atoms['phosphorus']) > 0:
            create_lipid_phosphorus_restraints(
                system,
                lipid_atoms['phosphorus'],
                solvated_mdtraj.xyz[0, :, :],
                lipid_p_restraint_k
            )
        
    map_unsolv_to_solv_list = find_known_atom_indices(solvated_mdtraj, unsolvated_mdtraj)
    unzipped_mapped_lists = list(zip(*map_unsolv_to_solv_list))
    unsolvated_atom_indices = unzipped_mapped_lists[0]
    solvated_atom_indices = unzipped_mapped_lists[1]
    assert len(unsolvated_atom_indices) == len(solvated_atom_indices), \
        "Number of mapped atoms between solvated and unsolvated structures should be the same"
    protein_backbone_solvated_indices, protein_backbone_unsolvated_indices \
        = select_both_from_mdtraj_and_indices(
            solvated_mdtraj, solvated_atom_indices, unsolvated_mdtraj, unsolvated_atom_indices, 
            "protein and backbone and not element H"
        )
    assert len(protein_backbone_unsolvated_indices) == len(protein_backbone_solvated_indices), \
        "Number of protein backbone atoms in solvated and unsolvated structures "\
        "should be the same"
    protein_sidechain_solvated_indices, protein_sidechain_unsolvated_indices \
        = select_both_from_mdtraj_and_indices(
            solvated_mdtraj, solvated_atom_indices, unsolvated_mdtraj, unsolvated_atom_indices, 
            "protein and sidechain and not element H"
        )
    assert len(protein_sidechain_unsolvated_indices) == len(protein_sidechain_solvated_indices), \
        "Number of protein sidechain atoms in solvated and unsolvated structures should be "\
        "the same"
    nonprotein_solvated_indices, nonprotein_unsolvated_indices \
        = select_both_from_mdtraj_and_indices(
            solvated_mdtraj, solvated_atom_indices, unsolvated_mdtraj, unsolvated_atom_indices, 
            "not protein and not element H"
        )
    assert len(nonprotein_unsolvated_indices) == len(nonprotein_solvated_indices), \
        "Number of protein sidechain atoms in solvated and unsolvated structures should be "\
        "the same"
    #print("nonprotein_solvated_indices:", nonprotein_solvated_indices)
    #print("nonprotein_unsolvated_indices:", nonprotein_unsolvated_indices)
    # Loop through the solvated and unsolvated atoms in pair using zip, and
    # assert that elements are the same, at least
    for solv_idx, unsolv_idx in zip(nonprotein_solvated_indices, nonprotein_unsolvated_indices):
        solv_atom = solvated_mdtraj.topology.atom(solv_idx)
        unsolv_atom = unsolvated_mdtraj.topology.atom(unsolv_idx)
        if solv_atom.element != unsolv_atom.element:
            print("ERROR: Non-matching elements found between solvated and unsolvated "
                  "structures for non-protein atom:")
            print(f"  Solvated atom index: {solv_idx}, name: {solv_atom.name}, "
                  f"element: {solv_atom.element}")
            print(f"  Unsolvated atom index: {unsolv_idx}, name: {unsolv_atom.name}, "
                  f"element: {unsolv_atom.element}")
        
    protein_backbone_restraint_k = DEFAULT_STAGE_1_K_BACKBONE
    if len(protein_backbone_solvated_indices) > 0:
        print(f"Imposing positional restraints on {len(protein_backbone_solvated_indices)} "
              "protein backbone atoms.")
        create_positional_restraints(
            system,
            unsolvated_mdtraj,
            protein_backbone_solvated_indices,
            protein_backbone_unsolvated_indices,
            protein_backbone_restraint_k,
            FORCE_CONST_NAME_BACKBONE
        )
    protein_sidechain_restraint_k = DEFAULT_STAGE_1_K_SIDECHAIN
    if len(protein_sidechain_solvated_indices) > 0:
        print(f"Imposing positional restraints on {len(protein_sidechain_solvated_indices)} "
              "protein sidechain atoms.")
        create_positional_restraints(
            system,
            unsolvated_mdtraj,
            protein_sidechain_solvated_indices,
            protein_sidechain_unsolvated_indices,
            protein_sidechain_restraint_k,
            FORCE_CONST_NAME_SIDECHAIN
        )
    if len(nonprotein_solvated_indices) > 0:
        print(f"Imposing positional restraints on {len(nonprotein_solvated_indices)} "
              "non-protein atoms.")
        create_positional_restraints(
            system,
            unsolvated_mdtraj,
            nonprotein_solvated_indices,
            nonprotein_unsolvated_indices,
            DEFAULT_STAGE_1_K_NONPROTEIN,
            FORCE_CONST_NAME_NONPROTEIN
        )
    
    box_vectors_eq1 = topology.getPeriodicBoxVectors()
    #print("box_vectors_eq1:", box_vectors_eq1)
    start_temp = DEFAULT_STAGE_1_START_TEMP
    time_step = stage_1_timestep
    simulation = create_new_simulation_objects(
        topology,
        system, 
        start_temp,
        temperature_freq=0.1,
        constant_pressure=True,
        target_pressure=target_pressure,
        pressure_frequency=25,
        time_step=time_step,
        cuda_index=cuda_index,
        use_membrane_barostat=is_membrane_system,
        surface_tension=surface_tension,
        random_seed=random_seed
    )

    integrator = simulation.context.getIntegrator()
    simulation.context.setPositions(solvated_pdb.positions)
    simulation.context.setPeriodicBoxVectors(*box_vectors_eq1)
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(
        start_temp, random_seed)
    #simulation.reporters.append(openmm_app.StateDataReporter(
    #     stdout, eq1_traj_interval, step=True,
    #     potentialEnergy=True, temperature=True, volume=True))
    if is_membrane_system:
        detailed_reporter_eq1 = DetailedEnergyReporter(
            open(eq_state_filename, 'w'), eq_state_interval, lipid_atoms['phosphorus'])
    else:
        detailed_reporter_eq1 = DetailedEnergyReporter(
            open(eq_state_filename, 'w'), eq_state_interval)
    simulation.reporters.append(detailed_reporter_eq1)
    #eq_pdb_reporter = openmm_app.PDBReporter(
    #    eq_traj_filename, eq_traj_interval)
    #simulation.reporters.append(eq_pdb_reporter)
    eq_dcd_reporter = openmm_app.DCDReporter(
        eq_traj_filename, eq_traj_interval)
    simulation.reporters.append(eq_dcd_reporter)
    temperature_ramp = np.arange(start_temp.value_in_unit(unit.kelvin),
                                 temperature.value_in_unit(unit.kelvin)+1,
                                 10.0)
    print("Starting eq1 simulation with temperature ramp.")
    steps_per_ramp = stage_1_steps_per_ramp
    for temp in temperature_ramp:
        print("Equilibrating at temperature: {:.2f} K".format(temp))
        #if temp == 190.0:
        #    print("Adding dcdreporter to identify numerical issue")
        #    dcd_reporter = openmm_app.DCDReporter("debug.dcd", 1)
        #    simulation.reporters.append(dcd_reporter)
        simulation.context.setVelocitiesToTemperature(
            temp*unit.kelvin, random_seed)
        integrator.setTemperature(temp*unit.kelvin)
        simulation.step(steps_per_ramp)
        
    print("Finished eq1 simulation.")
    
    simulation.context.setParameter(FORCE_CONST_NAME_BACKBONE, DEFAULT_STAGE_2_K_BACKBONE)
    simulation.context.setParameter(FORCE_CONST_NAME_SIDECHAIN, DEFAULT_STAGE_2_K_SIDECHAIN)
    simulation.context.setParameter(FORCE_CONST_NAME_NONPROTEIN, DEFAULT_STAGE_2_K_NONPROTEIN)
    if is_membrane_system:
        simulation.context.setParameter(FORCE_CONST_NAME_DIHE_PROPER, DEFAULT_STAGE_2_K_DIHE)
        simulation.context.setParameter(FORCE_CONST_NAME_DIHE_IMPROPER, DEFAULT_STAGE_2_K_DIHE)
        simulation.context.setParameter(FORCE_CONST_NAME_PHOS, DEFAULT_STAGE_2_K_HEAD)

    integrator.setStepSize(stage_2_timestep)
    # Equilibration 2: Run at constant P and T long enough to reach plateau in
    # energies and volumes.
    print(f"Starting eq2 simulation for {stage_2_total_steps} steps.")
    simulation.step(stage_2_total_steps)


    simulation.context.setParameter(FORCE_CONST_NAME_BACKBONE, DEFAULT_STAGE_3_K_BACKBONE)
    simulation.context.setParameter(FORCE_CONST_NAME_SIDECHAIN, DEFAULT_STAGE_3_K_SIDECHAIN)
    simulation.context.setParameter(FORCE_CONST_NAME_NONPROTEIN, DEFAULT_STAGE_3_K_NONPROTEIN)
    if is_membrane_system:
        simulation.context.setParameter(FORCE_CONST_NAME_DIHE_PROPER, DEFAULT_STAGE_3_K_DIHE)
        simulation.context.setParameter(FORCE_CONST_NAME_DIHE_IMPROPER, DEFAULT_STAGE_3_K_DIHE)
        simulation.context.setParameter(FORCE_CONST_NAME_PHOS, DEFAULT_STAGE_3_K_HEAD)
    integrator.setStepSize(stage_3_timestep)
    print("Starting eq3 simulation with temperature ramp.")
    steps_per_ramp = stage_3_steps_per_ramp
    for temp in temperature_ramp:
        print("Equilibrating at temperature: {:.2f} K".format(temp))
        simulation.context.setVelocitiesToTemperature(
            temp*unit.kelvin, random_seed)
        integrator.setTemperature(temp*unit.kelvin)
        simulation.step(steps_per_ramp)


    simulation.context.setParameter(FORCE_CONST_NAME_BACKBONE, DEFAULT_STAGE_4_K_BACKBONE)
    simulation.context.setParameter(FORCE_CONST_NAME_SIDECHAIN, DEFAULT_STAGE_4_K_SIDECHAIN)
    simulation.context.setParameter(FORCE_CONST_NAME_NONPROTEIN, DEFAULT_STAGE_4_K_NONPROTEIN)
    if is_membrane_system:
        simulation.context.setParameter(FORCE_CONST_NAME_DIHE_PROPER, DEFAULT_STAGE_4_K_DIHE)
        simulation.context.setParameter(FORCE_CONST_NAME_DIHE_IMPROPER, DEFAULT_STAGE_4_K_DIHE)
        simulation.context.setParameter(FORCE_CONST_NAME_PHOS, DEFAULT_STAGE_4_K_HEAD)
    integrator.setStepSize(stage_4_timestep)
    print(f"Starting eq4 simulation for {stage_4_total_steps} steps.")
    simulation.step(stage_4_total_steps)


    simulation.context.setParameter(FORCE_CONST_NAME_BACKBONE, DEFAULT_STAGE_5_K_BACKBONE)
    simulation.context.setParameter(FORCE_CONST_NAME_SIDECHAIN, DEFAULT_STAGE_5_K_SIDECHAIN)
    simulation.context.setParameter(FORCE_CONST_NAME_NONPROTEIN, DEFAULT_STAGE_5_K_NONPROTEIN)
    if is_membrane_system:
        simulation.context.setParameter(FORCE_CONST_NAME_DIHE_PROPER, DEFAULT_STAGE_5_K_DIHE)
        simulation.context.setParameter(FORCE_CONST_NAME_DIHE_IMPROPER, DEFAULT_STAGE_5_K_DIHE)
        simulation.context.setParameter(FORCE_CONST_NAME_PHOS, DEFAULT_STAGE_5_K_HEAD)
    integrator.setStepSize(stage_5_timestep)
    print(f"Starting eq5 simulation for {stage_5_total_steps} steps.")
    simulation.step(stage_5_total_steps)


    simulation.context.setParameter(FORCE_CONST_NAME_BACKBONE, DEFAULT_STAGE_6_K_BACKBONE)
    simulation.context.setParameter(FORCE_CONST_NAME_SIDECHAIN, DEFAULT_STAGE_6_K_SIDECHAIN)
    simulation.context.setParameter(FORCE_CONST_NAME_NONPROTEIN, DEFAULT_STAGE_6_K_NONPROTEIN)
    if is_membrane_system:
        simulation.context.setParameter(FORCE_CONST_NAME_DIHE_PROPER, DEFAULT_STAGE_6_K_DIHE)
        simulation.context.setParameter(FORCE_CONST_NAME_DIHE_IMPROPER, DEFAULT_STAGE_6_K_DIHE)
        simulation.context.setParameter(FORCE_CONST_NAME_PHOS, DEFAULT_STAGE_6_K_HEAD)
    integrator.setStepSize(stage_6_timestep)
    print(f"Starting eq6 simulation for {stage_6_total_steps} steps.")
    simulation.step(stage_6_total_steps)


    simulation.context.setParameter(FORCE_CONST_NAME_BACKBONE, DEFAULT_STAGE_7_K_BACKBONE)
    simulation.context.setParameter(FORCE_CONST_NAME_SIDECHAIN, DEFAULT_STAGE_7_K_SIDECHAIN)
    simulation.context.setParameter(FORCE_CONST_NAME_NONPROTEIN, DEFAULT_STAGE_7_K_NONPROTEIN)
    if is_membrane_system:
        simulation.context.setParameter(FORCE_CONST_NAME_DIHE_PROPER, DEFAULT_STAGE_7_K_DIHE)
        simulation.context.setParameter(FORCE_CONST_NAME_DIHE_IMPROPER, DEFAULT_STAGE_7_K_DIHE)
        simulation.context.setParameter(FORCE_CONST_NAME_PHOS, DEFAULT_STAGE_7_K_HEAD)
    integrator.setStepSize(stage_7_timestep)
    print(f"Starting eq7 simulation for {stage_7_total_steps} steps.")
    simulation.step(stage_7_total_steps)


    simulation.context.setParameter(FORCE_CONST_NAME_BACKBONE, DEFAULT_STAGE_8_K_BACKBONE)
    simulation.context.setParameter(FORCE_CONST_NAME_SIDECHAIN, DEFAULT_STAGE_8_K_SIDECHAIN)
    simulation.context.setParameter(FORCE_CONST_NAME_NONPROTEIN, DEFAULT_STAGE_8_K_NONPROTEIN)
    if is_membrane_system:
        simulation.context.setParameter(FORCE_CONST_NAME_DIHE_PROPER, DEFAULT_STAGE_8_K_DIHE)
        simulation.context.setParameter(FORCE_CONST_NAME_DIHE_IMPROPER, DEFAULT_STAGE_8_K_DIHE)
        simulation.context.setParameter(FORCE_CONST_NAME_PHOS, DEFAULT_STAGE_8_K_HEAD)
    integrator.setStepSize(stage_8_timestep)
    print(f"Starting eq8 simulation for {stage_8_total_steps} steps.")
    simulation.step(stage_8_total_steps)

    # Final structure output
    state = simulation.context.getState(getPositions=True, getVelocities=True,
                                        enforcePeriodicBox=True)
    box_vectors = state.getPeriodicBoxVectors()
    assert simulation.topology is not None
    simulation.topology.setPeriodicBoxVectors(box_vectors)
    #with openmm_app.PDBFile(output_pdb_file, simulation.topology) as output_pdb:
    #    output_pdb.writeFile(simulation.topology, state.getPositions(), file=stdout)
    print(f"Writing final equilibrated structure to {OUTPUT_PDB_FILE}.")
    with open(OUTPUT_PDB_FILE, "w") as out_file:
            openmm_app.PDBFile.writeHeader(simulation.topology, out_file)
            openmm_app.PDBFile.writeModel(simulation.topology, state.getPositions(), out_file)
            openmm_app.PDBFile.writeFooter(simulation.topology, out_file)
    
    equilibrated = mdtraj.load(OUTPUT_PDB_FILE)
    equilibrated.topology = solvated_mdtraj.topology
    #anchor_molecules = []
    #for index in receptor_indices:
    #    anchor_molecules.append(mdtraj_obj.topology.atom(index))

    #equilibrated.image_molecules(inplace=True, anchor_molecules=[set(anchor_molecules)])
    equilibrated.image_molecules(inplace=True)
    equilibrated.save(OUTPUT_PDB_FILE_IMAGED)
    equilibrated_pdb_filename = os.path.join(working_directory, OUTPUT_PDB_FILE_IMAGED)
    os.chdir(startdir)
    return equilibrated_pdb_filename