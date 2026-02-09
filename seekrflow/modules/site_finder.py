"""
Determine a good site definition within a protein structure for
a seekr calculation, based on a selection of small molecule ligand
indices.
"""

import numpy as np
import mdtraj

# If the number of site atoms falls significantly
TARGET_SITE_ATOMS = 20
TARGET_SITE_ATOMS_ENERGY_PENALTY = 0.1

# Monte-carlo temparture and energy penalty for distance from the ligand COM
MC_TEMP = 0.1
MC_DIST_ENERGY_PENALTY = 12.0

MAX_ITER = 10000
# The desired depth (in nm) of the center of the site behind the ligand
#  if the location of an opening is provided
DEEPER_THAN_LIGAND = 0.2

# If the selected site COM falls within this tolerance of the desired
#  location (in nm), the adjustments will stop.
TOLERANCE = 0.05

INNER_DISTANCE_CUTOFF = 0.5 # nm
OUTER_DISTANCE_CUTOFF = 0.8 # nm

def atom_selection_within_cutoff(
        traj: mdtraj.Trajectory, 
        ligand_indices: list[int], 
        distance_cutoff: float, 
        ligand_resname: str = ""
        ) -> list[int]:
    """
    Given a MDtraj structure, ligand indices, and a cutoff, return all
    alpha carbon indices within the cutoff distance of any ligand
    atoms
    """
    contact_prot_indices = []
    contact_resids = set()
    if ligand_resname == "":
        protein_indices = traj.topology.select(f"protein")
    else:
        protein_indices = traj.topology.select(f"protein and not resname {ligand_resname}")
    
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
        ca_index = traj.topology.select(ca_selection_str)[0]
        alpha_carbon_indices.append(ca_index)
    
    alpha_carbon_indices.sort()
    return alpha_carbon_indices

def traj_center_of_mass(traj):
    if traj.xyz.shape[1] == 1:
        trying_COM = traj.xyz[:,0,:]
    else:
        trying_COM = mdtraj.compute_center_of_mass(traj)
    return trying_COM

def try_ca_set(traj, trying_ca_set, desired_COM):
    """
    Given a MDtraj object and a set of alpha carbon indices, 
    compute the COM of those atoms and return the distance from the desired COM.
    """
    trying_traj = traj.atom_slice(trying_ca_set)
    trying_COM = traj_center_of_mass(trying_traj)
    trying_dist = np.linalg.norm(trying_COM - desired_COM)
    return trying_dist, trying_COM

def adjust_to_tolerance(traj, ligand_indices, current_atoms, possible_atoms,
                        site_opening_indices=None):
    """
    Adjust the atom selection until the site is within tolerance of the
    desired location.
    """
    lig_traj = traj.atom_slice(ligand_indices)
    lig_COM = traj_center_of_mass(lig_traj)
    current_traj = traj.atom_slice(current_atoms)
    current_COM = traj_center_of_mass(current_traj)
    if site_opening_indices is None:
        desired_COM = lig_COM
    else:
        opening_traj = traj.atom_slice(site_opening_indices)
        opening_COM = traj_center_of_mass(opening_traj)
        lig_to_opening = opening_COM - lig_COM
        desired_COM = lig_COM - DEEPER_THAN_LIGAND * lig_to_opening \
            / np.linalg.norm(lig_to_opening)
    
    # Find the alpha carbons that can be added
    indices_that_can_be_added = []
    for atom_index in possible_atoms:
        if atom_index not in current_atoms:
            indices_that_can_be_added.append(atom_index)
    
    new_atom_indices = current_atoms[:]
    found_good_index_set = False
    
    counter = 0
    current_dist = np.linalg.norm(current_COM - desired_COM)
    number_energy = TARGET_SITE_ATOMS_ENERGY_PENALTY * (len(new_atom_indices) - TARGET_SITE_ATOMS)**2
    dist_energy = MC_DIST_ENERGY_PENALTY * current_dist**2
    total_energy = number_energy + dist_energy
    while not found_good_index_set:
        #print("iteration:", counter, "distance to optimal:", current_dist)
        # Try adding one - choose a random index from indices_that_can_be_added
        trying_atom_index = new_atom_indices[0]
        if len(new_atom_indices) == len(possible_atoms):
            raise Exception("No more atoms to add.")
        
        while trying_atom_index in new_atom_indices:
            trying_atom_index = int(np.random.choice(indices_that_can_be_added))
        trying_atom_index_set = new_atom_indices[:] + [trying_atom_index]
        trying_atom_index_set.sort()
        trying_dist, trying_COM = try_ca_set(traj, trying_atom_index_set, desired_COM)
        
        number_energy = TARGET_SITE_ATOMS_ENERGY_PENALTY * (len(trying_atom_index_set) - TARGET_SITE_ATOMS)**2
        dist_energy = MC_DIST_ENERGY_PENALTY * trying_dist**2
        trying_energy = number_energy + dist_energy
        if trying_energy < total_energy:
            #print("adding atom (MC accepted):", trying_atom_index)
            new_atom_indices.append(trying_atom_index)
            indices_that_can_be_added.remove(trying_atom_index)
            total_energy = trying_energy
            current_dist = trying_dist
        else:
            # Run monte carlo acceptance criterion
            energy_diff = trying_energy - total_energy
            acceptance_prob = np.exp(-energy_diff / MC_TEMP)
            if np.random.rand() < acceptance_prob:
                #print("adding atom (MC accepted):", trying_atom_index)
                new_atom_indices.append(trying_atom_index)
                indices_that_can_be_added.remove(trying_atom_index)
                total_energy = trying_energy
                current_dist = trying_dist
        
        new_atom_indices.sort()
        if current_dist < TOLERANCE:
            #print("Close selection found by adding atom:", trying_atom_index)
            found_good_index_set = True
            break
        
        # Try removing one
        if len(new_atom_indices) == 0:
            raise Exception("No more atoms to remove.")
        
        popping_list_index = np.random.choice(len(new_atom_indices))
        popping_atom_index = new_atom_indices[popping_list_index]
        trying_atom_index_set = new_atom_indices[:]
        #print("trying_atom_index_set:", trying_atom_index_set)
        trying_atom_index_set.remove(popping_atom_index)
        trying_dist, trying_COM = try_ca_set(traj, trying_atom_index_set, desired_COM)
        
        number_energy = TARGET_SITE_ATOMS_ENERGY_PENALTY * (len(trying_atom_index_set) - TARGET_SITE_ATOMS)**2
        dist_energy = MC_DIST_ENERGY_PENALTY * trying_dist**2
        trying_energy = number_energy + dist_energy
        if trying_energy < total_energy:
            #print("removing atom (MC accepted):", popping_atom_index)
            new_atom_indices.remove(popping_atom_index)
            indices_that_can_be_added.append(popping_atom_index)
            total_energy = trying_energy
            current_dist = trying_dist
        else:
            # Run monte carlo acceptance criterion
            energy_diff = trying_energy - total_energy
            acceptance_prob = np.exp(-energy_diff / MC_TEMP)
            if np.random.rand() < acceptance_prob:
                #print("removing atom (MC accepted):", popping_atom_index)
                new_atom_indices.remove(popping_atom_index)
                indices_that_can_be_added.append(popping_atom_index)
                total_energy = trying_energy
                current_dist = trying_dist
        
        if current_dist < TOLERANCE:
            #print("Close selection found by removing atom:", popping_atom_index)
            found_good_index_set = True
            break

        counter += 1
        if counter > MAX_ITER:
            raise Exception("Maximum iterations exceeded.")
    
    #print("Found good ca set:", new_ca)
    final_dist, final_COM = try_ca_set(traj, new_atom_indices, desired_COM)
    #print("distance to optimal:", final_dist, "nm") 
    #print("new_atom_indices:", new_atom_indices)
    #print("final_dist:", final_dist, "nm")
    new_atom_indices = [int(index) for index in new_atom_indices]
    return new_atom_indices, final_COM

def site_finder_monte_carlo(
        pdb_filename: str, 
        ligand_indices: list[int],):
    """
    Given a PDB filename and ligand indices, perform a monte carlo search 
    to find a good site definition. The site definition will be a set of 
    alpha carbons within a certain distance of the ligand, and the COM of 
    those alpha carbons will be close to the desired location (either the 
    ligand COM or a location behind the ligand if an opening is provided).
    """
    traj = mdtraj.load(pdb_filename)
    current_atoms = atom_selection_within_cutoff(
        traj, ligand_indices, INNER_DISTANCE_CUTOFF)
    possible_atoms = atom_selection_within_cutoff(
        traj, ligand_indices, OUTER_DISTANCE_CUTOFF)
    new_atom_indices, final_COM = adjust_to_tolerance(
        traj, ligand_indices, current_atoms, possible_atoms)
    return new_atom_indices