"""
modules/base.py

Base routines for seekrflow.
"""

import os
import typing

from attrs import define, field, validators
import mdtraj
import parmed

PME_TOL = 2.5e-4
APBS_GRID_SPACING = 0.5
BENCHMARK_MIN_SIMULATION_LENGTH = 100000 # steps

@define
class Physical_attributes:
    """
    Physical attributes for the system.
    """
    temperature: float = field(
        default=300.0,
        validator=validators.instance_of(float),
        )
    pressure: float | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(float)),
        )
    ionic_strength: float = field(
        default=0.15,
        validator=validators.instance_of(float),
        )
    hydrogen_mass: float = field(
        default=1.008,
        validator=validators.instance_of(float),
        )

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

def assign_mbondi2_radii_to_parmed_structure(
        parmed_structure: parmed.Structure
        ) -> None:
    """
    Assign MBondi2 radii to a ParmEd structure.
    """
    # Assign mbondi2 radii to atoms
    for atom in parmed_structure.atoms:
        if atom.element == 1:
            for bonded_atom in atom.bond_partners:
                if bonded_atom.element == 7:
                    atom.solvent_radius = 1.3
                else:
                    atom.solvent_radius = 1.2
        elif atom.element == 6:
            atom.solvent_radius = 1.7
        elif atom.element == 7:
            atom.solvent_radius = 1.55
        elif atom.element == 8:
            atom.solvent_radius = 1.5
        elif atom.element == 9:
            atom.solvent_radius = 1.5
        elif atom.element == 14:
            atom.solvent_radius = 2.1
        elif atom.element == 15:
            atom.solvent_radius = 1.85
        elif atom.element == 16:
            atom.solvent_radius = 1.8
        elif atom.element == 17:
            atom.solvent_radius = 1.7
        else:
            atom.solvent_radius = 1.7
    return

def choose_only_receptor_atoms(
        complex_structure: parmed.Structure, 
        receptor_structure: parmed.Structure
        ) -> parmed.Structure:
    """
    If the complex_structure file contains the atomic charges and radii,
    conversely the receptor_structure file contains all atoms of the receptor,
    then select from complex_structure the atoms that will be present in the
    receptor_structure file and place them into a new parmed Structure
    """
    chain_ids = set()
    for atom in receptor_structure.atoms:
        if atom.residue.chain is not None:
            chain = atom.residue.chain
            chain_ids.add(chain)
    receptor = complex_structure[",".join(chain_ids), :, :]
    return receptor

def make_protein_openmm_and_mdtraj_top(protein_pdb_filename: str) -> typing.Tuple[
        "openmm.app.Topology",
        mdtraj.Topology,
        typing.List["openmm.unit.Quantity"]
        ]:
    """
    Given a PDB filename, return the OpenMM Topology, MDTraj Topology, and
    list of positions.
    """
    import openmm.app as openmm_app
    with open(protein_pdb_filename, 'r') as f:
        protein = openmm_app.PDBFile(f) 
    protein_topology = protein.topology
    protein_positions = protein.positions 
    protein_md_topology = mdtraj.Topology.from_openmm(protein_topology)
    return (
        protein_topology,
        protein_md_topology,
        protein_positions
    )

def make_ligand_openmm_and_mdtraj_top(
        ligand_sdf_filename: str,
        ligand_pdb_filename: str,
        param_directory: str,
        draw_ligand: bool = True,
        ligand_resname: str = ""
        ) -> typing.Tuple[
        "openmm.app.Topology",
        mdtraj.Topology,
        typing.List["openmm.unit.Quantity"],
        "openff.toolkit.topology.Molecule"
        ]:
    """
    Given a ligand SDF filename, return the OpenMM Topology, MDTraj Topology, and
    list of positions.
    """
    from rdkit import Chem
    from rdkit.Chem import rdForceFieldHelpers, AllChem
    from rdkit.Geometry import Point3D
    from rdkit.Chem import Draw
    import openmm.unit as unit
    import openmm.app as openmm_app
    from openff.toolkit.topology import Molecule
    from openff.units.openmm import to_openmm
    assert os.path.splitext(ligand_sdf_filename)[-1].lower() == ".sdf", \
        "Ligand SDF file must have a .sdf extension."
    suppl = Chem.SDMolSupplier(ligand_sdf_filename)
    mols = [x for x in suppl if x is not None]
    mol = mols[0]
    if draw_ligand:
        # Generate 2D coordinates for better visualization
        mol_2d = Chem.Mol(mol)
        AllChem.Compute2DCoords(mol_2d)
        img = Draw.MolToImage(mol_2d, size=(300, 300))
        img.save(os.path.join(param_directory,  "ligand.png"))
    ligand_pdb = parmed.load_file(str(ligand_pdb_filename), skip_bonds=True)
    ligand_positions_in_A_pdb = ligand_pdb.coordinates
    conf = mol.GetConformer()
    assert ligand_positions_in_A_pdb.shape[0] == mol.GetNumAtoms()
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, Point3D(ligand_positions_in_A_pdb[i,0],
                                        ligand_positions_in_A_pdb[i,1],
                                        ligand_positions_in_A_pdb[i,2]))
    
    rdForceFieldHelpers.MMFFOptimizeMolecule(mol)
    offmol = Molecule.from_rdkit(mol)
    ligand_positions = to_openmm(offmol.conformers[0])
    ligand_topology = offmol.to_topology().to_openmm()
    output_pdb_filename_ligand = os.path.join(param_directory, "ligand_with_connections.pdb")
    openmm_app.PDBFile.writeFile(ligand_topology, ligand_positions.value_in_unit(unit.nanometers), 
                                file=open(output_pdb_filename_ligand, 'w'))
    ligand_md_topology = mdtraj.Topology.from_openmm(ligand_topology)
    for atom in ligand_md_topology.atoms:
        if ligand_resname == "":
            atom.residue.name = "LIG"
        else:
            atom.residue.name = ligand_resname

    return (
        ligand_topology,
        ligand_md_topology,
        ligand_positions,
        offmol
    )
