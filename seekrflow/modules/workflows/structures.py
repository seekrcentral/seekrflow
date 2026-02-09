"""
modules/workflows/structures.py

Contain data structure classes used for workflows that can be
imported for any of the workflows.
"""
import os
import typing

from attrs import define, field, validators
import numpy as np
import mdtraj

import seekrflow.modules.base as base
import seekrflow.modules.parameters_topology_structures as parameters_topology_structures
import seekrflow.modules.parameterize_structures as parameterize_structures

MAX_MINIMIZATION_ITERATIONS = 0

@define
class MD_settings:
    """
    Settings for the MD calculation.
    """
    type: typing.Literal["MD"] = "MD"
    engine: str = field(
        default="openmm",
        validator=validators.in_(["openmm"]),
        )
    integrator: str = field(
        default="langevin",
        validator=validators.in_(["langevin", "langevinmiddle"]),
        )
    nonbonded_cutoff: float | None = field(
        default=0.9,
        validator=validators.instance_of(float),
        )
    friction: float = field(
        default=1.0,
        validator=validators.instance_of(float),
        )
    barostat_period: int | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(int)),
        )
    stepsize: float = field(
        default=0.002,
        validator=validators.instance_of(float),
        )
    
@define
class BD_settings:
    """
    Settings for the BD calculation.
    """
    type: typing.Literal["BD"] = "BD"
    engine: str = field(
        default="browndye",
        validator=validators.in_(["browndye", "sda"]),
        )
    binary_directory: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    auxiliary_directory: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    hydropro_directory: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    num_trajectories: int = field(
        default=10000,
        validator=validators.instance_of(int),
        )
    num_threads: int = field(
        default=1,
        validator=validators.instance_of(int),
        )
    
@define
class Workflow_base:
    """
    Base class for workflows.
    """
    type: typing.Literal["workflow_base"] = "workflow_base"

@define
class Solvated_system_for_md:
    """
    Information about the solvated system for MD simulations.
    """
    parameters_topology: parameters_topology_structures.Parameters_topology | None = None
    solvated_pdb: str = field(
        default="",
        validator=validators.instance_of(str),
        )

def regenerate_espaloma_system(
        system_generator: "SystemGenerator",
        solvated_topology: "openmm.app.Topology",
        solvated_positions: "openmm.unit.Quantity",
        parameterize_directory: str,
        ) -> "openmm.System":
    """
    Add espaloma parameterization to the solvated system.
    """
    import openmm.app as openmm_app
    from openff.toolkit.topology import Molecule
    output_pdb_basename = os.path.join(parameterize_directory, "complex")
    # Convert solvated topology to MDTraj format and identify protein chains
    mdtop = mdtraj.Topology.from_openmm(solvated_topology)
    chain_indices = [chain.index for chain in solvated_topology.chains()]
    protein_chain_indices \
        = [chain_index for chain_index in chain_indices if mdtop.select(
            f"protein and chainid == {chain_index}").any()]
    # Create new topology and copy chains
    new_solvated_topology = openmm_app.Topology()
    new_solvated_topology.setPeriodicBoxVectors(solvated_topology.getPeriodicBoxVectors())
    new_atoms = {}
    chain_counter = 0
    for chain in solvated_topology.chains():
        new_chain = new_solvated_topology.addChain(chain.id)
        if chain.index in protein_chain_indices:
            resname = f'XX{chain_counter:01d}'  # Assign unique residue name for each protein chain
            resid = "1"
            chain_counter += 1
            new_residue = new_solvated_topology.addResidue(resname, new_chain, resid)
        for residue in chain.residues():
            if residue.chain.index not in protein_chain_indices:
                new_residue = new_solvated_topology.addResidue(residue.name, new_chain, residue.id)
            for atom in residue.atoms():
                new_atom = new_solvated_topology.addAtom(atom.name, atom.element, new_residue, atom.id)
                new_atoms[atom] = new_atom
    for bond in solvated_topology.bonds():
        if bond[0] in new_atoms and bond[1] in new_atoms:
            new_solvated_topology.addBond(new_atoms[bond[0]], new_atoms[bond[1]])
    # Save the complex with ESPALOMA parameterization to PDB file
    complex_espaloma_filename = f"{output_pdb_basename}-espaloma-bad-resids.pdb"
    openmm_app.PDBFile.writeFile(new_solvated_topology, solvated_positions, 
                                 file=open(complex_espaloma_filename, 'w'))
    # Split protein chains into separate PDB files
    protein_espaloma_filenames = []
    for chain_index in protein_chain_indices:
        t = mdtraj.load_pdb(complex_espaloma_filename)
        indices = t.topology.select(f"chainid == {chain_index}")
        chain_espaloma_filename = f"{output_pdb_basename}-espaloma-{chain_index}.pdb"
        t.atom_slice(indices).save_pdb(chain_espaloma_filename)
        protein_espaloma_filenames.append(chain_espaloma_filename)

    # Load protein molecules and add to system generator template
    protein_molecules = [Molecule.from_file(protein_filename) for protein_filename in protein_espaloma_filenames]
    system_generator.template_generator.add_molecules(protein_molecules)
    # Create new solvated system
    new_solvated_system = system_generator.create_system(new_solvated_topology)
    return new_solvated_system

def parameterize_and_check_complex(
        complex_topology: "openmm.app.Topology",
        complex_positions: np.ndarray,
        offmol: typing.Optional["openff.toolkit.topology.Molecule"],
        param_dir: str,
        parameterizer: parameterize_structures.Parameterizer,
        physical_attributes: base.Physical_attributes,
        md_settings: MD_settings,
        ) -> tuple[str, str]:
    """
    Parameterize the complex structure for the given workflow, and 
    run minimizations and short sims to check the structure is stable.
    """
    import openmm
    import openmm.unit as unit
    import openmm.app as openmm_app
    from openmmforcefields.generators import SystemGenerator
    modeller = openmm_app.Modeller(complex_topology, complex_positions)
    nonsolvated_topology = modeller.getTopology()
    nonsolvated_positions = modeller.getPositions()
    forcefield_kwargs = {'removeCMMotion': True, 
                        'ewaldErrorTolerance': base.PME_TOL, 
                        'constraints': openmm_app.HBonds, 
                        'rigidWater': True, 
                        'hydrogenMass': physical_attributes.hydrogen_mass * unit.amu}
    periodic_forcefield_kwargs = {'nonbondedMethod': openmm_app.PME}
    if physical_attributes.pressure is not None:
        barostat = openmm.MonteCarloBarostat(
            physical_attributes.pressure * unit.atmosphere, 
            physical_attributes.temperature * unit.kelvin, 
            md_settings.barostat_period)
    else:
        barostat = None
    
    system_generator = SystemGenerator(
        #forcefields=FF_FILES, 
        forcefields=parameterizer.auxiliary_forcefield_files,
        forcefield_kwargs=forcefield_kwargs, 
        periodic_forcefield_kwargs=periodic_forcefield_kwargs, 
        barostat=barostat, 
        small_molecule_forcefield=parameterizer.forcefield, 
        molecules=offmol, 
        cache=None)
    nonsolvated_system = system_generator.create_system(nonsolvated_topology)
    integrator = openmm.LangevinMiddleIntegrator(
        physical_attributes.temperature * unit.kelvin, 
        md_settings.friction / unit.picosecond, 
        md_settings.stepsize * unit.picoseconds)
    simulation = openmm_app.Simulation(nonsolvated_topology, nonsolvated_system, 
                                    integrator)
    simulation.context.setPositions(nonsolvated_positions)
    #simulation.minimizeEnergy(maxIterations=MAX_MINIMIZATION_ITERATIONS)
    state = simulation.context.getState(getPositions=True, )
    #                                    enforcePeriodicBox=True) # No PBC yet
    nonsolvated_positions_minimized = state.getPositions()
    output_pdb_filename_nosolv = os.path.join(param_dir, "complex_no_solvent.pdb")
    openmm_app.PDBFile.writeFile(nonsolvated_topology, nonsolvated_positions_minimized, 
                                file=open(output_pdb_filename_nosolv, 'w'))
    
    # TODO: set up possibility for triclinic water box
    solv_modeller = openmm_app.Modeller(nonsolvated_topology, nonsolvated_positions_minimized)
    solv_modeller.addSolvent(
        system_generator.forcefield, 
        model=parameterizer.water_model, 
        padding=parameterizer.solvent_padding * unit.nanometers, 
        boxShape=parameterizer.box_shape,
        ionicStrength=physical_attributes.ionic_strength * unit.molar)
    solvated_topology = solv_modeller.getTopology()
    solvated_positions = solv_modeller.getPositions()
    solvated_system = system_generator.create_system(solvated_topology)

    # NOTE: the commented-out code below, I mass repartition the hydrogens of the
    #  solvent to be heavier. It turns out that this is not standard practice - 
    #  normally, the masses of hydrogen are left alone in solvent molecules with
    #  a rigid structure.
    """
    # Modify hydrogens as needed
    hydrogens_needing_mass_change = {atom.index: atom for atom in solvated_topology.atoms() 
                                    if atom.element.symbol == "H"}
    atom_indices_bonded_to_each_hydrogen = {}
    #for h_atom in hydrogens_needing_mass_change:
    #    bonded_atoms = [bonded_atom for bonded_atom in h_atom.bonds()]
    #    assert len(bonded_atoms) == 1, \
    #        "Hydrogen atom is bonded to multiple atoms, unexpected."
    #    atoms_bonded_to_each_hydrogen.append(bonded_atoms[0].index)
    for bond in solvated_topology.bonds():
        if bond[0].index in hydrogens_needing_mass_change:
            atom_indices_bonded_to_each_hydrogen[bond[0].index] = bond[1].index
        elif bond[1].index in hydrogens_needing_mass_change:
            atom_indices_bonded_to_each_hydrogen[bond[1].index] = bond[0].index

    for h_index in hydrogens_needing_mass_change.keys():
        new_H_mass = physical_attributes.hydrogen_mass * unit.amu
        old_H_mass = solvated_system.getParticleMass(h_index)
        mass_diff = new_H_mass - old_H_mass
        solvated_system.setParticleMass(h_index, new_H_mass)
        bonded_atom_index = atom_indices_bonded_to_each_hydrogen[h_index]
        old_bonded_atom_mass = solvated_system.getParticleMass(bonded_atom_index)
        new_bonded_atom_mass = old_bonded_atom_mass - mass_diff
        solvated_system.setParticleMass(bonded_atom_index, new_bonded_atom_mass)
    """
    #output_pdb_filename = f"{output_pdb_basename}.pdb"
    #openmm_app.PDBFile.writeFile(solvated_topology, solvated_positions, file=open(output_pdb_filename, 'w'))

    serialized_xml = os.path.join(param_dir, "complex.xml")
    if parameterizer.forcefield.startswith("espaloma"):
        new_system = regenerate_espaloma_system(
            system_generator,
            solvated_topology,
            solvated_positions,
            param_dir)
    else:
        new_system = solvated_system
    
    with open(serialized_xml, "w") as wf:
        xml = openmm.XmlSerializer.serialize(new_system)  # Serialize system
        wf.write(xml)
    
    integrator = openmm.LangevinMiddleIntegrator(
        physical_attributes.temperature * unit.kelvin, 
        md_settings.friction / unit.picosecond, 
        md_settings.stepsize * unit.picoseconds)
    simulation = openmm_app.Simulation(solvated_topology, solvated_system, integrator)
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    simulation.context.setPositions(solvated_positions)
    output_pdb_filename_solv = os.path.join(param_dir, "complex_solvent.pdb")
    openmm_app.PDBFile.writeFile(solvated_topology, solvated_positions, 
                                file=open(output_pdb_filename_solv, 'w'))
    simulation.minimizeEnergy(maxIterations=MAX_MINIMIZATION_ITERATIONS)
    simulation.step(1000)
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    solvated_positions_equilibrated = state.getPositions()
    output_pdb_filename = os.path.join(param_dir, "complex-equil.pdb")
    openmm_app.PDBFile.writeFile(solvated_topology, solvated_positions_equilibrated, 
                                 file=open(output_pdb_filename, 'w'))
    traj = mdtraj.load(output_pdb_filename)
    traj.topology = mdtraj.Topology.from_openmm(solvated_topology)
    traj.image_molecules(inplace=True)
    openmm_app.PDBFile.writeFile(solvated_topology, traj.xyz[0]*10.0, 
                                 file=open(output_pdb_filename, 'w'))
    
    return serialized_xml, output_pdb_filename