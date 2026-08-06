"""
modules/parameterize_algorithms.py

Core force-field generation algorithms used by the parameterize pipeline:
SystemGenerator (GAFF/OpenFF), espaloma protein remapping, and solvation.
"""

import os
import typing

import numpy as np
import mdtraj

import seekrflow.modules.base as base
import seekrflow.modules.parameterize_structures as parameterize_structures
import seekrflow.modules.workflows.scale_settings as scale_settings_module

# Default hydrogen mass (amu) when MD scale settings leave hydrogen_mass unset.
DEFAULT_HYDROGEN_MASS = 1.008
DEFAULT_BAROSTAT_PERIOD = 25


def regenerate_espaloma_system(
        system_generator: "SystemGenerator",
        solvated_topology: "openmm.app.Topology",
        solvated_positions: "openmm.unit.Quantity",
        parameterize_directory: str,
        ) -> "openmm.System":
    """
    Rebuild a solvated system so espaloma can parameterize protein chains
    as OpenFF Molecules (one residue per protein chain).
    """
    import openmm.app as openmm_app
    from openff.toolkit.topology import Molecule
    output_pdb_basename = os.path.join(parameterize_directory, "complex")
    mdtop = mdtraj.Topology.from_openmm(solvated_topology)
    chain_indices = [chain.index for chain in solvated_topology.chains()]
    protein_chain_indices = [
        chain_index for chain_index in chain_indices
        if mdtop.select(f"protein and chainid == {chain_index}").any()]
    new_solvated_topology = openmm_app.Topology()
    new_solvated_topology.setPeriodicBoxVectors(
        solvated_topology.getPeriodicBoxVectors())
    new_atoms = {}
    chain_counter = 0
    for chain in solvated_topology.chains():
        new_chain = new_solvated_topology.addChain(chain.id)
        if chain.index in protein_chain_indices:
            resname = f"XX{chain_counter:01d}"
            resid = "1"
            chain_counter += 1
            new_residue = new_solvated_topology.addResidue(
                resname, new_chain, resid)
        for residue in chain.residues():
            if residue.chain.index not in protein_chain_indices:
                new_residue = new_solvated_topology.addResidue(
                    residue.name, new_chain, residue.id)
            for atom in residue.atoms():
                new_atom = new_solvated_topology.addAtom(
                    atom.name, atom.element, new_residue, atom.id)
                new_atoms[atom] = new_atom
    for bond in solvated_topology.bonds():
        if bond[0] in new_atoms and bond[1] in new_atoms:
            new_solvated_topology.addBond(
                new_atoms[bond[0]], new_atoms[bond[1]])
    complex_espaloma_filename = f"{output_pdb_basename}-espaloma-bad-resids.pdb"
    openmm_app.PDBFile.writeFile(
        new_solvated_topology, solvated_positions,
        file=open(complex_espaloma_filename, "w"))
    protein_espaloma_filenames = []
    for chain_index in protein_chain_indices:
        t = mdtraj.load_pdb(complex_espaloma_filename)
        indices = t.topology.select(f"chainid == {chain_index}")
        chain_espaloma_filename = \
            f"{output_pdb_basename}-espaloma-{chain_index}.pdb"
        t.atom_slice(indices).save_pdb(chain_espaloma_filename)
        protein_espaloma_filenames.append(chain_espaloma_filename)

    protein_molecules = [
        Molecule.from_file(protein_filename)
        for protein_filename in protein_espaloma_filenames]
    system_generator.template_generator.add_molecules(protein_molecules)
    return system_generator.create_system(new_solvated_topology)


def parameterize_and_check_complex(
        complex_topology: "openmm.app.Topology",
        complex_positions: np.ndarray,
        offmols: typing.Optional[
            typing.Union[
                "openff.toolkit.topology.Molecule",
                list["openff.toolkit.topology.Molecule"],
            ]],
        param_dir: str,
        parameterizer: parameterize_structures.Parameterizer,
        physical_attributes: base.Physical_attributes,
        md_settings: scale_settings_module.Molecular_dynamics_scale_settings,
        ) -> tuple[str, str]:
    """
    Parameterize and solvate the complex; write canonical XML and PDB outputs.
    """
    import openmm
    import openmm.unit as unit
    import openmm.app as openmm_app
    from openmmforcefields.generators import SystemGenerator

    modeller = openmm_app.Modeller(complex_topology, complex_positions)
    nonsolvated_topology = modeller.getTopology()
    nonsolvated_positions = modeller.getPositions()
    hydrogen_mass = md_settings.hydrogen_mass
    if hydrogen_mass is None:
        hydrogen_mass = DEFAULT_HYDROGEN_MASS
    forcefield_kwargs = {
        "removeCMMotion": True,
        "ewaldErrorTolerance": base.PME_TOL,
        "constraints": openmm_app.HBonds,
        "rigidWater": True,
        "hydrogenMass": hydrogen_mass * unit.amu,
    }
    periodic_forcefield_kwargs = {"nonbondedMethod": openmm_app.PME}
    if physical_attributes.pressure is not None:
        barostat = openmm.MonteCarloBarostat(
            physical_attributes.pressure * unit.atmosphere,
            physical_attributes.temperature * unit.kelvin,
            DEFAULT_BAROSTAT_PERIOD)
    else:
        barostat = None

    system_generator = SystemGenerator(
        forcefields=parameterizer.forcefields,
        forcefield_kwargs=forcefield_kwargs,
        periodic_forcefield_kwargs=periodic_forcefield_kwargs,
        barostat=barostat,
        small_molecule_forcefield=parameterizer.small_molecule_forcefield,
        molecules=offmols,
        cache=None)
    nonsolvated_system = system_generator.create_system(nonsolvated_topology)
    friction = md_settings.friction_coefficient
    if friction is None:
        friction = 1.0
    integrator = openmm.LangevinMiddleIntegrator(
        physical_attributes.temperature * unit.kelvin,
        friction / unit.picosecond,
        md_settings.timestep * unit.picoseconds)
    simulation = openmm_app.Simulation(
        nonsolvated_topology, nonsolvated_system, integrator)
    simulation.context.setPositions(nonsolvated_positions)
    state = simulation.context.getState(getPositions=True)
    nonsolvated_positions_minimized = state.getPositions()
    output_pdb_filename_nosolv = os.path.join(
        param_dir, parameterize_structures.COMPLEX_NO_SOLVENT_PDB)
    openmm_app.PDBFile.writeFile(
        nonsolvated_topology, nonsolvated_positions_minimized,
        file=open(output_pdb_filename_nosolv, "w"))

    solv_modeller = openmm_app.Modeller(
        nonsolvated_topology, nonsolvated_positions_minimized)
    solv_modeller.addSolvent(
        system_generator.forcefield,
        model=parameterizer.water_model,
        padding=parameterizer.solvent_padding * unit.nanometers,
        boxShape=parameterizer.box_shape,
        ionicStrength=physical_attributes.ionic_strength * unit.molar)
    solvated_topology = solv_modeller.getTopology()
    solvated_positions = solv_modeller.getPositions()
    solvated_system = system_generator.create_system(solvated_topology)

    serialized_system_xml = os.path.join(
        param_dir, parameterize_structures.COMPLEX_SYSTEM_XML)
    if parameterizer.is_espaloma():
        new_system = regenerate_espaloma_system(
            system_generator,
            solvated_topology,
            solvated_positions,
            param_dir)
    else:
        new_system = solvated_system

    with open(serialized_system_xml, "w") as wf:
        wf.write(openmm.XmlSerializer.serialize(new_system))

    # Keep a Simulation context for periodic imaging of positions.
    integrator = openmm.LangevinMiddleIntegrator(
        physical_attributes.temperature * unit.kelvin,
        friction / unit.picosecond,
        md_settings.timestep * unit.picoseconds)
    simulation = openmm_app.Simulation(
        solvated_topology, solvated_system, integrator)
    simulation.context.setPositions(solvated_positions)
    output_pdb_filename = os.path.join(
        param_dir, parameterize_structures.COMPLEX_SOLVENT_PDB)
    openmm_app.PDBFile.writeFile(
        solvated_topology, solvated_positions,
        file=open(output_pdb_filename, "w"))
    traj = mdtraj.load(output_pdb_filename)
    traj.topology = mdtraj.Topology.from_openmm(solvated_topology)
    traj.image_molecules(inplace=True)
    openmm_app.PDBFile.writeFile(
        solvated_topology, traj.xyz[0] * 10.0,
        file=open(output_pdb_filename, "w"))

    return serialized_system_xml, output_pdb_filename
