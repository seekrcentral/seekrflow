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
        parametrize_directory: str,
        ) -> "openmm.System":
    """
    Add espaloma parameterization to the solvated system.
    """
    import openmm.app as openmm_app
    from openff.toolkit.topology import Molecule
    output_pdb_basename = os.path.join(parametrize_directory, "complex")
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

def parametrize_and_check_complex(
        complex_topology: "openmm.app.Topology",
        complex_positions: np.ndarray,
        offmol: "openff.toolkit.topology.Molecule",
        param_dir: str,
        parameterizer: parameterize_structures.Parameterizer,
        physical_attributes: base.Physical_attributes,
        md_settings: MD_settings,
        ) -> tuple[str, str]:
    """
    Parametrize the complex structure for the given workflow, and 
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
    
    # TODO: these need to be set in the seekrflow object
    #FF_FILES = [
    #    "amber/ff14SB.xml",
    #    "amber/tip3p_standard.xml",
    #    "amber/tip3p_HFE_multivalent.xml"
    #]
    
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
    simulation.minimizeEnergy(maxIterations=MAX_MINIMIZATION_ITERATIONS)
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
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
        ionicStrength=physical_attributes.ionic_strength * unit.molar)
    solvated_topology = solv_modeller.getTopology()
    solvated_positions = solv_modeller.getPositions()
    solvated_system = system_generator.create_system(solvated_topology)
    
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
    simulation.context.setPositions(solvated_positions)
    simulation.minimizeEnergy(maxIterations=MAX_MINIMIZATION_ITERATIONS)
    simulation.step(1000)
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    solvated_positions_equilibrated = state.getPositions()
    output_pdb_filename = os.path.join(param_dir, "complex-equil.pdb")
    openmm_app.PDBFile.writeFile(solvated_topology, solvated_positions_equilibrated, file=open(output_pdb_filename, 'w'))
    
    return serialized_xml, output_pdb_filename