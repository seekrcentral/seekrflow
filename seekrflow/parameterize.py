"""
parameterize.py

Parameterize a molecular complex for input to a seekr calculation.

Primary entry is via ``flow.py`` (instruction ``parameterize`` / ``any``),
the same pattern as prepare/run. The ``main()`` CLI here is a deprecated
thin convenience wrapper only.
"""

import os
import typing
import argparse
from shutil import copyfile

import mdtraj
import numpy as np
import parmed
import openmm.unit as unit

import seekrflow.modules.base as base
import seekrflow.modules.structures as structures
import seekrflow.modules.parameterize_structures as parameterize_structures
import seekrflow.modules.parameterize_algorithms as parameterize_algorithms
import seekrflow.modules.parameters_topology_structures \
    as parameters_topology_structures
import seekrflow.modules.workflows.components as components_module
import seekrflow.modules.workflows.scale_settings as scale_settings_module
import seekrflow.modules.client.structures as client_structures
import seekrflow.modules.client.start as client_start

PARAM_SEEKRFLOW_GLOB = "seekrflow_pre_param_*.json"
PARAM_SEEKRFLOW_BASE = "seekrflow_pre_param_{}.json"


def _md_system_is_populated(
        md_settings: scale_settings_module.Molecular_dynamics_scale_settings,
        ) -> bool:
    """
    True if the MD scale already has an explicit system definition.
    """
    if md_settings.system is None:
        return False
    system = md_settings.system
    if system.solvated_pdb:
        return True
    pt = system.parameters_topology
    if pt is None:
        return False
    # Any concrete topology path counts as populated.
    for attr in ("system_filename", "prmtop_filename", "psf_filename",
                 "top_filename", "gro_filename"):
        if getattr(pt, attr, ""):
            return True
    if getattr(pt, "built_in_forcefield_filenames", None) \
            or getattr(pt, "custom_forcefield_filenames", None) \
            or getattr(pt, "param_filename_list", None):
        return True
    return False


def _assert_parameterize_system_contract(
        seekrflow: structures.Seekrflow,
        ) -> None:
    """
    Parameterizer and an already-populated MD system are mutually exclusive.
    """
    if seekrflow.parameterizer is None:
        return
    md_settings = seekrflow.workflow.get_md_settings()
    if _md_system_is_populated(md_settings):
        raise Exception(
            "Seekrflow.parameterizer is set, but "
            "workflow.scale_settings molecular_dynamics.system is already "
            "populated. Leave system as null/blank when parameterizing, or "
            "set parameterizer to null when providing an existing system.")


def fill_md_system_from_parameterize(
        seekrflow: structures.Seekrflow,
        ) -> None:
    """
    If MD ``system`` is blank, fill it from canonical parameterize outputs
    under ``work/parameterize/``.
    """
    md_settings = seekrflow.workflow.get_md_settings()
    if _md_system_is_populated(md_settings):
        return
    param_dir = structures.PARAMETERIZE
    xml_rel = os.path.join(
        param_dir, parameterize_structures.COMPLEX_SYSTEM_XML)
    pdb_rel = os.path.join(
        param_dir, parameterize_structures.COMPLEX_SOLVENT_PDB)
    xml_abs = os.path.join(seekrflow.work_directory, xml_rel)
    pdb_abs = os.path.join(seekrflow.work_directory, pdb_rel)
    if not (os.path.exists(xml_abs) and os.path.exists(pdb_abs)):
        raise FileNotFoundError(
            "MD scale system is blank and parameterized outputs were not "
            f"found at {xml_abs} and {pdb_abs}. Run the parameterize "
            "instruction first, or set scale_settings.system explicitly.")
    md_settings.system = scale_settings_module.System_for_md(
        parameters_topology=parameters_topology_structures.Openmm_system(
            system_filename=xml_rel),
        solvated_pdb=pdb_rel,
    )


def _component_pdb_filename(component_name: str) -> str:
    return f"{component_name}.pdb"


def _component_sdf_filename(component_name: str) -> str:
    return f"{component_name}.sdf"


def make_small_molecule_sdf_from_pdb(
        ligand_pdb_filename: str,
        ligand_sdf_filename: str,
        pdb_to_sdf_settings: parameterize_structures.PDB_to_SDF_settings,
        ) -> None:
    """
    Convert a small-molecule PDB to SDF using ``pdb_to_sdf_settings``.

    Prefer an explicit ``Small_molecule_component.sdf_file`` when possible.
    Callers must pass non-None settings; ``Parameterizer.pdb_to_sdf_settings``
    of ``None`` means conversion is disabled.
    """
    pdb_to_sdf_settings.run(ligand_pdb_filename, ligand_sdf_filename)


def split_components(
        seekrflow: structures.Seekrflow,
        param_directory: str,
        complex_pdb_filename: str,
        ) -> dict[str, list[int]]:
    """
    Split the complex PDB into one PDB per Workflow component. Returns
    resolved component name -> atom indices (in the complex).
    """
    traj = mdtraj.load(complex_pdb_filename)
    full_structure = parmed.load_file(complex_pdb_filename, skip_bonds=True)
    resolved: dict[str, list[int]] = {}
    members = seekrflow.workflow.components.members
    assert len(members) > 0, \
        "Workflow.components.members must be non-empty to parameterize."
    for member in members:
        indices = member.resolve(traj)
        assert len(indices) > 0, \
            f"Component '{member.name}' resolved to no atoms."
        resolved[member.name] = indices
        # ParmEd uses 1-based atom numbers in selection strings.
        selection = "@" + ",".join(str(i + 1) for i in indices)
        piece = full_structure[selection]
        out_path = os.path.join(
            param_directory, _component_pdb_filename(member.name))
        piece.save(str(out_path), overwrite=True)
    return resolved


def _apply_protein_fixer_pdb2pqr(
        parameterizer: parameterize_structures.Parameterizer,
        protein_pdb: str,
        output_basename: str,
        ) -> str:
    """
    Optionally run PDBFixer and PDB2PQR on a protein PDB; return latest PDB.
    """
    latest = protein_pdb
    if parameterizer.pdb_fixer_settings is not None:
        print("PDB chains before fixer:",
              parameterize_structures.describe_pdb_chains(latest))
        fixed = output_basename + "_fixed.pdb"
        parameterizer.pdb_fixer_settings.run(latest, fixed)
        latest = fixed
    if parameterizer.pdb2pqr_settings is not None:
        pqr_out = output_basename + "_pdb2pqr.pqr"
        pdb_out = output_basename + "_pdb2pqr.pdb"
        parameterizer.pdb2pqr_settings.run(latest, pqr_out, pdb_out)
        latest = pdb_out
    return latest

def _add_missing_hydrogens(
        input_pdb: str,
        output_pdb: str,
        pH: float = 7.0,
        ) -> str:
    """
    Add hydrogens to a crystallographic-water PDB. Does not remove chains.
    """
    from pdbfixer import PDBFixer
    import openmm.app as openmm_app
    fixer = PDBFixer(filename=input_pdb)
    fixer.addMissingHydrogens(pH)
    with open(output_pdb, "w") as handle:
        openmm_app.PDBFile.writeFile(
            fixer.topology, fixer.positions, handle)
    return output_pdb

def _append_onto_chain(
        base_top,
        base_positions,
        add_top,
        add_positions,
        chain_index: int = 0,
        ):
    """
    Add residues onto an existing chain so their chainid stays chain_index.
    """
    if base_top is None:
        return add_top, add_positions
    chain = base_top.chain(chain_index)
    atom_map = {}
    for residue in add_top.residues:
        new_residue = base_top.add_residue(
            residue.name, chain, resSeq=residue.resSeq)
        for atom in residue.atoms:
            new_atom = base_top.add_atom(atom.name, atom.element, new_residue)
            atom_map[atom] = new_atom
    for atom1, atom2 in add_top.bonds:
        if atom1 in atom_map and atom2 in atom_map:
            base_top.add_bond(atom_map[atom1], atom_map[atom2])
    n_base = base_positions.shape[0]
    new_positions = np.zeros(
        [n_base + add_positions.shape[0], 3]) * unit.nanometers
    new_positions[:n_base, :] = base_positions
    new_positions[n_base:, :] = add_positions
    return base_top, new_positions

def _join_mdtraj_topologies(
        base_md_top: mdtraj.Topology | None,
        base_positions,
        add_md_top: mdtraj.Topology,
        add_positions,
        ):
    """
    Join an MDTraj topology/positions onto an accumulating complex.
    """
    if base_md_top is None:
        return add_md_top, add_positions
    n_base = base_md_top.n_atoms
    n_add = add_md_top.n_atoms
    joined = base_md_top.join(add_md_top)
    new_positions = np.zeros([n_base + n_add, 3]) * unit.nanometers
    new_positions[:n_base, :] = base_positions
    new_positions[n_base:, :] = add_positions
    return joined, new_positions

def _stamp_member_chain_ids(topology, members) -> None:
    """
    After a MDtraj complex is formed from component topologies, the chains 
    have empty IDs. Assign them here.
    """
    chains = list(topology.chains())
    if len(chains) != len(members):
        print(
            f"Not stamping chain ids: {len(chains)} chains for "
            f"{len(members)} members.")
        return
    for chain, member in zip(chains, members):
        chain_str = getattr(member.selector, "chain", None)
        if chain_str:
            chain.id = str(chain_str).strip()

def create_complex(
        seekrflow: structures.Seekrflow,
        param_directory: str,
        ) -> typing.Tuple[str, str]:
    """
    Build the solvated, parameterized complex from Workflow components.

    Non-small-molecule components (proteins, etc.) are fixed/protonated first
    and joined; small molecules are then converted via SDF and appended.
    """
    parameterizer = seekrflow.parameterizer
    assert parameterizer is not None
    members = seekrflow.workflow.components.members
    complex_md_topology = None
    complex_positions = None
    offmols = []
    
    solvent = [
        m for m in members
        if isinstance(m, components_module.Solvent_component)]
    non_small = [
        m for m in members
        if not isinstance(m,( 
            components_module.Small_molecule_component,
            components_module.Solvent_component))]
    small = [
        m for m in members
        if isinstance(m, components_module.Small_molecule_component)]

    for member in non_small:
        component_pdb = os.path.join(
            param_directory, _component_pdb_filename(member.name))
        assert os.path.exists(component_pdb), \
            f"Component PDB {component_pdb} does not exist."
        basename = os.path.join(param_directory, member.name)
        latest_pdb = _apply_protein_fixer_pdb2pqr(
            parameterizer, component_pdb, basename)
        _prot_top, prot_md_top, prot_pos = \
            base.make_protein_openmm_and_mdtraj_top(latest_pdb)
        complex_md_topology, complex_positions = _join_mdtraj_topologies(
            complex_md_topology, complex_positions, prot_md_top, prot_pos)
    
    for member in solvent:
        component_pdb = os.path.join(
            param_directory, _component_pdb_filename(member.name))
        assert os.path.exists(component_pdb), \
            f"Component PDB {component_pdb} does not exist."
        pH = 7.0
        if (parameterizer.pdb_fixer_settings is not None
                and parameterizer.pdb_fixer_settings.add_missing_hydrogens_pH
                is not None):
            pH = parameterizer.pdb_fixer_settings.add_missing_hydrogens_pH
        hydrogenated_pdb = os.path.join(
            param_directory, f"{member.name}_h.pdb")
        _add_missing_hydrogens(component_pdb, hydrogenated_pdb, pH=pH)
        _wat_top, wat_md_top, wat_pos = \
            base.make_protein_openmm_and_mdtraj_top(hydrogenated_pdb)
        complex_md_topology, complex_positions = _append_onto_chain(
            complex_md_topology, complex_positions, wat_md_top, wat_pos,
            chain_index=0)

    for member in small:
        component_pdb = os.path.join(
            param_directory, _component_pdb_filename(member.name))
        assert os.path.exists(component_pdb), \
            f"Component PDB {component_pdb} does not exist."
        sdf_dest = os.path.join(
            param_directory, _component_sdf_filename(member.name))
        if member.sdf_file: 
            assert os.path.exists(member.sdf_file), \
                f"SDF file {member.sdf_file} does not exist."
            if os.path.abspath(member.sdf_file) != os.path.abspath(sdf_dest):
                copyfile(member.sdf_file, sdf_dest)
        else:
            if parameterizer.pdb_to_sdf_settings is None:
                raise ValueError(
                    f"Small molecule '{member.name}' has no usable sdf_file "
                    f"and parameterizer.pdb_to_sdf_settings is None. Provide "
                    f"Small_molecule_component.sdf_file, or set "
                    f"pdb_to_sdf_settings (e.g. engine='rdkit') to convert "
                    f"from the split PDB."
                )
            make_small_molecule_sdf_from_pdb(
                component_pdb, sdf_dest,
                pdb_to_sdf_settings=parameterizer.pdb_to_sdf_settings)
        member.sdf_file = sdf_dest
        chain = getattr(member.selector, "chain", None)
        _lig_top, lig_md_top, lig_pos, offmol = \
            base.make_ligand_openmm_and_mdtraj_top(
                sdf_dest,
                component_pdb,
                param_directory,
                draw_ligand=True,
                ligand_resname=member.parameterize_resname(),
                chain=chain,
            )
        offmols.append(offmol)
        complex_md_topology, complex_positions = _join_mdtraj_topologies(
            complex_md_topology, complex_positions, lig_md_top, lig_pos)

    assert complex_md_topology is not None, \
        "No components produced a topology for parameterization."
    complex_topology = complex_md_topology.to_openmm()
    _stamp_member_chain_ids(
        complex_topology, list(non_small) + list(small))
    md_settings = seekrflow.workflow.get_md_settings()
    molecules_arg: typing.Any = None
    if len(offmols) == 1:
        molecules_arg = offmols[0]
    elif len(offmols) > 1:
        molecules_arg = offmols
    return parameterize_algorithms.parameterize_and_check_complex(
        complex_topology,
        complex_positions,
        molecules_arg,
        param_directory,
        parameterizer,
        seekrflow.physical_attributes,
        md_settings,
    )

def parameterize(
        batch: client_structures.Batch,
        ) -> typing.Tuple[typing.Optional[str], typing.Optional[str]]:
    """
    Parameterize a batch of seekrflow workflows. Writes canonical outputs 
    under work/parameterize/.
    """
    curdir = os.getcwd()
    if batch.batch_directory is not None:
        os.chdir(batch.batch_directory)
    for seekrflow, batch_system in zip(batch.create_seekrflow_objects(), batch.systems):
        _assert_parameterize_system_contract(seekrflow)
        if seekrflow.parameterizer is None:
            print("seekrflow.parameterizer is None; nothing to parameterize.")
            return None, None

        parameterizer = seekrflow.parameterizer
        assert parameterizer.complex_pdb_filename, \
            "parameterizer.complex_pdb_filename must be set."
        assert os.path.exists(parameterizer.complex_pdb_filename), \
            f"Complex PDB {parameterizer.complex_pdb_filename} does not exist."

        seekrflow.make_work_directory()
        work_param_dir = seekrflow.get_parameterize_directory()
        param_dir = structures.PARAMETERIZE
        
        pdb_src = parameterizer.complex_pdb_filename
        work_copy_pdb = os.path.join(
            seekrflow.work_directory, os.path.basename(pdb_src))
        copyfile(pdb_src, work_copy_pdb)
        work_copy = mdtraj.load(work_copy_pdb)
        non_h_indices = work_copy.topology.select("not element H")
        traj_noh = work_copy.atom_slice(non_h_indices)
        base_name, ext = os.path.splitext(work_copy_pdb)
        work_copy_pdb_noh = f"{base_name}_noh{ext}"
        print(f"Saving noh file at: {work_copy_pdb_noh}")
        traj_noh.save(work_copy_pdb_noh)
        parameterizer.complex_pdb_filename = os.path.basename(work_copy_pdb)

        # Copy user-provided SDFs into the parameterize directory up front.
        for member in seekrflow.workflow.components.members:
            if isinstance(member, components_module.Small_molecule_component) \
                    and member.sdf_file:
                assert os.path.exists(member.sdf_file), \
                    f"SDF file {member.sdf_file} does not exist."
                dest = os.path.join(
                    work_param_dir, _component_sdf_filename(member.name))
                copyfile(member.sdf_file, dest)
                member.sdf_file = os.path.join(
                    param_dir, _component_sdf_filename(member.name))

        os.chdir(seekrflow.work_directory)
        try:
            split_components(
                seekrflow, param_dir, parameterizer.complex_pdb_filename)
            serialized_system_xml, output_pdb_filename = create_complex(
                seekrflow, param_dir)
        finally:
            os.chdir(batch.batch_directory)
            
    os.chdir(curdir)
    return serialized_system_xml, output_pdb_filename


def main() -> None:
    argparser = argparse.ArgumentParser(
        description="Parameterize a molecular complex for input to a seekr "
        "calculation.")
    argparser.add_argument(
        "input_json", metavar="INPUT_JSON", type=str,
        help="Input JSON file containing the job specification.")
    
    args = vars(argparser.parse_args())
    input_json = args["input_json"]
    if not os.path.exists(input_json):
        raise FileNotFoundError(f"Input JSON file {input_json} does not exist.")
    #seekrflow = structures.load_seekrflow(input_json)
    #if seekrflow.parameterizer is None:
    #    seekrflow.parameterizer = parameterize_structures.Parameterizer()
    #seekrflow.make_work_directory(seekrflow.work_directory)
    if client_start.input_is_batch_file(input_json):
        batch = client_structures.Batch.load_batch_file(input_json)
    else:
        seekrflow = structures.load_seekrflow(input_json)
        work_dir = os.path.abspath(
            os.path.expanduser(seekrflow.work_directory))
        seekrflow.make_work_directory(work_dir)
        batch = client_start.make_batch_from_single_seekrflow(seekrflow)
    parameterize(batch)
    #structures.save_new_seekrflow(
    #    seekrflow, PARAM_SEEKRFLOW_GLOB, PARAM_SEEKRFLOW_BASE,
    #    directory=seekrflow.work_directory)

if __name__ == "__main__":
    main()
