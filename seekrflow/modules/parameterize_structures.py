"""
modules/parameterize_structures.py

Data structures and settings for seekrflow parameterization.
"""

import os
import typing

from attrs import define, field, validators, Factory
import openmm.app as openmm_app

# Canonical outputs written under work/parameterize/
COMPLEX_SYSTEM_XML = "complex_system.xml"
COMPLEX_SOLVENT_PDB = "complex_solvent.pdb"
COMPLEX_NO_SOLVENT_PDB = "complex_no_solvent.pdb"


def describe_pdb_chains(pdb_filename: str) -> list[dict]:
    """
    Inspect chains in a PDB for PDBFixer keep/remove decisions.

    Returns a list of dicts with keys: index, id, n_residues, n_atoms,
    residue_names (unique residue names in the chain). Does not modify
    the structure; callers decide what to keep.
    """
    traj = None
    try:
        import mdtraj
        traj = mdtraj.load(pdb_filename)
    except Exception:
        pdb = openmm_app.PDBFile(pdb_filename)
        chains_info = []
        for i, chain in enumerate(pdb.topology.chains()):
            residues = list(chain.residues())
            atoms = [a for r in residues for a in r.atoms()]
            resnames = sorted({r.name for r in residues})
            chains_info.append({
                "index": i,
                "id": chain.id,
                "n_residues": len(residues),
                "n_atoms": len(atoms),
                "residue_names": resnames,
            })
        return chains_info

    chains_info = []
    for chain in traj.topology.chains:
        indices = traj.topology.select(f"chainid {chain.index}")
        resnames = sorted({
            traj.topology.atom(int(i)).residue.name for i in indices
        }) if len(indices) else []
        n_residues = len({
            traj.topology.atom(int(i)).residue.index for i in indices
        }) if len(indices) else 0
        chains_info.append({
            "index": int(chain.index),
            "id": chain.chain_id if hasattr(chain, "chain_id") else str(chain.index),
            "n_residues": n_residues,
            "n_atoms": int(len(indices)),
            "residue_names": resnames,
        })
    return chains_info


@define
class PDBFixer_settings:
    """
    Settings for PDBFixer.

    Chain handling (checked in order):
      * ``keep_chain_indices`` — if set, keep only these zero-based chain
        indices (preferred for explicit control of symmetry mates / cofactors).
      * ``remove_extra_chains`` — if True and ``keep_chain_indices`` is None,
        keep only chain 0 (legacy default).

    Use :func:`describe_pdb_chains` to inspect a structure before choosing
    which chains to keep. Do not silently drop ambiguous cofactor chains.
    """
    remove_extra_chains: bool = field(
        default=True,
        validator=validators.instance_of(bool),
        )
    keep_chain_indices: list[int] | None = field(
        default=None,
        validator=validators.optional(
            validators.deep_iterable(
                member_validator=validators.instance_of(int),
                iterable_validator=validators.instance_of(list),
            )),
        )
    find_missing_residues: bool = field(
        default=True,
        validator=validators.instance_of(bool),
        )
    find_and_replace_nonstandard_residues: bool = field(
        default=True,
        validator=validators.instance_of(bool),
        )
    remove_heterogens: bool = field(
        default=False,
        validator=validators.instance_of(bool),
        )
    find_and_add_missing_atoms: bool = field(
        default=True,
        validator=validators.instance_of(bool),
        )
    add_missing_hydrogens_pH: float | None = None

    def run(
            self,
            input_pdb_filename: str,
            output_pdb_filename: str
            ) -> None:
        """
        Run PDBFixer on the given PDB file and write the fixed PDB.
        """
        from pdbfixer import PDBFixer
        fixer = PDBFixer(filename=input_pdb_filename)
        num_chains = len(list(fixer.topology.chains()))
        if self.keep_chain_indices is not None:
            assert len(self.keep_chain_indices) > 0, \
                "keep_chain_indices must be non-empty when set."
            for idx in self.keep_chain_indices:
                assert 0 <= idx < num_chains, \
                    f"keep_chain_indices entry {idx} is out of range "\
                    f"(structure has {num_chains} chains)."
            remove = [i for i in range(num_chains)
                      if i not in self.keep_chain_indices]
            if remove:
                fixer.removeChains(remove)
        elif self.remove_extra_chains and num_chains > 1:
            fixer.removeChains(range(1, num_chains))
        if self.find_missing_residues:
            fixer.findMissingResidues()
        if self.find_and_replace_nonstandard_residues:
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
        if self.remove_heterogens:
            fixer.removeHeterogens(keepWater=True)
        if self.find_and_add_missing_atoms:
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
        if self.add_missing_hydrogens_pH is not None:
            fixer.addMissingHydrogens(pH=self.add_missing_hydrogens_pH)
        openmm_app.PDBFile.writeFile(
            fixer.topology, fixer.positions, open(output_pdb_filename, "w"))
        return


@define
class PDB2PQR_settings:
    """
    Settings for PDB2PQR (receptor protonation for MD, not BD PQR output).
    """
    forcefield: str = field(
        default="PARSE",
        validator=validators.instance_of(str),
        )
    forcefield_output_format: str = field(
        default="AMBER",
        validator=validators.instance_of(str),
        )
    pH: float | None = field(
        default=7.0,
        validator=validators.optional(validators.instance_of(float))
        )

    def run(
            self,
            input_pdb_filename: str,
            output_pqr_filename: str,
            output_pdb_filename: str | None = None
            ) -> None:
        """
        Run PDB2PQR and write the resulting PQR (and optional PDB).
        """
        if output_pdb_filename is not None:
            output_pdb_string = f"--pdb-output {output_pdb_filename} "
        else:
            output_pdb_string = ""
        cmd = f"pdb2pqr --ff {self.forcefield} --ffout {self.forcefield_output_format} "\
        +f"{output_pdb_string} --with-ph {self.pH} --log-level CRITICAL --drop-water "\
        +f"--nodebump --noopt {input_pdb_filename} {output_pqr_filename}"
        print("running command:", cmd)
        os.system(cmd)
        assert os.path.exists(output_pqr_filename), \
            f"PDB2PQR output PQR file {output_pqr_filename} was not written. "\
            "A problem must have occurred"
        return


@define
class PDB_to_SDF_settings:
    """
    How to convert a split small-molecule PDB to SDF when no SDF was provided.

    PDB lacks reliable bond-order information, so conversion is best-effort.
    Prefer supplying ``Small_molecule_component.sdf_file`` when possible.
    When every small molecule already has an SDF, set
    ``Parameterizer.pdb_to_sdf_settings`` to ``None`` to disable conversion.

    Engines:
      * ``rdkit`` (default) — open-source; uses PDB connectivity and sanitizes.
      * ``openeye`` — optional; often more robust for tricky ligands, but we
        intend to phase it out in favor of RDKit / explicit SDFs.
    """
    engine: typing.Literal["rdkit", "openeye"] = field(
        default="rdkit",
        validator=validators.in_(["rdkit", "openeye"]),
        )

    def run(
            self,
            input_pdb_filename: str,
            output_sdf_filename: str,
            ) -> None:
        """
        Convert ``input_pdb_filename`` to ``output_sdf_filename``.
        """
        if self.engine == "rdkit":
            self._run_rdkit(input_pdb_filename, output_sdf_filename)
        elif self.engine == "openeye":
            self._run_openeye(input_pdb_filename, output_sdf_filename)
        else:
            raise ValueError(f"Unknown PDB→SDF engine: {self.engine}")

    @staticmethod
    def _run_rdkit(
            input_pdb_filename: str,
            output_sdf_filename: str,
            ) -> None:
        from rdkit import Chem
        # sanitize=False first so we can attempt bond-order perception ourselves
        # when PDB CONECT / element data is incomplete.
        mol = Chem.MolFromPDBFile(
            input_pdb_filename, removeHs=False, sanitize=False)
        assert mol is not None, \
            f"RDKit could not parse ligand PDB {input_pdb_filename}. "\
            "Provide Small_molecule_component.sdf_file explicitly, or try "\
            "pdb_to_sdf_settings.engine='openeye'."
        try:
            Chem.SanitizeMol(mol)
        except Exception as exc:
            raise RuntimeError(
                f"RDKit could not sanitize molecule from "
                f"{input_pdb_filename} ({exc}). PDB files often lack bond "
                "orders; supply an SDF via Small_molecule_component.sdf_file."
            ) from exc
        writer = Chem.SDWriter(output_sdf_filename)
        writer.write(mol)
        writer.close()
        assert os.path.exists(output_sdf_filename), \
            f"RDKit failed to write SDF {output_sdf_filename}."

    @staticmethod
    def _run_openeye(
            input_pdb_filename: str,
            output_sdf_filename: str,
            ) -> None:
        try:
            from openeye import oechem
        except ImportError as exc:
            raise ImportError(
                "pdb_to_sdf_settings.engine='openeye' requires the OpenEye "
                "toolkits. Install them, or use engine='rdkit', or provide "
                "Small_molecule_component.sdf_file."
            ) from exc
        ifs = oechem.oemolistream()
        ofs = oechem.oemolostream()
        if not ifs.open(str(input_pdb_filename)):
            oechem.OEThrow.Fatal(
                f"Unable to open the input PDB file: {input_pdb_filename}")
        if not ofs.open(output_sdf_filename):
            oechem.OEThrow.Fatal(
                f"Unable to create the output SDF file: {output_sdf_filename}")
        for mol in ifs.GetOEGraphMols():
            oechem.OEWriteMolecule(ofs, mol)
        ifs.close()
        ofs.close()
        assert os.path.exists(output_sdf_filename), \
            f"OpenEye failed to write SDF {output_sdf_filename}."


@define
class Parameterizer:
    """
    Global settings and inputs for parameterizing a molecular complex.

    Per-component chemistry inputs (e.g. SDF files) live on Workflow
    Components. This object holds the starting complex PDB, force-field /
    solvent settings, optional PDBFixer/PDB2PQR settings, and is the single
    Seekrflow attribute that enables parameterization when not None.

    When every small molecule already has ``sdf_file``, set
    ``pdb_to_sdf_settings`` to ``None`` to disable PDB→SDF conversion.
    Otherwise omit it (defaults to RDKit) or set ``engine`` explicitly.

    Force-field field names mirror ``openmmforcefields.SystemGenerator``:
      * ``forcefields`` — biopolymer / water / ion XML list (SystemGenerator
        ``forcefields=``), e.g. ``amber/ff14SB.xml``, ``amber/tip3p_standard.xml``.
      * ``small_molecule_forcefield`` — GAFF / OpenFF / espaloma name or path
        (SystemGenerator ``small_molecule_forcefield=``), e.g. ``gaff-2.11``,
        ``openff-2.0.0``, or ``espaloma-0.3.2.pt``.
    """
    complex_pdb_filename: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    # Maps to SystemGenerator(forcefields=...)
    forcefields: typing.List[str] = field(
        default=Factory(list),
        validator=validators.instance_of(list),
        )
    # Maps to SystemGenerator(small_molecule_forcefield=...)
    small_molecule_forcefield: str = field(
        default="gaff-2.11",
        validator=validators.instance_of(str),
        )
    water_model: str = field(
        default="tip3p",
        validator=validators.instance_of(str),
        )
    pdb_fixer_settings: PDBFixer_settings | None = field(
        default=Factory(PDBFixer_settings),
        validator=validators.optional(validators.instance_of(PDBFixer_settings)),
        )
    pdb2pqr_settings: PDB2PQR_settings | None = field(
        default=Factory(PDB2PQR_settings),
        validator=validators.optional(validators.instance_of(PDB2PQR_settings)),
        )
    # None disables PDB→SDF conversion; require Small_molecule_component.sdf_file.
    pdb_to_sdf_settings: PDB_to_SDF_settings | None = field(
        default=Factory(PDB_to_SDF_settings),
        validator=validators.optional(
            validators.instance_of(PDB_to_SDF_settings)),
        )
    solvent_padding: float = field(
        default=0.9,
        validator=validators.instance_of(float),
        )
    box_shape: str = field(
        default="octahedron",
        validator=validators.instance_of(str),
        )

    def is_espaloma(self) -> bool:
        """
        Whether this parameterizer should use the espaloma protein remapping path.
        """
        name = os.path.basename(self.small_molecule_forcefield).lower()
        return name.startswith("espaloma") or name.endswith(".pt")
