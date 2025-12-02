"""
modules/parameterize_structures.py

Contain data structure classes used for seekrflow parameterize-
related objects.
"""

import os
import typing

from attrs import define, field, validators, Factory
import openmm.app as openmm_app

@define
class PDBFixer_settings:
    """
    Settings for PDBFixer.
    """
    remove_extra_chains: bool = field(
        default=True,
        validator=validators.instance_of(bool),
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
        Run PDBFixer on the given PDB file and return the fixed PDB filename.
        """
        from pdbfixer import PDBFixer
        fixer = PDBFixer(filename=input_pdb_filename)
        numChains = len(list(fixer.topology.chains()))
        if self.remove_extra_chains:
            fixer.removeChains(range(1, numChains))
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
    Settings for PDB2PQR.
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
        Run PDB2PQR on the given PDB file and return the PQR and PDB resulting
        files with the hydrogens properly added.
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
class Parameterizer:
    """
    The parameterizer object contains the inputs needed to run the
    parameterization of a protein-ligand complex.
    """
    forcefield: str = field(
        #default="espaloma-0.3.2.pt",
        default="gaff-2.11",
        validator=validators.instance_of(str),
        # TODO: possible choices?
        )
    water_model: str = field(
        default="tip3p",
        validator=validators.instance_of(str),
        # TODO: possible choices?
        )
    auxiliary_forcefield_files: typing.List[str] = field(
        default=Factory(list),
        validator=validators.instance_of(list),
        )
    pdb_fixer_settings: PDBFixer_settings | None = field(
        default=Factory(PDBFixer_settings),
        validator=validators.optional(validators.instance_of(PDBFixer_settings)),
        )
    pdb2pqr_settings: PDB2PQR_settings | None = field(
        default=Factory(PDB2PQR_settings),
        validator=validators.optional(validators.instance_of(PDB2PQR_settings)),
        )
    solvent_padding: float | None = field(
        default=0.9,
        validator=validators.instance_of(float),
        )
    box_shape: str = field(
        default="octahedron",
        validator=validators.instance_of(str),
        )
    