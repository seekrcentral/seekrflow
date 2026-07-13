"""
modules/parameterize_workflow.py

The parameterize-side workflow object. This holds the state and methods used by
``parameterize.py`` and ``minimize_equilibrate.py`` to prepare a solvated,
parameterized molecular system (the input to the prepare-side ``Workflow``).

This code was relocated out of the now-removed per-system workflow classes
(``protein_ligand_seekr2`` etc.) so that parameterization is decoupled from the
prepare-side seekr3 ``Workflow``.
"""

import os
import typing

from attrs import define, field, validators, Factory
import parmed
import openmm
import openmm.unit as unit
import numpy as np

#import seekrtools.scripts.pqr_resid_for_each_atom as pqr_resid_for_each_atom

import seekrflow.modules.base as base
import seekrflow.modules.workflows.structures as workflow_structures
import seekrflow.modules.parameterize_structures as parameterize_structures
import seekrflow.modules.parameters_topology_structures \
    as parameters_topology_structures

LIGAND_PDB_FILENAME = "ligand.pdb"
RECEPTOR_PDB_FILENAME = "receptor.pdb"
DEFAULT_LIGAND_SDF_FILENAME = "ligand.sdf"


@define
class Parameterizer_information:
    """
    Information the parameterizer can use to prepare the system.
    """
    ligand_sdf_file: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    ligand_resname: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    receptor_ligand_pdb_filename: str = field(
        default="",
        validator=validators.instance_of(str),
        )


@define
class Parameterize_workflow:
    """
    Parameterize-side workflow holding the system definition and the methods
    used to split, parameterize, solvate, equilibrate, and emit PQR files.
    """
    type: typing.Literal["parameterize_workflow"] = "parameterize_workflow"
    solvated_system_for_md: workflow_structures.Solvated_system_for_md | None \
        = field(
        default=None,
        validator=validators.optional(
            validators.instance_of(workflow_structures.Solvated_system_for_md)),
        )
    ligand_indices: typing.List[int] = field(
        default=Factory(list),
        validator=validators.instance_of(list),
        )
    receptor_indices: typing.List[int] = field(
        default=Factory(list),
        validator=validators.instance_of(list),
        )
    receptor_pqr_filename_for_bd: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
        )
    ligand_pqr_filename_for_bd: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
        )
    parameterizer_information: Parameterizer_information = field(
        default=Factory(Parameterizer_information),
        validator=validators.instance_of(Parameterizer_information),
        )
    md_settings: workflow_structures.MD_settings | None = field(
        default=None,
        validator=validators.optional(
            validators.instance_of(workflow_structures.MD_settings)),
        )
    bd_settings: workflow_structures.BD_settings | None = field(
        default=None,
        validator=validators.optional(
            validators.instance_of(workflow_structures.BD_settings)),
        )

    def has_small_molecule_ligand(self) -> bool:
        """
        Return whether the workflow has a small molecule ligand.
        """
        return True

    def get_parameterizer_pdb_filename(self) -> str:
        """
        Get the filename of the PDB file to be used for parameterization.
        """
        assert self.parameterizer_information is not None, \
            "Parameterizer information is not set."
        assert self.parameterizer_information\
            .receptor_ligand_pdb_filename != "", \
            "parameterizer_information.receptor_ligand_pdb_filename "\
            "is not set, cannot parameterize."
        return self.parameterizer_information.receptor_ligand_pdb_filename

    def set_parameterizer_pdb_filename(self, filename: str) -> None:
        """
        Set the filename of the PDB file to be used for parameterization.
        """
        assert self.parameterizer_information is not None, \
            "Parameterizer information is not set."
        self.parameterizer_information\
            .receptor_ligand_pdb_filename = filename

    def get_parameterizer_ligand_pdb_filename(self) -> str:
        """
        Get the filename of the ligand PDB file used for parameterization.
        """
        return LIGAND_PDB_FILENAME

    def get_parameterizer_sdf_filename(self) -> str:
        """
        Get the filename of the SDF file to be used for parameterization.
        """
        assert self.parameterizer_information is not None, \
            "Parameterizer information is not set."
        if self.parameterizer_information.ligand_sdf_file != "":
            assert os.path.exists(
                self.parameterizer_information.ligand_sdf_file), \
                f"Ligand SDF file "\
                f"{self.parameterizer_information.ligand_sdf_file} "\
                "does not exist."
        return self.parameterizer_information.ligand_sdf_file

    def set_parameterizer_sdf_filename(self, filename: str) -> None:
        """
        Set the filename of the SDF file to be used for parameterization.
        """
        assert self.parameterizer_information is not None, \
            "Parameterizer information is not set."
        self.parameterizer_information.ligand_sdf_file = filename
        return

    def get_parameterizer_default_sdf_filename(self) -> str:
        """
        Get the default filename of the SDF file used for parameterization.
        """
        return DEFAULT_LIGAND_SDF_FILENAME

    def split_molecules(
            self,
            param_directory: str
            ) -> None:
        """
        Split the receptor-ligand complex into separate PDB files for the
        receptor and the ligand.
        """
        pdb_with_ligand = self.get_parameterizer_pdb_filename()
        full_structure = parmed.load_file(pdb_with_ligand, skip_bonds=True)
        receptor_filename = os.path.join(param_directory, RECEPTOR_PDB_FILENAME)
        ligand_filename = os.path.join(param_directory, LIGAND_PDB_FILENAME)
        assert len(self.ligand_indices) > 0, \
            "No ligand indices in seekrflow object."
        ligand_index_list: typing.List[int] = []
        resname = None
        for ligand_index in self.ligand_indices:
            # Parmed uses 1-based indexing for atom selections
            ligand_index_list.append(ligand_index + 1)
            if resname is None:
                resname = full_structure.atoms[ligand_index].residue.name
            else:
                assert resname == full_structure.atoms[ligand_index]\
                    .residue.name, \
                    "Ligand atoms belong to multiple residue names, "\
                    "cannot proceed."

        ligand_selection_str = f"@{','.join(map(str, ligand_index_list))}"
        receptor_selection_str = f"!{ligand_selection_str}"
        ligand_structure = full_structure[ligand_selection_str]
        ligand_structure.save(str(ligand_filename), overwrite=True)
        nonligand_structure = full_structure[receptor_selection_str]
        nonligand_structure.save(str(receptor_filename), overwrite=True)
        return

    def write_component_pqr_files(
            self,
            parmed_complex: parmed.Structure,
            param_directory: str
            ) -> None:
        """
        Write receptor and ligand PQR files (with MBondi2 radii) for use by the
        Brownian-dynamics engine.
        """
        assert self.parameterizer_information is not None, \
            "Parameterizer information is not set."

        # Assign radii
        base.assign_mbondi2_radii_to_parmed_structure(parmed_complex)

        # First, the ligand indices: find all residues the ligand atoms belong
        #  to and select the whole residue(s).
        ligand_resname_set = set()
        for ligand_index in self.ligand_indices:
            ligand_atom = parmed_complex.atoms[ligand_index]
            ligand_atom_resname = ligand_atom.residue.name
            ligand_resname_set.add(ligand_atom_resname)
        assert len(ligand_resname_set) >= 1, \
            "Ligand must have at least one residue."
        lig_selection_string = ":" + ",".join(ligand_resname_set)
        ligand_structure_all_resids = parmed_complex[lig_selection_string]

        # Choose only a single instance of the ligand in case there are
        #  multiple duplicates in the system.
        last_resnum = None
        last_chain = None
        lig_atom_number_list = []
        for j, atom in enumerate(ligand_structure_all_resids.atoms):
            if last_resnum is None:
                last_resnum = atom.residue.number
            if last_chain is None:
                last_chain = atom.residue.chain
            if (atom.residue.number != last_resnum) \
                    or (atom.residue.chain != last_chain):
                break
            lig_atom_number_list.append(j)
        ligand = ligand_structure_all_resids[:len(lig_atom_number_list)]

        # Next the receptor indices.
        not_ligand = parmed_complex[f"!{lig_selection_string}"]
        not_ligand_tmp_filename = os.path.join(
            param_directory, "not_ligand_tmp.pdb")
        not_ligand.save(not_ligand_tmp_filename, overwrite=True)
        fixer_settings = parameterize_structures.PDBFixer_settings(
            remove_extra_chains=True,
            find_missing_residues=False,
            find_and_replace_nonstandard_residues=False,
            find_and_add_missing_atoms=False,
            add_missing_hydrogens_pH=False,
        )
        receptor_tmp_filename = os.path.join(
            param_directory, "receptor_tmp.pdb")
        fixer_settings.run(not_ligand_tmp_filename, receptor_tmp_filename)
        full_receptor_pqr_filename = os.path.join(
            param_directory, "receptor.pqr")
        receptor_tmp_structure = parmed.load_file(receptor_tmp_filename)
        receptor = base.choose_only_receptor_atoms(
            parmed_complex, receptor_tmp_structure)

        full_ligand_pqr_filename_one_resid = os.path.join(
            param_directory, "ligand_one_resid.pqr")
        receptor.box = None  # Set box vectors to None for PQR
        ligand.box = None  # Set box vectors to None for PQR
        receptor.save(full_receptor_pqr_filename, overwrite=True)
        ligand.save(full_ligand_pqr_filename_one_resid, overwrite=True)

        full_ligand_pqr_filename = os.path.join(param_directory, "ligand.pqr")
        #pqr_resid_for_each_atom.pqr_resid_for_each_atom(
        #    full_ligand_pqr_filename_one_resid, full_ligand_pqr_filename)
        receptor_pqr_basename = os.path.basename(full_receptor_pqr_filename)
        ligand_pqr_basename = os.path.basename(full_ligand_pqr_filename)
        self.receptor_pqr_filename_for_bd = os.path.join(
            param_directory, receptor_pqr_basename)
        self.ligand_pqr_filename_for_bd = os.path.join(
            param_directory, ligand_pqr_basename)

    def create_complex(
            self,
            parameterizer: parameterize_structures.Parameterizer,
            physical_attributes: base.Physical_attributes,
            md_settings: workflow_structures.MD_settings,
            param_directory: str
            ) -> typing.Tuple[str, str]:
        """
        Create the solvated, parameterized complex structure.
        """
        receptor_filename = os.path.join(param_directory, RECEPTOR_PDB_FILENAME)
        ligand_filename = os.path.join(param_directory, LIGAND_PDB_FILENAME)
        assert os.path.exists(receptor_filename), \
            f"Receptor PDB file {receptor_filename} does not exist."
        assert os.path.exists(ligand_filename), \
            f"Ligand PDB file {ligand_filename} does not exist."
        apo_pdb_filename_no_ext = os.path.join(param_directory, "receptor")
        latest_pdb_filename = str(receptor_filename)

        if parameterizer.pdb_fixer_settings is not None:
            fixed_pdb_filename = str(apo_pdb_filename_no_ext) + "_fixed.pdb"
            parameterizer.pdb_fixer_settings.run(
                latest_pdb_filename, fixed_pdb_filename)
            latest_pdb_filename = fixed_pdb_filename
        if parameterizer.pdb2pqr_settings is not None:
            pdb2pqr_output_pqr_filename = str(apo_pdb_filename_no_ext) \
                + "_pdb2pqr.pqr"
            pdb2pqr_output_pdb_filename = str(apo_pdb_filename_no_ext) \
                + "_pdb2pqr.pdb"
            parameterizer.pdb2pqr_settings.run(
                latest_pdb_filename, pdb2pqr_output_pqr_filename,
                pdb2pqr_output_pdb_filename)
            latest_pdb_filename = pdb2pqr_output_pdb_filename

        # Load the receptor into openmm
        protein_topology, protein_md_topology, protein_positions = \
            base.make_protein_openmm_and_mdtraj_top(latest_pdb_filename)

        # Load the ligand
        ligand_topology, ligand_md_topology, ligand_positions, offmol = \
            base.make_ligand_openmm_and_mdtraj_top(
                self.parameterizer_information.ligand_sdf_file,
                ligand_filename,
                param_directory,
                draw_ligand=True,
                ligand_resname=self.parameterizer_information.ligand_resname
            )

        complex_md_topology = protein_md_topology.join(ligand_md_topology)
        complex_topology = complex_md_topology.to_openmm()
        # Ensure the number of atoms is consistent after merging
        n_atoms_total = complex_md_topology.n_atoms
        n_atoms_protein = protein_md_topology.n_atoms
        n_atoms_ligand = ligand_md_topology.n_atoms
        assert n_atoms_total == n_atoms_protein + n_atoms_ligand, \
            "Mismatch in atom numbers after merging."
        complex_positions = np.zeros([n_atoms_total, 3]) * unit.nanometers
        complex_positions[:n_atoms_protein, :] = protein_positions
        complex_positions[
            n_atoms_protein:n_atoms_protein + n_atoms_ligand, :] \
            = ligand_positions

        serialized_system_xml, output_pdb_filename \
            = workflow_structures.parameterize_and_check_complex(
                complex_topology,
                complex_positions,
                offmol,
                param_directory,
                parameterizer,
                physical_attributes,
                md_settings,
            )

        return serialized_system_xml, output_pdb_filename

    def minimize_equilibrate(
            self,
            physical_attributes: base.Physical_attributes,
            unsolvated_pdb_filename: str = "") -> str:
        """
        Minimize and equilibrate the solvated structure.
        """
        assert self.solvated_system_for_md is not None, \
            "Solvated system for MD is not set, cannot minimize and equilibrate."
        solvated_pdb_full_path = self.solvated_system_for_md.solvated_pdb
        complex_system_filename = self.solvated_system_for_md\
            .parameters_topology.system_filename
        with open(complex_system_filename, "r") as f:
            complex_system = openmm.XmlSerializer.deserialize(f.read())
        if unsolvated_pdb_filename == "":
            unsolvated_pdb_filename = self.parameterizer_information\
                .receptor_ligand_pdb_filename
        assert os.path.exists(unsolvated_pdb_filename), \
            f"Unsolvated PDB file {unsolvated_pdb_filename} does not exist."
        equilibrated_pdb_filename \
            = workflow_structures.minimize_and_equilibrate_complex(
                complex_system,
                solvated_pdb_full_path,
                unsolvated_pdb_filename,
                self.md_settings,
                physical_attributes,
                is_membrane_system=False,
            )
        return equilibrated_pdb_filename
