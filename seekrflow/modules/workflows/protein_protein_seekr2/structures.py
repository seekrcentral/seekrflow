"""
modules/workflows/protein_protein_seekr2/structures.py

Contain data structure classes used for seekrflow protein-protein
objects for use in seekr2. The classification of seekr2 implies
that the old HIDR program is used for MetaD/SMD.
"""

import os
import typing

from attrs import define, field, validators, Factory
import parmed
import mdtraj
import openmm.unit as unit
import numpy as np

import seekrflow.modules.base as base
import seekrflow.modules.workflows.structures as workflow_structures
import seekrflow.modules.parameterize_structures as parameterize_structures
import seekrflow.modules.parameters_topology_structures as parameters_topology_structures

PROTEIN1_PDB_FILENAME = "protein1.pdb"
PROTEIN2_PDB_FILENAME = "protein2.pdb"
PROTEIN1_SPLIT_PDB_FILENAME = "protein1_split.pdb"
PROTEIN2_SPLIT_PDB_FILENAME = "protein2_split.pdb"
PROTEIN1_NOH_SPLIT_PDB_FILENAME = "protein1_split_noh.pdb"
PROTEIN2_NOH_SPLIT_PDB_FILENAME = "protein2_split_noh.pdb"
SELECTION_NOH_FOR_PROTEIN_EXCEPTING_TERMINI = "(resname ACE NME) or (not element H)"


@define
class HIDR_settings_base:
    """
    Base class for HIDR settings.
    """
    type: typing.Literal["hidr_base"] = "hidr_base"


#class HIDR_settings_metaD(HIDR_settings_base):
@define
class HIDR_settings_metaD:
    """
    Settings for the HIDR calculation using metadynamics.
    """
    type: typing.Literal["hidr_metaD"] = "hidr_metaD"
    gaussian_height: float = field(
        default=1.0,
        validator=validators.instance_of(float),
        )
    gaussian_width: float = field(
        default=0.05,
        validator=validators.instance_of(float),
        )
    bias_factor: float = field(
        default=10.0,
        validator=validators.instance_of(float),
        )

#class HIDR_settings_SMD(HIDR_settings_base):
@define
class HIDR_settings_SMD:
    """
    Settings for the HIDR calculation using steered molecular dynamics.
    """
    type: typing.Literal["hidr_SMD"] = "hidr_SMD"
    restraint_force_constant: float = field(
        default=90000.0,
        validator=validators.instance_of(float),
        )
    translation_velocity: float = field(
        default=0.01,
        validator=validators.instance_of(float),
        )

#class MMVT_seekr_settings(Seekr_settings):
@define
class MMVT_seekr_settings:
    """
    Settings for the MMVT calculation.
    """
    type: typing.Literal["MMVT"] = "MMVT"
    md_output_interval: int = field(
        default=1000000,
        validator=validators.instance_of(int),
        )
    md_steps_per_anchor: int = field(
        default=100000000,
        validator=validators.instance_of(int),
        )
    cv_type: str = field(
        default="com_com_distance_rmsd",
        validator=validators.instance_of(str),
    )
    anchor_value_list: typing.List[typing.List[float]] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(list),
            iterable_validator=validators.instance_of(list),
        ))

@define
class Parameterizer_information:
    """
    Information the parameterizer can use to prepare the protein-protein
    seekr2 workflow.
    """
    protein_protein_pdb_filename: str = field(
        default="",
        validator=validators.instance_of(str),
        )

@define
class Protein_protein_seekr2_workflow:
    """
    Settings for the protein-protein seekr2 workflow.
    """
    type: typing.Literal["protein_protein_seekr2"] = "protein_protein_seekr2"
    solvated_system_for_md: workflow_structures.Solvated_system_for_md | None = field(
        default=None,
        validator=validators.optional(
            validators.instance_of(workflow_structures.Solvated_system_for_md)),
        )
    protein1_chain: str = field(
        default="A",
        validator=validators.instance_of(str),
        )
    protein2_chain: str = field(
        default="B",
        validator=validators.instance_of(str),
        )
    protein1_pqr_filename_for_bd: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
        )
    protein2_pqr_filename_for_bd: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
        )
    parameterizer_information: Parameterizer_information = field(
        default=Factory(Parameterizer_information),
        validator=validators.instance_of(Parameterizer_information),
        )
    hidr_settings: HIDR_settings_metaD | HIDR_settings_SMD = field(
        default=Factory(HIDR_settings_metaD),
        validator=validators.instance_of(HIDR_settings_metaD | HIDR_settings_SMD),
        )
    mmvt_settings: MMVT_seekr_settings = field(
        default=Factory(MMVT_seekr_settings),
        validator=validators.instance_of(MMVT_seekr_settings),
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
        Assuming no for this workflow type.
        """
        return False

    def get_parametrizer_pdb_filename(self) -> str:
        """
        Get the filename of the PDB file to be used for parameterization.
        """
        assert self.parameterizer_information is not None, \
            "Parameterizer information is not set."
        assert self.parameterizer_information\
            .protein_protein_pdb_filename != "", \
            "parameterizer_information.protein_protein_pdb_filename "\
            "is not set, cannot parameterize."
        return self.parameterizer_information.protein_protein_pdb_filename
    
    def set_parametrizer_pdb_filename(self, filename: str) -> None:
        """
        Set the filename of the PDB file to be used for parameterization.
        """
        assert self.parameterizer_information is not None, \
            "Parameterizer information is not set."
        self.parameterizer_information\
            .protein_protein_pdb_filename = filename
    
    def get_parametrizer_protein1_pdb_filename(self) -> str:
        """
        Get the filename of the protein1 PDB file to be used for parameterization.
        """
        return PROTEIN1_PDB_FILENAME

    def get_parametrizer_protein2_pdb_filename(self) -> str:
        """
        Get the filename of the protein2 PDB file to be used for parameterization.
        """
        return PROTEIN2_PDB_FILENAME
        
    def get_parametrizer_sdf_filename(self) -> str:
        """
        Not used in this type of workflow.
        """
        return ""
    
    #def set_parametrizer_sdf_filename(self, filename: str) -> None:
    #    """
    #    Set the filename of the SDF file to be used for parameterization.
    #    """
    #    assert self.parameterizer_information is not None, \
    #        "Parameterizer information is not set."
    #    self.parameterizer_information.ligand_sdf_file = filename
    #    return
    
    #def get_parametrizer_default_sdf_filename(self) -> str:
    #    """
    #    Get the default filename of the SDF file to be used for parameterization.
    #    """
    #    return DEFAULT_LIGAND_SDF_FILENAME

    def split_molecules(
            self,
            param_directory: str
            ) -> None:
        """
        Split the protein1-protein2 complex into separate PDB files.
        """
        pdb_with_ligand = self.get_parametrizer_pdb_filename()
        #full_structure = parmed.load_file(pdb_with_ligand, skip_bonds=True)
        full_structure = mdtraj.load(pdb_with_ligand)
        protein1_filename = os.path.join(param_directory, PROTEIN1_SPLIT_PDB_FILENAME)
        protein2_filename = os.path.join(param_directory, PROTEIN2_SPLIT_PDB_FILENAME)
        protein1_noh_filename = os.path.join(
            param_directory, PROTEIN1_NOH_SPLIT_PDB_FILENAME)
        protein2_noh_filename = os.path.join(
            param_directory, PROTEIN2_NOH_SPLIT_PDB_FILENAME)
        #protein1_selection_str = f"::{self.protein1_chain}"
        #protein2_selection_str = f"::{self.protein2_chain}"
        #protein1_structure = full_structure[protein1_selection_str]
        #protein2_structure = full_structure[protein2_selection_str]
        #protein1_structure = full_structure[self.protein1_chain,:,:]
        #protein2_structure = full_structure[self.protein2_chain,:,:]
        protein1_indices = []
        protein2_indices = []
        for i, chain in enumerate(full_structure.topology.chains):
            for atom in chain.atoms:
                if chain.chain_id == self.protein1_chain:
                    protein1_indices.append(atom.index)
                elif chain.chain_id == self.protein2_chain:
                    protein2_indices.append(atom.index)
        protein1_structure = full_structure.atom_slice(protein1_indices)
        protein2_structure = full_structure.atom_slice(protein2_indices)
        protein1_structure.save_pdb(protein1_filename)
        protein2_structure.save_pdb(protein2_filename)
        protein1_structure_noh_indices = protein1_structure.topology.select(
            SELECTION_NOH_FOR_PROTEIN_EXCEPTING_TERMINI)
        protein2_structure_noh_indices = protein2_structure.topology.select(
            SELECTION_NOH_FOR_PROTEIN_EXCEPTING_TERMINI)
        protein1_structure_noh = protein1_structure.atom_slice(protein1_structure_noh_indices)
        protein2_structure_noh = protein2_structure.atom_slice(protein2_structure_noh_indices)
        protein1_structure_noh.save_pdb(protein1_noh_filename)
        protein2_structure_noh.save_pdb(protein2_noh_filename)
        return

    def write_component_pqr_files(
            self, 
            parmed_complex: parmed.Structure,
            param_directory: str
            ) -> None:
        """
        Get the atom indices for the components in the system.

        In the case of this workflow, return the receptor and ligand
        atom indices.
        """
        assert self.parameterizer_information is not None, \
            "Parameterizer information is not set."
        
        # Assign radii
        base.assign_mbondi2_radii_to_parmed_structure(parmed_complex)

        #protein1_selection_str = f"::{self.protein1_chain}"
        #protein2_selection_str = f"::{self.protein2_chain}"
        #protein1_structure_all_resids = parmed_complex[protein1_selection_str]
        #protein2_structure_all_resids = parmed_complex[protein2_selection_str]

        protein1_structure_all_resids = parmed_complex[self.protein1_chain,:,:]
        protein2_structure_all_resids = parmed_complex[self.protein2_chain,:,:]

        protein1_tmp_filename = os.path.join(param_directory, "protein1_tmp.pdb")
        protein2_tmp_filename = os.path.join(param_directory, "protein2_tmp.pdb")
        protein1_structure_all_resids.save(protein1_tmp_filename, overwrite=True)
        protein2_structure_all_resids.save(protein2_tmp_filename, overwrite=True)
        fixer_settings = parameterize_structures.PDBFixer_settings(
            remove_extra_chains=True,
            find_missing_residues=False,
            find_and_replace_nonstandard_residues=False,
            find_and_add_missing_atoms=False,
            add_missing_hydrogens_pH=False,
        )
        protein1_filename = os.path.join(param_directory, "protein1.pdb")
        protein2_filename = os.path.join(param_directory, "protein2.pdb")
        fixer_settings.run(protein1_tmp_filename, protein1_filename)
        fixer_settings.run(protein2_tmp_filename, protein2_filename)
        full_protein1_pqr_filename = os.path.join(param_directory, "protein1.pqr")
        full_protein2_pqr_filename = os.path.join(param_directory, "protein2.pqr")
        protein1_structure = parmed.load_file(protein1_filename)
        protein2_structure = parmed.load_file(protein2_filename)
        protein1_structure.box = None  # Set box vectors to None for PQR
        protein2_structure.box = None  # Set box vectors to None for PQR
        protein1_structure.save(full_protein1_pqr_filename, overwrite=True)
        protein2_structure.save(full_protein2_pqr_filename, overwrite=True)
        protein1_pqr_basename = os.path.basename(full_protein1_pqr_filename)
        protein2_pqr_basename = os.path.basename(full_protein2_pqr_filename)
        self.protein1_pqr_filename_for_bd = os.path.join(
            param_directory, protein1_pqr_basename)
        self.protein2_pqr_filename_for_bd = os.path.join(
            param_directory, protein2_pqr_basename)
        return
        
    def create_complex(
            self,
            parameterizer: parameterize_structures.Parameterizer,
            physical_attributes: base.Physical_attributes,
            md_settings: workflow_structures.MD_settings,
            param_directory: str
            ) -> typing.Tuple[str, str]:
        """
        Create the complex structure.
        """
        #protein1_filename = os.path.join(param_directory, PROTEIN1_SPLIT_PDB_FILENAME)
        #protein2_filename = os.path.join(param_directory, PROTEIN2_SPLIT_PDB_FILENAME)
        protein1_filename = os.path.join(param_directory, PROTEIN1_NOH_SPLIT_PDB_FILENAME)
        protein2_filename = os.path.join(param_directory, PROTEIN2_NOH_SPLIT_PDB_FILENAME)
        assert os.path.exists(protein1_filename), \
            f"Protein 1 PDB file {protein1_filename} does not exist."
        assert os.path.exists(protein2_filename), \
            f"Protein 2 PDB file {protein2_filename} does not exist."
        protein1_pdb_filename_no_ext = os.path.join(param_directory, "protein1")
        latest_pdb_filename = str(protein1_filename)
        
        if parameterizer.pdb_fixer_settings is not None:
            fixed_pdb_filename = str(protein1_pdb_filename_no_ext) + "_fixed.pdb"
            parameterizer.pdb_fixer_settings.run(latest_pdb_filename, fixed_pdb_filename)
            latest_pdb_filename = fixed_pdb_filename
        if parameterizer.pdb2pqr_settings is not None:
            pdb2pqr_output_pqr_filename = str(protein1_pdb_filename_no_ext) + "_pdb2pqr.pqr"
            pdb2pqr_output_pdb_filename = str(protein1_pdb_filename_no_ext) + "_pdb2pqr.pdb"
            parameterizer.pdb2pqr_settings.run(latest_pdb_filename, pdb2pqr_output_pqr_filename,
                                        pdb2pqr_output_pdb_filename)
            latest_pdb_filename = pdb2pqr_output_pdb_filename

        protein1_topology, protein1_md_topology, protein1_positions = \
            base.make_protein_openmm_and_mdtraj_top(latest_pdb_filename)
        
        protein2_pdb_filename_no_ext = os.path.join(param_directory, "protein2")
        latest_pdb_filename = str(protein2_filename)

        if parameterizer.pdb_fixer_settings is not None:
            fixed_pdb_filename = str(protein2_pdb_filename_no_ext) + "_fixed.pdb"
            parameterizer.pdb_fixer_settings.run(latest_pdb_filename, fixed_pdb_filename)
            latest_pdb_filename = fixed_pdb_filename
        if parameterizer.pdb2pqr_settings is not None:
            pdb2pqr_output_pqr_filename = str(protein2_pdb_filename_no_ext) + "_pdb2pqr.pqr"
            pdb2pqr_output_pdb_filename = str(protein2_pdb_filename_no_ext) + "_pdb2pqr.pdb"
            parameterizer.pdb2pqr_settings.run(latest_pdb_filename, pdb2pqr_output_pqr_filename,
                                        pdb2pqr_output_pdb_filename)
            latest_pdb_filename = pdb2pqr_output_pdb_filename

        protein2_topology, protein2_md_topology, protein2_positions = \
            base.make_protein_openmm_and_mdtraj_top(latest_pdb_filename)

        complex_md_topology = protein1_md_topology.join(protein2_md_topology)
        complex_topology = complex_md_topology.to_openmm()
        # Ensure the number of atoms is consistent after merging
        n_atoms_total = complex_md_topology.n_atoms
        n_atoms_protein1 = protein1_md_topology.n_atoms
        n_atoms_protein2 = protein2_md_topology.n_atoms
        assert n_atoms_total == n_atoms_protein1 + n_atoms_protein2, \
            "Mismatch in atom numbers after merging."
        complex_positions = np.zeros([n_atoms_total, 3]) * unit.nanometers
        complex_positions[:n_atoms_protein1, :] = protein1_positions
        complex_positions[n_atoms_protein1:n_atoms_protein1 + n_atoms_protein2, :] \
            = protein2_positions
        
        serialized_xml, output_pdb_filename = workflow_structures.parametrize_and_check_complex(
            complex_topology,
            complex_positions,
            None,
            param_directory,
            parameterizer,
            physical_attributes,
            md_settings,
        )

        return serialized_xml, output_pdb_filename
    
    # TODO: fill out later
    def create_mmvt_protein_protein_com_com_rmsd_model_input_seekr2(
            self,
            root_directory: str,
            param_directory: str,
            physical_attributes: base.Physical_attributes,
            alpha_carbon_ligand_threshold: float = 0.6
            ) -> "seekr2_common_prepare.Model_input":
        """
        Create the input object for the MMVT model in SEEKR2.
        """
        import seekr2.modules.common_base as seekr2_base
        import seekr2.modules.common_cv as seekr2_common_cv
        import seekr2.modules.common_prepare as seekr2_common_prepare
        import seekrflow.modules.cvs as cvs

        # Copy starting PDB file over to the prepare directory
        model_input = seekr2_common_prepare.Model_input()
        model_input.calculation_type = "mmvt"
        model_input.calculation_settings = seekr2_common_prepare.MMVT_input_settings()
        model_input.calculation_settings.md_output_interval \
            = self.mmvt_settings.md_output_interval
        model_input.calculation_settings.md_steps_per_anchor \
            = self.mmvt_settings.md_steps_per_anchor
        model_input.temperature = physical_attributes.temperature
        if physical_attributes.pressure is None:
            model_input.pressure = 1.0
            model_input.ensemble = "nvt"
        else:
            model_input.pressure = physical_attributes.pressure
            model_input.ensemble = "npt"
        model_input.root_directory = root_directory
        model_input.md_program = "openmm"
        model_input.constraints = "HBonds"
        model_input.rigidWater = True
        if physical_attributes.hydrogen_mass == 1.008:
            model_input.hydrogenMass = None
        else:
            model_input.hydrogenMass = physical_attributes.hydrogen_mass
        model_input.timestep = self.md_settings.stepsize
        model_input.nonbonded_cutoff = self.md_settings.nonbonded_cutoff

        if self.solvated_system_for_md is None:
            starting_pdb_basename = "complex-equil.pdb"
            starting_pdb_full_path = os.path.join(param_directory, starting_pdb_basename)
        elif self.solvated_system_for_md.solvated_pdb == "":
            starting_pdb_basename = "complex-equil.pdb"
            starting_pdb_full_path = os.path.join(param_directory, starting_pdb_basename)
        else:
            starting_pdb_basename = os.path.basename(
                self.solvated_system_for_md.solvated_pdb)
            starting_pdb_full_path = self.solvated_system_for_md.solvated_pdb

        cv_input1 = seekr2_common_cv.Spherical_cv_input()
        receptor_atoms, ligand_atoms = cvs.get_receptor_ligand_com_com_selections(
            self, starting_pdb_full_path, alpha_carbon_ligand_threshold)
        if len(self.receptor_indices) == 0:
            self.receptor_indices = receptor_atoms
        if len(self.ligand_indices) == 0:
            self.ligand_indices = ligand_atoms
        cv_input1.group1 = self.receptor_indices
        cv_input1.group2 = self.ligand_indices
        cv_input1.input_anchors = []
        radius_list = self.mmvt_settings.anchor_radius_list
        assert len(radius_list) > 0, \
            "Anchor radius list is empty in the input JSON file."
        if self.solvated_system_for_md is None:
            self.solvated_system_for_md = workflow_structures.Solvated_system_for_md()
            self.solvated_system_for_md.parameters_topology \
                = parameters_topology_structures.Openmm_system()
            self.solvated_system_for_md.parameters_topology.system_filename \
                = os.path.join(param_directory, "complex.xml")
            self.solvated_system_for_md.solvated_pdb = starting_pdb_full_path

        for i, radius in enumerate(radius_list):
            input_anchor = seekr2_common_cv.Spherical_cv_anchor()
            input_anchor.radius = radius
            if self.solvated_system_for_md.parameters_topology.type == "Amber":
                amber_prmtop_filename = self.solvated_system_for_md\
                    .parameters_topology.prmtop_filename
                cvs.assign_amber_params(input_anchor, amber_prmtop_filename, "")
            elif self.solvated_system_for_md.parameters_topology.type == "Gromacs":
                raise NotImplementedError(
                    "Gromacs topology type is not implemented yet.")
            elif self.solvated_system_for_md.parameters_topology.type == "Charmm":
                charmm_psf_filename = self.solvated_system_for_md\
                    .parameters_topology.psf_filename
                charmm_ff_filenames = self.solvated_system_for_md\
                    .parameters_topology.param_filename_list
                cvs.assign_charmm_params(input_anchor, charmm_psf_filename, 
                                    charmm_ff_filenames, "")
            elif self.solvated_system_for_md.parameters_topology.type == "OpenMM_forcefield":
                built_in_ff_list = self.solvated_system_for_md.parameters_topology\
                    .built_in_forcefield_filenames
                custom_ff_list = self.solvated_system_for_md.parameters_topology\
                    .custom_forcefield_filenames
                cvs.assign_forcefield_params(input_anchor, built_in_ff_list, 
                                        custom_ff_list, "")
            elif self.solvated_system_for_md.parameters_topology.type == "OpenMM_system":
                system_filename = self.solvated_system_for_md\
                    .parameters_topology.system_filename
                cvs.assign_system_params(input_anchor, system_filename, "")
            else:
                raise Exception(
                    "Type not supported for seekrflow.md_parameters_topology: "\
                    f"{self.solvated_system_for_md.parameters_topology.type}")

            if i == 0:
                input_anchor.bound_state = True
            else:
                input_anchor.bound_state = False
                
            if i == len(radius_list)-1:
                input_anchor.bulk_anchor = True
            else:
                input_anchor.bulk_anchor = False
        
            cv_input1.input_anchors.append(input_anchor)
        
        model_input.cv_inputs = [cv_input1]
        if self.bd_settings is not None:
            # TODO: fill out from parameterize
            if self.receptor_pqr_filename_for_bd is None:
                self.receptor_pqr_filename_for_bd = os.path.join(
                    param_directory, "receptor.pqr")

            if self.ligand_pqr_filename_for_bd is None:
                self.ligand_pqr_filename_for_bd = os.path.join(
                    param_directory, "ligand.pqr")

            bd_rec_indices, bd_lig_indices = cvs.get_bd_receptor_ligand_selections(
                self.receptor_pqr_filename_for_bd, self.ligand_pqr_filename_for_bd,
                starting_pdb_full_path, self.receptor_indices, 
                self.ligand_indices)
            
            # TODO: Check for BD engine type - Browndye vs. SDA
            # If self.bd_settings.engine == "Browndye":
            model_input.browndye_settings_input \
                = seekr2_common_prepare.Browndye_settings_input()
            model_input.browndye_settings_input.binary_directory \
                = self.bd_settings.binary_directory
            model_input.browndye_settings_input.receptor_pqr_filename \
                = self.receptor_pqr_filename_for_bd
            model_input.browndye_settings_input.ligand_pqr_filename \
                = self.ligand_pqr_filename_for_bd
            model_input.browndye_settings_input.apbs_grid_spacing \
                = base.APBS_GRID_SPACING
            model_input.browndye_settings_input.receptor_indices = bd_rec_indices
            model_input.browndye_settings_input.ligand_indices = bd_lig_indices
            
            ion1 = seekr2_base.Ion()
            ion1.radius = 1.2
            ion1.charge = -1.0
            ion1.conc = physical_attributes.ionic_strength
            ion2 = seekr2_base.Ion()
            ion2.radius = 0.9
            ion2.charge = 1.0
            ion2.conc = physical_attributes.ionic_strength
            model_input.browndye_settings_input.ions = [ion1, ion2]
            model_input.browndye_settings_input.num_b_surface_trajectories \
                = self.bd_settings.num_trajectories
            model_input.browndye_settings_input.n_threads = 1
        else:
            model_input.browndye_settings_input = None

        return model_input
