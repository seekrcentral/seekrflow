"""
modules/workflows/protein_ligand_membrane_seekr2/structures.py

Contain data structure classes used for seekrflow protein-ligand
objects for use in seekr2. The classification of seekr2 implies
that the old HIDR program is used for MetaD/SMD.
"""

import os
import typing

from attrs import define, field, validators, Factory
import parmed
import openmm.unit as unit
import numpy as np

import seekrtools.scripts.pqr_resid_for_each_atom as pqr_resid_for_each_atom

import seekrflow.modules.base as base
import seekrflow.modules.workflows.structures as workflow_structures
import seekrflow.modules.parameterize_structures as parameterize_structures
import seekrflow.modules.parameters_topology_structures as parameters_topology_structures

def separate_protein_membrane_from_complex(parmed_structure):
    """
    Use graph connectivity to identify protein and membrane components.

    Args:
        parmed_structure: parmed.Structure object

    Returns:
        tuple: (protein_atom_indices, membrane_atom_indices)
    """
    import networkx as nx

    # Standard amino acid residue names (including common variants)
    protein_resnames = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "HIS", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "MET",
        "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        # Common terminal and modified residues
        "ACE", "NME", "NHE", "CT3", "NTER", "CTER"
    }

    # Common cholesterol residue names
    cholesterol_resnames = {"CHOL", "CHL1", "CHLS", "CHLE"}

    # Build connectivity graph
    graph = nx.Graph()

    # Add all atoms as nodes
    for i, atom in enumerate(parmed_structure.atoms):
        graph.add_node(i, atom=atom)

    # Add bonds as edges
    for bond in parmed_structure.bonds:
        atom1_idx = parmed_structure.atoms.index(bond.atom1)
        atom2_idx = parmed_structure.atoms.index(bond.atom2)
        graph.add_edge(atom1_idx, atom2_idx)

    # Find connected components (molecules)
    connected_components = list(nx.connected_components(graph))

    protein_atom_indices = []
    membrane_atom_indices = []

    for component in connected_components:
        component_atoms = [parmed_structure.atoms[i] for i in component]

        # Check if this component is a protein
        is_protein = False
        component_resnames = {atom.residue.name for atom in component_atoms}

        # If any standard amino acid residue is found, assume it's protein
        if component_resnames.intersection(protein_resnames):
            is_protein = True

        if is_protein:
            protein_atom_indices.extend(list(component))
            continue

        # Check if this component is a membrane lipid
        is_membrane = False

        # Method 1: Check for phosphorus (likely phospholipid)
        has_phosphorus = any(atom.element == 15 or atom.name == 'P'
                             for atom in component_atoms)  # Element 15 = P
        if has_phosphorus and len(component) > 10:
            # Must be substantial molecule
            is_membrane = True

        # Method 2: Check for cholesterol by residue name
        if component_resnames.intersection(cholesterol_resnames):
            is_membrane = True

        # Method 3: Check for cholesterol by composition (C27H46O)
        if not is_membrane and len(component) > 50:
            # Cholesterol has 74 atoms
            element_counts = {}
            for atom in component_atoms:
                if hasattr(atom, 'element_name'):
                    element = atom.element_name
                else:
                    element = atom.name[0]
                element_counts[element] = element_counts.get(element, 0) + 1

            # Check if composition matches cholesterol (allow some flexibility)
            c_count = element_counts.get('C', 0)
            h_count = element_counts.get('H', 0)
            o_count = element_counts.get('O', 0)

            cholesterol_composition = (
                25 <= c_count <= 29 and    # Allow some variation
                40 <= h_count <= 50 and
                o_count == 1 and
                len(element_counts) <= 3   # Should only have C, H, O
            )
            if cholesterol_composition:
                is_membrane = True

        if is_membrane:
            membrane_atom_indices.extend(list(component))

    return protein_atom_indices, membrane_atom_indices

def convert_hetatm_to_atom(pqr_filename: str) -> None:
    """
    Convert all HETATM records to ATOM records in a PQR file.
    
    Args:
        pqr_filename: Path to the PQR file to modify in-place
    """
    with open(pqr_filename, 'r') as f:
        lines = f.readlines()
    
    with open(pqr_filename, 'w') as f:
        for line in lines:
            if line.startswith("HETATM"):
                line = "ATOM  " + line[6:]
            f.write(line)

from seekrflow.modules.workflows.protein_ligand_seekr2.structures import \
    MMVT_seekr_settings, HIDR_settings_metaD, HIDR_settings_SMD

LIGAND_PDB_FILENAME = "ligand.pdb"
RECEPTOR_PDB_FILENAME = "receptor_membrane.pdb"
DEFAULT_LIGAND_SDF_FILENAME = "ligand.sdf"

@define
class Parameterizer_information:
    """
    Information the parameterizer can use to prepare the ligand-protein
    seekr2 workflow.
    """
    ligand_sdf_file: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    ligand_resname: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    receptor_ligand_membrane_pdb_filename: str = field(
        default="",
        validator=validators.instance_of(str),
        )

@define
class Protein_ligand_membrane_seekr2_workflow:
    """
    Settings for the protein-ligand-membrane seekr2 workflow.
    """
    type: typing.Literal["protein_ligand_membrane_seekr2"] = "protein_ligand_membrane_seekr2"
    solvated_system_for_md: workflow_structures.Solvated_system_for_md | None = field(
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
    parameterizer_information: Parameterizer_information | None = field(
        default=None,
        validator=validators.optional(
            validators.instance_of(Parameterizer_information)),
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
        Assuming yes for this workflow type.
        """
        return True

    def get_parameterizer_pdb_filename(self) -> str:
        """
        Get the filename of the PDB file to be used for parameterization.
        """
        assert self.parameterizer_information is not None, \
            "Parameterizer information is not set."
        assert self.parameterizer_information\
            .receptor_ligand_membrane_pdb_filename != "", \
            "parameterizer_information.receptor_ligand_membrane_pdb_filename "\
            "is not set, cannot parameterize."
        return self.parameterizer_information.receptor_ligand_membrane_pdb_filename

    def set_parameterizer_pdb_filename(self, filename: str) -> None:
        """
        Set the filename of the PDB file to be used for parameterization.
        """
        assert self.parameterizer_information is not None, \
            "Parameterizer information is not set."
        self.parameterizer_information\
            .receptor_ligand_membrane_pdb_filename = filename

    def get_parameterizer_ligand_pdb_filename(self) -> str:
        """
        Get the filename of the ligand PDB file to be used for parameterization.
        """
        return LIGAND_PDB_FILENAME
        
    def get_parameterizer_sdf_filename(self) -> str:
        """
        Get the filename of the SDF file to be used for parameterization.
        """
        assert self.parameterizer_information is not None, \
            "Parameterizer information is not set."
        if self.parameterizer_information.ligand_sdf_file != "":
            assert os.path.exists(self.parameterizer_information.ligand_sdf_file), \
            f"Ligand SDF file {self.parameterizer_information.ligand_sdf_file} does not exist."
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
        Get the default filename of the SDF file to be used for parameterization.
        """
        return DEFAULT_LIGAND_SDF_FILENAME

    def split_molecules(
            self,
            param_directory: str
            ) -> None:
        """
        Split the receptor-ligand complex into separate PDB files for the receptor and the ligand.
        """
        pdb_with_ligand = self.get_parameterizer_pdb_filename()
        full_structure = parmed.load_file(pdb_with_ligand, skip_bonds=True)
        receptor_filename = os.path.join(param_directory, RECEPTOR_PDB_FILENAME)
        ligand_filename = os.path.join(param_directory, LIGAND_PDB_FILENAME)
        assert len(self.ligand_indices) > 0, \
            "No ligand indices in seekrflow object."
        #ligand_serial_list:  typing.List[int] = []
        ligand_index_list: typing.List[int] = []
        resname = None
        for ligand_index in self.ligand_indices:
            # Parmed uses 1-based indexing for atom selections
            ligand_index_list.append(ligand_index+1)
            #ligand_serial = full_structure.atoms[ligand_index].number
            #ligand_serial_list.append(ligand_serial)
            if resname is None:
                resname = full_structure.atoms[ligand_index].residue.name
            else:
                assert resname == full_structure.atoms[ligand_index].residue.name, \
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
            ) -> list[int]:
        """
        Get the atom indices for the components in the system.

        In the case of this workflow, return the receptor and ligand
        atom indices.
        """
        assert self.parameterizer_information is not None, \
            "Parameterizer information is not set."
        
        # Assign radii
        base.assign_mbondi2_radii_to_parmed_structure(parmed_complex)
        
        # First, the ligand indices.
        # Find all residues that the ligand atoms belong to and select the whole
        #  residue(s).
        ligand_resname_set = set()
        for ligand_index in self.ligand_indices:
            ligand_atom = parmed_complex.atoms[ligand_index]
            ligand_atom_resname = ligand_atom.residue.name
            ligand_resname_set.add(ligand_atom_resname)
        assert len(ligand_resname_set) >= 1, \
            "Ligand must have at least one residue."
        lig_selection_string = ":"+",".join(ligand_resname_set)
        ligand_structure_all_resids = parmed_complex[lig_selection_string]

        # This is an attempt to choose only one single instance of the ligand
        # in case there are multiple duplicates in the system.
        last_resnum = None
        last_chain = None
        lig_atom_number_list = []
        for j, atom in enumerate(ligand_structure_all_resids.atoms):
            if last_resnum is None:
                last_resnum = atom.residue.number
            if last_chain is None:
                last_chain = atom.residue.chain
            if (atom.residue.number != last_resnum) or (atom.residue.chain != last_chain):
                break
            lig_atom_number_list.append(j)
        ligand = ligand_structure_all_resids[:len(lig_atom_number_list)]
        
        # Next the receptor indices.
        not_ligand = parmed_complex[f"!{lig_selection_string}"]
        not_ligand_tmp_filename = os.path.join(param_directory, "not_ligand_tmp.pdb")
        not_ligand.save(not_ligand_tmp_filename, overwrite=True)
        full_receptor_pqr_filename = os.path.join(param_directory, "receptor_membrane.pqr")
        not_ligand_structure = parmed.load_file(
            not_ligand_tmp_filename)
        
        # Use graph-based approach to separate protein and membrane components
        protein_indices, membrane_indices = separate_protein_membrane_from_complex(
            not_ligand_structure)

        # Create receptor structure (protein + membrane)
        all_receptor_indices = protein_indices + membrane_indices

        if all_receptor_indices:
            # Convert to atom numbers for selection
            receptor_atom_numbers = [
                not_ligand_structure.atoms[i].number
                for i in all_receptor_indices
            ]
            receptor_selection = "@" + ",".join(map(str, receptor_atom_numbers))
            receptor = not_ligand_structure[receptor_selection]
            print(f"Found receptor with {len(receptor.atoms)} atoms "
                  f"({len(protein_indices)} protein, "
                  f"{len(membrane_indices)} membrane)")
        else:
            # Fallback: use all non-ligand atoms except water/ions
            print("Warning: Could not identify protein/membrane components. "
                  "Using all non-water/ion atoms as receptor.")
            receptor = not_ligand_structure[":!WAT,HOH,NA,CL,K,MG,CA,ZN"]
        
        full_ligand_pqr_filename_one_resid = os.path.join(param_directory, "ligand_one_resid.pqr")
        receptor.box = None  # Set box vectors to None for PQR
        ligand.box = None  # Set box vectors to None for PQR
        base.assign_mbondi2_radii_to_parmed_structure(receptor)

        receptor.save(full_receptor_pqr_filename, overwrite=True)
        ligand.save(full_ligand_pqr_filename_one_resid, overwrite=True)
        # For each of these files, change any HETATM lines to ATOM lines
        convert_hetatm_to_atom(full_receptor_pqr_filename)
        convert_hetatm_to_atom(full_ligand_pqr_filename_one_resid)

        full_ligand_pqr_filename = os.path.join(param_directory, "ligand.pqr")
        pqr_resid_for_each_atom.pqr_resid_for_each_atom(
            full_ligand_pqr_filename_one_resid, full_ligand_pqr_filename)
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
        Create the complex structure.
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
            parameterizer.pdb_fixer_settings.run(latest_pdb_filename, fixed_pdb_filename)
            latest_pdb_filename = fixed_pdb_filename
        if parameterizer.pdb2pqr_settings is not None:
            pdb2pqr_output_pqr_filename = str(apo_pdb_filename_no_ext) + "_pdb2pqr.pqr"
            pdb2pqr_output_pdb_filename = str(apo_pdb_filename_no_ext) + "_pdb2pqr.pdb"
            parameterizer.pdb2pqr_settings.run(latest_pdb_filename, pdb2pqr_output_pqr_filename,
                                        pdb2pqr_output_pdb_filename)
            latest_pdb_filename = pdb2pqr_output_pdb_filename

        # TODO: need functions to parameterize proteins vs. small molecules - call the functions to
        # enforce DRY
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
        assert n_atoms_total == n_atoms_protein + n_atoms_ligand, "Mismatch in atom numbers after merging."
        complex_positions = np.zeros([n_atoms_total, 3]) * unit.nanometers
        complex_positions[:n_atoms_protein, :] = protein_positions
        complex_positions[n_atoms_protein:n_atoms_protein + n_atoms_ligand, :] = ligand_positions

        serialized_xml, output_pdb_filename = workflow_structures.parameterize_and_check_complex(
            complex_topology,
            complex_positions,
            offmol,
            param_directory,
            parameterizer,
            physical_attributes,
            md_settings,
        )

        return serialized_xml, output_pdb_filename
        
    def create_mmvt_receptor_ligand_com_com_model_input_seekr2(
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
        import seekr2.modules.common_sim_sda as seekr2_common_sim_sda
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
            #model_input.pressure = 1.0
            #model_input.ensemble = "nvt"
            raise Exception("Pressure must be set for membrane simulation.")
        else:
            model_input.pressure = physical_attributes.pressure
            model_input.ensemble = "npt_membrane"
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
                    param_directory, "receptor_membrane.pqr")
                
            if self.ligand_pqr_filename_for_bd is None:
                self.ligand_pqr_filename_for_bd = os.path.join(
                    param_directory, "ligand.pqr")
                
            bd_rec_indices, bd_lig_indices = cvs.get_bd_receptor_ligand_selections(
                self.receptor_pqr_filename_for_bd, self.ligand_pqr_filename_for_bd,
                starting_pdb_full_path, self.receptor_indices, 
                self.ligand_indices)
            
            model_input.sda_settings_input \
                = seekr2_common_prepare.SDA_settings_input()
            model_input.sda_settings_input.sda_bin_dir \
                = self.bd_settings.binary_directory
            model_input.sda_settings_input.sda_auxi_dir \
                = self.bd_settings.auxiliary_directory
            model_input.sda_settings_input.hydropro_dir \
                = self.bd_settings.hydropro_directory
            model_input.sda_settings_input.type_calculation = "sda_2proteins"
            model_input.sda_settings_input.receptor_indices = bd_rec_indices
            model_input.sda_settings_input.ligand_indices = bd_lig_indices
            model_input.sda_settings_input.geom_type = "box"
            if os.path.exists(starting_pdb_full_path):
                pdb_structure = parmed.load_file(starting_pdb_full_path, skip_bonds=True)
                coords = pdb_structure.get_coordinates()  # shape: (n_atoms, 3)
                xmin = np.min(coords[0, :, 0])
                xmax = np.max(coords[0, :, 0])
                ymin = np.min(coords[0, :, 1])
                ymax = np.max(coords[0, :, 1])
                zmin = np.min(coords[0, :, 2])
                zmax = np.max(coords[0, :, 2])
                xwidth = xmax - xmin
                ywidth = ymax - ymin
                zwidth = zmax - zmin
                padding = 1.0  # A
                model_input.sda_settings_input.xmin = -0.5*xwidth - padding
                model_input.sda_settings_input.xmax = 0.5*xwidth + padding
                model_input.sda_settings_input.ymin = -0.5*ywidth - padding
                model_input.sda_settings_input.ymax = 0.5*ywidth + padding
                model_input.sda_settings_input.zmin = -0.5*zwidth - padding
                model_input.sda_settings_input.zmax = 0.5*zwidth + padding
            else:
                raise Exception(
                    f"Starting PDB file {starting_pdb_full_path} does not exist - "
                    "cannot set bounding box coordinates for BD.")
                
            receptor_solute = seekr2_common_sim_sda.Solute()
            receptor_solute.type = "Protein"
            receptor_solute.apbs_grid_dime = 257
            receptor_solute.apbs_grid_spacing = base.APBS_GRID_SPACING
            receptor_solute.dielectric = 2.0
            receptor_solute.pqr_filename = self.receptor_pqr_filename_for_bd
            receptor_solute.solute_grid = seekr2_common_sim_sda.Solute_Grid()
            receptor_solute.solute_grid.nb_solute = 1
            receptor_solute.solute_grid.surface = "yes"
            ligand_solute = seekr2_common_sim_sda.Solute()
            ligand_solute.type = "Non_protein"
            ligand_solute.apbs_grid_dime = 129
            ligand_solute.apbs_grid_spacing = base.APBS_GRID_SPACING
            ligand_solute.dielectric = 2.0
            ligand_solute.pqr_filename = self.ligand_pqr_filename_for_bd
            ligand_solute.solute_grid = seekr2_common_sim_sda.Solute_Grid()
            ligand_solute.solute_grid.nb_solute = 1
            model_input.sda_settings_input.solutes = [receptor_solute, ligand_solute]
            
            atomic_parameters1 = seekr2_common_sim_sda.Atomic_parameters()
            atomic_parameters1.type = "H*"
            atomic_parameters1.resname = None
            atomic_parameters1.vdw = 1.20

            atomic_parameters2 = seekr2_common_sim_sda.Atomic_parameters()
            atomic_parameters2.type = "I"
            atomic_parameters2.resname = None
            atomic_parameters2.vdw = 1.98

            atomic_parameters3 = seekr2_common_sim_sda.Atomic_parameters()
            atomic_parameters3.type = "BR"
            atomic_parameters3.resname = None
            atomic_parameters3.vdw = 1.85

            atomic_parameters4 = seekr2_common_sim_sda.Atomic_parameters()
            atomic_parameters4.type = "CL"
            atomic_parameters4.resname = None
            atomic_parameters4.vdw = 1.75

            atomic_parameters5 = seekr2_common_sim_sda.Atomic_parameters()
            atomic_parameters5.type = "F"
            atomic_parameters5.resname = None
            atomic_parameters5.vdw = 1.47

            atomic_parameters6 = seekr2_common_sim_sda.Atomic_parameters()
            atomic_parameters6.type = "N*"
            atomic_parameters6.resname = "NTR"
            atomic_parameters6.test_charge = 1.0

            atomic_parameters7 = seekr2_common_sim_sda.Atomic_parameters()
            atomic_parameters7.type = "O*"
            atomic_parameters7.resname = "CTR"
            atomic_parameters7.test_charge = -0.5

            model_input.sda_settings_input.atoms = [
                atomic_parameters1,
                atomic_parameters2,
                atomic_parameters3,
                atomic_parameters4,
                atomic_parameters5,
                atomic_parameters6,
                atomic_parameters7
            ]

            ion1 = seekr2_base.Ion()
            ion1.radius = 1.2
            ion1.charge = -1.0
            ion1.conc = physical_attributes.ionic_strength
            ion2 = seekr2_base.Ion()
            ion2.radius = 0.9
            ion2.charge = 1.0
            ion2.conc = physical_attributes.ionic_strength
            model_input.sda_settings_input.ions = [ion1, ion2]
            model_input.sda_settings_input.num_b_surface_trajectories \
                = self.bd_settings.num_trajectories
        else:
            model_input.sda_settings_input = None

        return model_input
