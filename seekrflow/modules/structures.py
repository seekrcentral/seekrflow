"""
modules/structures.py

Contain data structure classes used for seekrflow parameters/inputs.
"""

import os
import json
import glob
import typing
import pathlib
from shutil import copyfile

from attrs import define, field, validators, Factory
import cattrs
from cattrs.strategies import include_subclasses
#import seekr.modules.engines.structures as seekr_engines_structures
import parmed
import openmm
import openmm.app as openmm_app

WORK = "work"
PARAMETRIZE = "parametrize"
ROOT = "root"
RUN = "run"

# ================== BEGIN COMMENT ===================
# TODO: remove these structures once seekr3 becomes the default
# and the old seekr program is no longer needed.
@define
class Amber_parameters_topology:
    type: typing.Literal["Amber"] = "Amber"
    prmtop_filename: str = field(default="",
                                validator=validators.instance_of(str))
    
    def same_parameters(
            self, 
            other: "Parameters_topology"
            ) -> bool:
        assert self.type == other.type, "Type mismatch."
        if self.prmtop_filename == other.prmtop_filename:
            return True
        else:
            return False
        
    def copy_files(self,
                   dest_directory: str,
                   new_md_parameters_topology: "Parameters_topology") -> None:
        """
        Copy the files for the initial 'stage' from this Anchor_input
        to the actual anchor.
        """
        assert self.type == new_md_parameters_topology.type, "Type mismatch."
        assert self.prmtop_filename != "", "prmtop_filename must be defined."
        prmtop_basename = os.path.basename(self.prmtop_filename)
        anchor_prmtop_full_path = os.path.join(dest_directory,
                                            prmtop_basename)
        copyfile(os.path.expanduser(self.prmtop_filename), anchor_prmtop_full_path)
        new_md_parameters_topology.prmtop_filename = prmtop_basename
        return
    
    def make_parmed(
            self,
            pdb_filename: None | str = None,
            directory: str = "."
            ) -> parmed.Structure:
        """
        Create a parmed structure from the Amber parameters.
        """
        assert self.prmtop_filename != "", "prmtop_filename must be defined."
        prmtop_full_path = os.path.join(directory, self.prmtop_filename)
        if pdb_filename is None:
            structure = parmed.load_file(prmtop_full_path)
        else:
            structure = parmed.load_file(prmtop_full_path, pdb_filename)
        return structure
    
@define
class Gromacs_parameters_topology:
    type: typing.Literal["Gromacs"] = "Gromacs"
    top_filename: str = field(default="",
                                validator=validators.instance_of(str))
    gro_filename: str = field(default="",
                                validator=validators.instance_of(str))
    
    def same_parameters(
            self, 
            other: "Parameters_topology"
            ) -> bool:
        assert self.type == other.type, "Type mismatch."
        if (self.top_filename == other.top_filename) \
                and (self.gro_filename == other.gro_filename):
            return True
        else:
            return False
    
    def copy_files(self,
                   dest_directory: str,
                   new_md_parameters_topology: "Parameters_topology") -> None:
        """
        Copy the files for the initial 'stage' from this Anchor_input
        to the actual anchor.
        """
        assert self.type == new_md_parameters_topology.type, "Type mismatch."
        assert self.top_filename != "", "top_filename must be defined."
        assert self.gro_filename != "", "gro_filename must be defined."
        top_basename = os.path.basename(self.top_filename)
        anchor_top_full_path = os.path.join(dest_directory,
                                            top_basename)
        copyfile(os.path.expanduser(self.top_filename), anchor_top_full_path)
        new_md_parameters_topology.top_filename = top_basename
        gro_basename = os.path.basename(self.gro_filename)
        anchor_gro_full_path = os.path.join(dest_directory,
                                            gro_basename)
        copyfile(os.path.expanduser(self.gro_filename), anchor_gro_full_path)
        new_md_parameters_topology.gro_filename = gro_basename
        return
    
    def make_parmed(
            self,
            pdb_filename: None | str = None,
            directory: str = "."
            ) -> parmed.Structure:
        """
        Create a parmed structure from the Gromacs parameters.
        """
        assert self.top_filename != "", "top_filename must be defined."
        top_full_path = os.path.join(directory, self.top_filename)
        gro_full_path = os.path.join(directory, self.gro_filename)
        if pdb_filename is None:
            structure = parmed.gromacs.GromacsTopologyFile(top_full_path)
            gmx_gro = parmed.gromacs.GromacsGroFile.parse(gro_full_path)
            structure.box = gmx_gro.box
            structure.positions = gmx_gro.positions
        else:
            structure = parmed.load_file(top_full_path, pdb_filename)
        return structure

@define
class Charmm_parameters_topology:
    type: typing.Literal["Charmm"] = "Charmm"
    psf_filename: str = field(default="",
                                validator=validators.instance_of(str))
    param_filename_list: list[str] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(str),
            iterable_validator=validators.instance_of(list),
        ))
    
    def same_parameters(
            self, 
            other: "Parameters_topology"
            ) -> bool:
        assert self.type == other.type, "Type mismatch."
        if (self.psf_filename == other.psf_filename) \
                and (self.param_filename_list == other.param_filename_list):
            return True
        else:
            return False
        
    def copy_files(self,
                   dest_directory: str,
                   new_md_parameters_topology: "Parameters_topology") -> None:
        """
        Copy the files for the initial 'stage' from this Anchor_input
        to the actual anchor.
        """
        assert self.type == new_md_parameters_topology.type, "Type mismatch."
        assert self.psf_filename != "", "psf_filename must be defined."
        assert len(self.param_filename_list) > 0, "param_filename_list must be defined."
        psf_basename = os.path.basename(self.psf_filename)
        anchor_psf_full_path = os.path.join(dest_directory,
                                            psf_basename)
        copyfile(os.path.expanduser(self.psf_filename), anchor_psf_full_path)
        new_md_parameters_topology.psf_filename = psf_basename
        new_md_parameters_topology.param_filename_list = []
        for param_filename in self.param_filename_list:
            param_basename = os.path.basename(param_filename)
            anchor_param_full_path = os.path.join(dest_directory,
                                                param_basename)
            copyfile(os.path.expanduser(param_filename), anchor_param_full_path)
            new_md_parameters_topology.param_filename_list.append(param_basename)
        return
    
    def make_parmed(
            self,
            pdb_filename: None | str = None,
            directory: str = "."
            ) -> parmed.Structure:
        """
        Create a parmed structure from the Charmm parameters.
        """
        assert self.psf_filename != "", "psf_filename must be defined."
        psf_full_path = os.path.join(directory, self.psf_filename)
        if pdb_filename is None:
            structure = parmed.load_file(psf_full_path)
        else:
            structure = parmed.load_file(psf_full_path, pdb_filename)
        return structure
    
@define
class Forcefield_parameters:
    type: typing.Literal["OpenMM_forcefield"] = "OpenMM_forcefield"
    built_in_forcefield_filenames: list[str] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(str),
            iterable_validator=validators.instance_of(list),
        ))
    custom_forcefield_filenames: list[str] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(str),
            iterable_validator=validators.instance_of(list),
        ))
    
    def same_parameters(
            self, 
            other: "Parameters_topology"
            ) -> bool:
        assert self.type == other.type, "Type mismatch."
        if (self.built_in_forcefield_filenames == other.built_in_forcefield_filenames) \
                and (self.custom_forcefield_filenames == other.custom_forcefield_filenames):
            return True
        else:
            return False
        
    def copy_files(self,
                   dest_directory: str,
                   new_md_parameters_topology: "Parameters_topology") -> None:
        """
        Copy the files for the initial 'stage' from this Anchor_input
        to the actual anchor.
        """
        assert self.type == new_md_parameters_topology.type, "Type mismatch."
        new_md_parameters_topology.built_in_forcefield_filenames = []
        for forcefield_filename in self.built_in_forcefield_filenames:
            forcefield_basename = os.path.basename(forcefield_filename)
            anchor_forcefield_full_path = os.path.join(dest_directory,
                                                    forcefield_basename)
            copyfile(os.path.expanduser(forcefield_filename), anchor_forcefield_full_path)
            new_md_parameters_topology.built_in_forcefield_filenames.append(forcefield_basename)

        new_md_parameters_topology.custom_forcefield_filenames = []
        for forcefield_filename in self.custom_forcefield_filenames:
            forcefield_basename = os.path.basename(forcefield_filename)
            anchor_forcefield_full_path = os.path.join(dest_directory,
                                                    forcefield_basename)
            copyfile(os.path.expanduser(forcefield_filename), anchor_forcefield_full_path)
            new_md_parameters_topology.custom_forcefield_filenames.append(forcefield_basename)

        return
    
    def make_parmed(
            self,
            pdb_filename: None | str = None,
            directory: str = "."
            ) -> parmed.Structure:
        """
        Create a parmed structure from the OpenMM XML parameters.
        """
        assert pdb_filename != "", "pdb_filename must be defined."
        structure = parmed.load_file(pdb_filename)
        return structure
    
@define
class Openmm_system:
    type: typing.Literal["OpenMM_system"] = "OpenMM_system"
    system_filename: str = field(default="",
                                validator=validators.instance_of(str))
    
    def same_parameters(
            self, 
            other: "Parameters_topology"
            ) -> bool:
        assert self.type == other.type, "Type mismatch."
        if self.system_filename == other.system_filename:
            return True
        else:
            return False
        
    def copy_files(self,
                   dest_directory: str,
                   new_md_parameters_topology: "Parameters_topology") -> None:
          """
          Copy the files for the initial 'stage' from this Anchor_input
          to the actual anchor.
          """
          assert self.type == new_md_parameters_topology.type, "Type mismatch."
          assert self.system_filename != "", "system_filename must be defined."
          system_basename = os.path.basename(self.system_filename)
          anchor_system_full_path = os.path.join(dest_directory,
                                              system_basename)
          copyfile(os.path.expanduser(self.system_filename), anchor_system_full_path)
          new_md_parameters_topology.system_filename = system_basename
          return
    
    def make_parmed(
            self,
            pdb_filename: None | str = None,
            directory: str = "."
            ) -> parmed.Structure:
        """
        Create a parmed structure from the OpenMM XML parameters.
        """
        assert pdb_filename != "", "pdb_filename must be defined."
        assert self.system_filename != "", "system_filename must be defined."
        full_system_filename = os.path.join(directory, self.system_filename)
        pdb = openmm_app.PDBFile(pdb_filename)
        with open(full_system_filename) as f:
            system = openmm.XmlSerializer.deserialize(f.read())
        structure = parmed.openmm.load_topology(
            pdb.topology, 
            system, 
            pdb.positions)
        return structure
        

Parameters_topology = typing.Union[
    Amber_parameters_topology, Gromacs_parameters_topology, Charmm_parameters_topology, 
    Forcefield_parameters, Openmm_system]
# ================== END COMMENT ===================

@define
class Parametrizer:
    """
    The parametrizer object contains the inputs needed to run the
    parametrization of a protein-ligand complex.
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
        default=True,
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
        default="AMBER",
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
        cmd = f"pdb2pqr --ff AMBER --ffout AMBER {output_pdb_string}"\
        +f"--with-ph {self.pH} --log-level CRITICAL --drop-water "\
        +f"--nodebump --noopt "\
        +f"{input_pdb_filename} {output_pqr_filename}"
        print("running command:", cmd)
        os.system(cmd)
        assert os.path.exists(output_pqr_filename), \
            f"PDB2PQR output PQR file {output_pqr_filename} was not written. "\
            "A problem must have occurred"
        return

@define
class BD_settings:
    """
    Settings for the BD calculation.
    """
    type: typing.Literal["BD"] = "BD"
    binary_directory: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    receptor_pqr_filename: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
        )
    ligand_pqr_filename: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
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
        default=10000,
        validator=validators.instance_of(int),
        )
    md_steps_per_anchor: int = field(
        default=1000000,
        validator=validators.instance_of(int),
        )
    anchor_radius_list: typing.List[float] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(float),
            iterable_validator=validators.instance_of(list),
        ))

@define
class Transfer_settings_base:
    """
    Base transfer settings class.
    """
    type: typing.Literal["base"] = "base"

@define
class Transfer_settings_globus(Transfer_settings_base):
    """
    Globus transfer settings for transferring files.
    """
    type: typing.Literal["globus"] = "globus"
    local_collection_id: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    remote_collection_id: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    

@define
class Resource_base:
    """
    Base class for resources.
    """
    type: typing.Literal["base"] = "base"
    

@define
class Resource_local(Resource_base):
    """
    Local resource for running the protocol.
    """
    #type: typing.Literal["local"] = "local"
    pass

@define
class Resource_remote_base(Resource_base):
    """
    Base class for remote resources.
    """
    """    type: typing.Literal["remote"] = "remote"
    name: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    remote_seekr2_directory: str = field(
        default="$HOME/seekr2/seekr2/",
        validator=validators.instance_of(str),
        )
    remote_seekrtools_directory: str = field(
        default="$HOME/seekrtools/seekrtools/",
        validator=validators.instance_of(str),
        )
    remote_working_directory: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    globus_compute_endpoint_id: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    transfer_settings: Transfer_settings_base = field(
        default=Factory(Transfer_settings_base),
        validator=validators.instance_of(Transfer_settings_base),
        )
    """
    pass

@define
class Resource_remote_slurm(Resource_remote_base):
    """
    Slurm resource for running the protocol.
    """
    type: typing.Literal["slurm_remote"] = "slurm_remote"
    name: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    remote_seekr2_directory: str = field(
        default="$HOME/seekr2/seekr2/",
        validator=validators.instance_of(str),
        )
    remote_seekrtools_directory: str = field(
        default="$HOME/seekrtools/seekrtools/",
        validator=validators.instance_of(str),
        )
    remote_working_directory: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    max_workers_per_node: int = field(
        default=1,
        validator=validators.instance_of(int),
        )
    partition: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    account: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    constraint: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
        )
    nodes_per_block: int = field(
        default=1,
        validator=validators.instance_of(int),
        )
    cores_per_node: int = field(
        default=1,
        validator=validators.instance_of(int),
        )
    memory_per_node: int = field(
        default=4,
        validator=validators.instance_of(int),
        )
    init_blocks: int | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(int)),
        )
    max_blocks: int | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(int)),
        )
    time_limit: str = field(
        default="00:30:00",
        validator=validators.instance_of(str),
        )
    scheduler_options: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    worker_init: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    globus_compute_endpoint_id: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    transfer_settings: Transfer_settings_globus = field(
        default=Factory(Transfer_settings_globus),
        validator=validators.instance_of(Transfer_settings_globus),
        )
    
@define
class Run_settings:
    """
    Settings for seekrflow runs.
    resource: the machine that a protocol will run on
    """
    resources: typing.List[Resource_remote_slurm] = field(
        default=Factory(list),)
    bd_stage_resource_name: str = field(
        default="local",
        validator=validators.instance_of(str),
        )
    hidr_stage_resource_name: str = field(
        default="local",
        validator=validators.instance_of(str),
        )
    seekr_stage_resource_name: str = field(
        default="local",
        validator=validators.instance_of(str),
        )
    allow_parsl_usage_tracking: bool = field(
        default=False,
        validator=validators.instance_of(bool),
        )
    
    def get_resource_by_name(
            self,
            resource_name: str
            ) -> Resource_base | None:
        """
        Get a resource by its name.
        """
        for resource in self.resources:
            if resource.name == resource_name:
                return resource
        raise ValueError(
            f"Resource with name '{resource_name}' not found in run_settings.resources.")
    
# TODO: make this a protein-ligand subclass of a superclass of generic seekrflows
@define
class Seekrflow:
    """
    All the inputs and parameters needed for a seekrflow calculation.
    """
    name: str = field(
        default="my_name",
        validator=validators.instance_of(str),
        )
    structure_version: str = field(default="1.0",
                                 validator=validators.instance_of(str))
    workflow_type: str = field(
        default="protein_ligand",
        validator=validators.instance_of(str),
    )
    
    # NOTE: although this will typically not be provided in the input file, I want that
    #  option to be available.
    receptor_ligand_pdb: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    ligand_sdf_file: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    ligand_resname: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    ligand_indices: typing.List[int] = field(
        default=Factory(list),
        validator=validators.instance_of(list),
        )
    receptor_selection: str = field(
        default="protein",
        validator=validators.instance_of(str),
        )
    work_directory: str = field(
        default="work",
        validator=validators.instance_of(str),
        )
    basename_output: str = field(
        default="complex",
        validator=validators.instance_of(str),
        )
    solvent_padding: float | None = field(
        default=0.9,
        validator=validators.instance_of(float),
        )
    ionic_strength: float = field(
        default=0.15,
        validator=validators.instance_of(float),
        )
    hmass: float = field(
        default=1.008,
        validator=validators.instance_of(float),
        )
    nonbonded_cutoff: float | None = field(
        default=0.9,
        validator=validators.instance_of(float),
        )
    temperature: float = field(
        default=300.0,
        validator=validators.instance_of(float),
        )
    friction: float = field(
        default=1.0,
        validator=validators.instance_of(float),
        )
    pressure: float | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(float)),
        )
    barostat_period: int | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(int)),
        )
    stepsize: float = field(
        default=0.002,
        validator=validators.instance_of(float),
        )
    parametrizer: Parametrizer | None = field(
        default=Factory(Parametrizer),
        validator=validators.optional(validators.instance_of(Parametrizer)),
        )
    pdb_fixer_settings: PDBFixer_settings | None = field(
        default=Factory(PDBFixer_settings),
        validator=validators.optional(validators.instance_of(PDBFixer_settings)),
        )
    pdb2pqr_settings: PDB2PQR_settings | None = field(
        default=Factory(PDB2PQR_settings),
        validator=validators.optional(validators.instance_of(PDB2PQR_settings)),
        )
    hidr_settings: HIDR_settings_metaD = field(
        default=Factory(HIDR_settings_metaD),
        validator=validators.instance_of(HIDR_settings_metaD),
        )
    #seekr_settings: Seekr_settings = field(
    seekr_settings: MMVT_seekr_settings = field(
        default=Factory(MMVT_seekr_settings),
        validator=validators.instance_of(MMVT_seekr_settings),
        )
    bd_settings: BD_settings | None = field(
        default=Factory(BD_settings),
        validator=validators.optional(validators.instance_of(BD_settings)),
        )
    # This is filled out by parametrize.py, but can be entered manually if desired.
    md_parameters_topology: Parameters_topology | None = None
    starting_pdb_filename: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    
    run_settings: Run_settings | None = None

    def save(
            self,
            filename: str
            ) -> None:
        """
        Save the Seekrflow object to a JSON file.
        """
        converter: cattrs.Converter = cattrs.Converter()
        # Make sure that interited data classes are unstructured as their
        #  subtypes.
        include_subclasses(HIDR_settings_base, converter)
        include_subclasses(Transfer_settings_base, converter)
        include_subclasses(Resource_base, converter)
        seekrflow_dict: dict = converter.unstructure(self)
        json_dump: str = json.dumps(seekrflow_dict, indent=4)
        with open(filename, "w") as file:
            file.write(json_dump)
        return
    
    def make_work_directory(
            self,
            work_directory: pathlib.Path | str | None = None
            ) -> None:
        """
        Make the work directory for the Seekrflow calculation.
        """
        if work_directory is not None:
            self.work_directory = str(work_directory)
        os.makedirs(self.work_directory, exist_ok=True)
        return

    def get_work_directory(self) -> pathlib.Path:
        """
        Get the directory where the preparation files are stored.
        """
        work_dir = pathlib.Path(self.work_directory)
        #os.makedirs(work_dir, exist_ok=True)
        return work_dir

    def get_parametrize_directory(self) -> pathlib.Path:
        """
        Get the directory where the preparation files are stored.
        """
        param_dir = pathlib.Path(self.work_directory) / PARAMETRIZE
        os.makedirs(param_dir, exist_ok=True)
        return param_dir

    def get_root_directory(self) -> pathlib.Path:
        """
        Get the root directory where the Seekrflow files are stored.
        """
        root_dir = pathlib.Path(self.work_directory) / ROOT
        os.makedirs(root_dir, exist_ok=True)
        return root_dir
    
    def get_run_directory(self) -> pathlib.Path:
        """
        Get the directory where the run files are stored.
        """
        run_dir = pathlib.Path(self.work_directory) / RUN
        os.makedirs(run_dir, exist_ok=True)
        return run_dir

def load_seekrflow(
        filename: str
    ) -> Seekrflow:
    """
    Load a Seekrflow object from a JSON file.
    """
    with open(filename, "r") as file:
        json_string: str = json.load(file)
    converter: cattrs.Converter = cattrs.Converter()
    seekrflow_obj: Seekrflow = converter.structure(json_string, Seekrflow)
    return seekrflow_obj

def save_new_seekrflow(
        seekrflow: Seekrflow, 
        seekrflow_glob: str, 
        seekrflow_base: str, 
        save_old_seekrflow=True,
        directory = "."):
    """
    Generate a new seekrflow file. The old seekrflow file(s) will be renamed with a 
    numerical index.
    
    """
    
    model_path = os.path.join(directory, "seekrflow.json")
    if os.path.exists(model_path) and save_old_seekrflow:
        # This is expected, because this old model was loaded
        full_model_glob = os.path.join(directory, seekrflow_glob)
        num_globs = len(glob.glob(full_model_glob))
        new_pre_model_filename = seekrflow_base.format(num_globs)
        new_pre_model_path = os.path.join(directory, 
                                          new_pre_model_filename)
        print("Renaming model.xml to {}".format(new_pre_model_filename))
        copyfile(model_path, new_pre_model_path)
        
    print("Saving new seekrflow.json")
    seekrflow.save(model_path)
    return
