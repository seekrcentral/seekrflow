"""
modules/structures.py

Contain data structure classes used for seekrflow parameters/inputs.
"""

import os
import json
import typing
from shutil import copyfile

from attrs import define, field, validators, Factory
import cattrs
#import seekr.modules.engines.structures as seekr_engines_structures
import openmm.app as openmm_app

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
        default="gaff-2.2.20",
        validator=validators.instance_of(str),
        # TODO: possible choices?
        )
    water_model: str = field(
        default="tip3p",
        validator=validators.instance_of(str),
        # TODO: possible choices?
        )
    intermediate_forcefield_files: typing.List[str] = field(
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
        +f"{input_pdb_filename} {output_pqr_filename}"
        print("running command:", cmd)
        os.system(cmd)
        assert os.path.exists(output_pqr_filename), \
            f"PDB2PQR output PQR file {output_pqr_filename} was not written. "\
            "A problem must have occurred"
        return


@define
class Seekrflow:
    """
    All the inputs and parameters needed for a seekrflow calculation.
    """
    structure_version: str = field(default="1.0",
                                 validator=validators.instance_of(str))
    workflow_type: str = field(
        default="ligand_protein",
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
    work_directory: str = field(
        default=".",
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
    temperature: float = field(
        default=300.0,
        validator=validators.instance_of(float),
        )
    friction: float = field(
        default=1.0,
        validator=validators.instance_of(float),
        )
    pressure: float | None = field(
        default=1.0,
        validator=validators.instance_of(float),
        )
    barostat_period: int = field(
        default=50,
        validator=validators.instance_of(int),
        )
    stepsize: float = field(
        default=0.002,
        validator=validators.instance_of(float),
        )
    parametrizer: Parametrizer = field(
        default=Factory(Parametrizer),
        validator=validators.instance_of(Parametrizer),
        )
    pdb_fixer_settings: PDBFixer_settings = field(
        default=Factory(PDBFixer_settings),
        validator=validators.instance_of(PDBFixer_settings),
        )
    pdb2pqr_settings: PDB2PQR_settings = field(
        default=Factory(PDB2PQR_settings),
        validator=validators.instance_of(PDB2PQR_settings),
        )
    # This is filled out by parametrize.py, but can be entered manually if desired.
    md_parameters_topology: Parameters_topology | None = None
    
    def save(
            self,
            filename: str
            ) -> None:
        """
        Save the Seekrflow object to a JSON file.
        """
        converter: cattrs.Converter = cattrs.Converter()
        seekrflow_dict: dict = converter.unstructure(self)
        json_dump: str = json.dumps(seekrflow_dict, indent=4)
        with open(filename, "w") as file:
            file.write(json_dump)
        return
    
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