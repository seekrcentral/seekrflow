"""
modules/parameters_topology_structures.py

Contain data structure classes used for seekrflow parameters
and topologies.

NOTE: These classes are TEMPORARY and will be replaced with
classes that exist in seekr>=3.
"""

import os
import typing
from shutil import copyfile

from attrs import define, field, validators, Factory

import parmed
import openmm
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
            structure = parmed.load_file(prmtop_full_path, xyz=pdb_filename)
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
            #structure = parmed.load_file(top_full_path)
            #coords_struct = parmed.load_file(pdb_filename, skip_bonds=True)
            #structure.coordinates = coords_struct.coordinates
            structure = parmed.load_file(top_full_path, xyz=pdb_filename)
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
            structure = parmed.load_file(psf_full_path)
            coords_struct = parmed.load_file(pdb_filename, skip_bonds=True)
            structure.coordinates = coords_struct.coordinates
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
            copyfile(os.path.expanduser(forcefield_filename), 
                     anchor_forcefield_full_path)
            new_md_parameters_topology.built_in_forcefield_filenames.append(
                forcefield_basename)

        new_md_parameters_topology.custom_forcefield_filenames = []
        for forcefield_filename in self.custom_forcefield_filenames:
            forcefield_basename = os.path.basename(forcefield_filename)
            anchor_forcefield_full_path = os.path.join(dest_directory,
                                                    forcefield_basename)
            copyfile(os.path.expanduser(forcefield_filename), 
                     anchor_forcefield_full_path)
            new_md_parameters_topology.custom_forcefield_filenames.append(
                forcefield_basename)

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
        full_system_filename = os.path.join(directory, pdb_filename)
        structure = parmed.load_file(full_system_filename)
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
          copyfile(os.path.expanduser(self.system_filename), 
                   anchor_system_full_path)
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
    Amber_parameters_topology, Gromacs_parameters_topology, \
    Charmm_parameters_topology, Forcefield_parameters, Openmm_system]
# ================== END COMMENT ===================