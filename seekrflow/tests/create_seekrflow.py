
"""
tests/create_seekrflow.py

Produce test seekflow objects
"""

import os

import seekrflow.modules.structures as structures

TEST_DIRECTORY = os.path.dirname(__file__)

def create_unparametrized_seekrflow(
        starting_structure: str,
        ligand_resname: str,
        anchor_radius_list: list[float],
        ff: str = "amber",
        ) -> structures.Seekrflow:
    """
    Create a testing seekrflow object before the parametrization step.
    """
    seekrflow = structures.Seekrflow()
    seekrflow.receptor_ligand_pdb = starting_structure
    seekrflow.ligand_resname = ligand_resname
    seekrflow.parametrizer = structures.Parametrizer()
    seekrflow.parametrizer.water_model = "tip3p"
    if ff == "amber":
        seekrflow.parametrizer.forcefield = "gaff-2.11"
        seekrflow.parametrizer.auxiliary_forcefield_files = [
            "amber/ff14SB.xml",
            "amber/tip3p_standard.xml",
            "amber/tip3p_HFE_multivalent.xml"
        ]
    else:
        seekrflow.parametrizer.forcefield = ff
        seekrflow.parametrizer.auxiliary_forcefield_files = [
            "amber/ff14SB.xml",
            "amber/tip3p_standard.xml",
            "amber/tip3p_HFE_multivalent.xml"
        ]
    seekrflow.seekr_settings = structures.MMVT_seekr_settings()
    seekrflow.seekr_settings.anchor_radius_list = anchor_radius_list
    seekrflow.hidr_settings = structures.HIDR_settings_metaD()
    seekrflow.run_settings = structures.Run_settings()
    seekrflow.run_settings.bd_stage_resource_name = "local"
    seekrflow.run_settings.hidr_stage_resource_name = "local"
    seekrflow.run_settings.seekr_stage_resource_name = "local"
    seekrflow.run_settings.allow_parsl_usage_tracking = False
    return seekrflow

def create_parametrized_seekrflow(
        md_parameters_topology: structures.Parameters_topology,
        starting_structure: str,
        receptor_pqr_filename: str,
        ligand_pqr_filename: str,
        ligand_resname: str,
        anchor_radius_list: list[float],
        ff: str = "amber",
        ) -> structures.Seekrflow:
    """
    Create a parametrized seekrflow object right before the prepare stage.
    """
    seekrflow = create_unparametrized_seekrflow(
        starting_structure=starting_structure,
        ligand_resname=ligand_resname,
        anchor_radius_list=anchor_radius_list,
        ff=ff
    )
    seekrflow.md_parameters_topology = md_parameters_topology
    seekrflow.starting_pdb_filename = starting_structure
    seekrflow.bd_settings.receptor_pqr_filename = receptor_pqr_filename
    seekrflow.bd_settings.ligand_pqr_filename = ligand_pqr_filename
    return seekrflow

def create_unparametrized_tryp_ben_seekrflow(
        ff: str = "amber",
        ) -> structures.Seekrflow:
    """
    Create a seekrflow object for the Tryp-Ben system without ff parameters.
    """
    starting_structure = os.path.join(TEST_DIRECTORY, "data", "trypsin_benzamidine.pdb")
    ligand_resname = "BEN"
    anchor_radius_list = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
    return create_unparametrized_seekrflow(
        starting_structure=starting_structure,
        ligand_resname=ligand_resname,
        anchor_radius_list=anchor_radius_list,
        ff=ff
    )

def create_parametrized_tryp_ben_openmm_xml_seekrflow(
        ff: str = "amber",
        ) -> structures.Seekrflow:
    """
    Create a seekrflow object for the Tryp-Ben system with ff parameters.
    """
    starting_pdb = os.path.join(TEST_DIRECTORY, "data", "tryp_ben_system_bound.pdb")
    starting_system_xml = os.path.join(TEST_DIRECTORY, "data", "tryp_ben_system.xml")
    receptor_pqr_filename = os.path.join(TEST_DIRECTORY, "data", "tryp_ben_receptor.pqr")
    ligand_pqr_filename = os.path.join(TEST_DIRECTORY, "data", "tryp_ben_ligand.pqr")
    md_parameters_topology = structures.Openmm_system()
    md_parameters_topology.system_filename = starting_system_xml
    ligand_resname = "BEN"
    anchor_radius_list = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
    return create_parametrized_seekrflow(
        md_parameters_topology=md_parameters_topology,
        starting_structure=starting_pdb,
        receptor_pqr_filename=receptor_pqr_filename,
        ligand_pqr_filename=ligand_pqr_filename,
        ligand_resname=ligand_resname,
        anchor_radius_list=anchor_radius_list,
        ff=ff
    )

def create_parametrized_host_guest_amber_seekrflow(
        ff: str = "amber",
        ) -> structures.Seekrflow:
    """
    Create a seekrflow object for the Tryp-Ben system with ff parameters.
    """
    starting_pdb = os.path.join(TEST_DIRECTORY, "data", "hostguest_at0.5.pdb")
    starting_system_parm7 = os.path.join(TEST_DIRECTORY, "data", "hostguest.parm7")
    receptor_pqr_filename = os.path.join(TEST_DIRECTORY, "data", "hostguest_receptor.pqr")
    ligand_pqr_filename = os.path.join(TEST_DIRECTORY, "data", "hostguest_ligand.pqr")
    md_parameters_topology = structures.Amber_parameters_topology()
    md_parameters_topology.prmtop_filename = starting_system_parm7
    ligand_resname = "BEN"
    anchor_radius_list = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
    return create_parametrized_seekrflow(
        md_parameters_topology=md_parameters_topology,
        starting_structure=starting_pdb,
        receptor_pqr_filename=receptor_pqr_filename,
        ligand_pqr_filename=ligand_pqr_filename,
        ligand_resname=ligand_resname,
        anchor_radius_list=anchor_radius_list,
        ff=ff
    )
