
"""
tests/create_seekrflow.py

Produce test seekflow objects
"""

import os

import seekrflow.modules.structures as structures
import seekrflow.modules.workflows.structures as workflow_structures
import seekrflow.modules.parameterize_structures as parameterizer_structures
import seekrflow.modules.parameters_topology_structures \
    as parameters_topology_structures
import seekrflow.modules.workflows.protein_ligand_seekr2.structures as protein_ligand_seekr2_structures

TEST_DIRECTORY = os.path.dirname(__file__)

def create_unparameterized_seekrflow(
        starting_structure: str,
        ligand_resname: str,
        anchor_radius_list: list[float],
        ff: str = "amber",
        ) -> structures.Seekrflow:
    """
    Create a testing seekrflow object before the parameterization step.
    """
    seekrflow = structures.Seekrflow()
    seekrflow.name = "test_seekrflow"
    workflow = protein_ligand_seekr2_structures.Protein_ligand_seekr2_workflow()
    workflow.parameterizer_information \
        = protein_ligand_seekr2_structures.Parameterizer_information()
    workflow.parameterizer_information.ligand_resname = ligand_resname
    workflow.parameterizer_information.receptor_ligand_pdb_filename = starting_structure
    workflow.hidr_settings = protein_ligand_seekr2_structures.HIDR_settings_metaD()
    workflow.mmvt_settings = protein_ligand_seekr2_structures.MMVT_seekr_settings()
    workflow.mmvt_settings.anchor_radius_list = anchor_radius_list
    workflow.md_settings = workflow_structures.MD_settings()
    seekrflow.workflow = workflow
    parameterizer = parameterizer_structures.Parameterizer()
    parameterizer.water_model = "tip3p"
    if ff == "amber":
        parameterizer.forcefield = "gaff-2.11"
        parameterizer.auxiliary_forcefield_files = [
            "amber/ff14SB.xml",
            "amber/tip3p_standard.xml",
            "amber/tip3p_HFE_multivalent.xml"
        ]
    else:
        parameterizer.forcefield = ff
        parameterizer.auxiliary_forcefield_files = [
            "amber/ff14SB.xml",
            "amber/tip3p_standard.xml",
            "amber/tip3p_HFE_multivalent.xml"
        ]
    seekrflow.parameterizer = parameterizer
    seekrflow.run_settings = structures.Run_settings()
    seekrflow.run_settings.placements = []
    return seekrflow

def create_parameterized_seekrflow(
        solvated_system_for_md: workflow_structures.Solvated_system_for_md,
        receptor_pqr_filename: str,
        ligand_pqr_filename: str,
        ligand_resname: str,
        anchor_radius_list: list[float],
        ff: str = "amber",
        ) -> structures.Seekrflow:
    """
    Create a parameterized seekrflow object right before the prepare stage.
    """
    seekrflow = create_unparameterized_seekrflow(
        starting_structure=solvated_system_for_md.solvated_pdb,
        ligand_resname=ligand_resname,
        anchor_radius_list=anchor_radius_list,
        ff=ff
    )
    seekrflow.workflow.solvated_system_for_md = solvated_system_for_md
    seekrflow.workflow.receptor_pqr_filename_for_bd = receptor_pqr_filename
    seekrflow.workflow.ligand_pqr_filename_for_bd = ligand_pqr_filename
    return seekrflow

def create_unparameterized_tryp_ben_seekrflow(
        ff: str = "amber",
        ) -> structures.Seekrflow:
    """
    Create a seekrflow object for the Tryp-Ben system without ff parameters.
    """
    starting_structure = os.path.join(TEST_DIRECTORY, "data", "trypsin_benzamidine.pdb")
    ligand_resname = "BEN"
    anchor_radius_list = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
    return create_unparameterized_seekrflow(
        starting_structure=starting_structure,
        ligand_resname=ligand_resname,
        anchor_radius_list=anchor_radius_list,
        ff=ff
    )

def create_parameterized_tryp_ben_openmm_xml_seekrflow(
        ff: str = "amber",
        ) -> structures.Seekrflow:
    """
    Create a seekrflow object for the Tryp-Ben system with ff parameters.
    """
    starting_pdb = os.path.join(TEST_DIRECTORY, "data", "tryp_ben_system_bound.pdb")
    starting_system_xml = os.path.join(TEST_DIRECTORY, "data", "tryp_ben_system.xml")
    receptor_pqr_filename = os.path.join(TEST_DIRECTORY, "data", "tryp_ben_receptor.pqr")
    ligand_pqr_filename = os.path.join(TEST_DIRECTORY, "data", "tryp_ben_ligand.pqr")
    md_parameters_topology = parameters_topology_structures.Openmm_system()
    md_parameters_topology.system_filename = starting_system_xml
    ligand_resname = "BEN"
    anchor_radius_list = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
    solvated_system_for_md = workflow_structures.Solvated_system_for_md()
    solvated_system_for_md.solvated_pdb = starting_pdb
    solvated_system_for_md.parameters_topology = md_parameters_topology
    return create_parameterized_seekrflow(
        solvated_system_for_md=solvated_system_for_md,
        receptor_pqr_filename=receptor_pqr_filename,
        ligand_pqr_filename=ligand_pqr_filename,
        ligand_resname=ligand_resname,
        anchor_radius_list=anchor_radius_list,
        ff=ff
    )

def create_parameterized_host_guest_amber_seekrflow(
        ff: str = "amber",
        ) -> structures.Seekrflow:
    """
    Create a seekrflow object for the Tryp-Ben system with ff parameters.
    """
    starting_pdb = os.path.join(TEST_DIRECTORY, "data", "hostguest_at0.5.pdb")
    starting_system_parm7 = os.path.join(TEST_DIRECTORY, "data", "hostguest.parm7")
    receptor_pqr_filename = os.path.join(TEST_DIRECTORY, "data", "hostguest_receptor.pqr")
    ligand_pqr_filename = os.path.join(TEST_DIRECTORY, "data", "hostguest_ligand.pqr")
    md_parameters_topology = parameters_topology_structures.Amber_parameters_topology()
    md_parameters_topology.prmtop_filename = starting_system_parm7
    ligand_resname = "APN"
    anchor_radius_list = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
    solvated_system_for_md = workflow_structures.Solvated_system_for_md()
    solvated_system_for_md.solvated_pdb = starting_pdb
    solvated_system_for_md.parameters_topology = md_parameters_topology
    flow = create_parameterized_seekrflow(
        solvated_system_for_md=solvated_system_for_md,
        receptor_pqr_filename=receptor_pqr_filename,
        ligand_pqr_filename=ligand_pqr_filename,
        ligand_resname=ligand_resname,
        anchor_radius_list=anchor_radius_list,
        ff=ff
    )
    flow.workflow.receptor_indices = list(range(147))
    flow.workflow.ligand_indices = list(range(147,162))
    return flow
