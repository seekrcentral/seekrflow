"""
Tests for seekrflow/tests/create_seekrflow.py
"""

import os

import mdtraj

import seekrflow.modules.parameters_topology_structures as pts
import seekrflow.modules.workflows.components as components_module
import seekrflow.modules.workflows.stage_procedures as stage_procedures
import seekrflow.tests.create_seekrflow as create_seekrflow


def _data(*names):
    return os.path.join(create_seekrflow.TEST_DIRECTORY, "data", *names)


def _seeding(flow):
    return next(
        procedure for procedure in flow.workflow.procedure.procedures
        if isinstance(procedure, stage_procedures.Seeding_stage_procedure))


def test_unparameterized_tryp_ben_wires_parameterizer_and_cv():
    flow = create_seekrflow.create_unparameterized_tryp_ben_seekrflow()
    assert flow.parameterizer is not None
    assert flow.workflow.get_md_settings().system is None
    binding_site = next(
        selector for selector in flow.workflow.components.group_selectors
        if isinstance(
            selector, components_module.Binding_site_group_selector))
    assert flow.workflow.cv_specs[0].group1_selection_name == binding_site.name
    member_names = {member.name for member in flow.workflow.components.members}
    bd = flow.workflow.get_scale_settings("brownian_dynamics")
    assert {
        molecule.component_name for molecule in bd.system.molecules
    } == member_names
    assert _seeding(flow).cv_names[0] == flow.workflow.cv_specs[0].name


def test_parameterized_tryp_ben_points_at_openmm_files():
    flow = create_seekrflow.create_parameterized_tryp_ben_openmm_xml_seekrflow()
    assert flow.parameterizer is None
    system = flow.workflow.get_md_settings().system
    assert isinstance(system.parameters_topology, pts.Openmm_system)
    solvated = _data("tryp_ben_system_bound.pdb")
    xml = _data("tryp_ben_system.xml")
    assert system.solvated_pdb == solvated
    assert os.path.isfile(solvated)
    assert system.parameters_topology.system_filename == xml
    assert os.path.isfile(xml)


def test_host_guest_uses_small_molecules_and_explicit_equilibration():
    flow = create_seekrflow.create_parameterized_host_guest_amber_seekrflow()
    members = flow.workflow.components.members
    assert all(
        isinstance(member, components_module.Small_molecule_component)
        for member in members)
    assert flow.workflow.components.group_selectors == []
    assert flow.workflow.cv_specs[0].group1_selection_name == "receptor"
    assert isinstance(
        flow.workflow.procedure.procedures[0],
        stage_procedures.Explicit_stage_procedure)
    assert flow.physical_attributes.ionic_strength == 0.0
    member_names = {member.name for member in members}
    bd = flow.workflow.get_scale_settings("brownian_dynamics")
    assert {
        molecule.component_name for molecule in bd.system.molecules
    } == member_names
    traj = mdtraj.load(_data("hostguest_at0.5.pdb"))
    host = flow.workflow.components.get_member("receptor").resolve(traj)
    guest = flow.workflow.components.get_member("ligand").resolve(traj)
    assert host
    assert guest
    assert set(host).isdisjoint(guest)
