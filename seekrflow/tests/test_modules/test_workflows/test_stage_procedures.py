"""
Tests for seekrflow/modules/workflows/stage_procedures.py
"""

import pytest

import seekrflow.modules.workflows.stage_procedures as stage_procedures
from seekrflow.tests.test_modules.test_workflows.test_membrane_utils import (
    write_popc_fragment,
)


class _FakeCtx:
    def __init__(self, root, md_pdb):
        self.root_directory = str(root)
        self.md_structure_filename = str(md_pdb)
        self.target_temperature = 300.0
        self.selections = {"molecular_dynamics": {}}


def _patch_protein_mapping(monkeypatch, n_atoms):
    indices = list(range(n_atoms))

    def _fake_map(ref_mdtraj, md_mdtraj, align_selection_str=""):
        return indices, indices

    monkeypatch.setattr(
        stage_procedures.utilities_module,
        "obtain_md_ref_structure_map",
        _fake_map,
    )
    monkeypatch.setattr(
        stage_procedures.utilities_module,
        "select_both_from_mdtraj_and_indices",
        lambda md_t, md_i, ref_t, ref_i, sel: ([], []),
    )


def _lipid_row(**overrides):
    row = {
        "name": "equil_custom",
        "run_minimization": False,
        "backbone_force_constant": 0.0,
        "sidechain_force_constant": 0.0,
        "nonprotein_force_constant": 0.0,
        "phosphorus_atoms_force_constant": 12.5,
        "cis_double_bonds_force_constant": 3.5,
        "glycerol_impropers_force_constant": 4.5,
        "number_of_steps": 10,
        "heat_ramp_start": 100.0,
        "heat_ramp_stride": 100.0,
        "heat_ramp_step_interval": 10,
        "timestep": 0.001,
        "ensemble": "NPT_membrane",
    }
    row.update(overrides)
    return row


def _expand(tmp_path, monkeypatch, row):
    pdb = tmp_path / "lipid.pdb"
    traj = write_popc_fragment(pdb)
    monkeypatch.setattr(stage_procedures.mdtraj, "load", lambda filename: traj)
    _patch_protein_mapping(monkeypatch, traj.n_atoms)
    proc = stage_procedures.Equilibration_membrane_stage_procedure(
        name="membrane_equil",
        align_selection_str="not element H",
        schedule=[row],
    )
    resolved = proc.expand(_FakeCtx(tmp_path, pdb))
    return resolved[0][1]


def test_membrane_expand_uses_schedule_row_for_lipid_restraints(
        tmp_path, monkeypatch):
    item = _expand(tmp_path, monkeypatch, _lipid_row())
    assert item.ensemble == "NPT_membrane"
    assert {restraint.type for restraint in item.restraints} == {"z", "torsion"}
    z_restraints = [
        restraint for restraint in item.restraints if restraint.type == "z"]
    assert len(z_restraints) == 1
    assert z_restraints[0].force_constant == 12.5
    torsions = {
        restraint.force_constant: restraint.torsion_values
        for restraint in item.restraints if restraint.type == "torsion"}
    assert torsions[3.5] == [0.0]
    assert torsions[4.5]
    temperatures = item.sampling.heat_ramp_temperatures
    assert temperatures[0] == 100.0
    assert temperatures[-1] == 300.0
    assert item.completion.number_of_steps == 10 * len(temperatures)


def test_membrane_expand_skips_zero_lipid_force_constants(tmp_path, monkeypatch):
    item = _expand(tmp_path, monkeypatch, _lipid_row(
        phosphorus_atoms_force_constant=0.0,
        cis_double_bonds_force_constant=0.0,
        glycerol_impropers_force_constant=0.0,
        heat_ramp_start=None,
    ))
    assert all(
        restraint.type not in ("z", "torsion") for restraint in item.restraints)


def host_guest_procedure():
    equil = stage_procedures.Explicit_stage_procedure(
        name="equilibration",
        items=[stage_procedures.MD_stage_item(name="equilibration")])
    seeding = stage_procedures.Seeding_stage_procedure(
        name="metadynamics", cv_names=["distance"])
    mmvt = stage_procedures.MMVT_stage_procedure(name="MMVT")
    bd = stage_procedures.BD_stage_procedure(name="bd")
    return stage_procedures.Composite_stage_procedure(
        name="full_procedure",
        procedures=[equil, seeding, mmvt, bd])


def test_host_guest_stage_addresses():
    address_map = stage_procedures.build_stage_address_map(host_guest_procedure())
    assert address_map["equilibration"] == (
        ["equilibration", "equilibration"], "equilibration")
    assert address_map["metadynamics_seed"] == (
        ["metadynamics", "sampling"], "sampling")
    assert address_map["metadynamics_assign_structures"] == (
        ["metadynamics", "logistic"], "logistic")
    assert address_map["MMVT"] == (["MMVT"], None)
    assert address_map["bd"] == (["bd"], None)


def test_nested_composite_stage_address():
    grandchild = stage_procedures.MMVT_stage_procedure(name="mmvt_inner")
    child = stage_procedures.Composite_stage_procedure(
        name="child", procedures=[grandchild])
    root = stage_procedures.Composite_stage_procedure(
        name="root", procedures=[child])
    address_map = stage_procedures.build_stage_address_map(root)
    assert address_map["mmvt_inner"] == (["child", "mmvt_inner"], None)


def test_duplicate_stage_names_raise():
    proc = stage_procedures.Composite_stage_procedure(
        name="root",
        procedures=[
            stage_procedures.BD_stage_procedure(name="bd"),
            stage_procedures.MMVT_stage_procedure(name="bd"),
        ])
    with pytest.raises(ValueError, match="Duplicate stage name"):
        stage_procedures.build_stage_address_map(proc)


def test_mmvt_and_seeding_and_ramd_policy_defaults():
    mmvt = stage_procedures.build_stage_policy_map(
        stage_procedures.MMVT_stage_procedure(name="MMVT"))
    assert mmvt["MMVT"].dispatch.dimensions == ["anchor"]
    assert mmvt["MMVT"].dispatch.group_size == 1

    seeding = stage_procedures.build_stage_policy_map(
        stage_procedures.Seeding_stage_procedure(
            name="metadynamics", cv_names=["distance"]))
    assert seeding["metadynamics_seed"].co_schedule_with is None
    assert seeding["metadynamics_assign_structures"].co_schedule_with == (
        "predecessor")

    ramd = stage_procedures.build_stage_policy_map(
        stage_procedures.RAMD_stage_procedure(name="ramd_proc"))
    assert ramd["ramd_proc_ramd"].dispatch.dimensions == ["swarm"]
    assert ramd["ramd_proc_ramd"].dispatch.group_size == 1


def test_equilibration_co_schedules_later_rows_with_predecessor():
    equil = stage_procedures.Equilibration_globular_stage_procedure(
        name="equilibration",
        schedule=[{"name": "heat"}, {"name": "relax"}],
    )
    policy = stage_procedures.build_stage_policy_map(equil)
    assert policy["heat"].co_schedule_with is None
    assert policy["relax"].co_schedule_with == "predecessor"
