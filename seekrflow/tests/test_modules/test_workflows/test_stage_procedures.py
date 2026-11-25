"""
Tests for seekrflow/modules/workflows/stage_procedures.py
"""

import numpy as np
import mdtraj
import pytest
from seekr.modules.base import mdtraj_center_of_mass

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


def _atom(name, element):
    return (name, element)


def _trajectory(chains, xyz):
    topology = mdtraj.Topology()
    for chain_id, residues in chains:
        chain = topology.add_chain(chain_id=chain_id)
        for resname, res_seq, atoms in residues:
            residue = topology.add_residue(resname, chain, resSeq=res_seq)
            for atom_name, element in atoms:
                topology.add_atom(atom_name, element, residue)
    return mdtraj.Trajectory(np.asarray(xyz, dtype=np.float32), topology)


def _globular_row(**overrides):
    row = {
        "name": "equil_custom",
        "run_minimization": False,
        "backbone_force_constant": 0.0,
        "sidechain_force_constant": 0.0,
        "nonprotein_force_constant": 0.0,
        "number_of_steps": 10,
        "timestep": 0.001,
        "ensemble": "NPT",
    }
    row.update(overrides)
    return row


def _expand_globular(tmp_path, monkeypatch, md_traj, schedule, com_traj=None):
    md_name = str(tmp_path / "solvated.pdb")
    com_name = str(tmp_path / "com_reference.pdb")

    def _load(filename):
        filename = str(filename)
        if filename == com_name:
            assert com_traj is not None
            return com_traj
        if filename == md_name:
            return md_traj
        raise AssertionError(filename)

    monkeypatch.setattr(stage_procedures.mdtraj, "load", _load)
    _patch_protein_mapping(monkeypatch, md_traj.n_atoms)
    proc = stage_procedures.Equilibration_globular_stage_procedure(
        name="equilibration",
        com_reference_structure_filename=com_name if com_traj is not None else None,
        schedule=schedule,
    )
    ctx = _FakeCtx(tmp_path, md_name)
    resolved = proc.expand(ctx)
    return [item for _, item in resolved], ctx


def _com_restraints(item):
    return [
        restraint for restraint in item.restraints
        if restraint.type == "center_of_mass"]


def test_globular_expand_restrains_every_atom_of_a_com_reference_residue(
        tmp_path, monkeypatch):
    carbon = mdtraj.element.carbon
    md_traj = _trajectory(
        [("A", [
            ("ALA", 1, [
                _atom("N", mdtraj.element.nitrogen),
                _atom("H", mdtraj.element.hydrogen),
                _atom("CA", carbon),
            ]),
            ("ALA", 2, [_atom("CA", carbon)]),
        ])],
        [[[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.2, 0.0, 0.0], [5.0, 0.0, 0.0]]],
    )
    com_traj = _trajectory(
        [("A", [("ALA", 1, [_atom("CX", carbon)])])],
        [[[1.5, -0.25, 0.3]]],
    )
    items, ctx = _expand_globular(
        tmp_path, monkeypatch, md_traj,
        [_globular_row(com_force_constant=12.5)], com_traj)
    restraints = _com_restraints(items[0])
    assert len(restraints) == 1
    assert restraints[0].force_constant == 12.5
    assert restraints[0].selection_names == ["equilibration_com_res_A_1_ALA"]
    selections = ctx.selections["molecular_dynamics"]
    assert selections["equilibration_com_res_A_1_ALA"] == [0, 1, 2]
    assert "equilibration_com_res_A_2_ALA" not in selections
    coords = np.loadtxt(restraints[0].coordinates_filename, ndmin=2)
    np.testing.assert_allclose(coords, [[1.5, -0.25, 0.3]])


def test_globular_expand_writes_mass_weighted_com_of_reference_atoms(
        tmp_path, monkeypatch):
    md_traj = _trajectory(
        [("A", [("ALA", 1, [_atom("CA", mdtraj.element.carbon)])])],
        [[[0.0, 0.0, 0.0]]],
    )
    com_traj = _trajectory(
        [("A", [("ALA", 1, [
            _atom("N", mdtraj.element.nitrogen),
            _atom("C", mdtraj.element.carbon),
        ])])],
        [[[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]]],
    )
    items, _ctx = _expand_globular(
        tmp_path, monkeypatch, md_traj,
        [_globular_row(com_force_constant=1.0)], com_traj)
    coords = np.loadtxt(_com_restraints(items[0])[0].coordinates_filename)
    expected = mdtraj_center_of_mass(com_traj)[0]
    np.testing.assert_allclose(np.atleast_2d(coords), np.atleast_2d(expected))


def test_globular_expand_stacks_one_com_block_per_reference_frame(
        tmp_path, monkeypatch):
    carbon = mdtraj.element.carbon
    md_traj = _trajectory(
        [("A", [
            ("GLY", 5, [
                _atom("N", mdtraj.element.nitrogen),
                _atom("CA", carbon),
            ]),
            ("ALA", 1, [_atom("CA", carbon)]),
        ])],
        [[[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [3.0, 0.0, 0.0]]],
    )
    com_xyz = [
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
    ]
    com_traj = _trajectory(
        [("A", [
            ("ALA", 1, [_atom("CX", carbon)]),
            ("GLY", 5, [_atom("CX", carbon)]),
        ])],
        com_xyz,
    )
    items, ctx = _expand_globular(
        tmp_path, monkeypatch, md_traj,
        [_globular_row(com_force_constant=4.0)], com_traj)
    restraint = _com_restraints(items[0])[0]
    assert restraint.selection_names == [
        "equilibration_com_res_A_1_ALA",
        "equilibration_com_res_A_5_GLY",
    ]
    selections = ctx.selections["molecular_dynamics"]
    assert selections["equilibration_com_res_A_1_ALA"] == [2]
    assert selections["equilibration_com_res_A_5_GLY"] == [0, 1]
    coords = np.loadtxt(restraint.coordinates_filename)
    np.testing.assert_allclose(coords, [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
    ])


def test_globular_expand_omits_com_restraint_when_force_constant_is_zero(
        tmp_path, monkeypatch):
    md_traj = _trajectory(
        [("A", [("ALA", 1, [_atom("CA", mdtraj.element.carbon)])])],
        [[[0.0, 0.0, 0.0]]],
    )
    md_name = str(tmp_path / "solvated.pdb")

    def _load(filename):
        assert str(filename) == md_name
        return md_traj

    monkeypatch.setattr(stage_procedures.mdtraj, "load", _load)
    _patch_protein_mapping(monkeypatch, md_traj.n_atoms)
    proc = stage_procedures.Equilibration_globular_stage_procedure(
        name="equilibration",
        com_reference_structure_filename=str(tmp_path / "missing_com.pdb"),
        schedule=[_globular_row(com_force_constant=0.0)],
    )
    ctx = _FakeCtx(tmp_path, md_name)
    item = proc.expand(ctx)[0][1]
    assert _com_restraints(item) == []
    assert all(
        "com_res" not in name
        for name in ctx.selections["molecular_dynamics"])


def test_globular_expand_requires_a_com_file_when_force_constant_is_positive(
        tmp_path):
    proc = stage_procedures.Equilibration_globular_stage_procedure(
        name="equilibration",
        schedule=[_globular_row(com_force_constant=5.0)],
    )
    with pytest.raises(ValueError, match="com_reference_structure_filename"):
        proc.expand(_FakeCtx(tmp_path, tmp_path / "solvated.pdb"))


def test_globular_expand_rejects_a_com_residue_missing_from_the_solvated_structure(
        tmp_path, monkeypatch):
    carbon = mdtraj.element.carbon
    md_traj = _trajectory(
        [("A", [("ALA", 1, [_atom("CA", carbon)])])],
        [[[0.0, 0.0, 0.0]]],
    )
    com_traj = _trajectory(
        [("A", [("ALA", 99, [_atom("CX", carbon)])])],
        [[[0.0, 0.0, 0.0]]],
    )
    with pytest.raises(ValueError, match="no solvated residue"):
        _expand_globular(
            tmp_path, monkeypatch, md_traj,
            [_globular_row(com_force_constant=5.0)], com_traj)


def test_membrane_expand_adds_com_restraint_from_the_reference_file(
        tmp_path, monkeypatch):
    carbon = mdtraj.element.carbon
    md_traj = _trajectory(
        [("A", [("ALA", 1, [_atom("CA", carbon), _atom("H", mdtraj.element.hydrogen)])])],
        [[[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]]],
    )
    com_traj = _trajectory(
        [("A", [("ALA", 1, [_atom("CX", carbon)])])],
        [[[0.4, 0.5, 0.6]]],
    )
    md_name = str(tmp_path / "solvated.pdb")
    com_name = str(tmp_path / "com_reference.pdb")

    def _load(filename):
        filename = str(filename)
        if filename == com_name:
            return com_traj
        if filename == md_name:
            return md_traj
        raise AssertionError(filename)

    monkeypatch.setattr(stage_procedures.mdtraj, "load", _load)
    _patch_protein_mapping(monkeypatch, md_traj.n_atoms)
    proc = stage_procedures.Equilibration_membrane_stage_procedure(
        name="membrane_equil",
        align_selection_str="not element H",
        com_reference_structure_filename=com_name,
        schedule=[_lipid_row(
            phosphorus_atoms_force_constant=0.0,
            cis_double_bonds_force_constant=0.0,
            glycerol_impropers_force_constant=0.0,
            com_force_constant=8.0,
            heat_ramp_start=None,
        )],
    )
    ctx = _FakeCtx(tmp_path, md_name)
    item = proc.expand(ctx)[0][1]
    restraints = _com_restraints(item)
    assert len(restraints) == 1
    assert restraints[0].force_constant == 8.0
    assert restraints[0].selection_names == ["membrane_equil_com_res_A_1_ALA"]
    assert ctx.selections["molecular_dynamics"][
        "membrane_equil_com_res_A_1_ALA"] == [0, 1]


def test_equilibration_co_schedules_later_rows_with_predecessor():
    equil = stage_procedures.Equilibration_globular_stage_procedure(
        name="equilibration",
        schedule=[{"name": "heat"}, {"name": "relax"}],
    )
    policy = stage_procedures.build_stage_policy_map(equil)
    assert policy["heat"].co_schedule_with is None
    assert policy["relax"].co_schedule_with == "predecessor"
