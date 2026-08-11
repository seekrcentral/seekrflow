"""
Tests for membrane equilibration lipid detection and stage expansion.
"""

import os
import tempfile
import textwrap

import numpy as np
import pytest
import mdtraj

import seekrflow.modules.workflows.membrane_utils as membrane_utils
import seekrflow.modules.workflows.stage_procedures as stage_procedures


def _write_tiny_lipid_pdb(path: str) -> None:
    """
    Minimal CHARMM-like POPC fragment: glycerol C1/C2/C3/O21, P, and a
    cis C=C (C22=C23) with flanking carbons for dihedral detection.
    """
    # Coordinates chosen so C1-C3-C2-O21 improper is near -120 deg.
    content = textwrap.dedent("""
    CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1
    ATOM      1  C1  POPC A   1      0.000   0.000   0.000  1.00  0.00           C
    ATOM      2  C2  POPC A   1      1.500   0.000   0.000  1.00  0.00           C
    ATOM      3  C3  POPC A   1      2.000   1.400   0.000  1.00  0.00           C
    ATOM      4  O21 POPC A   1      2.000  -0.700  -1.200  1.00  0.00           O
    ATOM      5  P   POPC A   1     -1.500   0.000   0.000  1.00  0.00           P
    ATOM      6  C21 POPC A   1      3.200  -0.700  -1.200  1.00  0.00           C
    ATOM      7  C22 POPC A   1      4.500  -0.700  -1.200  1.00  0.00           C
    ATOM      8  C23 POPC A   1      5.700  -0.700  -1.200  1.00  0.00           C
    ATOM      9  C24 POPC A   1      7.000  -0.700  -1.200  1.00  0.00           C
    ATOM     10  H22 POPC A   1      4.500   0.300  -1.200  1.00  0.00           H
    ATOM     11  H23 POPC A   1      5.700   0.300  -1.200  1.00  0.00           H
    CONECT    1    2
    CONECT    2    1    3    4
    CONECT    3    2
    CONECT    4    2    6
    CONECT    5    1
    CONECT    6    4    7
    CONECT    7    6    8   10
    CONECT    8    7    9   11
    CONECT    9    8
    CONECT   10    7
    CONECT   11    8
    END
    """).strip() + "\n"
    with open(path, "w") as f:
        f.write(content)


def test_default_membrane_schedule_has_eight_stages_no_headgroup():
    assert len(stage_procedures.DEFAULT_MEMBRANE_EQUIL_SCHEDULE) == 8
    for row in stage_procedures.DEFAULT_MEMBRANE_EQUIL_SCHEDULE:
        assert "lipid_headgroup_force_constant" not in row
        assert row["ensemble"] == "NPT_membrane"
    # Doubled torsion defaults vs legacy draft (1046 -> 2092 at stage 1)
    assert stage_procedures.DEFAULT_MEMBRANE_EQUIL_SCHEDULE[0][
        "cis_double_bonds_force_constant"] == 2092.0
    assert stage_procedures.DEFAULT_MEMBRANE_EQUIL_SCHEDULE[0][
        "glycerol_impropers_force_constant"] == 2092.0


def test_equilibration_membrane_in_stage_procedure_union():
    assert stage_procedures.Equilibration_membrane_stage_procedure in \
        stage_procedures.Stage_procedure.__args__


def test_glycerol_improper_theta0_prefers_nearest_120():
    # Synthetic positions for (C1, C3, C2, O21) near -120 deg
    positions = np.zeros((4, 3))
    positions[0] = [0.0, 0.0, 0.0]       # C1
    positions[1] = [1.0, 1.0, 0.0]       # C3
    positions[2] = [1.0, 0.0, 0.0]       # C2
    positions[3] = [1.0, -0.5, -0.866]   # O21 ~ -120-ish
    theta0 = membrane_utils.glycerol_improper_theta0(
        positions, (0, 1, 2, 3))
    assert abs(theta0) == pytest.approx(
        membrane_utils.GLYCEROL_IMPROPER_THETA0_POS, abs=1e-9) \
        or abs(theta0) == pytest.approx(
            abs(membrane_utils.GLYCEROL_IMPROPER_THETA0_NEG), abs=1e-9)
    assert abs(abs(theta0) - abs(
        membrane_utils.GLYCEROL_IMPROPER_THETA0_NEG)) < 1e-9


def test_identify_lipid_atoms_finds_phosphorus_and_glycerol():
    with tempfile.TemporaryDirectory() as tmpdir:
        pdb = os.path.join(tmpdir, "lipid.pdb")
        _write_tiny_lipid_pdb(pdb)
        traj = mdtraj.load(pdb)
        # mdtraj may not keep CONECT as bonds for all formats; ensure bonds
        if traj.topology.n_bonds == 0:
            pytest.skip("PDB CONECT bonds not loaded by mdtraj in this env")
        result = membrane_utils.identify_lipid_atoms(traj)
        assert len(result["phosphorus"]) >= 1
        assert len(result["glycerol_impropers"]) >= 1
        improper = result["glycerol_impropers"][0]
        assert len(improper) == 4


class _FakeCtx:
    def __init__(self, root, md_pdb):
        self.root_directory = root
        self.md_structure_filename = md_pdb
        self.target_temperature = 300.0
        self.selections = {"molecular_dynamics": {}}


def test_membrane_expand_emits_z_and_torsion_restraints(monkeypatch):
    with tempfile.TemporaryDirectory() as tmpdir:
        pdb = os.path.join(tmpdir, "lipid.pdb")
        _write_tiny_lipid_pdb(pdb)
        traj = mdtraj.load(pdb)
        if traj.topology.n_bonds == 0:
            pytest.skip("PDB CONECT bonds not loaded by mdtraj in this env")

        n = traj.n_atoms
        indices = list(range(n))

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

        proc = stage_procedures.Equilibration_membrane_stage_procedure(
            name="membrane_equil",
            align_selection_str="not element H",
            schedule=[dict(row) for row in
                      stage_procedures.DEFAULT_MEMBRANE_EQUIL_SCHEDULE[:1]],
        )
        ctx = _FakeCtx(tmpdir, pdb)
        resolved = proc.expand(ctx, input_stage_name=stage_procedures.INITIAL)
        assert len(resolved) == 1
        item = resolved[0][1]
        restraint_types = {r.type for r in item.restraints}
        assert "z" in restraint_types
        assert "torsion" in restraint_types
        assert item.ensemble == "NPT_membrane"
        z_specs = [r for r in item.restraints if r.type == "z"]
        assert z_specs[0].force_constant == 2092.0
        torsion_specs = [r for r in item.restraints if r.type == "torsion"]
        assert len(torsion_specs) >= 1
        flat_vals = [v for spec in torsion_specs for v in spec.torsion_values]
        assert any(abs(v) < 1e-12 for v in flat_vals) or any(
            abs(abs(v) - abs(membrane_utils.GLYCEROL_IMPROPER_THETA0_NEG)) < 1e-5
            for v in flat_vals)


def test_membrane_expand_skips_zero_force_constants(monkeypatch):
    with tempfile.TemporaryDirectory() as tmpdir:
        pdb = os.path.join(tmpdir, "lipid.pdb")
        _write_tiny_lipid_pdb(pdb)
        traj = mdtraj.load(pdb)
        if traj.topology.n_bonds == 0:
            pytest.skip("PDB CONECT bonds not loaded by mdtraj in this env")
        n = traj.n_atoms
        indices = list(range(n))

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
        row = dict(stage_procedures.DEFAULT_MEMBRANE_EQUIL_SCHEDULE[-1])
        assert row["phosphorus_atoms_force_constant"] == 0.0
        proc = stage_procedures.Equilibration_membrane_stage_procedure(
            name="membrane_equil",
            align_selection_str="not element H",
            schedule=[row],
        )
        ctx = _FakeCtx(tmpdir, pdb)
        resolved = proc.expand(ctx)
        item = resolved[0][1]
        assert all(r.type not in ("z", "torsion") for r in item.restraints)
