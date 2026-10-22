"""
Tests for seekrflow/modules/pqr.py
"""

import pytest
import parmed

import seekrflow.modules.pqr as pqr


def test_write_pqr_atoms_truncates_names_and_round_trips(tmp_path):
    atoms = [
        pqr.Pqr_atom(
            name="H27C",
            resname="CHL1",
            resid=360,
            x=5.504,
            y=-17.135,
            z=-1.666,
            charge=0.09,
            radius=1.2,
        ),
        pqr.Pqr_atom(
            name="O",
            resname="SER",
            resid=240,
            x=8.543,
            y=-5.519,
            z=-30.648,
            charge=-0.51,
            radius=1.5,
        ),
        pqr.Pqr_atom(
            name="H27CX",
            resname="ALA",
            resid=1,
            x=0.0,
            y=0.0,
            z=0.0,
            charge=0.1,
            radius=1.0,
        ),
    ]
    path = tmp_path / "test.pqr"
    pqr.write_pqr_atoms(atoms, str(path))
    lines = path.read_text().splitlines()
    assert len(lines) == 3
    assert "HETATM" not in path.read_text()
    tokens = lines[0].split()
    assert tokens[0] == "ATOM"
    assert tokens[3] == "CHL"
    assert tokens[4] == "360"
    parsed = pqr.read_pqr_atoms(str(path))
    assert parsed[0].resname == "CHL"
    assert parsed[0].name == "H27C"
    assert parsed[1].resname == "SER"
    assert parsed[2].name == "H27C"
    assert parsed[0].x == pytest.approx(5.504, abs=1e-3)
    assert parsed[0].charge == pytest.approx(0.09, abs=1e-4)


def test_write_pqr_atoms_rejects_empty():
    with pytest.raises(ValueError, match="no atoms"):
        pqr.write_pqr_atoms([], "unused.pqr")


def test_read_pqr_atoms_ignores_chain_and_element(tmp_path):
    path = tmp_path / "dialects.pqr"
    path.write_text(
        "REMARK ignored\n"
        "ATOM      1  N   ALA A  10      11.104  13.207  10.000 -0.4157 1.5500\n"
        "ATOM      2  CA  SER B  11      12.000  13.000  10.000  0.5973 1.7000"
        "           C\n"
    )
    parsed = pqr.read_pqr_atoms(str(path))
    assert len(parsed) == 2
    assert parsed[0].resname == "ALA"
    assert parsed[0].resid == 10
    assert parsed[0].x == pytest.approx(11.104, abs=1e-3)
    assert parsed[0].charge == pytest.approx(-0.4157, abs=1e-4)
    assert parsed[1].resname == "SER"
    assert parsed[1].resid == 11
    assert parsed[1].x == pytest.approx(12.0, abs=1e-3)
    assert parsed[1].charge == pytest.approx(0.5973, abs=1e-4)


def _cholesterol_fragment():
    struct = parmed.Structure()
    specs = (
        ("C1", 0.1, 6, 0.0),
        ("H1", 0.05, 1, 1.0),
        ("O1", -0.15, 8, 2.0),
    )
    for name, charge, znum, x in specs:
        atom = parmed.Atom(name=name, atomic_number=znum)
        atom.charge = charge
        atom.xx, atom.xy, atom.xz = x, x + 0.5, x + 1.0
        struct.add_atom(atom, "CHL1", 7)
    struct.atoms[0].solvent_radius = 1.25
    struct.atoms[1].solvent_radius = None
    struct.atoms[1].rmin = 0.8
    struct.atoms[2].solvent_radius = None
    struct.atoms[2].rmin = 0.0
    return struct


def test_structure_to_pqr_atoms_uses_radii_and_index_order():
    struct = _cholesterol_fragment()
    sliced = pqr.structure_to_pqr_atoms(struct, [2, 0])
    assert [atom.name for atom in sliced] == ["O1", "C1"]
    assert sliced[0].x == pytest.approx(2.0)
    assert sliced[0].y == pytest.approx(2.5)
    assert sliced[0].z == pytest.approx(3.0)
    assert sliced[0].charge == pytest.approx(-0.15)
    assert sliced[1].charge == pytest.approx(0.1)
    assert [atom.resid for atom in sliced] == [7, 7]
    assert sliced[0].radius == pytest.approx(0.0)
    assert sliced[1].radius == pytest.approx(1.25)
    rmin_only = pqr.structure_to_pqr_atoms(struct, [1])
    assert rmin_only[0].radius == pytest.approx(0.8)
    renumbered = pqr.structure_to_pqr_atoms(
        struct, [2, 0], one_resid_per_atom=True)
    assert [atom.resid for atom in renumbered] == [1, 2]
    with pytest.raises(IndexError):
        pqr.structure_to_pqr_atoms(struct, [3])


def test_write_pqr_from_parmed_copies_assigned_radii(tmp_path):
    struct = _cholesterol_fragment()
    path = tmp_path / "from_parmed.pqr"
    assert pqr.write_pqr_from_parmed(
        struct, [2, 0], str(path)) == str(path)
    parsed = pqr.read_pqr_atoms(str(path))
    assert [atom.radius for atom in parsed] == pytest.approx([
        struct.atoms[2].solvent_radius,
        struct.atoms[0].solvent_radius,
    ])
    with pytest.raises(AssertionError, match="no atoms"):
        pqr.write_pqr_from_parmed(struct, [], str(tmp_path / "empty.pqr"))
