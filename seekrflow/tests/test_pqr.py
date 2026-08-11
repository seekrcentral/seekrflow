"""Tests for strict classic PQR write/read helpers."""

import os

import pytest

import seekrflow.modules.pqr as pqr


def test_write_pqr_atoms_truncates_resname_and_uses_atom_only(tmp_path):
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
    ]
    path = tmp_path / "test.pqr"
    pqr.write_pqr_atoms(atoms, str(path))
    text = path.read_text().splitlines()
    assert len(text) == 2
    assert text[0].startswith("ATOM")
    assert "HETATM" not in path.read_text()
    # No chain-id token between resname and resid for first record.
    tokens0 = text[0].split()
    assert tokens0[0] == "ATOM"
    assert tokens0[3] == "CHL"  # truncated
    assert tokens0[4] == "360"
    parsed = pqr.read_pqr_atoms(str(path))
    assert len(parsed) == 2
    assert parsed[0].resname == "CHL"
    assert parsed[0].name.strip() == "H27C"
    assert parsed[1].resname == "SER"
    assert parsed[0].x == pytest.approx(5.504, abs=1e-3)
    assert parsed[0].charge == pytest.approx(0.09, abs=1e-4)


def test_write_pqr_atoms_rejects_empty():
    with pytest.raises(ValueError, match="no atoms"):
        pqr.write_pqr_atoms([], "unused.pqr")


def test_structure_to_pqr_atoms_and_write_from_parmed(tmp_path):
    parmed = pytest.importorskip("parmed")
    struct = parmed.Structure()
    for i, (name, charge, znum) in enumerate(
            (("C1", 0.1, 6), ("H1", 0.05, 1), ("O1", -0.15, 8))):
        atom = parmed.Atom(name=name, atomic_number=znum)
        atom.charge = charge
        atom.xx, atom.xy, atom.xz = float(i), float(i) + 0.5, float(i) + 1.0
        struct.add_atom(atom, "CHL1", 7)

    path = tmp_path / "from_parmed.pqr"
    out = pqr.write_pqr_from_parmed(
        struct, [0, 1, 2], str(path), one_resid_per_atom=False)
    assert out == str(path)
    parsed = pqr.read_pqr_atoms(str(path))
    assert len(parsed) == 3
    assert all(a.resname == "CHL" for a in parsed)
    assert [a.resid for a in parsed] == [7, 7, 7]
    assert os.path.isfile(path)

    path2 = tmp_path / "one_resid.pqr"
    pqr.write_pqr_from_parmed(
        struct, [0, 1, 2], str(path2), one_resid_per_atom=True)
    parsed2 = pqr.read_pqr_atoms(str(path2))
    assert [a.resid for a in parsed2] == [1, 2, 3]
