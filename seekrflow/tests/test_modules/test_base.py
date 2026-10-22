"""
Tests for seekrflow/modules/base.py
"""

import textwrap

import pytest

import seekrflow.modules.base as base


def test_initialize_ref_indices_parses_nonnegative_list():
    assert base.initialize_ref_indices("0, 1,2") == [0, 1, 2]


def test_initialize_ref_indices_rejects_empty_noninteger_and_negative():
    with pytest.raises(AssertionError):
        base.initialize_ref_indices("")
    with pytest.raises(ValueError):
        base.initialize_ref_indices("1,2,abc,4")
    with pytest.raises(AssertionError):
        base.initialize_ref_indices("-1,0,1")


def _complex_pdb(tmp_path):
    """Two LIG residues, one hydrogen on the first, and an H-only residue."""
    contents = textwrap.dedent(
        """\
        ATOM      1  N   LIG A   1      11.104  13.207  10.000  1.00  0.00           N
        ATOM      2  C   LIG A   1      12.000  13.000  10.000  1.00  0.00           C
        ATOM      3  O   LIG A   1      13.000  13.000  10.000  1.00  0.00           O
        ATOM      4  H   LIG A   1      11.104  13.207  11.000  1.00  0.00           H
        ATOM      5  N   LIG A   2      21.104  13.207  10.000  1.00  0.00           N
        ATOM      6  C   LIG A   2      22.000  13.000  10.000  1.00  0.00           C
        ATOM      7  CA  PRO A   3      14.000  13.000  10.000  1.00  0.00           C
        ATOM      8  H1  HYD A   4      15.000  13.000  10.000  1.00  0.00           H
        ATOM      9  H2  HYD A   4      16.000  13.000  10.000  1.00  0.00           H
        TER
        END
        """
    )
    pdb_path = tmp_path / "complex.pdb"
    pdb_path.write_text(contents)
    return str(pdb_path)


def test_get_ligand_indices_selects_every_residue_and_optional_hydrogens(tmp_path):
    pdb_path = _complex_pdb(tmp_path)
    assert base.get_ligand_indices(pdb_path, "LIG") == [0, 1, 2, 4, 5]
    assert base.get_ligand_indices(pdb_path, "LIG", include_H=True) == [
        0, 1, 2, 3, 4, 5]


def test_get_ligand_indices_rejects_empty_selections(tmp_path):
    pdb_path = _complex_pdb(tmp_path)
    with pytest.raises(AssertionError):
        base.get_ligand_indices(pdb_path, "BEN")
    with pytest.raises(AssertionError):
        base.get_ligand_indices(pdb_path, "HYD", include_H=False)
