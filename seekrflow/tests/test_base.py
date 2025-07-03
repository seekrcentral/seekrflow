"""
test_base.py
"""

import pytest
import tempfile
import textwrap

import mdtraj
import numpy as np

import seekrflow.modules.base as base

def test_initialize_ref_indices_basic():
    result = base.initialize_ref_indices("1,2,3,4")
    assert result == [1, 2, 3, 4]

def test_initialize_ref_indices_single_value():
    result = base.initialize_ref_indices("42")
    assert result == [42]

def test_initialize_ref_indices_with_spaces():
    result = base.initialize_ref_indices(" 5, 6 ,7 ")
    assert result == [5, 6, 7]

def test_initialize_ref_indices_empty_string():
    with pytest.raises(AssertionError):
        base.initialize_ref_indices("")

def _write_pdb(contents):
    """Helper to write a temporary PDB file and return its path."""
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", mode="w")
    tmp.write(contents)
    tmp.close()
    return tmp.name

def test_get_ligand_indices_basic():
    pdb = textwrap.dedent("""
    ATOM      1  N   LIG A   1      11.104  13.207  10.000  1.00  0.00           N
    ATOM      2  C   LIG A   1      12.000  13.000  10.000  1.00  0.00           C
    ATOM      3  O   LIG A   1      13.000  13.000  10.000  1.00  0.00           O
    ATOM      4  H   LIG A   1      11.104  13.207  11.000  1.00  0.00           H
    ATOM      5  CA  PRO A   2      14.000  13.000  10.000  1.00  0.00           C
    TER
    END
    """)
    pdb_path = _write_pdb(pdb)
    indices = base.get_ligand_indices(pdb_path, "LIG")
    # Should exclude H (index 3) by default
    assert indices == [0, 1, 2]

def test_get_ligand_indices_include_H():
    pdb = textwrap.dedent("""
    ATOM      1  N   LIG A   1      11.104  13.207  10.000  1.00  0.00           N
    ATOM      2  C   LIG A   1      12.000  13.000  10.000  1.00  0.00           C
    ATOM      3  O   LIG A   1      13.000  13.000  10.000  1.00  0.00           O
    ATOM      4  H   LIG A   1      11.104  13.207  11.000  1.00  0.00           H
    ATOM      5  CA  PRO A   2      14.000  13.000  10.000  1.00  0.00           C
    TER
    END
    """)
    pdb_path = _write_pdb(pdb)
    indices = base.get_ligand_indices(pdb_path, "LIG", include_H=True)
    # Should include all LIG atoms (indices 0,1,2,3)
    assert indices == [0, 1, 2, 3]

import pytest

def test_get_ligand_indices_not_found():
    pdb = textwrap.dedent("""
    ATOM      1  N   PRO A   1      11.104  13.207  10.000  1.00  0.00           N
    ATOM      2  C   PRO A   1      12.000  13.000  10.000  1.00  0.00           C
    TER
    END
    """)
    pdb_path = _write_pdb(pdb)
    with pytest.raises(AssertionError):
        base.get_ligand_indices(pdb_path, "LIG")

def test_initialize_ref_indices_invalid_integer():
    """Test that invalid integers raise ValueError"""
    with pytest.raises(ValueError):
        base.initialize_ref_indices("1,2,abc,4")

def test_initialize_ref_indices_negative_numbers():
    """Test that negative numbers are rejected"""
    with pytest.raises(AssertionError):
        base.initialize_ref_indices("-1,0,1")

def test_initialize_ref_indices_single_negative():
    """Test that a single negative number is rejected"""
    with pytest.raises(AssertionError):
        base.initialize_ref_indices("-5")

def test_initialize_ref_indices_mixed_negative():
    """Test that mixed positive and negative numbers are rejected"""
    with pytest.raises(AssertionError):
        base.initialize_ref_indices("1,2,-3,4")

def test_initialize_ref_indices_zero_allowed():
    """Test that zero is allowed (valid array index)"""
    result = base.initialize_ref_indices("0,1,2")
    assert result == [0, 1, 2]

def test_get_ligand_indices_multiple_residues():
    """Test with multiple residues of the same name"""
    pdb = textwrap.dedent("""
    ATOM      1  N   LIG A   1      11.104  13.207  10.000  1.00  0.00           N
    ATOM      2  C   LIG A   1      12.000  13.000  10.000  1.00  0.00           C
    ATOM      3  N   LIG A   2      21.104  13.207  10.000  1.00  0.00           N
    ATOM      4  C   LIG A   2      22.000  13.000  10.000  1.00  0.00           C
    ATOM      5  CA  PRO A   3      14.000  13.000  10.000  1.00  0.00           C
    TER
    END
    """)
    pdb_path = _write_pdb(pdb)
    indices = base.get_ligand_indices(pdb_path, "LIG")
    # Should find all LIG atoms from both residues
    assert indices == [0, 1, 2, 3]

def test_get_ligand_indices_only_hydrogens():
    """Test ligand with only hydrogen atoms"""
    pdb = textwrap.dedent("""
    ATOM      1  H1  LIG A   1      11.104  13.207  10.000  1.00  0.00           H
    ATOM      2  H2  LIG A   1      12.000  13.000  10.000  1.00  0.00           H
    ATOM      3  CA  PRO A   2      14.000  13.000  10.000  1.00  0.00           C
    TER
    END
    """)
    pdb_path = _write_pdb(pdb)
    
    # Without include_H, should find no atoms and raise assertion
    with pytest.raises(AssertionError):
        base.get_ligand_indices(pdb_path, "LIG", include_H=False)
    
    # With include_H, should find hydrogen atoms
    indices = base.get_ligand_indices(pdb_path, "LIG", include_H=True)
    assert indices == [0, 1]

def test_get_ligand_indices_invalid_file():
    """Test behavior with invalid PDB file"""
    with pytest.raises(Exception):  # mdtraj will raise an exception
        base.get_ligand_indices("/nonexistent/file.pdb", "LIG")
