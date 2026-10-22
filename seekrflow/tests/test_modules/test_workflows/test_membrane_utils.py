"""
Tests for seekrflow/modules/workflows/membrane_utils.py
"""

import textwrap

import mdtraj
import numpy as np

import seekrflow.modules.workflows.membrane_utils as membrane_utils


_POPC_FRAGMENT = textwrap.dedent("""
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

_POPC_BONDS = (
    (0, 1), (1, 2), (1, 3), (3, 5), (4, 0),
    (5, 6), (6, 7), (6, 9), (7, 8), (7, 10),
)


def write_popc_fragment(path):
    path.write_text(_POPC_FRAGMENT)
    traj = mdtraj.load(str(path))
    if traj.topology.n_bonds == 0:
        for i, j in _POPC_BONDS:
            traj.topology.add_bond(
                traj.topology.atom(i), traj.topology.atom(j))
    return traj


def test_identify_lipid_atoms_finds_popc_groups(tmp_path):
    traj = write_popc_fragment(tmp_path / "lipid.pdb")
    result = membrane_utils.identify_lipid_atoms(traj)
    assert len(result["phosphorus"]) == 1
    assert traj.topology.atom(result["phosphorus"][0]).name == "P"
    assert len(result["glycerol_impropers"]) == 1
    improper_names = [
        traj.topology.atom(index).name
        for index in result["glycerol_impropers"][0]]
    assert improper_names == ["C1", "C3", "C2", "O21"]
    assert len(result["cis_double_bonds"]) == 1
    cis_names = [
        traj.topology.atom(index).name
        for index in result["cis_double_bonds"][0]]
    assert cis_names == ["C21", "C22", "C23", "C24"]


def test_identify_lipid_atoms_empty_without_phosphorus(tmp_path):
    pdb = tmp_path / "carbon.pdb"
    pdb.write_text(textwrap.dedent(
        """\
        ATOM      1  C   XXX A   1       0.000   0.000   0.000  1.00  0.00           C
        END
        """
    ))
    result = membrane_utils.identify_lipid_atoms(mdtraj.load(str(pdb)))
    assert result["phosphorus"] == []
    assert result["cis_double_bonds"] == []
    assert result["glycerol_impropers"] == []


def test_glycerol_improper_theta0_picks_nearest_sign():
    positions = np.zeros((4, 3))
    positions[0] = [0.0, 0.0, 0.0]
    positions[1] = [1.0, 1.0, 0.0]
    positions[2] = [1.0, 0.0, 0.0]
    positions[3] = [1.0, -0.5, 0.866]
    improper = (0, 1, 2, 3)
    assert membrane_utils.glycerol_improper_theta0(
        positions, improper) == membrane_utils.GLYCEROL_IMPROPER_THETA0_NEG
    flipped = positions.copy()
    flipped[3, 2] *= -1
    assert membrane_utils.glycerol_improper_theta0(
        flipped, improper) == membrane_utils.GLYCEROL_IMPROPER_THETA0_POS
