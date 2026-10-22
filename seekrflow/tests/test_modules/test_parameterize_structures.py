"""
Tests for seekrflow/modules/parameterize_structures.py
"""

import textwrap

import pytest

import seekrflow.modules.parameterize_structures as parameterize_structures


_TWO_CHAIN_PDB = textwrap.dedent("""
    ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N
    ATOM      2  CA  ALA A   1      1.500   0.000   0.000  1.00  0.00           C
    ATOM      3  O   LIG B   1      5.000   0.000   0.000  1.00  0.00           O
    END
    """).strip() + "\n"

_METHANE_PDB = textwrap.dedent("""
    ATOM      1  C   MET A   1      0.000   0.000   0.000  1.00  0.00           C
    ATOM      2  H1  MET A   1      1.090   0.000   0.000  1.00  0.00           H
    ATOM      3  H2  MET A   1     -0.363   1.028   0.000  1.00  0.00           H
    ATOM      4  H3  MET A   1     -0.363  -0.514   0.890  1.00  0.00           H
    ATOM      5  H4  MET A   1     -0.363  -0.514  -0.890  1.00  0.00           H
    CONECT    1    2    3    4    5
    CONECT    2    1
    CONECT    3    1
    CONECT    4    1
    CONECT    5    1
    END
    """).strip() + "\n"


def test_describe_pdb_chains_counts_each_chain(tmp_path):
    pdb = tmp_path / "two_chain.pdb"
    pdb.write_text(_TWO_CHAIN_PDB)
    chains = parameterize_structures.describe_pdb_chains(str(pdb))
    by_index = {chain["index"]: chain for chain in chains}
    assert by_index[0]["n_atoms"] == 2
    assert by_index[0]["residue_names"] == ["ALA"]
    assert by_index[1]["n_atoms"] == 1
    assert by_index[1]["residue_names"] == ["LIG"]


def test_pdb_to_sdf_rdkit_writes_methane(tmp_path):
    pytest.importorskip("rdkit")
    from rdkit import Chem

    pdb = tmp_path / "methane.pdb"
    sdf = tmp_path / "methane.sdf"
    pdb.write_text(_METHANE_PDB)
    parameterize_structures.PDB_to_SDF_settings(engine="rdkit").run(
        str(pdb), str(sdf))
    mol = Chem.SDMolSupplier(str(sdf), removeHs=False)[0]
    assert mol is not None
    assert mol.GetNumAtoms() == 5
