"""
Tests for seekrflow/modules/workflows/utilities.py
"""

import numpy as np
import mdtraj
import pytest

import seekrflow.modules.workflows.utilities as utilities


def _trajectory(chains, xyz):
    """
    ``chains`` is a list of ``(chain_id, residues)``. Each residue is
    ``(resname, res_seq, [(atom_name, element), ...])``.
    """
    topology = mdtraj.Topology()
    for chain_id, residues in chains:
        chain = topology.add_chain(chain_id=chain_id)
        for resname, res_seq, atoms in residues:
            residue = topology.add_residue(resname, chain, resSeq=res_seq)
            for atom_name, element in atoms:
                topology.add_atom(atom_name, element, residue)
    return mdtraj.Trajectory(np.asarray(xyz, dtype=np.float32), topology)


def test_com_reference_matches_chain_and_resseq():
    md_traj = _trajectory(
        [(None, [("ALA", 4, [("CA", mdtraj.element.carbon)])])],
        [[[0.0, 0.0, 0.0]]],
    )
    com_traj = _trajectory(
        [(None, [("ALA", 4, [("CX", mdtraj.element.carbon)])])],
        [[[1.0, 2.0, 3.0]]],
    )
    groups, coords = utilities.com_reference_groups(com_traj, md_traj)
    assert len(groups) == 1
    assert groups[0].chain_id == ""
    assert groups[0].res_seq == 4
    assert groups[0].atom_indices == [0]
    np.testing.assert_allclose(coords, [[1.0, 2.0, 3.0]])


def test_com_reference_rejects_ambiguous_chain_and_resseq():
    md_traj = _trajectory(
        [
            ("A", [("ALA", 1, [("CA", mdtraj.element.carbon)])]),
            ("A", [("GLY", 1, [("CA", mdtraj.element.carbon)])]),
        ],
        [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]],
    )
    com_traj = _trajectory(
        [("A", [("ALA", 1, [("CX", mdtraj.element.carbon)])])],
        [[[0.0, 0.0, 0.0]]],
    )
    with pytest.raises(ValueError, match="matches 2 solvated residues"):
        utilities.com_reference_groups(com_traj, md_traj)
