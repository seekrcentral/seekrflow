"""
Tests for seekrflow/modules/cvs.py
"""

import os
import textwrap

import mdtraj
import parmed
import pytest

import seekrflow.modules.cvs as cvs
import seekrflow.modules.workflows.components as components_module
import seekrflow.tests.create_seekrflow as create_seekrflow


def _data_file(name):
    return os.path.join(
        create_seekrflow.TEST_DIRECTORY,
        "data",
        name)


def _cutoff_pdb(tmp_path):
    """Ligand, a close CA, a far CA, a CA-less protein residue, and a non-protein residue."""
    contents = textwrap.dedent(
        """\
        ATOM      1  C   LIG A   1      10.000  10.000  10.000  1.00  0.00           C
        ATOM      2  CA  ALA A   2      12.000  10.000  10.000  1.00  0.00           C
        ATOM      3  CA  ALA A   3      18.000  10.000  10.000  1.00  0.00           C
        ATOM      4  N   GLY A   4      11.000  10.000  10.000  1.00  0.00           N
        ATOM      5  C   HST A   5      11.500  10.000  10.000  1.00  0.00           C
        ATOM      6  H   HST A   5      11.600  10.000  10.000  1.00  0.00           H
        TER
        END
        """
    )
    pdb_path = tmp_path / "cutoff.pdb"
    pdb_path.write_text(contents)
    return mdtraj.load(str(pdb_path))


def test_alpha_carbon_cutoff_returns_contact_cas(tmp_path):
    traj = _cutoff_pdb(tmp_path)
    ligand = [int(i) for i in traj.topology.select("resname LIG")]
    close_ca = int(traj.topology.select("resname ALA and resid 1")[0])
    far_ca = int(traj.topology.select("resname ALA and resid 2")[0])
    no_ca = int(traj.topology.select("resname GLY")[0])

    near = cvs.alpha_carbon_selection_within_cutoff(
        traj, ligand, distance_cutoff=0.3)
    assert near == [close_ca]
    assert no_ca not in near
    assert cvs.alpha_carbon_selection_within_cutoff(
        traj, ligand, distance_cutoff=1.0) == [close_ca, far_ca]
    assert cvs.alpha_carbon_selection_within_cutoff(
        traj, ligand, distance_cutoff=0.0) == []
    assert cvs.alpha_carbon_selection_within_cutoff(
        traj, [], distance_cutoff=1.0) == []


def test_alpha_carbon_cutoff_falls_back_when_residue_has_no_ca(tmp_path):
    traj = _cutoff_pdb(tmp_path)
    ligand = [int(i) for i in traj.topology.select("resname LIG")]
    heavy = int(traj.topology.select("resname HST and not element H")[0])
    hydrogen = int(traj.topology.select("resname HST and element H")[0])
    selected = cvs.alpha_carbon_selection_within_cutoff(
        traj,
        ligand,
        distance_cutoff=0.3,
        receptor_selection="resname HST")
    assert selected == [heavy]
    assert hydrogen not in selected


def _trypsin_seekrflow(resname):
    return create_seekrflow.create_unparameterized_seekrflow(
        _data_file("trypsin_benzamidine.pdb"), resname, 12, [0.5, 1.0])


def test_com_com_selections_follow_ligand_component():
    pdb_file = _data_file("trypsin_benzamidine.pdb")
    traj = mdtraj.load(pdb_file)
    seekrflow = _trypsin_seekrflow("BEN")
    receptor_sel, ligand_sel = cvs.get_receptor_ligand_com_com_selections(
        seekrflow.workflow, pdb_file)
    expected_ligand = [int(i) for i in traj.topology.select("resname BEN")]
    assert ligand_sel == expected_ligand
    assert receptor_sel
    assert receptor_sel == sorted(receptor_sel)
    for index in receptor_sel:
        assert traj.topology.atom(index).name == "CA"

    ligand_indices = expected_ligand[:4]
    ligand = seekrflow.workflow.components.get_member("ligand")
    ligand.selector = components_module.Selector_by_indices(
        indices=ligand_indices)
    _, index_ligand_sel = cvs.get_receptor_ligand_com_com_selections(
        seekrflow.workflow, pdb_file)
    assert index_ligand_sel == ligand_indices


def test_com_com_selections_reject_missing_ligand_atoms():
    pdb_file = _data_file("trypsin_benzamidine.pdb")
    missing_resname = _trypsin_seekrflow("XXX")
    with pytest.raises(AssertionError, match="resolved to no atoms"):
        cvs.get_receptor_ligand_com_com_selections(
            missing_resname.workflow, pdb_file)

    protein_only = _trypsin_seekrflow("BEN")
    protein_only.workflow.components.members = [
        member for member in protein_only.workflow.components.members
        if not isinstance(member, components_module.Small_molecule_component)]
    with pytest.raises(AssertionError, match="Small_molecule_component"):
        cvs.get_receptor_ligand_com_com_selections(
            protein_only.workflow, pdb_file)


def test_com_com_selections_forward_the_distance_threshold():
    pdb_file = _data_file("trypsin_benzamidine.pdb")
    seekrflow = _trypsin_seekrflow("BEN")
    receptor_small, _ = cvs.get_receptor_ligand_com_com_selections(
        seekrflow.workflow, pdb_file, alpha_carbon_ligand_threshold=0.3)
    receptor_large, _ = cvs.get_receptor_ligand_com_com_selections(
        seekrflow.workflow, pdb_file, alpha_carbon_ligand_threshold=0.8)
    assert len(receptor_large) >= len(receptor_small)


def _bd_paths():
    return (
        _data_file("tryp_ben_receptor.pqr"),
        _data_file("tryp_ben_ligand.pqr"),
        _data_file("tryp_ben_system_bound.pdb"),
    )


def test_bd_selections_map_complex_atom_names_onto_pqr_indices():
    receptor_pqr, ligand_pqr, pdb_file = _bd_paths()
    traj = mdtraj.load(pdb_file)
    md_receptor = [int(i) for i in traj.topology.select("protein and name CA")][:5]
    md_ligand = [int(i) for i in traj.topology.select("resname BEN")][:3]
    bd_receptor, bd_ligand = cvs.get_bd_receptor_ligand_selections(
        receptor_pqr, ligand_pqr, pdb_file, md_receptor, md_ligand)
    assert len(bd_receptor) == len(md_receptor)
    assert len(bd_ligand) == len(md_ligand)
    receptor_struct = parmed.load_file(receptor_pqr)
    ligand_struct = parmed.load_file(ligand_pqr)
    for md_index, bd_index in zip(md_receptor, bd_receptor):
        assert receptor_struct.atoms[bd_index].name.strip() \
            == traj.topology.atom(md_index).name
    for md_index, bd_index in zip(md_ligand, bd_ligand):
        assert ligand_struct.atoms[bd_index].name.strip() \
            == traj.topology.atom(md_index).name


def test_bd_selections_reject_unmatched_atoms_and_empty_input():
    receptor_pqr, ligand_pqr, pdb_file = _bd_paths()
    traj = mdtraj.load(pdb_file)
    ca_index = int(traj.topology.select("protein and name CA")[0])
    ben_index = int(traj.topology.select("resname BEN")[0])
    with pytest.raises(AssertionError, match="No match found for receptor atom"):
        cvs.get_bd_receptor_ligand_selections(
            receptor_pqr, ligand_pqr, pdb_file, [ben_index], [ben_index])
    with pytest.raises(AssertionError, match="No match found for ligand atom"):
        cvs.get_bd_receptor_ligand_selections(
            receptor_pqr, ligand_pqr, pdb_file, [ca_index], [ca_index])
    assert cvs.get_bd_receptor_ligand_selections(
        receptor_pqr, ligand_pqr, pdb_file, [], []) == ([], [])
    with pytest.raises(AssertionError):
        cvs.get_bd_receptor_ligand_selections(
            "nonexistent_receptor.pqr",
            ligand_pqr,
            pdb_file,
            [ca_index],
            [ben_index])
