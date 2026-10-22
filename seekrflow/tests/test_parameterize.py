"""
Tests for seekrflow/parameterize.py
"""

import os
import textwrap

import mdtraj
import numpy as np
import openmm.unit as unit
import parmed
import pytest

import seekrflow.modules.client.structures as client_structures
import seekrflow.modules.parameterize_structures as parameterize_structures
import seekrflow.modules.parameters_topology_structures as pts
import seekrflow.modules.structures as structures
import seekrflow.modules.workflows.scale_settings as scale_settings
import seekrflow.parameterize as parameterize
import seekrflow.tests.create_seekrflow as create_seekrflow


_COMPLEX_PDB = textwrap.dedent("""
    ATOM      1  N   ALA A   1      -6.060   2.624   4.140  1.00  0.00           N
    ATOM      2  CA  ALA A   1      -5.287   1.688   3.297  1.00  0.00           C
    ATOM      3  C   ALA A   1      -4.066   1.131   4.042  1.00  0.00           C
    ATOM      4  O   ALA A   1      -3.287   1.688   4.642  1.00  0.00           O
    ATOM      5  CB  ALA A   1      -4.444   0.853   2.525  1.00  0.00           C
    ATOM      6  N   LIG B   1      11.104  13.207  10.000  1.00  0.00           N
    ATOM      7  C   LIG B   1      12.000  13.000  10.000  1.00  0.00           C
    ATOM      8  O   LIG B   1      13.000  13.000  10.000  1.00  0.00           O
    ATOM      9  H   LIG B   1      11.104  13.207  11.000  1.00  0.00           H
    TER
    END
    """).strip() + "\n"


class _FakeOpenMM:
    def chains(self):
        return []


class _FakeTop:
    def __init__(self, n_atoms):
        self.n_atoms = n_atoms

    def join(self, other):
        return _FakeTop(self.n_atoms + other.n_atoms)

    def to_openmm(self):
        return _FakeOpenMM()


def _write_complex(tmp_path):
    pdb = tmp_path / "complex.pdb"
    pdb.write_text(_COMPLEX_PDB)
    return pdb


def _flow(pdb, resname="LIG"):
    return create_seekrflow.create_unparameterized_seekrflow(
        str(pdb), resname, 12, [0.5, 1.0])


def _batch(flow, batch_directory=None):
    batch = client_structures.Batch(
        batch_directory=batch_directory,
        template={},
        systems=[client_structures.Batch_system(name=flow.name)],
    )
    batch._existing_seekrflow = flow
    return batch


def _disable_structure_prep(flow):
    flow.parameterizer.pdb_fixer_settings = None
    flow.parameterizer.pdb2pqr_settings = None
    flow.parameterizer.pdb_to_sdf_settings = None


def test_split_components_writes_receptor_and_ligand(tmp_path):
    pdb = _write_complex(tmp_path)
    param = tmp_path / "parameterize"
    param.mkdir()
    resolved = parameterize.split_components(_flow(pdb), str(param), str(pdb))
    receptor = parmed.load_file(str(param / "receptor.pdb"))
    ligand = parmed.load_file(str(param / "ligand.pdb"))
    assert len(receptor.atoms) == 5
    assert len(ligand.atoms) == 4
    assert set(resolved) == {"receptor", "ligand"}


def test_split_components_empty_ligand_selection(tmp_path):
    pdb = _write_complex(tmp_path)
    param = tmp_path / "parameterize"
    param.mkdir()
    with pytest.raises(AssertionError, match="ligand.*resolved to no atoms"):
        parameterize.split_components(_flow(pdb, "XXX"), str(param), str(pdb))


def test_parameterize_rejects_populated_system(tmp_path):
    pdb = _write_complex(tmp_path)
    flow = _flow(pdb)
    flow.workflow.get_md_settings().system = scale_settings.System_for_md(
        parameters_topology=pts.Openmm_system(system_filename="already.xml"),
        solvated_pdb="already.pdb",
    )
    with pytest.raises(Exception, match="already populated"):
        parameterize.parameterize(_batch(flow))


def test_parameterize_none_returns_without_complex_pdb(tmp_path):
    flow = _flow(tmp_path / "missing.pdb")
    flow.parameterizer = None
    assert parameterize.parameterize(_batch(flow)) == (None, None)


def test_parameterize_copies_pdb_strips_hydrogens_and_sdf(tmp_path, monkeypatch):
    pdb = _write_complex(tmp_path)
    sdf = tmp_path / "ligand_src.sdf"
    sdf.write_text("mock sdf\n")
    work = tmp_path / "work"
    flow = _flow(pdb)
    flow.make_work_directory(work)
    ligand = flow.workflow.components.get_member("ligand")
    ligand.sdf_file = str(sdf)
    batch_dir = tmp_path / "batch"
    batch_dir.mkdir()
    batch = _batch(flow, str(batch_dir))
    monkeypatch.setattr(
        client_structures.Batch, "create_seekrflow_objects",
        lambda self: [flow])
    monkeypatch.setattr(parameterize, "split_components", lambda *args: {})
    monkeypatch.setattr(
        parameterize, "create_complex",
        lambda *args: ("complex_system.xml", "complex_solvent.pdb"))
    cwd = os.getcwd()
    result = parameterize.parameterize(batch)
    assert os.getcwd() == cwd
    assert result == ("complex_system.xml", "complex_solvent.pdb")
    assert (work / "complex.pdb").read_text() == pdb.read_text()
    noh = mdtraj.load(str(work / "complex_noh.pdb"))
    assert noh.n_atoms == 8
    assert len(noh.topology.select("element H")) == 0
    assert "mock sdf" in (work / "parameterize" / "ligand.sdf").read_text()
    assert ligand.sdf_file == os.path.join(
        structures.PARAMETERIZE, "ligand.sdf")
    assert flow.parameterizer.complex_pdb_filename == "complex.pdb"


def test_fill_md_system_from_parameterize_outputs(tmp_path):
    pdb = _write_complex(tmp_path)
    flow = _flow(pdb)
    work = tmp_path / "work"
    flow.make_work_directory(work)
    with pytest.raises(FileNotFoundError):
        parameterize.fill_md_system_from_parameterize(flow)

    param = work / structures.PARAMETERIZE
    param.mkdir()
    (param / parameterize_structures.COMPLEX_SYSTEM_XML).write_text("<System/>")
    (param / parameterize_structures.COMPLEX_SOLVENT_PDB).write_text("END\n")
    parameterize.fill_md_system_from_parameterize(flow)
    system = flow.workflow.get_md_settings().system
    assert system.solvated_pdb == os.path.join(
        structures.PARAMETERIZE, parameterize_structures.COMPLEX_SOLVENT_PDB)
    assert system.parameters_topology.system_filename == os.path.join(
        structures.PARAMETERIZE, parameterize_structures.COMPLEX_SYSTEM_XML)

    populated = _flow(pdb)
    populated.make_work_directory(tmp_path / "work2")
    populated.workflow.get_md_settings().system = scale_settings.System_for_md(
        parameters_topology=pts.Openmm_system(system_filename="already.xml"),
        solvated_pdb="already.pdb",
    )
    parameterize.fill_md_system_from_parameterize(populated)
    kept = populated.workflow.get_md_settings().system
    assert kept.solvated_pdb == "already.pdb"
    assert kept.parameters_topology.system_filename == "already.xml"


def test_create_complex_missing_sdf_with_null_settings_raises(tmp_path):
    pdb = _write_complex(tmp_path)
    flow = _flow(pdb)
    _disable_structure_prep(flow)
    flow.workflow.components.get_member("ligand").sdf_file = ""
    param = tmp_path / "parameterize"
    param.mkdir()
    parameterize.split_components(flow, str(param), str(pdb))
    with pytest.raises(ValueError, match="pdb_to_sdf_settings is None"):
        parameterize.create_complex(flow, str(param))


def test_create_complex_provided_sdf_skips_conversion(tmp_path, monkeypatch):
    pdb = _write_complex(tmp_path)
    sdf = tmp_path / "ligand_src.sdf"
    sdf.write_text("mock sdf\n")
    flow = _flow(pdb)
    _disable_structure_prep(flow)
    flow.workflow.components.get_member("ligand").sdf_file = str(sdf)
    param = tmp_path / "parameterize"
    param.mkdir()
    parameterize.split_components(flow, str(param), str(pdb))

    def boom(*_args, **_kwargs):
        raise AssertionError("PDB to SDF conversion should not run")

    def fake_protein(_pdb_filename):
        return None, _FakeTop(5), np.zeros([5, 3]) * unit.nanometers

    def fake_ligand(
            ligand_sdf_filename, ligand_pdb_filename, param_directory,
            draw_ligand=True, ligand_resname="LIG", chain=None):
        with open(ligand_sdf_filename) as handle:
            assert "mock sdf" in handle.read()
        return (
            None,
            _FakeTop(4),
            np.zeros([4, 3]) * unit.nanometers,
            "fake_offmol",
        )

    monkeypatch.setattr(parameterize, "make_small_molecule_sdf_from_pdb", boom)
    monkeypatch.setattr(
        parameterize.base, "make_protein_openmm_and_mdtraj_top", fake_protein)
    monkeypatch.setattr(
        parameterize.base, "make_ligand_openmm_and_mdtraj_top", fake_ligand)
    monkeypatch.setattr(
        parameterize.parameterize_algorithms,
        "parameterize_and_check_complex",
        lambda *args, **kwargs: ("complex_system.xml", "complex_solvent.pdb"),
    )
    assert parameterize.create_complex(flow, str(param)) == (
        "complex_system.xml", "complex_solvent.pdb")
    assert "mock sdf" in (param / "ligand.sdf").read_text()
