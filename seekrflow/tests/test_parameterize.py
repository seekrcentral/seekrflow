"""
Test the automated parameterization scripts.
"""

import pytest
import tempfile
import os
import pathlib
import textwrap
import parmed

import seekrflow.modules.base as base
import seekrflow.modules.structures as structures
import seekrflow.modules.parameterize_structures as parameterize_structures
import seekrflow.tests.create_seekrflow as create_seekrflow
import seekrflow.parameterize as parameterize


@pytest.mark.needs_openff
def test_parameterize_amber_tryp_ben(tryp_ben_seekrflow_amber_unparameterized):
    pytest.importorskip("openff.toolkit")
    flow = tryp_ben_seekrflow_amber_unparameterized
    flow.make_work_directory(flow.work_directory)
    system_filename, positions_filename = parameterize.parameterize(flow)
    assert system_filename is not None
    assert positions_filename is not None
    assert os.path.exists(os.path.join(
        flow.work_directory, structures.PARAMETERIZE,
        parameterize_structures.COMPLEX_SYSTEM_XML))
    assert os.path.exists(os.path.join(
        flow.work_directory, structures.PARAMETERIZE,
        parameterize_structures.COMPLEX_SOLVENT_PDB))


@pytest.mark.needs_espaloma
@pytest.mark.needs_openff
@pytest.mark.skip
def test_parameterize_espaloma_tryp_ben(
        tryp_ben_seekrflow_espaloma_unparameterized):
    flow = tryp_ben_seekrflow_espaloma_unparameterized
    flow.make_work_directory(flow.work_directory)
    system_filename, positions_filename = parameterize.parameterize(flow)
    assert system_filename is not None
    assert positions_filename is not None


class TestParameterizeHelperFunctions:
    """Test individual helper functions in parameterize module"""

    def _create_test_pdb(self):
        """Helper to create a test PDB file with receptor and ligand"""
        pdb_content = textwrap.dedent("""
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
        """).strip()
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", mode="w")
        tmp.write(pdb_content)
        tmp.close()
        return tmp.name

    def test_split_components_basic(self):
        """Test basic component splitting functionality"""
        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_file = self._create_test_pdb()
            flow = create_seekrflow.create_unparameterized_seekrflow(
                pdb_file, "LIG", [0.5, 1.0]
            )
            flow.parameterizer.complex_pdb_filename = pdb_file
            flow.make_work_directory(pathlib.Path(tmpdir))
            os.chdir(flow.work_directory)
            os.makedirs(structures.PARAMETERIZE, exist_ok=True)
            resolved = parameterize.split_components(
                flow, structures.PARAMETERIZE, pdb_file)
            receptor_file = os.path.join(
                structures.PARAMETERIZE, "receptor.pdb")
            ligand_file = os.path.join(structures.PARAMETERIZE, "ligand.pdb")
            assert os.path.exists(receptor_file)
            assert os.path.exists(ligand_file)
            receptor_structure = parmed.load_file(receptor_file)
            ligand_structure = parmed.load_file(ligand_file)
            assert len(ligand_structure.atoms) == 4
            assert len(receptor_structure.atoms) == 5
            assert "ligand" in resolved
            os.unlink(pdb_file)

    def test_split_components_empty_ligand_selection(self):
        """Test split function when ligand selector matches nothing"""
        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_file = self._create_test_pdb()
            flow = create_seekrflow.create_unparameterized_seekrflow(
                pdb_file, "XXX", [0.5, 1.0]
            )
            flow.parameterizer.complex_pdb_filename = pdb_file
            flow.make_work_directory(pathlib.Path(tmpdir))
            os.chdir(flow.work_directory)
            os.makedirs(structures.PARAMETERIZE, exist_ok=True)
            with pytest.raises(AssertionError, match="resolved to no atoms"):
                parameterize.split_components(
                    flow, structures.PARAMETERIZE, pdb_file)
            os.unlink(pdb_file)

    def test_parameterize_rejects_populated_system(self):
        """parameterizer + populated MD system must raise."""
        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_file = self._create_test_pdb()
            flow = create_seekrflow.create_unparameterized_seekrflow(
                pdb_file, "LIG", [0.5, 1.0]
            )
            flow.make_work_directory(pathlib.Path(tmpdir))
            import seekrflow.modules.workflows.scale_settings as scale_settings
            import seekrflow.modules.parameters_topology_structures as pts
            flow.workflow.get_md_settings().system = \
                scale_settings.System_for_md(
                    parameters_topology=pts.Openmm_system(
                        system_filename="already.xml"),
                    solvated_pdb="already.pdb",
                )
            with pytest.raises(Exception, match="already populated"):
                parameterize.parameterize(flow)
            os.unlink(pdb_file)


class TestParameterizeIntegration:
    """Integration tests for parameterize file handling"""

    @pytest.mark.needs_openff
    def test_parameterize_file_copying(self):
        """Test that input files are properly copied to work directory"""
        with tempfile.TemporaryDirectory() as tmpdir:
            test_data_dir = os.path.join(os.path.dirname(__file__), "data")
            test_pdb = os.path.join(test_data_dir, "trypsin_benzamidine.pdb")
            if not os.path.exists(test_pdb):
                pytest.skip(f"Test file {test_pdb} not found")

            original_sdf = os.path.join(tmpdir, "ligand.sdf")
            with open(original_sdf, "w") as f:
                f.write("mock sdf content\n")

            flow = create_seekrflow.create_unparameterized_seekrflow(
                test_pdb, "BEN", [0.5, 1.0]
            )
            ligand = flow.workflow.components.get_member("ligand")
            ligand.sdf_file = original_sdf
            work_dir = os.path.join(tmpdir, "work")
            flow.make_work_directory(pathlib.Path(work_dir))

            from shutil import copyfile
            work_copy_pdb = os.path.join(work_dir, os.path.basename(test_pdb))
            work_copy_sdf = os.path.join(work_dir, "parameterize", "ligand.sdf")
            os.makedirs(os.path.join(work_dir, "parameterize"), exist_ok=True)
            copyfile(test_pdb, work_copy_pdb)
            copyfile(original_sdf, work_copy_sdf)
            assert os.path.exists(work_copy_pdb)
            assert os.path.exists(work_copy_sdf)
            with open(test_pdb, "r") as f1, open(work_copy_pdb, "r") as f2:
                assert f1.read() == f2.read()


def test_describe_pdb_chains():
    test_data_dir = os.path.join(os.path.dirname(__file__), "data")
    test_pdb = os.path.join(test_data_dir, "trypsin_benzamidine.pdb")
    if not os.path.exists(test_pdb):
        pytest.skip(f"Test file {test_pdb} not found")
    chains = parameterize_structures.describe_pdb_chains(test_pdb)
    assert len(chains) >= 1
    assert "index" in chains[0]
    assert "n_atoms" in chains[0]


def test_pdb_to_sdf_rdkit_engine():
    """RDKit PDB→SDF conversion on a tiny ligand fragment."""
    pdb_content = textwrap.dedent("""
    ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C  
    ATOM      2  O   LIG A   1       1.200   0.000   0.000  1.00  0.00           O  
    CONECT    1    2
    END
    """).strip()
    with tempfile.TemporaryDirectory() as tmpdir:
        pdb_path = os.path.join(tmpdir, "lig.pdb")
        sdf_path = os.path.join(tmpdir, "lig.sdf")
        with open(pdb_path, "w") as f:
            f.write(pdb_content)
        settings = parameterize_structures.PDB_to_SDF_settings(engine="rdkit")
        settings.run(pdb_path, sdf_path)
        assert os.path.exists(sdf_path)
        assert os.path.getsize(sdf_path) > 0


def test_pdb_to_sdf_settings_none_allowed():
    """Parameterizer accepts pdb_to_sdf_settings=None (disable conversion)."""
    param = parameterize_structures.Parameterizer(
        complex_pdb_filename="x.pdb",
        pdb_to_sdf_settings=None,
    )
    assert param.pdb_to_sdf_settings is None


def test_missing_sdf_with_null_settings_raises():
    """No sdf_file + pdb_to_sdf_settings=None must raise a clear error."""
    pdb_content = textwrap.dedent("""
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
    """).strip()
    with tempfile.TemporaryDirectory() as tmpdir:
        pdb_file = os.path.join(tmpdir, "complex.pdb")
        with open(pdb_file, "w") as f:
            f.write(pdb_content)
        flow = create_seekrflow.create_unparameterized_seekrflow(
            pdb_file, "LIG", [0.5, 1.0]
        )
        flow.parameterizer.complex_pdb_filename = pdb_file
        flow.parameterizer.pdb_to_sdf_settings = None
        flow.parameterizer.pdb_fixer_settings = None
        flow.parameterizer.pdb2pqr_settings = None
        ligand = flow.workflow.components.get_member("ligand")
        ligand.sdf_file = ""
        flow.make_work_directory(pathlib.Path(tmpdir) / "work")
        os.chdir(flow.work_directory)
        os.makedirs(structures.PARAMETERIZE, exist_ok=True)
        parameterize.split_components(
            flow, structures.PARAMETERIZE, pdb_file)
        with pytest.raises(ValueError, match="pdb_to_sdf_settings is None"):
            parameterize.create_complex(flow, structures.PARAMETERIZE)


def test_provided_sdf_with_null_settings_skips_conversion(monkeypatch):
    """Explicit sdf_file + pdb_to_sdf_settings=None skips PDB→SDF conversion."""
    pdb_content = textwrap.dedent("""
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
    """).strip()
    with tempfile.TemporaryDirectory() as tmpdir:
        pdb_file = os.path.join(tmpdir, "complex.pdb")
        sdf_file = os.path.join(tmpdir, "ligand.sdf")
        with open(pdb_file, "w") as f:
            f.write(pdb_content)
        with open(sdf_file, "w") as f:
            f.write("mock sdf\n")
        flow = create_seekrflow.create_unparameterized_seekrflow(
            pdb_file, "LIG", [0.5, 1.0]
        )
        flow.parameterizer.complex_pdb_filename = pdb_file
        flow.parameterizer.pdb_to_sdf_settings = None
        flow.parameterizer.pdb_fixer_settings = None
        flow.parameterizer.pdb2pqr_settings = None
        ligand = flow.workflow.components.get_member("ligand")
        ligand.sdf_file = sdf_file
        flow.make_work_directory(pathlib.Path(tmpdir) / "work")
        os.chdir(flow.work_directory)
        os.makedirs(structures.PARAMETERIZE, exist_ok=True)
        parameterize.split_components(
            flow, structures.PARAMETERIZE, pdb_file)

        def boom(*_args, **_kwargs):
            raise AssertionError("PDB→SDF conversion should not run")

        monkeypatch.setattr(
            parameterize, "make_small_molecule_sdf_from_pdb", boom)

        class _FakeTop:
            def __init__(self, n_atoms):
                self.n_atoms = n_atoms

            def join(self, other):
                return _FakeTop(self.n_atoms + other.n_atoms)

            def to_openmm(self):
                return "fake_openmm_top"

        import numpy as np
        import openmm.unit as unit

        def fake_protein(pdb_filename):
            return None, _FakeTop(5), np.zeros([5, 3]) * unit.nanometers

        def fake_ligand(
                ligand_sdf_filename, ligand_pdb_filename, param_directory,
                draw_ligand=True, ligand_resname="LIG"):
            assert os.path.exists(ligand_sdf_filename)
            with open(ligand_sdf_filename) as f:
                assert "mock sdf" in f.read()
            return (
                None,
                _FakeTop(4),
                np.zeros([4, 3]) * unit.nanometers,
                "fake_offmol",
            )

        monkeypatch.setattr(
            base, "make_protein_openmm_and_mdtraj_top", fake_protein)
        monkeypatch.setattr(
            base, "make_ligand_openmm_and_mdtraj_top", fake_ligand)
        monkeypatch.setattr(
            parameterize.parameterize_algorithms,
            "parameterize_and_check_complex",
            lambda *a, **k: ("complex_system.xml", "complex_solvent.pdb"),
        )
        out = parameterize.create_complex(flow, structures.PARAMETERIZE)
        assert out == ("complex_system.xml", "complex_solvent.pdb")
        assert os.path.exists(
            os.path.join(structures.PARAMETERIZE, "ligand.sdf"))
