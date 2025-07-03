"""
test_structures.py
"""

import os
import tempfile
import json
import pytest
import pathlib

import seekrflow.modules.structures as structures

def test_load_save_seekrflow(tryp_ben_seekrflow_system_xml_parametrized):
    myflow = tryp_ben_seekrflow_system_xml_parametrized
    myflow.make_work_directory(myflow.work_directory)
    seekrflow_filename = os.path.join(myflow.work_directory, "seekrflow.json")
    myflow.save(seekrflow_filename)
    myflow2 = structures.load_seekrflow(seekrflow_filename)

def test_seekrflow_creation_with_defaults():
    """Test that Seekrflow can be created with default values"""
    flow = structures.Seekrflow()
    assert flow.name == "my_name"
    assert flow.structure_version == "1.0"
    assert flow.workflow_type == "protein_ligand"
    assert flow.work_directory == "work"
    assert flow.basename_output == "complex"
    assert flow.ionic_strength == 0.15
    assert len(flow.ligand_indices) == 0

def test_make_work_directory():
    """Test work directory creation"""
    flow = structures.Seekrflow()
    with tempfile.TemporaryDirectory() as tmpdir:
        work_path = pathlib.Path(tmpdir) / "test_work"
        flow.make_work_directory(work_path)
        
        assert flow.work_directory == str(work_path)
        assert os.path.exists(work_path)
        assert os.path.isdir(work_path)

def test_get_directory_methods():
    """Test the various get_*_directory methods"""
    flow = structures.Seekrflow()
    with tempfile.TemporaryDirectory() as tmpdir:
        work_path = pathlib.Path(tmpdir) / "test_work"
        flow.make_work_directory(work_path)
        
        # Test get_work_directory
        work_dir = flow.get_work_directory()
        assert work_dir == work_path
        
        # Test get_parametrize_directory
        param_dir = flow.get_parametrize_directory()
        assert param_dir == work_path / "parametrize"
        assert os.path.exists(param_dir)  # Should be created
        
        # Test get_root_directory
        root_dir = flow.get_root_directory()
        assert root_dir == work_path / "root"
        assert os.path.exists(root_dir)  # Should be created
        
        # Test get_run_directory
        run_dir = flow.get_run_directory()
        assert run_dir == work_path / "run"
        assert os.path.exists(run_dir)  # Should be created

def test_save_load_roundtrip():
    """Test that save/load preserves all data correctly"""
    # Create a Seekrflow with non-default values
    flow = structures.Seekrflow()
    flow.name = "test_flow"
    flow.receptor_ligand_pdb = "test.pdb"
    flow.ligand_resname = "LIG"
    flow.ligand_indices = [1, 2, 3, 4]
    flow.ionic_strength = 0.2
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Save
        save_path = os.path.join(tmpdir, "test_flow.json")
        flow.save(save_path)
        
        # Verify file was created and is valid JSON
        assert os.path.exists(save_path)
        with open(save_path, 'r') as f:
            data = json.load(f)
        
        # Load and compare
        loaded_flow = structures.load_seekrflow(save_path)
        assert loaded_flow.name == "test_flow"
        assert loaded_flow.receptor_ligand_pdb == "test.pdb"
        assert loaded_flow.ligand_resname == "LIG"
        assert loaded_flow.ligand_indices == [1, 2, 3, 4]
        assert loaded_flow.ionic_strength == 0.2

def test_load_seekrflow_invalid_file():
    """Test load_seekrflow with invalid file"""
    with pytest.raises(FileNotFoundError):
        structures.load_seekrflow("nonexistent_file.json")

def test_load_seekrflow_invalid_json():
    """Test load_seekrflow with invalid JSON"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        f.write("invalid json content {{{")
        f.flush()
        
        with pytest.raises(json.JSONDecodeError):
            structures.load_seekrflow(f.name)
        
        os.unlink(f.name)

def test_save_seekrflow_invalid_path():
    """Test save with invalid directory path"""
    flow = structures.Seekrflow()
    with pytest.raises(FileNotFoundError):
        flow.save("/nonexistent/directory/file.json")