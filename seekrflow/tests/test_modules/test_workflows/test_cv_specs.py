"""
Tests for seekrflow/modules/workflows/cv_specs.py
"""

import pytest

import seekrflow.modules.workflows.cv_specs as cv_specs


def test_resolve_cv_names_single_cv_uses_spec_name():
    specs = [cv_specs.Com_com_distance_CV_spec(name="dist")]
    assert cv_specs.resolve_cv_names(["dist"], specs) == ["dist"]
    assert cv_specs.resolve_cv_names(["dist.distance"], specs) == ["dist"]


def test_resolve_cv_names_combo_expands_roles_in_request_order():
    specs = [cv_specs.Com_com_distance_rmsd_CV_spec(name="combo")]
    assert cv_specs.resolve_cv_names(["combo"], specs) == [
        "combo_distance", "combo_rmsd"]
    assert cv_specs.resolve_cv_names(["combo.rmsd"], specs) == ["combo_rmsd"]
    assert cv_specs.resolve_cv_names(
        ["combo.rmsd", "combo.distance"], specs) == [
        "combo_rmsd", "combo_distance"]


def test_resolve_cv_names_rejects_bad_references():
    specs = [cv_specs.Com_com_distance_rmsd_CV_spec(name="combo")]
    with pytest.raises(ValueError, match="Empty CV name"):
        cv_specs.resolve_cv_names([""], specs)
    with pytest.raises(ValueError, match="Duplicate CV_spec name"):
        cv_specs.resolve_cv_names(
            ["combo"],
            [specs[0], cv_specs.Com_com_distance_CV_spec(name="combo")])
    with pytest.raises(ValueError, match="not defined"):
        cv_specs.resolve_cv_names(["missing"], specs)
    with pytest.raises(ValueError, match="not defined"):
        cv_specs.resolve_cv_names(["missing.distance"], specs)
    with pytest.raises(ValueError, match="no role"):
        cv_specs.resolve_cv_names(["combo.foo"], specs)
    with pytest.raises(ValueError, match="Invalid CV reference"):
        cv_specs.resolve_cv_names(["combo.a.b"], specs)
