"""Tests for CV_spec generated names and spec.role resolution."""

import pytest

import seekrflow.modules.workflows.cv_specs as cv_specs
import seekrflow.modules.workflows.prepare as prepare
import seekrflow.modules.workflows.stage_procedures as stage_procedures


def test_com_com_distance_generated_cvs():
    spec = cv_specs.Com_com_distance_CV_spec(name="dist")
    assert spec.generated_cvs() == [("distance", "dist")]


def test_com_com_distance_rmsd_generated_cvs():
    spec = cv_specs.Com_com_distance_rmsd_CV_spec(name="combo")
    assert spec.generated_cvs() == [
        ("distance", "combo_distance"),
        ("rmsd", "combo_rmsd"),
    ]


def test_resolve_cv_names_expands_multi_cv_spec():
    specs = [cv_specs.Com_com_distance_rmsd_CV_spec(name="combo")]
    assert cv_specs.resolve_cv_names(["combo"], specs) == [
        "combo_distance", "combo_rmsd"]


def test_resolve_cv_names_spec_role_subset():
    specs = [cv_specs.Com_com_distance_rmsd_CV_spec(name="combo")]
    assert cv_specs.resolve_cv_names(["combo.rmsd"], specs) == ["combo_rmsd"]
    assert cv_specs.resolve_cv_names(
        ["combo.distance", "combo.rmsd"], specs) == [
        "combo_distance", "combo_rmsd"]


def test_resolve_cv_names_single_cv_identity():
    specs = [cv_specs.Com_com_distance_CV_spec(name="dist")]
    assert cv_specs.resolve_cv_names(["dist"], specs) == ["dist"]
    assert cv_specs.resolve_cv_names(["dist.distance"], specs) == ["dist"]


def test_resolve_cv_names_unknown_spec_or_role():
    specs = [cv_specs.Com_com_distance_rmsd_CV_spec(name="combo")]
    with pytest.raises(ValueError, match="not defined"):
        cv_specs.resolve_cv_names(["missing"], specs)
    with pytest.raises(ValueError, match="no role"):
        cv_specs.resolve_cv_names(["combo.foo"], specs)
    with pytest.raises(ValueError, match="Invalid CV reference"):
        cv_specs.resolve_cv_names(["combo.a.b"], specs)


def test_map_sampling_metadynamics_resolves_cv_names():
    specs = [cv_specs.Com_com_distance_rmsd_CV_spec(name="combo")]
    sampling = stage_procedures.Metadynamics_sampling_spec(
        cv_names=["combo"])
    mapped = prepare._map_sampling(sampling, specs)
    assert mapped.cv_names == ["combo_distance", "combo_rmsd"]

    sampling_subset = stage_procedures.Metadynamics_sampling_spec(
        cv_names=["combo.distance"])
    mapped_subset = prepare._map_sampling(sampling_subset, specs)
    assert mapped_subset.cv_names == ["combo_distance"]
