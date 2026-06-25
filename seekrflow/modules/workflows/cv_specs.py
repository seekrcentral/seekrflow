"""
modules/workflows/cv_specs.py

Collective-variable specifications. A CV_spec is a seekrflow-native, JSON-
serializable description of one or more collective variables. It references the
atom groups it needs by NAME (those names are produced by ``components.py``).

``modules/workflows/prepare.py`` converts each CV_spec into the corresponding
seekr3 CV input object(s).
"""

import typing

from attrs import define, field, validators


@define
class CV_spec_base:
    """
    Base class for collective-variable specifications.
    """
    type: typing.Literal["base"] = "base"
    name: str = field(
        default="cv",
        validator=validators.instance_of(str),
        )


@define
class Com_com_distance_CV_spec(CV_spec_base):
    """
    Distance between the centers of mass of two named atom groups.
    Maps to a single seekr3 ``Distance_cv_input``.
    """
    type: typing.Literal["com_com_distance"] = "com_com_distance"
    group1_selection_name: str = field(
        default="receptor_site",
        validator=validators.instance_of(str),
        )
    group2_selection_name: str = field(
        default="ligand",
        validator=validators.instance_of(str),
        )
    min_value: float = field(
        default=0.0,
        validator=validators.instance_of(float),
        )
    max_value: float = field(
        default=3.0,
        validator=validators.instance_of(float),
        )


@define
class Com_com_distance_rmsd_CV_spec(CV_spec_base):
    """
    A center-of-mass distance collective variable plus an RMSD collective
    variable. Maps to a seekr3 ``Distance_cv_input`` AND an ``RMSD_cv_input``.
    """
    type: typing.Literal["com_com_distance_rmsd"] = "com_com_distance_rmsd"
    group1_selection_name: str = field(
        default="receptor_site",
        validator=validators.instance_of(str),
        )
    group2_selection_name: str = field(
        default="ligand",
        validator=validators.instance_of(str),
        )
    min_value: float = field(
        default=0.0,
        validator=validators.instance_of(float),
        )
    max_value: float = field(
        default=3.0,
        validator=validators.instance_of(float),
        )
    rmsd_group_selection_name: str = field(
        default="ligand",
        validator=validators.instance_of(str),
        )
    rmsd_align_selection_name: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
        )
    rmsd_ref_structure: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    rmsd_min_value: float = field(
        default=0.0,
        validator=validators.instance_of(float),
        )
    rmsd_max_value: float = field(
        default=1.0,
        validator=validators.instance_of(float),
        )


CV_spec = Com_com_distance_CV_spec | Com_com_distance_rmsd_CV_spec
