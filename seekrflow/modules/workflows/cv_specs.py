"""
modules/workflows/cv_specs.py

Collective-variable specifications. A CV_spec is a seekrflow-native, JSON-
serializable description of one or more collective variables. It references the
atom groups it needs by NAME (those names are produced by ``components.py``).

``modules/workflows/prepare.py`` converts each CV_spec into the corresponding
seekr3 CV input object(s).

Sampling procedures reference CV_specs by name (or ``spec.role`` to select one
generated CV). Use :func:`resolve_cv_names` at prepare time to expand those
references into the concrete seekr CV names produced by :meth:`generated_cvs`.
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

    def generated_cvs(self) -> list[tuple[str, str]]:
        """
        Return ``(role, cv_name)`` pairs for the seekr CVs this spec produces,
        in order. ``role`` is the stable sub-name used in ``spec.role``
        references; ``cv_name`` is the concrete name written into the model.
        """
        raise NotImplementedError


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

    def generated_cvs(self) -> list[tuple[str, str]]:
        return [("distance", self.name)]


@define
class Com_com_distance_rmsd_CV_spec(CV_spec_base):
    """
    A center-of-mass distance collective variable plus an RMSD collective
    variable. Maps to a seekr3 ``Distance_cv_input`` AND an ``RMSD_cv_input``.

    Roles for ``spec.role`` references: ``distance``, ``rmsd``.
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

    def generated_cvs(self) -> list[tuple[str, str]]:
        return [
            ("distance", f"{self.name}_distance"),
            ("rmsd", f"{self.name}_rmsd"),
        ]


CV_spec = Com_com_distance_CV_spec | Com_com_distance_rmsd_CV_spec


def resolve_cv_names(
        requested: list[str],
        cv_specs: list[CV_spec_base],
        ) -> list[str]:
    """
    Expand sampling ``cv_names`` entries into concrete seekr CV names.

    Each entry is either:
      * ``\"spec\"`` — all CVs produced by that CV_spec, or
      * ``\"spec.role\"`` — one generated CV (e.g. ``\"foo.rmsd\"``).
    """
    by_name = {spec.name: spec for spec in cv_specs}
    if len(by_name) != len(cv_specs):
        dupes = [name for name in by_name if
                 sum(1 for s in cv_specs if s.name == name) > 1]
        raise ValueError(
            f"Duplicate CV_spec name(s): {', '.join(sorted(set(dupes)))}.")

    resolved: list[str] = []
    for entry in requested:
        if not entry:
            raise ValueError("Empty CV name reference is not allowed.")
        if "." in entry:
            spec_name, role = entry.split(".", 1)
            if not spec_name or not role or "." in role:
                raise ValueError(
                    f"Invalid CV reference '{entry}'. Expected 'spec' or "
                    f"'spec.role' with a single role segment.")
            if spec_name not in by_name:
                raise ValueError(
                    f"Sampling method references CV_spec '{spec_name}' that "
                    f"is not defined in the workflow's CV specs.")
            role_map = dict(by_name[spec_name].generated_cvs())
            if role not in role_map:
                known = ", ".join(sorted(role_map))
                raise ValueError(
                    f"CV_spec '{spec_name}' has no role '{role}'. "
                    f"Known roles: {known}.")
            resolved.append(role_map[role])
        else:
            if entry not in by_name:
                raise ValueError(
                    f"Sampling method references CV_spec '{entry}' that is "
                    f"not defined in the workflow's CV specs.")
            resolved.extend(name for _, name in by_name[entry].generated_cvs())
    return resolved
