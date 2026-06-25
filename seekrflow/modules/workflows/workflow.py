"""
modules/workflows/workflow.py

The single, generic seekrflow ``Workflow``. It replaces the old per-system
workflow classes (protein-ligand, protein-protein, protein-ligand-membrane),
which had heavy duplicated boilerplate. A workflow is assembled from reusable,
orthogonal parts:

  * ``components``     -- named pieces of the system + derived atom groups,
  * ``cv_specs``       -- collective-variable specifications,
  * ``anchor_spec``    -- where the milestoning anchors sit (or ``None``),
  * ``procedure``      -- the ordered stage procedure to run,
  * ``scale_settings`` -- per-scale (MD / BD / ...) settings.

This class is DATA ONLY. The logic that turns a Workflow into a seekr3
``Seekr_input`` lives in ``modules/workflows/prepare.py``.
"""

import typing

from attrs import define, field, validators, Factory

import seekrflow.modules.workflows.components as components_module
import seekrflow.modules.workflows.cv_specs as cv_specs_module
import seekrflow.modules.workflows.anchor_specs as anchor_specs_module
import seekrflow.modules.workflows.scale_settings as scale_settings_module
import seekrflow.modules.workflows.stage_procedures as stage_procedures_module


@define
class Workflow:
    """
    A generic, composable seekrflow workflow.
    """
    type: typing.Literal["workflow"] = "workflow"
    # NOTE: we are assuming that the old workflow approach was structure
    # version 1.0, in case we need any backwards compatibility.
    structure_version: str = field(default="1.1",
                                 validator=validators.instance_of(str))
    components: components_module.Components = field(
        default=Factory(components_module.Components),
        validator=validators.instance_of(components_module.Components),
        )
    cv_specs: list[cv_specs_module.CV_spec_base] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(
                cv_specs_module.CV_spec_base),
            iterable_validator=validators.instance_of(list),
        ))
    anchor_spec: anchor_specs_module.Anchor_spec_base | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(
            anchor_specs_module.Anchor_spec_base)),
        )
    procedure: stage_procedures_module.Stage_procedure_base \
        = field(default=Factory(stage_procedures_module\
            .Composite_stage_procedure),
        validator=validators.instance_of(
            stage_procedures_module.Stage_procedure_base),
        )
    scale_settings: list[scale_settings_module.Scale_settings_base] = field(
        default=Factory(
            lambda: [scale_settings_module.Molecular_dynamics_scale_settings()]),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(
                scale_settings_module.Scale_settings_base),
            iterable_validator=validators.instance_of(list),
        ))
    plugins: list[str] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(str),
            iterable_validator=validators.instance_of(list),
        ))

    def get_scale_settings(
            self,
            scale_type: str,
            ) -> scale_settings_module.Scale_settings_base:
        """
        Return the settings object for the given scale type
        (e.g. ``"molecular_dynamics"`` or ``"brownian_dynamics"``).
        """
        for settings in self.scale_settings:
            if settings.get_scale_type() == scale_type:
                return settings
        raise KeyError(f"No scale settings for scale type {scale_type}.")

    def has_scale(self, scale_type: str) -> bool:
        return any(
            settings.get_scale_type() == scale_type
            for settings in self.scale_settings)

    def has_bd(self) -> bool:
        return self.has_scale("brownian_dynamics")

    def get_md_settings(
            self,
            ) -> scale_settings_module.Molecular_dynamics_scale_settings:
        return self.get_scale_settings("molecular_dynamics")

    def has_small_molecule(self) -> bool:
        return any(
            isinstance(member,
                       components_module.Small_molecule_component)
            for member in self.components.members)
