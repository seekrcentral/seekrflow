"""
modules/workflows/anchor_specs.py

Anchor (milestoning) specifications. An anchor_spec is a seekrflow-native,
JSON-serializable description of how the Voronoi anchors are placed along a
single collective variable. Placement is computed automatically from the
collective variable's value range (and, for graded placement, the receptor
surface inferred from the starting structure).

``modules/workflows/prepare.py`` converts an anchor_spec into a seekr3
``Anchors_as_voronoi_cells_anchor_scheme`` (with one ``Voronoi_cell`` per
anchor value). State points (the bound state) are derived separately by
``prepare.py`` from the starting structure; the bulk state is simply the
outermost anchor.
"""

import typing

from attrs import define, field, validators


@define
class Anchor_spec_base:
    """
    Base class for anchor specifications. Subclasses place anchors along a
    single collective variable, given that variable's value range.
    """
    type: typing.Literal["base"] = "base"

    def get_anchor_value_vectors(
            self,
            min_value: float,
            max_value: float,
            surface_value: float | None = None,
            ) -> list[list[float]]:
        """
        Return the anchor values as a list of one-component vectors (one
        component per collective variable; a single collective variable is
        currently supported).

        Parameters
        ----------
        min_value, max_value :
            The collective variable's value range, in the collective
            variable's own units (e.g. nanometers for a distance).
        surface_value :
            For graded placement, the collective-variable value at which the
            spacing transitions from fine to coarse (e.g. the receptor surface
            inferred from the starting structure). Ignored by uniform
            placement.
        """
        raise NotImplementedError


@define
class Uniform_anchor_spec(Anchor_spec_base):
    """
    Evenly-spaced anchors along a single collective variable. ``n_anchors``
    anchors are placed at the centers of ``n_anchors`` equal-width cells
    spanning ``[min_value, max_value]``.
    """
    type: typing.Literal["uniform"] = "uniform"
    n_anchors: int = field(
        default=10,
        validator=[validators.instance_of(int), validators.gt(0)],
        )

    def get_anchor_value_vectors(
            self,
            min_value: float,
            max_value: float,
            surface_value: float | None = None,
            ) -> list[list[float]]:
        width = (max_value - min_value) / self.n_anchors
        return [
            [round(min_value + (i + 0.5) * width, 6)]
            for i in range(self.n_anchors)]


@define
class Graded_anchor_spec(Anchor_spec_base):
    """
    Anchors that are finely spaced near the receptor and more widely spaced
    out in the solvent. The fine-to-coarse transition occurs at the receptor
    surface: the closest approach of solvent to the receptor center (computed
    by ``prepare.py`` from the starting structure) plus ``surface_margin``.

    All spacings and distances are given in NANOMETERS.
    """
    type: typing.Literal["graded"] = "graded"
    inner_spacing: float = field(
        default=0.075,
        validator=[validators.instance_of(float), validators.gt(0.0)],
        )
    outer_spacing: float = field(
        default=0.15,
        validator=[validators.instance_of(float), validators.gt(0.0)],
        )
    surface_margin: float = field(
        default=0.3,
        validator=validators.instance_of(float),
        )
    #max_distance: float = field(
    #    default=1.5,
    #    validator=[validators.instance_of(float), validators.gt(0.0)],
    #    )

    def get_anchor_value_vectors(
            self,
            min_value: float,
            max_value: float,
            surface_value: float | None = None,
            ) -> list[list[float]]:
        if surface_value is None:
            raise ValueError(
                "Graded_anchor_spec requires a surface crossover value "
                "(derived from the starting structure).")
        inner = self.inner_spacing
        outer = self.outer_spacing
        surface = surface_value + self.surface_margin
        #max_distance = self.max_distance
        max_distance = max_value
        positions: list[float] = []
        value = min_value + inner / 2.0
        while value < surface and value <= max_distance:
            positions.append(round(value, 6))
            value += inner
        while value <= max_distance:
            positions.append(round(value, 6))
            value += outer
        return [[position] for position in positions]


Anchor_spec = Uniform_anchor_spec | Graded_anchor_spec
