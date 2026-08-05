"""
modules/workflows/components.py

A "component" is a named, chemically-meaningful piece of the molecular system
(a small-molecule ligand, a receptor protein, a protein partner, a membrane,
...). Components own the logic that turns a molecular structure into concrete
atom selections (lists of indices). They also own optional PQR emission used by
the Brownian-dynamics scale.

Selections produced here feed three downstream consumers (all by NAME):
  * collective variables (CV_specs) reference selection names,
  * restraints reference selection names,
  * seekr3 engine settings carry the resolved Selection list.

This module deliberately contains NO seekr3 imports; it produces plain
``{name: [atom_index, ...]}`` mappings. ``modules/workflows/prepare.py`` wraps
those into seekr3 ``Selection`` objects.
"""

import typing

from attrs import define, field, validators, Factory

import mdtraj

import seekrflow.modules.cvs as cvs
import seekrflow.modules.site_finder as site_finder


# ======================================================================
#  Selectors: turn a structure into a list of atom indices.
# ======================================================================

@define
class Selector_base:
    """
    Base class for atom selectors.
    """
    type: typing.Literal["base"] = "base"
    resolved_indices: list[int] = field(
        init=False,
        default=Factory(list))

    def resolve(self, traj: mdtraj.Trajectory) -> list[int]:
        raise NotImplementedError


@define
class Selector_by_indices(Selector_base):
    """
    Select atoms by an explicit list of (zero-based) atom indices.
    """
    type: typing.Literal["indices"] = "indices"
    indices: list[int] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(int),
            iterable_validator=validators.instance_of(list),
        ))
    resolved_indices: list[int] = field(
        init=False,
        default=Factory(list))

    def resolve(self, traj: mdtraj.Trajectory) -> list[int]:
        self.resolved_indices = [int(i) for i in self.indices]
        return self.resolved_indices


@define
class Selector_by_resname(Selector_base):
    """
    Select all atoms belonging to a given residue name (e.g. a ligand).
    """
    type: typing.Literal["resname"] = "resname"
    resname: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    resolved_indices: list[int] = field(
        init=False,
        default=Factory(list))
    include_hydrogens: bool = field(
        default=True,
        validator=validators.instance_of(bool),
        )

    def resolve(self, traj: mdtraj.Trajectory) -> list[int]:
        assert self.resname != "", "resname must be defined."
        if self.include_hydrogens:
            indices = traj.topology.select(f"resname '{self.resname}'")
        else:
            indices = traj.topology.select(f"resname '{self.resname}' and not element H")
        self.resolved_indices = [int(i) for i in indices]
        return self.resolved_indices


@define
class Selector_mdtraj(Selector_base):
    """
    Select atoms with an arbitrary MDTraj atom-selection string.
    """
    type: typing.Literal["mdtraj"] = "mdtraj"
    selection_string: str = field(
        default="protein",
        validator=validators.instance_of(str),
        )
    resolved_indices: list[int] = field(
        init=False,
        default=Factory(list))

    def resolve(self, traj: mdtraj.Trajectory) -> list[int]:
        indices = traj.topology.select(self.selection_string)
        self.resolved_indices = [int(i) for i in indices]
        return self.resolved_indices


@define
class Selector_by_chain(Selector_base):
    """
    Select all atoms belonging to a given (zero-based) chain index.
    """
    type: typing.Literal["chain"] = "chain"
    chain_index: int = field(
        default=0,
        validator=validators.instance_of(int),
        )
    resolved_indices: list[int] = field(
        init=False,
        default=Factory(list))

    def resolve(self, traj: mdtraj.Trajectory) -> list[int]:
        indices = traj.topology.select(f"chainid {self.chain_index}")
        self.resolved_indices = [int(i) for i in indices]
        return self.resolved_indices


Selector = (
    Selector_by_indices | Selector_by_resname | Selector_mdtraj
    | Selector_by_chain
)


# ======================================================================
#  Components: named pieces of the system.
# ======================================================================

@define
class Component_base:
    """
    Base class for a named component of the molecular system.
    """
    type: typing.Literal["base"] = "base"
    name: str = field(
        default="component",
        validator=validators.instance_of(str),
        )
    selector: Selector = field(
        default=Factory(Selector_mdtraj),
        validator=validators.instance_of(Selector_base),
        )
    resolved_indices: list[int] = field(
        init=False,
        default=Factory(list))
    # TODO: this one is also extraneous because this file is already
    # optionally supplied by BD_molecule.
    #pqr_filename_for_bd: str | None = field(
    #    default=None,
    #    validator=validators.optional(validators.instance_of(str)),
    #    )

    @property
    def one_resid_per_atom_pqr(self) -> bool:
        """
        Whether this component's Brownian-dynamics PQR should place each atom
        in its own residue. Determined by the component's chemistry (only small
        molecules need this), so it is derived rather than user-set.
        """
        return False

    def resolve(self, traj: mdtraj.Trajectory) -> list[int]:
        """
        Resolve this component's atom indices and cache them.
        """
        self.resolved_indices = self.selector.resolve(traj)
        return self.resolved_indices
        
@define
class Protein_component(Component_base):
    """
    A protein, whether it acts as the receptor or as a binding partner. Uses
    standard biomolecular force-field parameterization and a single PQR for
    Brownian dynamics.
    """
    type: typing.Literal["protein"] = "protein"


@define
class Small_molecule_component(Component_base):
    """
    A small molecule, e.g. a drug-like ligand or a small-molecule host such as
    a cyclodextrin. Needs special parameterization (GAFF/OpenFF/Espaloma) and a
    per-atom-resid PQR for Brownian dynamics.

    Optional parameterize inputs (ignored when Seekrflow.parameterizer is None):
      * ``sdf_file`` — path to an SDF used for OpenFF/GAFF/Espaloma. If empty,
        parameterize may convert the split PDB via
        ``Parameterizer.pdb_to_sdf_settings`` (RDKit by default). If you already
        have an SDF, set this path and optionally set
        ``pdb_to_sdf_settings`` to ``None``.
      * ``resname`` — residue name to assign when building the OpenFF molecule;
        if empty, falls back to the selector's resname (when applicable) or ``LIG``.
    """
    type: typing.Literal["small_molecule"] = "small_molecule"
    sdf_file: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    resname: str = field(
        default="",
        validator=validators.instance_of(str),
        )

    @property
    def one_resid_per_atom_pqr(self) -> bool:
        return True

    def parameterize_resname(self) -> str:
        """
        Residue name to use when building OpenFF/topology for parameterization.
        """
        if self.resname:
            return self.resname
        if isinstance(self.selector, Selector_by_resname) and self.selector.resname:
            return self.selector.resname
        return "LIG"


@define
class Membrane_component(Component_base):
    """
    A lipid membrane.
    """
    type: typing.Literal["membrane"] = "membrane"


# ======================================================================
#  Group selectors: derive CV/restraint groups from components.
# ======================================================================

@define
class Group_selector_base:
    """
    Base class for a named group selector. Group selectors run AFTER components
    are resolved and may reference resolved component indices.
    """
    type: typing.Literal["base"] = "base"
    name: str = field(
        default="group",
        validator=validators.instance_of(str),
        )

    def resolve(
            self,
            traj: mdtraj.Trajectory,
            member_indices: dict[str, list[int]],
            ) -> list[int]:
        raise NotImplementedError


@define
class Alpha_carbon_site_group_selector(Group_selector_base):
    """
    Select the alpha carbons of receptor residues within ``threshold`` (nm) of a
    reference component (typically the ligand). Defines a binding-site group for
    a center-of-mass collective variable.
    """
    type: typing.Literal["alpha_carbon_site"] = "alpha_carbon_site"
    reference_component_name: str = field(
        default="ligand",
        validator=validators.instance_of(str),
        )
    threshold: float = field(
        default=0.6,
        validator=validators.instance_of(float),
        )
    receptor_selection: str = field(
        default="protein",
        validator=validators.instance_of(str),
        )

    def resolve(
            self,
            traj: mdtraj.Trajectory,
            member_indices: dict[str, list[int]],
            ) -> list[int]:
        assert self.reference_component_name in member_indices, \
            f"Reference component {self.reference_component_name} not found."
        reference_indices = member_indices[self.reference_component_name]
        return cvs.alpha_carbon_selection_within_cutoff(
            traj,
            reference_indices,
            self.threshold,
            receptor_selection=self.receptor_selection,
        )


@define
class Component_atoms_group_selector(Group_selector_base):
    """
    A group that is simply (a copy of) a component's atoms.
    """
    type: typing.Literal["component_atoms"] = "component_atoms"
    component_name: str = field(
        default="ligand",
        validator=validators.instance_of(str),
        )

    def resolve(
            self,
            traj: mdtraj.Trajectory,
            member_indices: dict[str, list[int]],
            ) -> list[int]:
        assert self.component_name in member_indices, \
            f"Component {self.component_name} not found."
        return list(member_indices[self.component_name])


@define
class Binding_site_group_selector(Group_selector_base):
    """
    Select a protein binding-site group via the Monte-Carlo site finder
    (``modules/site_finder.py``). Starting from the receptor atoms near a
    reference component (typically the ligand), it searches for a compact set
    of alpha carbons whose center of mass sits in the pocket, giving a stable
    group for a center-of-mass collective variable.

    Unlike :class:`Alpha_carbon_site_group_selector` (a plain distance cutoff),
    this performs a stochastic optimization, so a ``random_seed`` is used for
    reproducibility and the result is cached in ``resolved_indices`` (frozen
    into the saved input) the first time it is computed. Only meaningful when
    the receptor is a protein.
    """
    type: typing.Literal["binding_site"] = "binding_site"
    reference_component_name: str = field(
        default="ligand",
        validator=validators.instance_of(str),
        )
    inner_cutoff: float = field(
        default=0.5,
        validator=validators.instance_of(float),
        )
    outer_cutoff: float = field(
        default=0.8,
        validator=validators.instance_of(float),
        )
    random_seed: int= field(
        default=0, validator=validators.instance_of(int),
        )
    resolved_indices: list[int] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(int),
            iterable_validator=validators.instance_of(list),
        ))

    def resolve(
            self,
            traj: mdtraj.Trajectory,
            member_indices: dict[str, list[int]],
            ) -> list[int]:
        # The Monte-Carlo search is expensive and stochastic, so run it once
        #  and reuse the frozen result on subsequent resolves.
        if self.resolved_indices:
            return list(self.resolved_indices)
        assert self.reference_component_name in member_indices, \
            f"Reference component {self.reference_component_name} not found."
        ligand_indices = member_indices[self.reference_component_name]
        self.resolved_indices = site_finder.site_finder_monte_carlo_from_traj(
            traj,
            ligand_indices,
            inner_cutoff=self.inner_cutoff,
            outer_cutoff=self.outer_cutoff,
            random_seed=self.random_seed,
        )
        return list(self.resolved_indices)


Group_selector = (
    Alpha_carbon_site_group_selector | Component_atoms_group_selector
    | Binding_site_group_selector
)


# ======================================================================
#  Components container.
# ======================================================================

@define
class Components:
    """
    The full collection of components and derived group selectors for a
    Workflow.
    """
    members: list[Component_base] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(Component_base),
            iterable_validator=validators.instance_of(list),
        ))
    group_selectors: list[Group_selector_base] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(Group_selector_base),
            iterable_validator=validators.instance_of(list),
        ))

    def get_member(self, name: str) -> Component_base:
        for member in self.members:
            if member.name == name:
                return member
        raise KeyError(f"No component named {name}.")

    def bd_molecule_md_atoms(
            self, bd_molecule: typing.Any) -> list[int]:
        """
        Return the molecular-dynamics (topology) atom indices that make up a
        single Brownian-dynamics rigid body, as the ascending-sorted union of
        the molecule's referenced components' resolved indices.

        Components must already be resolved (this is guaranteed when called
        after :meth:`resolve_all_selections`).
        """
        atoms: list[int] = []
        member = self.get_member(bd_molecule.component_name)
        atoms.extend(member.resolved_indices)
        return sorted(set(atoms))

    def resolve_all_selections(
            self,
            scale_settings: list,
            ) -> dict[str, dict[str, list[int]]]:
        """
        Resolve selections for every scale. Returns a nested mapping::

            { scale_type: { selection_name: [atom_index, ...], ... }, ... }

        The molecular-dynamics scale resolves components and group selectors
        against its solvated structure, in molecular-dynamics (topology) atom
        space. The Brownian-dynamics scale, if present, remaps each selection
        onto the per-molecule, PQR-local index space of the rigid body that
        contains it.
        """
        md_settings = None
        for settings in scale_settings:
            if settings.get_scale_type() == "molecular_dynamics":
                md_settings = settings
                break
        assert md_settings is not None, \
            "A molecular_dynamics scale is required to resolve selections."
        assert md_settings.system is not None, \
            "molecular_dynamics.system must be set (or filled from "\
            "parameterize outputs) before resolving selections."
        assert md_settings.system.solvated_pdb, \
            "molecular_dynamics.system.solvated_pdb must be set."
        result: dict[str, dict[str, list[int]]] = {}
        md_selections = self._resolve_md_selections(
            md_settings.system.solvated_pdb)
        result["molecular_dynamics"] = md_selections
        for settings in scale_settings:
            if settings.get_scale_type() == "brownian_dynamics":
                result["brownian_dynamics"] = self._resolve_bd_selections(
                    settings, md_selections)
        return result

    def _resolve_md_selections(
            self,
            pdb_filename: str,
            ) -> dict[str, list[int]]:
        """
        Resolve every component and group selector against the given structure.
        Returns an ordered ``{name: [atom_index, ...]}`` mapping (members first,
        then group selectors), in molecular-dynamics atom space.
        """
        traj = mdtraj.load(pdb_filename)
        selections: dict[str, list[int]] = {}
        for member in self.members:
            selections[member.name] = member.resolve(traj)
            assert len(selections[member.name]) > 0, \
                f"Member '{member.name}' produced an empty selection."
        for group_selector in self.group_selectors:
            selections[group_selector.name] = group_selector.resolve(
                traj, selections)
            assert len(selections[group_selector.name]) > 0, \
                f"Group selector '{group_selector.name}' produced an empty selection."
        return selections

    def _resolve_bd_selections(
            self,
            bd_settings: typing.Any,
            md_selections: dict[str, list[int]],
            ) -> dict[str, dict[str, list[int]]]:
        """
        Remap molecular-dynamics selections onto each Brownian-dynamics
        molecule's PQR-local index space. Returns
        ``{ molecule_name: { selection_name: [local_index, ...], ... }, ... }``.

        A selection is attached to a molecule only if all of its atoms are
        contained in that molecule. Selections that fall entirely outside the
        molecule, or that span more than one molecule (e.g. the synthesized
        molecular-dynamics ``solute`` group covering several rigid bodies), are
        simply omitted from this molecule's Brownian-dynamics selections.
        """
        bd_result: dict[str, dict[str, list[int]]] = {}
        for molecule in bd_settings.system.molecules:
            molecule_atoms = self.bd_molecule_md_atoms(molecule)
            molecule_atom_set = set(molecule_atoms)
            atom_position = {atom: i for i, atom in enumerate(molecule_atoms)}
            local_selections: dict[str, list[int]] = {}
            for name, indices in md_selections.items():
                if not indices:
                    continue
                index_set = set(indices)
                if index_set <= molecule_atom_set:
                    local_selections[name] = [
                        atom_position[atom] for atom in indices]
            bd_result[molecule.name] = local_selections
        return bd_result
