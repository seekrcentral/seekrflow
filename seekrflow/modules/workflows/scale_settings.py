"""
modules/workflows/scale_settings.py

Per-scale (molecular dynamics, Brownian dynamics, ...) settings for a
generic seekrflow Workflow. These are pure data classes; the logic that
maps them onto seekr3 engine-settings objects lives in
``modules/workflows/prepare.py``.

A Workflow holds a ``list[Scale_settings_base]`` discriminated by ``.type``.
Use ``Workflow.get_scale_settings(scale_type)`` to retrieve a specific one.
"""

import typing

from attrs import define, field, validators, Factory

import seekrflow.modules.parameters_topology_structures \
    as parameters_topology_structures

# The accepted parameters/topology container types (seekrflow-native; mapped to
#  their seekr3 equivalents at prepare time).
Parameters_topology = (
    parameters_topology_structures.Amber_parameters_topology
    | parameters_topology_structures.Gromacs_parameters_topology
    | parameters_topology_structures.Charmm_parameters_topology
    | parameters_topology_structures.Forcefield_parameters
    | parameters_topology_structures.Openmm_system
)


@define
class System_for_md:
    """
    The solvated molecular system used for the molecular-dynamics scale: the
    parameter/topology definition plus a solvated starting structure (PDB).
    """
    parameters_topology: Parameters_topology = field(
        default=Factory(
            parameters_topology_structures.Amber_parameters_topology),
        validator=validators.instance_of(
            parameters_topology_structures.Amber_parameters_topology
            | parameters_topology_structures.Gromacs_parameters_topology
            | parameters_topology_structures.Charmm_parameters_topology
            | parameters_topology_structures.Forcefield_parameters
            | parameters_topology_structures.Openmm_system),
        )
    solvated_pdb: str = field(
        default="",
        validator=validators.instance_of(str),
        )


@define
class BD_molecule:
    """
    One rigid body for the Brownian-dynamics scale, corresponding to a single
    PQR file.

    A molecule may aggregate more than one component (for example a receptor
    protein together with a membrane, or a receptor that hosts two binding
    sites), so it references components by a list of names rather than a single
    name. Atoms are taken in ascending molecular-dynamics (topology) order
    across all referenced components, which is the same order a sliced or
    auto-generated PQR will follow.

    An empty ``pqr_filename`` means "auto-generate this PQR from the
    molecular-dynamics system at prepare time".
    """
    name: str = field(
        default="molecule",
        validator=validators.instance_of(str),
        )
    # NOTE entire molecule
    component_name: str = field(
        default="NOT ASSIGNED",
        validator=validators.instance_of(str),
        )
    # Not yet needed - the CV_spec handles this automatically.
    #site_component_names: list[str] = field(
    #    default=Factory(list),
    #    validator=validators.deep_iterable(
    #        member_validator=validators.instance_of(str),
    #        iterable_validator=validators.instance_of(list),
    #    ))
    pqr_filename: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    # Optional, lightweight hint (e.g. "receptor" or "ligand") for downstream
    # tooling. The operative receptor/ligand distinction is the ORDER of
    # ``System_for_bd.molecules`` (first molecule is the receptor), not this.
    role: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
        )


@define
class System_for_bd:
    """
    The set of rigid bodies (one PQR each) used by the Brownian-dynamics scale.
    """
    molecules: list[BD_molecule] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(BD_molecule),
            iterable_validator=validators.instance_of(list),
        ))


@define
class Scale_settings_base:
    """
    Base class for per-scale settings.
    """
    type: typing.Literal["base"] = "base"
    #scale_type: str = field(
    #    default="base",
    #    validator=validators.instance_of(str),
    #    )
    def get_scale_type(self) -> str:
        raise NotImplementedError(
            "get_scale_type() must be implemented by subclasses of "\
            "Scale_settings_base.")


@define
class Molecular_dynamics_scale_settings(Scale_settings_base):
    """
    OpenMM molecular-dynamics scale settings. Maps to a seekr3
    ``Openmm_settings_input``.
    """
    type: typing.Literal["molecular_dynamics"] = "molecular_dynamics"
    system: System_for_md = field(
        default=Factory(System_for_md),
        validator=validators.instance_of(System_for_md),
        )
    nonbonded_method: str = field(
        default="pme",
        validator=validators.instance_of(str),
        )
    nonbonded_cutoff: float = field(
        default=0.9,
        validator=validators.instance_of(float),
        )
    constraints: str = field(
        default="hbonds",
        validator=validators.instance_of(str),
        )
    hydrogen_mass: float | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(float)),
        )
    timestep: float = field(
        default=0.002,
        validator=validators.instance_of(float),
        )
    friction_coefficient: float | None = field(
        default=1.0,
        validator=validators.optional(validators.instance_of(float)),
        )
    integrator_type: str = field(
        default="langevin",
        validator=validators.instance_of(str),
        )
    platform_type: str = field(
        default="cuda",
        validator=validators.instance_of(str),
        )

    def get_scale_type(self) -> str:
        return "molecular_dynamics"

@define
class Brownian_dynamics_scale_settings(Scale_settings_base):
    """
    Brownian-dynamics scale settings. Maps to a seekr3 BD engine-settings
    object (Browndye2 or SDA).
    """
    type: typing.Literal["brownian_dynamics"] = "brownian_dynamics"
    engine: str = field(
        default="browndye",
        validator=validators.in_(["browndye", "sda"]),
        )
    system: System_for_bd = field(
        default=Factory(System_for_bd),
        validator=validators.instance_of(System_for_bd),
        )
    binary_directory: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    
    def get_scale_type(self) -> str:
        return "brownian_dynamics"
