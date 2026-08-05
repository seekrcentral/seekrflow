"""
modules/structures.py

Contain data structure classes used for seekrflow parameters/inputs.
"""

import os
import json
import glob
import typing
import pathlib
from shutil import copyfile

from attrs import define, field, validators, Factory
import cattrs
from cattrs.strategies import include_subclasses

import seekrflow.modules.base as base
import seekrflow.modules.workflows.workflow as workflow_module
import seekrflow.modules.workflows.components as components_module
import seekrflow.modules.workflows.cv_specs as cv_specs_module
import seekrflow.modules.workflows.anchor_specs as anchor_specs_module
import seekrflow.modules.workflows.scale_settings as scale_settings_module
import seekrflow.modules.workflows.stage_procedures as stage_procedures_module
import seekrflow.modules.parameterize_structures as parameterizer_structures


WORK = "work"
PARAMETERIZE = "parameterize"
ROOT = "root"
RUN = "run"
INPUT_TARBALL_NAME = "model_input.tar.gz" # For Cloud transfers


def _register_workflow_subclasses(converter: cattrs.Converter) -> None:
    """
    Register all polymorphic base classes used by the generic Workflow so that
    cattrs (un)structures them as their concrete subtypes.

    Registration order matters: ``include_subclasses`` eagerly bakes structure
    hooks for each subclass's fields at registration time, so any union that is
    NESTED inside another must be registered BEFORE its container. The order
    below therefore proceeds from the most deeply nested unions outward.
    """
    # Selectors are nested inside components.
    include_subclasses(components_module.Selector_base, converter)
    include_subclasses(components_module.Component_base, converter)
    include_subclasses(components_module.Group_selector_base, converter)
    include_subclasses(cv_specs_module.CV_spec_base, converter)
    include_subclasses(anchor_specs_module.Anchor_spec_base, converter)
    include_subclasses(scale_settings_module.Scale_settings_base, converter)
    # Sampling / restraint / completion specs are nested inside stage items,
    #  which are in turn nested inside stage procedures.
    include_subclasses(stage_procedures_module.Sampling_spec_base, converter)
    include_subclasses(stage_procedures_module.Restraint_spec_base, converter)
    include_subclasses(stage_procedures_module.Completion_spec_base, converter)
    include_subclasses(stage_procedures_module.Stage_item_base, converter)
    include_subclasses(stage_procedures_module.Stage_procedure_base, converter)
    return

@define
class Remote_interface_base:
    """
    Base remote interface settings class.
    """
    type: typing.Literal["base"] = "base"

@define
class Remote_interface_ssh(Remote_interface_base):
    """
    SSH remote interface settings class.
    """
    type: typing.Literal["ssh"] = "ssh"
    hostname: str = field(
        default="",
        validator=validators.instance_of(str),
    )
    username: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
    )
    password: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
    )
    port: int = field(
        default=22,
        validator=validators.optional(validators.instance_of(int)),
    )
    private_key_filename: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
    )
    private_key_passphrase: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
    )

@define
class Remote_interface_globus_compute_sdk(Remote_interface_base):
    """
    Globus Compute SDK remote interface settings class.
    """
    type: typing.Literal["globus_compute_sdk"] = "globus_compute_sdk"
    endpoint_id: str = field(
        default="",
        validator=validators.instance_of(str),
    )

@define
class Transfer_settings_base:
    """
    Base transfer settings class.
    """
    type: typing.Literal["base"] = "base"

@define
class Transfer_settings_globus(Transfer_settings_base):
    """
    Globus transfer settings for transferring files.
    """
    type: typing.Literal["globus"] = "globus"
    local_collection_id: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    remote_collection_id: str = field(
        default="",
        validator=validators.instance_of(str),
        )

@define
class Transfer_settings_rsync(Transfer_settings_base):
    """
    Rsync transfer settings for transferring files.
    """
    type: typing.Literal["rsync"] = "rsync"
    hostname: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    username: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
        )
    password: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
        )
    port: int = field(
        default=22,
        validator=validators.instance_of(int)
        )
    private_key_filename: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
        )
    private_key_passphrase: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
        )

@define
class Transfer_settings_aws_s3(Transfer_settings_base):
    """
    Amazon web services cloud file transfer utility.
    """
    type: typing.Literal["aws_s3"] = "aws_s3"
    bucket: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    prefix: str = field(
        default="",
        validator=validators.instance_of(str),
        )

    def get_uri(self):
        return f"s3://{self.bucket}/{self.prefix}"

    def get_input_tarball_name(self):
        return INPUT_TARBALL_NAME

@define
class Resource_base:
    """
    Base class for resources.
    """
    type: typing.Literal["base"] = "base"
    

@define
class Resource_local(Resource_base):
    """
    Local resource for running the protocol.
    """
    pass

@define
class Resource_remote_base(Resource_base):
    """
    Base class for remote resources.
    """
    pass

@define
class Resource_remote_slurm(Resource_remote_base):
    """
    Slurm resource for running the protocol.
    """
    type: typing.Literal["slurm_remote"] = "slurm_remote"
    name: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    remote_working_directory: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    partition: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    account: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    constraint: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
        )
    nodes_per_block: int = field(
        default=1,
        validator=validators.instance_of(int),
        )
    cpus_per_task: int = field(
        default=1,
        validator=validators.instance_of(int),
        )
    memory_per_node: int = field(
        default=4,
        validator=validators.instance_of(int),
        )
    time_limit: str = field(
        default="00:30:00",
        validator=validators.instance_of(str),
        )
    scheduler_options: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    worker_init: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    mps: int = field(
        default=1,
        validator=validators.instance_of(int),
    )
    remote_interface: Remote_interface_globus_compute_sdk | Remote_interface_ssh = field(
        default=Factory(Remote_interface_globus_compute_sdk),
    )
    transfer_settings: Transfer_settings_rsync | Transfer_settings_globus = field(
        default=Factory(Transfer_settings_globus),
    )

@define
class Resource_remote_pbs(Resource_remote_base):
    """
    PBS/Torque resource for running the protocol.
    """
    type: typing.Literal["pbs_remote"] = "pbs_remote"
    name: str = field(
        default="",
        validator=validators.instance_of(str),
    )
    remote_working_directory: str = field(
        default="",
        validator=validators.instance_of(str),
    )
    queue: str = field(  # PBS uses "queue" instead of "partition"
        default="",
        validator=validators.instance_of(str),
    )
    account: str = field(
        default="",
        validator=validators.instance_of(str),
    )
    resource_list: str | None = field(  # PBS -l options (e.g., "qos=high")
        default=None,
        validator=validators.optional(validators.instance_of(str)),
    )
    nodes_per_block: int = field(
        default=1,
        validator=validators.instance_of(int),
    )
    cpus_per_task: int = field(  # Keep same name as SLURM for consistency; convert to ppn in qsub
        default=1,
        validator=validators.instance_of(int),
    )
    memory_per_node: int = field(
        default=4,
        validator=validators.instance_of(int),
    )
    time_limit: str = field(  # Format: HH:MM:SS (walltime)
        default="00:30:00",
        validator=validators.instance_of(str),
    )
    scheduler_options: str = field(  # Additional PBS options
        default="",
        validator=validators.instance_of(str),
    )
    worker_init: str = field(
        default="",
        validator=validators.instance_of(str),
    )
    mps: int = field(
        default=1,
        validator=[validators.instance_of(int),
                            validators.ge(1)],
    )
    remote_interface: Remote_interface_globus_compute_sdk | Remote_interface_ssh = field(
        default=Factory(Remote_interface_globus_compute_sdk),
    )
    transfer_settings: Transfer_settings_rsync | Transfer_settings_globus = field(
        default=Factory(Transfer_settings_globus),
    )

@define
class Resource_cloud_base(Resource_base):
    """
    Base class for cloud resources.
    """
    pass

@define
class Resource_cloud_aws(Resource_cloud_base):
    """
    AWS resource for running the protocol.
    """
    type: typing.Literal["aws_cloud"] = "aws_cloud"
    name: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    # AWS region: such as 'us-west-2'
    region: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    # Your AWS account id
    account_id: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    # The name of your compute environment
    compute_env_name: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    # The name of your job queue
    job_queue_name: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    # Maximum Batch wall-clock (seconds). Upper bound for this resource;
    # may later be tightened dynamically from benchmarks.
    job_timeout_seconds: int = field(
        default=86400,
        validator=[validators.instance_of(int), validators.gt(0)],
        )

    n_gpus: int = field(
        default=1,
        validator=[validators.instance_of(int),
                   validators.ge(0)]
        )
    n_vcpus: int = field(
        default=4,
        validator=[validators.instance_of(int),
                   validators.ge(0)]
        )
    memory_mb: int = field(
        default=14000,
        validator=[validators.instance_of(int),
                   validators.ge(0)]
        )
    # Max concurrent GPU-sharing processes per Batch task (Placement clamp).
    mps: int = field(
        default=1,
        validator=[validators.instance_of(int), validators.ge(1)],
        )
    transfer_settings: Transfer_settings_aws_s3 = field(
        default=Factory(Transfer_settings_aws_s3),
        )

    def get_seekr_image_uri(self):
        return f"{self.account_id}.dkr.ecr.{self.region}.amazonaws.com/seekr-engines-bd:latest"

@define
class Time_policy_base:
    """
    Base class for Placement walltime policies.
    """
    type: typing.Literal["base"] = "base"


@define
class Time_policy_fixed(Time_policy_base):
    """
    Always request a fixed walltime.

    If ``time_limit`` is unset, resolve uses the Placement ``time_limit``
    override if present, otherwise the Resource cap.
    """
    type: typing.Literal["fixed"] = "fixed"
    time_limit: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)))


@define
class Time_policy_adaptive(Time_policy_base):
    """
    Estimate walltime from remaining work and a performance rate.

    ``estimated_performance`` units depend on stage scale (docstring only):
    molecular dynamics → ns/day; Brownian dynamics → trajectories/day.
    Unset ``min_time_limit`` / ``max_time_limit`` fall back to backend
    defaults (1 h remote / 10 min AWS) and the Resource cap, respectively.
    """
    type: typing.Literal["adaptive"] = "adaptive"
    estimated_performance: float | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(float)))
    safety_factor: float = field(
        default=0.2,
        validator=validators.instance_of(float))
    min_time_limit: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)))
    max_time_limit: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)))


@define
class Placement:
    """
    Deployment policy for stages whose address path begins with ``target``.

    Compute fields (``cpus``, ``memory_mb``, ``time_limit``, ``mps``) are
    backend-agnostic overrides. When unset, ``resolve_stage_execution`` fills
    them from the selected Resource's native defaults.

    Unset ``time_policy`` is treated as adaptive at resolve time.
    """
    target: list[str] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(str),
            iterable_validator=validators.instance_of(list),
        ))
    resource: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)))
    dispatch: stage_procedures_module.Dispatch | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(
            stage_procedures_module.Dispatch)))
    co_schedule_with: str | None = field(
        default=None,
        validator=validators.optional(
            validators.in_(["predecessor", "successor"])))
    cpus: int | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(int)))
    memory_mb: int | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(int)))
    time_limit: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)))
    mps: int | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(int)))
    time_policy: Time_policy_base | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(Time_policy_base)))


@define
class Resolved_execution:
    """Fully-resolved, ready-to-run policy for one stage."""
    stage_name: str
    resource_name: str
    resource: Resource_base | None
    dispatch: stage_procedures_module.Dispatch
    co_schedule_with: str | None
    cpus: int | None
    memory_mb: int | None
    time_limit: str | None  # HH:MM:SS; None when adaptive
    mps: int | None
    time_policy: Time_policy_base = field(
        default=Factory(Time_policy_adaptive),
        validator=validators.instance_of(Time_policy_base))


def time_limit_to_seconds(time_limit: str) -> int:
    """Parse ``HH:MM:SS`` (optional ``D-`` day prefix) to integer seconds."""
    time_str = (time_limit or "").strip()
    if not time_str:
        raise ValueError("time_limit must be a non-empty HH:MM:SS string")
    days = 0
    if "-" in time_str:
        day_part, time_str = time_str.split("-", 1)
        days = int(day_part)
    parts = time_str.split(":")
    if len(parts) != 3:
        raise ValueError(
            f"time_limit must be HH:MM:SS, got {time_limit!r}")
    hours, minutes, seconds = (int(p) for p in parts)
    return days * 86400 + hours * 3600 + minutes * 60 + seconds


def seconds_to_time_limit(seconds: int) -> str:
    """Format non-negative seconds as ``HH:MM:SS`` (hours may exceed 24)."""
    if seconds < 0:
        raise ValueError(f"seconds must be >= 0, got {seconds}")
    hours, rem = divmod(int(seconds), 3600)
    minutes, secs = divmod(rem, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def _resource_compute_defaults(
        resource: Resource_base | None,
        ) -> dict:
    """
    Backend-native Resource fields → agnostic compute defaults for resolve.
    """
    if resource is None:
        return {
            "cpus": None,
            "memory_mb": None,
            "time_limit": None,
            "mps": 1,
        }
    if isinstance(resource, (Resource_remote_slurm, Resource_remote_pbs)):
        return {
            "cpus": resource.cpus_per_task,
            "memory_mb": resource.memory_per_node,
            "time_limit": resource.time_limit,
            "mps": resource.mps,
        }
    if isinstance(resource, Resource_cloud_aws):
        return {
            "cpus": resource.n_vcpus,
            "memory_mb": resource.memory_mb,
            "time_limit": seconds_to_time_limit(resource.job_timeout_seconds),
            "mps": resource.mps,
        }
    return {
        "cpus": None,
        "memory_mb": None,
        "time_limit": None,
        "mps": 1,
    }


def resource_supports_arrays(resource: Resource_base | None) -> bool:
    return isinstance(
        resource,
        (Resource_remote_slurm, Resource_remote_pbs, Resource_cloud_aws),
    )


def _dispatch_uses_array_spread(dispatch: stage_procedures_module.Dispatch) -> bool:
    return bool(dispatch.dimensions)


def _co_schedule_host_name(
        stage_name: str,
        co_schedule_with: str,
        stage_index: int,
        model_stages: list,
        stage_names: list[str],
        ) -> str:
    stage = model_stages[stage_index]
    if co_schedule_with == "predecessor":
        parent_one_based = getattr(stage, "input_stage_index", 0)
        if parent_one_based <= 0:
            raise ValueError(
                f"Stage {stage_name!r} has co_schedule_with='predecessor' "
                f"but has no predecessor in the model chain.")
        return model_stages[parent_one_based - 1].name
    for idx, other in enumerate(model_stages):
        if getattr(other, "input_stage_index", 0) - 1 == stage_index:
            return other.name
    raise ValueError(
        f"Stage {stage_name!r} has co_schedule_with='successor' "
        f"but has no successor in the model chain.")


def _placement_targets_match(address: list[str], target: list[str]) -> bool:
    if len(target) > len(address):
        return False
    return address[:len(target)] == target


def _collect_matching_placements(
        address: list[str],
        placements: list[Placement],
        ) -> list[Placement]:
    seen_targets: set[tuple[str, ...]] = set()
    matching: list[Placement] = []
    for placement in placements:
        target_key = tuple(placement.target)
        if target_key in seen_targets:
            raise ValueError(
                f"Duplicate placement target {list(target_key)!r}.")
        seen_targets.add(target_key)
        if _placement_targets_match(address, placement.target):
            matching.append(placement)
    # Broadest (shortest target) first so more specific placements override.
    # An empty target ([]) is the broadest of all and acts as the default.
    matching.sort(key=lambda p: len(p.target))
    return matching


def _apply_placement_fields(
        accumulator: dict,
        placement: Placement,
        ) -> None:
    if placement.resource is not None:
        accumulator["resource"] = placement.resource
    if placement.co_schedule_with is not None:
        accumulator["co_schedule_with"] = placement.co_schedule_with
    if placement.cpus is not None:
        accumulator["cpus"] = placement.cpus
    if placement.memory_mb is not None:
        accumulator["memory_mb"] = placement.memory_mb
    if placement.time_limit is not None:
        accumulator["time_limit"] = placement.time_limit
    if placement.mps is not None:
        accumulator["mps"] = placement.mps
    if placement.time_policy is not None:
        accumulator["time_policy"] = placement.time_policy
    if placement.dispatch is not None:
        base = accumulator.get(
            "dispatch", stage_procedures_module.Dispatch())
        accumulator["dispatch"] = base.merged_with(placement.dispatch)


def _resource_cap_time_limit(resource: Resource_base | None) -> str | None:
    """Resource walltime cap as HH:MM:SS, or None for local/unknown."""
    return _resource_compute_defaults(resource).get("time_limit")


def _ensure_time_limit_within_cap(
        label: str,
        time_limit: str | None,
        resource_cap: str | None,
        resource_name: str,
        ) -> None:
    if time_limit is None or resource_cap is None:
        return
    if time_limit_to_seconds(time_limit) > time_limit_to_seconds(resource_cap):
        raise ValueError(
            f"{label} time_limit {time_limit!r} exceeds resource "
            f"{resource_name!r} cap {resource_cap!r}.")


def _resolve_time_policy_fields(
        *,
        stage_name: str,
        resource_name: str,
        resource: Resource_base | None,
        placement_time_limit: str | None,
        time_policy: Time_policy_base | None,
        ) -> tuple[str | None, Time_policy_base]:
    """
    Resolve Placement time_policy + time_limit against the Resource cap.

    Returns ``(resolved_time_limit, concrete_policy)``. For adaptive policies
    ``resolved_time_limit`` is None (submit estimates walltime).
    """
    policy: Time_policy_base = (
        time_policy if time_policy is not None else Time_policy_adaptive())
    resource_cap = _resource_cap_time_limit(resource)
    _ensure_time_limit_within_cap(
        f"Placement for stage {stage_name!r}",
        placement_time_limit,
        resource_cap,
        resource_name,
    )

    if isinstance(policy, Time_policy_fixed):
        _ensure_time_limit_within_cap(
            f"Time_policy_fixed for stage {stage_name!r}",
            policy.time_limit,
            resource_cap,
            resource_name,
        )
        effective = policy.time_limit or placement_time_limit or resource_cap
        _ensure_time_limit_within_cap(
            f"Resolved fixed walltime for stage {stage_name!r}",
            effective,
            resource_cap,
            resource_name,
        )
        return effective, policy

    if isinstance(policy, Time_policy_adaptive):
        _ensure_time_limit_within_cap(
            f"Time_policy_adaptive.max_time_limit for stage {stage_name!r}",
            policy.max_time_limit,
            resource_cap,
            resource_name,
        )
        _ensure_time_limit_within_cap(
            f"Time_policy_adaptive.min_time_limit for stage {stage_name!r}",
            policy.min_time_limit,
            resource_cap,
            resource_name,
        )
        max_tl = policy.max_time_limit or placement_time_limit
        _ensure_time_limit_within_cap(
            f"Adaptive max walltime for stage {stage_name!r}",
            max_tl,
            resource_cap,
            resource_name,
        )
        if (policy.min_time_limit is not None and max_tl is not None
                and time_limit_to_seconds(policy.min_time_limit)
                > time_limit_to_seconds(max_tl)):
            raise ValueError(
                f"Stage {stage_name!r} adaptive min_time_limit "
                f"{policy.min_time_limit!r} exceeds max {max_tl!r}.")
        normalized = Time_policy_adaptive(
            estimated_performance=policy.estimated_performance,
            safety_factor=policy.safety_factor,
            min_time_limit=policy.min_time_limit,
            max_time_limit=max_tl,
        )
        return None, normalized

    # Unknown / base policy: treat as adaptive defaults.
    return None, Time_policy_adaptive(max_time_limit=placement_time_limit)


@define
class Run_settings:
    """
    Settings for seekrflow runs.
    resource: the machine that a protocol will run on
    """
    resources: typing.List[
        Resource_remote_slurm | Resource_remote_pbs | Resource_cloud_aws
    ] = field(
        default=Factory(list),)
    # A placement with an empty ``target`` ([]) matches every stage and acts
    #  as the default; more specific (longer) targets override it.
    placements: list[Placement] = field(
        default=Factory(list),
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(Placement),
            iterable_validator=validators.instance_of(list),
        ))
    
    
    def get_resource_by_name(
            self,
            resource_name: str
            ) -> Resource_base | None:
        """
        Get a resource by its name.
        """
        if resource_name == "local":
            return None
        for resource in self.resources:
            if resource.name == resource_name:
                return resource
        raise ValueError(
            f"Resource with name '{resource_name}' not found in "
            f"run_settings.resources.")

    def set_placement_resource(
            self,
            target: list[str],
            resource_name: str,
            ) -> None:
        """
        Assign a resource to a placement target, upserting an existing entry.
        """
        self.get_resource_by_name(resource_name)
        for placement in self.placements:
            if placement.target == target:
                placement.resource = resource_name
                return
        self.placements.append(Placement(target=list(target), resource=resource_name))

    def resolve_stage_execution(
            self,
            stage_name: str,
            procedure: stage_procedures_module.Stage_procedure_base,
            ) -> Resolved_execution:
        """
        Resolve the full execution policy for a stage from procedure defaults,
        placements (longest-prefix wins), and resource defaults.
        """
        address_map = stage_procedures_module.build_stage_address_map(procedure)
        if stage_name not in address_map:
            raise ValueError(
                f"Stage '{stage_name}' is not produced by any procedure.")
        address, _role = address_map[stage_name]

        accumulator: dict = {}
        policy_map = stage_procedures_module.build_stage_policy_map(procedure)
        emitted = policy_map[stage_name]
        accumulator["dispatch"] = emitted.dispatch
        if emitted.co_schedule_with is not None:
            accumulator["co_schedule_with"] = emitted.co_schedule_with

        for placement in _collect_matching_placements(
                address, self.placements):
            _apply_placement_fields(accumulator, placement)

        resource_name = accumulator.get("resource", "local")
        resource = self.get_resource_by_name(resource_name)
        defaults = _resource_compute_defaults(resource)

        cpus = accumulator.get("cpus", defaults["cpus"])
        memory_mb = accumulator.get("memory_mb", defaults["memory_mb"])
        mps = accumulator.get("mps", defaults["mps"])
        time_limit, time_policy = _resolve_time_policy_fields(
            stage_name=stage_name,
            resource_name=resource_name,
            resource=resource,
            placement_time_limit=accumulator.get("time_limit"),
            time_policy=accumulator.get("time_policy"),
        )

        dispatch = accumulator.get(
            "dispatch", stage_procedures_module.Dispatch())
        if resource is None or not resource_supports_arrays(resource):
            dispatch = stage_procedures_module.Dispatch(
                dimensions=dispatch.dimensions,
                group_size=None,
                concurrency=dispatch.concurrency,
            )
        max_concurrency = mps if mps is not None else 1
        if dispatch.concurrency > max_concurrency:
            raise ValueError(
                f"Stage '{stage_name}' has dispatch.concurrency={dispatch.concurrency} "
                f"but resource '{resource_name}' has mps={max_concurrency}. "
                f"Reduce dispatch.concurrency to <= {max_concurrency}.")

        return Resolved_execution(
            stage_name=stage_name,
            resource_name=resource_name,
            resource=resource,
            dispatch=dispatch,
            co_schedule_with=accumulator.get("co_schedule_with"),
            cpus=cpus,
            memory_mb=memory_mb,
            time_limit=time_limit,
            mps=mps,
            time_policy=time_policy,
        )

    def get_stage_resource(
            self,
            stage_name: str,
            procedure: stage_procedures_module.Stage_procedure_base,
            ) -> Resource_base | None:
        """
        Get the resource for a given stage (None means local).
        """
        return self.resolve_stage_execution(
            stage_name, procedure).resource


def validate_run_settings(
        seekrflow: "Seekrflow",
        model: typing.Any | None = None,
        ) -> None:
    """
    Validate placement targets, resource references, and co-scheduling rules.
    If model is None, it's merely a check placements and resources.
    """
    procedure = seekrflow.workflow.procedure
    address_map = stage_procedures_module.build_stage_address_map(procedure)
    all_addresses = [path for path, _name in address_map.values()]
    def _target_is_valid(target: list[str]) -> bool:
        # Assert that the target procedure/child names in the placements
        # have valid procedure/children among the stages.
        return any(
            _placement_targets_match(address, target)
            for address in all_addresses)

    seen_targets: set[tuple[str, ...]] = set()
    for placement in seekrflow.run_settings.placements:
        target_key = tuple(placement.target)
        if target_key in seen_targets:
            raise ValueError(
                f"Duplicate placement target {list(target_key)!r}.")
        seen_targets.add(target_key)
        if len(placement.target) > 0 and not _target_is_valid(placement.target):
            valid = sorted({tuple(p) for p in all_addresses})
            raise ValueError(
                f"Unknown placement target {placement.target!r}. "
                f"Valid stage address paths include: "
                f"{[list(p) for p in valid]}.")
        if placement.resource is not None:
            resource = seekrflow.run_settings.get_resource_by_name(
                placement.resource)
            resource_cap = _resource_cap_time_limit(resource)
            _ensure_time_limit_within_cap(
                f"Placement target {placement.target!r}",
                placement.time_limit,
                resource_cap,
                placement.resource,
            )
            if isinstance(placement.time_policy, Time_policy_fixed):
                _ensure_time_limit_within_cap(
                    f"Time_policy_fixed for target {placement.target!r}",
                    placement.time_policy.time_limit,
                    resource_cap,
                    placement.resource,
                )
            elif isinstance(placement.time_policy, Time_policy_adaptive):
                _ensure_time_limit_within_cap(
                    f"Time_policy_adaptive.max for target {placement.target!r}",
                    placement.time_policy.max_time_limit,
                    resource_cap,
                    placement.resource,
                )
                _ensure_time_limit_within_cap(
                    f"Time_policy_adaptive.min for target {placement.target!r}",
                    placement.time_policy.min_time_limit,
                    resource_cap,
                    placement.resource,
                )
                if placement.time_policy.estimated_performance is not None:
                    _validate_estimated_performance_scope(
                        placement, address_map, model)

    if model is None:
        return

    stage_names = [stage.name for stage in model.stages]
    for stage_name in stage_names:
        resolved = seekrflow.run_settings.resolve_stage_execution(
            stage_name, procedure)
        if resolved.co_schedule_with is None:
            continue
        stage_index = stage_names.index(stage_name)
        neighbor_name = _co_schedule_host_name(
            stage_name,
            resolved.co_schedule_with,
            stage_index,
            model.stages,
            stage_names,
        )
        neighbor_resolved = seekrflow.run_settings.resolve_stage_execution(
            neighbor_name, procedure)
        if neighbor_resolved.resource_name != resolved.resource_name:
            raise ValueError(
                f"Stage {stage_name!r} co_schedule_with "
                f"{resolved.co_schedule_with!r} requires the same resource "
                f"as neighbor {neighbor_name!r}, but "
                f"{resolved.resource_name!r} != "
                f"{neighbor_resolved.resource_name!r}.")
        # TODO: this is a problem! We need to be able to chain co-scheduled
        # stages into each other.
        #if neighbor_resolved.co_schedule_with is not None:
        #    #raise ValueError(
        #    #    f"Stage {neighbor_name!r} cannot host co-scheduled stage "
        #    #    f"{stage_name!r} because it is itself co-scheduled into "
        #    #    f"another stage.")
        #    print("Chaining multiple co-scheduled stages together - I want this to work!")
        if _dispatch_uses_array_spread(resolved.dispatch):
            raise ValueError(
                f"Stage {stage_name!r} cannot be co-scheduled: "
                f"dispatch.dimensions={resolved.dispatch.dimensions!r} requires "
                f"array spreading and cannot be fused into a neighbor job.")
        if _dispatch_uses_array_spread(neighbor_resolved.dispatch):
            raise ValueError(
                f"Host stage {neighbor_name!r} cannot co-schedule "
                f"{stage_name!r}: dispatch.dimensions="
                f"{neighbor_resolved.dispatch.dimensions!r} requires array "
                f"spreading.")


def _stage_scale_kind_from_model(
        stage: typing.Any,
        ) -> str | None:
    """Return ``\"md\"``, ``\"bd\"``, or None for non-countable scales."""
    scale_type = getattr(stage, "scale_type", None)
    if scale_type == "molecular_dynamics":
        return "md"
    if scale_type == "brownian_dynamics":
        return "bd"
    return None


def _validate_estimated_performance_scope(
        placement: Placement,
        address_map: dict,
        model: typing.Any | None,
        ) -> None:
    """
    Refuse estimated_performance when one Placement matches both MD and BD.
    """
    if model is None:
        return
    kinds: set[str] = set()
    for stage in model.stages:
        address_info = address_map.get(stage.name)
        if address_info is None:
            continue
        address, _role = address_info
        if not _placement_targets_match(address, placement.target):
            continue
        kind = _stage_scale_kind_from_model(stage)
        if kind is not None:
            kinds.add(kind)
    if "md" in kinds and "bd" in kinds:
        raise ValueError(
            f"Placement target {placement.target!r} sets "
            f"estimated_performance but matches both MD and BD stages. "
            f"Narrow the target or omit estimated_performance.")


@define
class Seekrflow:
    """
    All the inputs and parameters needed for a seekrflow calculation.
    """
    name: str = field(
        default="my_name",
        validator=validators.instance_of(str),
        )
    # NOTE: structure version 1.0 had a globus_compute_endpoint_id directly
    #  within the slurm_remote resource object. This was deprecated in 1.1
    #  in order to support different remote_interfaces, such as SSH, by
    #  including a remote_interface attribute instead with the settings 
    #  specific to that remote resource. In addition, some attributes were
    #  removed from the resource objects since remote workflows no longer
    #  employ parsl. I did not attempt to implement backwards compatibility
    #  to structure v1.0 because seekrflow was not yet released.
    # NOTE: structure version 1.2 implemented an entirely different structure
    #  framework. Versions after 1.1 include a workflow object with many nested
    #  attribute objects. Additionally, several other structures are nested as
    #  attributes within the parameterize attribute. Backwards compatibility to
    #  structure version 1.1 was not implemented because seekrflow was not yet 
    #  released. Removed parsl-related attributes.
    # NOTE: structure version 1.4 implemented the generalized, composable
    #  Workflow (single prepare-side Workflow object replacing the per-system
    #  workflow classes) and split the parameterize-side state into a separate
    #  parameterize_workflow object. Added random_seed to physical_attributes.
    # Backwards compatibility to earlier versions was not implemented because 
    # not needed.
    # NOTE: structure version 1.5 collapses parameterization into a single
    #  Parameterizer (removed parameterize_workflow), adds optional
    #  per-component param inputs (e.g. sdf_file), and allows blank MD
    #  scale_settings.system when parameterizing (prepare auto-fills from
    #  work/parameterize/). Backwards compatibility not implemented.
    structure_version: str = field(default="1.5",
                                 validator=validators.instance_of(str))
    workflow: workflow_module.Workflow = field(
        default=Factory(workflow_module.Workflow),
        validator=validators.instance_of(workflow_module.Workflow),
        )
    physical_attributes: base.Physical_attributes = field(
        default=Factory(base.Physical_attributes),
        validator=validators.instance_of(base.Physical_attributes),
        )
    work_directory: str | None = field(
        default="work",
        validator=validators.optional(validators.instance_of(str))
    )
    # Basename of the SEEKR model directory under work_directory.
    # Default (None or empty) resolves to structures.ROOT ("root").
    # Example: "root_1D" -> work/root_1D. Not a full path.
    root_directory: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str))
    )
    parameterizer: parameterizer_structures.Parameterizer | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(
            parameterizer_structures.Parameterizer)),
    )
    run_settings: Run_settings = field(
        default=Factory(Run_settings),
        validator=validators.instance_of(Run_settings),
    )

    def save(
            self,
            filename: str
            ) -> None:
        """
        Save the Seekrflow object to a JSON file.
        """
        converter: cattrs.Converter = cattrs.Converter()
        # Make sure that interited data classes are unstructured as their
        #  subtypes.
        _register_workflow_subclasses(converter)
        include_subclasses(Transfer_settings_base, converter)
        include_subclasses(Resource_base, converter)
        include_subclasses(Time_policy_base, converter)
        seekrflow_dict: dict = converter.unstructure(self)
        json_dump: str = json.dumps(seekrflow_dict, indent=4)
        with open(filename, "w") as file:
            file.write(json_dump)
        return
    
    def make_work_directory(
            self,
            work_directory: pathlib.Path | str | None = None
            ) -> None:
        """
        Make the work directory for the Seekrflow calculation.
        """
        if work_directory is not None:
            self.work_directory = str(work_directory)
        os.makedirs(self.work_directory, exist_ok=True)
        return

    def get_work_directory(self) -> pathlib.Path:
        """
        Get the directory where the preparation files are stored.
        """
        work_dir = pathlib.Path(self.work_directory)
        #os.makedirs(work_dir, exist_ok=True)
        return work_dir

    def get_parameterize_directory(self) -> pathlib.Path:
        """
        Get the directory where the preparation files are stored.
        """
        param_dir = pathlib.Path(self.work_directory) / PARAMETERIZE
        os.makedirs(param_dir, exist_ok=True)
        return param_dir

    def get_root_directory(self) -> pathlib.Path:
        """
        Get the SEEKR model root directory under the work directory.

        Uses ``root_directory`` as a subdirectory name within ``work_directory``
        when set; otherwise defaults to ``structures.ROOT`` (``"root"``).
        Must be a single leaf name under work — not an absolute or nested path.
        """
        name = self.root_directory if self.root_directory else ROOT
        name = str(name).rstrip("/\\")
        assert name and not os.path.isabs(name), \
            f"root_directory must be a subdirectory name under work, "\
            f"not an absolute path: {name!r}"
        assert os.path.basename(name) == name, \
            f"root_directory must be a single subdirectory name under work, "\
            f"not a nested path: {name!r}"
        root_dir = pathlib.Path(self.work_directory) / name
        os.makedirs(root_dir, exist_ok=True)
        return root_dir
    
    def get_run_directory(self) -> pathlib.Path:
        """
        Get the directory where the run files are stored.
        """
        run_dir = pathlib.Path(self.work_directory) / RUN
        os.makedirs(run_dir, exist_ok=True)
        return run_dir
    
def try_to_load_resources_json(
        seekrflow: "Seekrflow",
        ) -> None:
    """
    Try to load a resources JSON file and add the resources to the seekrflow
    run_settings.
    """
    home_seekrflow_resources_filename = os.path.expanduser(
        "~/.seekrflow_resources.json")
    if seekrflow.work_directory is not None:
        work_seekrflow_resources_filename = os.path.join(
            seekrflow.work_directory, ".seekrflow_resources.json")
    else:
        work_seekrflow_resources_filename = None
    if seekrflow.root_directory is not None:
        root_seekrflow_resources_filename = os.path.join(
            seekrflow.root_directory, ".seekrflow_resources.json")
    else:
        root_seekrflow_resources_filename = None

    json_str: dict | None = None
    if root_seekrflow_resources_filename is not None:
        if os.path.exists(root_seekrflow_resources_filename):
            with open(root_seekrflow_resources_filename, "r") as file:
                json_str = json.load(file)

    if json_str is None and work_seekrflow_resources_filename is not None:
        if os.path.exists(work_seekrflow_resources_filename):
            with open(work_seekrflow_resources_filename, "r") as file:
                json_str = json.load(file)

    if json_str is None:
        if os.path.exists(home_seekrflow_resources_filename):
            with open(home_seekrflow_resources_filename, "r") as file:
                json_str = json.load(file)

    if json_str is None:
        return

    converter: cattrs.Converter = cattrs.Converter()
    include_subclasses(Resource_base, converter)
    run_settings_obj: Run_settings = converter.structure(json_str, Run_settings)
    if seekrflow.run_settings is None:
        seekrflow.run_settings = Run_settings()
    for resource in run_settings_obj.resources:
        if resource.name in [r.name for r in seekrflow.run_settings.resources]:
            continue
        seekrflow.run_settings.resources.append(resource)
    return

def load_seekrflow(
        filename: str
    ) -> Seekrflow:
    """
    Load a Seekrflow object from a JSON file.
    """
    with open(filename, "r") as file:
        json_string: str = json.load(file)
    converter: cattrs.Converter = cattrs.Converter()
    _register_workflow_subclasses(converter)
    include_subclasses(Transfer_settings_base, converter)
    include_subclasses(Resource_base, converter)
    include_subclasses(Time_policy_base, converter)
    seekrflow_obj: Seekrflow = converter.structure(json_string, Seekrflow)
    try_to_load_resources_json(seekrflow_obj)
    return seekrflow_obj

def save_new_seekrflow(
        seekrflow: Seekrflow, 
        seekrflow_glob: str, 
        seekrflow_base: str, 
        save_old_seekrflow=True,
        directory = "."):
    """
    Generate a new seekrflow file. The old seekrflow file(s) will be renamed with a 
    numerical index.
    """
    
    model_path = os.path.join(directory, "seekrflow.json")
    if os.path.exists(model_path) and save_old_seekrflow:
        # This is expected, because this old model was loaded
        full_model_glob = os.path.join(directory, seekrflow_glob)
        num_globs = len(glob.glob(full_model_glob))
        new_pre_model_filename = seekrflow_base.format(num_globs)
        new_pre_model_path = os.path.join(directory, 
                                          new_pre_model_filename)
        print("Renaming model.xml to {}".format(new_pre_model_filename))
        copyfile(model_path, new_pre_model_path)
        
    print("Saving new seekrflow.json")
    seekrflow.save(model_path)
    return

def assign_default_parameterizer_settings(
        seekrflow: Seekrflow
        ) -> None:
    """
    Assign a default Parameterizer to the seekrflow structure.
    """
    seekrflow.parameterizer = parameterizer_structures.Parameterizer()
    seekrflow.physical_attributes = base.Physical_attributes()
    return
