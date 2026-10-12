"""
modules/client/validation.py

Validate the run settings of a seekrflow object.
"""

import typing

import seekr.modules.structures as seekr_structures

import seekrflow.modules.structures as structures
import seekrflow.modules.workflows.stage_procedures as stage_procedures_module

def _placement_targets_match(address: list[str], target: list[str]) -> bool:
    """
    Assert that the target procedure/child names in the placements
    have valid procedure/children among the stages.
    """
    if len(target) > len(address):
        return False
    return address[:len(target)] == target

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
        resource: structures.Resource_base | None,
        ) -> dict:
    """
    Backend-native Resource fields → agnostic compute defaults for resolve.
    """
    if isinstance(resource, (structures.Resource_direct)):
        return {
            "cpus": resource.cpus_per_task,
            "memory_mb": None,
            "time_limit": resource.max_time_limit,
            "mps": resource.mps,
        }
    if isinstance(resource, (structures.Resource_remote_slurm, 
                             structures.Resource_remote_pbs)):
        return {
            "cpus": resource.cpus_per_task,
            "memory_mb": resource.memory_per_node,
            "time_limit": resource.max_time_limit,
            "mps": resource.mps,
        }
    if isinstance(resource, structures.Resource_cloud_aws):
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

def _resource_cap_time_limit(resource: structures.Resource_base | None) -> str | None:
    """
    Resource walltime cap as HH:MM:SS, or None for local/unknown.
    """
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

def _stage_scale_kind_from_model(
        stage: typing.Any,
        ) -> str | None:
    """
    Return md, bd, or None for non-countable scales.
    """
    scale_type = getattr(stage, "scale_type", None)
    if scale_type == "molecular_dynamics":
        return "md"
    if scale_type == "brownian_dynamics":
        return "bd"
    return None

def _validate_estimated_performance_scope(
        placement: structures.Placement,
        address_map: dict,
        model: seekr_structures.Seekr_model | None,
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

def _co_schedule_host_name(
        stage_name: str,
        co_schedule_with: str,
        stage_index: int,
        model_stages: list,
        ) -> str:
    """
    Determine the name of the stage that a given stage is co-scheduled with.
    """
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

def validate_run_settings(
        seekrflow: structures.Seekrflow,
        model: seekr_structures.Seekr_model | None = None,
        ) -> None:
    """
    Validate placement targets, resource references, and co-scheduling rules.
    If model is None, it's merely a check placements and resources.
    """
    procedure = seekrflow.workflow.procedure
    address_map = stage_procedures_module.build_stage_address_map(procedure)
    all_addresses = [path for path, _name in address_map.values()]
    def _target_is_valid(target: list[str]) -> bool:
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
            if isinstance(placement.time_policy, structures.Time_policy_fixed):
                _ensure_time_limit_within_cap(
                    f"Time_policy_fixed for target {placement.target!r}",
                    placement.time_policy.time_limit,
                    resource_cap,
                    placement.resource,
                )
            elif isinstance(placement.time_policy, structures.Time_policy_adaptive):
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
        resolved = resolve_stage_execution(
            seekrflow.run_settings, stage_name, procedure)
        if resolved.co_schedule_with is None:
            continue
        stage_index = stage_names.index(stage_name)
        neighbor_name = _co_schedule_host_name(
            stage_name,
            resolved.co_schedule_with,
            stage_index,
            model.stages,
        )
        neighbor_resolved = resolve_stage_execution(
            seekrflow.run_settings, neighbor_name, procedure)
        if neighbor_resolved.resource.name != resolved.resource.name:
            raise ValueError(
                f"Stage {stage_name!r} co_schedule_with "
                f"{resolved.co_schedule_with!r} requires the same resource "
                f"as neighbor {neighbor_name!r}, but "
                f"{resolved.resource.name!r} != "
                f"{neighbor_resolved.resource.name!r}.")
        if resolved.dispatch.dimensions != neighbor_resolved.dispatch.dimensions:
            raise ValueError(
                f"Host stage {neighbor_name!r} cannot co-schedule "
                f"{stage_name!r}: dispatch.dimensions="
                f"{neighbor_resolved.dispatch.dimensions!r} does not equal"
                f"{resolved.dispatch.dimensions!r}.")

def _collect_matching_placements(
        address: list[str],
        placements: list[structures.Placement],
        ) -> list[structures.Placement]:
    seen_targets: set[tuple[str, ...]] = set()
    matching: list[structures.Placement] = []
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
        placement: structures.Placement,
        ) -> None:
    if placement.resource is not None:
        accumulator["resource"] = placement.resource
    if placement.co_schedule_with is not None:
        accumulator["co_schedule_with"] = placement.co_schedule_with
    if placement.cpus is not None:
        accumulator["cpus"] = placement.cpus
    if placement.memory_mb is not None:
        accumulator["memory_mb"] = placement.memory_mb
    if placement.mps is not None:
        accumulator["mps"] = placement.mps
    if placement.dispatch is not None:
        base = accumulator.get(
            "dispatch", stage_procedures_module.Dispatch())
        accumulator["dispatch"] = base.merged_with(placement.dispatch)

def _resolve_time_policy_fields(
        *,
        stage_name: str,
        resource_name: str,
        resource: structures.Resource_base | None,
        time_policy: structures.Time_policy_base | None,
        ) -> tuple[str | None, structures.Time_policy_base]:
    """
    Resolve Placement time_policy + time_limit against the Resource cap.

    Returns ``(resolved_time_limit, concrete_policy)``. For adaptive policies
    ``resolved_time_limit`` is None (submit estimates walltime).
    """
    policy: structures.Time_policy_base = (
        time_policy if time_policy is not None else structures.Time_policy_adaptive())
    resource_cap = _resource_cap_time_limit(resource)
    if isinstance(policy, structures.Time_policy_fixed):
        _ensure_time_limit_within_cap(
            f"Time_policy_fixed for stage {stage_name!r}",
            policy.time_limit,
            resource_cap,
            resource_name,
        )
        effective = policy.time_limit or resource_cap
        _ensure_time_limit_within_cap(
            f"Resolved fixed walltime for stage {stage_name!r}",
            effective,
            resource_cap,
            resource_name,
        )
        return effective, policy

    if isinstance(policy, structures.Time_policy_adaptive):
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
        if (policy.min_time_limit is not None and policy.max_time_limit is not None
                and time_limit_to_seconds(policy.min_time_limit)
                > time_limit_to_seconds(policy.max_time_limit)):
            raise ValueError(
                f"Stage {stage_name!r} adaptive min_time_limit "
                f"{policy.min_time_limit!r} exceeds max {policy.max_time_limit!r}.")
        return None, policy

    # Unknown / base policy: treat as adaptive defaults.
    return None, structures.Time_policy_adaptive(max_time_limit=resource_cap)

def resource_supports_arrays(resource: structures.Resource_base | None) -> bool:
    return isinstance(
        resource,
        (structures.Resource_remote_slurm, structures.Resource_remote_pbs, \
            structures.Resource_cloud_aws),
    )

def resolve_stage_execution(
        run_settings: structures.Run_settings,
        stage_name: str,
        procedure: stage_procedures_module.Stage_procedure_base,
        ) -> structures.Resolved_execution:
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
            address, run_settings.placements):
        _apply_placement_fields(accumulator, placement)

    resource_name = accumulator.get("resource", "local")
    resource = run_settings.get_resource_by_name(resource_name)
    defaults = _resource_compute_defaults(resource)

    cpus = accumulator.get("cpus", defaults["cpus"])
    memory_mb = accumulator.get("memory_mb", defaults["memory_mb"])
    mps = accumulator.get("mps", defaults["mps"])
    time_limit, time_policy = _resolve_time_policy_fields(
        stage_name=stage_name,
        resource_name=resource_name,
        resource=resource,
        time_policy=accumulator.get("time_policy"),
    )

    dispatch = accumulator.get(
        "dispatch", stage_procedures_module.Dispatch())
    if not resource_supports_arrays(resource):
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

    return structures.Resolved_execution(
        stage_name=stage_name,
        resource=resource,
        dispatch=dispatch,
        co_schedule_with=accumulator.get("co_schedule_with"),
        cpus=cpus,
        memory_mb=memory_mb,
        time_limit=time_limit,
        mps=mps,
        time_policy=time_policy,
    )
