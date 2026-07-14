"""
Pure logic for co_schedule_with job-script fusion.
"""
from __future__ import annotations

import typing

if typing.TYPE_CHECKING:
    from seekrflow.modules.structures import Resource_remote_base

def combine_fused_commands(command_strings: list[str]) -> str:
    """Join per-stage remote commands, short-circuiting on failure."""
    fused_commands = " && ".join(f"({cmd})" for cmd in command_strings)
    return fused_commands


def _resource_supports_fusion(
        resource: "Resource_remote_base | None",
        ) -> bool:
    return resource is not None and resource.type in (
        "slurm_remote", "pbs_remote")


def populate_fusion_map(
        stage_workflows: list,
        benchmark_stage: str | None = None,
        ) -> None:
    """
    Wire co-scheduled stages into host fused_before / fused_after lists.
    Remote-only; local co_schedule_with is a no-op.
    """
    # stage_workflow_index_map: key - stage name, value - stage_workflow index
    stage_workflow_index_map = {}
    # upstream_stage_map: key - stage name, value - upstream stage name
    upstream_stage_map = {}
    # host_map: key - stage_name, value - host stage name
    host_map = {}

    for i, sw in enumerate(stage_workflows):
        if sw.resolved_execution is not None:
            sw.co_schedule_with = sw.resolved_execution.co_schedule_with
        else:
            sw.co_schedule_with = None
        stage_workflow_index_map[sw.stage.name] = i

    for sw in stage_workflows:
        co = sw.co_schedule_with
        if not co:
            continue
        #if sw.benchmark_mode or sw.stage.name == benchmark_stage:
        #    continue
        if sw.resource_name == "local":
            continue
        if not _resource_supports_fusion(sw.resource):
            continue
        #if sw.fused_before or sw.fused_after:
        #    raise ValueError(
        #        f"Stage {sw.stage.name!r} cannot be co-scheduled because it "
        #        f"already hosts fused stages.")

        if co == "predecessor":
            upstream_idx = getattr(sw.stage, "input_stage_index", 0)
            if upstream_idx <= 0:
                continue
            upstream_sw = next(
                (candidate for candidate in stage_workflows
                 if candidate.stage.index == upstream_idx),
                None,
            )
            if upstream_sw is None:
                continue
            upstream_stage_map[sw.stage.name] = upstream_sw.stage.name

        else:
            fused_idx = stage_workflows.index(sw)
            upstream_sw = next(
                (candidate for candidate in stage_workflows
                 if getattr(candidate.stage, "input_stage_index", 0) - 1
                 == fused_idx),
                None,
            )
            if upstream_sw is None:
                continue
            upstream_stage_map[sw.stage.name] = upstream_sw.stage.name

        for stage_name in upstream_stage_map:
            visited_indices_backstop = set()
            starting_index = stage_workflow_index_map[stage_name]
            visited_indices_backstop.add(starting_index)
            cur_upstream_name = upstream_stage_map[stage_name]
            counter = 0
            while cur_upstream_name is not None:
                host_map[stage_name] = cur_upstream_name
                next_index = stage_workflow_index_map[cur_upstream_name]
                next_name = stage_workflows[next_index].stage.name
                cur_upstream_name = upstream_stage_map.get(next_name, None)
                if next_index in visited_indices_backstop:
                    raise Exception("Infinite loop detected. Is there a cycle in the "\
                                    "co-schedule definition?")
                visited_indices_backstop.add(next_index)

    for stage_name, host_name in host_map.items():
        stage_index = stage_workflow_index_map[stage_name]
        sw = stage_workflows[stage_index]
        sw.fusion_host = host_name
        co = sw.co_schedule_with
        host_index = stage_workflow_index_map[host_name]
        host_sw = stage_workflows[host_index]
        if co == "predecessor":
            host_sw.fused_after.append(sw.stage.name)
        else:
            host_sw.fused_before.append(sw.stage.name)
        
def fusion_dependencies_satisfied(
        stage_workflow,
        stage_by_name: dict,
        ) -> bool:
    """Fused stages wait for the host job to be submitted, not host completion."""
    host_name = stage_workflow.fusion_host
    if host_name is None:
        return True
    host_sw = stage_by_name.get(host_name)
    if host_sw is None:
        return False
    if host_sw.state in ("started", "completed"):
        return True
    return host_sw.manager_status in {
        "running", "queued", "running/queued",
    }


def skips_remote_submit(stage_workflow) -> bool:
    return (
        stage_workflow.fusion_host is not None
        and stage_workflow.resource_name != "local"
    )


def is_fusion_host(stage_workflow) -> bool:
    """Remote stage that owns a fused job (submits the combined command)."""
    return (
        stage_workflow.fusion_host is None
        and stage_workflow.resource_name != "local"
        and bool(stage_workflow.fused_before or stage_workflow.fused_after)
    )


def is_fusion_member(stage_workflow) -> bool:
    """Remote stage fused into a host job (monitor-only submit)."""
    return (
        stage_workflow.fusion_host is not None
        and stage_workflow.resource_name != "local"
    )


def is_in_fused_set(stage_workflow) -> bool:
    return is_fusion_host(stage_workflow) or is_fusion_member(stage_workflow)


def fusion_host_name(stage_workflow) -> str | None:
    """Name of the root host for a fused set, or None if standalone."""
    if is_fusion_member(stage_workflow):
        return stage_workflow.fusion_host
    if is_fusion_host(stage_workflow):
        return stage_workflow.stage.name
    return None


def fused_set_members(
        host_stage_workflow: typing.Any,
        ) -> list[str]:
    return (
        list(host_stage_workflow.fused_before)
        + [host_stage_workflow.stage.name]
        + list(host_stage_workflow.fused_after)
    )


def _stage_counts_for_progress(sw: typing.Any) -> bool:
    """Exclude cheap logistic / co-scheduled stages from set progress."""
    if getattr(sw.stage, "scale_type", None) == "logistic":
        return False
    if sw.co_schedule_with is not None:
        return False
    return True


def fused_set_completed(
        host_stage_workflow_name: str,
        stage_by_name: dict,
        ) -> bool:
    host_stage_workflow = stage_by_name[host_stage_workflow_name]
    for stage_name in fused_set_members(host_stage_workflow):
        sw = stage_by_name[stage_name]
        if sw.state != "completed":
            return False
    return True


def fused_set_progress(
        host_stage_workflow_name: str,
        stage_by_name: dict,
        ) -> float:
    total_progress = 0.0
    normalize = 0.0
    host_stage_workflow = stage_by_name[host_stage_workflow_name]
    for stage_name in fused_set_members(host_stage_workflow):
        sw = stage_by_name[stage_name]
        if _stage_counts_for_progress(sw):
            total_progress += sw.progress
            normalize += 1.0
    if normalize == 0.0:
        return 0.0
    return total_progress / normalize


def mark_fused_set_completed(
        host_stage_workflow: typing.Any,
        stage_by_name: dict,
        ) -> None:
    """Mark every member of a fused set completed (probe short-circuit)."""
    for stage_name in fused_set_members(host_stage_workflow):
        sw = stage_by_name[stage_name]
        sw.state = "completed"
        sw.progress = 1.0


def classify_fused_probe_status(
        member_statuses: dict[str, dict | None],
        host_stage_name: str,
        ) -> str:
    """
    Classify a one-shot probe over all fused-set members.

    Returns ``completed``, ``reattach``, or ``submit``.
    """
    for status in member_statuses.values():
        if status is None:
            return "submit"
        stage_status = status.get("stage_status") or {}
        if stage_status.get("state") != "completed":
            break
    else:
        return "completed"

    host_status = member_statuses.get(host_stage_name)
    if host_status is not None:
        jobs = (host_status.get("manager_status") or {}).get("jobs") or []
        if jobs:
            return "reattach"
    return "submit"