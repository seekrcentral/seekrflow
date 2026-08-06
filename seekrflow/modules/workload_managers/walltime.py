"""
Adaptive walltime estimation for remote and cloud job submits.

Estimates HH:MM:SS (or seconds) from remaining countable work and a
performance rate, with pad / min / max clamping. Used by SLURM, PBS, and AWS.
"""

from __future__ import annotations

import math
import typing

import seekrflow.modules.structures as structures

SECONDS_PER_DAY = 86400.0
MIN_REMOTE_SECONDS = 3600  # 1 hour (SLURM/PBS round early finish up)
MIN_AWS_SECONDS = 600  # 10 minutes
DEFAULT_SAFETY_FACTOR = 0.2

WorkKind = typing.Literal["ns", "trajectories"]
RoundMode = typing.Literal["hour", "minute"]


def parse_elapsed_to_seconds(elapsed: str | float | int | None) -> float | None:
    """Parse SLURM/PBS-style elapsed (``HH:MM:SS`` or ``D-HH:MM:SS``) to seconds."""
    if elapsed is None:
        return None
    if isinstance(elapsed, (int, float)):
        value = float(elapsed)
        return value if value > 0.0 else None
    text = str(elapsed).strip()
    if not text:
        return None
    try:
        return float(structures.time_limit_to_seconds(text))
    except ValueError:
        return None


def steps_to_ns(steps: float, timestep_ps: float) -> float:
    """Convert MD steps × timestep (picoseconds) to nanoseconds."""
    return float(steps) * float(timestep_ps) * 1.0e-3


def extract_work_snapshot(
        partition_progress: dict | None,
        ) -> tuple[WorkKind, dict[str, float]] | None:
    """
    Snapshot countable work units from a seekr progress partition map.

    Returns ``(kind, {unit_key: amount})`` where amount is current steps
    (MD) or trajectories completed (BD). Returns None if criteria are not
    countable or the map is empty/unusable.
    """
    if not isinstance(partition_progress, dict) or not partition_progress:
        return None

    md_amounts: dict[str, float] = {}
    bd_amounts: dict[str, float] = {}
    saw_non_countable = False

    for part_key, part in partition_progress.items():
        if not isinstance(part, dict):
            continue
        criteria = part.get("completion_criteria")
        swarms = part.get("swarms")
        if isinstance(swarms, dict) and swarms:
            for swarm_key, swarm in swarms.items():
                if not isinstance(swarm, dict):
                    continue
                swarm_criteria = swarm.get("completion_criteria", criteria)
                unit_key = f"{part_key}:{swarm_key}"
                if swarm_criteria == "number_of_steps":
                    current = swarm.get("current_step")
                    if isinstance(current, (int, float)):
                        md_amounts[unit_key] = float(current)
                elif swarm_criteria == "number_of_trajectories":
                    done = swarm.get(
                        "trajectories_completed",
                        swarm.get("current_step"),
                    )
                    if isinstance(done, (int, float)):
                        bd_amounts[unit_key] = float(done)
                else:
                    saw_non_countable = True
            continue

        if criteria == "number_of_steps":
            current = part.get("current_step")
            if isinstance(current, (int, float)):
                md_amounts[str(part_key)] = float(current)
        elif criteria == "number_of_trajectories":
            done = part.get("trajectories_completed")
            if isinstance(done, (int, float)):
                bd_amounts[str(part_key)] = float(done)
        elif criteria is not None:
            saw_non_countable = True

    if md_amounts and not bd_amounts:
        return "ns", md_amounts
    if bd_amounts and not md_amounts:
        return "trajectories", bd_amounts
    if saw_non_countable or (md_amounts and bd_amounts):
        return None
    return None


def extract_remaining_work(
        partition_progress: dict | None,
        *,
        timestep_ps: float | None = None,
        fallback_total_steps: float | None = None,
        fallback_trajectories_needed: float | None = None,
        ) -> tuple[WorkKind, float] | None:
    """
    Remaining work for walltime estimation (max over incomplete units).

    MD → ``(\"ns\", remaining_ns)`` using ``timestep_ps`` (picoseconds).
    BD → ``(\"trajectories\", remaining_trajs)``.
    Non-countable criteria → None.

    When progress lists ``number_of_steps`` / ``number_of_trajectories`` but
    has not yet populated step/traj counts (common before swarms exist),
    ``fallback_total_steps`` / ``fallback_trajectories_needed`` supply the
    target from the stage's completion criteria (current assumed 0).
    """
    if not isinstance(partition_progress, dict) or not partition_progress:
        if (fallback_total_steps is not None and timestep_ps is not None
                and timestep_ps > 0.0 and fallback_total_steps > 0.0):
            return "ns", steps_to_ns(float(fallback_total_steps), timestep_ps)
        if (fallback_trajectories_needed is not None
                and fallback_trajectories_needed > 0.0):
            return "trajectories", float(fallback_trajectories_needed)
        return None

    md_remaining: list[float] = []
    bd_remaining: list[float] = []
    saw_steps_criteria = False
    saw_traj_criteria = False
    saw_non_countable = False

    def _md_remaining_from(node: dict) -> float | None:
        if node.get("completion_criteria") != "number_of_steps":
            return None
        current = node.get("current_step")
        total = node.get("total_steps")
        if not isinstance(current, (int, float)) or not isinstance(
                total, (int, float)):
            return None
        if timestep_ps is None or timestep_ps <= 0.0:
            return None
        return steps_to_ns(max(0.0, float(total) - float(current)), timestep_ps)

    def _bd_remaining_from(node: dict) -> float | None:
        if node.get("completion_criteria") != "number_of_trajectories":
            return None
        needed = node.get("trajectories_needed")
        done = node.get("trajectories_completed")
        if not isinstance(needed, (int, float)) or not isinstance(
                done, (int, float)):
            return None
        return max(0.0, float(needed) - float(done))

    for _part_key, part in partition_progress.items():
        if not isinstance(part, dict):
            continue
        criteria = part.get("completion_criteria")
        swarms = part.get("swarms")
        if isinstance(swarms, dict) and swarms:
            for _swarm_key, swarm in swarms.items():
                if not isinstance(swarm, dict):
                    continue
                if swarm.get("attained", False):
                    continue
                swarm_criteria = swarm.get("completion_criteria", criteria)
                if swarm_criteria == "number_of_steps":
                    saw_steps_criteria = True
                    rem = _md_remaining_from({
                        **swarm,
                        "completion_criteria": "number_of_steps",
                    })
                    if rem is not None:
                        md_remaining.append(rem)
                    elif timestep_ps is None:
                        return None
                elif swarm_criteria == "number_of_trajectories":
                    saw_traj_criteria = True
                    rem = _bd_remaining_from({
                        **swarm,
                        "completion_criteria": "number_of_trajectories",
                    })
                    if rem is not None:
                        bd_remaining.append(rem)
                else:
                    saw_non_countable = True
            continue

        if part.get("attained", False):
            continue
        if criteria == "number_of_steps":
            saw_steps_criteria = True
            rem = _md_remaining_from(part)
            if rem is not None:
                md_remaining.append(rem)
            elif timestep_ps is None:
                return None
        elif criteria == "number_of_trajectories":
            saw_traj_criteria = True
            rem = _bd_remaining_from(part)
            if rem is not None:
                bd_remaining.append(rem)
        elif criteria is not None:
            saw_non_countable = True

    if saw_non_countable and not md_remaining and not bd_remaining:
        return None
    if md_remaining and not bd_remaining:
        return "ns", max(md_remaining)
    if bd_remaining and not md_remaining:
        return "trajectories", max(bd_remaining)

    # Criteria known countable, but step/traj counts not populated yet
    # (e.g. equil stage waiting on predecessor starting structures).
    if saw_steps_criteria and fallback_total_steps is not None:
        if timestep_ps is None or timestep_ps <= 0.0:
            return None
        return "ns", steps_to_ns(float(fallback_total_steps), timestep_ps)
    if saw_traj_criteria and fallback_trajectories_needed is not None:
        return "trajectories", float(fallback_trajectories_needed)
    if (not partition_progress and fallback_total_steps is not None
            and timestep_ps is not None and timestep_ps > 0.0):
        return "ns", steps_to_ns(float(fallback_total_steps), timestep_ps)
    return None


def fallback_targets_from_stage(stage: typing.Any) -> tuple[
        float | None, float | None, str | None]:
    """
    ``(total_steps, trajectories_needed, criteria_type)`` from a seekr stage.
    """
    cc = getattr(stage, "completion_criteria", None)
    if cc is None:
        return None, None, None
    criteria_type = getattr(cc, "type", None)
    total_steps = None
    trajs = None
    if criteria_type == "number_of_steps":
        raw = getattr(cc, "number_of_steps", None)
        if isinstance(raw, (int, float)) and raw > 0:
            total_steps = float(raw)
    elif criteria_type == "number_of_trajectories":
        raw = getattr(cc, "number_of_trajectories", None)
        if raw is None:
            raw = getattr(cc, "trajectories_needed", None)
        if isinstance(raw, (int, float)) and raw > 0:
            trajs = float(raw)
    return total_steps, trajs, criteria_type


def is_open_ended_criteria(criteria_type: str | None) -> bool:
    """True when walltime cannot be estimated from remaining work."""
    if criteria_type is None:
        return False
    return criteria_type not in (
        "number_of_steps",
        "number_of_trajectories",
    )


def work_done_delta(
        current: dict[str, float],
        baseline: dict[str, float] | None,
        ) -> float:
    """
    Total work advanced since baseline (sum of per-unit deltas).

    Missing baseline keys are treated as 0. Units only in baseline are ignored.
    """
    if not current:
        return 0.0
    baseline = baseline or {}
    total = 0.0
    for key, value in current.items():
        prev = float(baseline.get(key, 0.0))
        total += max(0.0, float(value) - prev)
    return total


def measure_rate(
        elapsed_seconds: float,
        work_done: float,
        kind: WorkKind,
        *,
        timestep_ps: float | None = None,
        ) -> float | None:
    """
    Performance rate from a prior attempt: ns/day or trajectories/day.

    For MD snapshots stored as steps, pass ``timestep_ps`` so ``work_done``
    (steps) is converted to ns before computing the rate.
    """
    if elapsed_seconds <= 0.0 or work_done <= 0.0:
        return None
    if kind == "ns":
        if timestep_ps is not None and timestep_ps > 0.0:
            work_ns = steps_to_ns(work_done, timestep_ps)
        else:
            work_ns = work_done
        if work_ns <= 0.0:
            return None
        return work_ns / (elapsed_seconds / SECONDS_PER_DAY)
    return work_done / (elapsed_seconds / SECONDS_PER_DAY)


def estimate_seconds_from_rate(
        remaining_amount: float,
        rate_per_day: float,
        *,
        safety_factor: float = DEFAULT_SAFETY_FACTOR,
        min_seconds: int,
        max_seconds: int,
        round_mode: RoundMode,
        finalize: bool = True,
        ) -> int:
    """
    Pad an estimate from remaining / rate.

    When ``finalize`` is True (single-stage jobs), also round and clamp to
    ``[min_seconds, max_seconds]``. When False (fused members), return the
    padded raw seconds so the caller can sum then finalize once.
    """
    if remaining_amount <= 0.0 or rate_per_day <= 0.0:
        return int(max_seconds) if finalize else 0
    raw = (remaining_amount / rate_per_day) * SECONDS_PER_DAY
    padded = raw * (1.0 + max(0.0, safety_factor))
    if not finalize:
        return max(0, int(math.ceil(padded)))
    if round_mode == "hour":
        rounded = int(math.ceil(padded / 3600.0) * 3600)
    else:
        rounded = int(math.ceil(padded / 60.0) * 60)
    return int(min(max(rounded, min_seconds), max_seconds))


def finalize_walltime_seconds(
        seconds: float | int,
        *,
        min_seconds: int,
        max_seconds: int,
        round_mode: RoundMode,
        ) -> int:
    """Round and clamp a (possibly fused) walltime total once."""
    value = max(0.0, float(seconds))
    if round_mode == "hour":
        rounded = int(math.ceil(value / 3600.0) * 3600) if value > 0 else 0
    else:
        rounded = int(math.ceil(value / 60.0) * 60) if value > 0 else 0
    if rounded <= 0:
        rounded = int(min_seconds)
    return int(min(max(rounded, min_seconds), max_seconds))


def backend_min_seconds(resource: structures.Resource_base | None) -> int:
    if isinstance(resource, structures.Resource_cloud_aws):
        return MIN_AWS_SECONDS
    return MIN_REMOTE_SECONDS


def backend_round_mode(resource: structures.Resource_base | None) -> RoundMode:
    if isinstance(resource, structures.Resource_cloud_aws):
        return "minute"
    return "hour"


def resource_max_seconds(resource: structures.Resource_base | None) -> int | None:
    cap = structures._resource_cap_time_limit(resource)
    if cap is None:
        return None
    return structures.time_limit_to_seconds(cap)


def resolve_adaptive_bounds(
        policy: structures.Time_policy_adaptive,
        resource: structures.Resource_base | None,
        ) -> tuple[int, int]:
    """Return ``(min_seconds, max_seconds)`` for an adaptive policy."""
    resource_cap = resource_max_seconds(resource)
    if resource_cap is None:
        resource_cap = backend_min_seconds(resource) * 48

    if policy.max_time_limit is not None:
        max_s = structures.time_limit_to_seconds(policy.max_time_limit)
    else:
        max_s = resource_cap

    if policy.min_time_limit is not None:
        min_s = structures.time_limit_to_seconds(policy.min_time_limit)
    else:
        min_s = backend_min_seconds(resource)

    min_s = min(min_s, max_s)
    return min_s, max_s


def resolve_rate(
        *,
        kind: WorkKind,
        estimated_performance: float | None,
        elapsed_seconds: float | None,
        work_done: float | None,
        timestep_ps: float | None = None,
        ) -> float | None:
    """
    Prefer measured same-stage rate; else Placement estimated_performance.
    """
    if elapsed_seconds is not None and work_done is not None:
        measured = measure_rate(
            elapsed_seconds,
            work_done,
            kind,
            timestep_ps=timestep_ps if kind == "ns" else None,
        )
        if measured is not None:
            return measured
    if estimated_performance is not None and estimated_performance > 0.0:
        return float(estimated_performance)
    return None


def estimate_stage_walltime_seconds(
        *,
        resolved: structures.Resolved_execution,
        partition_progress: dict | None,
        elapsed_seconds: float | None = None,
        work_baseline: dict[str, float] | None = None,
        timestep_ps: float | None = None,
        fallback_total_steps: float | None = None,
        fallback_trajectories_needed: float | None = None,
        finalize: bool = True,
        ) -> int | None:
    """
    Walltime seconds for one stage under its resolved time_policy.

    Fixed policies return the resolved ``time_limit``. Adaptive policies
    estimate from remaining work when countable; otherwise use max.

    Set ``finalize=False`` when summing fused members so min/round/max are
    applied once on the fused total instead of per stage.

    Returns ``None`` when ``finalize=False`` and the stage cannot contribute a
    countable estimate (unknown remaining work, or remaining known but no
    performance rate). Callers must treat ``None`` as “request full max”,
    never as zero seconds.
    """
    resource = resolved.resource
    policy = resolved.time_policy

    if isinstance(policy, structures.Time_policy_fixed):
        if resolved.time_limit is not None:
            return structures.time_limit_to_seconds(resolved.time_limit)
        cap = resource_max_seconds(resource)
        return int(cap if cap is not None else backend_min_seconds(resource))

    if not isinstance(policy, structures.Time_policy_adaptive):
        policy = structures.Time_policy_adaptive()

    min_s, max_s = resolve_adaptive_bounds(policy, resource)
    remaining = extract_remaining_work(
        partition_progress,
        timestep_ps=timestep_ps,
        fallback_total_steps=fallback_total_steps,
        fallback_trajectories_needed=fallback_trajectories_needed,
    )
    if remaining is None:
        return max_s if finalize else None

    kind, remaining_amount = remaining
    if remaining_amount <= 0.0:
        return min_s if finalize else 0

    snapshot = extract_work_snapshot(partition_progress)
    work_done = None
    if snapshot is not None and snapshot[0] == kind:
        work_done = work_done_delta(snapshot[1], work_baseline)

    rate = resolve_rate(
        kind=kind,
        estimated_performance=policy.estimated_performance,
        elapsed_seconds=elapsed_seconds,
        work_done=work_done,
        timestep_ps=timestep_ps if kind == "ns" else None,
    )
    if rate is None:
        # Do not return 0 under finalize=False: fused submit would treat that
        # as a successful tiny estimate and clamp to the backend minimum.
        return max_s if finalize else None

    return estimate_seconds_from_rate(
        remaining_amount,
        rate,
        safety_factor=policy.safety_factor,
        min_seconds=min_s,
        max_seconds=max_s,
        round_mode=backend_round_mode(resource),
        finalize=finalize,
    )


def partition_progress_from_status_payload(
        status_payload: dict | None,
        stage_name: str,
        ) -> dict | None:
    """
    Extract the per-partition progress map for ``stage_name`` from a seekr
    ``status(..., progress)`` JSON payload.
    """
    if not isinstance(status_payload, dict):
        return None
    root = status_payload.get("progress", status_payload)
    if not isinstance(root, dict):
        return None
    entry = None
    if stage_name in root and isinstance(root[stage_name], dict):
        entry = root[stage_name]
    else:
        for _key, value in root.items():
            if not isinstance(value, dict):
                continue
            if (value.get("stage_name") == stage_name
                    or value.get("name") == stage_name):
                entry = value
                break
    if entry is None:
        return None
    partition_map = entry.get("progress")
    if isinstance(partition_map, dict):
        return partition_map
    return None


def load_runner_walltime_state(
        model_directory: str,
        stage_name: str,
        ) -> dict:
    """
    Load walltime-related fields from ``.aws_runner`` / ``.slurm_runner`` /
    ``.pbs_runner`` latest state files, if present.
    """
    import json
    import os

    for dirname in (".aws_runner", ".slurm_runner", ".pbs_runner"):
        path = os.path.join(
            model_directory, dirname, f"{stage_name}_latest.json")
        if not os.path.isfile(path):
            continue
        try:
            with open(path, "r") as handle:
                data = json.load(handle)
            if isinstance(data, dict):
                return data
        except Exception:
            continue
    pending = os.path.join(
        model_directory, ".seekrflow_walltime", f"{stage_name}.json")
    if os.path.isfile(pending):
        try:
            with open(pending, "r") as handle:
                data = json.load(handle)
            if isinstance(data, dict):
                return data
        except Exception:
            pass
    return {}


def molecular_dynamics_timestep_ps(model: typing.Any) -> float | None:
    """Best-effort MD timestep in picoseconds from a seekr model."""
    value = None
    for name in ("get_timestep_by_type", "get_timestep_by_scale"):
        getter = getattr(model, name, None)
        if getter is None:
            continue
        try:
            value = getter("molecular_dynamics")
        except Exception:
            value = None
        if value is not None:
            break
    if value is None:
        return None
    try:
        return float(value)
    except Exception:
        pass
    raw = getattr(value, "_value", None)
    if raw is not None:
        try:
            return float(raw)
        except Exception:
            return None
    return None


def estimate_submit_time_limit(
        *,
        host_resolved: structures.Resolved_execution,
        stage_inputs: list[dict],
        ) -> str:
    """
    Walltime for a (possibly fused) submit.

    ``stage_inputs`` entries: ``resolved``, optional ``partition_progress``,
    ``elapsed_seconds``, ``work_baseline``, ``timestep_ps``, ``scale_type``.

    Adaptive fused rules:
    - logistic stages contribute 0
    - open-ended MD/BD criteria (e.g. ``reporter_progress``) force the
      **entire** fused job to the adaptive max
    - countable MD/BD stages (including when step counts are not yet in
      progress but exist on the stage completion criteria) are summed as
      padded raw seconds; min / round / max are applied **once** to the sum
    """
    policy = host_resolved.time_policy
    if isinstance(policy, structures.Time_policy_fixed):
        seconds = estimate_stage_walltime_seconds(
            resolved=host_resolved,
            partition_progress=None,
        )
        return structures.seconds_to_time_limit(seconds)

    if not isinstance(policy, structures.Time_policy_adaptive):
        policy = structures.Time_policy_adaptive()

    min_s, max_s = resolve_adaptive_bounds(policy, host_resolved.resource)
    total = 0
    any_estimate = False
    for entry in stage_inputs:
        resolved = entry["resolved"]
        scale_type = entry.get("scale_type")
        if scale_type == "logistic":
            continue

        stage_policy = resolved.time_policy
        if isinstance(stage_policy, structures.Time_policy_fixed):
            total += estimate_stage_walltime_seconds(
                resolved=resolved,
                partition_progress=None,
            )
            any_estimate = True
            continue

        criteria_type = entry.get("completion_criteria_type")
        if is_open_ended_criteria(criteria_type):
            print(
                f"[seekr-time] stage {resolved.stage_name}: open-ended "
                f"completion criteria {criteria_type!r}; requesting full max "
                f"{structures.seconds_to_time_limit(max_s)}"
            )
            return structures.seconds_to_time_limit(max_s)

        remaining = extract_remaining_work(
            entry.get("partition_progress"),
            timestep_ps=entry.get("timestep_ps"),
            fallback_total_steps=entry.get("fallback_total_steps"),
            fallback_trajectories_needed=entry.get(
                "fallback_trajectories_needed"),
        )
        if remaining is None:
            print(
                f"[seekr-time] stage {resolved.stage_name}: cannot estimate "
                f"remaining work; requesting full max "
                f"{structures.seconds_to_time_limit(max_s)}"
            )
            return structures.seconds_to_time_limit(max_s)

        seconds = estimate_stage_walltime_seconds(
            resolved=resolved,
            partition_progress=entry.get("partition_progress"),
            elapsed_seconds=entry.get("elapsed_seconds"),
            work_baseline=entry.get("work_baseline"),
            timestep_ps=entry.get("timestep_ps"),
            fallback_total_steps=entry.get("fallback_total_steps"),
            fallback_trajectories_needed=entry.get(
                "fallback_trajectories_needed"),
            finalize=False,
        )
        if seconds is None:
            hint = getattr(
                resolved.time_policy, "estimated_performance", None)
            print(
                f"[seekr-time] stage {resolved.stage_name}: remaining="
                f"{remaining[1]:.4g} {remaining[0]}, but no performance rate "
                f"(estimated_performance={hint!r}; no measured rate); "
                f"requesting full max "
                f"{structures.seconds_to_time_limit(max_s)}"
            )
            return structures.seconds_to_time_limit(max_s)
        print(
            f"[seekr-time] stage {resolved.stage_name}: remaining="
            f"{remaining[1]:.4g} {remaining[0]}, "
            f"raw_contribution={structures.seconds_to_time_limit(seconds)}"
        )
        total += seconds
        any_estimate = True

    if not any_estimate:
        print(
            f"[seekr-time] no countable stage estimates; requesting full max "
            f"{structures.seconds_to_time_limit(max_s)}"
        )
        return structures.seconds_to_time_limit(max_s)

    finalized = finalize_walltime_seconds(
        total,
        min_seconds=min_s,
        max_seconds=max_s,
        round_mode=backend_round_mode(host_resolved.resource),
    )
    print(
        f"[seekr-time] fused/job total raw="
        f"{structures.seconds_to_time_limit(int(total))} -> "
        f"{structures.seconds_to_time_limit(finalized)} "
        f"(min/round/max applied)"
    )
    return structures.seconds_to_time_limit(finalized)
