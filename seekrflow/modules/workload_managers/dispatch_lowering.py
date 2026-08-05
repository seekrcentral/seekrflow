"""
Lower ``Resolved_execution.dispatch`` to scheduler array specs and seekr commands.

Pure logic (no scheduler I/O) lives here for testability.
"""
from __future__ import annotations

import json
import math
import os
import typing
from dataclasses import dataclass


@dataclass(frozen=True)
class StageUnitCounts:
    """Counts used to size a stage launch."""
    scope: str
    num_anchors: int | None
    num_swarms: int


@dataclass(frozen=True)
class RunUnit:
    """One addressable launch unit (anchor and/or swarm)."""
    anchor: int | str
    swarm_id: int | None


@dataclass(frozen=True)
class LoweredDispatch:
    """Concrete launch plan for one stage."""
    units: tuple[RunUnit, ...]
    group_size: int
    concurrency: int
    n_members: int
    use_array: bool


def extract_info_entry(
        message: dict,
        stage_name: str | None = None,
        ) -> dict:
    """
    Return one ``info`` entry from a seekr ``status(..., instruction='info')``
    payload. Match by ``stage_name`` when provided; otherwise return the sole
    entry or the first.
    """
    info = message.get("info") or {}
    if not isinstance(info, dict) or len(info) == 0:
        raise ValueError(
            f"seekr info instruction returned no data: {message!r}")
    if stage_name is not None:
        for entry in info.values():
            if isinstance(entry, dict) and entry.get("stage_name") == stage_name:
                return entry
        raise ValueError(
            f"No info entry for stage {stage_name!r} in {list(info.keys())!r}.")
    if len(info) == 1:
        entry = next(iter(info.values()))
        if isinstance(entry, dict):
            return entry
    raise ValueError(
        f"Multiple info entries and no stage_name to match: {list(info.keys())!r}.")


def _coerce_optional_int(value: typing.Any) -> int | None:
    if value is None:
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, str) and value.upper() in ("N/A", "NA", "NONE", ""):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def resolve_unit_counts(
        input_info: dict,
        launching_info: dict | None,
        *,
        input_is_initial: bool,
        ) -> StageUnitCounts:
    """
    Merge ``info`` results for the input stage and (when needed) the launching
    stage into counts for array sizing.
    """
    num_swarms = _coerce_optional_int(input_info.get("num_swarms"))
    if num_swarms is None:
        num_swarms = 0
    if num_swarms == 0:
        input_name = input_info.get("stage_name", "<input>")
        raise ValueError(
            f"Input stage {input_name!r} has num_swarms=0 — it has not "
            f"produced structures yet; cannot launch the downstream stage.")

    if input_is_initial:
        if launching_info is None:
            raise ValueError(
                "Launching-stage info is required when input_stage_index is 0.")
        scope = str(launching_info.get("scope", ""))
        num_anchors = _coerce_optional_int(launching_info.get("num_anchors"))
    else:
        scope = str(input_info.get("scope", ""))
        num_anchors = _coerce_optional_int(input_info.get("num_anchors"))

    if scope.upper() == "N/A":
        scope = "unpartitioned"

    return StageUnitCounts(
        scope=scope,
        num_anchors=num_anchors,
        num_swarms=num_swarms,
    )


def validate_dispatch_dimensions(
        dimensions: list[str] | None,
        counts: StageUnitCounts,
        ) -> None:
    """Raise if ``dimensions`` are incompatible with ``counts``."""
    if not dimensions:
        return
    if "anchor" in dimensions:
        if counts.scope != "partitioned":
            raise ValueError(
                f"Dispatch dimension 'anchor' requires a partitioned scope, "
                f"got {counts.scope!r}.")
        if counts.num_anchors is None or counts.num_anchors <= 0:
            raise ValueError(
                "Dispatch dimension 'anchor' requires a positive num_anchors.")
    if "swarm" in dimensions:
        if counts.num_swarms <= 1:
            raise ValueError(
                f"Dispatch dimension 'swarm' requires num_swarms > 1, "
                f"got {counts.num_swarms}.")


def build_run_units(
        dimensions: list[str] | None,
        counts: StageUnitCounts,
        ) -> list[RunUnit]:
    """
    Build the flat unit list from dispatch dimensions and launch-time counts.
    """
    validate_dispatch_dimensions(dimensions, counts)
    if not dimensions:
        return [RunUnit(anchor="any", swarm_id=None)]

    if dimensions == ["anchor"]:
        return [
            RunUnit(anchor=anchor, swarm_id=None)
            for anchor in range(counts.num_anchors or 0)
        ]

    if dimensions == ["swarm"]:
        return [
            RunUnit(anchor="any", swarm_id=swarm)
            for swarm in range(counts.num_swarms)
        ]

    if dimensions == ["anchor", "swarm"]:
        units: list[RunUnit] = []
        for anchor in range(counts.num_anchors or 0):
            for swarm in range(counts.num_swarms):
                units.append(RunUnit(anchor=anchor, swarm_id=swarm))
        return units

    raise ValueError(f"Unsupported dispatch dimensions: {dimensions!r}.")


def member_count(n_units: int, group_size: int | None) -> int:
    if n_units <= 0:
        return 0
    if group_size is None:
        return 1
    return math.ceil(n_units / group_size)


def member_unit_slice(
        member_index: int,
        n_units: int,
        group_size: int,
        ) -> slice:
    start = member_index * group_size
    end = min(start + group_size, n_units)
    return slice(start, end)


def wave_slices(n_units_in_member: int, concurrency: int) -> list[slice]:
    """Return slice objects for each concurrency wave within one member."""
    if concurrency < 1:
        concurrency = 1
    waves: list[slice] = []
    for start in range(0, n_units_in_member, concurrency):
        waves.append(slice(start, min(start + concurrency, n_units_in_member)))
    return waves


def partition_entry_is_complete(entry: dict | None) -> bool:
    """True when a seekr per-partition progress entry is done."""
    if not isinstance(entry, dict):
        return False
    if entry.get("attained"):
        return True
    prog = entry.get("progress")
    if isinstance(prog, (int, float)) and float(prog) >= 1.0:
        return True
    return False


def lookup_unit_progress_entry(
        unit: RunUnit,
        progress_map: dict,
        ) -> dict | None:
    """
    Find the per-partition progress dict for ``unit``.

    Seekr keys partitions by anchor and/or swarm id (int or str). Missing
    entries mean the unit has not started and is treated as incomplete.
    """
    if not isinstance(progress_map, dict):
        return None

    def _get(key: typing.Any) -> dict | None:
        if key in progress_map and isinstance(progress_map[key], dict):
            return progress_map[key]
        return None

    if unit.swarm_id is None:
        if unit.anchor == "any":
            return None
        return _get(unit.anchor) or _get(str(unit.anchor))

    if unit.anchor == "any":
        return _get(unit.swarm_id) or _get(str(unit.swarm_id))

    # anchor + swarm: nested map or compound key
    for ak in (unit.anchor, str(unit.anchor)):
        nested = _get(ak)
        if nested is None:
            continue
        for sk in (unit.swarm_id, str(unit.swarm_id)):
            if sk in nested and isinstance(nested[sk], dict):
                return nested[sk]
        # Nested entry may itself be the partition (no swarm nesting).
        if "progress" in nested or "attained" in nested:
            return nested

    for key in (
            f"{unit.anchor}_{unit.swarm_id}",
            f"{unit.anchor}:{unit.swarm_id}",
            (unit.anchor, unit.swarm_id),
            ):
        found = _get(key)
        if found is not None:
            return found
    return None


def filter_incomplete_units(
        units: list[RunUnit] | tuple[RunUnit, ...],
        progress_map: dict | None,
        *,
        stage_finished: bool = False,
        force_overwrite: bool = False,
        ) -> list[RunUnit]:
    """
    Drop units whose partition is already attained / complete.

    When ``force_overwrite`` is set, all units are kept. When the stage is
    finished (or every partition is complete), returns an empty list.
    Missing progress entries are treated as incomplete.
    """
    if force_overwrite:
        return list(units)
    if stage_finished:
        return []
    if not progress_map:
        return list(units)
    return [
        unit for unit in units
        if not partition_entry_is_complete(
            lookup_unit_progress_entry(unit, progress_map))
    ]


def lower_dispatch_from_units(
        units: list[RunUnit] | tuple[RunUnit, ...],
        group_size: int | None,
        concurrency: int,
        ) -> LoweredDispatch:
    """Lower an explicit (possibly filtered) unit list to an array plan."""
    units_t = tuple(units)
    n_units = len(units_t)
    if n_units == 0:
        return LoweredDispatch(
            units=(),
            group_size=1,
            concurrency=max(1, concurrency),
            n_members=0,
            use_array=False,
        )
    effective_group = group_size if group_size is not None else n_units
    if effective_group < 1:
        effective_group = max(n_units, 1)
    n_members = member_count(n_units, group_size)
    use_array = n_members > 1
    return LoweredDispatch(
        units=units_t,
        group_size=effective_group,
        concurrency=max(1, concurrency),
        n_members=n_members,
        use_array=use_array,
    )


def lower_dispatch(
        dimensions: list[str] | None,
        group_size: int | None,
        concurrency: int,
        counts: StageUnitCounts,
        ) -> LoweredDispatch:
    """
    Combine resolved dispatch fields with launch-time counts into a plan.
    """
    return lower_dispatch_from_units(
        build_run_units(dimensions, counts),
        group_size,
        concurrency,
    )


def fetch_stage_progress_local(
        model: typing.Any,
        launching_stage: typing.Any,
        status_fn: typing.Callable[..., dict] | None = None,
        ) -> tuple[bool, dict]:
    """
    Call seekr ``status(..., instruction='progress')`` for ``launching_stage``.

    Returns ``(stage_finished, partition_progress_map)``. The map is the
    per-anchor/swarm ``progress`` dict under the stage entry (may be empty).
    """
    if status_fn is None:
        import seekr.status as seekr_status
        status_fn = seekr_status.status

    stage_arg = getattr(launching_stage, "name", launching_stage)
    msg = status_fn(
        model, "progress", stage_arg=stage_arg, print_json=True)
    progress_root = msg.get("progress", {})
    if not isinstance(progress_root, dict):
        return False, {}

    stage_progress: dict = {}
    stage_index = getattr(launching_stage, "index", None)
    stage_name = getattr(launching_stage, "name", None)
    for key, value in progress_root.items():
        if not isinstance(value, dict):
            continue
        if stage_index is not None and (
                key == stage_index or key == str(stage_index)):
            stage_progress = value
            break
        if stage_name is not None and value.get("stage_name") == stage_name:
            stage_progress = value
            break
    if not stage_progress and len(progress_root) == 1:
        sole = next(iter(progress_root.values()))
        if isinstance(sole, dict):
            stage_progress = sole

    finished = bool(stage_progress.get("finished", False))
    partition_map = stage_progress.get("progress", {})
    if not isinstance(partition_map, dict):
        partition_map = {}
    return finished, partition_map


def array_member_indices(lowered: LoweredDispatch) -> list[int] | None:
    if not lowered.use_array:
        return None
    return list(range(lowered.n_members))


def needs_unit_enumeration(dimensions: list[str] | None) -> bool:
    return bool(dimensions)


def build_info_fetch_snippet(
        model_filename: str,
        launching_stage_name: str,
        input_stage_index: int,
        launching_stage_index: int,
        ) -> str:
    """
    Python snippet (no shebang) that fetches ``info`` for unit enumeration and
    prints ``__SEEKR_INFO__`` + JSON on stdout.
    """
    lines = [
        "import json, sys",
        "try:",
        "    import seekr.modules.structures as S",
        "    import seekr.status as st",
        f"    model = S.load_model({model_filename!r})",
        f"    launching_name = {launching_stage_name!r}",
        f"    input_index = {input_stage_index}",
        f"    launching_index = {launching_stage_index}",
        "    def _entry(msg, stage_name=None):",
        "        info = msg.get('info') or {}",
        "        if stage_name is not None:",
        "            for v in info.values():",
        "                if isinstance(v, dict) and v.get('stage_name') == stage_name:",
        "                    return v",
        "            raise KeyError(stage_name)",
        "        return next(iter(info.values()))",
        "    input_msg = st.status(model, instruction='info',",
        "        stage_arg=input_index if input_index else 0, print_json=True)",
        "    input_info = _entry(input_msg)",
        "    launching_info = None",
        "    if input_index == 0:",
        "        launch_msg = st.status(model, instruction='info',",
        "            stage_arg=launching_index, print_json=True)",
        "        launching_info = _entry(launch_msg, launching_name)",
        "    payload = {",
        "        'ok': True,",
        "        'input_info': input_info,",
        "        'launching_info': launching_info,",
        "        'input_is_initial': input_index == 0,",
        "    }",
        "    print('__SEEKR_INFO__' + json.dumps(payload))",
        "except Exception as e:",
        "    import traceback",
        "    print('__SEEKR_INFO__' + json.dumps({",
        "        'ok': False, 'error': str(e), 'traceback': traceback.format_exc()}))",
    ]
    return "\n".join(lines)


def parse_info_fetch_output(stdout: str) -> StageUnitCounts:
    marker = "__SEEKR_INFO__"
    payload = None
    for line in stdout.splitlines():
        line = line.strip()
        if line.startswith(marker):
            payload = json.loads(line[len(marker):])
    if payload is None or not payload.get("ok"):
        raise RuntimeError(
            f"seekr info fetch failed: {payload!r}")
    return resolve_unit_counts(
        payload["input_info"],
        payload.get("launching_info"),
        input_is_initial=payload.get("input_is_initial", False),
    )


def fetch_unit_counts_local(
        model: typing.Any,
        launching_stage: typing.Any,
        status_fn: typing.Callable[..., dict] | None = None,
        ) -> StageUnitCounts:
    """
    Call seekr ``status(..., instruction='info')`` for the input stage (and
    launching stage when the input is initial).
    """
    if status_fn is None:
        import seekr.status as seekr_status
        status_fn = seekr_status.status

    input_index = getattr(launching_stage, "input_stage_index", 0)
    input_is_initial = input_index == 0

    if input_is_initial:
        input_msg = status_fn(
            model, instruction="info", stage_arg=0, print_json=True)
        launching_msg = status_fn(
            model,
            instruction="info",
            stage_arg=getattr(launching_stage, "index", launching_stage.name),
            print_json=True,
        )
        input_info = extract_info_entry(input_msg)
        launching_info = extract_info_entry(
            launching_msg, stage_name=launching_stage.name)
        return resolve_unit_counts(
            input_info, launching_info, input_is_initial=True)

    input_msg = status_fn(
        model, instruction="info", stage_arg=input_index, print_json=True)
    input_info = extract_info_entry(input_msg)
    return resolve_unit_counts(
        input_info, None, input_is_initial=False)


def build_remote_command_string(
        stage_name: str,
        lowered: LoweredDispatch,
        *,
        force_overwrite: bool,
        benchmark: bool,
        output_filename: str,
        ) -> str:
    """
    Build the remote ``sbatch --wrap`` / PBS command for one stage launch.
    """
    units_payload = [[u.anchor, u.swarm_id] for u in lowered.units]
    if len(lowered.units) == 1 and not lowered.use_array:
        unit = lowered.units[0]
        anchor_repr = repr(unit.anchor)
        swarm_repr = "None" if unit.swarm_id is None else repr(unit.swarm_id)
        return (
            "python -c \"import seekr.modules.structures as structures; "
            "import seekr.run as seekr_run; "
            "model = structures.load_model('model.json'); "
            f"seekr_run.run(model, {stage_name!r}, {anchor_repr}, None, "
            f"{force_overwrite!r}, {swarm_repr}, benchmark={benchmark!r})\""
            f" > {output_filename}"
        )

    units_json = json.dumps(units_payload)
    driver_lines = [
        "import json, os, subprocess, sys",
        "import seekr.modules.structures as structures",
        "import seekr.run as seekr_run",
        f"units = json.loads({units_json!r})",
        f"group_size = {lowered.group_size}",
        f"concurrency = {lowered.concurrency}",
        f"stage_name = {stage_name!r}",
        f"force_overwrite = {force_overwrite!r}",
        f"benchmark = {benchmark!r}",
        "member = int(os.environ.get('AWS_BATCH_JOB_ARRAY_INDEX', "
        "os.environ.get('SLURM_ARRAY_TASK_ID', "
        "os.environ.get('PBS_ARRAY_INDEX', '0'))) or 0)",
        "start = member * group_size",
        "batch = units[start:min(start + group_size, len(units))]",
        "model = structures.load_model('model.json')",
        "for wave_start in range(0, len(batch), concurrency):",
        "    procs = []",
        "    for u in batch[wave_start:wave_start + concurrency]:",
        "        anchor, swarm = u[0], u[1]",
        "        code = (",
        "            'import seekr.modules.structures as structures; '",
        "            'import seekr.run as seekr_run; '",
        "            \"model = structures.load_model('model.json'); \"",
        "            f'seekr_run.run(model, ' + repr(stage_name) + ', ' + repr(anchor)",
        "            + ', None, ' + repr(force_overwrite) + ', ' + repr(swarm)",
        "            + ', benchmark=' + repr(benchmark) + ')'",
        "        )",
        "        procs.append(subprocess.Popen([sys.executable, '-c', code]))",
        "    for p in procs:",
        "        rc = p.wait()",
        "        if rc != 0:",
        "            sys.exit(rc)",
    ]
    driver = "\n".join(driver_lines)
    escaped = driver.replace("\\", "\\\\").replace('"', '\\"')
    if lowered.use_array:
        # Avoid interleaved writes when array members share a cwd.
        root, ext = os.path.splitext(output_filename)
        if not ext:
            ext = ".out"
        array_idx = (
            "${AWS_BATCH_JOB_ARRAY_INDEX:-"
            "${SLURM_ARRAY_TASK_ID:-"
            "${PBS_ARRAY_INDEX:-0}}}"
        )
        return f'python -c "{escaped}" > {root}_{array_idx}{ext}'
    return f'python -c "{escaped}" > {output_filename}'
