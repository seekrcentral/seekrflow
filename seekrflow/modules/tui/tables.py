"""
modules/tui/tables.py

Turn the client's status file into the rows and columns of the TUI's three
table levels: systems, stages of one system, and details of one stage.

To add a column, append a Column to SYSTEM_COLUMNS or STAGE_COLUMNS with a
header, a width (0 fills the remaining width), a function that turns a row
dict into text, and optionally a color: one color name for the whole column,
or a dict from cell text to color name (other text stays white).
"""

import textwrap
from typing import Callable, NamedTuple

class Column(NamedTuple):
    header: str
    width: int
    text: Callable[[dict], str]
    color: str | dict[str, str] | None = None

    def color_of(self, text: str) -> str:
        if isinstance(self.color, dict):
            return self.color.get(text, "white")
        return self.color or "white"

MANAGER_COLORS = {"queued": "blue", "queued/running": "blue", "running": "green",
                  "failed": "red", "pushing": "aqua", "pulling": "aqua"}
SEMAPHORE_COLORS = {"go": "green", "wait": "yellow", "stop": "red"}
STAGE_STATE_COLORS = {"unstarted": "yellow", "started": "green",
                      "completed": "blue", "error": "red"}

def stage_workflow_key(stage_workflow: dict) -> str:
    """
    The control file names a stage workflow by its first stage.
    """
    return next(iter(stage_workflow["stages"]))

def error_text(stage_workflow: dict, stage: dict | None = None) -> str:
    """
    The most relevant error message, or "" if there is none.
    """
    text = stage_workflow.get("last_error") or stage_workflow.get("transfer_error")
    if not text and stage is not None and stage.get("stage_state") == "error":
        text = "error"
    return text or ""

def progress(stage: dict) -> float:
    """
    Mean progress over all anchors and swarms of a stage.
    """
    if stage.get("stage_state") == "completed":
        return 1.0
    values = list((stage.get("progress") or {}).values())
    return sum(values) / len(values) if values else 0.0

def stage_rows(system_status: dict) -> list[dict]:
    """
    One row per stage, joined with the stage workflow that runs it.
    """
    rows = [
        {"stage": stage, "sw": stage_workflow}
        for stage_workflow in (system_status.get("stage_workflows") or {}).values()
        for stage in stage_workflow["stages"].values()
    ]
    return sorted(rows, key=lambda row: row["stage"]["index"])

def system_rows(status: dict) -> list[dict]:
    """
    One row per system, with its stages and its currently active stage: the
    first started stage, otherwise the first one not completed.
    """
    rows = []
    for name, system_status in status.items():
        stages = stage_rows(system_status)
        active = next(
            (r for r in stages if r["stage"]["stage_state"] == "started"),
            next((r for r in stages if r["stage"]["stage_state"] != "completed"),
                 None))
        errors = [error_text(r["sw"], r["stage"]) for r in stages]
        rows.append({"name": name, "status": system_status, "stages": stages,
                     "active": active, "error": next(filter(None, errors), "")})
    return rows

def run_unit_label(key: str) -> str:
    """
    Turn a progress key like "anchor_3_swarm_None" into "anchor 3". An anchor
    of "any" means one run covers all anchors; a swarm of None means no swarms.
    """
    _, anchor, _, swarm = key.split("_")
    parts = [f"anchor {anchor}"] if anchor != "any" else []
    parts += [f"swarm {swarm}"] if swarm != "None" else []
    return ", ".join(parts) or "all anchors"

def detail_rows(row: dict, width: int) -> list[dict]:
    """
    Every field of one stage and its stage workflow, with long values
    wrapped to width.
    """
    stage, stage_workflow = row["stage"], row["sw"]
    fields = [(k, v) for k, v in stage.items() if k != "progress"]
    fields.append(("stage names", ", ".join(stage_workflow["stages"])))
    fields += [(k, v) for k, v in stage_workflow.items() if k != "stages"]
    fields += [(f"progress %: {run_unit_label(k)}", f"{100 * v:.2f}")
               for k, v in (stage.get("progress") or {}).items()]
    rows = []
    for key, value in fields:
        lines = [wrapped for line in str(value).splitlines() or [""]
                 for wrapped in textwrap.wrap(line, max(width, 10)) or [""]]
        rows += [{"key": key if i == 0 else "", "value": line}
                 for i, line in enumerate(lines)]
    return rows

SYSTEM_COLUMNS = [
    Column("System", 20, lambda r: r["name"]),
    Column("Stage", 16, lambda r:
           r["active"]["stage"]["name"] if r["active"] else "(done)"),
    Column("Manager", 15, lambda r:
           r["active"]["sw"]["manager_status"] if r["active"] else "",
           MANAGER_COLORS),
    Column("Done", 7, lambda r: "%d/%d" % (
        sum(s["stage"]["stage_state"] == "completed" for s in r["stages"]),
        len(r["stages"]))),
    Column("Semaphore", 10, lambda r: "/".join(sorted(
        {s["stage"]["semaphore"] for s in r["stages"]})), SEMAPHORE_COLORS),
    Column("Detach", 7, lambda r: "yes" if r["status"].get("detached_requested") else "",
           "red"),
    Column("Error", 0, lambda r: r["error"], "red"),
]

STAGE_COLUMNS = [
    Column("Stage", 18, lambda r: r["stage"]["name"]),
    Column("State", 10, lambda r: r["stage"]["stage_state"], STAGE_STATE_COLORS),
    Column("Progress %", 10, lambda r: f"{100 * progress(r['stage']):.2f}"),
    Column("Manager", 15, lambda r: r["sw"]["manager_status"], MANAGER_COLORS),
    Column("Semaphore", 10, lambda r: r["stage"]["semaphore"], SEMAPHORE_COLORS),
    Column("Resource", 14, lambda r: r["sw"]["resource_name"]),
    Column("Error", 0, lambda r: error_text(r["sw"], r["stage"]), "red"),
]

DETAIL_COLUMNS = [
    Column("Field", 32, lambda r: r["key"]),
    Column("Value", 0, lambda r: r["value"]),
]
