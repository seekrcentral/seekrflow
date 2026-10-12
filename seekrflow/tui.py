"""
tui.py

Terminal user interface for a client.py session. Displays the status file
(client.py -o) as nested tables of systems, stages, and stage details, and
writes semaphore, transfer, and detach commands into the control file
(client.py -c).

Keys apply to everything in view: all systems at the top level, the focused
system at the stage level, and the focused stage's StageWorkflow at the
detail level. Detach always applies to the whole session.
"""

import os
import json
import time
import curses
import argparse

import seekrflow.modules.tui.control as control
import seekrflow.modules.tui.tables as tables

LEVEL_COLUMNS = [tables.SYSTEM_COLUMNS, tables.STAGE_COLUMNS, tables.DETAIL_COLUMNS]
SEMAPHORE_KEYS = {"g": "go", "w": "wait", "s": "stop"}
CONFIRM_PROMPTS = {"s": "Stop (permanent for this session)?",
                   "d": "Detach the whole session?"}
ENTER_KEYS = {curses.KEY_ENTER, 10, 13}
ESCAPE_KEY = 27
# Colors as (xterm 256-color number, 8-color fallback)
PALETTE = {"white": (252, curses.COLOR_WHITE), "blue": (75, curses.COLOR_BLUE),
           "green": (114, curses.COLOR_GREEN), "yellow": (221, curses.COLOR_YELLOW),
           "red": (203, curses.COLOR_RED), "aqua": (80, curses.COLOR_CYAN)}
BACKGROUND, HIGHLIGHT = 233, 238  # off-black and dark grey (256-color only)
HELP = ("Up/Down move  Enter open  Esc back  g go  w wait  s stop  "
        "t transfer  d detach session  q quit")

def color_attributes() -> dict[tuple[str, bool], int]:
    """
    The curses attribute for each (color name, highlighted) pair.
    """
    if not curses.has_colors():
        return {(name, highlighted): curses.A_REVERSE if highlighted else 0
                for name in PALETTE for highlighted in (False, True)}
    curses.start_color()
    rich = curses.COLORS >= 256
    attributes = {}
    for i, (name, (color_256, color_8)) in enumerate(PALETTE.items()):
        for highlighted in (False, True):
            pair = 1 + 2 * i + highlighted
            background = HIGHLIGHT if highlighted else BACKGROUND
            curses.init_pair(pair, color_256 if rich else color_8,
                             background if rich else curses.COLOR_BLACK)
            attributes[name, highlighted] = curses.color_pair(pair) | (
                curses.A_REVERSE if highlighted and not rich else 0)
    return attributes

class Tui:
    def __init__(self, control_file: str, status_file: str, interval: float) -> None:
        self.control_file = control_file
        self.status_file = status_file
        self.timeout_ms = int(interval * 1000)
        self.status: dict = {}
        self.path: list[str] = []  # focused system name, then stage name
        self.cursor = [0, 0, 0]
        self.message = ""

    def refresh_status(self) -> None:
        try:
            with open(self.status_file, "r") as f:
                self.status = json.load(f)
        except (OSError, ValueError) as error:
            self.message = f"Cannot read status file: {error}"
        else:
            if self.message.startswith("Cannot read status file:"):
                self.message = ""

    def rows(self, width: int) -> list[dict]:
        """
        The rows of the table at the current level.
        """
        if not self.path:
            return tables.system_rows(self.status)
        stages = tables.stage_rows(self.status.get(self.path[0], {}))
        if len(self.path) == 1:
            return stages
        row = next((r for r in stages if r["stage"]["name"] == self.path[1]), None)
        return tables.detail_rows(row, width) if row else []

    def target(self) -> tuple[str | None, str | None]:
        """
        The (system, stage workflow) that commands apply to; None means all.
        """
        if len(self.path) < 2:
            return (self.path[0] if self.path else None), None
        stages = tables.stage_rows(self.status.get(self.path[0], {}))
        row = next(r for r in stages if r["stage"]["name"] == self.path[1])
        return self.path[0], tables.stage_workflow_key(row["sw"])

    def command(self, key: str) -> None:
        system, stage_workflow = self.target()
        scope = " / ".join(filter(None, [system, stage_workflow])) or "all systems"
        if key in SEMAPHORE_KEYS:
            control.set_stage_workflows(self.control_file, "semaphore",
                SEMAPHORE_KEYS[key], system, stage_workflow)
            self.message = f"Set {SEMAPHORE_KEYS[key]} for {scope}"
        elif key == "t":
            control.set_stage_workflows(self.control_file, "transfer", True,
                system, stage_workflow)
            self.message = f"Requested transfer for {scope}"
        elif key == "d":
            control.detach_session(self.control_file)
            self.message = "Requested detach for the whole session"

    def draw(self, screen: curses.window) -> list[dict]:
        screen.erase()
        height, width = screen.getmaxyx()
        level = len(self.path)
        columns = LEVEL_COLUMNS[level]
        rows = self.rows(width - sum(c.width + 1 for c in columns[:-1]))
        self.cursor[level] = max(0, min(self.cursor[level], len(rows) - 1))
        try:
            age = f"{time.time() - os.path.getmtime(self.status_file):.0f}s ago"
        except OSError:
            age = "never"
        title = " > ".join(["seekrflow"] + self.path) + f"   (status updated {age})"
        screen.addnstr(0, 0, title, width - 1, curses.A_BOLD)
        body = height - 4
        top = max(0, self.cursor[level] - body + 1)
        for y, row in enumerate([None] + rows[top:top + body]):
            x = 0
            selected = row is not None and top + y - 1 == self.cursor[level]
            for column in columns:
                cell_width = column.width or max(1, width - x - 1)
                text = column.header if row is None else column.text(row)
                attr = self.colors[column.color_of(text), selected] if row \
                    else self.colors["white", False] | curses.A_UNDERLINE
                if x < width - 1:
                    screen.addnstr(y + 1, x, f"{text:<{cell_width}} ",
                                   min(cell_width + 1, width - x - 1), attr)
                x += cell_width + 1
        screen.addnstr(height - 2, 0, self.message, width - 1)
        screen.addnstr(height - 1, 0, HELP, width - 1, curses.A_DIM)
        return rows

    def confirm(self, screen: curses.window, prompt: str) -> bool:
        height, width = screen.getmaxyx()
        screen.move(height - 2, 0)
        screen.clrtoeol()
        screen.addnstr(height - 2, 0, prompt + " [y/N]", width - 1, curses.A_BOLD)
        screen.timeout(-1)
        answer = screen.getch()
        screen.timeout(self.timeout_ms)
        return answer == ord("y")

    def run(self, screen: curses.window) -> None:
        curses.curs_set(0)
        self.colors = color_attributes()
        screen.bkgd(" ", self.colors["white", False])
        screen.timeout(self.timeout_ms)
        while True:
            self.refresh_status()
            rows = self.draw(screen)
            key = screen.getch()
            level = len(self.path)
            name = chr(key) if 0 <= key < 256 else ""
            if name == "q":
                return
            elif key == curses.KEY_UP:
                self.cursor[level] -= 1
            elif key == curses.KEY_DOWN:
                self.cursor[level] += 1
            elif key in ENTER_KEYS and level < 2 and rows:
                row = rows[self.cursor[level]]
                self.path.append(row["name"] if level == 0 else row["stage"]["name"])
                self.cursor[level + 1] = 0
            elif key == ESCAPE_KEY and level > 0:
                self.path.pop()
            elif name in SEMAPHORE_KEYS or name in {"t", "d"}:
                if name in CONFIRM_PROMPTS \
                        and not self.confirm(screen, CONFIRM_PROMPTS[name]):
                    self.message = "Cancelled"
                    continue
                try:
                    self.command(name)
                except (OSError, ValueError, KeyError, StopIteration) as error:
                    self.message = f"Cannot write control file: {error!r}"

def main():
    parser = argparse.ArgumentParser(
        description="Terminal user interface for monitoring and controlling "
        "a client.py session.")
    parser.add_argument(
        "control_file", metavar="CONTROL_FILE", type=str,
        help="The control JSON file passed to client.py with -c.")
    parser.add_argument(
        "status_file", metavar="STATUS_FILE", type=str,
        help="The status (output) JSON file passed to client.py with -o.")
    parser.add_argument(
        "-i", "--interval", dest="interval", metavar="SECONDS", type=float,
        default=2.0, help="How often to re-read the status file. Default: 2.")
    args = parser.parse_args()
    os.environ.setdefault("ESCDELAY", "25")
    curses.wrapper(Tui(args.control_file, args.status_file, args.interval).run)

if __name__ == "__main__":
    main()
