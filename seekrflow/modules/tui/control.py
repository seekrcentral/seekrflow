"""
modules/tui/control.py

Write user commands into the client's control file, using the same lock
file and atomic rename as modules/client/run.py so neither side clobbers
the other's changes.
"""

import os
import json
import fcntl
from typing import Callable

def update(
        control_file: str,
        change: Callable[[dict], None],
        ) -> None:
    """
    Read the control file, apply change() to it, and write it back.
    """
    with open(control_file + ".lock", "a+") as lockf:
        fcntl.flock(lockf, fcntl.LOCK_EX)
        try:
            with open(control_file, "r") as f:
                control_dict = json.load(f)
            change(control_dict)
            tmp_file = control_file + f".{os.getpid()}.tmp"
            with open(tmp_file, "w") as f:
                json.dump(control_dict, f, indent=2)
            os.rename(tmp_file, control_file)
        finally:
            fcntl.flock(lockf, fcntl.LOCK_UN)

def set_stage_workflows(
        control_file: str,
        field: str,
        value: object,
        system: str | None = None,
        stage_workflow: str | None = None,
        ) -> None:
    """
    Set a stage workflow field ("semaphore" or "transfer"). A system or
    stage_workflow of None means all of them.
    """
    def change(control_dict: dict) -> None:
        for name in [system] if system else control_dict:
            stage_workflows = control_dict[name]["stage_workflows"]
            for key in [stage_workflow] if stage_workflow else stage_workflows:
                stage_workflows[key][field] = value
    update(control_file, change)

def detach_session(control_file: str) -> None:
    """
    Detach every system, which lets the client exit while jobs keep running.
    """
    def change(control_dict: dict) -> None:
        for system_control in control_dict.values():
            system_control["detach"] = True
    update(control_file, change)
