"""
modules/job/structures.py

Define the structures for the job.py script, including:
- JobSpec
- StageInfo
- StageState
- StageBundleState
- TelemetryState
"""

import os
import time
import json
import pathlib
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, asdict

@dataclass
class StageInfo:
    internal_id: int
    stage_index: int
    anchor_indices: Optional[List[int]]
    swarm_indices: Optional[List[str]]
    filename: str
    starting_timestamp: float

    def save(
            self, 
            root_dir_path: pathlib.Path | None = None, 
            path: pathlib.Path | None = None) -> None:
        if path is None:
            if self.filename is None:
                raise ValueError("StageInfo filename is required.")
            self.filename = pathlib.Path(self.filename).name
            path_basename = pathlib.Path(self.filename)
            path = root_dir_path / ".stage_info_states" / path_basename
        else:
            path_basename = path.name
            self.filename = str(path_basename)
        # Creates root_dir/.stage_states; fails if root_dir is missing
        path.parent.mkdir(exist_ok=True)
        path_tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
        with open(path_tmp, "w") as f:
            json.dump(asdict(self), f, indent=2)
        os.rename(path_tmp, path)
        
    @staticmethod
    def load(path: pathlib.Path) -> "StageInfo":
        with open(path, "r") as f:
            return StageInfo(**json.load(f))

@dataclass
class StageState:
    internal_id: int
    state: str = "unstarted" #'unstarted', 'started', 'completed', 'error', 'unknown'
    finished: bool = False
    stage_anchor_swarm_progress_list: Optional[List[Tuple[int, int | str, int | None, float]]] = None
    stage_anchor_swarm_starting_step_list: Optional[List[Tuple[int, int | str, int | None, int]]] = None
    stage_anchor_swarm_current_step_list: Optional[List[Tuple[int, int | str, int | None, int]]] = None
    stage_anchor_swarm_total_steps_list: Optional[List[Tuple[int, int | str, int | None, int]]] = None
    stage_anchor_swarm_time_of_first_progress_list: Optional[List[Tuple[int, int | str, int | None, float]]] = None
    stage_anchor_swarm_time_of_last_progress_list: Optional[List[Tuple[int, int | str, int | None, float]]] = None
    status_info: Optional[dict] = None
    notes: Optional[str] = None
    filename: Optional[str] = None
    stage_info_filename: Optional[str] = None
    latest_timestamp: Optional[float] = None

    def save(
            self, 
            root_dir_path: pathlib.Path | None = None, 
            path: pathlib.Path | None = None) -> None:
        if path is None:
            if self.filename is None:
                raise ValueError("StageState filename is required.")
            self.filename = pathlib.Path(self.filename).name
            path_basename = pathlib.Path(self.filename)
            path = root_dir_path / ".stage_info_states" / path_basename
        else:
            path_basename = path.name
            self.filename = str(path_basename)
        # Creates root_dir/.stage_states; fails if root_dir is missing
        path.parent.mkdir(exist_ok=True)
        path_tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
        self.latest_timestamp = time.time()
        with open(path_tmp, "w") as f:
            json.dump(asdict(self), f, indent=2)
        os.rename(path_tmp, path)

    @staticmethod
    def load(path: pathlib.Path) -> "StageState":
        with open(path, "r") as f:
            return StageState(**json.load(f))

@dataclass(frozen=True)
class RunUnit:
    """One addressable launch unit (anchor and/or swarm)."""
    anchor: int | str
    swarm_id: int | None

@dataclass
class StageSpec:
    stage_index: int
    stage_name: str
    force_overwrite: bool
    benchmark: bool

    def to_dict(self) -> dict:
        return asdict(self)
    @classmethod
    def from_dict(cls, data: dict) -> "StageSpec":
        return cls(**data)

    def save(self, path) -> None:
        path = pathlib.Path(path)
        path_tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
        with open(path_tmp, "w") as f:
            json.dump(self.to_dict(), f, indent=2)
        os.rename(path_tmp, path)

    @classmethod
    def load(cls, path) -> "JobSpec":
        with open(path) as f:
            return cls.from_dict(json.load(f))

@dataclass
class JobSpec:
    internal_id: int
    remote_root_dir: str
    status_write_interval: float
    telemetry_poll_interval: Optional[float]
    stage_specs: List[StageSpec]
    run_units: List[RunUnit]
    concurrency: int # number of run units to launch in parallel
    array_index: int # index of this run unit in the group
    schema_version: str = "1.0.0"

    def to_dict(self) -> dict:
        return asdict(self)
        
    @classmethod
    def from_dict(cls, data: dict) -> "JobSpec":
        stage_specs = [
            StageSpec(**spec) if isinstance(spec, dict) else spec
            for spec in data["stage_specs"]
        ]
        run_units = [
            RunUnit(**unit) if isinstance(unit, dict) else unit
            for unit in data["run_units"]
        ]
        return cls(**{**data, "stage_specs": stage_specs, "run_units": run_units})

    def save(self, path) -> None:
        path = pathlib.Path(path)
        path_tmp = path.with_suffix(".tmp")
        with open(path_tmp, "w") as f:
            json.dump(self.to_dict(), f, indent=2)
        os.rename(path_tmp, path)
        
    @classmethod
    def load(cls, path) -> "JobSpec":
        with open(path) as f:
            return cls.from_dict(json.load(f))

    def get_stage_state(self, index: int) -> StageState:
        stage_spec = self.stage_specs[index]
        stage_index = stage_spec.stage_index
        root_dir_path = pathlib.Path(self.remote_root_dir)
        stage_info_basename = f"stage_info_{self.internal_id}_{stage_index}.json"
        stage_info_path = root_dir_path / ".stage_info_states" / stage_info_basename
        stage_state_basename = f"stage_state_{self.internal_id}_{stage_index}_{self.array_index}.json"
        stage_state_path = root_dir_path / ".stage_states" / stage_state_basename
        if stage_state_path.exists():
            stage_state = StageState.load(stage_state_path)
        
        else:
            stage_state = StageState(internal_id=self.internal_id,
                                     stage_info_filename=str(stage_info_path))

        return stage_state

    def get_stage_info(self, index: int) -> StageInfo:
        stage_spec = self.stage_specs[index]
        stage_index = stage_spec.stage_index
        root_dir_path = pathlib.Path(self.remote_root_dir)
        stage_info_basename = f"stage_info_{self.internal_id}_{stage_index}.json"
        stage_info_path = root_dir_path / ".stage_info_states" / stage_info_basename
        stage_state_basename = f"stage_state_{self.internal_id}_{stage_index}_{self.array_index}.json"
        stage_state_path = root_dir_path / ".stage_states" / stage_state_basename
        
        if stage_info_path.exists():
            stage_info = StageInfo.load(stage_info_path)
        else:
            stage_info = StageInfo(
                internal_id=self.internal_id, 
                stage_index=stage_index,
                anchor_indices=[run_unit.anchor for run_unit in self.run_units],
                swarm_indices=[run_unit.swarm_id for run_unit in self.run_units],
                filename=str(stage_info_path),
                starting_timestamp=time.time())
        return stage_info