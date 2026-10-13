"""
modules/client/structures.py

Define structures for the client module, including:
- RunSession
- SystemRun - one model.json
- StagePlan - names, dependencies, resource, force, benchmark

"""

import os
import json
import copy
import typing
import asyncio
from typing import List, Dict

from attrs import define, field, validators
import cattrs
import seekr.modules.structures as seekr_structures

import seekrflow.modules.structures as structures
Seekrflow = structures.Seekrflow

DEFAULT_BACKGROUND_POLL_INTERVAL = 300.0
DEFAULT_TELEMETRY_POLL_INTERVAL = 30.0
DEFAULT_FOCUSED_POLL_INTERVAL = 5.0
DEFAULT_MAX_CONCURRENT_LOCAL_RUNS = 1
SEEKRFLOW_JSON_NAME = "seekrflow.json"

@define
class SystemRun:
    """
    An individual system run.
    """
    name: str = field(validator=validators.instance_of(str))
    seekrflow: Seekrflow = field(validator=validators.instance_of(Seekrflow))
    model: seekr_structures.Seekr_model = field(validator=validators.instance_of(seekr_structures.Seekr_model))
    semaphore_dict: Dict[str, str] = field(validator=validators.instance_of(dict))
    force_rerun_stages: List[str] = field(validator=validators.instance_of(list))
    benchmark_stage: str | None = field(validator=validators.optional(validators.instance_of(str))) # TODO: move to stage plan?
    telemetry_poll_interval: float = field(
        default=DEFAULT_TELEMETRY_POLL_INTERVAL,
        validator=validators.and_(
            validators.instance_of(float),
            validators.gt(0.0),
        ),
    )

@define
class RunSession:
    """
    The entire set of systems that are being run within this session of the 
    client, including all settings.
    """
    systemrun_objects: list[SystemRun] = field(
        validator=validators.instance_of(list),
        factory=list,
    )
    batch_directory: str | None = field(
        default=None, validator=validators.optional(
            validators.instance_of(str)))
    backend: typing.Any = field(default=None)
    workflow_engine: typing.Any = field(default=None)
    output_file: str = field(
        default="", validator=validators.instance_of(str))
    control_file: str = field(
        default="", validator=validators.instance_of(str))
    background_poll_interval: float = field(
        default=DEFAULT_BACKGROUND_POLL_INTERVAL,
        validator=validators.and_(
            validators.instance_of(float),
            validators.gt(0.0),
        ),
    )
    focused_poll_interval: float = field(
        default=DEFAULT_FOCUSED_POLL_INTERVAL,
        validator=validators.and_(
            validators.instance_of(float),
            validators.gt(0.0),
        ),
    )
    telemetry_poll_interval: float = field(
        default=DEFAULT_TELEMETRY_POLL_INTERVAL,
        validator=validators.and_(
            validators.instance_of(float),
            validators.gt(0.0),
        ),
    )
    limit_semaphores_by_resource: Dict[str, asyncio.Semaphore] = field(
        default={},
        validator=validators.instance_of(dict),
    )

# TODO: put this somewhere else? Used by Batch objects which are used in
#  prepare and analyze stages as well as the client run
def deep_merge(
        base: typing.Any, 
        override: typing.Any
        ) -> typing.Any:
    """
    Deep-merge override onto a deep copy of base.

    Dicts merge recursively. Lists merge element-by-index: dict elements
    are deep-merged; other element types are replaced. Extra override
    elements are appended; trailing base elements (omitted indices) are kept.
    Scalars and other types are replaced by override.
    """
    if override is None:
        return copy.deepcopy(base)
    if isinstance(base, dict) and isinstance(override, dict):
        result = copy.deepcopy(base)
        for key, value in override.items():
            if key in result:
                result[key] = deep_merge(result[key], value)
            else:
                result[key] = copy.deepcopy(value)
        return result
    if isinstance(base, list) and isinstance(override, list):
        result: list = []
        shared = min(len(base), len(override))
        for i in range(shared):
            if isinstance(base[i], dict) and isinstance(override[i], dict):
                result.append(deep_merge(base[i], override[i]))
            else:
                result.append(copy.deepcopy(override[i]))
        if len(override) > shared:
            result.extend(copy.deepcopy(override[shared:]))
        elif len(base) > shared:
            result.extend(copy.deepcopy(base[shared:]))
        return result
    return copy.deepcopy(override)

# TODO: put this somewhere else? this will be used by prepare and analyze stages
def absolutize_existing_paths(
        obj: typing.Any,
        base_dir: str,
        ) -> typing.Any:
    """
    Recursively turn relative path strings that exist under base_dir into
    absolute paths so children can chdir into work directories safely.
    """
    if isinstance(obj, dict):
        return {
            key: absolutize_existing_paths(value, base_dir)
            for key, value in obj.items()
        }
    if isinstance(obj, list):
        return [absolutize_existing_paths(value, base_dir) for value in obj]
    if isinstance(obj, str) and obj and not os.path.isabs(obj):
        candidate = os.path.normpath(os.path.join(base_dir, obj))
        if os.path.exists(candidate):
            return os.path.abspath(candidate)
    return obj

# TODO: put this somewhere else? this will be used by prepare and analyze stages
#  as well as the client run
@define
class Batch_system:
    """
    One system entry in a batch.
    """
    name: str = field(validator=validators.instance_of(str))
    skip: bool = field(default=False, validator=validators.instance_of(bool))
    overrides: dict = field(factory=dict, validator=validators.instance_of(dict))


# TODO: put this somewhere else? this will be used by prepare and analyze stages
#  as well as the client run
@define
class Batch:
    """
    For run stages: the set of seekrflow objects to be run, as well as 
    batch-level options and settings.
    """
    batch_directory: str | None = field(
        validator=validators.optional(validators.instance_of(str)))
    # Path string or inline seekrflow dict (validated at materialize time).
    template: typing.Any = field()
    systems: list[Batch_system] = field(factory=list)
    prepare_concurrency: int = field(
        default=1,
        validator=validators.and_(
            validators.instance_of(int),
            validators.ge(1),
        ),
    )
    background_poll_interval: float = field(
        default=DEFAULT_BACKGROUND_POLL_INTERVAL,
        validator=validators.and_(
            validators.instance_of(float),
            validators.gt(0.0),
        ),
    )
    focused_poll_interval: float = field(
        default=DEFAULT_FOCUSED_POLL_INTERVAL,
        validator=validators.and_(
            validators.instance_of(float),
            validators.gt(0.0),
        ),
    )
    telemetry_poll_interval: float = field(
        default=DEFAULT_TELEMETRY_POLL_INTERVAL,
        validator=validators.and_(
            validators.instance_of(float),
            validators.gt(0.0),
        ),
    )
    # TODO: implement
    batch_analyses: list[typing.Any] = field(
        factory=list,
    )

    # Populated after load; not part of JSON input - used to resolve relative paths
    # in inputs
    _source_path: str | None = field(default=None, eq=False, repr=False)
    _source_dir: str | None = field(default=None, eq=False, repr=False)
    _existing_seekrflow: Seekrflow | None = field(default=None, eq=False, repr=False)

    def resolve_directory(self) -> str:
        """
        Resolve the batch directory.
        """
        assert self.batch_directory is not None, "Batch directory must be set"
        path = os.path.abspath(self.batch_directory)
        os.makedirs(path, exist_ok=True)
        return path

    @classmethod
    def load_batch_file(
        cls, 
        batch_json: str) -> "Batch":
        """
        Create a Batch from a batch file.
        """
        with open(batch_json, "r") as f:
            data = json.load(f)
        converter = cattrs.Converter()
        batch: "Batch" = converter.structure(data, cls)
        batch._source_path = os.path.abspath(batch_json)
        batch._source_dir = os.path.dirname(batch._source_path)
        names = [s.name for s in batch.systems]
        if len(names) != len(set(names)):
            raise ValueError(f"Duplicate system names in batch: {names}")
        
        # Take care of directory
        batch.batch_directory = batch.resolve_directory()
        return batch

    def _load_template_dict(self) -> dict:
        if isinstance(self.template, dict):
            return copy.deepcopy(self.template)
        if not isinstance(self.template, str):
            raise TypeError(
                f"batch template must be a path string or dict, got "
                f"{type(self.template)}")
        template_path = self._template_path_resolved()
        with open(template_path, "r") as f:
            template_dict = json.load(f)
        return template_dict

    def _template_path_resolved(self) -> str | None:
        if isinstance(self.template, dict):
            return None
        if not isinstance(self.template, str):
            return None
        template_path = self.template
        if not os.path.isabs(template_path) and self._source_dir:
            template_path = os.path.join(self._source_dir, template_path)
        return os.path.abspath(template_path)

    def system_work_directory(
            self, 
            system_name: str
            ) -> str:
        """
        Obtain the work directory for a single system in this batch.
        """
        return os.path.join(
            self.resolve_directory(), f"work_{system_name}")

    def materialize_one(self, system: Batch_system) -> Seekrflow:
        """
        produce the seekrflow object for a single system in this batch
        """
        template_dict = self._load_template_dict()
        template_path = self._template_path_resolved()
        path_base = (
            os.path.dirname(template_path) if template_path
            else (self._source_dir or os.getcwd())
        )
        merged = deep_merge(template_dict, system.overrides or {})
        merged = absolutize_existing_paths(merged, path_base)
        work_dir = self.system_work_directory(system.name)
        os.makedirs(work_dir, exist_ok=True)
        #os.makedirs(os.path.join(work_dir, "logs"), exist_ok=True)
        merged["name"] = system.name
        merged["work_directory"] = work_dir
        json_path = os.path.join(work_dir, SEEKRFLOW_JSON_NAME)
        with open(json_path, "w") as f:
            json.dump(merged, f, indent=4)
        seekrflow = structures.load_seekrflow(json_path)
        return seekrflow

    def materialize_all(self) -> List[Seekrflow]:
        """
        produce the list of seekrflow objects for all systems in this batch
        """
        return [self.materialize_one(system) for system in self.systems]
    
    def create_seekrflow_objects(self) -> List[Seekrflow]:
        """
        Create a list of seekrflow objects for all systems in this batch.
        """
        if self.batch_directory is not None:
            seekrflow_objects=self.materialize_all()
        else:
            seekrflow_objects=[self._existing_seekrflow] \
                if self._existing_seekrflow else []
        return seekrflow_objects
    
    def create_session(
            self,
            output_file: str | None,
            control_file: str | None,
            ) -> RunSession:
        """
        Create a RunSession from this Batch.
        """
        seekrflow_objects = self.create_seekrflow_objects()
        systemrun_objects = []
        curdir = os.getcwd()
        for seekrflow in seekrflow_objects:
            if seekrflow.work_directory is not None:
                seekrflow.work_directory = os.path.abspath(seekrflow.work_directory)
                os.chdir(seekrflow.work_directory)
            root_directory = str(seekrflow.get_root_directory())
            model_filename = os.path.join(root_directory, "model.json")
            if not os.path.isfile(model_filename):
                raise FileNotFoundError(
                    f"No prepared model for system {seekrflow.name!r}: "
                    f"{model_filename} does not exist. Run prepare for this "
                    f"batch before client.py run.")
            model = seekr_structures.load_model(model_filename)
            systemrun_objects.append(SystemRun(
                name=seekrflow.name,
                seekrflow=seekrflow,
                model=model,
                semaphore_dict={},
                force_rerun_stages=[],
                benchmark_stage=None,
                telemetry_poll_interval=self.telemetry_poll_interval,
                #stage_plans=[],
            ))
            os.chdir(curdir)

        return RunSession(
            systemrun_objects=systemrun_objects,
            output_file=output_file,
            control_file=control_file,
            batch_directory=self.batch_directory,
            background_poll_interval=self.background_poll_interval,
            focused_poll_interval=self.focused_poll_interval,
            telemetry_poll_interval=self.telemetry_poll_interval,
        )