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
import seekrflow.modules.workflows.structures as workflow_structures
import seekrflow.modules.workflows.workflow as workflow_module
import seekrflow.modules.workflows.components as components_module
import seekrflow.modules.workflows.cv_specs as cv_specs_module
import seekrflow.modules.workflows.anchor_specs as anchor_specs_module
import seekrflow.modules.workflows.scale_settings as scale_settings_module
import seekrflow.modules.workflows.stage_procedures as stage_procedures_module
import seekrflow.modules.parameterize_workflow as parameterize_workflow_module
import seekrflow.modules.parameterize_structures as parameterizer_structures
import seekrflow.modules.parameters_topology_structures \
    as parameters_topology_structures
import seekrflow.modules.site_finder as site_finder


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
    # Job definition name
    # TODO: Is this chosen automatically?
    job_def_generic: str = field(
        default="",
        validator=validators.instance_of(str),
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
    transfer_settings: Transfer_settings_aws_s3 = field(
        default=Factory(Transfer_settings_aws_s3),
        )

    def get_seekr_image_uri(self):
        return f"{self.account_id}.dkr.ecr.{self.region}.amazonaws.com/seekr-engines-bd:latest"

@define
class Placement:
    """
    Deployment policy for stages whose address path begins with ``target``.
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
    cpus_per_task: int | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(int)))
    memory_per_node: int | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(int)))
    time_limit: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)))
    mps: int | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(int)))
    # time_policy: reserved for Phase 8+ (adaptive walltime).


@define
class Resolved_execution:
    """Fully-resolved, ready-to-run policy for one stage."""
    stage_name: str
    resource_name: str
    resource: Resource_base | None
    dispatch: stage_procedures_module.Dispatch
    co_schedule_with: str | None
    cpus_per_task: int | None
    memory_per_node: int | None
    time_limit: str | None
    mps: int | None


def resource_supports_arrays(resource: Resource_base | None) -> bool:
    # TODO: handle this for cloud resources.
    return isinstance(resource, (Resource_remote_slurm, Resource_remote_pbs))


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
    if placement.cpus_per_task is not None:
        accumulator["cpus_per_task"] = placement.cpus_per_task
    if placement.memory_per_node is not None:
        accumulator["memory_per_node"] = placement.memory_per_node
    if placement.time_limit is not None:
        accumulator["time_limit"] = placement.time_limit
    if placement.mps is not None:
        accumulator["mps"] = placement.mps
    if placement.dispatch is not None:
        base = accumulator.get(
            "dispatch", stage_procedures_module.Dispatch())
        accumulator["dispatch"] = base.merged_with(placement.dispatch)


@define
class Run_settings:
    """
    Settings for seekrflow runs.
    resource: the machine that a protocol will run on
    """
    resources: typing.List[
        Resource_remote_slurm | Resource_remote_pbs
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

        cpus_per_task = accumulator.get("cpus_per_task")
        memory_per_node = accumulator.get("memory_per_node")
        time_limit = accumulator.get("time_limit")
        mps = accumulator.get("mps")
        if resource is not None:
            if cpus_per_task is None:
                cpus_per_task = resource.cpus_per_task
            if memory_per_node is None:
                memory_per_node = resource.memory_per_node
            if time_limit is None:
                time_limit = resource.time_limit
            if mps is None:
                mps = resource.mps

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
            cpus_per_task=cpus_per_task,
            memory_per_node=memory_per_node,
            time_limit=time_limit,
            mps=mps,
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
            seekrflow.run_settings.get_resource_by_name(placement.resource)

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
    structure_version: str = field(default="1.4",
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
    root_directory: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str))
    )
    # TODO: resolve these two
    parameterize_workflow: \
        parameterize_workflow_module.Parameterize_workflow | None = field(
        default=Factory(parameterize_workflow_module.Parameterize_workflow),
        validator=validators.optional(validators.instance_of(
            parameterize_workflow_module.Parameterize_workflow)),
        )
    parameterizer: parameterizer_structures.Parameterizer | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(parameterizer_structures.Parameterizer)),
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
        Get the root directory where the Seekrflow files are stored.
        """
        root_dir = pathlib.Path(self.work_directory) / ROOT
        os.makedirs(root_dir, exist_ok=True)
        return root_dir
    
    def get_run_directory(self) -> pathlib.Path:
        """
        Get the directory where the run files are stored.
        """
        run_dir = pathlib.Path(self.work_directory) / RUN
        os.makedirs(run_dir, exist_ok=True)
        return run_dir
    
    def handle_ligand_indices(
            self,
            ligand_indices: str,
            ligand_resname: str,
            pdb_filename: str
            ) -> None:
        """
        Handle the ligand indices string and set the ligand_indices attribute.
        """
        if self.parameterize_workflow.has_small_molecule_ligand():
            # If the ligand resname is not provided, use the one from the seekrflow input
            if ligand_resname == "":
                ligand_resname = self.parameterize_workflow.parameterizer_information.ligand_resname
            else:
                self.parameterize_workflow.parameterizer_information.ligand_resname = ligand_resname
            
            # if the ligand indices are provided, use them preferentially
            if ligand_indices != "":
                ligand_indices = base.initialize_ref_indices(ligand_indices)
            else:
                if len(self.parameterize_workflow.ligand_indices) > 0:
                    ligand_indices = self.parameterize_workflow.ligand_indices
                elif ligand_resname != "":
                    ligand_indices = base.get_ligand_indices(pdb_filename, ligand_resname)
                else:
                    # TODO: implement some automated way to identify the ligand molecule
                    # in a molecular complex?
                    ligand_indices = []
            if len(ligand_indices) > 0:
                self.parameterize_workflow.ligand_indices = ligand_indices
            else:
                if len(self.parameterize_workflow.ligand_indices) == 0:
                    if self.parameterize_workflow.parameterizer_information.ligand_resname != "":
                        self.parameterize_workflow.ligand_indices = base.get_ligand_indices(
                            self.parameterize_workflow.parameterizer_information\
                                .receptor_ligand_pdb_filename, 
                            self.parameterize_workflow.parameterizer_information.ligand_resname)
                    else:
                        raise Exception("No ligand indices provided and no ligand "
                                        "residue name specified.")
        return
    
    def handle_receptor_indices(
            self,
            pdb_filename: str
            ) -> None:
        """
        Handle the receptor indices string and set the receptor_indices attribute.
        """
        if self.parameterize_workflow.has_small_molecule_ligand():
            ligand_indices = self.parameterize_workflow.ligand_indices
            assert len(ligand_indices) > 0, \
                "Cannot determine receptor indices without ligand " \
                "indices."
            receptor_indices = site_finder.site_finder_monte_carlo(
                pdb_filename,
                ligand_indices
            )
            self.parameterize_workflow.receptor_indices = receptor_indices

        return
    
    def assign_seekrflow_parameter_topology_files(
            self,
            parameter_topology_files: list[str]
            ) -> None:
        """
        Assign parameter and topology files to seekrflow structure.
        """
        assert parameter_topology_files is not None
        # First, we must try to discern the types of parameter files provided
        file_extensions = [os.path.splitext(f)[-1].lower() for f in parameter_topology_files]
        if ".prmtop" in file_extensions or ".parm7" in file_extensions:
            # Then assume AMBER - find the prmtop or parm7 file and assign it
            for f in parameter_topology_files:
                ext = os.path.splitext(f)[-1].lower()
                if ext == ".prmtop" or ext == ".parm7":
                    self.parameterize_workflow.solvated_system_for_md.parameters_topology \
                        = parameters_topology_structures.Amber_parameters_topology()
                    self.parameterize_workflow.solvated_system_for_md.parameters_topology\
                        .prmtop_filename = f
                    break
        elif ".psf" in file_extensions:
            # Then assume CHARMM - find the psf file and assign it
            other_parameter_files = []
            for f in parameter_topology_files:
                ext = os.path.splitext(f)[-1].lower()
                if ext == ".psf":
                    self.parameterize_workflow.solvated_system_for_md.parameters_topology \
                        = parameters_topology_structures.Charmm_parameters_topology()
                    self.parameterize_workflow.solvated_system_for_md.parameters_topology\
                        .psf_filename = f
                else:
                    other_parameter_files.append(f)
            self.parameterize_workflow.solvated_system_for_md.parameters_topology\
                .param_filename_list = other_parameter_files
        elif ".top" in file_extensions or ".gro" in file_extensions:
            raise NotImplementedError(
                "GROMACS parameter and topology files are not yet supported.")
        elif ".xml" in file_extensions:
            # Then assume OpenMM XML file
            for f in parameter_topology_files:
                ext = os.path.splitext(f)[-1].lower()
                if ext == ".xml":
                    self.parameterize_workflow.solvated_system_for_md.parameters_topology \
                        = parameters_topology_structures.Openmm_system()
                    self.parameterize_workflow.solvated_system_for_md.parameters_topology\
                        .system_filename = f
                    break
        else:
            raise ValueError(
                "Could not discern parameter and topology file types from "
                "provided files. Supported types include AMBER (.prmtop, .parm7), "
                "CHARMM (.psf), GROMACS (.top, .itp), and OpenMM (.xml).")
        return

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

def assign_default_prepare_settings(
        seekrflow: Seekrflow
        ) -> None:
    """
    Assign default parameterize-side settings to the seekrflow structure.
    """
    seekrflow.parameterize_workflow = \
        parameterize_workflow_module.Parameterize_workflow()
    seekrflow.parameterize_workflow.parameterizer_information = \
        parameterize_workflow_module.Parameterizer_information()
    seekrflow.parameterize_workflow.md_settings = \
        workflow_structures.MD_settings()
    seekrflow.parameterize_workflow.bd_settings = \
        workflow_structures.BD_settings()
    seekrflow.physical_attributes = base.Physical_attributes()
    return
