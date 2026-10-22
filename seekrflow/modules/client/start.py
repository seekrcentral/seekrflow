"""
modules/client/start.py:

Tasks to perform at the start of the client run.
"""

import json
from typing import Optional, List, Dict

import seekr.modules.structures as seekr_structures

import seekrflow.modules.structures as structures
import seekrflow.modules.client.structures as client_structures
import seekrflow.modules.client.validation as client_validation
import seekrflow.modules.transfer.base as transfer_base

def input_is_batch_file(input_json: str) -> bool:
    with open(input_json, "r") as f:
        data = json.load(f)
    if "batch_directory" in data:
        return True
    return False # It's probably a single seekrflow file

def make_batch_from_single_seekrflow(
        seekrflow_object: structures.Seekrflow,
        ) -> client_structures.Batch:
    seekrflow_object_dict = seekrflow_object.to_dict()
    batch_system = client_structures.Batch_system(
        name=seekrflow_object.name,
        skip=False,
        overrides={},
    )
    batch = client_structures.Batch(
        batch_directory=None,
        template=seekrflow_object_dict,
        systems=[batch_system],
        prepare_concurrency=1,
        max_concurrent_local_runs=client_structures.DEFAULT_MAX_CONCURRENT_LOCAL_RUNS,
        telemetry_poll_interval=client_structures.DEFAULT_TELEMETRY_POLL_INTERVAL,
    )
    batch._existing_seekrflow = seekrflow_object
    return batch

def _targets_overlap(a: tuple[str, str], b: tuple[str, str]) -> bool:
    def part(x: str, y: str) -> bool:
        return x == "*" or y == "*" or x == y
    return part(a[0], b[0]) and part(a[1], b[1])

def handle_semaphore(
        model: seekr_structures.Seekr_model,
        semaphore_dict: Dict[tuple[str, str], str],
        system_name: str,
        ) -> Dict[str, str]:
    """
    Handle the semaphore.
    """
    all_stage_names = [s.name for s in model.stages]
    this_system_semaphore_dict = {}
    for (sys, stage), value in semaphore_dict.items():
        if value not in ["go", "wait", "stop"]:
            raise ValueError(
                f"Invalid semaphore value: {value}. Must be go, wait, or stop")
        if sys == "*" or sys == system_name:
            if stage == "*":
                for name in all_stage_names:
                    this_system_semaphore_dict[name] = value
            else:
                if stage not in all_stage_names:
                    raise ValueError(
                        f"Unknown stage {stage!r} in --semaphore for {system_name!r}. "
                        f"Available: {all_stage_names}")
                this_system_semaphore_dict[(stage)] = value
    
    return this_system_semaphore_dict
            
def handle_benchmark_stage(
        model: seekr_structures.Seekr_model,
        benchmark_target: tuple[str, str],
        this_system_semaphore_dict: Dict[str, str],
        system_name: str,
        ) -> str | None:
    """
    Handle the benchmark stage.
    """
    system, stage = benchmark_target
    if not (system == "*" or system == system_name):
        return None
    if stage == "*":
        raise ValueError(
            f"Cannot benchmark all stages; pass one STAGE or SYSTEM:STAGE.")
    if this_system_semaphore_dict.get(stage) == "stop":
        raise ValueError(
            f"Conflicting options: --benchmark {system}:{stage} and "
            f"--semaphore {system}:{stage}:stop. "
            "Cannot benchmark a stopped stage.")    
    all_stage_names = [s.name for s in model.stages]
    if stage not in all_stage_names:
        raise ValueError(
            f"Unknown benchmark_stage {stage!r}. "
            f"Available stages: {all_stage_names}")
    return stage

def handle_force_rerun(
        model: seekr_structures.Seekr_model,
        force_targets: Optional[List[tuple[str, str]]],
        this_system_semaphore_dict: Dict[str, str],
        system_name: str,
        ) -> List:
    """
    Handle the force rerun stage.
    """
    
    all_stage_names = [s.name for s in model.stages]
    if force_targets is None:
        force_rerun_stages: set[str] = set()
    elif len(force_targets) == 0:
        force_rerun_stages = set(all_stage_names)
    else:
        force_rerun_stages: set[str] = set()
        for fsys, fstage in force_targets:
            if not (fsys == system_name or fsys == "*"):
                continue
            if fstage == "*":
                force_rerun_stages.update(all_stage_names)
            elif fstage not in all_stage_names:
                raise ValueError(
                    f"Unknown stage {fstage!r} in --force_rerun for {system_name!r}. "
                    f"Available: {all_stage_names}")
            else:
                force_rerun_stages.add(fstage)

    if force_targets is not None:
        if len(force_targets) == 0:
            force_targets = [("*", "*")]
        for fsys, fstage in force_targets:
            for sstage, value in this_system_semaphore_dict.items():
                if value != "stop":
                    continue
                if _targets_overlap((fsys, fstage), (system_name, sstage)):
                    raise ValueError(f"Conflicting options: --force_rerun {fsys}:{fstage} and "
                                     f"--semaphore {system_name}:{sstage}:stop.")
    if force_rerun_stages:
        print(f"Force-rerun requested for stages: "
              f"{sorted(force_rerun_stages)}")
    return sorted(force_rerun_stages)

def handle_transfer(
        model: seekr_structures.Seekr_model,
        transfer_targets: Optional[List[tuple[str, str]]],
        system_name: str,
        ) -> List[str]:
    """
    Handle the transfer stage.
    """
    
    all_stage_names = [s.name for s in model.stages]
    if transfer_targets is None:
        transfer_stages: set[str] = set()
    elif len(transfer_targets) == 0:
        transfer_stages = set(all_stage_names)
    else:
        transfer_stages: set[str] = set()
        for fsys, fstage in transfer_targets:
            if not (fsys == system_name or fsys == "*"):
                continue
            if fstage == "*":
                transfer_stages.update(all_stage_names)
            elif fstage not in all_stage_names:
                raise ValueError(
                    f"Unknown stage {fstage!r} in --force_rerun for {system_name!r}. "
                    f"Available: {all_stage_names}")
            else:
                transfer_stages.add(fstage)

    if transfer_stages:
        print(f"Transfer requested for stages: "
              f"{sorted(transfer_stages)}")
    return sorted(transfer_stages)

def transfer_unique_stage_resources(
        seekrflow: structures.Seekrflow,
        stage_names: list[str],
        backwards: bool
        ) -> list[str]:
    """
    Transfer the stage resources to its remote resource.
    """
    root_directory = str(seekrflow.get_root_directory())
    resources_by_name: dict[str, structures.Resource_remote_base] = {}
    for stage_name in stage_names:
        try:
            resource = client_validation.resolve_stage_execution(
                seekrflow.run_settings,
                stage_name,
                seekrflow.workflow.procedure,
            ).resource
        except ValueError:
            resource = None
        if resource is None:
            continue
        resources_by_name[resource.name] = resource

    # TODO: check if the remote resource already has the files?
    for resource in resources_by_name.values():
        transfer_base.transfer_files_to_from_remote_resource(
            seekrflow.name,
            resource,
            root_directory,
            backwards=backwards,
            seekrflow_name=seekrflow.name,
        )
    return list(resources_by_name.keys())