"""
modules/job/job.py

This script will be run from within the SLURM/PBS/AWS/local etc. job. 
It will take a JobSpec JSON file as input and be responsible for:
- Running the seekr calculation
- Monitoring the progress of the calculation - updating a state json file.
- Provide a potential launch point for Dragon HPC or AI.

The orchestration of these tasks is handled by radical asyncflow.
"""

import os
import time
import pathlib
import asyncio
import argparse
from concurrent.futures import ProcessPoolExecutor

from radical.asyncflow import WorkflowEngine
from rhapsody.backends import ConcurrentExecutionBackend #, DragonExecutionBackend?

try: 
    from .structures import JobSpec
except ImportError:
    from structures import JobSpec

async def status_job_writer_loop(
        job_spec: JobSpec,
        stop_event: asyncio.Event,
        ) -> None:
    import seekr.modules.structures as structures
    import seekr.status as seekr_status
    model_filename = os.path.join(job_spec.remote_root_dir, "model.json")
    assert os.path.exists(model_filename), f"Model file does not exist: {model_filename}"
    model = structures.load_model(str(model_filename))
    instruction = "progress"
    stage_states = []
    stage_infos = []
    stage_state_paths = []
    for i, stage_spec in enumerate(job_spec.stage_specs):
        stage_state_path = pathlib.Path(job_spec.remote_root_dir) \
            / ".stage_states" \
            / f"stage_state_{job_spec.internal_id}_{stage_spec.stage_index}_{job_spec.array_index}.json"
        stage_state_path.parent.mkdir(exist_ok=True)
        stage_state_paths.append(stage_state_path)
        stage_state = job_spec.get_stage_state(i)
        stage_info = job_spec.get_stage_info(i)
        stage_states.append(stage_state)
        stage_infos.append(stage_info)

    while True:
        # Write the status snaphots
        for i, stage_spec in enumerate(job_spec.stage_specs):
            stage_info = stage_infos[i]
            stage_state = stage_states[i]
            start_by_key = {
                (stage_id, anchor_id, swarm_id): step
                for stage_id, anchor_id, swarm_id, step in (
                    stage_state.stage_anchor_swarm_starting_step_list or []
                )
            }

            first_by_key = {
                (stage_id, anchor_id, swarm_id): stamp
                for stage_id, anchor_id, swarm_id, stamp in (
                    stage_state.stage_anchor_swarm_time_of_first_progress_list or []
                )
            }
            stage_state_path = stage_state_paths[i]
            stage_anchor_swarm_progress_list = []
            stage_anchor_swarm_current_step_list = []
            stage_anchor_swarm_total_steps_list = []
            stage_anchor_swarm_time_of_last_progress_list = []
            max_progress = 0.0
            stage_this_job_finished = True
            scope = "unpartitioned"
            for model_stage in model.stages:
                if getattr(model_stage, "index", None) == stage_info.stage_index:
                    scope = getattr(model_stage, "scope", "unpartitioned")
                    break

            for run_unit in job_spec.run_units:
                message_dict = seekr_status.status(
                        model,
                        instruction,
                        stage_arg=stage_info.stage_index,
                        anchor_arg=run_unit.anchor,
                        swarm_id=run_unit.swarm_id,
                        print_json=True,
                    )
                progress = message_dict.get("progress") or {}
                stage_state_dict = progress.get(stage_info.stage_index) \
                    or progress.get(str(stage_info.stage_index)) or {}
                if not stage_state_dict:
                    stage_this_job_finished = False
                    continue
                anchors = stage_state_dict.get("anchors") or {}
                if run_unit.anchor == "any" and scope == "partitioned":
                    entries = [
                        (int(key) if str(key).isdigit() else key, value)
                        for key, value in anchors.items()
                        if key != "unpartitioned" and isinstance(value, dict)
                    ]
                    if not entries:
                        stage_this_job_finished = False
                        continue
                else:
                    anchor_key = (
                        "unpartitioned" if run_unit.anchor == "any"
                        else str(run_unit.anchor)
                    )
                    entries = [(run_unit.anchor, anchors.get(anchor_key) or {})]


                for anchor, stage_anchor_dict in entries:
                    key = (stage_info.stage_index, anchor, run_unit.swarm_id)
                    if run_unit.swarm_id is None:
                        stage_swarm_dict = stage_anchor_dict
                    else:
                        stage_swarm_dict = (stage_anchor_dict.get("swarms") or {}).get(
                            str(run_unit.swarm_id)) or {}
                    if not stage_swarm_dict or (
                            "progress" not in stage_swarm_dict
                            and "finished" not in stage_swarm_dict):
                        stage_this_job_finished = False
                        continue
                    # Now we have the information specific to this stage/anchor/swarm
                    swarm_attained = stage_swarm_dict.get("finished", False)
                    if not swarm_attained:
                        stage_this_job_finished = False
                    swarm_progress = stage_swarm_dict.get("progress", 0.0)
                    swarm_current_step = stage_swarm_dict.get("current_step", None)
                    swarm_total_steps = stage_swarm_dict.get("total_steps", None)
                    if swarm_current_step is None or swarm_total_steps is None:
                        # Handle the BD situation
                        swarm_current_step = stage_swarm_dict.get("trajectories_completed")
                        swarm_total_steps = stage_swarm_dict.get("trajectories_needed")
                    stage_anchor_swarm_progress \
                        = (stage_info.stage_index, anchor, run_unit.swarm_id, swarm_progress)
                    if swarm_current_step is not None and swarm_total_steps is not None:
                        # Then this completion criteria can know steps in advance
                        if (swarm_current_step >= 0) and (key not in start_by_key):
                            start_by_key[key] = swarm_current_step
                        
                        stage_anchor_swarm_current_step \
                            = (stage_info.stage_index, anchor, run_unit.swarm_id, swarm_current_step)
                        stage_anchor_swarm_total_steps \
                            = (stage_info.stage_index, anchor, run_unit.swarm_id, swarm_total_steps)
                        stage_anchor_swarm_current_step_list.append(stage_anchor_swarm_current_step)
                        stage_anchor_swarm_total_steps_list.append(stage_anchor_swarm_total_steps)
                    
                    if swarm_progress > 0.0:
                        
                        if key not in first_by_key:
                            first_by_key[key] = time.time()

                        stage_anchor_swarm_time_of_last_progress \
                            = (stage_info.stage_index, anchor, run_unit.swarm_id, time.time())
                        stage_anchor_swarm_time_of_last_progress_list.append(
                            stage_anchor_swarm_time_of_last_progress)

                    stage_anchor_swarm_progress_list.append(stage_anchor_swarm_progress)
                    max_progress = max(max_progress, swarm_progress)

            stage_state.finished = stage_this_job_finished
            if stage_this_job_finished:
                stage_state.state = "completed"
            elif max_progress > 0.0:
                stage_state.state = "started"
            else:
                stage_state.state = "unstarted"

            stage_state.stage_anchor_swarm_progress_list = stage_anchor_swarm_progress_list
            
            stage_state.stage_anchor_swarm_starting_step_list = [
                (stage_id, anchor, swarm, step)
                for (stage_id, anchor, swarm), step in start_by_key.items()
            ]
            
            if len(stage_anchor_swarm_current_step_list) > 0:
                stage_state.stage_anchor_swarm_current_step_list = stage_anchor_swarm_current_step_list
            else:
                stage_state.stage_anchor_swarm_current_step_list = None
            if len(stage_anchor_swarm_total_steps_list) > 0:
                stage_state.stage_anchor_swarm_total_steps_list = stage_anchor_swarm_total_steps_list
            else:
                stage_state.stage_anchor_swarm_total_steps_list = None
            
            stage_state.stage_anchor_swarm_time_of_first_progress_list = [
                (stage_id, anchor, swarm, stamp)
                for (stage_id, anchor, swarm), stamp in first_by_key.items()
            ]
            
            if len(stage_anchor_swarm_time_of_last_progress_list) > 0:
                stage_state.stage_anchor_swarm_time_of_last_progress_list = stage_anchor_swarm_time_of_last_progress_list
            else:
                stage_state.stage_anchor_swarm_time_of_last_progress_list = None
            stage_state.save(path=stage_state_path)
        # wait for the next write interval
        if stop_event.is_set():
            break
        try:
            await asyncio.wait_for(
                stop_event.wait(), job_spec.status_write_interval)
        except asyncio.TimeoutError:
            pass
    return

def get_output_file_basename(
        stage_name: str,
        anchor: int | str, 
        swarm_id: int | None, 
        ) -> str:
    anchor_str = ""
    if anchor != "any":
        anchor_str = f"anchor_{anchor}_"

    swarm_id_str = ""
    if swarm_id is not None:
        swarm_id_str = f"swarm_id_{swarm_id}_"

    return f"{stage_name}_{anchor_str}{swarm_id_str}run.out"

async def do_job(
        job_spec: JobSpec,
        ) -> None:
    backend = await ConcurrentExecutionBackend(ProcessPoolExecutor())
    flow = await WorkflowEngine.create(backend=backend)

    if job_spec.telemetry_poll_interval is not None:
        telemetry_dir = pathlib.Path(job_spec.remote_root_dir) \
            / ".telemetry" / f"{job_spec.internal_id}_{job_spec.array_index}"
        telemetry = await flow.start_telemetry(
            resource_poll_interval=job_spec.telemetry_poll_interval,
            checkpoint_interval=job_spec.telemetry_poll_interval,
            checkpoint_path=telemetry_dir)
    else:
        telemetry = None

    @flow.function_task
    async def run_stage(*args):
        """
        Run a single stage of the job.
        """
        import sys
        import subprocess
        list_index = args[0]
        stage_spec = job_spec.stage_specs[list_index]
        stage_index = stage_spec.stage_index
        
        model_filename = os.path.join(job_spec.remote_root_dir, "model.json")
        assert os.path.exists(model_filename), f"Model file does not exist: {model_filename}"
        group_size = len(job_spec.run_units)
        for group_index in range(0, group_size, job_spec.concurrency):
            group_run_units = job_spec.run_units[group_index:group_index+job_spec.concurrency]
            procs = []
            for run_unit in group_run_units:
                output_basename = get_output_file_basename(
                    stage_spec.stage_name, run_unit.anchor, run_unit.swarm_id)
                output_file_path = pathlib.Path(job_spec.remote_root_dir) / "logs" \
                    / output_basename
                output_file_path.parent.mkdir(exist_ok=True)
                code = (
                    f"import seekr.modules.structures as structures;"
                    f"import seekr.run as seekr_run;"
                    f"model = structures.load_model('{model_filename}');"
                    f"seekr_run.run("
                    f"model, {stage_index}, {run_unit.anchor!r}, None, "
                    f"force_overwrite={stage_spec.force_overwrite}, " 
                    f"swarm_id={run_unit.swarm_id!r}, "
                    f"benchmark={stage_spec.benchmark}"
                    f")"
                )
                with open(output_file_path, "w") as f:
                    proc = subprocess.Popen(
                        [sys.executable, "-c", code],
                        stdout=f,
                        stderr=subprocess.STDOUT,
                    )
                procs.append(proc)
            
            for proc in procs:
                rc = proc.wait()
                if rc != 0:
                    raise RuntimeError(
                        f"seekr exited with code {rc}; output is {output_file_path}")
                
        await asyncio.sleep(1)
        return

    # Launch all the stage runs
    stage_tasks = []
    last_stage_task = None
    for list_index, stage_spec in enumerate(job_spec.stage_specs):
        stage_task = run_stage(list_index, last_stage_task)
        stage_tasks.append(stage_task)
        last_stage_task = stage_task
    # Create the monitor loop
    stop_event = asyncio.Event()
    status_writer_task = asyncio.create_task(
        status_job_writer_loop(job_spec, stop_event))
    try:
        await asyncio.gather(*stage_tasks)
    finally:
        stop_event.set()
        await status_writer_task
    
    #await finish_all_stages_successful()
    if telemetry is not None:
        await telemetry.stop()
    await flow.shutdown()
    return # The full stage state?

def main():
    parser = argparse.ArgumentParser(
        description="Script running interior to a remote/cloud Seekrflow job.")
    parser.add_argument(
        "job_spec", help="The path to the JobSpec JSON file.")
    args = parser.parse_args()
    job_spec = args.job_spec

    job_spec_path = pathlib.Path(job_spec)
    assert job_spec_path.exists(), f"JobSpec file does not exist: {job_spec_path}"
    job_spec = JobSpec.load(job_spec_path)
    asyncio.run(do_job(job_spec))

if __name__ == '__main__':
    main()