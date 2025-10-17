"""
modules/seekr_run.py

Using the seekr API as well as parsl, run stages of the seekr calculation.
"""

import os
import sys
import uuid
import shutil
import typing
import select
import datetime
from collections import defaultdict

import parsl
from parsl.app.app import python_app
from parsl.providers import LocalProvider
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
import seekr2.modules.common_base as seekr2_base
import seekr2.run as seekr2_run

import seekrflow.modules.structures as structures
import seekrflow.modules.transfer.globus as transfer_globus
import seekrflow.modules.transfer.rsync as transfer_rsync
import seekrflow.modules.workload_managers.slurm as workload_slurm
import seekrflow.modules.remote_interfaces.globus_compute_sdk as remote_globus
import seekrflow.modules.remote_interfaces.ssh as remote_ssh

def bd_finished(
        model: seekr2_base.Model,
        ) -> bool: 
    """
    Check if the BD stage has finished. This is done by reading the
    results XML files and seeing how many BD steps have elapsed.
    """
    bd_milestone_info_to_run = seekr2_run.choose_next_simulation_browndye2(
            model, "any_bd", None, False, None)
    if len(bd_milestone_info_to_run) > 0:
        return False
    else:
        return True

def hidr_finished(
        model: seekr2_base.Model
        ) -> bool:
    """
    Check if the Hidr stage has finished. This is done by reading the
    model files and seeing if all anchors have starting PDB files.
    """
    for alpha, anchor in enumerate(model.anchors):
        if anchor.bulkstate:
            continue
        anchor_pdb_filename = seekr2_base.get_anchor_pdb_filename(anchor)
        if anchor_pdb_filename == "":
            return False
        
    return True

def seekr_anchors_to_run(
        model: seekr2_base.Model
        ) -> bool:
    """
    Check if SEEKR has finished.
    """
    md_info_to_run = seekr2_run.choose_next_simulation_openmm(
        model, "any_md", None, None, None, None, False, False, None)
    anchor_indices_to_run = []
    for md_info in md_info_to_run:
        anchor_indices_to_run.append(md_info[2])
    return anchor_indices_to_run

@python_app(executors=["bd_executor"])
def run_bd_locally(inputs):
    import seekr2.modules.common_base as seekr2_base
    import seekr2.run as seekr2_run
    model_filename = inputs[0]
    n_threads = inputs[1]
    model = seekr2_base.load_model(model_filename)
    seekr2_run.run(model, "any_bd", n_threads=n_threads)
    return

@python_app(executors=["hidr_executor"])
def run_hidr_metaD_locally(inputs):
    import os
    import seekr2.modules.common_base as seekr2_base
    import seekrtools.hidr.hidr as seekr2_hidr
    model_filename = inputs[0]
    input_pdb_file = inputs[1]
    model = seekr2_base.load_model(model_filename)
    os.chdir(os.path.abspath(model.anchor_rootdir))
    # TODO: cuda_device_index?
    seekr2_hidr.hidr(model, "any", pdb_files=[input_pdb_file], mode="MetaD")
    return model.openmm_settings.cuda_platform_settings.cuda_device_index

@python_app(executors=["seekr_executor"])
def run_seekr_locally(inputs):
    import seekr2.modules.common_base as seekr2_base
    import seekr2.run as seekr2_run
    model_filename = inputs[0]
    anchor_index = inputs[1]
    model = seekr2_base.load_model(model_filename)
    seekr2_run.run(model, f"{anchor_index}")

def assign_local_executor(
        label: str,
        cores_per_worker: int = 1
        ) -> HighThroughputExecutor:
    """
    Assign a local executor for parsl.
    """
    return HighThroughputExecutor(
        label=label,
        worker_debug=True,
        cores_per_worker=cores_per_worker,
        provider=LocalProvider(
            init_blocks=1,
            max_blocks=1,
        ),
    )

def transfer_files_to_from_remote_resource(
        remote_root_directory_name: str,
        resource: structures.Resource_remote_base,
        local_directory: str,
        backwards: bool = False,
        ) -> None:
    
    if resource.transfer_settings.type == "globus":
        remote_path = os.path.join(
            resource.remote_working_directory, remote_root_directory_name)
        local_collection_id = resource.transfer_settings.local_collection_id
        remote_collection_id = resource.transfer_settings.remote_collection_id
        transfer_globus.transfer_files_with_globus(
            remote_root_directory_name, local_directory, remote_path, 
            local_collection_id, remote_collection_id, backwards=backwards)
                                    
        # TODO: confirm that the files were transferred successfully before continuing.
        # To prevent later errors if the file transfer failed quietly (as has happened).
        # NOTE: Implemented checksums: Fixed ??
        
    elif resource.transfer_settings.type == "rsync":
        remote_path = os.path.join(
            resource.remote_working_directory, remote_root_directory_name)
        remote_hostname = resource.transfer_settings.remote_hostname
        remote_username = resource.transfer_settings.remote_username
        remote_password = resource.transfer_settings.remote_password
        port = resource.transfer_settings.port
        private_key_filename = resource.transfer_settings.private_key_filename
        private_key_passphrase = resource.transfer_settings.private_key_passphrase
        transfer_rsync.transfer_files_with_rsync(
            local_directory, remote_path, remote_hostname, remote_username, 
            remote_password, port, private_key_filename=private_key_filename,
            private_key_passphrase=private_key_passphrase, 
            backwards=backwards)
    else:
        raise NotImplementedError(
            "Only rsync and globus transfers are implemented.")
    return

def assign_remote_workflow(
        name: str,
        destination_path: str,
        resource: structures.Resource_remote_base,
        usage_tracking: int,
        command_string: str,
        model_filename: str,
        workflow_type: str,
        indices: list | None = None,
        resubmit_on_timeout: bool = True
    ) -> typing.Any:
    
    log_directory = os.path.join(destination_path, "logs")
    if resource.type == "slurm_remote":
        slurm_runner_dir = os.path.join(destination_path, ".slurm_runner")
        stop_flag = os.path.join(slurm_runner_dir, f"stop_{uuid.uuid4()}")
        slurm_args = [
            destination_path,
            log_directory,
            resource.partition,
            resource.account,
            resource.constraint,
            name,
            resource.cores_per_node, # TODO: replace with cpus_per_task
            resource.memory_per_node,
            resubmit_on_timeout,
            stop_flag,
            resource.time_limit,
            resource.scheduler_options,
            resource.worker_init,
            usage_tracking,
            command_string,
            indices,
            model_filename,
            workflow_type
        ]
        run_workflow = workload_slurm.slurm_remote_run_workflow
        cancel_workflow = workload_slurm.slurm_remote_cancel_workflow
        status_workflow = workload_slurm.slurm_remote_status_workflow
    else:
        raise NotImplementedError(
            f"Resource type {resource.type} is not implemented.")

    return run_workflow, cancel_workflow, status_workflow, slurm_args
            
def run_model(
        seekrflow: structures.Seekrflow,
        transfer_before: bool = False,
        transfer_from_host_only: bool = False,
        force_overwrite: bool = False,
        ) -> None:
    """
    Run the SEEKR calculation using parsl.
    """
    seekrflow.work_directory = os.path.abspath(seekrflow.work_directory)
    curdir = os.getcwd()
    os.chdir(seekrflow.work_directory)
    source_directory = os.path.join(seekrflow.work_directory, structures.ROOT)
    # Load the model
    model_filename = os.path.join(source_directory, "model.xml")
    model = seekr2_base.load_model(model_filename)
    
    if seekrflow.bd_settings is not None:
        bd_n_threads = seekrflow.bd_settings.num_threads
    
    if seekrflow.run_settings.allow_parsl_usage_tracking:
        usage_tracking = 3
    else:
        usage_tracking = 0

    transferred_files = defaultdict(bool)
    local_executors = []

    anchors_to_run = seekr_anchors_to_run(model)

    if (seekrflow.run_settings.bd_stage_resource_name == "local") \
            and (model.using_bd()):
        bd_local_executor = assign_local_executor(
            label="bd_executor",
            cores_per_worker=bd_n_threads)
        local_executors.append(bd_local_executor)
        bd_running_locally = True
    else:
        if model.using_bd():
            if not bd_finished(model):
                resource: structures.Resource_remote_base \
                    = seekrflow.run_settings.get_resource_by_name(
                        seekrflow.run_settings.bd_stage_resource_name)
                if not transferred_files[seekrflow.run_settings.bd_stage_resource_name] \
                        and not no_transfer_before:
                    transfer_files_to_from_remote_resource(
                        seekrflow.name, resource, source_directory)
                    transferred_files[seekrflow.run_settings.hidr_stage_resource_name] = True

                destination_path = os.path.join(
                    resource.remote_working_directory, seekrflow.name)
                destination_model_filename = os.path.join(destination_path, "model.xml")
                if force_overwrite:
                    force_str = "-f "
                else:
                    force_str = ""
                command_string = f"python {resource.remote_seekr2_directory}/run.py "\
                    f"any_bd model.xml -N {bd_n_threads} {force_str}&> bd_run.out" 
                run_workflow, cancel_workflow, status_workflow, args \
                    = assign_remote_workflow(
                        seekrflow.name, destination_path, resource, usage_tracking, 
                        command_string, destination_model_filename, workflow_type="bd")
                if resource.remote_interface.type == "globus_compute_sdk":
                    endpoint = resource.remote_interface.endpoint_id
                    transfer_back = remote_globus\
                        .submit_remote_workflow_with_globus_compute(
                            run_workflow, cancel_workflow, status_workflow, endpoint, 
                            args)
                elif resource.remote_interface.type == "ssh":
                    hostname = resource.remote_interface.hostname
                    username = resource.remote_interface.username
                    password = resource.remote_interface.password
                    port = resource.remote_interface.port
                    #private_key_filename = resource.remote_interface.private_key_filename
                    #private_key_passphrase = resource.remote_interface.private_key_passphrase
                    transfer_back = remote_ssh\
                        .submit_remote_workflow_with_ssh(
                            run_workflow, cancel_workflow, status_workflow, args, 
                            hostname, username, password, port)
                else:
                    raise NotImplementedError(
                        f"Remote interface type not implemented: "\
                        f"{resource.remote_interface.type}")
                    
                if transfer_back:
                    transfer_files_to_from_remote_resource(
                        seekrflow.name, resource, source_directory, backwards=True)

        bd_running_locally = False

    if seekrflow.run_settings.hidr_stage_resource_name == "local":
        hidr_local_executor = assign_local_executor(
            label="hidr_executor",
            cores_per_worker=1)
        local_executors.append(hidr_local_executor)
        hidr_running_locally = True
    else:
        if not hidr_finished(model):
            resource: structures.Resource_remote_base \
                = seekrflow.run_settings.get_resource_by_name(
                    seekrflow.run_settings.hidr_stage_resource_name)
            if not transferred_files[seekrflow.run_settings.hidr_stage_resource_name] \
                    and not no_transfer_before:
                transfer_files_to_from_remote_resource(
                    seekrflow.name, resource, source_directory)
                transferred_files[seekrflow.run_settings.hidr_stage_resource_name] = True

            destination_path = os.path.join(
                resource.remote_working_directory, seekrflow.name)
            destination_model_filename = os.path.join(destination_path, "model.xml")
            input_pdb_filename = os.path.basename(seekrflow.starting_pdb_filename)
            if seekrflow.hidr_settings.type == "hidr_metaD":
                gaussian_height = seekrflow.hidr_settings.gaussian_height
                gaussian_width = seekrflow.hidr_settings.gaussian_width
                bias_factor = seekrflow.hidr_settings.bias_factor
                command_string \
                    = f"python {resource.remote_seekrtools_directory}/hidr/hidr.py any "\
                      f"model.xml -M metadyn -p {input_pdb_filename} "\
                      f"-H {gaussian_height} -w {gaussian_width} -b {bias_factor} -c 0 "\
                      f"&> hidr_run.out"

            elif seekrflow.hidr_settings.type == "hidr_SMD":
                restraint_force_constant = seekrflow.hidr_settings.restraint_force_constant
                translation_velocity = seekrflow.hidr_settings.translation_velocity
                command_string \
                    = f"python {resource.remote_seekrtools_directory}/hidr/hidr.py any"\
                      f" model.xml -M SMD -p {input_pdb_filename} "\
                      f"-k {restraint_force_constant} -v {translation_velocity} -c 0 "\
                      f"&> hidr_run.out"

            else:
                raise NotImplementedError(
                    f"HIDR type {seekrflow.hidr_settings.type} is not implemented.")
            run_workflow, cancel_workflow, status_workflow, args = assign_remote_workflow(
                seekrflow.name, destination_path, resource, usage_tracking, command_string,
                destination_model_filename, workflow_type="hidr")
            if resource.remote_interface.type == "globus_compute_sdk":
                endpoint = resource.remote_interface.endpoint_id
                transfer_back = remote_globus\
                    .submit_remote_workflow_with_globus_compute(
                        run_workflow, cancel_workflow, status_workflow, endpoint, 
                        args)
            elif resource.remote_interface.type == "ssh":
                hostname = resource.remote_interface.hostname
                username = resource.remote_interface.username
                password = resource.remote_interface.password
                port = resource.remote_interface.port
                #private_key_filename = resource.remote_interface.private_key_filename
                #private_key_passphrase = resource.remote_interface.private_key_passphrase
                transfer_back = remote_ssh\
                    .submit_remote_workflow_with_ssh(
                        run_workflow, cancel_workflow, status_workflow, args, 
                        hostname, username, password, port)
            else:
                raise NotImplementedError(
                    f"Remote interface type not implemented: "\
                    f"{resource.remote_interface.type}")
                
            if transfer_back:
                transfer_files_to_from_remote_resource(
                    seekrflow.name, resource, source_directory, backwards=True)
            # Re-load the model for SEEKR
            model = seekr2_base.load_model(model_filename)
            anchors_to_run = seekr_anchors_to_run(model)
            
        hidr_running_locally = False
    
    if seekrflow.run_settings.seekr_stage_resource_name == "local":
        seekr_local_executor = assign_local_executor(
            label="seekr_executor",
            cores_per_worker=1)
        local_executors.append(seekr_local_executor)
        seekr_running_locally = True
    else:
        if len(anchors_to_run) > 0:
            resource: structures.Resource_remote_base = seekrflow.run_settings\
                .get_resource_by_name(
                    seekrflow.run_settings.seekr_stage_resource_name)
            if not transferred_files[seekrflow.run_settings.seekr_stage_resource_name] \
                    and not no_transfer_before:
                transfer_files_to_from_remote_resource(
                    seekrflow.name, resource, source_directory)
                transferred_files[seekrflow.run_settings.seekr_stage_resource_name] = True

            destination_path = os.path.join(
                resource.remote_working_directory, seekrflow.name)
            destination_model_filename = os.path.join(destination_path, "model.xml")
            if force_overwrite:
                force_str = "-f "
            else:
                force_str = ""
            command_string = f"python {resource.remote_seekr2_directory}/run.py "\
                f"{{index}} model.xml {force_str} -c 0 &> seekr_{{index}}_run.out" 
            run_workflow, cancel_workflow, status_workflow, args = assign_remote_workflow(
                seekrflow.name, destination_path, resource, usage_tracking, command_string,
                destination_model_filename, workflow_type="seekr", indices=anchors_to_run)
            if resource.remote_interface.type == "globus_compute_sdk":
                endpoint = resource.remote_interface.endpoint_id
                transfer_back = remote_globus\
                    .submit_remote_workflow_with_globus_compute(
                        run_workflow, cancel_workflow, status_workflow, endpoint, 
                        args)
            elif resource.remote_interface.type == "ssh":
                hostname = resource.remote_interface.hostname
                username = resource.remote_interface.username
                password = resource.remote_interface.password
                port = resource.remote_interface.port
                #private_key_filename = resource.remote_interface.private_key_filename
                #private_key_passphrase = resource.remote_interface.private_key_passphrase
                transfer_back = remote_ssh\
                    .submit_remote_workflow_with_ssh(
                        run_workflow, cancel_workflow, status_workflow, args, 
                        hostname, username, password, port)
            else:
                raise NotImplementedError(
                    f"Remote interface type not implemented: "\
                    f"{resource.remote_interface.type}")

            if transfer_back:
                transfer_files_to_from_remote_resource(
                    seekrflow.name, resource, source_directory, backwards=True)

        seekr_running_locally = False
    
    if len(local_executors) > 0:
        multi_site_config = Config(
            executors=local_executors,
            strategy=None,
            usage_tracking=usage_tracking
        )

        parsl.load(multi_site_config)

        # Run BD - if needed
        if model.using_bd():
            if bd_running_locally and not bd_finished(model):
                bd_future = run_bd_locally(inputs=[model_filename, bd_n_threads])
                print("BD simulations launched")
                
        # Run HIDR if needed
        if hidr_running_locally and not hidr_finished(model):
            # TODO: if using metaD, SMD, ...
            input_pdb_file = os.path.basename(seekrflow.starting_pdb_filename)
            hidr_future = run_hidr_metaD_locally(inputs=[model_filename, input_pdb_file])
            # Wait until this is done to run SEEKR
            print("HIDR simulations launched")
            hidr_future.result()
            # Re-load the model for SEEKR
            model = seekr2_base.load_model(model_filename)
            anchors_to_run = seekr_anchors_to_run(model)

        # Run SEEKR if needed and HIDR is finished
        if seekr_running_locally and len(anchors_to_run) > 0:
            for anchor_index in anchors_to_run:
                seekr_future = run_seekr_locally(inputs=[model_filename, str(anchor_index)])
                print("anchor_index:", anchor_index)
                # TODO: this probably must be done for local job where there's just one GPU
                # but for the cluster, we will want to spin them off all at once.
                seekr_future.result()

        if model.using_bd():
            if bd_running_locally and not bd_finished(model):
                print("Waiting for BD simulations to finish...")
                bd_future.result()            

    os.chdir(curdir)
