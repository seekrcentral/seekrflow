"""
modules/seekr_run.py

Using the seekr API as well as parsl, run stages of the seekr calculation.
"""

import os
import typing
from collections import defaultdict

import parsl
from parsl.app.app import python_app
from parsl.providers import LocalProvider, SlurmProvider
from parsl.launchers import SrunLauncher
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.addresses import address_by_hostname

import seekr2.modules.common_base as seekr2_base
import seekr2.run as seekr2_run

import seekrflow.modules.structures as structures

GLOBUS_SEEKRFLOW_APP_CLIENT_ID = "683ab038-1578-4520-bfb4-57de7411102f"
GLOBUS_TRANSFER_RESOURCE_SERVER = "transfer.api.globus.org"

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
    Check if the Hidr stage has finished. This is done by reading the model files
    and seeing if all anchors have starting PDB files.
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
def run_bd(inputs):
    import seekr2.modules.common_base as seekr2_base
    import seekr2.run as seekr2_run
    model_filename = inputs[0]
    n_threads = inputs[1]
    model = seekr2_base.load_model(model_filename)
    seekr2_run.run(model, "any_bd", n_threads=n_threads)
    return

@python_app(executors=["hidr_executor"])
def run_hidr_metaD(inputs):
    import os
    import seekr2.modules.common_base as seekr2_base
    import seekrtools.hidr.hidr as seekr2_hidr
    model_filename = inputs[0]
    input_pdb_file = inputs[1]
    model = seekr2_base.load_model(model_filename)
    os.chdir(os.path.abspath(model.anchor_rootdir))
    seekr2_hidr.hidr(model, "any", pdb_files=[input_pdb_file], mode="MetaD") #TODO: cuda_device_index?
    return model.openmm_settings.cuda_platform_settings.cuda_device_index

@python_app(executors=["seekr_executor"])
def run_seekr(inputs):
    import seekr2.modules.common_base as seekr2_base
    import seekr2.run as seekr2_run
    model_filename = inputs[0]
    anchor_index = inputs[1]
    model = seekr2_base.load_model(model_filename)
    seekr2_run.run(model, f"{anchor_index}")

def transfer_files_with_globus(
        name: str,
        local_path: str,
        remote_path: str,
        local_collection_id: str,
        remote_collection_id: str,
        backwards: bool = False,
        ) -> None:
    """
    Transfer files to remote system using Globus.
    """
    from globus_sdk import TransferClient, NativeAppAuthClient, TransferData, RefreshTokenAuthorizer, UserApp
    from globus_sdk.scopes import TransferScopes, MutableScope
    from globus_sdk.tokenstorage import SimpleJSONFileAdapter

    # Attempt to transfer files
    token_file = os.path.expanduser("~/.globus_tokens.json")
    transfer_scope = MutableScope(TransferScopes.all)
    
    # Setup token storage and auth client
    file_adapter = SimpleJSONFileAdapter(token_file)
    client = NativeAppAuthClient(GLOBUS_SEEKRFLOW_APP_CLIENT_ID)        
    client.oauth2_start_flow(requested_scopes=transfer_scope, refresh_tokens=True)

    # Get or refresh tokens
    tokens = file_adapter.get_token_data(GLOBUS_TRANSFER_RESOURCE_SERVER)
    if not tokens:
        authorize_url = client.oauth2_get_authorize_url()
        print(f"Please go to this URL and login: {authorize_url}")
        get_input = getattr(__builtins__, "raw_input", input)
        auth_code = get_input(
            "Please enter the code you get after login here: ").strip()
        token_response = client.oauth2_exchange_code_for_tokens(auth_code)
        file_adapter.store(token_response)
        tokens = token_response.by_resource_server[GLOBUS_TRANSFER_RESOURCE_SERVER]

    transfer_refresh_token = tokens['refresh_token']
    transfer_access_token = tokens['access_token']
    expires_at_sec = tokens['expires_at_seconds']

    # Ensure trailing slashes for directories
    if not local_path.endswith("/"):
        local_path += "/"
    if not remote_path.endswith("/"):
        remote_path += "/"

    # Create RefreshTokenAuthorizer
    authorizer = RefreshTokenAuthorizer(
        refresh_token=transfer_refresh_token,
        auth_client=client,
        access_token=transfer_access_token, 
        expires_at=expires_at_sec
    )

    # ========== TRANSFER CLIENT SETUP ==========
    app = UserApp("seekrflow", client_id=GLOBUS_SEEKRFLOW_APP_CLIENT_ID)
    tc = TransferClient(app=app).add_app_data_access_scope(remote_collection_id)

    if backwards:
        # Transferring from remote to local
        source_path = remote_path
        destination_path = local_path
        source_collection_id = remote_collection_id
        destination_collection_id = local_collection_id

    else:
        # Transferring from local to remote
        source_path = local_path
        destination_path = remote_path
        source_collection_id = local_collection_id
        destination_collection_id = remote_collection_id

    tdata = TransferData(
        tc,
        source_collection_id,
        destination_collection_id,
        label=f"{name} files transfer",
        sync_level="checksum",
        encrypt_data=True,
    )
    tdata.add_item(source_path, destination_path, recursive=True)

    submit_result = tc.submit_transfer(tdata)
    print(f"Filetree transfer task submitted to Globus with ID: {submit_result['task_id']}")
    while not tc.task_wait(submit_result['task_id'], timeout=60):
        pass

    print("Transfer complete!")
    return

def slurm_remote_workflow(args):
    import os
    import sys
    import parsl
    from parsl.app.app import bash_app
    from parsl.providers import SlurmProvider
    from parsl.launchers import SrunLauncher
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor

    def bd_finished_wkflow(
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

    def hidr_finished_wkflow(
            model: seekr2_base.Model
            ) -> bool:
        """
        Check if the Hidr stage has finished. This is done by reading the model files
        and seeing if all anchors have starting PDB files.
        """
        for alpha, anchor in enumerate(model.anchors):
            if anchor.bulkstate:
                continue
            anchor_pdb_filename = seekr2_base.get_anchor_pdb_filename(anchor)
            if anchor_pdb_filename == "":
                return False
            
        return True

    def seekr_anchors_to_run_wkflow(
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

    working_dir = args[0]
    max_workers_per_node = args[1]
    partition = args[2]
    account = args[3]
    constraint = args[4]
    nodes_per_block = args[5]
    cores_per_node = args[6]
    mem_per_node = args[7]
    init_blocks = args[8]
    max_blocks = args[9]
    time_limit = args[10]
    scheduler_options = args[11]
    worker_init = args[12]
    usage_tracking = args[13]
    command_string = args[14]
    indices = args[15]
    model_filename = args[16]
    workflow_type = args[17]

    my_executor = HighThroughputExecutor(
        label="my_executor",
        working_dir=working_dir,
        #worker_debug=False,
        max_workers_per_node=max_workers_per_node,
        provider=SlurmProvider(
            partition=partition,
            account=account,
            constraint=constraint,
            nodes_per_block=nodes_per_block,
            cores_per_node=cores_per_node,
            mem_per_node=mem_per_node,
            init_blocks=init_blocks,
            max_blocks=max_blocks, # TODO: make this configurable
            parallelism=1.0,
            walltime=time_limit,
            scheduler_options=scheduler_options,
            worker_init=worker_init,
            exclusive=False,
            launcher=SrunLauncher(),
            
        ),
    )

    parsl.clear()
    my_config = Config(
        executors=[my_executor],
        strategy="simple",
        usage_tracking=usage_tracking
    )

    parsl.load(my_config)

    # TODO: to be implemented in version 3.0
    #import seekr2.modules.common_base as seekr2_base
    #import seekr2.run as seekr2_run
    #model = seekr2_base.load_model(model_filename)

    @bash_app
    def run_workflow(inputs):
        cmd_str = inputs[0]
        return cmd_str
        
    if indices is None:
        job_already_finished = False
        # TODO: to be implemented in version 3.0
        #if workflow_type == "bd":
        #    job_already_finished = bd_finished_wkflow(model)
        #elif workflow_type == "hidr":
        #    job_already_finished = hidr_finished_wkflow(model)
        #elif workflow_type == "seekr":
        #    anchors_to_do = seekr_anchors_to_run_wkflow(model)
        #    if len(anchors_to_do) == 0:
        #        job_already_finished = True
        #else:
        #    raise NotImplementedError(
        #        f"Workflow type '{workflow_type}' is not implemented."
        #    )
        if not job_already_finished:
            full_command_string = f"""cd {working_dir};
        {command_string}"""
            future = run_workflow(inputs=[full_command_string])
            future.result()
    else:
        # TODO: to be implemented in version 3.0
        #if workflow_type == "seekr":
        #    anchors_to_do = seekr_anchors_to_run_wkflow(model)
        #else:
        #    raise Exception(
        #        f"Workflow type '{workflow_type}' is not implemented for indices."
        #    )
        futures = []
        for index in indices:
            #if index in anchors_to_do:
            full_command_string = f"""cd {working_dir};
    {command_string.format(index=index)}"""
            #future = run_workflow(inputs=[full_command_string])
            #futures.append(future)
            futures.append(run_workflow(inputs=[full_command_string]))

        #for future in futures:
        #    future.result()
        outputs = [i.result() for i in futures]

    return

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
        remote_path = os.path.join(resource.remote_working_directory, remote_root_directory_name)
        local_collection_id = resource.transfer_settings.local_collection_id
        remote_collection_id = resource.transfer_settings.remote_collection_id
        transfer_files_with_globus(remote_root_directory_name, local_directory, remote_path, 
                                    local_collection_id, remote_collection_id, backwards=backwards)
    else:
        raise NotImplementedError(
            "Only globus transfer is implemented.")
    return

def assign_remote_workflow(
        destination_path: str,
        resource: structures.Resource_remote_base,
        usage_tracking: int,
        command_string: str,
        model_filename: str,
        workflow_type: str,
        num_blocks: int,
        indices: list | None = None,
    ) -> typing.Any:
    
    if resource.init_blocks is None:
        init_blocks = num_blocks
    else:
        init_blocks = resource.init_blocks
    
    if resource.max_blocks is None:
        max_blocks = num_blocks
    else:
        max_blocks = resource.max_blocks

    if resource.type == "slurm_remote":
        slurm_args = [
            destination_path,
            resource.max_workers_per_node,
            resource.partition,
            resource.account,
            resource.constraint,
            resource.nodes_per_block,
            resource.cores_per_node,
            resource.memory_per_node,
            init_blocks,
            max_blocks,
            resource.time_limit,
            resource.scheduler_options,
            resource.worker_init,
            usage_tracking,
            command_string,
            indices,
            model_filename,
            workflow_type
        ]
        workflow = slurm_remote_workflow
    else:
        raise NotImplementedError(
            f"Resource type {resource.type} is not implemented.")
    
    return workflow, slurm_args

def submit_remote_workflow(
        workflow: typing.Any,
        endpoint: str,
        args: tuple,
):
    from globus_compute_sdk import Client, Executor
    from globus_compute_sdk.serialize import ComputeSerializer, CombinedCode
    with Executor(endpoint) as gcx:
        gcx.serializer = ComputeSerializer(strategy_code=CombinedCode())
        future = gcx.submit(workflow, args)
        future.result()

def run_model(
        seekrflow: structures.Seekrflow,
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
    # TODO: implement a argument entry for this.
    force_overwrite = False

    if seekrflow.bd_settings is not None:
        bd_n_threads = seekrflow.bd_settings.num_threads
    
    if seekrflow.run_settings.allow_parsl_usage_tracking:
        usage_tracking = 3
    else:
        usage_tracking = 0

    transferred_files = defaultdict(bool)
    local_executors = []

    anchors_to_run = seekr_anchors_to_run(model)

    if seekrflow.run_settings.bd_stage_resource_name == "local":
        bd_local_executor = assign_local_executor(
            label="bd_executor",
            cores_per_worker=bd_n_threads)
        local_executors.append(bd_local_executor)
        bd_running_locally = True
    else:
        if model.using_bd():
            if not bd_finished(model):
                resource: structures.Resource_remote_base = seekrflow.run_settings.get_resource_by_name(
                    seekrflow.run_settings.bd_stage_resource_name)
                if not transferred_files[seekrflow.run_settings.bd_stage_resource_name]:
                    transfer_files_to_from_remote_resource(seekrflow.name, resource, source_directory)
                    transferred_files[seekrflow.run_settings.hidr_stage_resource_name] = True

                destination_path = os.path.join(resource.remote_working_directory, seekrflow.name)
                destination_model_filename = os.path.join(destination_path, "model.xml")
                endpoint = resource.globus_compute_endpoint_id
                if force_overwrite:
                    force_str = "-f "
                else:
                    force_str = ""
                command_string = f"python {resource.remote_seekr2_directory}/run.py any_bd model.xml "\
                    f"-N {bd_n_threads} {force_str}&> bd_run.out" 
                workflow, args = assign_remote_workflow(destination_path, resource, usage_tracking, command_string,
                                                         destination_model_filename, workflow_type="bd", num_blocks=1)
                submit_remote_workflow(workflow, endpoint, args)
                transfer_files_to_from_remote_resource(seekrflow.name, resource, source_directory, backwards=True)
        
        bd_running_locally = False

    if seekrflow.run_settings.hidr_stage_resource_name == "local":
        hidr_local_executor = assign_local_executor(
            label="hidr_executor",
            cores_per_worker=1)
        local_executors.append(hidr_local_executor)
        hidr_running_locally = True
    else:
        if not hidr_finished(model):
            resource: structures.Resource_remote_base = seekrflow.run_settings.get_resource_by_name(
                seekrflow.run_settings.hidr_stage_resource_name)
            if not transferred_files[seekrflow.run_settings.hidr_stage_resource_name]:
                transfer_files_to_from_remote_resource(seekrflow.name, resource, source_directory)
                transferred_files[seekrflow.run_settings.hidr_stage_resource_name] = True

            destination_path = os.path.join(resource.remote_working_directory, seekrflow.name)
            destination_model_filename = os.path.join(destination_path, "model.xml")
            endpoint = resource.globus_compute_endpoint_id
            input_pdb_filename = os.path.basename(seekrflow.starting_pdb_filename)
            if seekrflow.hidr_settings.type == "hidr_metaD":
                gaussian_height = seekrflow.hidr_settings.gaussian_height
                gaussian_width = seekrflow.hidr_settings.gaussian_width
                bias_factor = seekrflow.hidr_settings.bias_factor
                # TODO: once seekrtools has the lastest changes pushed - restore this.
                #f"-M metadyn -p {input_pdb_filename} -w {gaussian_width} -H {gaussian_height} "\
                command_string = f"python {resource.remote_seekrtools_directory}/hidr/hidr.py any model.xml "\
                    f"-M metadyn -p {input_pdb_filename} -H {gaussian_height} "\
                    f"-b {bias_factor} -c 0 &> hidr_run.out" 
                
            elif seekrflow.hidr_settings.type == "hidr_SMD":
                restraint_force_constant = seekrflow.hidr_settings.restraint_force_constant
                translation_velocity = seekrflow.hidr_settings.translation_velocity
                command_string = f"python {resource.remote_seekrtools_directory}/hidr/hidr.py any model.xml "\
                    f"-M SMD -p {input_pdb_filename} -k {restraint_force_constant} -v {translation_velocity} "\
                    f"-c 0 &> hidr_run.out" 
                
            else:
                raise NotImplementedError(
                    f"HIDR type {seekrflow.hidr_settings.type} is not implemented.")
            workflow, args = assign_remote_workflow(destination_path, resource, usage_tracking, command_string,
                                                     destination_model_filename, workflow_type="hidr", num_blocks=1)
            submit_remote_workflow(workflow, endpoint, args)
            transfer_files_to_from_remote_resource(seekrflow.name, resource, source_directory, backwards=True)
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
            resource: structures.Resource_remote_base = seekrflow.run_settings.get_resource_by_name(
                seekrflow.run_settings.seekr_stage_resource_name)
            if not transferred_files[seekrflow.run_settings.seekr_stage_resource_name]:
                transfer_files_to_from_remote_resource(seekrflow.name, resource, source_directory)
                transferred_files[seekrflow.run_settings.seekr_stage_resource_name] = True

            destination_path = os.path.join(resource.remote_working_directory, seekrflow.name)
            destination_model_filename = os.path.join(destination_path, "model.xml")
            endpoint = resource.globus_compute_endpoint_id
            if force_overwrite:
                force_str = "-f "
            else:
                force_str = ""
            command_string = f"python {resource.remote_seekr2_directory}/run.py {{index}} model.xml "\
                f"{force_str} -c 0 &> seekr_{{index}}_run.out" 
            workflow, args = assign_remote_workflow(destination_path, resource, usage_tracking, command_string,
                                                     destination_model_filename, workflow_type="seekr", 
                                                     num_blocks=len(anchors_to_run), indices=anchors_to_run)
            submit_remote_workflow(workflow, endpoint, args)
            transfer_files_to_from_remote_resource(seekrflow.name, resource, source_directory, backwards=True)
        
        seekr_running_locally = False
    
    if len(local_executors) > 0:
        multi_site_config = Config(
            #executors=[bd_local_executor, hidr_local_executor, seekr_local_executor],
            executors=local_executors,
            strategy=None,
            usage_tracking=usage_tracking
        )

        parsl.load(multi_site_config)

        # Run BD - if needed
        if model.using_bd():
            if bd_running_locally and not bd_finished(model):
                bd_future = run_bd(inputs=[model_filename, bd_n_threads])
                print("BD simulations launched")
                
        # Run HIDR if needed
        if hidr_running_locally and not hidr_finished(model):
            # TODO: if using metaD, SMD, ...
            input_pdb_file = os.path.basename(seekrflow.starting_pdb_filename)
            hidr_future = run_hidr_metaD(inputs=[model_filename, input_pdb_file])
            # Wait until this is done to run SEEKR
            print("HIDR simulations launched")
            hidr_future.result()
            # Re-load the model for SEEKR
            model = seekr2_base.load_model(model_filename)
            anchors_to_run = seekr_anchors_to_run(model)

        # Run SEEKR if needed and HIDR is finished
        if seekr_running_locally and len(anchors_to_run) > 0:
            for anchor_index in anchors_to_run:
                seekr_future = run_seekr(inputs=[model_filename, str(anchor_index)])
                print("anchor_index:", anchor_index)
                # TODO: this probably must be done for local job where there's just one GPU
                # but for the cluster, we will want to spin them off all at once.
                seekr_future.result()

        if model.using_bd():
            if bd_running_locally and not bd_finished(model):
                print("Waiting for BD simulations to finish...")
                bd_future.result()            

    os.chdir(curdir)