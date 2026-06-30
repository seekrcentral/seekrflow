"""
modules/transfer/globus.py

Handle file transfers via Globus, which can be useful for file transfers
where password entry and two-factor authentication would otherwise be
burdensome.
"""

GLOBUS_SEEKRFLOW_APP_CLIENT_ID = "683ab038-1578-4520-bfb4-57de7411102f"
GLOBUS_TRANSFER_RESOURCE_SERVER = "transfer.api.globus.org"

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
    from globus_sdk import TransferClient, TransferData, UserApp, \
        GlobusAppConfig

    # Ensure trailing slashes for directories
    if not local_path.endswith("/"):
        local_path += "/"
    if not remote_path.endswith("/"):
        remote_path += "/"

    # ========== TRANSFER CLIENT SETUP ==========
    # The GlobusApp framework handles token storage (in ~/.globus/app/),
    # the login flow, and authorization automatically. When a required
    # token is missing, it prints the login URL and prompts for the code.
    #
    # Force the command-line login flow so it uses the redirect URI that is
    # registered on the native ("Thick Client") app
    # (https://auth.globus.org/v2/web/auth-code). The default would
    # auto-select a local-server flow that redirects to http://localhost:<port>,
    # which is not a registered redirect URI and causes a
    # "Mismatching redirect URI" error.
    config = GlobusAppConfig(
        login_flow_manager="command-line",
        request_refresh_tokens=True,
    )
    app = UserApp(
        "seekrflow",
        client_id=GLOBUS_SEEKRFLOW_APP_CLIENT_ID,
        config=config,
    )
    tc = TransferClient(app=app)
    tc.add_app_data_access_scope(remote_collection_id)

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
        source_collection_id,
        destination_collection_id,
        label=f"{name} files transfer",
        sync_level="mtime",
        verify_checksum=True,
        encrypt_data=True,
    )
    tdata.add_item(source_path, destination_path, recursive=True)

    submit_result = tc.submit_transfer(tdata)
    print(f"Filetree transfer task submitted to Globus with ID: "
          f"{submit_result['task_id']}")
    while not tc.task_wait(submit_result['task_id'], timeout=60):
        pass

    print("Transfer complete!")
    return