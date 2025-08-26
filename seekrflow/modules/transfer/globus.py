"""
modules/transfer/globus.py

Handle file transfers via Globus, which can be useful for file transfers
where password entry and two-factor authentication would otherwise be
burdensome.
"""

import os

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
    from globus_sdk import TransferClient, NativeAppAuthClient, \
        TransferData, RefreshTokenAuthorizer, UserApp
    from globus_sdk.scopes import TransferScopes, MutableScope
    from globus_sdk.tokenstorage import SimpleJSONFileAdapter

    # Attempt to transfer files
    token_file = os.path.expanduser("~/.globus_tokens.json")
    transfer_scope = MutableScope(TransferScopes.all)
    
    # Setup token storage and auth client
    file_adapter = SimpleJSONFileAdapter(token_file)
    client = NativeAppAuthClient(GLOBUS_SEEKRFLOW_APP_CLIENT_ID)        
    client.oauth2_start_flow(requested_scopes=transfer_scope, 
                             refresh_tokens=True)

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