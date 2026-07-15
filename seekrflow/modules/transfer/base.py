"""
modules/transfer/base.py

Base module for transfer methods in seekrflow - common functions to 
all transfer methods.
"""

import os

from seekrflow.modules import structures
import seekrflow.modules.transfer.globus as transfer_globus
import seekrflow.modules.transfer.rsync as transfer_rsync
import seekrflow.modules.transfer.aws_s3 as transfer_aws_s3

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
        
    elif resource.transfer_settings.type == "rsync":
        remote_path = os.path.join(
            resource.remote_working_directory, remote_root_directory_name)
        remote_hostname = resource.transfer_settings.hostname
        remote_username = resource.transfer_settings.username
        remote_password = resource.transfer_settings.password
        port = resource.transfer_settings.port
        private_key_filename = resource.transfer_settings.private_key_filename
        private_key_passphrase = resource.transfer_settings.private_key_passphrase
        transfer_rsync.transfer_files_with_rsync(
            local_directory, remote_path, remote_hostname, remote_username, 
            remote_password, port, private_key_filename=private_key_filename,
            private_key_passphrase=private_key_passphrase, 
            backwards=backwards)
    
    elif resource.transfer_settings.type == "aws_s3":
        s3_uri = resource.transfer_settings.get_uri()
        input_tarball_name = resource.transfer_settings.get_input_tarball_name()
        transfer_aws_s3.transfer_files_with_aws_s3(
            local_directory, s3_uri, input_tarball_name, resource.region, 
            backwards=backwards)

    else:
        raise NotImplementedError(
            "Only rsync and globus transfers are implemented.")
    return