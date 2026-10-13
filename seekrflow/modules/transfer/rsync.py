"""
modules/transfer/rsync.py

Handle file transfers via rsync, which is a simpler, more basic
way to transfer files. However, this can be annoying if, say, two-
factor authentication is always required.
"""

import os
import six
import shlex
import fabric
import signal
import subprocess

_TRANSFER_PIDS: set[int] = set()

def cancel_running_transfers() -> None:
    """
    Cancel all running transfer processes.
    """
    for pid in list(_TRANSFER_PIDS):
        try:
            os.killpg(pid, signal.SIGTERM)
        except OSError:
            pass
    

class Connection(fabric.Connection):
    def rsync(self, source, target, exclude=(),  rsync_opts="", ssh_opts="",
              delete=False, strict_host_keys=True, backwards=False):
        if isinstance(exclude, six.string_types):
            exclude = [exclude]
        # Double-backslash-escape
        replace_bs = lambda s: str(s).replace('"', '\\\\"')
        exclusions = " ".join((f'--exclude "{replace_bs(s)}"' for s in exclude))
        # Honor SSH key(s)
        keys = self.connect_kwargs.get("key_filename", [])
        if isinstance(keys, six.string_types):
            keys = [keys]
        key_string = "-i " + " -i ".join(keys) if keys else ""
        port_string = f"-p {self.port}"

        # Workaround for our old keys:
        ssh_opts += '-o "PubkeyAcceptedKeyTypes +ssh-dss"'
        # Strict host key checking
        disable_keys = "-o StrictHostKeyChecking=no"
        if not strict_host_keys and disable_keys not in ssh_opts:
            ssh_opts += " {}".format(disable_keys)
        # Remote shell (SSH) options
        rsh_parts = [key_string, port_string, ssh_opts]
        rsh_string = "--rsh='ssh {}'".format(" ".join(rsh_parts)) \
            if any(rsh_parts) else ""
        # Set up options part of string
        options = f'{"--delete" if delete else ""}{exclusions} ' \
                  f'-pthrvz {rsync_opts} {rsh_string}'
        remote_prefix = f"{self.user}@{self.host}"
        if self.host.count(":") > 1:
            remote_prefix = f"[{remote_prefix}]"

        upload_direction = f"{source} {remote_prefix}:{target}" \
            if not backwards else f"{remote_prefix}:{source} {target}"
        #return self.local(f"rsync {options} {upload_direction}")
        """ # DEBUG
        cmd = f"rsync {options} {upload_direction}"
        proc = subprocess.Popen(cmd, shell=True, start_new_session=True)
        _TRANSFER_PIDS.add(proc.pid)
        try:
            rc = proc.wait()
        finally:
            _TRANSFER_PIDS.discard(proc.pid)
        if rc != 0:
            raise RuntimeError(f"rsync exited with {rc}")
        """
        cmd = f"rsync {options} {upload_direction}"
        print(f"rsync command: {cmd}")
        proc = subprocess.Popen(
            cmd, shell=True, start_new_session=True,
            stderr=subprocess.PIPE, text=True)
        _TRANSFER_PIDS.add(proc.pid)
        try:
            _, err = proc.communicate()
        finally:
            _TRANSFER_PIDS.discard(proc.pid)
        if proc.returncode != 0:
            raise RuntimeError(
                f"rsync exited with {proc.returncode}\n{err}")

def transfer_files_with_rsync(
        local_path: str, 
        remote_path: str, 
        remote_hostname: str, 
        remote_username: str | None = None, 
        remote_password: str | None = None, 
        port: int | None = 22, 
        private_key_filename: str | None = None, 
        private_key_passphrase: str | None = None, 
        backwards: bool = False
        ) -> None:
    
    print(f"Filetree transfer task submitted to host: '{remote_hostname}' with rsync")
    fabric_kwargs = {}
    if remote_password:
        fabric_kwargs['password'] = remote_password
    if private_key_filename:
        fabric_kwargs['key_filename'] = private_key_filename
    if private_key_passphrase:
        fabric_kwargs['passphrase'] = private_key_passphrase
    c = Connection(host=remote_hostname, user=remote_username, 
                   port=port, connect_kwargs=fabric_kwargs)
    if not local_path.endswith("/"):
        local_path += "/"
    if not remote_path.endswith("/"):
        remote_path += "/"
    # Ensure the remote destination directory exists before rsync upload.
    # DEBUG
    #c.run(f"mkdir -p {shlex.quote(remote_path)}", hide=True)
    print(
        f"ssh mkdir: {remote_username}@{remote_hostname}:{port} {remote_path}")
    try:
        c.run(f"mkdir -p {shlex.quote(remote_path)}", hide=True)
    except Exception as e:
        raise RuntimeError(
            f"ssh mkdir failed ({type(e).__name__}): {e}"
        ) from e
    
    if backwards:
        c.rsync(source=remote_path, target=local_path, rsync_opts="-uc", 
                backwards=True)

    else:
        c.rsync(source=local_path, target=remote_path, rsync_opts="-uc", 
                backwards=False)

    print("Transfer complete!")
    return