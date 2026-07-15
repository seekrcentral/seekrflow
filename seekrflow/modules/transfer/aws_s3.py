"""
modules/transfer/aws_s3.py

Handle file transfers via AWS S3 - the required method for Amazon web services.
"""

import os
import tempfile

def transfer_files_with_aws_s3(
        local_path: str,
        s3_uri: str,
        input_tarball_name: str,
        region: str,
        backwards: bool,
        ) -> None:
    """Tar the local model directory and upload it as a single S3 object.

    NOTE: 'aws s3 sync' silently drops empty directories because S3 has no
    real directories -- only objects (keys). The prepared SEEKR calculation
    relies on empty dirs (e.g. .../stage_1_equilibration/ where run.log is
    written), so we package the tree as a tarball, which preserves empty
    directories, and extract it inside the container.
    """
    if backwards:
        """Sync the (now updated) S3 prefix back down for local verification."""
        #print(f"\n=== Downloading results {S3_URI} -> {LOCAL_RESULTS_DIR} ===")
        os.makedirs(local_path, exist_ok=True)
        _run(["aws", "s3", "sync", S3_URI, local_path,
            "--region", region, "--only-show-errors"])
        # The input tarball lives in the same S3 prefix, so the sync pulls it back
        # down too. It is not a result -- drop it from the local results directory.
        downloaded_tarball = os.path.join(local_path, input_tarball_name)
        if os.path.isfile(downloaded_tarball):
            os.remove(downloaded_tarball)

    else:
        #print(f"\n=== Uploading {local_path} -> {s3_uri}/{input_tarball_name} ===")
        tarball_path = os.path.join(tempfile.gettempdir(), input_tarball_name)
        # -C so the archive contains the directory's contents at the top level
        # (model.json, anchor_*/, ...) rather than the local_path name itself.
        _run(["tar", "czf", tarball_path, "-C", local_path, "."])
        _run(["aws", "s3", "cp", tarball_path, f"{s3_uri}/{input_tarball_name}",
            "--region", region, "--only-show-errors"])
        os.remove(tarball_path)

    
