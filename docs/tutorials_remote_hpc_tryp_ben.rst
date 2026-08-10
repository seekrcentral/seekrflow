Running Trypsin-Benzamidine on a Remote HPC System
==================================================

This tutorial walks you through setting up and running our trypsin-benzamidine 
example system on a remote HPC system using Globus.

Prerequisites
-------------

- seekrflow installed and working (see :doc:`getting_started`) along with all required, recommended, and optional dependencies, including:
  - seekr
  - Browndye2
  - OpenEye Toolkits
  - PDBFixer
  - PDB2PQR
  - openmmforcefields
  - Globus Endpoints (optional)
  - Globus Compute SDK
- Access to a remote HPC system, where you can submit jobs and manage resources.

.. note::

    It is assumed that you have completed the 
    :doc:`Trypsin-Benzamidine System: Basic Protein-Ligand Setup<tutorials_basic_tryp_ben>` 
    tutorial for this to work properly, at least the parameterization step.

Step 3.1: Justification
-----------------------
Full seekr calculations almost always require a power GPU cluster or 
supercomputer, yet transferring files to and from a remote system, as well 
as managing job submissions with SLURM/PBS scripts can be slow and cumbersome. 
This tutorial shows how to use seekrflow to streamline this process, although
care must be taken to ensure that the remote system is configured correctly, 
and that the Globus endpoints are set up properly (unless one wants to use 
SSH to control one's remote job). So make sure that all dependences are 
installed on both the local and remote machines, and that the remote machine 
is set up as defined in :doc:`getting_started`.

Step 3.2: Prepare the Configuration JSON Files
----------------------------------------------
An example configuration file is provided in the directory, named 
```seekrflow_remote.json``` . 
If one opens this file, one will see a large section filled out
titled *run_settings*. This section contains settings, which are used to
transfer files to and from the remote system, as well as to submit jobs to 
the remote system, either 'delta' or 'anvil'. One will need to configure 
these settings to their own system - they will be different for everyone, and
it's impossible for me to anticipate the changes you will make, so you will 
need to be proactive and resourceful in order to get this tutorial completed. 
Yet, if you can complete this tutorial, you should be all set for running 
your own systems on HPC. Let us consider some of the settings that
one will likely need to modify in order to complete this tutorial.

- "type": This should be set to either "slurm_remote" or "pbs_remote", 
  depending on the job scheduler used by your HPC system.

- "name": Choose any name for your resource, which will be referenced by the 
  "placements" fields lower in the configuration file.

- "remote_working_directory": This is the directory on the remote system 
  where the seekr workflow will be copied into and run. Typically, HPC 
  resources have a "scratch", "work", or "projects" directory where intensive 
  read/write operations can be performed. Make sure that you have write 
  permissions to this directory.

- "partition": This is the partition on the remote system where the jobs will 
  be submitted. This is often something like "gpu" or "compute". Check with 
  your HPC documentation to find the correct partition name.
  
.. note::

    The attribute is called "partition" when using SLURM, but is renamed 
    "queue" when using PBS.

- "account": This is the account name that you were assigned for job 
  submissions on the remote system. You should check with any online 
  portal or HPC documentation to find the correct account name.

- "nodes_per_block": This is the number of nodes that will be requested per 
  "block", and will probably usually be kept at 1.

- "cpus_per_task": This is the number of cores that will be requested per 
  task. This should be set to the number of cores that you would like to 
  use for each job. If you use MPS, it should at least be equal to the number 
  of concurrent tasks running per GPU. Note that seekrflow is designed to 
  request shared resources, so this should not exceed the proportional 
  number of cores that you would like to use for a shared job (using a 
  single GPU, for instance). Consult your HPC documentation to find the 
  correct number of cores to request for your jobs.

- "memory_per_node": This is the amount of memory that will be requested 
  per node. This should be set to the amount of memory that you would like 
  to use for each job. Note that seekrflow is designed to request shared 
  resources, so this should not exceed the proportional amount of memory 
  that you would like to use for a shared job (using a single GPU, for 
  instance). Consult your HPC documentation to find the correct amount of 
  memory to request for your jobs.

- "time_limit": This is the maximum amount of time that the job will be 
  allowed to run on the remote system. This should be set according to your HPC
  documentation. Note that this is not the same as the total simulation time, 
  but rather the maximum time that the job will be allowed to run before it is 
  killed. Example: "time_limit": "24:00:00" would be 24 hours.

- "scheduler_options": These are settings that will be passed to the job 
  scheduler when submitting jobs. This can include things like job names, 
  output files, error files, etc. Most importantly, this line will probably 
  be used to assign GPU settings. Consult your HPC documentation to find 
  the correct settings for your system.

- "worker_init": These settings define which commands will be run upon 
  the creation of a new Globus compute "worker". This might include the 
  loading of important modules, setting environment variables, or
  activating a conda/mamba environment. This will depend on which HPC 
  resource is being used. One should consult the HPC documentation, and 
  probably experiment with debug/test job submissions in order to find the 
  correct settings for their system.

- "remote_interface":
  
  - "type": This should be set to "globus_compute_sdk" to use Globus 
    Compute for remote job handling. At present, Globus Compute, SSH, 
    and head_shell remote job handling are supported.
  - "endpoint_id": This is the Globus Compute endpoint ID that will be used 
    to submit jobs to the remote system. This should be set to the endpoint 
    ID that you created for your remote system. You can find this ID by 
    running ```globus-compute-endpoint list``` in the terminal of the 
    remote HPC resource.

- "transfer_settings":
  
  - "type": This should be set to "globus" to use Globus for file transfers. 
    At present, Globus, rsync, and head_shell transfers are supported.
  - "local_collection_id": The Globus collection (endpoint) ID for the local 
     system. One can find it by getting globus_connect_personal running on 
     one's own machine, and then using the Globus web portal "Collections" 
     page to find its UUID.
  - "remote_collection_id": The Globus collection (endpoint) ID for the remote 
     system. One can find it by searching for the Globus collection UUID of 
     the HPC resource in the Globus web portal "Collections" page.

- "placements": A list of placement objects. Each placement has a 
  ``target`` path (a list of procedure and/or role names, e.g. 
  ``["metadynamics", "logistic"]``), an optional ``resource`` name, 
  optional ``dispatch`` settings (``dimensions``, ``group_size``, 
  ``concurrency`` for array/MPS control), optional ``co_schedule_with`` 
  (``"predecessor"`` or ``"successor"`` to fold a stage into a neighbor's 
  allocation), and optional per-stage scheduler overrides. Time policies
  may also be set to control how much time is requested for each stage to run.
  Resource names should match the ``resources`` section  (use ``"local"`` 
  for the local computer). Longer matching target paths override shorter 
  ones. A placement with an empty target (``{"target": [], "resource": "local"}``) matches every
  stage and acts as the default.

.. note::

    It is also possible to use rsync to transfer files between local and 
    remote machines. In that case, one would use the "rsync" type within 
    the "transfer_settings" attribute. The *rsync* has many optional 
    attributes, including "hostname", "username", "password", "port", 
    "private_key_filename", and "private_key_passphrase", which can be 
    used to transfer files to and from the remote machine without using 
    Globus Endpoints. However, Globus can be helpful to avoid password and 
    two-factor authentication on certain HPC systems, which might make 
    using rsync cumbersome. However, if one has set up one's 
    ~/.ssh/config file automatically with usernames, passwords, and/or 
    private keys, then only the "hostname" alias will be needed within 
    the "transfer_settings" attribute.
    
.. note::

    Similarly to rsync, it is possible to use SSH to control your remote job,
    instead of Globus Compute. In that case, one would use the "ssh" type 
    within the "remote_interface" attribute. The *ssh* type, like rsync, 
    has many optional attributes, including "hostname", "username", "password",
    "port", "private_key_filename", and "private_key_passphrase", which can 
    be used to control one's remote job without using Globus Compute SDK. 
    However, Globus Compute SDK can be helpful to avoid password and 
    two-factor authentication on certain HPC systems, which might make using 
    SSH cumbersome. However, if one has set up one's ~/.ssh/config file 
    automatically with usernames, passwords, and/or private keys, then only 
    the "hostname" alias will be needed within the "transfer_settings" attribute.

Step 3.3: Run Seekrflow on a Remote HPC System
----------------------------------------------
Once these settings are configured, one can run the seekrflow workflow on 
the remote HPC system in the same way as before:

.. code-block:: bash

    mamba activate SEEKR
    python ~/seekrflow/seekrflow/flow.py prepare -i seekrflow_remote.json
    python ~/seekrflow/seekrflow/flow.py run -i seekrflow_remote.json
    python ~/seekr/seekr/analyze.py work_hpc/root/model.xml

The job will probably take quite a long time to run, depending on the 
resources available on the remote system, as well as the backlog in the 
remote job queue. However, the BD simulations should still be run remotely 
and synchronously with the rest of the jobs. In this configuration, HIDR 
will be run on the remote system first, and then the seekr anchor calculations 
will be run synchronously with each other. All file transfers should
be automatically handled to and from the remote resource.
