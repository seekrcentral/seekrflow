"""
modules/structures.py

Contain data structure classes used for seekrflow parameters/inputs.
"""

import os
import json
import glob
import typing
import pathlib
from shutil import copyfile

from attrs import define, field, validators, Factory
import cattrs
from cattrs.strategies import include_subclasses

import seekrflow.modules.base as base
import seekrflow.modules.workflows.protein_ligand_seekr2.structures \
    as protein_ligand_seekr2_structures
import seekrflow.modules.workflows.protein_ligand_membrane_seekr2.structures \
    as protein_ligand_membrane_seekr2_structures
import seekrflow.modules.parameterize_structures as parameterizer_structures

WORK = "work"
PARAMETERIZE = "parameterize"
ROOT = "root"
RUN = "run"

@define
class Remote_interface_base:
    """
    Base remote interface settings class.
    """
    type: typing.Literal["base"] = "base"

@define
class Remote_interface_ssh(Remote_interface_base):
    """
    SSH remote interface settings class.
    """
    type: typing.Literal["ssh"] = "ssh"
    hostname: str = field(
        default="",
        validator=validators.instance_of(str),
    )
    username: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
    )
    password: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
    )
    port: int = field(
        default=22,
        validator=validators.optional(validators.instance_of(int)),
    )
    private_key_filename: str = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
    )
    private_key_passphrase: str = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
    )

@define
class Remote_interface_globus_compute_sdk(Remote_interface_base):
    """
    Globus Compute SDK remote interface settings class.
    """
    type: typing.Literal["globus_compute_sdk"] = "globus_compute_sdk"
    endpoint_id: str = field(
        default="",
        validator=validators.instance_of(str),
    )

@define
class Transfer_settings_base:
    """
    Base transfer settings class.
    """
    type: typing.Literal["base"] = "base"

@define
class Transfer_settings_globus(Transfer_settings_base):
    """
    Globus transfer settings for transferring files.
    """
    type: typing.Literal["globus"] = "globus"
    local_collection_id: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    remote_collection_id: str = field(
        default="",
        validator=validators.instance_of(str),
        )

@define
class Transfer_settings_rsync(Transfer_settings_base):
    """
    Rsync transfer settings for transferring files.
    """
    type: typing.Literal["rsync"] = "rsync"
    hostname: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    username: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
        )
    password: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
        )
    port: int = field(
        default=22,
        validator=validators.instance_of(int)
        )
    private_key_filename: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
        )
    private_key_passphrase: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
        )

@define
class Resource_base:
    """
    Base class for resources.
    """
    type: typing.Literal["base"] = "base"
    

@define
class Resource_local(Resource_base):
    """
    Local resource for running the protocol.
    """
    pass

@define
class Resource_remote_base(Resource_base):
    """
    Base class for remote resources.
    """
    pass

@define
class Resource_remote_slurm(Resource_remote_base):
    """
    Slurm resource for running the protocol.
    """
    type: typing.Literal["slurm_remote"] = "slurm_remote"
    name: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    remote_seekr2_directory: str = field(
        default="$HOME/seekr2/seekr2/",
        validator=validators.instance_of(str),
        )
    remote_seekrtools_directory: str = field(
        default="$HOME/seekrtools/seekrtools/",
        validator=validators.instance_of(str),
        )
    remote_working_directory: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    partition: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    account: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    constraint: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
        )
    nodes_per_block: int = field(
        default=1,
        validator=validators.instance_of(int),
        )
    cores_per_node: int = field(
        default=1,
        validator=validators.instance_of(int),
        )
    memory_per_node: int = field(
        default=4,
        validator=validators.instance_of(int),
        )
    time_limit: str = field(
        default="00:30:00",
        validator=validators.instance_of(str),
        )
    scheduler_options: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    worker_init: str = field(
        default="",
        validator=validators.instance_of(str),
        )
    remote_interface: Remote_interface_globus_compute_sdk | Remote_interface_ssh = field(
        default=Factory(Remote_interface_globus_compute_sdk),
    )
    transfer_settings: Transfer_settings_rsync | Transfer_settings_globus = field(
        default=Factory(Transfer_settings_globus),
    )
    
@define
class Run_settings:
    """
    Settings for seekrflow runs.
    resource: the machine that a protocol will run on
    """
    resources: typing.List[Resource_remote_slurm] = field(
        default=Factory(list),)
    bd_stage_resource_name: str = field(
        default="local",
        validator=validators.instance_of(str),
        )
    hidr_stage_resource_name: str = field(
        default="local",
        validator=validators.instance_of(str),
        )
    seekr_stage_resource_name: str = field(
        default="local",
        validator=validators.instance_of(str),
        )
    
    def get_resource_by_name(
            self,
            resource_name: str
            ) -> Resource_base | None:
        """
        Get a resource by its name.
        """
        for resource in self.resources:
            if resource.name == resource_name:
                return resource
        raise ValueError(
            f"Resource with name '{resource_name}' not found in run_settings.resources.")

@define
class Seekrflow:
    """
    All the inputs and parameters needed for a seekrflow calculation.
    """
    name: str = field(
        default="my_name",
        validator=validators.instance_of(str),
        )
    # NOTE: structure version 1.0 had a globus_compute_endpoint_id directly
    #  within the slurm_remote resource object. This was deprecated in 1.1
    #  in order to support different remote_interfaces, such as SSH, by
    #  including a remote_interface attribute instead with the settings 
    #  specific to that remote resource. In addition, some attributes were
    #  removed from the resource objects since remote workflows no longer
    #  employ parsl. I did not attempt to implement backwards compatibility
    #  to structure v1.0 because seekrflow was not yet released.
    # NOTE: structure version 1.2 implemented an entirely different structure
    #  framework. Versions after 1.1 include a workflow object with many nested
    #  attribute objects. Additionally, several other structures are nested as
    #  attributes within the parameterize attribute. Backwards compatibility to
    #  structure version 1.1 was not implemented because seekrflow was not yet 
    #  released. Removed parsl-related attributes.
    structure_version: str = field(default="1.2",
                                 validator=validators.instance_of(str))
    workflow: protein_ligand_seekr2_structures.Protein_ligand_seekr2_workflow \
         | protein_ligand_membrane_seekr2_structures.Protein_ligand_membrane_seekr2_workflow \
        = field(
        default=Factory(protein_ligand_seekr2_structures.Protein_ligand_seekr2_workflow),
        validator=validators.instance_of(
            protein_ligand_seekr2_structures.Protein_ligand_seekr2_workflow \
            | protein_ligand_membrane_seekr2_structures.Protein_ligand_membrane_seekr2_workflow
            ),
        )
    physical_attributes: base.Physical_attributes = field(
        default=Factory(base.Physical_attributes),
        validator=validators.instance_of(base.Physical_attributes),
        )
    work_directory: str = field(
        default="work",
        validator=validators.instance_of(str),
        )
    parameterizer: parameterizer_structures.Parameterizer | None = field(
        default=Factory(parameterizer_structures.Parameterizer),
        validator=validators.optional(validators.instance_of(parameterizer_structures.Parameterizer)),
        )
    run_settings: Run_settings | None = None

    def save(
            self,
            filename: str
            ) -> None:
        """
        Save the Seekrflow object to a JSON file.
        """
        converter: cattrs.Converter = cattrs.Converter()
        # Make sure that interited data classes are unstructured as their
        #  subtypes.
        include_subclasses(protein_ligand_seekr2_structures.HIDR_settings_base, converter)
        include_subclasses(Transfer_settings_base, converter)
        include_subclasses(Resource_base, converter)
        seekrflow_dict: dict = converter.unstructure(self)
        json_dump: str = json.dumps(seekrflow_dict, indent=4)
        with open(filename, "w") as file:
            file.write(json_dump)
        return
    
    def make_work_directory(
            self,
            work_directory: pathlib.Path | str | None = None
            ) -> None:
        """
        Make the work directory for the Seekrflow calculation.
        """
        if work_directory is not None:
            self.work_directory = str(work_directory)
        os.makedirs(self.work_directory, exist_ok=True)
        return

    def get_work_directory(self) -> pathlib.Path:
        """
        Get the directory where the preparation files are stored.
        """
        work_dir = pathlib.Path(self.work_directory)
        #os.makedirs(work_dir, exist_ok=True)
        return work_dir

    def get_parameterize_directory(self) -> pathlib.Path:
        """
        Get the directory where the preparation files are stored.
        """
        param_dir = pathlib.Path(self.work_directory) / PARAMETERIZE
        os.makedirs(param_dir, exist_ok=True)
        return param_dir

    def get_root_directory(self) -> pathlib.Path:
        """
        Get the root directory where the Seekrflow files are stored.
        """
        root_dir = pathlib.Path(self.work_directory) / ROOT
        os.makedirs(root_dir, exist_ok=True)
        return root_dir
    
    def get_run_directory(self) -> pathlib.Path:
        """
        Get the directory where the run files are stored.
        """
        run_dir = pathlib.Path(self.work_directory) / RUN
        os.makedirs(run_dir, exist_ok=True)
        return run_dir

def load_seekrflow(
        filename: str
    ) -> Seekrflow:
    """
    Load a Seekrflow object from a JSON file.
    """
    with open(filename, "r") as file:
        json_string: str = json.load(file)
    converter: cattrs.Converter = cattrs.Converter()
    seekrflow_obj: Seekrflow = converter.structure(json_string, Seekrflow)
    return seekrflow_obj

def save_new_seekrflow(
        seekrflow: Seekrflow, 
        seekrflow_glob: str, 
        seekrflow_base: str, 
        save_old_seekrflow=True,
        directory = "."):
    """
    Generate a new seekrflow file. The old seekrflow file(s) will be renamed with a 
    numerical index.
    """
    
    model_path = os.path.join(directory, "seekrflow.json")
    if os.path.exists(model_path) and save_old_seekrflow:
        # This is expected, because this old model was loaded
        full_model_glob = os.path.join(directory, seekrflow_glob)
        num_globs = len(glob.glob(full_model_glob))
        new_pre_model_filename = seekrflow_base.format(num_globs)
        new_pre_model_path = os.path.join(directory, 
                                          new_pre_model_filename)
        print("Renaming model.xml to {}".format(new_pre_model_filename))
        copyfile(model_path, new_pre_model_path)
        
    print("Saving new seekrflow.json")
    seekrflow.save(model_path)
    return
