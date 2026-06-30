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
import seekrflow.modules.workflows.structures as workflow_structures
import seekrflow.modules.workflows.workflow as workflow_module
import seekrflow.modules.workflows.components as components_module
import seekrflow.modules.workflows.cv_specs as cv_specs_module
import seekrflow.modules.workflows.anchor_specs as anchor_specs_module
import seekrflow.modules.workflows.scale_settings as scale_settings_module
import seekrflow.modules.workflows.stage_procedures as stage_procedures_module
import seekrflow.modules.parameterize_workflow as parameterize_workflow_module
import seekrflow.modules.parameterize_structures as parameterizer_structures
import seekrflow.modules.parameters_topology_structures \
    as parameters_topology_structures
import seekrflow.modules.site_finder as site_finder


WORK = "work"
PARAMETERIZE = "parameterize"
ROOT = "root"
RUN = "run"


def _register_workflow_subclasses(converter: cattrs.Converter) -> None:
    """
    Register all polymorphic base classes used by the generic Workflow so that
    cattrs (un)structures them as their concrete subtypes.

    Registration order matters: ``include_subclasses`` eagerly bakes structure
    hooks for each subclass's fields at registration time, so any union that is
    NESTED inside another must be registered BEFORE its container. The order
    below therefore proceeds from the most deeply nested unions outward.
    """
    # Selectors are nested inside components.
    include_subclasses(components_module.Selector_base, converter)
    include_subclasses(components_module.Component_base, converter)
    include_subclasses(components_module.Group_selector_base, converter)
    include_subclasses(cv_specs_module.CV_spec_base, converter)
    include_subclasses(anchor_specs_module.Anchor_spec_base, converter)
    include_subclasses(scale_settings_module.Scale_settings_base, converter)
    # Sampling / restraint / completion specs are nested inside stage items,
    #  which are in turn nested inside stage procedures.
    include_subclasses(stage_procedures_module.Sampling_spec_base, converter)
    include_subclasses(stage_procedures_module.Restraint_spec_base, converter)
    include_subclasses(stage_procedures_module.Completion_spec_base, converter)
    include_subclasses(stage_procedures_module.Stage_item_base, converter)
    include_subclasses(stage_procedures_module.Stage_procedure_base, converter)
    return

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
    private_key_filename: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)),
    )
    private_key_passphrase: str | None = field(
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
    cpus_per_task: int = field(
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
    mps: int = field(
        default=1,
        validator=validators.instance_of(int),
    )
    remote_interface: Remote_interface_globus_compute_sdk | Remote_interface_ssh = field(
        default=Factory(Remote_interface_globus_compute_sdk),
    )
    transfer_settings: Transfer_settings_rsync | Transfer_settings_globus = field(
        default=Factory(Transfer_settings_globus),
    )

@define
class Resource_remote_pbs(Resource_remote_base):
    """
    PBS/Torque resource for running the protocol.
    """
    type: typing.Literal["pbs_remote"] = "pbs_remote"
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
    queue: str = field(  # PBS uses "queue" instead of "partition"
        default="",
        validator=validators.instance_of(str),
    )
    account: str = field(
        default="",
        validator=validators.instance_of(str),
    )
    resource_list: str | None = field(  # PBS -l options (e.g., "qos=high")
        default=None,
        validator=validators.optional(validators.instance_of(str)),
    )
    nodes_per_block: int = field(
        default=1,
        validator=validators.instance_of(int),
    )
    cpus_per_task: int = field(  # Keep same name as SLURM for consistency; convert to ppn in qsub
        default=1,
        validator=validators.instance_of(int),
    )
    memory_per_node: int = field(
        default=4,
        validator=validators.instance_of(int),
    )
    time_limit: str = field(  # Format: HH:MM:SS (walltime)
        default="00:30:00",
        validator=validators.instance_of(str),
    )
    scheduler_options: str = field(  # Additional PBS options
        default="",
        validator=validators.instance_of(str),
    )
    worker_init: str = field(
        default="",
        validator=validators.instance_of(str),
    )
    mps: int = field(
        default=1,
        validator=[validators.instance_of(int),
                            validators.ge(1)],
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
    resources: typing.List[
        Resource_remote_slurm | Resource_remote_pbs
    ] = field(
        default=Factory(list),)
    # Maps a user-authored PROCEDURE name to the name of the resource that runs
    #  it. Users author procedures (not stages) in the input, and each
    #  procedure expands into one or more auto-named stages; resources are
    #  therefore assigned per procedure and resolved to stages on demand via
    #  the procedure's ``get_stage_names``.
    procedure_resource_name: dict[str, str] = field(
        default=Factory(dict),
        validator=validators.instance_of(dict),
        )
    
    
    def get_resource_by_name(
            self,
            resource_name: str
            ) -> Resource_base | None:
        """
        Get a resource by its name.
        """
        if resource_name == "local":
            return None
        for resource in self.resources:
            if resource.name == resource_name:
                return resource
        raise ValueError(
            f"Resource with name '{resource_name}' not found in run_settings.resources.")

    def get_procedure_resource(
            self,
            procedure_name: str,
            ) -> Resource_base | None:
        """
        Get the resource assigned to a given procedure (None means local).
        """
        if procedure_name not in self.procedure_resource_name:
            raise ValueError(
                f"No resource assigned to procedure '{procedure_name}'.")
        resource_name = self.procedure_resource_name[procedure_name]
        return self.get_resource_by_name(resource_name)

    def set_procedure_resource(
            self,
            procedure_name: str,
            resource_name: str,
            ) -> None:
        """
        Assign a resource name to a procedure after validating that the
        resource exists (or is local).
        """
        # Validate resource name early for clearer CLI/runtime errors.
        self.get_resource_by_name(resource_name)
        self.procedure_resource_name[procedure_name] = resource_name
        return

    def get_stage_resource_name(
            self,
            stage_name: str,
            procedure: stage_procedures_module.Stage_procedure_base,
            ) -> str:
        """
        Resolve the resource name configured for a given stage by finding the
        procedure that produces the stage and returning its assigned resource.

        ``procedure`` is the workflow's top-level stage procedure (usually a
        ``Composite_stage_procedure`` whose children are the user-authored
        procedures).
        """
        stage_to_procedure = stage_procedures_module\
            .build_stage_to_procedure_map(procedure)
        if stage_name not in stage_to_procedure:
            raise ValueError(
                f"Stage '{stage_name}' is not produced by any procedure.")
        procedure_name = stage_to_procedure[stage_name]
        if procedure_name not in self.procedure_resource_name:
            raise ValueError(
                f"No resource assigned to procedure '{procedure_name}' "
                f"(which produces stage '{stage_name}').")
        return self.procedure_resource_name[procedure_name]

    def get_stage_resource(
            self,
            stage_name: str,
            procedure: stage_procedures_module.Stage_procedure_base,
            ) -> Resource_base | None:
        """
        Get the resource for a given stage by resolving the procedure that
        produces it (None means local).
        """
        resource_name = self.get_stage_resource_name(stage_name, procedure)
        return self.get_resource_by_name(resource_name)

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
    # NOTE: structure version 1.4 implemented the generalized, composable
    #  Workflow (single prepare-side Workflow object replacing the per-system
    #  workflow classes) and split the parameterize-side state into a separate
    #  parameterize_workflow object. Added random_seed to physical_attributes.
    # Backwards compatibility to earlier versions was not implemented because 
    # not needed.
    structure_version: str = field(default="1.4",
                                 validator=validators.instance_of(str))
    workflow: workflow_module.Workflow = field(
        default=Factory(workflow_module.Workflow),
        validator=validators.instance_of(workflow_module.Workflow),
        )
    physical_attributes: base.Physical_attributes = field(
        default=Factory(base.Physical_attributes),
        validator=validators.instance_of(base.Physical_attributes),
        )
    work_directory: str | None = field(
        default="work",
        validator=validators.optional(validators.instance_of(str))
    )
    root_directory: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str))
    )
    # TODO: resolve these two
    parameterize_workflow: \
        parameterize_workflow_module.Parameterize_workflow | None = field(
        default=Factory(parameterize_workflow_module.Parameterize_workflow),
        validator=validators.optional(validators.instance_of(
            parameterize_workflow_module.Parameterize_workflow)),
        )
    parameterizer: parameterizer_structures.Parameterizer | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(parameterizer_structures.Parameterizer)),
    )
    run_settings: Run_settings = field(
        default=Factory(Run_settings),
        validator=validators.instance_of(Run_settings),
    )

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
        _register_workflow_subclasses(converter)
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
    
    def handle_ligand_indices(
            self,
            ligand_indices: str,
            ligand_resname: str,
            pdb_filename: str
            ) -> None:
        """
        Handle the ligand indices string and set the ligand_indices attribute.
        """
        if self.parameterize_workflow.has_small_molecule_ligand():
            # If the ligand resname is not provided, use the one from the seekrflow input
            if ligand_resname == "":
                ligand_resname = self.parameterize_workflow.parameterizer_information.ligand_resname
            else:
                self.parameterize_workflow.parameterizer_information.ligand_resname = ligand_resname
            
            # if the ligand indices are provided, use them preferentially
            if ligand_indices != "":
                ligand_indices = base.initialize_ref_indices(ligand_indices)
            else:
                if len(self.parameterize_workflow.ligand_indices) > 0:
                    ligand_indices = self.parameterize_workflow.ligand_indices
                elif ligand_resname != "":
                    ligand_indices = base.get_ligand_indices(pdb_filename, ligand_resname)
                else:
                    # TODO: implement some automated way to identify the ligand molecule
                    # in a molecular complex?
                    ligand_indices = []
            if len(ligand_indices) > 0:
                self.parameterize_workflow.ligand_indices = ligand_indices
            else:
                if len(self.parameterize_workflow.ligand_indices) == 0:
                    if self.parameterize_workflow.parameterizer_information.ligand_resname != "":
                        self.parameterize_workflow.ligand_indices = base.get_ligand_indices(
                            self.parameterize_workflow.parameterizer_information\
                                .receptor_ligand_pdb_filename, 
                            self.parameterize_workflow.parameterizer_information.ligand_resname)
                    else:
                        raise Exception("No ligand indices provided and no ligand "
                                        "residue name specified.")
        return
    
    def handle_receptor_indices(
            self,
            pdb_filename: str
            ) -> None:
        """
        Handle the receptor indices string and set the receptor_indices attribute.
        """
        if self.parameterize_workflow.has_small_molecule_ligand():
            ligand_indices = self.parameterize_workflow.ligand_indices
            assert len(ligand_indices) > 0, \
                "Cannot determine receptor indices without ligand " \
                "indices."
            receptor_indices = site_finder.site_finder_monte_carlo(
                pdb_filename,
                ligand_indices
            )
            self.parameterize_workflow.receptor_indices = receptor_indices

        return
    
    def assign_seekrflow_parameter_topology_files(
            self,
            parameter_topology_files: list[str]
            ) -> None:
        """
        Assign parameter and topology files to seekrflow structure.
        """
        assert parameter_topology_files is not None
        # First, we must try to discern the types of parameter files provided
        file_extensions = [os.path.splitext(f)[-1].lower() for f in parameter_topology_files]
        if ".prmtop" in file_extensions or ".parm7" in file_extensions:
            # Then assume AMBER - find the prmtop or parm7 file and assign it
            for f in parameter_topology_files:
                ext = os.path.splitext(f)[-1].lower()
                if ext == ".prmtop" or ext == ".parm7":
                    self.parameterize_workflow.solvated_system_for_md.parameters_topology \
                        = parameters_topology_structures.Amber_parameters_topology()
                    self.parameterize_workflow.solvated_system_for_md.parameters_topology\
                        .prmtop_filename = f
                    break
        elif ".psf" in file_extensions:
            # Then assume CHARMM - find the psf file and assign it
            other_parameter_files = []
            for f in parameter_topology_files:
                ext = os.path.splitext(f)[-1].lower()
                if ext == ".psf":
                    self.parameterize_workflow.solvated_system_for_md.parameters_topology \
                        = parameters_topology_structures.Charmm_parameters_topology()
                    self.parameterize_workflow.solvated_system_for_md.parameters_topology\
                        .psf_filename = f
                else:
                    other_parameter_files.append(f)
            self.parameterize_workflow.solvated_system_for_md.parameters_topology\
                .param_filename_list = other_parameter_files
        elif ".top" in file_extensions or ".gro" in file_extensions:
            raise NotImplementedError(
                "GROMACS parameter and topology files are not yet supported.")
        elif ".xml" in file_extensions:
            # Then assume OpenMM XML file
            for f in parameter_topology_files:
                ext = os.path.splitext(f)[-1].lower()
                if ext == ".xml":
                    self.parameterize_workflow.solvated_system_for_md.parameters_topology \
                        = parameters_topology_structures.Openmm_system()
                    self.parameterize_workflow.solvated_system_for_md.parameters_topology\
                        .system_filename = f
                    break
        else:
            raise ValueError(
                "Could not discern parameter and topology file types from "
                "provided files. Supported types include AMBER (.prmtop, .parm7), "
                "CHARMM (.psf), GROMACS (.top, .itp), and OpenMM (.xml).")
        return

def try_to_load_resources_json(
        seekrflow: "Seekrflow",
        ) -> None:
    """
    Try to load a resources JSON file and add the resources to the seekrflow
    run_settings.
    """
    home_seekrflow_resources_filename = os.path.expanduser(
        "~/.seekrflow_resources.json")
    if seekrflow.work_directory is not None:
        work_seekrflow_resources_filename = os.path.join(
            seekrflow.work_directory, ".seekrflow_resources.json")
    else:
        work_seekrflow_resources_filename = None
    if seekrflow.root_directory is not None:
        root_seekrflow_resources_filename = os.path.join(
            seekrflow.root_directory, ".seekrflow_resources.json")
    else:
        root_seekrflow_resources_filename = None

    json_str: dict | None = None
    if root_seekrflow_resources_filename is not None:
        if os.path.exists(root_seekrflow_resources_filename):
            with open(root_seekrflow_resources_filename, "r") as file:
                json_str = json.load(file)

    if json_str is None and work_seekrflow_resources_filename is not None:
        if os.path.exists(work_seekrflow_resources_filename):
            with open(work_seekrflow_resources_filename, "r") as file:
                json_str = json.load(file)

    if json_str is None:
        if os.path.exists(home_seekrflow_resources_filename):
            with open(home_seekrflow_resources_filename, "r") as file:
                json_str = json.load(file)

    if json_str is None:
        return

    converter: cattrs.Converter = cattrs.Converter()
    include_subclasses(Resource_base, converter)
    run_settings_obj: Run_settings = converter.structure(json_str, Run_settings)
    if seekrflow.run_settings is None:
        seekrflow.run_settings = Run_settings()
    for resource in run_settings_obj.resources:
        if resource.name in [r.name for r in seekrflow.run_settings.resources]:
            continue
        seekrflow.run_settings.resources.append(resource)
    return

def load_seekrflow(
        filename: str
    ) -> Seekrflow:
    """
    Load a Seekrflow object from a JSON file.
    """
    with open(filename, "r") as file:
        json_string: str = json.load(file)
    converter: cattrs.Converter = cattrs.Converter()
    _register_workflow_subclasses(converter)
    include_subclasses(Transfer_settings_base, converter)
    include_subclasses(Resource_base, converter)
    seekrflow_obj: Seekrflow = converter.structure(json_string, Seekrflow)
    try_to_load_resources_json(seekrflow_obj)
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

def assign_default_prepare_settings(
        seekrflow: Seekrflow
        ) -> None:
    """
    Assign default parameterize-side settings to the seekrflow structure.
    """
    seekrflow.parameterize_workflow = \
        parameterize_workflow_module.Parameterize_workflow()
    seekrflow.parameterize_workflow.parameterizer_information = \
        parameterize_workflow_module.Parameterizer_information()
    seekrflow.parameterize_workflow.md_settings = \
        workflow_structures.MD_settings()
    seekrflow.parameterize_workflow.bd_settings = \
        workflow_structures.BD_settings()
    seekrflow.physical_attributes = base.Physical_attributes()
    return
