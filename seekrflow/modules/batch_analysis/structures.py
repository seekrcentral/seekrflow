"""
Batch analysis objects.

The analyzers are run for each system in the batch, and their results are
collected and aggregated into a single result object.
"""

import os
import csv

from attrs import define, field, validators
import numpy as np
import scipy.stats as stats
from typing import List, Optional, Literal
import matplotlib.pyplot as plt

@define
class Data_point:
    """
    A data point: value, error, and unit, for either reference dataset 
    or batch results.
    """
    value: float = field(validator=validators.instance_of(float))
    error: float = field(validator=validators.instance_of(float))
    # TODO: use openmm.units?
    units: str = field(validator=validators.instance_of(str))

@define
class Reference_dataset:
    """
    Reference data that could be compared to the batch results. For example,
    experimental dissociation rate constants.
    """
    dataset_name: str = field(validator=validators.instance_of(str))
    field_name: str = field(validator=validators.instance_of(str))

@define
class Reference_dataset_from_csv_file(Reference_dataset):
    """
    Reference data that could be compared to the batch results. For example,
    experimental dissociation rate constants.
    """
    dataset_name: str = field(validator=validators.instance_of(str))
    field_name: str = field(validator=validators.instance_of(str))
    csv_filename: str = field(validator=validators.instance_of(str))
    column_name_value: str = field(validator=validators.instance_of(str))
    column_value_unit: str = field(validator=validators.instance_of(str))
    type: Literal["csv"] = "csv"
    column_name_error: str | None = field(default=None, validator=validators.optional(
        validators.instance_of(str)))
    rows: None | list[int] = field(
        default=None,
        validator=validators.optional(validators.deep_iterable(
            member_validator=validators.instance_of(int),
            iterable_validator=validators.instance_of(list),
        )),
    )
    _data_points: dict[int, Data_point] = field(
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(Data_point),
            iterable_validator=validators.instance_of(dict),
        ),
        factory=dict, eq=False, repr=False,
    )
    
    def __attrs_post_init__(self):
        self._data_points = {}
        with open(self.csv_filename, newline="") as file:
            reader = csv.DictReader(file)
            fieldnames = reader.fieldnames or []
            if self.column_name_value not in fieldnames:
                raise ValueError(
                    f"Column {self.column_name_value} not found in "
                    f"file {self.csv_filename}. Available columns: {fieldnames}.")
            if self.column_name_error and self.column_name_error not in fieldnames:
                raise ValueError(
                    f"Column {self.column_name_error} not found in "
                    f"file {self.csv_filename}. Available columns: {fieldnames}.")
            selected_rows = set(self.rows) if self.rows is not None else None
            for row_number, row in enumerate(reader, start=1):
                if selected_rows is not None and row_number not in selected_rows:
                    continue
                raw_value = row[self.column_name_value]
                if raw_value is None or raw_value.strip() == "":
                    continue
                if self.column_name_error:
                    raw_error = row[self.column_name_error]
                    error = float(raw_error) if raw_error and raw_error.strip() else 0.0
                else:
                    error = 0.0
                self._data_points[row_number] = Data_point(
                    value=float(raw_value),
                    error=error,
                    units=self.column_value_unit)

@define
class Scatter_plot:
    """
    A scatter plot batch analysis object comparing to a reference dataset.
    """
    title: str = field(validator=validators.instance_of(str))
    x_axis_label: str = field(validator=validators.instance_of(str))
    x_axis_unit: str = field(validator=validators.instance_of(str))
    y_axis_label: str = field(validator=validators.instance_of(str))
    y_axis_unit: str = field(validator=validators.instance_of(str))
    type: Literal["scatter"] = "scatter"
    scale: str = field(default="linear", validator=validators.in_(['linear', 'log']))

    def make_plot(
            self,
            data_points: dict[int, Data_point],
            reference_data_points: dict[int, Data_point],
            directory_full_path: str
            ) -> None:
        """
        Make the plot for this analysis series.
        """
        plot_basename = f"{self.title}.png"
        plot_full_path = os.path.join(directory_full_path, plot_basename)
        shared_indices = sorted(
            set(data_points) & set(reference_data_points))
        fig, ax = plt.subplots(figsize=(7, 7))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        for index in shared_indices:
            reference = reference_data_points[index]
            computed = data_points[index]
            ax.errorbar(
                reference.value, computed.value,
                xerr=reference.error, yerr=computed.error,
                fmt="o", color="#3c708f", ecolor="k",
                capsize=4, markersize=7)
            ax.annotate(
                str(index), (reference.value, computed.value),
                textcoords="offset points", xytext=(5, 5), fontsize=9)
        ax.set_xscale(self.scale)
        ax.set_yscale(self.scale)
        ax.set_xlabel(rf"{self.x_axis_label} ({self.x_axis_unit})", fontsize=16)
        ax.set_ylabel(rf"{self.y_axis_label} ({self.y_axis_unit})", fontsize=16)
        ax.set_title(self.title, fontsize=18)
        plt.tight_layout()
        plt.savefig(plot_full_path)
        plt.close()

@define
class Analysis_series:
    """
    A base class of data series to be analyzed: could be batch or reference.
    """

@define
class Analysis_batch_reference_series(Analysis_series):
    """
    A comparison of a batch series to a reference dataset.
    """
    analysis_field_name: str = field(validator=validators.instance_of(str))
    reference_dataset_name: str = field(validator=validators.instance_of(str))
    plots: list[Scatter_plot] = field(
        validator=validators.instance_of(list), factory=list)
    type: Literal["batch_reference"] = "batch_reference"
    _data_points: dict[int, Data_point] = field(
        validator=validators.deep_mapping(
            key_validator=validators.instance_of(int),
            value_validator=validators.instance_of(Data_point),
            mapping_validator=validators.instance_of(dict),
        ),
        factory=dict, eq=False, repr=False,
    )
    _reference_data_points: dict[int, Data_point] = field(
        validator=validators.deep_mapping(
            key_validator=validators.instance_of(int),
            value_validator=validators.instance_of(Data_point),
            mapping_validator=validators.instance_of(dict),
        ),
        factory=dict, eq=False, repr=False,
    )

    def fill_data_points(self, all_analysis_jsons, reference_dataset):
        self._data_points = {}
        for index, analysis_json in enumerate(all_analysis_jsons, start=1):
            raw = None
            for stage_key in sorted(analysis_json, key=lambda key: int(key)):
                results = analysis_json[stage_key].get("results", {})
                if self.analysis_field_name in results:
                    raw = results[self.analysis_field_name]
            if raw is None or raw.get("value") is None:
                continue
            error = raw.get("error")
            self._data_points[index] = Data_point(
                value=float(raw["value"]),
                error=float(error) if error is not None else 0.0,
                units=str(raw.get("units")))
        self._reference_data_points = reference_dataset._data_points
        
    def make_plots(
            self,
            directory_full_path: str
            ) -> None:
        """
        Make the plots for this analysis series.
        """
        for plot in self.plots:
            plot.make_plot(
                self._data_points, self._reference_data_points, 
                directory_full_path)
        
    def make_results(
            self,
            directory_full_path: str
            ) -> dict:
        """
        Make the results for this analysis series.

        Compute R, Spearman's rho, and Kendall's tau.
        """
        shared_indices = sorted(
            set(self._data_points) & set(self._reference_data_points))
        comp_values = []
        reference_values = []
        for index in shared_indices:
            reference = self._reference_data_points[index]
            computed = self._data_points[index]
            comp_values.append(computed.value)
            reference_values.append(reference.value)
        log_comp = np.log10(comp_values)
        log_exp = np.log10(reference_values)
        kendall_tau, _ = stats.kendalltau(log_comp, log_exp)
        spearman_rho, _ = stats.spearmanr(log_comp, log_exp)
        pearson_r, _ = stats.pearsonr(log_comp, log_exp)
        mae = np.mean(np.abs(log_comp - log_exp))
        results = {
            "kendall_tau": kendall_tau,
            "spearman_rho": spearman_rho,
            "pearson_r": pearson_r,
            "MeanAbsLog10Err": mae
        }
        return results

@define
class Bar_plot:
    """
    A bar plot batch analysis object comparing data points to each other.
    """
    title: str = field(validator=validators.instance_of(str))
    x_axis_label: str = field(validator=validators.instance_of(str))
    y_axis_label: str = field(validator=validators.instance_of(str))
    y_axis_unit: str = field(validator=validators.instance_of(str))
    type: Literal["bar"] = "bar"
    scale: str = field(default="linear", validator=validators.in_(['linear', 'log']))

    def make_plot(
            self,
            data_points: list[Data_point],
            directory_full_path: str
            ) -> None:
        """
        Make the plot for this analysis series.
        """
        plot_basename = f"{self.title}.png"
        plot_full_path = os.path.join(directory_full_path, plot_basename)
        sorted_indices = sorted(data_points.keys(), key=lambda x: data_points[x].value)
        value_list_sorted = [data_points[index].value for index in sorted_indices]
        error_list_sorted = [data_points[index].error for index in sorted_indices]
        width = 0.6
        capsize = 4
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        #sorted_dirs = [os.path.basename(os.path.normpath(dirs[index])) for index in sorted_indices]
        sorted_dirs = [str(index+1) for index in sorted_indices]
        bars1 = ax.bar(sorted_dirs, value_list_sorted, width=width,
                        yerr=error_list_sorted, color="#3c708f", ecolor="k", capsize=capsize,
                        label="Tau-RAMD") # color="tab:gray", color="#3c708f"
        ax.legend(prop={"size": 22})
        plt.xticks(rotation=25, fontsize=14, ha="right")
        plt.yticks(fontsize=14)
        plt.ylabel("$Tau$ (ps)", fontsize=24)
        plt.xlabel("System", fontsize=24)
        plt.title("Tau-RAMD Unbinding Times", fontsize=28)
        plt.yscale("log") #, nonpositive="mask")
        plt.tight_layout()
        plt.savefig(plot_full_path)
        plt.close()

@define
class Analysis_batch_series(Analysis_series):
    """
    An analysis of the batch series by itself.
    """
    analysis_field_name: str = field(validator=validators.instance_of(str))
    #data_series_names: Optional[List[str]] = field(
    #    default=None, validator=validators.optional(
    #    validators.instance_of(List[str])))
    plots: list[Bar_plot] = field(validator=validators.instance_of(
        list), factory=list)
    type: Literal["batch"] = "batch"
    _data_points: dict[int, Data_point] = field(
        validator=validators.deep_iterable(
            member_validator=validators.instance_of(Data_point),
            iterable_validator=validators.instance_of(dict),
        ),
        factory=dict, eq=False, repr=False,
    )

    def fill_data_points(self, all_analysis_jsons: list[dict]):
        """
        Fill the data points for this analysis series.
        """
        self._data_points = {}
        for i, analysis_json in enumerate(all_analysis_jsons):
            for stage_key in sorted(analysis_json.keys(), key=lambda x: int(x)):
                stage_json = analysis_json[stage_key]
                if self.analysis_field_name in stage_json["results"]:
                    data_point = stage_json["results"][self.analysis_field_name]
                    self._data_points[int(i)] = Data_point(
                        value=data_point["value"],
                        error=data_point["error"],
                        units=data_point["units"])

    def make_plots(
            self,
            directory_full_path: str
            ) -> None:
        """
        Make the plots for this analysis series.
        """
        plots = []
        for plot in self.plots:
            plot.make_plot(self._data_points, directory_full_path)
        return plots
    
    def make_results(
            self,
            directory_full_path: str
            ) -> dict:
        """
        Make the results for this analysis series.

        None for this analysis series type.
        """
        return {}
    
@define
class Batch_analysis:
    """
    A batch analysis object.
    """
    output_subdirectory: str = field(validator=validators.instance_of(str))
    reference_datasets: List[Reference_dataset] = field(
        validator=validators.instance_of(list), factory=list)
    analysis_series: List[Analysis_series] = field(
        validator=validators.instance_of(list),factory=list)
    
    def get_reference_dataset(self, reference_dataset_name: str) -> Reference_dataset:
        """
        Get a reference dataset by name.
        """
        for reference_dataset in self.reference_datasets:
            if reference_dataset.dataset_name == reference_dataset_name:
                return reference_dataset
        raise ValueError(f"Reference dataset {reference_dataset_name} not found.")