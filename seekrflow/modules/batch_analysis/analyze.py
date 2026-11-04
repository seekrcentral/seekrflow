"""
Analyze batch data.
"""

import os
import json

import seekr.analyze as seekr_analyze

import seekrflow.modules.client.structures as client_structures

def analyze_batch(batch: client_structures.Batch):
    """
    Analyze batch data.
    """
    os.makedirs(batch.batch_analysis.output_subdirectory, exist_ok=True)
    current_directory = os.path.abspath(os.getcwd())
    output_directory = os.path.abspath(batch.batch_analysis.output_subdirectory)
    seekrflow_objects = batch.create_seekrflow_objects()
    model_objects = batch.create_model_objects(seekrflow_objects)
    all_analysis_jsons = []
    for seekrflow, model in zip(
           seekrflow_objects, model_objects):
        os.chdir(seekrflow.get_root_directory())
        analysis_json = seekr_analyze.analyze_all(
            model, stage_arg=None, print_json=True)
        all_analysis_jsons.append(analysis_json)
        os.chdir(current_directory)
    
    all_results = {}
    for series in batch.batch_analysis.analysis_series:
        if series.type == "batch_reference":
            reference = batch.batch_analysis.get_reference_dataset(
                series.reference_dataset_name)
            series.fill_data_points(all_analysis_jsons, reference)
        else:
            series.fill_data_points(all_analysis_jsons)
        series.make_plots(output_directory)
        results = series.make_results(output_directory)
        all_results[series.type] = results
    
    with open(os.path.join(output_directory, "results.json"), "w") as file:
            json.dump(all_results, file, indent=4)