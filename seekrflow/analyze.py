"""
Analyze batch data.
"""

import argparse

import seekrflow.modules.client.structures as client_structures
import seekrflow.modules.client.start as client_start
from seekrflow.modules.batch_analysis.structures import Batch_analysis
import seekrflow.modules.batch_analysis.analyze as batch_analyze

def analyze(batch: client_structures.Batch):
    """
    Analyze batch data.
    """
    batch_analyze.analyze_batch(batch)

def main():
    parser = argparse.ArgumentParser(
        description="Prepare a batch of seekrflow workflows")
    parser.add_argument(
        "input_json", metavar="INPUT_JSON", type=str,
        help="Input JSON file containing the job specification.")
    parser.add_argument(
        "-s", "--skip_checks", dest="skip_checks", default=False,
        help="By default, pre-simulation checks will be run after the "
        "preparation is complete, and if the checks fail, the SEEKR "
        "model will not be saved. This argument bypasses those "
        "checks and allows the model to be generated anyways.",
        action="store_true")

    args = vars(parser.parse_args())
    input_json = args.get("input_json")
    skip_checks = args.get("skip_checks")

    assert client_start.input_is_batch_file(input_json), \
        "Input JSON file must be a batch file."
    batch = client_structures.Batch.load_batch_file(input_json)
    assert batch.batch_analysis is not None, \
        "Batch analysis must be specified in the batch file."
    analyze(batch)

if __name__ == "__main__":
    main()