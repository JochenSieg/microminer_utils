"""
Compares mutation pairs in MicroMiner results with known mutations of a mutation dataset
and reports the number of correctly retrieved mutation structure pairs.
"""

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import List

import helper
from helper.data_operations import read_microminer_csv, merge_results_for_pair_eval

logger = logging.getLogger(__name__)


def run_mutation_checking(
    csv_input: List[Path], dataset_names: List[str], outdir: Path, backward: bool
) -> None:
    """Checks if known mutations are in MicroMiner output.

    :param csv_input: List of MicroMiner result files.
    :param dataset_names: Name of mutation data sets.
    :param outdir: Directory to write results.
    :param backward: Whether to consider wild-type in the mutation data sets as the query
                     or hit in the MicroMiner results.
    :return: None
    """
    # gather all 'resultStatistic.csv' in the csv_input list (including recursive read of dirs)
    files = []
    for path in csv_input:
        files.extend(
            [p for p in helper.utils.scantree(path) if p.name == "resultStatistic.csv"]
        )
    if len(files) == 0:
        print(
            "Error: No resultStatistic.csv in input (and not in subdirs of any input dir)."
        )
        sys.exit(1)
    logger.info(f"Gathered {len(files)} input resultStatistic.csv files.")

    df_res = read_microminer_csv(files)

    if df_res.shape[0] == 0:
        print("Error: resultStatistic.csv input files are empty.")
        sys.exit(1)
    logger.info(f"Gathered {df_res.shape[0]} input structure pairs.")

    merged_file_suffix = "_eval.tsv"
    report_file = outdir / "eval_report.txt"
    if backward:
        merged_file_suffix = "_eval_backward.tsv"
        report_file = outdir / "eval_report_backwards.txt"

    dataset_collection = helper.get_dataset_collection()

    for dataset_name in dataset_names:
        dataset = dataset_collection.get_dataset(dataset_name)

        logger.info(f"Evaluating known mutations of {dataset.name}")
        df_ref = dataset.read_single_mutations(pdb_mutant_only=True)

        df_anno, df_not_found, df_merged = merge_results_for_pair_eval(
            df_ref, df_res, backward=backward
        )
        logger.info(
            f"MicroMiner found {df_anno.shape[0]} of {df_merged.shape[0]} mutations in"
            f" {dataset_name}" + (" (backward)" if backward else "")
        )

        df_merged.to_csv(
            outdir / f"{dataset_name}{merged_file_suffix}",
            sep="\t",
            header=True,
            index=False,
        )
        not_found_report_file = outdir / f"{dataset_name}_not_found{merged_file_suffix}"
        df_not_found.to_csv(not_found_report_file, sep="\t", header=True, index=False)
        with open(report_file, "a+") as f:
            f.write(
                f"{dataset_name}{merged_file_suffix}"
                f"\t{df_anno.shape[0]}"
                f"\t{df_merged.shape[0]}"
                f"\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description="""
        Compares mutation pairs in MicroMiner results with known mutations of a mutation dataset.
        """
    )
    dataset_collection = helper.get_dataset_collection()
    supported_datasets = [
        dataset.name
        for dataset in dataset_collection.get_mutation_datasets_with_structure_pairs()
    ]
    parser.add_argument(
        "--csv",
        "-c",
        required=True,
        type=str,
        nargs="+",
        help="resultStatistic.csv file or dir in which resultStatistic.csv files"
        " will be searched recursively.",
    )
    parser.add_argument(
        "--dataset",
        "-d",
        required=True,
        type=str.lower,
        choices=supported_datasets,
        nargs="+",
        help="Data set name. Data location is inferred from config.ini. Will be"
        " used as ground truth of mutation structure pairs.",
    )
    parser.add_argument(
        "--outdir", "-o", default=os.getcwd(), type=str, help="Path to output directory"
    )
    parser.add_argument(
        "--backward",
        "-b",
        default=False,
        action="store_true",
        help="Invert wild-type and mutant",
    )

    args = parser.parse_args()

    csv_input = args.csv
    dataset_names = args.dataset
    outdir = Path(args.outdir)
    backward = args.backward

    if not outdir.is_dir():
        print("Error: Specified output directory does not exist or is not a directory.")
        sys.exit(1)

    logging.basicConfig(
        filename=str((outdir / "log.log").absolute()),
        level=logging.INFO,
        format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
    )
    logger.info(f'Starting scripts: {" ".join(sys.argv)}')

    csv_input = [Path(_) for _ in csv_input]
    run_mutation_checking(csv_input, dataset_names, outdir, backward)


if __name__ == "__main__":
    main()
