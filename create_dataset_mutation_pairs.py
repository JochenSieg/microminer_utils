import os
from pathlib import Path
import argparse
import sys
import logging

import helper
from helper.data_operations import make_pair_parameter_table

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="""
        Creates a "dataset" CSV file of protein mutation structure pairs from pre-defined datasets. 
        """
    )
    dataset_collection = helper.get_dataset_collection()
    supported_datasets = [
        name
        for name in dataset_collection.get_dataset_names()
        if dataset_collection.get_dataset(name).has_structure_pairs
    ]

    parser.add_argument(
        "--dataset",
        "-d",
        required=False,
        type=str.lower,
        choices=supported_datasets,
        nargs="+",
        default=supported_datasets,
        help="Data set name. Data location is inferred from config.ini",
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

    dataset_names = args.dataset
    outdir = Path(args.outdir)
    backward = args.backward

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s"
    )
    logger.info(f'Starting scripts: {" ".join(sys.argv)}')

    if not outdir.is_dir():
        print("Error: Specified output directory does not exist or is not a directory.")
        sys.exit(1)

    for dataset_name in dataset_names:
        dataset = dataset_collection.get_dataset(dataset_name)

        logger.info(f"Creating pair datasets for: {dataset.name}")

        df = make_pair_parameter_table(dataset, backward)

        outfile_name = f"{dataset.name}_pair.tsv"
        if backward:
            outfile_name = f"{dataset.name}_pair_backward.tsv"
        df.to_csv(outdir / outfile_name, sep="\t", index=False)


if __name__ == "__main__":
    main()
