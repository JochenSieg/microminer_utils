import argparse
import logging
import os
import sys
from pathlib import Path

import helper
from helper import BAD_PDBIDS
from helper import constants
from helper.constants import one_2_three_dict
from helper.data_operations import read_microminer_csv
from helper.utils import scantree

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="""
        Annotate a pre-defined mutation dataset with mutant structures searched by MicroMiner.
        This script just matches mutations in a mutation dataset like ProTherm with MicroMiner
        results. MicroMiner is not executed by this script.
        """
    )
    dataset_collection = helper.get_dataset_collection()
    supported_datasets = [
        dataset.name
        for dataset in dataset_collection
        if dataset_collection.is_mutation_dataset(dataset)
    ]
    parser.add_argument(
        "--dataset",
        "-d",
        required=True,
        type=str.lower,
        choices=supported_datasets,
        nargs="+",
        help="Data set name. Data location is inferred from config.ini",
    )
    parser.add_argument(
        "--mm_resultdir",
        "-m",
        required=True,
        type=str,
        help="Path to directory of MicroMiner results. Is searched recursively"
        "for MicroMiner result CSV files.",
    )
    parser.add_argument(
        "--outdir", "-o", default=os.getcwd(), type=str, help="Path to output directory"
    )
    parser.add_argument(
        "--ids_only",
        default=False,
        action="store_true",
        help="Only write PDB IDs of annotated structures instead of the"
        " whole mutation datasets with structure annotations.",
    )

    args = parser.parse_args()

    dataset_names = args.dataset
    mm_resultdir = Path(args.mm_resultdir)
    outdir = Path(args.outdir)
    ids_only = args.ids_only

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s"
    )
    logger.info(f'Starting scripts: {" ".join(sys.argv)}')

    if not mm_resultdir.is_dir():
        print(
            "Error: Specified MicroMiner result directory does not exist or is not a directory."
        )
        sys.exit(1)

    if not outdir.is_dir():
        print("Error: Specified output directory does not exist or is not a directory.")
        sys.exit(1)

    mm_result_file_paths = [
        path for path in scantree(mm_resultdir) if path.name == "resultStatistic.csv"
    ]
    df_mm = read_microminer_csv(mm_result_file_paths)

    logger.info(f"Collected {df_mm.shape[0]} MicroMiner hits from disk")

    b4 = df_mm.shape[0]

    # filter MM results similarity measures: drop all hits with too low global sequence identity
    # df_mm.drop(df_mm[df_mm["fullSeqId"] < 0.4].index, inplace=True)
    # drop exact duplicate hits if there are any.
    df_mm.drop_duplicates(
        inplace=True,
        subset=[
            constants.MM_QUERY_NAME,
            constants.MM_QUERY_POS,
            constants.MM_QUERY_AA,
            constants.MM_QUERY_CHAIN,
            constants.MM_HIT_NAME,
            constants.MM_HIT_POS,
            constants.MM_HIT_AA,
            constants.MM_HIT_CHAIN,
        ],
    )
    if df_mm.shape[0] != b4:
        logger.info(f"Filtered {b4 - df_mm.shape[0]} MicroMiner hits")
    assert (
        df_mm.drop_duplicates(
            subset=[
                constants.MM_QUERY_NAME,
                constants.MM_QUERY_POS,
                constants.MM_QUERY_CHAIN,
                constants.MM_HIT_NAME,
                constants.MM_HIT_POS,
                constants.MM_HIT_CHAIN,
            ]
        ).shape[0]
        == df_mm.shape[0]
    )

    mm_merge_cols = [
        constants.MM_QUERY_NAME,
        constants.MM_QUERY_AA,
        constants.MM_QUERY_POS,
        constants.MM_HIT_AA,
    ]
    mm_merge_cols_with_chain = mm_merge_cols + [constants.MM_QUERY_CHAIN]

    mut_merge_cols = [
        constants.WILD_COL,
        constants.WILD_AA,
        constants.WILD_SEQ_NUM,
        constants.MUT_AA,
    ]
    mut_merge_cols_with_chain = mut_merge_cols + [constants.WILD_CHAIN]

    dataset_collection = helper.get_dataset_collection()

    for dataset_name in dataset_names:
        logger.info(f"Annotating mutant structures for: {dataset_name}")

        dataset = dataset_collection.get_dataset(dataset_name)

        df_mut = dataset.read_single_mutations(pdb_mutant_only=False)

        # some mutation data set do not provide chain info in the wild-type protein
        left_on = mut_merge_cols
        right_on = mm_merge_cols
        if constants.WILD_CHAIN in df_mut.columns:
            left_on = mut_merge_cols_with_chain
            right_on = mm_merge_cols_with_chain

        # drop PDB Ids that are obsolete.
        df_mut = df_mut[~df_mut[helper.WILD_COL].isin(BAD_PDBIDS)]

        # convert mutation dataset 1-letter code to 3-letter code for merging
        df_mut[constants.WILD_AA] = df_mut[constants.WILD_AA].apply(
            lambda aa: one_2_three_dict[aa] if aa in one_2_three_dict else aa
        )
        df_mut[constants.MUT_AA] = df_mut[constants.MUT_AA].apply(
            lambda aa: one_2_three_dict[aa] if aa in one_2_three_dict else aa
        )

        # join mutation data with MicroMiner results.
        df_mut = df_mut.merge(df_mm, left_on=left_on, right_on=right_on, how="left")

        if ids_only:
            # only write structure annotations for the mutant without mutation data
            outfile_path = outdir / f"{dataset_name}_annotated_ids.tsv"
            df_mut[left_on + df_mm.columns.tolist()].to_csv(
                outfile_path, sep="\t", index=False, header=True
            )
        else:
            # write the full data sets with structure annotations for the mutant
            outfile_path = outdir / f"{dataset_name}_annotated.tsv"
            df_mut.to_csv(outfile_path, sep="\t", index=False, header=True)


if __name__ == "__main__":
    main()
