import argparse
import logging
import os
import sys
from pathlib import Path

import pandas as pd

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
        Compute statistics for the structure annotation calculated by MicroMiner.
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

    args = parser.parse_args()

    dataset_names = args.dataset
    mm_resultdir = Path(args.mm_resultdir)
    outdir = Path(args.outdir)

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

    logger.info(f"Collected {df_mm.shape[0]} MicroMiner rows from disk")

    b4 = df_mm.shape[0]

    # filter MM results similarity measures: drop all hits with too low global sequence identity
    df_mm.drop(df_mm[df_mm["fullSeqId"] < 0.4].index, inplace=True)
    # drop exact duplicates
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
        logger.info(f"Filtered {b4 - df_mm.shape[0]} rows")
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

    # read problematic mutations for statistics
    col_dtypes_problematic = {
        "pdb_wild": str,
        "pdb_mutant": str,
        "wild_aa": str,
        "mutation_aa": str,
        "seq_num": str,
    }
    df_problematic = pd.read_csv(
        helper.CONFIG["DATA"]["PROBLEMATIC_SINGLE_MUTATIONS"],
        sep="\t",
        dtype=col_dtypes_problematic,
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
    report_data = {
        "dataset": [],
        "original_muts_mutants": [],
        "corrected_muts_mutants": [],
        "mm_muts_mutants": [],
        "all_muts_wild_pdb": [],
    }
    for dataset_name in dataset_names:
        logger.info(
            f"Calculating statistics for mutant structure annotations for: {dataset_name}"
        )

        dataset = dataset_collection.get_dataset(dataset_name)

        df_mut = dataset.read_single_mutations(pdb_mutant_only=False)

        # some mutation data set do not provide chain info in the wild-type protein
        left_on = mut_merge_cols
        right_on = mm_merge_cols
        if constants.WILD_CHAIN in df_mut.columns:
            left_on = mut_merge_cols_with_chain
            right_on = mm_merge_cols_with_chain

        # drop PDB Ids are obsolete.
        df_mut = df_mut[~df_mut[helper.WILD_COL].isin(BAD_PDBIDS)]

        # convert mutation dataset 1-letter code to 3-letter code for merging
        df_mut[constants.WILD_AA] = df_mut[constants.WILD_AA].apply(
            lambda aa: one_2_three_dict[aa] if aa in one_2_three_dict else aa
        )
        df_mut[constants.MUT_AA] = df_mut[constants.MUT_AA].apply(
            lambda aa: one_2_three_dict[aa] if aa in one_2_three_dict else aa
        )

        # left join will likely increase the number of rows in df_mut since often
        # there are multiple structure hits for the same mutation row.
        df_mut = df_mut.merge(df_mm, left_on=left_on, right_on=right_on, how="left")

        # CALC STATISTICS:

        # remove entries where the same mutation is measured multiple times
        df_mut_u = df_mut.drop_duplicates(left_on)
        nof_all_wild_pdb_mutations = df_mut_u.shape[0]

        # calc dataset mutation with at least one MM hit
        nof_annotated_mutations = (~df_mut_u[constants.MM_QUERY_NAME].isnull()).sum()

        nof_mutant_structures_original = 0
        nof_mutant_structures_original_corrected = 0
        # calc number of mutations with mutant structure originally in the mutation dataset
        df_mut_pdbmut = dataset.read_single_mutations(pdb_mutant_only=True)
        if constants.MUTANT_COL in df_mut_pdbmut.columns:
            df_mut_pdbmut = df_mut_pdbmut[
                ~df_mut_pdbmut[constants.WILD_COL].isin(BAD_PDBIDS)
            ]
            df_mut_pdbmut = df_mut_pdbmut[
                ~df_mut_pdbmut[constants.MUTANT_COL].isin(BAD_PDBIDS)
            ]

            left_on_mutant_structures = left_on + [constants.MUTANT_COL]
            nof_mutant_structures_original = df_mut_pdbmut.drop_duplicates(
                left_on_mutant_structures
            ).shape[0]
            # correct original mutation count
            key_mut = [
                constants.WILD_COL,
                constants.MUTANT_COL,
                constants.WILD_AA,
                constants.MUT_AA,
                constants.WILD_SEQ_NUM,
            ]
            key_problematic = [
                "pdb_wild",
                "pdb_mutant",
                "wild_aa",
                "mutation_aa",
                "seq_num",
            ]
            if constants.WILD_CHAIN in df_mut_pdbmut.columns:
                key_mut = key_mut + [constants.WILD_CHAIN]
                key_problematic = key_problematic + ["wild_chain"]

            df_merged = df_mut_pdbmut.merge(
                df_problematic,
                left_on=key_mut,
                right_on=key_problematic,
                how="left",
                indicator=True,
            )
            df_mut_pdbmut_corrected = df_merged[df_merged["_merge"] == "left_only"]
            nof_mutant_structures_original_corrected = (
                df_mut_pdbmut_corrected.drop_duplicates(
                    left_on_mutant_structures
                ).shape[0]
            )

        # write statistic numbers
        logger.info(
            f"{dataset_name}: mutations with mutant structure: "
            f"original={nof_mutant_structures_original}"
            f" corrected={nof_mutant_structures_original_corrected}"
            f" MicroMiner={nof_annotated_mutations}"
            f"; all_mutation_wild_pdb={nof_all_wild_pdb_mutations}"
        )
        report_data["dataset"].append(dataset_name)
        report_data["original_muts_mutants"].append(nof_mutant_structures_original)
        report_data["corrected_muts_mutants"].append(
            nof_mutant_structures_original_corrected
        )
        report_data["mm_muts_mutants"].append(nof_annotated_mutations)
        report_data["all_muts_wild_pdb"].append(nof_all_wild_pdb_mutations)

        report_file = outdir / "annotation_statistics.txt"
        df_report = pd.DataFrame(report_data)
        df_report.to_csv(report_file, sep="\t", header=True, index=False)


if __name__ == "__main__":
    main()
