import logging
from pathlib import Path
from typing import List, Tuple

import pandas as pd
from pandas.errors import EmptyDataError

import helper
from helper import BAD_PDBIDS
from helper.constants import (
    MM_QUERY_NAME,
    MM_HIT_NAME,
    MM_QUERY_POS,
    MM_HIT_POS,
    MM_QUERY_AA,
    MM_HIT_AA,
    WILD_COL,
    MUTANT_COL,
    one_2_three_dict,
    WILD_AA,
    MUT_AA,
    WILD_SEQ_NUM,
    MM_QUERY_CHAIN,
    WILD_CHAIN,
    MM_HIT_CHAIN,
)
from helper.datasets.dataset import Dataset
from helper.datasets.scope import read_scope
from helper.datasets.utils import get_pdb_file_path

logger = logging.getLogger(__name__)
dataset_collection = helper.get_dataset_collection()


def make_search_parameter_table(dataset: Dataset, backward: bool) -> pd.DataFrame:
    """Generates a DataFrame for MicroMiner search from a dataset

    :param dataset: The dataset instance.
    :param backward: Whether to do a backward search of a mutation. Does only work for datasets
                     with structures for both wild-type and mutant.
    :return: Dataframe containing two columns: an identifier and filepath to structure file.
    """
    pdb_mirror = "standard"
    if dataset.has_custom_structure_files:
        pdb_mirror = dataset.name

    df = None
    id_col = None
    path_col = "structure_path"
    if dataset_collection.is_mutation_dataset(dataset):
        # pdb_mutant_only must be False. Can be set to True for faster dev cycles
        df = dataset.read_single_mutations(pdb_mutant_only=False)

        # drop PDB Ids that are Calpha only or obsolete. We do not need those for mutation
        # evaluation experiments
        df = df[~df[helper.WILD_COL].isin(BAD_PDBIDS)]

        if backward:
            if helper.MUTANT_COL not in df.columns:
                raise ValueError(
                    f"backward mode not possible. Dataset has no"
                    f"mutant structure annotation: {dataset.name}"
                )

            pdb_mirror = "standard"  # backward mode always with standard PDB files.
            df = df[~df[helper.MUTANT_COL].isin(BAD_PDBIDS)]

            df[path_col] = df[helper.MUTANT_COL].apply(
                lambda val: get_pdb_file_path(
                    val, mirror=pdb_mirror, allow_obsolete=True
                )
            )
            id_col = helper.MUTANT_COL
        else:
            df[path_col] = df[helper.WILD_COL].apply(
                lambda val: get_pdb_file_path(
                    val, mirror=pdb_mirror, allow_obsolete=True
                )
            )
            id_col = helper.WILD_COL
    else:
        if dataset.name == "scope":
            df = read_scope()
            id_col = "sid"
        elif dataset.name == "pdb" or dataset.name == "afdb":
            df = dataset.read()
            id_col = "id"
        else:
            raise NotImplementedError("dataset not known")
        if path_col not in df.columns:
            df[path_col] = df[id_col].apply(
                lambda val: get_pdb_file_path(
                    val, mirror=pdb_mirror, allow_obsolete=True
                )
            )

    df = df[[id_col, path_col]]

    na_mask = df.isna().any(axis=1)
    if na_mask.any():
        # remove missing values, e.g. when only CIF file is available, but we want PDB
        logger.info(f"Dropping {na_mask.sum()} missing values (of {df.shape[0]} total)")
        df = df[na_mask]

    # remove duplicate ROWS, i.e. when ID and PATH are identical. Would lead to redundant
    # computations. Duplicates are frequent in mutation datasets, e.g. the same mutation
    # was measured multiple times experimentally under different conditions.
    df_out = df.drop_duplicates()

    df = df_out.rename(columns={id_col: "id"})

    return df


def make_pair_parameter_table(dataset: Dataset, backward: bool) -> pd.DataFrame:
    """Generates a DataFrame for MicroMiner pair mode from a dataset.

    Note that the dataset must contain structures for both wild-type and mutant.

    :param dataset: The dataset
    :param backward: Whether to do a backward search of a mutation.
    :return: Dataframe containing two columns: an identifier and filepath to structure file.
    """
    if not dataset.has_structure_pairs:
        raise ValueError(
            f"dataset has no structures for both wild-type and mutant: {dataset.name}"
        )

    pdb_mirror1 = "standard"
    pdb_mirror2 = "standard"
    if dataset.has_custom_structure_files:
        pdb_mirror1 = dataset.name

    path_col1 = "structure_path1"
    path_col2 = "structure_path2"

    id_col1 = helper.WILD_COL
    id_col2 = helper.MUTANT_COL

    if not dataset_collection.is_mutation_dataset(dataset):
        raise NotImplementedError("")

    df = dataset.read_single_mutations(pdb_mutant_only=True)

    df = df[~df[helper.WILD_COL].isin(BAD_PDBIDS)]
    df = df[~df[helper.MUTANT_COL].isin(BAD_PDBIDS)]

    if backward:
        if helper.MUTANT_COL not in df.columns:
            raise ValueError(
                f"backward mode not possible. Dataset has no"
                f"mutant structure annotation: {dataset.name}"
            )
        # backward mode currently always with standard PDB files.
        pdb_mirror2, pdb_mirror1 = pdb_mirror1, pdb_mirror2
        id_col1, id_col2 = id_col2, id_col1

    df[path_col1] = df[id_col1].apply(
        lambda val: get_pdb_file_path(val, mirror=pdb_mirror1, allow_obsolete=True)
    )
    df[path_col2] = df[id_col2].apply(
        lambda val: get_pdb_file_path(val, mirror=pdb_mirror2, allow_obsolete=True)
    )
    df = df[[id_col1, path_col1, id_col2, path_col2]]

    na_mask = df.isna().any(axis=1)
    if na_mask.any():
        # remove missing values, e.g. when only CIF file is available, but we want PDB
        logger.info(f"Dropping {na_mask.sum()} missing values (of {df.shape[0]} total)")
        df = df[na_mask]

    # remove duplicate ROWS, i.e. when ID and PATH are identical. Would lead to redundant
    # computations.
    df = df.drop_duplicates()

    df = df.rename(columns={id_col1: "id1", id_col2: "id2"})
    return df


def read_microminer_csv(files: List[Path]) -> pd.DataFrame:
    """Reads result CSVs of a MicroMiner to a single dataframe.

     This function is convenient because sometimes we need to enforce dtypes.
     For example PDB IDs are sometimes interpreted as a float given in scientific notation
     or residue positions can or can not contain insertion code (iCode) or start with a minus
     or chain identifiers are '1' and interpreted as int.
    :param files: List of file paths to MicroMiner CSV files.
    :return: A single dataframe containing the content of all input CSV files
            (duplicate entries are removed).
    """
    # We use strings to guarantee consistent and correct behaviour when matching.
    col_dtypes = {
        MM_QUERY_NAME: str,
        MM_HIT_NAME: str,
        MM_QUERY_POS: str,
        MM_HIT_POS: str,
        MM_QUERY_CHAIN: str,
        MM_HIT_CHAIN: str,
    }

    def df_gen(files_list):
        for filepath in files_list:
            try:
                yield pd.read_csv(filepath, sep="\t", header=0, dtype=col_dtypes)
            except EmptyDataError:
                print("Warning: file is empty: ", filepath)

    return pd.concat(df_gen(files), ignore_index=True).drop_duplicates()


def merge_results_for_pair_eval(
    df_dataset: pd.DataFrame, df_microminer: pd.DataFrame, backward: bool
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Merges a mutation dataset table with a MicroMiner result table for evaluation.

    Wild-type columns are matched with MM query while mutant name and mutant AA are matched with
    MM hit columns to check if MicroMiner could retrieve the wild-type/mutant structure pair
    from the mutation dataset. Use backward parameter to invert matching.

    # If a chain ID is available in the mutation dataset it will be used for matching. Note that,
    # many mutation dataset contain no chain or inconsistent chain information (especially for
    # the backward case) in these cases chain ID will not be used for matching.

    :param df_dataset: The mutation dataset dataframe.
    :param df_microminer: The MicroMiner result dataframe.
    :param backward: Whether to match backward (match wild-type columns with MM hit columns
                     and mutant columns with MM query columns)
    :return: 3 dataframes: mutations found by MicroMiner, mutations not found, a dataframe of all
             mutation annotated with MicroMiner results.
    """
    df_ref = df_dataset
    df_mm = df_microminer

    assert MUTANT_COL in df_ref.columns, f"Need {MUTANT_COL} column"

    right_on = [MM_QUERY_NAME, MM_HIT_NAME, MM_QUERY_AA, MM_HIT_AA, MM_QUERY_POS]
    if backward:
        right_on = [MM_QUERY_NAME, MM_HIT_NAME, MM_QUERY_AA, MM_HIT_AA, MM_HIT_POS]

    ref_key_cols = [WILD_COL, WILD_AA, WILD_SEQ_NUM, MUT_AA, MUTANT_COL]

    if WILD_CHAIN in df_ref.columns:
        # Unfortunately, some mutation dataset hold no or inconsistent chain ID annotations.
        # We only use the chain ID if we trust the annotation. We trust the annotation
        # when a WILD_CHAIN column is available.
        right_on.append(MM_HIT_CHAIN if backward else MM_QUERY_CHAIN)
        ref_key_cols.append(WILD_CHAIN)

    # We can drop duplicates in MicroMiner results here. Duplicates are not necessary
    # to check for successful retrieval of a known wild-type/mutant structure pair.
    # In addition, probably all duplicates present are multiple matches in homo-meric structures.
    df_mm = df_mm.drop_duplicates(right_on)

    df_ref = df_ref[~df_ref[WILD_COL].isin(BAD_PDBIDS)]
    df_ref = df_ref[~df_ref[MUTANT_COL].isin(BAD_PDBIDS)]

    # drop duplicates. There are often duplicates, e.g. because of multiple ddG measurements
    df_ref = df_ref.drop_duplicates(ref_key_cols)

    from_one_to_three = lambda aa: one_2_three_dict[aa] if len(aa) != 3 else aa
    df_ref["wild_aa3"] = df_ref[WILD_AA].apply(from_one_to_three)
    df_ref["mutant_aa3"] = df_ref[MUT_AA].apply(from_one_to_three)

    left_on = [WILD_COL, MUTANT_COL, "wild_aa3", "mutant_aa3", WILD_SEQ_NUM]
    if backward:
        left_on = [MUTANT_COL, WILD_COL, "mutant_aa3", "wild_aa3", WILD_SEQ_NUM]

    if WILD_CHAIN in df_ref.columns:
        left_on.append(WILD_CHAIN)

    df_merged = df_ref.merge(df_mm, left_on=left_on, right_on=right_on, how="left")
    df_anno = df_merged.dropna(subset=right_on)
    df_not_found = df_ref.merge(
        df_mm, left_on=left_on, right_on=right_on, how="left", indicator=True
    )
    df_not_found.query('_merge == "left_only"', inplace=True)
    df_not_found.drop("_merge", axis=1, inplace=True)

    # clean up
    df_anno = df_anno.drop(["wild_aa3", "mutant_aa3"], axis=1)
    df_not_found.drop(["wild_aa3", "mutant_aa3"], axis=1, inplace=True)
    df_merged = df_merged.drop(["wild_aa3", "mutant_aa3"], axis=1)

    return df_anno, df_not_found, df_merged
