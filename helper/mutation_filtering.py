from typing import List, Set

import numpy as np
import pandas as pd

from .constants import WILD_AA, MUT_AA, physicochemical_groups, atom_count_grouping


def are_in_same_group(aa1: str, aa2: str, groups: List[Set[str]]) -> bool:
    """Test whether to amino acids are in the same group.

    :param aa1: First amino acid code.
    :param aa2: Second amino acid code.
    :param groups: Groups.
    :return: Treu if both amino acids are in the same group. False otherwise.
    """
    for g in groups:
        if aa1 in g:
            return aa2 in g
    return False


def mutations_are_silent(
    df: pd.DataFrame, aa_col1: str = WILD_AA, aa_col2: str = MUT_AA
) -> np.ndarray:
    """Get silent mutations.

    Returns a mask for the table rows. True values indicating mutation is silent
    a.k.a. mutant side chain has same properties as wild-type. False values mean
    a change there is a change of property.
    :param df: The table.
    :param aa_col1: Column name of first residue.
    :param aa_col2: Column name of second residue.
    :return: boolean mask indicating silent mutations.
    """

    is_silent_list = []
    for row in df.itertuples(index=False):
        is_silent_list.append(
            are_in_same_group(
                getattr(row, aa_col1), getattr(row, aa_col2), physicochemical_groups
            )
        )
    return np.array(is_silent_list)


def have_same_size(df: pd.DataFrame) -> list:
    is_diff = []
    for row in df.itertuples(index=False):
        is_diff.append(
            atom_count_grouping[getattr(row, WILD_AA)]
            != atom_count_grouping[getattr(row, MUT_AA)]
        )
    return is_diff


def is_proline_mutation(df: pd.DataFrame) -> list:
    """Return boolean mask indicating whether row contains a proline involving mutation.

    :param df: The table.
    :return: boolean mask indicating whether wild AA or mutant AA is proline.
    """
    return (df["queryAA"] == "PRO" or df["hitAA"] == "PRO").tolist()
