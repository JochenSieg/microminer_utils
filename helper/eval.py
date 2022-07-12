import pandas as pd
from pathlib import Path

from . import constants
from .utils import get_switched_column_names


def check_row_in_result(
        df: pd.DataFrame, row, backward) -> bool:
    """

    :param df: The resultStatistic table from the mutation search tool.
    :param row:
    :param backward:
    :return:
    """

    in_aa = getattr(row, constants.WILD_AA)
    out_aa = getattr(row, constants.MUT_AA)

    in_aa_3letter = constants.one_2_three_dict[in_aa] if len(in_aa) != 3 else in_aa
    out_aa_3letter = constants.one_2_three_dict[out_aa] if len(out_aa) != 3 else out_aa
    seq_num_int = int(getattr(row, constants.WILD_SEQ_NUM))

    query = getattr(row, constants.WILD_COL)
    mutant = getattr(row, constants.MUTANT_COL)

    # print('looking for ', query, mutant, in_aa_3letter, seq_num_int, out_aa_3letter)

    # if backward:
    #     # in the backward case just swap values
    #     in_aa_3letter, out_aa_3letter = out_aa_3letter, in_aa_3letter
    #     query, mutant = mutant, query

    # run select query to check if protherm mutation was found
    df_this = df.query(f'{constants.MICROMINER_QUERY_NAME} == @query and \
                        {constants.MICROMINER_HIT_NAME} == @mutant and \
                        {constants.MICROMINER_QUERY_AA} == @in_aa_3letter and \
                        {constants.MICROMINER_HIT_AA} == @out_aa_3letter')

    if backward:
        # seq_num usually describes ONlY the position in the wild_pdb.
        # For the reverse case we therefore do not know the seq_num for the wild type (the original mutant).
        # Therefore, we just check if the mutant seq num is the seq num we expect.
        df_this = df_this.query(f'{constants.MICROMINER_HIT_POS} == @seq_num_int')
    else:
        df_this = df_this.query(f'{constants.MICROMINER_QUERY_POS} == @seq_num_int')

    return df_this.shape[0] > 0


def run_evaluation(df: pd.DataFrame, experiment_policy, outdir: Path, backward: bool) -> list:
    """
    Deals with the computation of the mutation search experiment.
    :param df: Mutation data table.
    :param experiment_policy: Policy to handle experiment computation.
    :param outdir: Directory to write experiment output to.
    :param backward: Whether use mutant PDB as search query instead of wild-type.
    :return: bool mask indicating whether each row in df was found by the mutation search tool.
    """
    if backward:
        query_col, target_col = constants.MUTANT_COL, constants.WILD_COL
    else:
        target_col, query_col = constants.MUTANT_COL, constants.WILD_COL

    was_found = [None for _ in range(df.shape[0])]

    # use grouping to avoid duplicate computation.
    for name, df_group in df.groupby(experiment_policy.get_groupby_cols(query_col, target_col)):

        # get PDB ids
        query_pdbid = df_group[query_col].iloc[0]
        target_pdbid = df_group[target_col].iloc[0]

        this_outdir = outdir / query_pdbid
        this_outdir.mkdir(parents=True, exist_ok=True)

        # call the abstract function that executes the mutation search tool with some set of arguments
        experiment_policy.compute(query_pdbid, target_pdbid, this_outdir)

        results_csv_path = this_outdir / 'resultStatistic.csv'
        if not results_csv_path.is_file():
            raise FileNotFoundError('resultStatistic.csv was not generated! '
                                    'Something wrong with the executable?')

        df_resultStatistic = pd.read_csv(results_csv_path, sep='\t', header=0)

        for row in df_group.itertuples():
            row_pos_in_table = df.index.get_loc(row.Index)
            assert type(row_pos_in_table) == int  # there is something wrong if this is not int
            was_found[row_pos_in_table] = check_row_in_result(
                df_resultStatistic,
                row,
                backward)

    return was_found


def run_evaluation2(df: pd.DataFrame, experiment_policy, outdir: Path, backward: bool) -> list:
    """
    Deals with the computation of the mutation search experiment.
    :param df: Mutation data table.
    :param experiment_policy: Policy to handle experiment computation.
    :param outdir: Directory to write experiment output to.
    :param backward: Whether use mutant PDB as search query instead of wild-type.
    :return: bool mask indicating whether each row in df was found by the mutation search tool.
    """
    if backward:
        # exchange wild-type and mutant column
        df.columns = get_switched_column_names(df, constants.WILD_COL, constants.MUTANT_COL)
        df.columns = get_switched_column_names(df, constants.WILD_AA, constants.MUT_AA)

    experiment_policy.compute(df, outdir)

    was_found = [None for _ in range(df.shape[0])]

    # use grouping to avoid duplicate computation.
    for name, df_group in df.groupby(experiment_policy.get_groupby_cols(constants.WILD_COL,
                                                                        constants.MUTANT_COL)):

        # get PDB ids
        query_pdbid = df_group[constants.WILD_COL].iloc[0]
        target_pdbid = df_group[constants.MUTANT_COL].iloc[0]

        results_csv_path = experiment_policy.make_resultstatistic_path(query_pdbid, target_pdbid, outdir)
        if not results_csv_path.is_file():
            raise FileNotFoundError('resultStatistic.csv was not found!')

        df_resultStatistic = pd.read_csv(results_csv_path, sep='\t', header=0)

        for row in df_group.itertuples():
            row_pos_in_table = df.index.get_loc(row.Index)
            assert (type(row_pos_in_table) == int)  # there is something wrong if this is not int
            was_found[row_pos_in_table] = check_row_in_result(
                df_resultStatistic,
                row,
                backward)

    return was_found


# def eval_pair_tool2(df: pd.DataFrame, outdir: Path, cpus: int = 1, backward: bool = False, pdb_mirror='standard') -> list:
#     """
#
#     """
#
#     class PairPolicy:
#         @staticmethod
#         def compute(df, this_outdir):
#             # pdb_mirror is query only. The query might come from a specifically prepared mirror
#             run_pairs(df, this_outdir, cpus=cpus, pdb_mirror=pdb_mirror)
#
#         @staticmethod
#         def get_groupby_cols(query_col, target_col):
#             return [query_col, target_col]
#
#         @staticmethod
#         def make_resultstatistic_path(query_pdbid, target_pdbid, this_outdir):
#             return this_outdir / '{}_{}'.format(query_pdbid, target_pdbid) / 'resultStatistic.csv'
#
#     return run_evaluation2(df, PairPolicy, outdir, backward)


# def eval_screening_tool2(df: pd.DataFrame, outdir: Path, cpus: int = 1, backward: bool = False,
#                          pdb_mirror='standard') -> list:
#     """
#
#     """
#     class ScreenPolicy:
#         @staticmethod
#         def compute(df, this_outdir):
#             run_screens(df, this_outdir, cpus=cpus, pdb_mirror=pdb_mirror)
#
#         @staticmethod
#         def get_groupby_cols(query_col, _):
#             return [query_col]
#
#         @staticmethod
#         def make_resultstatistic_path(query_pdbid, target_pdbid, this_outdir):
#             return this_outdir / '{}'.format(query_pdbid) / 'resultStatistic.csv'
#
#     return run_evaluation2(df, ScreenPolicy, outdir, backward)
#
#
# def eval_pair_tool(df: pd.DataFrame, outdir: Path, backward: bool = False, pdb_mirror='standard') -> list:
#     """
#
#     """
#     class PairPolicy:
#         @staticmethod
#         def compute(query_pdbid, target_pdbid, this_outdir):
#             call_mutscreen_pairs(query_pdbid, target_pdbid, this_outdir,
#                                  pdb_mirror1=pdb_mirror,  # the query might come from a specifically prepared mirror
#                                  pdb_mirror2='standard')
#
#         @staticmethod
#         def get_groupby_cols(query_col, target_col):
#             return [query_col, target_col]
#
#     return run_evaluation(df, PairPolicy, outdir, backward)
#
#
# def eval_screening_tool(df: pd.DataFrame, outdir: Path, backward: bool = False, pdb_mirror='standard'):
#     """
#
#     """
#     class ScreenPolicy:
#         @staticmethod
#         def compute(query_pdbid, _, this_outdir):
#             call_mutscreen_screen(query_pdbid, this_outdir, pdb_mirror=pdb_mirror)
#
#         @staticmethod
#         def get_groupby_cols(query_col, _):
#             return [query_col]
#
#     return run_evaluation(df, ScreenPolicy, outdir, backward)
