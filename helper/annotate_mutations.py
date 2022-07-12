"""Functionality to annotate mutations without PDB mutant with experimental "potential mutant" structure."""
from pathlib import Path
import pandas as pd
import multiprocessing
from .cmdl_calls import call_mutscreen_screen
from .constants import one_2_three_dict, WILD_COL, WILD_SEQ_NUM, WILD_AA, MUT_AA, WILD_CHAIN
from .execution_handling import run_screens


def rm_non_standard_aas(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove non standard amino acids as wild-type and mutant. For example X.
    :param df: table to filter
    :return: filtered table
    """
    df = df[df[WILD_AA].isin(one_2_three_dict) & df[MUT_AA].isin(one_2_three_dict)]
    return df


def annotate_mutations(df: pd.DataFrame, outdir: Path, use_chainid=False, pdb_mirror='standard') -> pd.DataFrame:
    """

    :param df:
    :param outdir:
    :param use_chainid: Whether to use chain identifier to identfy wild-type mutation position. Some data sets do not
                        provide this information and an error will be thrown.
    :param pdb_mirror: The pdb_mirror to use to retrieve wild-type protein structures. Options: standard, skempi2
    :return:
    """

    b4 = df.shape[0]
    df = rm_non_standard_aas(df)
    print('Removed {} non-standard Amino Acids in mutations'.format(b4 - df.shape[0]))

    print('Starting to annotate PDB mutants for {} wild-type proteins from {} DB mutations.'.format(
        len(df[WILD_COL].unique()), df.shape[0]))

    for query_pdbid, df_group in df.groupby([WILD_COL]):

        this_outdir = outdir / query_pdbid
        this_outdir.mkdir(parents=True, exist_ok=True)
 
        try:
            call_mutscreen_screen(query_pdbid, this_outdir, pdb_mirror)
        except ValueError as e:
            # dont abort if a single screen fails. Failed screens are documented in the
            # log with this print. You can grep for it:
            print(e)

    result_dfs = []
    for query_pdbid, df_group in df.groupby([WILD_COL]):

        this_outdir = outdir / query_pdbid
        this_outdir.mkdir(parents=True, exist_ok=True)

        results_csv_path = this_outdir / 'resultStatistic.csv'
        if not results_csv_path.is_file():
            print('Warning: results for {} missing at {}'.format(query_pdbid, results_csv_path.absolute()))
            continue
            # raise FileNotFoundError('resultStatistic.csv for {} was not generated! '.format(query_pdbid))
        result_dfs.append(pd.read_csv(results_csv_path, sep='\t', header=0))

    df_all_results = pd.concat(result_dfs)

    df['wildAA_3letter'] = df[WILD_AA].apply(lambda aa: one_2_three_dict[aa] if len(aa) == 1 else aa)
    df['mutantAA_3letter'] = df[MUT_AA].apply(lambda aa: one_2_three_dict[aa] if len(aa) == 1 else aa)

    left_cols = [WILD_COL, 'wildAA_3letter', WILD_SEQ_NUM, 'mutantAA_3letter']
    right_cols = ['wildName', 'wildAA', 'wildPos', 'mutantAA']
    if use_chainid:
        left_cols.append(WILD_CHAIN)
        right_cols.append('wildChain')

    df_annotated = df.merge(df_all_results, left_on=left_cols, right_on=right_cols, how='left')
    # df_annotated.drop('siteId', axis=1, inplace=True)
    # df_annotated.drop('mutationId', axis=1, inplace=True)
    df_annotated.drop('wildAA_3letter', axis=1, inplace=True)
    df_annotated.drop('mutantAA_3letter', axis=1, inplace=True)
    print()
    return df_annotated


def annotate_mutations2(df: pd.DataFrame, outdir: Path, use_chainid=False, cpus: int = 1,
                        pdb_mirror='standard') -> pd.DataFrame:
    """

    :param df:
    :param outdir:
    :param use_chainid: Whether to use chain identifier to identfy wild-type mutation position. Some data sets do not
                        provide this information and an error will be thrown.
    :param pdb_mirror: The pdb_mirror to use to retrieve wild-type protein structures. Options: standard, skempi2
    :return:
    """

    b4 = df.shape[0]
    df = rm_non_standard_aas(df)
    print('Removed {} non-standard Amino Acids in mutations'.format(b4 - df.shape[0]))

    print('Starting to annotate PDB mutants for {} wild-type proteins from {} DB mutations.'.format(
        len(df[WILD_COL].unique()), df.shape[0]))

    run_screens(df, outdir, cpus=cpus, pdb_mirror=pdb_mirror, raise_error=False)
    # for query_pdbid, df_group in df.groupby([WILD_COL]):
    #
    #     this_outdir = outdir / query_pdbid
    #     this_outdir.mkdir(parents=True, exist_ok=True)
    #
    #     try:
    #         call_mutscreen_screen(query_pdbid, this_outdir, pdb_mirror)
    #     except ValueError as e:
    #         # dont abort if a single screen fails. Failed screens are documented in the
    #         # log with this print. You can grep for it:
    #         print(e)

    result_dfs = []
    for query_pdbid, df_group in df.groupby([WILD_COL]):

        this_outdir = outdir / query_pdbid
        this_outdir.mkdir(parents=True, exist_ok=True)

        results_csv_path = this_outdir / 'resultStatistic.csv'
        if not results_csv_path.is_file():
            print('Warning: results for {} missing at {}'.format(query_pdbid, results_csv_path.absolute()))
            continue
            # raise FileNotFoundError('resultStatistic.csv for {} was not generated! '.format(query_pdbid))
        result_dfs.append(pd.read_csv(results_csv_path, sep='\t', header=0))

    df_all_results = pd.concat(result_dfs)

    df['wildAA_3letter'] = df[WILD_AA].apply(lambda aa: one_2_three_dict[aa] if len(aa) == 1 else aa)
    df['mutantAA_3letter'] = df[MUT_AA].apply(lambda aa: one_2_three_dict[aa] if len(aa) == 1 else aa)

    left_cols = [WILD_COL, 'wildAA_3letter', WILD_SEQ_NUM, 'mutantAA_3letter']
    right_cols = ['wildName', 'wildAA', 'wildPos', 'mutantAA']
    if use_chainid:
        left_cols.append(WILD_CHAIN)
        right_cols.append('wildChain')

    df_annotated = df.merge(df_all_results, left_on=left_cols, right_on=right_cols, how='left')
    # df_annotated.drop('siteId', axis=1, inplace=True)
    # df_annotated.drop('mutationId', axis=1, inplace=True)
    df_annotated.drop('wildAA_3letter', axis=1, inplace=True)
    df_annotated.drop('mutantAA_3letter', axis=1, inplace=True)
    print()
    return df_annotated
