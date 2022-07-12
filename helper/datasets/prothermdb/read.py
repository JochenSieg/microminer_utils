import pandas as pd
from pathlib import Path
import helper
import logging

logger = logging.getLogger(__name__)

PROTHERMDB_PDB_WILD_COL = 'PDB_wild'
PROTHERMDB_PDB_MUTATION = 'PDB_Chain_Mutation'

PDB_REGEX = r'^[0-9][A-Z0-9]{3}$'
MUTATION_REGEX = r'^(?:[0-9][A-Z0-9]{3}_[A-Z|0-9]:[A-Z]\d+[A-Z],?\s?)+$'  # example: 1IRO_A:V24I 1IRO_A:I33L
SINGLE_MUTATION_REGEX = r'^[0-9][A-Z0-9]{3}_[A-Z|0-9]:[A-Z]\d+[A-Z]$'
SEQ_NUM_REGEX = r'^\d+$'


# def split_pdb_mutation_string(df: pd.DataFrame) -> pd.DataFrame:
#     _original = pd.options.mode.chained_assignment
#     pd.options.mode.chained_assignment = None  # avoid SettingWithCopyWarning
#     _split = df[PROTHERMDB_PDB_MUTATION].str.split('_')
#     df.loc[:, helper.WILD_COL] = _split.str[0]
#     wild_chain = _split.str[1].str.split(':').str[0]
#     mutation_code = _split.str[1].str.split(':').str[1]
#     df.loc[:, helper.WILD_CHAIN] = wild_chain
#     df.loc[:, helper.WILD_AA] = mutation_code.str[0]
#     df.loc[:, helper.WILD_SEQ_NUM] = mutation_code.str[1:-1].astype(int)
#     df.loc[:, helper.MUT_AA] = mutation_code.str[-1]
#     pd.options.mode.chained_assignment = _original
#     return df
#
#
# def filter_prothermdb(df: pd.DataFrame):
#
#     b4 = df.shape[0]
#     interesting_fields = [PROTHERMDB_PDB_WILD_COL, PROTHERMDB_PDB_MUTATION]
#     df_pdb_mutations = df[df[interesting_fields].notnull().all(1)]
#     print(log_prefix + '{} of {} mutations are non-NaN for both PDB wild ID and mutation specification.'.format(
#         df_pdb_mutations.shape[0], b4))
#
#     # normalize wild PDB id and mutation string to upper case
#     df[helper.WILD_COL] = df[PROTHERMDB_PDB_WILD_COL].str.upper()
#     df[PROTHERMDB_PDB_MUTATION] = df[PROTHERMDB_PDB_MUTATION].str.upper()
#
#     # drop wild type protein PDB ids that are not matching PDB ID pattern
#     b4 = df.shape[0]
#     df = df[df[PROTHERMDB_PDB_WILD_COL].str.contains(PDB_REGEX)]
#     print('ProThermDB: After removing invalid wild-type PDB IDs {} / {} mutations remain'.format(df.shape[0], b4))
#
#     # drop mutation string not matching mutation string pattern
#     # Entries that not match are: 1YRI:A_L42A; 1YRI:A_L61A; 1YRI:A_L69A; 1YRI:A_L75A; -; 1AXB:A_G236S
#     b4 = df.shape[0]
#     df = df.dropna(subset=[PROTHERMDB_PDB_MUTATION])
#     df = df[df[PROTHERMDB_PDB_MUTATION].str.contains(MUTATION_REGEX)]
#     print('ProThermDB: After removing invalid mutation strings {} / {} mutations remain'.format(df.shape[0], b4))
#
#     b4 = df.shape[0]
#     df_single = df[df[PROTHERMDB_PDB_MUTATION].str.split(' ').str.len() == 1]
#     df_single = split_pdb_mutation_string(df_single)
#     print('ProThermDB: Extracted {} / {} SINGLE mutations'.format(df.shape[0], b4))
#
#     return df_single


def read_prothermdb_single(path: Path, pdb_mutant_only: bool) -> pd.DataFrame:
    """
    """
    # return filter_prothermdb(pd.read_csv(path, sep='\t'))

    log_prefix = 'ProThermDB: '

    df = pd.read_csv(path, sep='\t')

    logger.info(log_prefix + '{} data points / mutations in plain data set'.format(df.shape[0]))

    b4 = df.shape[0]
    interesting_fields = [PROTHERMDB_PDB_WILD_COL, PROTHERMDB_PDB_MUTATION]
    df_pdb_mutations = df[df[interesting_fields].notnull().all(1)]
    logger.info(log_prefix + '{} of {} mutations are non-NaN for both PDB wild ID and mutation specification.'.format(
        df_pdb_mutations.shape[0], b4))

    # validate wild type identifiers
    b4 = df_pdb_mutations.shape[0]
    df_pdb_mutations = df_pdb_mutations[df_pdb_mutations[PROTHERMDB_PDB_WILD_COL].str.upper().str.contains(PDB_REGEX)]
    logger.info(log_prefix + '{} of {} mutations have valid WILD PDB ID.'.format(df_pdb_mutations.shape[0], b4))

    # validate mutation identifiers
    b4 = df_pdb_mutations.shape[0]
    df_pdb_mutations = df_pdb_mutations[
        df_pdb_mutations[PROTHERMDB_PDB_MUTATION].str.upper().str.contains(MUTATION_REGEX)]
    logger.info(log_prefix + '{} of {} mutations have valid MUTATION identifier.'.format(df_pdb_mutations.shape[0], b4))

    df_single = df_pdb_mutations[
        df_pdb_mutations[PROTHERMDB_PDB_MUTATION].str.upper().str.strip().str.contains(SINGLE_MUTATION_REGEX)]
    logger.info(log_prefix + '{} of {} mutations are single mutations'.format(df_single.shape[0],
                                                                        df_pdb_mutations.shape[0]))

    _original = pd.options.mode.chained_assignment
    pd.options.mode.chained_assignment = None  # avoid SettingWithCopyWarning
    _split = df_single[PROTHERMDB_PDB_MUTATION].str.split('_')
    df_single[helper.WILD_COL] = _split.str[0].str.upper()
    wild_chain = _split.str[1].str.split(':').str[0]
    mutation_code = _split.str[1].str.split(':').str[1].str.upper()
    df_single[helper.WILD_CHAIN] = wild_chain
    df_single[helper.WILD_AA] = mutation_code.str[0].str.upper()
    df_single[helper.WILD_SEQ_NUM] = mutation_code.str[1:-1].astype(str)
    df_single[helper.MUT_AA] = mutation_code.str[-1].str.upper()
    pd.options.mode.chained_assignment = _original

    # ensure the index's integrity
    df_single.reset_index(drop=True, inplace=True)

    logger.info(log_prefix + 'Read {} mutations, {} WILD PDB ids, {} MUTANT PDB ids (single)'.format(
        df_single.shape[0], df_single[helper.WILD_COL].unique().__len__(), 0
    ))

    return df_single