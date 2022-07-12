import pandas as pd
from pathlib import Path
import json
import helper
from helper import utils
import logging

logger = logging.getLogger(__name__)

THERMOMUTDB_PDB_WILD = 'PDB_wild'
THERMOMUTDB_PDB_MUTANT = 'pdb_mutant'
THERMOMUTDB_MUTATION = 'mutation_code'
THERMOMUTDB_WILD_CHAIN = 'mutated_chain'  # Like 98% are chain A and 7 are of chain "unsigned" which is problematic.

SINGLE_MUTATION_REGEX = r'^[A-Z]\d+[A-Z]$'
MUTATION_REGEX = r'^(?:[A-Z]\d+[A-Z],?\s?)+$'
PDB_REGEX = r'^[0-9][A-Z0-9]{3}$'
SEQ_NUM_REGEX = r'^\d+$'


def read_thermomutdb_json(json_path: Path):
    with open(json_path, 'r') as f:
        return json.load(f)


def read_thermomutdb_single(path: Path, pdb_mutant_only: bool):
    """

    :param path:
    :param pdb_mutant_only:
    :return:
    """

    log_prefix = 'ThermoMutDB: '

    df = pd.DataFrame(read_thermomutdb_json(path))

    logger.info(log_prefix + '{} data points / mutations in plain data set'.format(df.shape[0]))

    b4 = df.shape[0]
    interesting_fields = [THERMOMUTDB_PDB_WILD, THERMOMUTDB_MUTATION]
    df_pdb_mutations = df[df[interesting_fields].notnull().all(1)]
    logger.info(log_prefix + '{} of {} mutations are non-NaN for both PDB wild ID and mutation specification.'.format(
        df_pdb_mutations.shape[0], b4))

    # validate wild type identifiers
    b4 = df_pdb_mutations.shape[0]
    df_pdb_mutations[helper.WILD_COL] = df_pdb_mutations[THERMOMUTDB_PDB_WILD].str.upper().str.strip()
    df_pdb_mutations = df_pdb_mutations[df_pdb_mutations[helper.WILD_COL].str.contains(PDB_REGEX)]
    logger.info(log_prefix + '{} of {} mutations have valid WILD PDB ID.'.format(df_pdb_mutations.shape[0], b4))

    # validate mutation identifiers
    b4 = df_pdb_mutations.shape[0]
    df_pdb_mutations = df_pdb_mutations[df_pdb_mutations[THERMOMUTDB_MUTATION].str.upper().str.contains(MUTATION_REGEX)]
    logger.info(log_prefix + '{} of {} mutations have valid MUTATION identifier.'.format(df_pdb_mutations.shape[0], b4))

    if pdb_mutant_only:
        b4 = df_pdb_mutations.shape[0]
        df_pdb_mutations = df_pdb_mutations.dropna(subset=[THERMOMUTDB_PDB_MUTANT])

        df_pdb_mutations[helper.MUTANT_COL] = df_pdb_mutations[THERMOMUTDB_PDB_MUTANT].str.upper()

        # explode PDBid of mutant to different rows, if there are multiple PDBids in a single row.
        df_pdb_mutations = df_pdb_mutations.assign(pdb_mutant=df_pdb_mutations[helper.MUTANT_COL].str.split(','))
        df_pdb_mutations = df_pdb_mutations.explode(helper.MUTANT_COL)
        df_pdb_mutations = df_pdb_mutations.assign(pdb_mutant=df_pdb_mutations[helper.MUTANT_COL].str.split('/'))
        df_pdb_mutations = df_pdb_mutations.explode(helper.MUTANT_COL)

        df_pdb_mutations[helper.MUTANT_COL] = df_pdb_mutations[helper.MUTANT_COL].str.strip()

        df_pdb_mutations = df_pdb_mutations[df_pdb_mutations[helper.MUTANT_COL].str.contains(PDB_REGEX)]
        # df_pdb_mutations = df_pdb_mutations.query('`{}` != `{}`'.format(PLATINUM_PDB_MUTANT, PLATINUM_PDB_WILD))
        logger.info(log_prefix + '{} of {} mutations have valid MUTANT PDB ID.'.format(df_pdb_mutations.shape[0], b4))

    df_single = df_pdb_mutations[
        df_pdb_mutations[THERMOMUTDB_MUTATION].str.upper().str.strip().str.contains(SINGLE_MUTATION_REGEX)]
    logger.info(log_prefix + '{} of {} mutations are single mutations'.format(df_single.shape[0], df_pdb_mutations.shape[0]))

    # if needed you can get multiple mutations with:
    # df_multiple = df[~df[PLATINUM_MUTATION].str.upper().str.contains(SINGLE_MUTATION_REGEX)]

    pd.options.mode.chained_assignment = None
    df_single[helper.WILD_AA] = df_single[THERMOMUTDB_MUTATION].str[0]
    df_single[helper.WILD_SEQ_NUM] = df_single[THERMOMUTDB_MUTATION].str[1:-1].astype(str)
    df_single[helper.MUT_AA] = df_single[THERMOMUTDB_MUTATION].str[-1]
    df_single[helper.WILD_CHAIN] = df_single[THERMOMUTDB_WILD_CHAIN].str.strip()
    pd.options.mode.chained_assignment = 'warn'

    # ensure the index's integrity
    df_single.reset_index(drop=True, inplace=True)

    logger.info(log_prefix + 'Read {} mutations, {} WILD PDB ids, {} MUTANT PDB ids (single)'.format(
        df_single.shape[0], df_single[helper.WILD_COL].unique().__len__(),
        df_single[helper.MUTANT_COL].unique().__len__() if pdb_mutant_only else 0
    ))

    return df_single

