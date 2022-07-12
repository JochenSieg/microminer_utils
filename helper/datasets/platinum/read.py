import pandas as pd
from pathlib import Path
import json
import helper
import logging

logger = logging.getLogger(__name__)

"""
Many columns seem very similar but I could not find docu describing their differences. For example
its not clear what 'MUT.WT_PDB' describes. Sometimes MUT.MT_PDB equals MUT.WT_PDB.

TODO: das Paper mal ueberfliegen. Vielleicht ist es wichtig bei welchr PDB id die affinity gemessen wurde.

The chain column is well defined, no missing values.
"""
PLATINUM_PDB_WILD = 'AFFIN.PDB_ID'  # die aber auch teilweise: 'MUT.WT_PDB'
PLATINUM_PDB_MUTANT = 'MUT.MT_PDB'
PLATINUM_MUTATION = 'MUTATION'
PLATINUM_WILD_CHAIN = 'AFFIN.CHAIN'

SINGLE_MUTATION_REGEX = r'^[A-Z]\d+[A-Z]$'
MUTATION_REGEX = r'^(?:[A-Z]\d+[A-Z]/?)+$'
PDB_REGEX = r'^[0-9][A-Z0-9]{3}$'


def read_platinum_cleaned_pdb(pdbid: str) -> Path:
    pdb_dir = Path(helper.CONFIG['DATA']['PLATINUM_PDBS'])
    filename = '{}.pdb'.format(pdbid.upper())
    path = pdb_dir / filename

    if not path.is_file():
        raise FileNotFoundError('Could not find PLATINUM cleaned PDB file: {}'.format(path))
    return path


def read_platinum_single(path: Path, pdb_mutant_only: bool) -> pd.DataFrame:
    """
    Read all valid single mutations in the PLATINUM data set to table.
    NOTE: Mutations are NOT necessarily unique.
    :param path: Path to data set file.
    :param pdb_mutant_only: Only return mutations for which both wild-type and
                            mutant have a PDB structure annotated.
    :return: Table of valid single mutations.
    """
    log_prefix = 'Platinum: '

    df = pd.read_csv(path, sep=',')

    logger.info(log_prefix + '{} data points / mutations in plain data set'.format(df.shape[0]))

    b4 = df.shape[0]
    interesting_fields = [PLATINUM_PDB_WILD, PLATINUM_MUTATION]
    df_pdb_mutations = df[df[interesting_fields].notnull().all(1)]
    logger.info(log_prefix + '{} of {} mutations are non-NaN for both PDB wild ID and mutation specification.'.format(
        df_pdb_mutations.shape[0], b4))

    # validate wild type identifiers
    b4 = df_pdb_mutations.shape[0]
    df_pdb_mutations = df_pdb_mutations[df_pdb_mutations[PLATINUM_PDB_WILD].str.contains(PDB_REGEX)]
    logger.info(log_prefix + '{} of {} mutations have valid WILD PDB ID.'.format(df_pdb_mutations.shape[0], b4))

    # validate mutation identifiers
    b4 = df_pdb_mutations.shape[0]
    df_pdb_mutations = df_pdb_mutations[df_pdb_mutations[PLATINUM_MUTATION].str.contains(MUTATION_REGEX)]
    logger.info(log_prefix + '{} of {} mutations have valid MUTATION identifier.'.format(df_pdb_mutations.shape[0], b4))

    if pdb_mutant_only:
        b4 = df_pdb_mutations.shape[0]
        df_pdb_mutations = df_pdb_mutations[df_pdb_mutations[PLATINUM_PDB_MUTANT].str.contains(PDB_REGEX)]
        df_pdb_mutations = df_pdb_mutations.query('`{}` != `{}`'.format(PLATINUM_PDB_MUTANT, PLATINUM_PDB_WILD))
        logger.info(log_prefix + '{} of {} mutations have valid MUTANT PDB ID.'.format(df_pdb_mutations.shape[0], b4))
        df_pdb_mutations[helper.MUTANT_COL] = df_pdb_mutations[PLATINUM_PDB_MUTANT].str.upper()

    df_single = df_pdb_mutations[df_pdb_mutations[PLATINUM_MUTATION].str.upper().str.contains(SINGLE_MUTATION_REGEX)]
    logger.info(log_prefix + '{} of {} mutations are single mutations'.format(df_single.shape[0], df_pdb_mutations.shape[0]))

    # if needed you can get multiple mutations with:
    # df_multiple = df[~df[PLATINUM_MUTATION].str.upper().str.contains(SINGLE_MUTATION_REGEX)]

    pd.options.mode.chained_assignment = None
    df_single[helper.WILD_COL] = df_single[PLATINUM_PDB_WILD].str[:4].str.upper()
    df_single[helper.WILD_AA] = df_single[PLATINUM_MUTATION].str[0:1].str.upper()
    df_single[helper.WILD_SEQ_NUM] = df_single[PLATINUM_MUTATION].str[1:-1].astype(str)
    df_single[helper.MUT_AA] = df_single[PLATINUM_MUTATION].str[-1:].str.upper()
    pd.options.mode.chained_assignment = 'warn'

    # ensure the index's integrity
    df_single.reset_index(drop=True, inplace=True)

    logger.info(log_prefix + 'Read {} mutations, {} WILD PDB ids, {} MUTANT PDB ids (single)'.format(
        df_single.shape[0], df_single[helper.WILD_COL].unique().__len__(),
        df_single[helper.MUTANT_COL].unique().__len__() if pdb_mutant_only else 0
    ))

    return df_single