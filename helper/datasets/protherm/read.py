import pandas as pd
from pathlib import Path
import helper
import logging

logger = logging.getLogger(__name__)


PROTHERM_PDB_WILD = 'pdb_wild'
PROTHERM_PDB_MUTANT = 'pdb_mutant'
PROTHERM_MUTATION = 'mutation'

SINGLE_MUTATION_REGEX = r'^[A-Z]\s\d+\s[A-Z]$'
MUTATION_REGEX = r'^(?:[A-Z]\s\d+\s[A-Z],?\s?)+$'
PDB_REGEX = r'^[0-9][A-Z0-9]{3}$'


def read_raw_protherm_csv(csv_path: Path) -> pd.DataFrame:
    """
    Simply reads the raw ProTherm CSV into a pandas DataFrame.
    Automatically applying using correct separator and encoding.
    :param csv_path: Path to CSV.
    :return: DataFrame of raw ProtTherm csv.
    """
    return pd.read_csv(csv_path, sep=',', header=0, encoding='latin-1')


def read_protherm_single(path: Path, pdb_mutant_only: bool):
    """
    Read all valid single mutations in the ProTherm data set to table.
    NOTE: Mutations are NOT necessarily unique.
    :param path: Path to data set file.
    :param pdb_mutant_only: Only return mutations for which both wild-type and
                            mutant have a PDB structure annotated.
    :return: Table of valid single mutations.
    """

    log_prefix = 'ProTherm: '

    df = read_raw_protherm_csv(path)

    logger.info(log_prefix + '{} data points / mutations in plain data set'.format(df.shape[0]))

    # problems with ProTherm format:
    # 1. pdb_mutant column has list entries like: 1OUA, 1LOZ
    # 2. mutation column has list entries like: C 77 A, C 95 A  (this means the proteins differ at least more than
    #    one mutation).
    # 3. mutation column has entries like: I 3 C (S-H) and another I 3 C (S-S). this is not a 'mutation' we are
    #    looking for.
    # 4. mutation column has entries like: Y 30 F (PDB: Y 32 F; PIR: Y 32 F).
    # 5. mutation column has entry with just the string "wild*". seems like a formatting error.
    # 6. Almost no mutation has a valid value for the field mutated_chain. Meaning there is no wild chain info.

    # keep only rows that have all the data we need.
    b4 = df.shape[0]
    interesting_fields = [PROTHERM_PDB_WILD, PROTHERM_MUTATION]
    df_pdb_mutations = df[df[interesting_fields].notnull().all(1)]
    logger.info(log_prefix + '{} of {} mutations are non-NaN for both PDB wild ID and mutation specification.'.format(
        df_pdb_mutations.shape[0], b4))

    # normalize relevant fields. Write them into new columns to avoid modifying the original data set.
    df_pdb_mutations[helper.WILD_COL] = df_pdb_mutations[PROTHERM_PDB_WILD].str.upper()
    # df_pdb_mutations['_this_normalized_mutation'] = df_pdb_mutations[PROTHERM_MUTATION].str.upper()

    # validate wild type identifiers
    b4 = df_pdb_mutations.shape[0]
    df_pdb_mutations = df_pdb_mutations[df_pdb_mutations[helper.WILD_COL].str.contains(PDB_REGEX)]
    logger.info(log_prefix + '{} of {} mutations have valid WILD PDB ID.'.format(df_pdb_mutations.shape[0], b4))

    # validate mutation identifiers
    b4 = df_pdb_mutations.shape[0]
    df_pdb_mutations = df_pdb_mutations[df_pdb_mutations[PROTHERM_MUTATION].str.upper().str.contains(MUTATION_REGEX)]
    logger.info(log_prefix + '{} of {} mutations have valid MUTATION identifier.'.format(df_pdb_mutations.shape[0], b4))

    if pdb_mutant_only:
        df_pdb_mutations = df_pdb_mutations.dropna(subset=[PROTHERM_PDB_MUTANT])
        df_pdb_mutations[helper.MUTANT_COL] = df_pdb_mutations[PROTHERM_PDB_MUTANT].str.upper()

        # 166H seems to have never existed in the PDB.
        strange_mutant_pdbids = ['166H']
        df_pdb_mutations = df_pdb_mutations[~df_pdb_mutations[helper.MUTANT_COL].isin(strange_mutant_pdbids)]

        # explode PDBid of mutant to different rows, if there are multiple PDBids in a single row.
        b4 = df_pdb_mutations.shape[0]
        df_pdb_mutations = df_pdb_mutations.assign(pdb_mutant=df_pdb_mutations[helper.MUTANT_COL].str.split(','))
        df_pdb_mutations = df_pdb_mutations.explode(helper.MUTANT_COL)
        logger.info(log_prefix + '{} mutations (from {}) after explosion MUTANT PDB ID lists.'.format(
            df_pdb_mutations.shape[0], b4))

        df_pdb_mutations[helper.MUTANT_COL] = df_pdb_mutations[helper.MUTANT_COL].str.strip()

        b4 = df_pdb_mutations.shape[0]
        df_pdb_mutations = df_pdb_mutations[df_pdb_mutations[helper.MUTANT_COL].str.contains(PDB_REGEX)]
        logger.info(log_prefix + '{} of {} mutations have valid MUTANT PDB ID.'.format(df_pdb_mutations.shape[0], b4))

    df_single = df_pdb_mutations[df_pdb_mutations[PROTHERM_MUTATION].str.upper().str.contains(SINGLE_MUTATION_REGEX)]
    logger.info(log_prefix + '{} of {} mutations are single mutations'.format(df_single.shape[0], df_pdb_mutations.shape[0]))

    pd.options.mode.chained_assignment = None
    df_single[helper.WILD_AA] = df_single[PROTHERM_MUTATION].str[0:1].str.upper()
    df_single[helper.WILD_SEQ_NUM] = df_single[PROTHERM_MUTATION].str[2:-2].astype(str)
    df_single[helper.MUT_AA] = df_single[PROTHERM_MUTATION].str[-1:].str.upper()
    pd.options.mode.chained_assignment = 'warn'

    # ensure the index's integrity
    df_single.reset_index(drop=True, inplace=True)

    logger.info(log_prefix + 'Read {} mutations, {} WILD PDB ids, {} MUTANT PDB ids (single)'.format(
        df_single.shape[0], df_single[helper.WILD_COL].unique().__len__(),
        df_single[helper.MUTANT_COL].unique().__len__() if pdb_mutant_only else 0
    ))

    return df_single
