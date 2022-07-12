"""Reads the mutation dataset from
Shanthirabalan S, Chomilier J, Carpentier M. Structural effects of point mutations in proteins.
Proteins. 2018 Aug;86(8):853-867. doi: 10.1002/prot.25499.
"""
import logging
from pathlib import Path

import pandas as pd

import helper

logger = logging.getLogger(__name__)


def read_shanthirabalan_single(path: Path, pdb_mutant_only: bool = True):
    log_prefix = 'Shanthirabalan: '

    df = pd.read_csv(path, sep='\t', dtype={helper.WILD_SEQ_NUM: str, helper.MUT_SEQ_NUM: str})

    logger.info(log_prefix + '{} data points / mutations in plain data set'.format(df.shape[0]))

    # standardize names
    df[helper.WILD_COL] = df[helper.WILD_COL].str.upper()
    df[helper.MUTANT_COL] = df[helper.MUTANT_COL].str.upper()
    df[helper.WILD_SEQ_NUM] = df[helper.WILD_SEQ_NUM].str.strip()
    df[helper.MUT_SEQ_NUM] = df[helper.MUT_SEQ_NUM].str.strip()

    b4 = df.shape[0]
    # drop wrong pairs that are not related in similarity measures
    df.query('tmalign_tmscore1 > 0.5 and tmalign_tmscore2 > 0.5 and tmalign_seqid > 0.8',
             inplace=True)

    logger.info(log_prefix + '{} of {} protein pairs meet the similarity criteria .'.format(
        df.shape[0], b4))

    b4 = df.shape[0]
    # drop pairs where no reasonable mutation position could be retrieved
    # this is the case if both proteins are actually identical, there are
    # multiple mutations or indels.
    df.query(
        f'`{helper.WILD_SEQ_NUM}` == `{helper.WILD_SEQ_NUM}` and `{helper.WILD_AA}` != "-" and '
        f'`{helper.MUT_SEQ_NUM}` == `{helper.MUT_SEQ_NUM}` and `{helper.MUT_AA}` != "-" ',
        inplace=True)

    logger.info(log_prefix + '{} of {} have a reasonable mutation position.'.format(
        df.shape[0], b4))

    return df
