import pandas as pd
from pathlib import Path
import json
import helper
import logging

logger = logging.getLogger(__name__)


def read_pisces() -> pd.DataFrame:
    log_prefix = 'PISCES: '
    df = pd.read_csv(helper.CONFIG['DATA']['PISCES'], sep='\s+', header=0)
    df['pdb'] = df['PDBchain'].str[:4]
    logger.info(log_prefix + f'read {df.shape[0]} chains, which are {len(df["pdb"].unique())}'
                             f' PDB IDs.')
    return df
