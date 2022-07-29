import logging
from pathlib import Path

import pandas as pd

from helper import WILD_COL, WILD_AA, WILD_CHAIN, WILD_SEQ_NUM, MUT_AA
from helper.datasets.dataset import MutationDataset

logger = logging.getLogger(__name__)


class FireProtDB(MutationDataset):
    """Wrapper for FireProtDB
    Publication https://doi.org/10.1093/nar/gkaa981
    """

    name = "fireprotdb"
    has_custom_structure_files = False
    has_structure_pairs = False

    PDB_WILD_COL = "pdb_id"
    CHAIN_COL = "chain"
    POSITION_COL = "position"
    WILD_AA_COL = "wild_type"
    MUT_AA_COL = "mutation"

    RELEVANT_COLS = [PDB_WILD_COL, CHAIN_COL, POSITION_COL, WILD_AA_COL, MUT_AA_COL]

    PDB_REGEX = r"^[0-9][A-Z0-9]{3}$"
    CHAIN_REGEX = r"^[A-Za-z0-9]$"
    AA_REGEX = r"^[ATGCDEFHIKLMNPQRSVWY]$"
    SEQ_NUM_REGEX = r"^-?\d+[A-Za-z]?$"

    def __init__(self, file_path: Path):
        """Create a new instance.

        :param file_path: Path to FireProtDB flat file.
        """
        self.file_path = file_path

    def read(self) -> pd.DataFrame:
        """Read the data sets raw data.

        :return: The flat file as table.
        """
        return pd.read_csv(
            self.file_path,
            sep=",",
            header=0,
            dtype={c: str for c in FireProtDB.RELEVANT_COLS},
            low_memory=False,
        )

    def read_single_mutations(self, pdb_mutant_only: bool) -> pd.DataFrame:
        """Read all single mutations.

        :param pdb_mutant_only: Read only single mutations that have
        :return:
        """
        if pdb_mutant_only:
            # does not contain mutant structures
            return pd.DataFrame()
        df = self.read()

        logger.info(
            f"{FireProtDB.name}: "
            + "{} data points / mutations in plain data set".format(df.shape[0])
        )
        b4 = df.shape[0]
        df = df[~df[FireProtDB.RELEVANT_COLS].isnull().any(axis=1)]
        if b4 > df.shape[0]:
            logger.info(
                f"{FireProtDB.name}: dropped {b4 - df.shape[0]} rows with nan in relevant "
                f"columns"
            )

        # there can be a list of PDB ids in a row like <PDBid1>|<PDBid2>
        df[FireProtDB.PDB_WILD_COL] = (
            df[FireProtDB.PDB_WILD_COL].str.split("|").tolist()
        )
        df = df.explode(FireProtDB.PDB_WILD_COL)
        df[FireProtDB.PDB_WILD_COL] = df[FireProtDB.PDB_WILD_COL].str.upper()
        df = df.drop_duplicates(FireProtDB.RELEVANT_COLS)
        if (
            df[FireProtDB.WILD_AA_COL].values == df[FireProtDB.MUT_AA_COL].values
        ).any():
            raise ValueError("Bad mutation data")

        df = df[df[FireProtDB.PDB_WILD_COL].str.contains(FireProtDB.PDB_REGEX)]
        df = df[df[FireProtDB.CHAIN_COL].str.contains(FireProtDB.CHAIN_REGEX)]
        df = df[df[FireProtDB.WILD_AA_COL].str.contains(FireProtDB.AA_REGEX)]
        df = df[df[FireProtDB.MUT_AA_COL].str.contains(FireProtDB.AA_REGEX)]
        df = df[df[FireProtDB.POSITION_COL].str.contains(FireProtDB.SEQ_NUM_REGEX)]

        df = df.rename(
            columns={
                FireProtDB.PDB_WILD_COL: WILD_COL,
                FireProtDB.CHAIN_COL: WILD_CHAIN,
                FireProtDB.WILD_AA_COL: WILD_AA,
                FireProtDB.MUT_AA_COL: MUT_AA,
                FireProtDB.POSITION_COL: WILD_SEQ_NUM,
            }
        )

        # ensure the index's integrity
        df.reset_index(drop=True, inplace=True)

        logger.info(
            f"{FireProtDB.name}: Read {df.shape[0]} mutations,"
            f" {len(df[WILD_COL].unique())} WILD PDB ids,"
            f" {0} MUTANT PDB ids (single)"
        )

        return df
