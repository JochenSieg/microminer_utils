import logging
from pathlib import Path

import pandas as pd

import helper
from helper.datasets.dataset import MutationDataset

logger = logging.getLogger(__name__)


class Shanthirabalan(MutationDataset):
    """Wrapper for the data set by Shanthirabalan et al..

    Publication https://doi.org/10.1002/prot.25499

    """

    name = "shanthirabalan"
    has_custom_structure_files = False
    has_structure_pairs = True

    def __init__(self, file_path: Path):
        """Construct new instance.

        :param file_path: Path to flat file.
        """
        self.file_path = file_path

    def read(self) -> pd.DataFrame:
        """Read the raw data file to table.

        :return: Table of the raw data.
        """
        return pd.read_csv(
            self.file_path,
            sep="\t",
            dtype={helper.WILD_SEQ_NUM: str, helper.MUT_SEQ_NUM: str},
        )

    def read_single_mutations(self, pdb_mutant_only: bool = True):
        """Read all valid single mutations in the Shanthirabalan et al. data set to table.

        :param pdb_mutant_only: Only return mutations for which both wild-type and
                                mutant have a PDB structure annotated.
        :return: Table of valid single mutations.
        """
        df = self.read()

        logger.info(
            f"{Shanthirabalan.name}: "
            + "{} data points / mutations in plain data set".format(df.shape[0])
        )

        # standardize names
        df[helper.WILD_COL] = df[helper.WILD_COL].str.upper()
        df[helper.MUTANT_COL] = df[helper.MUTANT_COL].str.upper()
        df[helper.WILD_SEQ_NUM] = df[helper.WILD_SEQ_NUM].str.strip()
        df[helper.MUT_SEQ_NUM] = df[helper.MUT_SEQ_NUM].str.strip()

        b4 = df.shape[0]
        # drop wrong pairs that are not related in similarity measures
        df.query(
            "tmalign_tmscore1 > 0.5"
            "and tmalign_tmscore2 > 0.5"
            "and tmalign_seqid > 0.8"
            " and has_non_terminal_indels == False",
            inplace=True,
        )

        logger.info(
            f"{Shanthirabalan.name}: "
            + "{} of {} protein pairs meet the similarity criteria .".format(
                df.shape[0], b4
            )
        )

        b4 = df.shape[0]
        # drop pairs where no reasonable mutation position could be retrieved
        # this is the case if both proteins are actually identical, there are
        # multiple mutations or indels.
        df.query(
            f'`{helper.WILD_SEQ_NUM}` == `{helper.WILD_SEQ_NUM}` and `{helper.WILD_AA}` != "-" and '
            f'`{helper.MUT_SEQ_NUM}` == `{helper.MUT_SEQ_NUM}` and `{helper.MUT_AA}` != "-" ',
            inplace=True,
        )

        logger.info(
            f"{Shanthirabalan.name}: "
            + "{} of {} have a reasonable mutation position.".format(df.shape[0], b4)
        )

        return df
