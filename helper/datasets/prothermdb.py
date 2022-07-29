import logging
from pathlib import Path

import pandas as pd

import helper
from helper.datasets.dataset import MutationDataset

logger = logging.getLogger(__name__)


class ProThermDB(MutationDataset):
    """Wrapper for the ProTherm data set.

    Publication https://doi.org/10.1093/nar/gkaa1035

    """

    name = "prothermdb"
    has_custom_structure_files = False
    has_structure_pairs = False

    PDB_WILD_COL = "PDB_wild"
    PDB_MUTATION = "PDB_Chain_Mutation"

    PDB_REGEX = r"^[0-9][A-Z0-9]{3}$"
    MUTATION_REGEX = r"^(?:[0-9][A-Z0-9]{3}_[A-Z|0-9]:[A-Z]\d+[A-Z],?\s?)+$"  # example: 1IRO_A:V24I 1IRO_A:I33L
    SINGLE_MUTATION_REGEX = r"^[0-9][A-Z0-9]{3}_[A-Z|0-9]:[A-Z]\d+[A-Z]$"
    SEQ_NUM_REGEX = r"^-?\d+[A-Za-z]?$"

    def __init__(self, file_path: Path):
        """Construct a new instance.

        :param file_path: Path to flat file.
        """
        self.file_path = file_path

    def read(self) -> pd.DataFrame:
        """Read the raw data set.

        :return: Raw data set as table.
        """
        return pd.read_csv(self.file_path, sep="\t")

    def read_single_mutations(self, pdb_mutant_only: bool) -> pd.DataFrame:
        """Read all valid single mutations in the ProThermDB data set to table.

        :param pdb_mutant_only: Only return mutations for which both wild-type and
                                mutant have a PDB structure annotated.
        :return: Table of valid single mutations.
        """
        if pdb_mutant_only:
            # does not contain mutant structures
            return pd.DataFrame()
        df = self.read()
        logger.info(
            f"{ProThermDB.name}: "
            + "{} data points / mutations in plain data set".format(df.shape[0])
        )

        b4 = df.shape[0]
        interesting_fields = [ProThermDB.PDB_WILD_COL, ProThermDB.PDB_MUTATION]
        df_pdb_mutations = df[df[interesting_fields].notnull().all(1)]
        logger.info(
            f"{ProThermDB.name}: "
            + "{} of {} mutations are non-NaN for both PDB wild ID and mutation.".format(
                df_pdb_mutations.shape[0], b4
            )
        )

        # validate wild type identifiers
        b4 = df_pdb_mutations.shape[0]
        df_pdb_mutations = df_pdb_mutations[
            df_pdb_mutations[ProThermDB.PDB_WILD_COL]
            .str.upper()
            .str.contains(ProThermDB.PDB_REGEX)
        ]
        logger.info(
            f"{ProThermDB.name}: "
            + "{} of {} mutations have valid WILD PDB ID.".format(
                df_pdb_mutations.shape[0], b4
            )
        )

        # validate mutation identifiers
        b4 = df_pdb_mutations.shape[0]
        df_pdb_mutations = df_pdb_mutations[
            df_pdb_mutations[ProThermDB.PDB_MUTATION]
            .str.upper()
            .str.contains(ProThermDB.MUTATION_REGEX)
        ]
        logger.info(
            f"{ProThermDB.name}: "
            + "{} of {} mutations have valid MUTATION identifier.".format(
                df_pdb_mutations.shape[0], b4
            )
        )

        df_single = df_pdb_mutations[
            df_pdb_mutations[ProThermDB.PDB_MUTATION]
            .str.upper()
            .str.strip()
            .str.contains(ProThermDB.SINGLE_MUTATION_REGEX)
        ]
        logger.info(
            f"{ProThermDB.name}: "
            + "{} of {} mutations are single mutations".format(
                df_single.shape[0], df_pdb_mutations.shape[0]
            )
        )

        _original = pd.options.mode.chained_assignment
        pd.options.mode.chained_assignment = None  # avoid SettingWithCopyWarning
        _split = df_single[ProThermDB.PDB_MUTATION].str.split("_")
        df_single[helper.WILD_COL] = _split.str[0].str.upper()
        wild_chain = _split.str[1].str.split(":").str[0]
        mutation_code = _split.str[1].str.split(":").str[1].str.upper()
        df_single[helper.WILD_CHAIN] = wild_chain
        df_single[helper.WILD_AA] = mutation_code.str[0].str.upper()
        df_single[helper.WILD_SEQ_NUM] = mutation_code.str[1:-1].astype(str)
        df_single[helper.MUT_AA] = mutation_code.str[-1].str.upper()
        pd.options.mode.chained_assignment = _original

        # ensure the index's integrity
        df_single.reset_index(drop=True, inplace=True)

        logger.info(
            f"{ProThermDB.name}: "
            + "Read {} mutations, {} WILD PDB ids, {} MUTANT PDB ids (single)".format(
                df_single.shape[0], df_single[helper.WILD_COL].unique().__len__(), 0
            )
        )

        return df_single
