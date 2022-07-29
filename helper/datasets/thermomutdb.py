import json
import logging
from pathlib import Path

import pandas as pd

import helper
from helper.datasets.dataset import MutationDataset

logger = logging.getLogger(__name__)


class ThermoMutDB(MutationDataset):
    """Wrapper for the ThermoMutDB data set.

       Publication https://doi.org/10.1093/nar/gkaa925

    Chain identifiers are little strange in this data set. Most are chain A and some are "unsigned".
    """

    name = "thermomutdb"
    has_custom_structure_files = False
    has_structure_pairs = True

    THERMOMUTDB_PDB_WILD = "PDB_wild"
    THERMOMUTDB_PDB_MUTANT = "pdb_mutant"
    THERMOMUTDB_MUTATION = "mutation_code"
    THERMOMUTDB_WILD_CHAIN = "mutated_chain"

    SINGLE_MUTATION_REGEX = r"^[A-Z]\d+[A-Z]$"
    MUTATION_REGEX = r"^(?:[A-Z]\d+[A-Z],?\s?)+$"
    PDB_REGEX = r"^[0-9][A-Z0-9]{3}$"
    SEQ_NUM_REGEX = r"^-?\d+[A-Za-z]?$"
    PDB_CHAIN_ID_REGEX = r"^[A-Za-z0-9]$"  # all valid PDB format chain IDs

    def __init__(self, file_path: Path):
        """Construct new instance.

        :param file_path: Path to flat file.
        """
        self.file_path = file_path

    def read(self) -> pd.DataFrame:
        """Read the raw data file to table.

        :return: Table of the raw data.
        """
        with open(self.file_path, "r") as f:
            return pd.DataFrame(json.load(f))

    def read_single_mutations(self, pdb_mutant_only: bool):
        """Read all valid single mutations in the ThermoMutDB data set to table.

        :param pdb_mutant_only: Whether to read only mutations with a structure for the mutant.
        :return: Table of valid single mutations.
        """
        df = self.read()
        logger.info(
            f"{ThermoMutDB.name}: "
            + "{} data points / mutations in plain data set".format(df.shape[0])
        )

        b4 = df.shape[0]
        interesting_fields = [
            ThermoMutDB.THERMOMUTDB_PDB_WILD,
            ThermoMutDB.THERMOMUTDB_MUTATION,
        ]
        df_pdb_mutations = df[
            df[interesting_fields].notnull().all(1)
        ].copy()  # copy to avoid verbose warnings
        logger.info(
            f"{ThermoMutDB.name}: "
            + "{} of {} mutations are non-NaN for both PDB wild ID and mutation specification.".format(
                df_pdb_mutations.shape[0], b4
            )
        )

        # validate wild type identifiers
        b4 = df_pdb_mutations.shape[0]
        df_pdb_mutations.loc[:, helper.WILD_COL] = (
            df_pdb_mutations[ThermoMutDB.THERMOMUTDB_PDB_WILD].str.upper().str.strip()
        )
        df_pdb_mutations = df_pdb_mutations[
            df_pdb_mutations[helper.WILD_COL].str.contains(ThermoMutDB.PDB_REGEX)
        ]
        logger.info(
            f"{ThermoMutDB.name}: "
            + "{} of {} mutations have valid WILD PDB ID.".format(
                df_pdb_mutations.shape[0], b4
            )
        )

        # validate mutation identifiers
        b4 = df_pdb_mutations.shape[0]
        df_pdb_mutations = df_pdb_mutations[
            df_pdb_mutations[ThermoMutDB.THERMOMUTDB_MUTATION]
            .str.upper()
            .str.contains(ThermoMutDB.MUTATION_REGEX)
        ]
        logger.info(
            f"{ThermoMutDB.name}: "
            + "{} of {} mutations have valid MUTATION identifier.".format(
                df_pdb_mutations.shape[0], b4
            )
        )

        if pdb_mutant_only:
            b4 = df_pdb_mutations.shape[0]
            df_pdb_mutations = df_pdb_mutations.dropna(
                subset=[ThermoMutDB.THERMOMUTDB_PDB_MUTANT]
            )

            df_pdb_mutations[helper.MUTANT_COL] = df_pdb_mutations[
                ThermoMutDB.THERMOMUTDB_PDB_MUTANT
            ].str.upper()

            # explode PDBid of mutant to different rows, if multiple PDBids in a single row.
            df_pdb_mutations.loc[:, helper.MUTANT_COL] = (
                df_pdb_mutations[helper.MUTANT_COL].str.strip().str.split(r"\s*,\s*")
            )
            df_pdb_mutations = df_pdb_mutations.explode(helper.MUTANT_COL)
            df_pdb_mutations.loc[:, helper.MUTANT_COL] = (
                df_pdb_mutations[helper.MUTANT_COL].str.strip().str.split("\s*/\s*")
            )
            df_pdb_mutations = df_pdb_mutations.explode(helper.MUTANT_COL)

            df_pdb_mutations[helper.MUTANT_COL] = df_pdb_mutations[
                helper.MUTANT_COL
            ].str.strip()

            df_pdb_mutations = df_pdb_mutations[
                df_pdb_mutations[helper.MUTANT_COL].str.contains(ThermoMutDB.PDB_REGEX)
            ]
            logger.info(
                f"{ThermoMutDB.name}: "
                + "{} of {} mutations have valid MUTANT PDB ID.".format(
                    df_pdb_mutations.shape[0], b4
                )
            )

        df_single = df_pdb_mutations[
            df_pdb_mutations[ThermoMutDB.THERMOMUTDB_MUTATION]
            .str.upper()
            .str.strip()
            .str.contains(ThermoMutDB.SINGLE_MUTATION_REGEX)
        ].copy()  # copy to avoid verbose warnings
        logger.info(
            f"{ThermoMutDB.name}: "
            + "{} of {} mutations are single mutations".format(
                df_single.shape[0], df_pdb_mutations.shape[0]
            )
        )

        df_single.loc[:, helper.WILD_AA] = df_single[
            ThermoMutDB.THERMOMUTDB_MUTATION
        ].str[0]
        df_single.loc[:, helper.WILD_SEQ_NUM] = (
            df_single[ThermoMutDB.THERMOMUTDB_MUTATION].str[1:-1].astype(str)
        )
        df_single.loc[:, helper.MUT_AA] = df_single[
            ThermoMutDB.THERMOMUTDB_MUTATION
        ].str[-1]
        # Comment out the following line(s) to not use the chain IDs of ThermoMutDB.
        # Currently, it makes no difference in the MM result (15 Nov 2022).
        df_single = df_single[
            df_single[ThermoMutDB.THERMOMUTDB_WILD_CHAIN]
            .str.strip()
            .str.contains(ThermoMutDB.PDB_CHAIN_ID_REGEX)
        ]
        df_single.loc[:, helper.WILD_CHAIN] = df_single[
            ThermoMutDB.THERMOMUTDB_WILD_CHAIN
        ].str.strip()

        # ensure the index's integrity
        df_single.reset_index(drop=True, inplace=True)

        logger.info(
            f"{ThermoMutDB.name}: "
            + "Read {} mutations, {} WILD PDB ids, {} MUTANT PDB ids (single)".format(
                df_single.shape[0],
                df_single[helper.WILD_COL].unique().__len__(),
                df_single[helper.MUTANT_COL].unique().__len__()
                if pdb_mutant_only
                else 0,
            )
        )

        return df_single
