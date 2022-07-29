import logging
from pathlib import Path

import pandas as pd

import helper
from helper.datasets.dataset import MutationDataset

logger = logging.getLogger(__name__)


class ProTherm(MutationDataset):
    """Wrapper for the ProTherm data set.

    Publication https://doi.org/10.1093/nar/gkj103

    Problems to handle with ProTherm format:
        1. pdb_mutant column has list entries like: 1OUA, 1LOZ
        2. mutation column has list entries like: C 77 A, C 95 A  (this means the proteins differ
           at least more than one mutation).
        3. mutation column has entries like: I 3 C (S-H) and another I 3 C (S-S). This is not a
           'mutation' we are looking for.
        4. mutation column has entries like: Y 30 F (PDB: Y 32 F; PIR: Y 32 F).
        5. mutation column has entry with just the string "wild*".
        6. Almost no mutation has a valid value for the field mutated_chain. Meaning there is no
           wild chain info.

    """

    name = "protherm"
    has_custom_structure_files = False
    has_structure_pairs = True

    PROTHERM_PDB_WILD = "pdb_wild"
    PROTHERM_PDB_MUTANT = "pdb_mutant"
    PROTHERM_MUTATION = "mutation"

    SINGLE_MUTATION_REGEX = r"^[A-Z]\s\d+\s[A-Z]$"
    MUTATION_REGEX = r"^(?:[A-Z]\s\d+\s[A-Z],?\s?)+$"
    PDB_REGEX = r"^[0-9][A-Z0-9]{3}$"

    def __init__(self, file_path: Path):
        """Construct new instance.

        :param file_path: Path to flat file.
        """
        self.file_path = file_path

    def read(self) -> pd.DataFrame:
        """Read the raw data set.

        :return: Raw data set as table.
        """
        return pd.read_csv(
            self.file_path, sep=",", header=0, encoding="latin-1", low_memory=False
        )

    def read_single_mutations(self, pdb_mutant_only: bool):
        """Read all valid single mutations in the ProTherm data set to table.

        :param pdb_mutant_only: Only return mutations for which both wild-type and
                                mutant have a PDB structure annotated.
        :return: Table of valid single mutations.
        """
        df = self.read()

        logger.info(
            ProTherm.name
            + ": "
            + "{} data points / mutations in plain data set".format(df.shape[0])
        )

        # keep only rows that have all the data we need.
        b4 = df.shape[0]
        interesting_fields = [ProTherm.PROTHERM_PDB_WILD, ProTherm.PROTHERM_MUTATION]
        df_pdb_mutations = df[
            df[interesting_fields].notnull().all(1)
        ].copy()  # copy to avoid verbose warnings
        logger.info(
            f"{ProTherm.name}: "
            + "{} of {} mutations are non-NaN for both PDB wild ID and mutation specification.".format(
                df_pdb_mutations.shape[0], b4
            )
        )

        # normalize relevant fields.
        df_pdb_mutations.loc[:, helper.WILD_COL] = df_pdb_mutations[
            ProTherm.PROTHERM_PDB_WILD
        ].str.upper()

        # validate wild type identifiers
        b4 = df_pdb_mutations.shape[0]
        df_pdb_mutations = df_pdb_mutations[
            df_pdb_mutations[helper.WILD_COL].str.contains(ProTherm.PDB_REGEX)
        ]
        logger.info(
            f"{ProTherm.name}: "
            + "{} of {} mutations have valid WILD PDB ID.".format(
                df_pdb_mutations.shape[0], b4
            )
        )

        # validate mutation identifiers
        b4 = df_pdb_mutations.shape[0]
        df_pdb_mutations = df_pdb_mutations[
            df_pdb_mutations[ProTherm.PROTHERM_MUTATION]
            .str.upper()
            .str.contains(ProTherm.MUTATION_REGEX)
        ]
        logger.info(
            f"{ProTherm.name}: "
            + "{} of {} mutations have valid MUTATION identifier.".format(
                df_pdb_mutations.shape[0], b4
            )
        )

        if pdb_mutant_only:
            df_pdb_mutations = df_pdb_mutations.dropna(
                subset=[ProTherm.PROTHERM_PDB_MUTANT]
            )
            df_pdb_mutations.loc[:, helper.MUTANT_COL] = df_pdb_mutations[
                ProTherm.PROTHERM_PDB_MUTANT
            ].str.upper()

            # 166H seems to have never existed in the PDB.
            strange_mutant_pdbids = ["166H"]
            df_pdb_mutations = df_pdb_mutations[
                ~df_pdb_mutations[helper.MUTANT_COL].isin(strange_mutant_pdbids)
            ]

            # explode PDBid of mutant to different rows
            b4 = df_pdb_mutations.shape[0]
            df_pdb_mutations.loc[:, helper.MUTANT_COL] = df_pdb_mutations[
                helper.MUTANT_COL
            ].str.split("\s*,\s*")
            df_pdb_mutations = df_pdb_mutations.explode(helper.MUTANT_COL).copy()
            logger.info(
                f"{ProTherm.name}: "
                + "{} mutations (from {}) after explosion MUTANT PDB ID lists.".format(
                    df_pdb_mutations.shape[0], b4
                )
            )

            df_pdb_mutations.loc[:, helper.MUTANT_COL] = df_pdb_mutations[
                helper.MUTANT_COL
            ].str.strip()

            b4 = df_pdb_mutations.shape[0]
            df_pdb_mutations = df_pdb_mutations[
                df_pdb_mutations[helper.MUTANT_COL].str.contains(ProTherm.PDB_REGEX)
            ]
            logger.info(
                f"{ProTherm.name}: "
                + "{} of {} mutations have valid MUTANT PDB ID.".format(
                    df_pdb_mutations.shape[0], b4
                )
            )

        df_single = df_pdb_mutations[
            df_pdb_mutations[ProTherm.PROTHERM_MUTATION]
            .str.upper()
            .str.contains(ProTherm.SINGLE_MUTATION_REGEX)
        ].copy()  # copy to avoid verbose warnings
        logger.info(
            f"{ProTherm.name}: "
            + "{} of {} mutations are single mutations".format(
                df_single.shape[0], df_pdb_mutations.shape[0]
            )
        )

        df_single.loc[:, helper.WILD_AA] = (
            df_single[ProTherm.PROTHERM_MUTATION].str[0:1].str.upper()
        )
        df_single.loc[:, helper.WILD_SEQ_NUM] = (
            df_single[ProTherm.PROTHERM_MUTATION].str[2:-2].astype(str)
        )
        df_single.loc[:, helper.MUT_AA] = (
            df_single[ProTherm.PROTHERM_MUTATION].str[-1:].str.upper()
        )

        # ensure the index's integrity
        df_single = df_single.reset_index(drop=True)

        logger.info(
            f"{ProTherm.name}: "
            + "Read {} mutations, {} WILD PDB ids, {} MUTANT PDB ids (single)".format(
                df_single.shape[0],
                len(df_single[helper.WILD_COL].unique()),
                len(df_single[helper.MUTANT_COL].unique()) if pdb_mutant_only else 0,
            )
        )

        return df_single
