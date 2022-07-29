import logging
from pathlib import Path

import pandas as pd

import helper
from helper.datasets.dataset import MutationDataset

logger = logging.getLogger(__name__)


class Platinum(MutationDataset):
    """Wrapper for the Platinum data set.

    Publication https://doi.org/10.1093/nar/gku966

    Many columns in the Platinum flat file seem very similar, but I could not find docu
    describing their differences. For example, it's not clear what 'MUT.WT_PDB' describes.
    Sometimes MUT.MT_PDB equals MUT.WT_PDB.

    The chain column is well-curated, no missing values (in contrast to other data sets).
    """

    name = "platinum"
    has_custom_structure_files = True
    has_structure_pairs = True

    PLATINUM_PDB_WILD = (
        "AFFIN.PDB_ID"  # seems to give the wild PDB code, but 'MUT.WT_PDB' also
    )
    PLATINUM_PDB_MUTANT = "MUT.MT_PDB"
    PLATINUM_MUTATION = "MUTATION"
    PLATINUM_WILD_CHAIN = "AFFIN.CHAIN"

    SINGLE_MUTATION_REGEX = r"^[A-Z]\d+[A-Z]$"
    MUTATION_REGEX = r"^(?:[A-Z]\d+[A-Z]/?)+$"
    PDB_REGEX = r"^[0-9][A-Z0-9]{3}$"

    def __init__(self, file_path: Path):
        """Construct a new instance.

        :param file_path: Path to flat file.
        """
        self.file_path = file_path

    def read(self) -> pd.DataFrame:
        """Read the raw data set.

        :return: Table of the raw data set.
        """
        return pd.read_csv(self.file_path, sep=",")

    def read_single_mutations(self, pdb_mutant_only: bool) -> pd.DataFrame:
        """Read all valid single mutations in the PLATINUM data set to table.

        :param pdb_mutant_only: Only return mutations for which both wild-type and
                                mutant have a PDB structure annotated.
        :return: Table of valid single mutations.
        """
        df = self.read()

        logger.info(
            f"{Platinum.name}: "
            + "{} data points / mutations in plain data set".format(df.shape[0])
        )

        b4 = df.shape[0]
        interesting_fields = [Platinum.PLATINUM_PDB_WILD, Platinum.PLATINUM_MUTATION]
        df_pdb_mutations = df[
            df[interesting_fields].notnull().all(1)
        ].copy()  # copy to avoid verbose warnings
        logger.info(
            f"{Platinum.name}: {df_pdb_mutations.shape[0]} of {b4} mutations are non-NaN "
            f"for both PDB wild ID and mutation specification."
        )

        # validate wild type identifiers
        b4 = df_pdb_mutations.shape[0]
        df_pdb_mutations = df_pdb_mutations[
            df_pdb_mutations[Platinum.PLATINUM_PDB_WILD].str.contains(
                Platinum.PDB_REGEX
            )
        ]
        logger.info(
            f"{Platinum.name}: {df_pdb_mutations.shape[0]} of {b4} mutations have"
            f" valid WILD PDB ID."
        )

        # validate mutation identifiers
        b4 = df_pdb_mutations.shape[0]
        df_pdb_mutations = df_pdb_mutations[
            df_pdb_mutations[Platinum.PLATINUM_MUTATION].str.contains(
                Platinum.MUTATION_REGEX
            )
        ]
        logger.info(
            f"{Platinum.name}: {df_pdb_mutations.shape[0]} of {b4} mutations have valid MUTATION"
            f" identifier."
        )

        if pdb_mutant_only:
            b4 = df_pdb_mutations.shape[0]
            df_pdb_mutations = df_pdb_mutations[
                df_pdb_mutations[Platinum.PLATINUM_PDB_MUTANT].str.contains(
                    Platinum.PDB_REGEX
                )
            ]
            df_pdb_mutations = df_pdb_mutations.query(
                "`{}` != `{}`".format(
                    Platinum.PLATINUM_PDB_MUTANT, Platinum.PLATINUM_PDB_WILD
                )
            ).copy()
            logger.info(
                f"{Platinum.name}: "
                + "{} of {} mutations have valid MUTANT PDB ID.".format(
                    df_pdb_mutations.shape[0], b4
                )
            )
            df_pdb_mutations.loc[:, helper.MUTANT_COL] = df_pdb_mutations[
                Platinum.PLATINUM_PDB_MUTANT
            ].str.upper()

        df_single = df_pdb_mutations[
            df_pdb_mutations[Platinum.PLATINUM_MUTATION]
            .str.upper()
            .str.contains(Platinum.SINGLE_MUTATION_REGEX)
        ]
        logger.info(
            f"{Platinum.name}: "
            + "{} of {} mutations are single mutations".format(
                df_single.shape[0], df_pdb_mutations.shape[0]
            )
        )

        # if needed you can get multiple mutations with:
        # df_multiple = df[~df[PLATINUM_MUTATION].str.upper().str.contains(SINGLE_MUTATION_REGEX)]

        pd.options.mode.chained_assignment = None
        df_single[helper.WILD_COL] = (
            df_single[Platinum.PLATINUM_PDB_WILD].str[:4].str.upper()
        )
        df_single[helper.WILD_CHAIN] = df_single[
            Platinum.PLATINUM_WILD_CHAIN
        ].str.strip()
        df_single[helper.WILD_AA] = (
            df_single[Platinum.PLATINUM_MUTATION].str[0:1].str.upper()
        )
        df_single[helper.WILD_SEQ_NUM] = (
            df_single[Platinum.PLATINUM_MUTATION].str[1:-1].astype(str)
        )
        df_single[helper.MUT_AA] = (
            df_single[Platinum.PLATINUM_MUTATION].str[-1:].str.upper()
        )
        pd.options.mode.chained_assignment = "warn"

        # ensure the index's integrity
        df_single.reset_index(drop=True, inplace=True)

        logger.info(
            f"{Platinum.name}: "
            + "Read {} mutations, {} WILD PDB ids, {} MUTANT PDB ids (single)".format(
                df_single.shape[0],
                df_single[helper.WILD_COL].unique().__len__(),
                df_single[helper.MUTANT_COL].unique().__len__()
                if pdb_mutant_only
                else 0,
            )
        )

        return df_single

    @staticmethod
    def read_platinum_cleaned_pdb(pdbid: str) -> Path:
        """Read custom PDB files of the platinum data set.

        :param pdbid: PDB id for the structure.
        :return: Path to structure.
        """
        pdb_dir = Path(helper.CONFIG["DATA"]["PLATINUM_PDBS"])
        filename = "{}.pdb".format(pdbid.upper())
        path = pdb_dir / filename

        if not path.is_file():
            raise FileNotFoundError(
                "Could not find PLATINUM cleaned PDB file: {}".format(path)
            )
        return path
