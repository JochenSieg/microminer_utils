import logging
from pathlib import Path

import pandas as pd

import helper
from helper.datasets.dataset import MutationDataset

logger = logging.getLogger(__name__)


class SKEMPI2(MutationDataset):
    """Wrapper for the data set by Shanthirabalan et al..

    Publication https://doi.org/10.1093/bioinformatics/bty635

    From the SKEMPI website FAQ (https://life.bsc.es/pid/skempi2/info/faq_and_help):
    "1) The PDB entry for the complex, followed by the chain identifiers for the two subunits.
    The first chain(s) correspond to protein 1 (column 10) and the second chain(s) correspond
    to protein 2 (column 11). Following this link will lead you to the relevant page in the
    protein databank.

    2) The mutation(s) corresponding to the residue numbering found in the protein database.
    The first character is the one letter amino acid code for the original residue, the
    second character is the chain identifier, the third to penultimate characters indicate
    the residue number, followed by the residue insertion code where applicable, and the
    final character indicates the mutant amino acid. Where multiple mutations are present,
    they are separated by commas."
    """

    name = "skempi2"
    has_custom_structure_files = True
    has_structure_pairs = False

    SKEMPI2_PDB_WILD_COL = "#Pdb"
    SKEMPI2_MUTATION = "Mutation(s)_cleaned"

    # PDB_REGEX = r'^[0-9][A-Z0-9]{3}$'
    SKEMPI_PDB_REGEX = r"^[0-9][A-Z0-9]{3}_[A-Z0-9]+_[A-Z0-9]+$"
    MUTATION_REGEX = r"^(?:[A-Z][A-Z0-9]\d+[A-Z],?)+$"
    SINGLE_MUTATION_REGEX = r"^[A-Z][A-Z0-9]\d+[A-Z]$"
    SEQ_NUM_REGEX = r"^-?\d+[A-Za-z]?$"

    def __init__(self, file_path: Path):
        """Construct new instance.

        :param file_path: Path to flat file.
        """
        self.file_path = file_path

    def read(self) -> pd.DataFrame:
        """Read the raw data file to table.

        :return: Table of the raw data.
        """
        return pd.read_csv(self.file_path, sep=";")

    def read_single_mutations(self, pdb_mutant_only: bool) -> pd.DataFrame:
        """Read all valid single mutations in the SKEMPI2.0 data set to table.

        :param pdb_mutant_only: Only for interface compatibility. SKEMPI does not have PDB mutants.
        :return: Table of valid single mutations.
        """
        df = self.read()
        logger.info(
            f"{SKEMPI2.name}: "
            + "{} data points / mutations in plain data set".format(df.shape[0])
        )

        b4 = df.shape[0]
        interesting_fields = [SKEMPI2.SKEMPI2_PDB_WILD_COL, SKEMPI2.SKEMPI2_MUTATION]
        df_pdb_mutations = df[df[interesting_fields].notnull().all(1)]
        logger.info(
            f"{SKEMPI2.name}: "
            + "{} of {} mutations are non-NaN for both PDB wild ID and mutation.".format(
                df_pdb_mutations.shape[0], b4
            )
        )

        # validate wild type identifiers
        b4 = df_pdb_mutations.shape[0]
        df_pdb_mutations = df_pdb_mutations[
            df_pdb_mutations[SKEMPI2.SKEMPI2_PDB_WILD_COL].str.contains(
                SKEMPI2.SKEMPI_PDB_REGEX
            )
        ]
        logger.info(
            f"{SKEMPI2.name}: "
            + "{} of {} mutations have valid WILD PDB ID.".format(
                df_pdb_mutations.shape[0], b4
            )
        )

        # validate mutation identifiers
        b4 = df_pdb_mutations.shape[0]
        df_pdb_mutations = df_pdb_mutations[
            df_pdb_mutations[SKEMPI2.SKEMPI2_MUTATION].str.contains(
                SKEMPI2.MUTATION_REGEX
            )
        ]
        logger.info(
            f"{SKEMPI2.name}: "
            + "{} of {} mutations have valid MUTATION identifier.".format(
                df_pdb_mutations.shape[0], b4
            )
        )

        df_single = df_pdb_mutations[
            df_pdb_mutations[SKEMPI2.SKEMPI2_MUTATION].str.contains(
                SKEMPI2.SINGLE_MUTATION_REGEX
            )
        ]
        logger.info(
            f"{SKEMPI2.name}: "
            + "{} of {} mutations are single mutations".format(
                df_single.shape[0], df_pdb_mutations.shape[0]
            )
        )

        pd.options.mode.chained_assignment = None
        df_single[helper.WILD_COL] = (
            df_single[SKEMPI2.SKEMPI2_PDB_WILD_COL].str[:4].str.upper()
        )
        df_single[helper.WILD_AA] = (
            df_single[SKEMPI2.SKEMPI2_MUTATION].str[0:1].str.upper()
        )
        df_single[helper.WILD_CHAIN] = df_single[SKEMPI2.SKEMPI2_MUTATION].str[1:2]
        df_single[helper.WILD_SEQ_NUM] = (
            df_single[SKEMPI2.SKEMPI2_MUTATION].str[2:-1].astype(str)
        )
        df_single[helper.MUT_AA] = (
            df_single[SKEMPI2.SKEMPI2_MUTATION].str[-1:].str.upper()
        )
        pd.options.mode.chained_assignment = "warn"

        # ensure the index's integrity
        df_single.reset_index(drop=True, inplace=True)

        logger.info(
            f"{SKEMPI2.name}: "
            + "Read {} mutations, {} WILD PDB ids, {} MUTANT PDB ids (single)".format(
                df_single.shape[0], df_single[helper.WILD_COL].unique().__len__(), 0
            )
        )

        return df_single

    @staticmethod
    def read_skempi_cleaned_pdb(pdbid: str) -> Path:
        """Read custom PDB files of the skempi data set.

        :param pdbid: PDB id for the structure.
        :return: Path to structure.
        """
        skempi_dir = Path(helper.CONFIG["DATA"]["SKEMPI2_PDBS"])
        filename = "{}.pdb".format(pdbid.upper())
        path = skempi_dir / filename

        if not path.is_file():
            raise FileNotFoundError(
                "Could not find SKEMPI2.0 cleaned PDB file: {}".format(path)
            )
        return path
