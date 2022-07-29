import logging
from pathlib import Path

import pandas as pd

from helper.datasets.dataset import Dataset
from helper.utils import scantree

logger = logging.getLogger(__name__)


class AFDB(Dataset):
    """Wrapper to read AlphaFold database structures"""

    name = "afdb"
    has_structure_pairs = False
    has_custom_structure_files = False

    def __init__(self, afdb_dir: Path):
        """Generate a new AFDB wrapper.

        :param afdb_dir: The directory with structures.
        """
        self.afdb_dir = afdb_dir
        self.file_extension = "pdb.gz"

    def read(self) -> pd.DataFrame:
        """Collect all structures of AFDB to a table (might get large).

        :return: A dataframe containing the AFDB structure for querying
        """
        # get all AFDB paths .
        structure_paths = [
            path
            for path in scantree(self.afdb_dir)
            if path.name.endswith(self.file_extension)
        ]

        # AF PDB files are name like this:
        #   AF-<UniprotID>-<fragmentID>-model_v<version>.pdb.gz
        #   we use AF-<UniprotID>-<fragmentID>-model_v<version> as ID
        ids = [p.stem.split(".")[0] for p in structure_paths]
        logger.info("Finished collecting AFDB IDs from local AFDB mirror")
        return pd.DataFrame({"id": ids, "structure_path": structure_paths})
