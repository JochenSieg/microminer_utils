from pathlib import Path

import logging
import pandas as pd

from helper import CONFIG
from helper.datasets.dataset import Dataset
from helper.utils import scantree

logger = logging.getLogger(__name__)


class PDB(Dataset):
    """Wrapper to read PDB database structures"""

    name = "pdb"
    has_structure_pairs = False
    has_custom_structure_files = False

    def __init__(self, pdb_dir: Path):
        """Generate a new PDB wrapper.

        :param pdb_dir: The directory with structures.
        """
        self.pdb_dir = pdb_dir

    def read(self) -> pd.DataFrame:
        """Collect all structures of the PDB to a table.

        :return: A dataframe containing the PDB structure for querying.
        """
        # get all PDB IDs. There could be a better way, but we scan the PDB mirror on
        # the filesystem and parse PDBids from file names.
        file_prefix = CONFIG["PDB_FILE_INFO"]["PREFIX"]
        file_suffix = CONFIG["PDB_FILE_INFO"]["SUFFIX"]
        file_paths = [path for path in scantree(self.pdb_dir)]
        pdbids, structure_paths = zip(
            *[
                (p.name[len(file_prefix) : -len(file_suffix)].upper(), p)
                for p in file_paths
                if p.name.startswith(file_prefix) and p.name.endswith(file_suffix)
            ]
        )
        logger.info("Finished collecting PDB IDs from local PDB mirror")
        return pd.DataFrame({"id": pdbids, "structure_path": structure_paths})
