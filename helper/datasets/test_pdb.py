import unittest
import tempfile
from pathlib import Path

import pandas as pd

from helper import CONFIG
from helper.datasets import PDB


class PDBTests(unittest.TestCase):
    """Test handling PDB data set"""

    def test_read(self):
        """Test to read single mutations from PDB"""

        with tempfile.TemporaryDirectory() as pdb_dir:
            pdb_dir_path = Path(pdb_dir)
            for pdbid in ["1g9v", "1xg1", "2rn2"]:
                filename = (
                    f"{CONFIG['PDB_FILE_INFO']['PREFIX']}{pdbid}"
                    f"{CONFIG['PDB_FILE_INFO']['SUFFIX']}"
                )
                with open(pdb_dir_path / filename, "w") as f:
                    pass

            exp_df = pd.DataFrame([["1G9V"], ["1XG1"], ["2RN2"]], columns=["id"])

            dataset = PDB(Path(pdb_dir))
            df = dataset.read()
            self.assertTrue(df[["id"]].equals(exp_df))
