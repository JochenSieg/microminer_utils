import unittest
import tempfile
from pathlib import Path

import pandas as pd

from helper.datasets import AFDB


class AFDBTests(unittest.TestCase):
    """Test handling AlphaFoldDB data set"""

    def test_read(self):
        """Test reading file from alphafold database"""

        test_files = [
            "AF-A0A075B6H5-F1-model_v1.pdb.gz",
            "AF-A0A0G2KTJ5-F1-model_v1.pdb.gz",
            "AF-Q9XJ29-F1-model_v1.pdb.gz",
        ]
        exp_df = pd.DataFrame(
            [
                ["AF-A0A075B6H5-F1-model_v1"],
                ["AF-A0A0G2KTJ5-F1-model_v1"],
                ["AF-Q9XJ29-F1-model_v1"],
            ],
            columns=["id"],
        )
        exp_df["path"] = test_files

        with tempfile.TemporaryDirectory() as pdb_dir:
            afdb_dir_path = Path(pdb_dir)

            for filename in test_files:
                with open(afdb_dir_path / filename, "w") as f:
                    pass

            dataset = AFDB(Path(pdb_dir))
            df = dataset.read()
            self.assertTrue(df["id"].equals(exp_df["id"]))
