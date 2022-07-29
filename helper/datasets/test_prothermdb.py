import tempfile
import unittest
from pathlib import Path

import pandas as pd

from helper import WILD_COL, WILD_AA, WILD_CHAIN, WILD_SEQ_NUM, MUTANT_COL, MUT_AA
from helper.datasets.prothermdb import ProThermDB


class ProThermDBTests(unittest.TestCase):
    """Test handling ProThermDB data set"""

    def test_read_single_mutations(self):
        """Test to read single mutations from ProThermDB CSV"""
        with tempfile.NamedTemporaryFile(mode="wt") as prothermdb_csv:
            prothermdb_csv.write(
                "\n".join(
                    [
                        "NO\tPROTEIN UniProt_ID\tMUTATION\tSOURCE\tPDB_wild\tPDB_Chain_Mutation",
                        "1\t\t\tP99S (Based on UniProt and PDB)\t\t2KQ6\t2kq6_B:P99S",
                        "3\t\t\tC86A C77T K9A (Based on UniProt)\t\t1RBB\t1rbb_A:C86A 1rbb_A:C77T"
                        " 1rbb_A:K9A",
                        "4\t\t\twild-type\t\t2KQ6\t-",  # no mutation
                    ]
                )
            )
            prothermdb_csv.flush()

            dataset = ProThermDB(Path(prothermdb_csv.name))

            df = dataset.read_single_mutations(pdb_mutant_only=False)
            self.assertEqual(df.shape[0], 1)
            self.assertFalse(MUTANT_COL in df.columns)  # mutant PDB column omitted
            self.assertEqual(df[WILD_COL].iloc[0], "2KQ6")
            self.assertEqual(df[WILD_CHAIN].iloc[0], "B")
            self.assertEqual(df[WILD_AA].iloc[0], "P")
            self.assertEqual(df[WILD_SEQ_NUM].iloc[0], "99")
            self.assertEqual(df[MUT_AA].iloc[0], "S")

            df = dataset.read_single_mutations(pdb_mutant_only=True)
            self.assertTrue(df.equals(pd.DataFrame()))
