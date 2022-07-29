import unittest
import tempfile
from pathlib import Path

from helper.datasets.protherm import ProTherm
from helper import WILD_COL, WILD_AA, WILD_CHAIN, WILD_SEQ_NUM, MUTANT_COL, MUT_AA


class ProThermTests(unittest.TestCase):
    """Test handling ProTherm data set"""

    def test_read_single_mutations(self):
        """Test to read single mutations from ProTherm CSV"""
        with tempfile.NamedTemporaryFile(mode="wt") as protherm_csv:
            protherm_csv.write(
                "\n".join(
                    [
                        "pt_no,col1,protein,source,length,mol_weight,pir_id,swissprot_id,e_c_number,"
                        "pmd_no,pdb_wild,pdb_mutant,mutation,mutated_chain",
                        "1,,,,,,,,,,1BP2,,wild,",  # bad line
                        "2,,,,,,,,,,1BP2,,I 42 A,",  # single mutation
                        "2,,,,,,,,,,1G9V,2RN2,H 48 N,",  # mutation with PDB mutant structure
                        '2,,,,,,,,,,2LZM,169L,"E 128 A, V 131 A, N 132 A,",',  # multiple mutations
                    ]
                )
            )
            protherm_csv.flush()

            dataset = ProTherm(Path(protherm_csv.name))

            df = dataset.read_single_mutations(pdb_mutant_only=False)
            self.assertEqual(df.shape[0], 2)
            self.assertFalse(MUTANT_COL in df.columns)  # mutant PDB column omitted
            self.assertFalse(WILD_CHAIN in df.columns)
            self.assertEqual(df[WILD_COL].iloc[0], "1BP2")
            self.assertEqual(df[WILD_AA].iloc[0], "I")
            self.assertEqual(df[WILD_SEQ_NUM].iloc[0], "42")
            self.assertEqual(df[MUT_AA].iloc[0], "A")
            self.assertEqual(df[WILD_COL].iloc[1], "1G9V")
            self.assertEqual(df[WILD_AA].iloc[1], "H")
            self.assertEqual(df[WILD_SEQ_NUM].iloc[1], "48")
            self.assertEqual(df[MUT_AA].iloc[1], "N")

            df = dataset.read_single_mutations(pdb_mutant_only=True)
            self.assertEqual(df.shape[0], 1)
            self.assertFalse(WILD_CHAIN in df.columns)
            self.assertEqual(df[WILD_COL].iloc[0], "1G9V")
            self.assertEqual(df[MUTANT_COL].iloc[0], "2RN2")
            self.assertEqual(df[WILD_AA].iloc[0], "H")
            self.assertEqual(df[WILD_SEQ_NUM].iloc[0], "48")
            self.assertEqual(df[MUT_AA].iloc[0], "N")
