import unittest
import tempfile
from pathlib import Path

from helper import WILD_COL, WILD_AA, WILD_CHAIN, WILD_SEQ_NUM, MUTANT_COL, MUT_AA
from helper.datasets.skempi2 import SKEMPI2


class SKEMPI2Tests(unittest.TestCase):
    """Test handling SKEMPI2 data set"""

    def test_read_single_mutations(self):
        """Test to read single mutations from SKEMPI2 CSV"""
        with tempfile.NamedTemporaryFile(mode="wt") as skempi2_csv:
            skempi2_csv.write(
                "\n".join(
                    [
                        "#Pdb;Mutation(s)_PDB;Mutation(s)_cleaned;iMutation_Location(s);Hold_out_type;"
                        "Hold_out_proteins;Affinity_mut (M);Affinity_mut_parsed;Affinity_wt (M);"
                        "Affinity_wt_parsed;Reference;Protein 1;Protein 2;Temperature;"
                        "kon_mut (M^(-1)s^(-1));kon_mut_parsed;kon_wt (M^(-1)s^(-1));"
                        "kon_wt_parsed;koff_mut (s^(-1));koff_mut_parsed;koff_wt (s^(-1));"
                        "koff_wt_parsed;dH_mut (kcal mol^(-1));dH_wt (kcal mol^(-1));"
                        "dS_mut (cal mol^(-1) K^(-1));dS_wt (cal mol^(-1) K^(-1));Notes;Method;"
                        "SKEMPI version",
                        "1CSB_E_I;II46G;II40G;;;;;;;;;;;;;;;;;;;;;;;;;1",
                        "1CSC_E_I;TE41A,IE55Y,SI99A;TE41A,IE55Y,SI99A;;;;;;;;;;;;;;;;;;;;;;;;;;1",
                    ]
                )
            )
            skempi2_csv.flush()

            dataset = SKEMPI2(Path(skempi2_csv.name))

            for mut_only in [True, False]:
                df = dataset.read_single_mutations(pdb_mutant_only=mut_only)
                self.assertEqual(df.shape[0], 1)
                self.assertFalse(MUTANT_COL in df.columns)  # mutant PDB column omitted
                self.assertEqual(df[WILD_COL].iloc[0], "1CSB")
                self.assertEqual(df[WILD_CHAIN].iloc[0], "I")
                self.assertEqual(df[WILD_AA].iloc[0], "I")
                self.assertEqual(df[WILD_SEQ_NUM].iloc[0], "40")
                self.assertEqual(df[MUT_AA].iloc[0], "G")
