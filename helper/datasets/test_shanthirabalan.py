import unittest
import tempfile
from pathlib import Path

from helper import WILD_COL, WILD_AA, WILD_CHAIN, WILD_SEQ_NUM, MUTANT_COL, MUT_AA
from helper.datasets.shanthirabalan import Shanthirabalan


class ShanthirabalanTests(unittest.TestCase):
    """Test handling Shanthirabalan et al. data set"""

    def test_read_single_mutations(self):
        """Test to read single mutations from Shantirabalan et al. CSV"""
        with tempfile.NamedTemporaryFile(mode="wt") as shanthirabalan_csv:
            shanthirabalan_csv.write(
                "\n".join(
                    [
                        "wild_pdb\twild_chain\trms\tmut_pdb\tmut_chain\twild_aa\twild_seq_num\tmut_aa\t"
                        "mut_seq_num\ttmalign_rmsd\ttmalign_seqid\ttmalign_tmscore1\ttmalign_tmscore2\t"
                        "tmalign_seqlen1\ttmalign_seqlen2\ttmalign_alignlen\thas_non_terminal_indels",
                        "2d1c\tM\t0.24\t4dfa\tM\tT\t21\tL\t21\t0.3\t1.0\t0.9\t0.9\t10\t10\t10\tFalse",
                        # very low TMscores -> should be filtered out
                        "2d1d\tM\t0.24\t4dfb\tM\tT\t21\tL\t21\t0.3\t0.1\t0.1\t0.1\t10\t10\t10\tFalse",
                    ]
                )
            )
            shanthirabalan_csv.flush()

            dataset = Shanthirabalan(Path(shanthirabalan_csv.name))

            df = dataset.read_single_mutations(pdb_mutant_only=False)
            self.assertEqual(df.shape[0], 1)
            self.assertEqual(df[WILD_COL].iloc[0], "2D1C")
            # self.assertEqual(df[WILD_CHAIN].iloc[0], 'M')
            self.assertEqual(df[WILD_AA].iloc[0], "T")
            self.assertEqual(df[WILD_SEQ_NUM].iloc[0], "21")
            self.assertEqual(df[MUTANT_COL].iloc[0], "4DFA")
            self.assertEqual(df[MUT_AA].iloc[0], "L")
