import unittest
import tempfile
from pathlib import Path

from helper.datasets.platinum import Platinum
from helper import WILD_COL, WILD_AA, WILD_CHAIN, WILD_SEQ_NUM, MUTANT_COL, MUT_AA


class PlatinumTests(unittest.TestCase):
    """Test handling Platinum data set"""

    def test_read_single_mutations(self):
        """Test to read single mutations from Platinum CSV"""
        with tempfile.NamedTemporaryFile(mode="wt") as platinum_csv:
            platinum_csv.write(
                "\n".join(
                    [
                        "MUTATION,AFFIN.CHAIN,AFFIN.LIG_ID,AFFIN.LIG_NAME,AFFIN.PDB_ID,MUT.UNIPROT,"
                        "AFFIN.EXPTAL_METHOD,AFFIN.TEMPERATURE,AFFIN.PH,AFFIN.K_WT,AFFIN.K_MT,"
                        "AFFIN.DELTA_K,AFFIN.FOLD_CHANGE,AFFIN.UNIT,MUT.MT_PDB,MUT.WT_PDB,"
                        "MUT.IN_BINDING_SITE,MUT.DISTANCE_TO_LIG,MUT.DOI,MUT.PMID,"
                        "MUT.IS_SINGLE_POINT,MUT.MCSM_PPI,MUT.DUET_STABILITY,MUT.RESIDUE_RSA,"
                        "MUT.RESIDUE_SST,MUT.RESIDUE_HBOND,PROT.MOLECULE_NAME,PROT.ORGANISM,"
                        "PROT.RESOLUTION,PROT.R_VALUE,PROT.R_FREE,PROT.PH,PROT.TEMPERATURE,"
                        "PROT.EXPTAL_METHOD,PROT.STOICHIOMETRY,LIG.PDBE,LIG.CANONICAL_SMILES,"
                        "LIG.POLYMER_TYPE,LIG.TYPE_DESCRIPTION,LIG.MOL_WEIGHT,LIG.ATOM_COUNT,"
                        "LIG.LOGP,LIG.FORMULA,LIG.NUM_ACCEPTORS,LIG.NUM_DONORS,LIG.NUM_HETEROATOMS,"
                        "LIG.NUM_ROTATABLE_BONDS,LIG.NUM_RINGS,LIG.LABUTE_ASA,LIG.TPSA,"
                        "PROT.PROTEIN_CLASS",
                        ",".join(
                            [
                                "I89A",
                                "A",
                                "",
                                "",
                                "2IEM",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "2EIO",
                                "2IEM",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                            ]
                        ),
                        ",".join(
                            [
                                "I1V",
                                "A",
                                "",
                                "",
                                "1G9V",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "NO",
                                "1G9V",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                            ]
                        ),
                        ",".join(
                            [
                                "D31N/N87D/L90M",
                                "A",
                                "",
                                "",
                                "1G9V",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "NO",
                                "1G9V",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                            ]
                        ),
                    ]
                )
            )
            platinum_csv.flush()

            dataset = Platinum(Path(platinum_csv.name))

            df = dataset.read_single_mutations(pdb_mutant_only=False)
            self.assertEqual(df.shape[0], 2)
            self.assertFalse(MUTANT_COL in df.columns)  # mutant PDB column omitted
            self.assertEqual(df[WILD_CHAIN].iloc[0], "A")
            self.assertEqual(df[WILD_COL].iloc[0], "2IEM")
            self.assertEqual(df[WILD_AA].iloc[0], "I")
            self.assertEqual(df[WILD_SEQ_NUM].iloc[0], "89")
            self.assertEqual(df[MUT_AA].iloc[0], "A")
            self.assertEqual(df[WILD_COL].iloc[1], "1G9V")
            self.assertEqual(df[WILD_AA].iloc[1], "I")
            self.assertEqual(df[WILD_SEQ_NUM].iloc[1], "1")
            self.assertEqual(df[MUT_AA].iloc[1], "V")

            df = dataset.read_single_mutations(pdb_mutant_only=True)
            self.assertEqual(df.shape[0], 1)
            self.assertEqual(df[WILD_CHAIN].iloc[0], "A")
            self.assertEqual(df[WILD_COL].iloc[0], "2IEM")
            self.assertEqual(df[MUTANT_COL].iloc[0], "2EIO")
            self.assertEqual(df[WILD_AA].iloc[0], "I")
            self.assertEqual(df[WILD_SEQ_NUM].iloc[0], "89")
            self.assertEqual(df[MUT_AA].iloc[0], "A")
