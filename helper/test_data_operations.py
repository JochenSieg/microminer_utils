import unittest
import tempfile
from pathlib import Path

import pandas as pd
import numpy as np

from helper import WILD_COL, WILD_AA, WILD_CHAIN, WILD_SEQ_NUM, MUTANT_COL, MUT_AA
from helper.constants import (
    MM_QUERY_NAME,
    MM_QUERY_AA,
    MM_QUERY_CHAIN,
    MM_QUERY_POS,
    MM_HIT_NAME,
    MM_HIT_AA,
    MM_HIT_CHAIN,
    MM_HIT_POS,
    CONFIG,
)
from helper.data_operations import (
    make_search_parameter_table,
    make_pair_parameter_table,
    read_microminer_csv,
    merge_results_for_pair_eval,
)
from helper.datasets import PDB
from helper.datasets.dataset import MockDataset, MockMutationDataset


class DataOperationsTests(unittest.TestCase):
    """Test data operations"""

    def test_make_search_parameter_table(self):
        """test correct creation of search parameter table"""

        # empty dataset
        dataset = MockDataset()
        self.assertRaises(
            NotImplementedError, make_search_parameter_table, dataset, False
        )

        dataset = MockMutationDataset()
        dataset.read_single_mutations = lambda pdb_mutant_only: pd.DataFrame(
            {
                WILD_COL: ["1G9V", "2rn2"],
                MUTANT_COL: ["2RN2", "1g9v"],
                WILD_AA: ["H", "N"],
                WILD_SEQ_NUM: ["48", "48"],
                MUT_AA: ["N", "H"],
            }
        )

        df = make_search_parameter_table(dataset, backward=False)
        self.assertEqual(df.shape, (2, 2))
        self.assertTrue(df.columns.tolist() == ["id", "structure_path"])
        self.assertEqual(df["id"].iloc[0], "1G9V")
        self.assertTrue(Path(df["structure_path"].iloc[0]).is_file())
        self.assertEqual(df["id"].iloc[1], "2rn2")
        self.assertTrue(Path(df["structure_path"].iloc[1]).is_file())

        path1 = df["structure_path"].iloc[0]

        df = make_search_parameter_table(dataset, backward=True)
        self.assertEqual(df.shape, (2, 2))
        self.assertTrue(df.columns.tolist() == ["id", "structure_path"])
        self.assertEqual(df["id"].iloc[0], "2RN2")
        self.assertTrue(Path(df["structure_path"].iloc[0]).is_file())
        self.assertNotEqual(path1, df["structure_path"].iloc[0])
        self.assertEqual(df["id"].iloc[1], "1g9v")
        self.assertTrue(Path(df["structure_path"].iloc[1]).is_file())

        dataset = MockMutationDataset()
        dataset.read_single_mutations = lambda pdb_mutant_only: pd.DataFrame(
            {
                WILD_COL: ["1G9V"],
                WILD_AA: ["H"],
                WILD_SEQ_NUM: ["48"],
                MUT_AA: ["N"],
            }
        )
        dataset.has_structure_pairs = False
        self.assertRaises(ValueError, make_search_parameter_table, dataset, True)

        # test PDB as query DB
        with tempfile.TemporaryDirectory() as pdb_dir:
            pdb_dir_path = Path(pdb_dir)
            for pdbid in ["1g9v", "1xg1", "2rn2"]:
                filename = (
                    f"{CONFIG['PDB_FILE_INFO']['PREFIX']}{pdbid}"
                    f"{CONFIG['PDB_FILE_INFO']['SUFFIX']}"
                )
                with open(pdb_dir_path / filename, "w") as f:
                    pass

            dataset = PDB(pdb_dir_path)
            df = make_search_parameter_table(dataset, backward=False)
            self.assertEqual(df.columns.tolist(), ["id", "structure_path"])
            self.assertEqual(df["id"].tolist(), ["1G9V", "1XG1", "2RN2"])
            for i in range(df.shape[0]):
                self.assertTrue(Path(df["structure_path"].iloc[i]).is_file())

    def test_make_pair_parameter_table(self):
        """test correct creation of pair parameter table"""

        dataset = MockMutationDataset()
        dataset.read_single_mutations = lambda pdb_mutant_only: pd.DataFrame(
            {
                WILD_COL: ["1G9V", "2rn2"],
                MUTANT_COL: ["2RN2", "1g9v"],
                WILD_AA: ["H", "N"],
                WILD_SEQ_NUM: ["48", "48"],
                MUT_AA: ["N", "H"],
            }
        )

        df = make_pair_parameter_table(dataset, backward=False)
        self.assertEqual(df.shape, (2, 4))
        self.assertTrue(
            df.columns.tolist() == ["id1", "structure_path1", "id2", "structure_path2"]
        )
        self.assertEqual(df["id1"].iloc[0], "1G9V")
        self.assertEqual(df["id2"].iloc[0], "2RN2")
        self.assertTrue(Path(df["structure_path1"].iloc[0]).is_file())
        self.assertTrue(Path(df["structure_path2"].iloc[0]).is_file())
        self.assertNotEqual(
            df["structure_path1"].iloc[0], df["structure_path2"].iloc[0]
        )
        self.assertEqual(df["id1"].iloc[1], "2rn2")
        self.assertEqual(df["id2"].iloc[1], "1g9v")
        self.assertTrue(Path(df["structure_path1"].iloc[1]).is_file())
        self.assertTrue(Path(df["structure_path2"].iloc[1]).is_file())
        self.assertNotEqual(
            df["structure_path1"].iloc[1], df["structure_path2"].iloc[1]
        )

        df = make_pair_parameter_table(dataset, backward=True)
        self.assertEqual(df.shape, (2, 4))
        self.assertTrue(
            df.columns.tolist() == ["id1", "structure_path1", "id2", "structure_path2"]
        )
        self.assertEqual(df["id1"].iloc[0], "2RN2")
        self.assertEqual(df["id2"].iloc[0], "1G9V")
        self.assertTrue(Path(df["structure_path1"].iloc[0]).is_file())
        self.assertTrue(Path(df["structure_path2"].iloc[0]).is_file())
        self.assertNotEqual(
            df["structure_path1"].iloc[0], df["structure_path2"].iloc[0]
        )
        self.assertEqual(df["id1"].iloc[1], "1g9v")
        self.assertEqual(df["id2"].iloc[1], "2rn2")
        self.assertTrue(Path(df["structure_path1"].iloc[1]).is_file())
        self.assertTrue(Path(df["structure_path2"].iloc[1]).is_file())
        self.assertNotEqual(
            df["structure_path1"].iloc[1], df["structure_path2"].iloc[1]
        )

    def test_read_microminer_csv(self):
        with tempfile.NamedTemporaryFile(
            mode="wt"
        ) as mm_csv1, tempfile.NamedTemporaryFile(mode="wt") as mm_csv2:
            df1 = pd.DataFrame(
                {
                    MM_QUERY_NAME: ["1G9V"],
                    MM_QUERY_AA: ["ALA"],
                    MM_QUERY_CHAIN: ["A"],
                    MM_QUERY_POS: ["23"],
                    MM_HIT_NAME: ["2RN2"],
                    MM_HIT_AA: ["VAL"],
                    MM_HIT_CHAIN: ["B"],
                    MM_HIT_POS: ["-32a"],
                }
            )
            df1.to_csv(Path(mm_csv1.name), sep="\t", index=False)
            mm_csv1.flush()
            df2 = pd.DataFrame(
                {
                    MM_QUERY_NAME: ["1E23", "1G9V"],
                    MM_QUERY_AA: ["ILE", "ALA"],
                    MM_QUERY_CHAIN: ["C", "A"],
                    MM_QUERY_POS: ["88", "23"],
                    MM_HIT_NAME: ["8ABC", "2RN2"],
                    MM_HIT_AA: ["TYR", "VAL"],
                    MM_HIT_CHAIN: ["B", "B"],
                    MM_HIT_POS: ["99", "-32a"],
                }
            )
            df2.to_csv(Path(mm_csv2.name), sep="\t", index=False)
            mm_csv2.flush()

            df = read_microminer_csv([Path(mm_csv1.name), Path(mm_csv2.name)])
            self.assertEqual(df.shape, (2, df1.shape[1]))
            self.assertEqual(df[MM_QUERY_NAME].iloc[0], "1G9V")
            self.assertEqual(df[MM_QUERY_NAME].iloc[1], "1E23")
            self.assertEqual(df[MM_HIT_POS].iloc[0], "-32a")
            self.assertEqual(df[MM_HIT_POS].iloc[1], "99")

    def test_merge_results_for_pair_eval(self):
        # setup mutation dataset table
        df_dataset = pd.DataFrame(
            {
                WILD_COL: ["1G9V", "1E23", "1G9V"],
                MUTANT_COL: ["2RN2", "7ABC", "2RN2"],
                WILD_AA: ["H", "I", "H"],
                WILD_SEQ_NUM: ["48", "88", "48"],
                MUT_AA: ["N", "Y", "N"],
            }
        )

        df_mm = pd.DataFrame(
            {
                MM_QUERY_NAME: ["1E23", "1G9V", "2RN2", "1G9V", "1G9V"],
                MM_QUERY_AA: ["ILE", "HIS", "ASN", "HIS", "HIS"],
                MM_QUERY_CHAIN: ["C", "A", "B", "A", "A"],
                MM_QUERY_POS: ["88", "48", "123", "48", "48"],
                MM_HIT_NAME: ["8ABC", "2RN2", "1G9V", "2RN2", "1NOP"],
                MM_HIT_AA: ["TYR", "ASN", "HIS", "ASN", "PRO"],
                MM_HIT_CHAIN: ["B", "B", "A", "D", "Y"],
                MM_HIT_POS: ["99", "-32a", "48", "-32a", "272"],
                "some_col": [
                    "some_val1",
                    "some_val2",
                    "some_val3",
                    "some_val4",
                    "some_val5",
                ],
            }
        )

        # forward matching
        df_anno, df_not_found, df_merged = merge_results_for_pair_eval(
            df_dataset, df_mm, backward=False
        )

        # compare df of successful annotations with expected df
        exp_df_anno = pd.DataFrame(
            [pd.concat([df_dataset.iloc[0], df_mm.iloc[1]], axis=0)]
        )
        exp_df_anno.index = [0]
        self.assertTrue(df_anno.equals(exp_df_anno))

        # compare df of failed annotations with expected df
        exp_df_not_found = pd.DataFrame(
            [
                pd.concat(
                    [
                        df_dataset.iloc[1],
                        pd.Series(
                            np.full(len(df_mm.columns), np.nan), index=df_mm.columns
                        ),
                    ],
                    axis=0,
                )
            ]
        )
        exp_df_not_found = exp_df_not_found.astype(
            {col: object for col in df_mm.columns}
        )
        exp_df_not_found.index = [1]
        self.assertTrue(df_not_found.equals(exp_df_not_found))

        # compare df of all (successful + failed) with expected
        exp_df_merged = pd.DataFrame(
            [
                pd.concat([df_dataset.iloc[0], df_mm.iloc[1]]),
                pd.concat(
                    [
                        df_dataset.iloc[1],
                        pd.Series(
                            np.full(len(df_mm.columns), np.nan), index=df_mm.columns
                        ),
                    ],
                    axis=0,
                ),
            ]
        )
        exp_df_merged = exp_df_merged.astype({col: object for col in df_mm.columns})
        exp_df_merged.index = [0, 1]
        self.assertTrue(df_merged.equals(exp_df_merged))

        # backward matching
        df_anno, df_not_found, df_merged = merge_results_for_pair_eval(
            df_dataset, df_mm, backward=True
        )

        # compare df of successful annotations with expected df
        exp_df_anno = pd.DataFrame(
            [pd.concat([df_dataset.iloc[0], df_mm.iloc[2]], axis=0)]
        )
        exp_df_anno.index = [0]
        self.assertTrue(df_anno.equals(exp_df_anno))

        # compare df of failed annotations with expected df. Must be same as forward
        self.assertTrue(df_not_found.equals(exp_df_not_found))

        # compare df of all (successful + failed) with expected
        exp_df_merged = pd.DataFrame(
            [
                pd.concat([df_dataset.iloc[0], df_mm.iloc[2]]),
                pd.concat(
                    [
                        df_dataset.iloc[1],
                        pd.Series(
                            np.full(len(df_mm.columns), np.nan), index=df_mm.columns
                        ),
                    ],
                    axis=0,
                ),
            ]
        )
        exp_df_merged = exp_df_merged.astype({col: object for col in df_mm.columns})
        exp_df_merged.index = [0, 1]
        self.assertTrue(df_merged.equals(exp_df_merged))

        # test correct chain handling, if chain ids are available in mutation data set

        # with chain id and forward
        df_dataset[WILD_CHAIN] = ["A", "X", "A"]
        df_anno, df_not_found, df_merged = merge_results_for_pair_eval(
            df_dataset, df_mm, backward=False
        )

        # compare df of successful annotations with expected df
        exp_df_anno = pd.DataFrame(
            [pd.concat([df_dataset.iloc[0], df_mm.iloc[1]], axis=0)]
        )
        exp_df_anno.index = [0]
        self.assertTrue(df_anno.equals(exp_df_anno))

        # with chain id and backwards
        df_anno, df_not_found, df_merged = merge_results_for_pair_eval(
            df_dataset, df_mm, backward=True
        )

        # compare df of successful annotations with expected df
        exp_df_anno = pd.DataFrame(
            [pd.concat([df_dataset.iloc[0], df_mm.iloc[2]], axis=0)]
        )
        exp_df_anno.index = [0]
        self.assertTrue(df_anno.equals(exp_df_anno))

        df_dataset[WILD_CHAIN] = ["NotAChainID", "NotAChainID", "NotAChainID"]

        df_anno, df_not_found, df_merged = merge_results_for_pair_eval(
            df_dataset, df_mm, backward=False
        )
        self.assertEqual(df_anno.shape[0], 0)
