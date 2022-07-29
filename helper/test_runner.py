import tempfile
import unittest
from pathlib import Path

import pandas as pd

from helper import WILD_COL, WILD_AA, WILD_SEQ_NUM, MUTANT_COL, MUT_AA
from helper.data_operations import make_search_parameter_table
from helper.datasets.dataset import MockMutationDataset
from helper.runners import MicroMinerSearch


class RunnersTests(unittest.TestCase):
    """Test runners"""

    def test_MicroMinerSearch(self):
        """Test MM search runner"""

        dataset = MockMutationDataset()
        dataset.read_single_mutations = lambda pdb_mutant_only: pd.DataFrame(
            {
                WILD_COL: ["2rn2", "3loe"],
                MUTANT_COL: ["1law", "2rn1"],
                WILD_AA: ["V", "A"],
                WILD_SEQ_NUM: ["74", "42"],
                MUT_AA: ["ILE", "PRO"],
            }
        )

        df = make_search_parameter_table(dataset, backward=False)

        with tempfile.TemporaryDirectory() as outdir, tempfile.NamedTemporaryFile(
            mode="wt"
        ) as dataset_file:
            outdir_path = Path(outdir)

            df.to_csv(dataset_file, sep="\t", index=False)
            dataset_file.flush()

            runner = MicroMinerSearch(
                cpus=1, raise_error=True, mm_mode="single_mutation", mm_repr="monomer"
            )
            info_dict_list = runner.run(Path(dataset_file.name), outdir=outdir_path)

            self.assertIsNotNone(info_dict_list)
            self.assertEqual(len(info_dict_list), 2)
            for info_dict in info_dict_list:
                self.assertGreater(len(info_dict), 0)

            files = list(outdir_path.glob("2rn2/*resultStatistic.csv"))
            self.assertGreater(len(files), 0)
            files = list(outdir_path.glob("3loe/*resultStatistic.csv"))
            self.assertGreater(len(files), 0)

        # parallelized
        with tempfile.TemporaryDirectory() as outdir, tempfile.NamedTemporaryFile(
            mode="wt"
        ) as dataset_file:
            outdir_path = Path(outdir)

            df.to_csv(dataset_file, sep="\t", index=False)
            dataset_file.flush()

            runner = MicroMinerSearch(
                cpus=2, raise_error=True, mm_mode="single_mutation", mm_repr="monomer"
            )
            info_dict_list = runner.run(Path(dataset_file.name), outdir=outdir_path)

            self.assertIsNotNone(info_dict_list)
            self.assertEqual(len(info_dict_list), 2)
            for info_dict in info_dict_list:
                self.assertGreater(len(info_dict), 0)

            files = list(outdir_path.glob("2rn2/*resultStatistic.csv"))
            self.assertGreater(len(files), 0)
            files = list(outdir_path.glob("3loe/*resultStatistic.csv"))
            self.assertGreater(len(files), 0)
