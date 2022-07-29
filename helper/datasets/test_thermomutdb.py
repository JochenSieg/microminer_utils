import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from helper import WILD_COL, WILD_AA, WILD_CHAIN, WILD_SEQ_NUM, MUTANT_COL, MUT_AA
from helper.datasets.thermomutdb import ThermoMutDB


class ThermoMutDBTests(unittest.TestCase):
    """Test handling ThermoMutDB data set"""

    def test_read_single_mutations(self):
        """Test to read single mutations from ThermoMutDB CSV"""
        with tempfile.NamedTemporaryFile(mode="wt") as thermomutdb_csv:
            data = [
                {
                    "pos": 0,
                    "pdbs_template": "null",
                    "PDBs_template": "null",
                    "id": 1234,
                    "PDB_wild": "3DRI",
                    "mut_count": 0,
                    "pdb_mutant": "null",
                    "mutation_type": "Single",
                    "mutation_code": "L63C",
                    "mutation_based": "PDB",
                    "mutated_chain": "A",
                },
                {
                    "pos": 0,
                    "pdbs_template": "null",
                    "PDBs_template": "null",
                    "id": 1234,
                    "PDB_wild": "2DRA",
                    "mut_count": 0,
                    "pdb_mutant": "2RN2",
                    "mutation_type": "Single",
                    "mutation_code": "L123C",
                    "mutation_based": "PDB",
                    "mutated_chain": "B",
                },
                {
                    "pos": "null",
                    "pdbs_template": "null",
                    "PDBs_template": "null",
                    "id": 4321,
                    "PDB_wild": "1EEN",
                    "mut_count": 1,
                    "pdb_mutant": "null",
                    "mutation_type": "Multiple",
                    "mutation_code": "I42A, V70I",
                    "mutated_chain": "A",
                },
                {
                    "pos": 0,
                    "PDB_wild": "1PGA",
                    "pdb_mutant": "2KLK, 2RMM",
                    "mutation_type": "Single",
                    "mutation_code": "A34H",
                    "mutation_based": "unsigned",
                    "rsa": 0.01,
                    "effect": "destabilizing",
                    "mutated_chain": "A",
                    "res_depth": 3.05,
                    "temperature": 310.15,
                },
                {
                    # this entry should be filtered out because of the invalid chain id.
                    # this is a real case.
                    "pos": "null",
                    "PDB_wild": "2RN5",
                    "pdb_mutant": "1KVB",
                    "mutation_type": "Single",
                    "mutation_code": "D134N",
                    "mutated_chain": "unsigned",
                },
            ]
            json.dump(data, thermomutdb_csv)
            thermomutdb_csv.flush()

            dataset = ThermoMutDB(Path(thermomutdb_csv.name))

            exp_df = pd.DataFrame(
                [
                    ["3DRI", "A", "63", "L", "C"],
                    ["2DRA", "B", "123", "L", "C"],
                    ["1PGA", "A", "34", "A", "H"],
                ],
                columns=[WILD_COL, WILD_CHAIN, WILD_SEQ_NUM, WILD_AA, MUT_AA],
            )

            df = dataset.read_single_mutations(pdb_mutant_only=False)
            self.assertTrue(
                df[[WILD_COL, WILD_CHAIN, WILD_SEQ_NUM, WILD_AA, MUT_AA]].equals(exp_df)
            )

            exp_df = pd.DataFrame(
                [
                    ["2DRA", "B", "123", "L", "C", "2RN2"],
                    ["1PGA", "A", "34", "A", "H", "2KLK"],
                    ["1PGA", "A", "34", "A", "H", "2RMM"],
                ],
                columns=[
                    WILD_COL,
                    WILD_CHAIN,
                    WILD_SEQ_NUM,
                    WILD_AA,
                    MUT_AA,
                    MUTANT_COL,
                ],
            )

            df = dataset.read_single_mutations(pdb_mutant_only=True)
            self.assertTrue(
                df[
                    [WILD_COL, WILD_CHAIN, WILD_SEQ_NUM, WILD_AA, MUT_AA, MUTANT_COL]
                ].equals(exp_df)
            )
