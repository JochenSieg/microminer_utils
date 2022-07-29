import tempfile
import unittest
from pathlib import Path

import pandas as pd

from helper import WILD_COL, WILD_AA, WILD_CHAIN, WILD_SEQ_NUM, MUT_AA
from helper.datasets.fireprotdb import FireProtDB


class FireProtDBTests(unittest.TestCase):
    """Test handling FireProtDB data set"""

    def test_read_single_mutations(self):
        """Test to read single mutations from FireProtDB CSV"""
        with tempfile.NamedTemporaryFile(mode="wt") as fireprotdb_csv:
            fireprotdb_csv.write(
                "\n".join(
                    [
                        "experiment_id,protein_name,uniprot_id,pdb_id,chain,position,wild_type,mutation,"
                        "ddG,dTm,is_curated,type,derived_type,interpro_families,conservation,"
                        "is_essential,correlated_positions,is_back_to_consensus,secondary_structure,"
                        "asa,is_in_catalytic_pocket,is_in_tunnel_bottleneck,b_factor,method,"
                        "method_details,technique,technique_details,pH,tm,notes,publication_doi,"
                        "publication_pubmed,hsw_job_id,datasets,sequence",
                        "LL000123,,,1CWQ,A,254,P,F,,,,,,,2,,,,E,,true,true,,,,,,,42.42,,,,,,",
                        "LL000321,,,2RN2|2RN2,B,77,L,R,,,,,,,2,,,,E,,true,true,,,,,,,,,,,,,",
                        "LL000322,,,2RN2|1G9V|1PX1,Z,55,S,R,,,,,,,2,,,,E,,true,true,,,,,,,,,,,,,",
                        "LL000323,,,,B,55,S,R,,,,,,,2,,,,E,,true,true,,,,,,,,,,,,,",  # missing PDBid
                    ]
                )
            )
            fireprotdb_csv.flush()

            exp_df = pd.DataFrame(
                [
                    ["1CWQ", "A", "254", "P", "F"],
                    ["2RN2", "B", "77", "L", "R"],
                    ["2RN2", "Z", "55", "S", "R"],
                    ["1G9V", "Z", "55", "S", "R"],
                    ["1PX1", "Z", "55", "S", "R"],
                ],
                columns=[WILD_COL, WILD_CHAIN, WILD_SEQ_NUM, WILD_AA, MUT_AA],
            )

            dataset = FireProtDB(Path(fireprotdb_csv.name))

            df = dataset.read_single_mutations(pdb_mutant_only=False)
            self.assertTrue(
                df[[WILD_COL, WILD_CHAIN, WILD_SEQ_NUM, WILD_AA, MUT_AA]].equals(exp_df)
            )
            df = dataset.read_single_mutations(pdb_mutant_only=True)
            self.assertTrue(df.equals(pd.DataFrame()))
