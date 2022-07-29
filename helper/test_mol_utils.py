import unittest
import tempfile

from pathlib import Path

import pandas as pd

from helper.mol_utils import read_plddt_values


class MolUtilsTests(unittest.TestCase):
    """Test mol utils"""

    pdb_str = """ATOM    101  N   ASP A  16      -3.771   1.435 -17.204  1.00 94.55           N  
ATOM    102  CA  ASP A  16      -2.456   0.874 -16.935  1.00 94.55           C  
ATOM    103  C   ASP A  16      -1.828   1.573 -15.728  1.00 94.55           C  
ATOM    104  CB  ASP A  16      -1.565   0.957 -18.184  1.00 94.55           C  
ATOM    105  O   ASP A  16      -1.354   2.712 -15.796  1.00 94.55           O  
ATOM    106  CG  ASP A  16      -0.241   0.213 -17.991  1.00 94.55           C  
ATOM    107  OD1 ASP A  16       0.078  -0.196 -16.846  1.00 94.55           O  
ATOM    108  OD2 ASP A  16       0.533   0.103 -18.959  1.00 94.55           O  
ATOM    109  N   ALA A  17      -1.788   0.861 -14.602  1.00 93.02           N  
ATOM    110  CA  ALA A  17      -1.170   1.354 -13.378  1.00 93.02           C  
ATOM    111  C   ALA A  17       0.344   1.593 -13.535  1.00 93.02           C  
ATOM    112  CB  ALA A  17      -1.471   0.364 -12.249  1.00 93.02           C  
ATOM    113  O   ALA A  17       0.883   2.457 -12.843  1.00 93.02           O  
ATOM    114  N   SER A  18       1.021   0.889 -14.454  1.00 94.50           N  
ATOM    115  CA  SER A  18       2.436   1.130 -14.773  1.00 94.50           C  
ATOM    116  C   SER A  18       2.656   2.488 -15.436  1.00 94.50           C  
ATOM    117  CB  SER A  18       3.007   0.023 -15.672  1.00 94.50           C  
ATOM    118  O   SER A  18       3.715   3.081 -15.274  1.00 94.50           O  
ATOM    119  OG  SER A  18       2.717   0.196 -17.041  1.00 94.50           O  """

    def test_read_plddt(self):
        """Test read pLDDT from PDB file"""

        with tempfile.NamedTemporaryFile(mode="wt") as file:
            path = Path(file.name)
            file.write(MolUtilsTests.pdb_str)
            file.flush()

            # add column with name which is the tmp files basename
            exp_df = pd.DataFrame(
                {
                    "name": path.stem.split(".")[0],
                    "aa": ["ASP", "ALA", "SER"],
                    "chain": ["A", "A", "A"],
                    "pos": ["16", "17", "18"],
                    "plddt": ["94.55", "93.02", "94.50"],
                }
            )

            df = read_plddt_values(path)
            self.assertTrue(df.equals(exp_df))
