
# TODO aktuell nutze ich dieses Skript nicht da es zu einer kombinatorischen Explosion der Paare fuehrt. Unpraktisch.

import os
from pathlib import Path
import argparse
import sys
import logging
from itertools import product
import pandas as pd

logger = logging.getLogger(__name__)

PATH_COL1 = 'structure_path1'
PATH_COL2 = 'structure_path2'
ID_COL1 = 'id1'
ID_COL2 = 'id2'


def combine_to_all_against_all_pairs(csv_path1: Path, csv_path2: Path, outfile: Path):
    usecols = ['id', 'structure_path']
    df1 = pd.read_csv(csv_path1, sep='\t', header=0, usecols=usecols)
    df2 = pd.read_csv(csv_path2, sep='\t', header=0, usecols=usecols)

    df = pd.DataFrame(product(df1['id'], df2['id']), columns=[ID_COL1, ID_COL2])
    df[PATH_COL1], df[PATH_COL2] = zip(product(df1['structure_path'], df2['structure_path']))
    df.to_csv(outfile, sep='\t', header=True, index=False)


def main():
    parser = argparse.ArgumentParser(
        description=
        """
        Creates a "dataset" CSV file of protein structure pair files.
        
        Input: Two CSV files of with id and	structure_path columns each.
        Output: A single CSV file of id1, structure_path1, id2, structure_path2 pairs.

        Can be used for TMAlign.
        """)

    parser.add_argument('--csv1', required=True, type=str,
                        help='First CSV file.')
    parser.add_argument('--csv2', required=True, type=str,
                        help='Second CSV file.')
    parser.add_argument('--outfile', '-o', type=str, required=True,
                        help='Path to output file.')

    args = parser.parse_args()

    csv_path1 = args.csv1
    csv_path2 = args.csv2
    outfile = Path(args.outfile)

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    logger.info(f'Starting scripts: {" ".join(sys.argv)}')

    combine_to_all_against_all_pairs(csv_path1, csv_path2, outfile)


if __name__ == "__main__":
    main()
