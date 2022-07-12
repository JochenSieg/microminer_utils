"""
Merge CSV files by columns. More or less a wrapper around pandas merge functionality.
"""

import os
from pathlib import Path
import argparse
import sys
import logging
import functools
import pandas as pd

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description=
        """
        Merge CSV files by columns. 
        """)

    parser.add_argument('--csv', '-c', required=True, type=str, nargs='+',
                        help='List of CSV files to merge. CSVs are merged in the given order.')
    parser.add_argument('--outfile', '-o', required=True, type=str,
                        help='Path to output file')
    parser.add_argument('--cols', default=['id1', 'id2'], type=str, nargs='+',
                        help='Column names to merge CSVs on.')
    parser.add_argument('--logfile', default=Path(os.getcwd()) / 'log.log', type=Path,
                        help='Path to output file')

    args = parser.parse_args()

    csv_files = args.csv
    outfile = Path(args.outfile)
    col_names = args.cols
    logfile = Path(args.logfile)

    csv_files = [Path(_) for _ in csv_files]

    logging.basicConfig(filename=str(logfile.resolve()),
                        level=logging.INFO,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    logger.info(f'Starting scripts: {" ".join(sys.argv)}')

    df_list = [pd.read_csv(p, sep='\t', header=0) for p in csv_files]
    df_merged = functools.reduce(lambda left, right: pd.merge(left, right, on=[col_names],
                                                              how='left'), df_list)
    df_merged.to_csv(outfile, sep='\t', header=True, index=False)


if __name__ == "__main__":
    main()
