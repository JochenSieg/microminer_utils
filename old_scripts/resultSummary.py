"""
THIS IS WORK IN PROGRESS
"""

import logging
import os
from pathlib import Path
import argparse
import pandas as pd
import sys


parser = argparse.ArgumentParser(
    description=
    """ """)

parser.add_argument('--resultstatistics', '-r', nargs='+', required=True, type=str,
                    help='Path(s) to resultStatistic files.')
parser.add_argument('--outdir', '-o', default=os.getcwd(), type=str, help='Path to output directory')

args = parser.parse_args()

input_files = [Path(f) for f in args.resultstatistics]
outdir = Path(args.outdir)

for f in input_files:
    if not f.is_file():
        print('Error: Specified input file does not exist.')
        sys.exit(1)

if not outdir.is_dir():
    print('Error: Specified output directory does not exist or is not a directory.')
    sys.exit(1)


dfs = []
for input_file in input_files:
    df = pd.read_csv(input_file, header=0, sep='\t')

    if df.shape[0] < 1:
        print('Error: Empty file: {}'.format(input_file))
        sys.exit(1)


    #sort by wild protein infile id
    # df.sort_values(by=['wildPos'], inplace=True)

    # drop duplicates: make the found aa types for each wild position unique
    df.drop_duplicates(['wildName', 'wildAA', 'wildChain', 'wildPos', 'mutantAA'], inplace=True)

    dfs.append(df)

    # print(df)
df = pd.concat(dfs)

with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    print(df)


print('Nof found mutations: ', df.shape[0])

g = df.groupby(['wildName'])
print('Mean Nof mutations per protein: {} (std {})'.format(g.size(), 0))

# hier noch andere sachen zaehlen. Aber zu erst den Vergleich zwischenwieviele Mutationen in PRothemr sind und
# wieviele ich finde.


