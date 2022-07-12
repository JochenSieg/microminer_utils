"""
Calculate statistics over several resultStatistic.csv file from MutScreen.

For example counts of:
- mean number of found mutations per wild-type query
"""
import argparse
import pandas as pd
import numpy as np
import sys
from pathlib import Path

import helper
from helper import calc_tmalign
from helper.utils import scantree, rename_to_lib_standard


parser = argparse.ArgumentParser(
    description=
    """
    Calculate statistics for all resultStatistic.csv files in a directory tree.
    """)

parser.add_argument('--dir', '-d', required=True, type=str, help='Path to directory from which to read '
                                                                 'resultStatistic.csv files recursively.')
parser.add_argument('--outfile', '-o', required=True, type=str, help='Output CSV')
parser.add_argument('--cpus', '-c', default=1, type=int, help='number of processes to use for parallel execution')

args = parser.parse_args()

input_dir = Path(args.dir)
outfile = Path(args.outfile)
cpus = args.cpus

if not input_dir.is_dir():
    print('Error: Specified input dir does not exist.')
    sys.exit(1)

cpus = max(cpus, 1)

df_list = []
for p in scantree(input_dir):
    if p.is_file() and p.name == 'resultStatistic.csv':
        df_list.append(pd.read_csv(p, sep='\t'))

df = pd.concat(df_list)

if df.shape[0] == 0:
    print('Error: No resultStatistic.csv in input dir (and not in subdirs of input dir).')
    sys.exit(1)

# prepare for TMalign on the chains with the mutation

# df['wild_with_chain'] = df['wildName'].str.upper() + '_' + df['wildChain'].astype(str)
# df['mutant_with_chain'] = df['mutantName'].str.upper() + '_' + df['mutantChain'].astype(str)
#
# df = rename_to_lib_standard(df, wild='wild_with_chain', mutant='mutant_with_chain')
#
# print('Read {} ({} unique) wild/mutant structure pairs from {} CSV files from input dir: {}'.format(
#     df.shape[0], df.drop_duplicates([helper.WILD_COL, helper.MUTANT_COL]).shape[0], len(df_list), input_dir))
#
# df = calc_tmalign(df, cpus=cpus, use_monomers=True)

# calculate statistics
wildquery_groupby = df.groupby(['wildName'])

list_nof_group = []  # just count mutations per wild-type (including duplicates i.e. mutation to the same amino acid)
list_nof_group_pre_pos_u = []  # count the number of unique mutations
list_nof_group_pre_pos = []  # count the number of wild positions for which at least 1 mutation was found.

for name, group in wildquery_groupby:
    list_nof_group.append(group.shape[0])

    group_pre_pos_u = group.drop_duplicates(['wildChain', 'wildPos', 'wildAA', 'mutantAA'])
    list_nof_group_pre_pos_u.append(group_pre_pos_u.shape[0])

    group_pre_pos = group.drop_duplicates(['wildChain', 'wildPos', 'wildAA'])
    list_nof_group_pre_pos.append(group_pre_pos.shape[0])

print('Number of query structures: {}'.format(len(df['wildName'].unique())))
print('Number of mutations found (number of rows): {}'.format(df['wildName'].shape[0]))

print('ALL MUTATIONS: abs number: {}, mean: {}, median: {}, std: {}'.format(
    np.sum(list_nof_group), np.mean(list_nof_group), np.median(list_nof_group), np.std(list_nof_group)))
print('MUTATIONS FOR WILD POSITIONS (UNIQ MUTANT AA): abs number: {}, mean: {}, median: {}, std: {}'.format(
    np.sum(list_nof_group_pre_pos_u), np.mean(list_nof_group_pre_pos_u), np.median(list_nof_group_pre_pos_u),
    np.std(list_nof_group_pre_pos_u)))
print('WILD POSITIONS WITH MUTATIONS: abs number: {}, mean: {}, median: {}, std: {}'.format(
    np.sum(list_nof_group_pre_pos), np.mean(list_nof_group_pre_pos), np.median(list_nof_group_pre_pos), np.std(list_nof_group_pre_pos)))

# df.to_csv(outfile, sep='\t', index=False)
#
# nof_dissimilar = df.query('tm_score1 < 0.5 or tm_score2 < 0.5').shape[0]
# print(nof_dissimilar, '/', df.shape, ' have TM-score < 0.5 (are likely not related at all)')
