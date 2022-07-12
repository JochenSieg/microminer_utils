"""
Calculates all global similarities with TM-Align on MUTSCREEN result using the chain pairs having the mutation.
For a given directory all resultStatistic files are extracted and combined to a single table. On this
table all wild/mutant structure pairs are subjected to TM-Align alignment. The goal is to compare the
local similarity of the micro-environments with the global similarity.
"""
import argparse
import pandas as pd
import sys
from pathlib import Path

import helper
from helper import calc_tmalign
from helper.utils import scantree, rename_to_lib_standard, get_monomer_pdb_file_path


parser = argparse.ArgumentParser(
    description=
    """
    Calc TM-Align for all wild/mutant PDB pairs in all resultStatistic.csv files in a given directory (also in subdirs).
    """)

parser.add_argument('--dir', '-d', required=True, type=str, help='Path to directory.')
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
    if p.is_file():
        if p.name == 'resultStatistic.csv':
            df_list.append(pd.read_csv(p, sep='\t'))

df = pd.concat(df_list)

if df.shape[0] == 0:
    print('Error: No resultStatistic.csv in input dir (and not in subdirs of input dir).')
    sys.exit(1)

df['wild_with_chain'] = df['wildName'].str.upper() + '_' + df['wildChain'].astype(str)
df['mutant_with_chain'] = df['mutantName'].str.upper() + '_' + df['mutantChain'].astype(str)

df = rename_to_lib_standard(df, wild='wild_with_chain', mutant='mutant_with_chain')

print('Read {} ({} unique) wild/mutant structure pairs from {} CSV files from input dir: {}'.format(
    df.shape[0], df.drop_duplicates([helper.WILD_COL, helper.MUTANT_COL]).shape[0], len(df_list), input_dir))

# TODO since I do not rebuild the sienaDB everytime the PDB mirror changes there are cases that are now obsolete
# all_pdbids = set(list(df[helper.WILD_COL].unique()) + list(df[helper.MUTANT_COL].unique()))
# missing_pdbids = [_.stem.upper() for _ in scantree(helper.CONFIG['DATA']['MONOMER_DIR'])
#                   if _.suffixes[-1] == '.pdb' and not _.stem.upper() in all_pdbids]
# if len(missing_pdbids) > 0:
#     df = df[~df[helper.WILD_COL].isin(missing_pdbids)]
#     df = df[~df[helper.MUTANT_COL].isin(missing_pdbids)]
#     print('Warning: Removed {} PDB codes not existing in PDB mirror'.format(len(missing_pdbids)))

df = calc_tmalign(df, cpus=cpus, use_monomers=True)
df.to_csv(outfile, sep='\t', index=False)

nof_dissimilar = df.query('tm_score1 < 0.5 or tm_score2 < 0.5').shape[0]
print(nof_dissimilar, '/', df.shape, ' have TM-score < 0.5 (are likely not related at all)')
