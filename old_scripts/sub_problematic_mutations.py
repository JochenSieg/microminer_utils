import argparse
import pandas as pd
import sys
from pathlib import Path
from helper import CONFIG, WILD_COL, MUTANT_COL, MUT_AA, WILD_SEQ_NUM, WILD_AA, calc_tmalign


parser = argparse.ArgumentParser(
    description=
    """
    Substract manually review problematic/erroneous mutations from mutation CSV.
    This script removes all rows from the input CSV that are in the problematic_mutations.csv.
    """)

parser.add_argument('--csv', '-c', required=True, type=str, help='Path to the csv file.')
parser.add_argument('--outfile', '-o', required=True, type=str, help='Output CSV')

args = parser.parse_args()

input_file = Path(args.csv)
outfile = Path(args.outfile)

if not input_file.is_file():
    print('Error: Specified input file does not exist.')
    sys.exit(1)

prob_csv_path = Path(CONFIG['DATA']['PROBLEMATIC_MUTATIONS'])
if not prob_csv_path.is_file():
    print('Error: Can not find CSV file with problematic mutations at {}'.format(prob_csv_path))
    sys.exit(1)

relevant_cols = [WILD_COL, MUTANT_COL, MUT_AA, WILD_SEQ_NUM, WILD_AA]


df_mutations = pd.read_csv(input_file, sep=None, header=0, engine='python')
# using "slow" python engine to be able to read any csv file with sep=None

if not all(_ in df_mutations.columns for _ in relevant_cols):
    print('Error: Input CSV does not contain necessary columns: {}'.format(relevant_cols))
    sys.exit(1)

df_problem = pd.read_csv(prob_csv_path, sep='\t')

if not all(_ in df_problem.columns for _ in relevant_cols):
    print('Error: CSV with problematic mutations does not contain'
          ' necessary columns: {}'.format(relevant_cols))
    sys.exit(1)

# drop duplicates if they are any to avoid problems in merging
df_problem.drop_duplicates(subset=relevant_cols, inplace=True)

df_all = df_mutations.merge(df_problem, on=relevant_cols,
                            how='left', indicator=True)

# assert that all the PROBLEMATIC/WRONG mutations are also not predicted by us.
if 'mutscreen_pair_found' in df_all.columns:
    assert(df_all.query('_merge == "both"')['mutscreen_pair_found'].sum() == 0)
if 'mutscreen_screen_found' in df_all.columns:
    assert(df_all.query('_merge == "both"')['mutscreen_screen_found'].sum() == 0)

# get all rows in df_mutations without the rows that are in df_problematic
df_sub = df_all[df_all['_merge'] == 'left_only']

df_sub.to_csv(outfile, sep='\t', index=False)

print('{} of {} were "problematic" mutations'.format(
    df_mutations.shape[0] - df_sub.shape[0], df_mutations.shape[0]))
print('Resulting Table has {} mutations'.format(df_sub.shape[0]))

if 'tm_score1' in df_sub.columns and 'tm_score2' in df_sub.columns:
    print('\tfrom that {} / {} have TM-score < 0.5'.format(
        df_sub.query('tm_score1 < 0.5 or tm_score2 < 0.5').shape[0],
        df_sub.shape[0]))
    print('\t removing mutations with TM-score < 0.5...')
    df_sub = df_sub.query('tm_score1 >= 0.5 and tm_score2 >= 0.5')

if 'mutscreen_pair_found' in df_all.columns:
    print('MUTSCREEN-PAIR found {} of {} mutations'.format(
        df_sub['mutscreen_pair_found'].sum(), df_sub.shape[0]))
if 'mutscreen_screen_found' in df_all.columns:
    print('MUTSCREEN-SCREEN found {} of {} mutations'.format(
        df_sub['mutscreen_screen_found'].sum(), df_sub.shape[0]))
