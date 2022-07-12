import argparse
import pandas as pd
import sys
from pathlib import Path

import helper
from helper import calc_tmalign


def has_needed_columns(df:pd.DataFrame) -> bool:
    """
    Checks whether needed columns are present. The two columns that specify the
    PDB ids of the protein pair to align are needed.
    :param df: the table
    :return: true if has needed columns. False otherwise.
    """
    return helper.WILD_COL in df.columns and helper.MUTANT_COL in df.columns


def get_convertable_column_names(df: pd.DataFrame) -> tuple:
    """
    Check whether there are columns present that can be converted to the
    needed columns.
    :param df: the table
    :return:
    """
    if 'wildName' in df.columns and 'mutantName' in df.columns:
        return 'wildName', 'mutantName'
    return tuple()


def convert_column_names_if_needed(df: pd.DataFrame) -> pd.DataFrame:
    if has_needed_columns(df):
        return df
    conv_col_names = get_convertable_column_names(df)
    if conv_col_names:
        return df.rename({conv_col_names[0]: helper.WILD_COL, conv_col_names[1]: helper.MUTANT_COL},
                         axis=1, errors='raise')
    return pd.DataFrame()


parser = argparse.ArgumentParser(
    description=
    """
    Calculates TMalign 3D structure alignment between wild and mutant structures
    and annotated the TM-Scores in a CSV.
    This script simply adds some columns about protein structure/sequence similarity
    based on TM-Align.
    """)

parser.add_argument('--csv', '-c', required=True, type=str, help='Path to the csv file.')
parser.add_argument('--outfile', '-o', required=True, type=str, help='Output CSV')

args = parser.parse_args()

input_file = Path(args.csv)
outfile = Path(args.outfile)

if not input_file.is_file():
    print('Error: Specified input file does not exist.')
    sys.exit(1)

df_mutations = pd.read_csv(input_file, sep=None, header=0, engine='python')
# using "slow" python engine to be able to read any csv file with sep=None

df_mutations = convert_column_names_if_needed(df_mutations)
if df_mutations.shape[0] == 0:
    print('Error: Provided table empty or invalid format')
    sys.exit(1)

df = calc_tmalign(df_mutations)
df.to_csv(outfile, sep='\t', index=False)

nof_dissimilar = df.query('tm_score1 < 0.5 or tm_score2 < 0.5').shape[0]
print(nof_dissimilar, '/', df.shape, ' have TM-score < 0.5 (are likely not related at all)')
