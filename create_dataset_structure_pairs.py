import os
from pathlib import Path
import argparse
import sys
import logging

import pandas as pd

import helper
from helper import BAD_PDBIDS, SUPPORTED_DATASETS, MUTATION_DATASETS, \
    DATASETS_WITH_CUSTOM_STRUCTUREFILES, MUTATION_DATASETS_WITH_STRUCTURE_PAIRS
from helper.datasets import read_mutation_dataset_single
from helper.utils import scantree

logger = logging.getLogger(__name__)

PATH_COL1 = 'structure_path1'
PATH_COL2 = 'structure_path2'
ID_COL1 = 'id1'
ID_COL2 = 'id2'


def make_from_resultStatistics(csv_input, backward, outfile, pdb_mirror1, pdb_mirror2):
    # gather all 'resultStatistic.csv' in the csv_input list (including recursive read of dirs)
    files = []
    for path in csv_input:
        files.extend([p for p in helper.utils.scantree(path) if p.name == 'resultStatistic.csv'])
    if len(files) == 0:
        print('Error: No resultStatistic.csv in input (and not in subdirs of any input dir).')
        sys.exit(1)
    logger.info(f'Gathered {len(files)} input resultStatistic.csv files.')

    # collect all csv files in a single dataframe.
    # Careful this could get large. Maybe this needs to be processed in chunks in the future if we
    # move to millions/billions of structure pairs.
    interesting_cols = ['queryName', 'hitName']
    df = pd.concat([pd.read_csv(p, sep='\t', header=0, usecols=interesting_cols).drop_duplicates(
        subset=interesting_cols) for p in files]).drop_duplicates(
        subset=interesting_cols)
    if df.shape[0] == 0:
        print('Error: resultStatistic.csv input files are empty.')
        sys.exit(1)
    logger.info(f'Gathered {df.shape[0]} input structure pairs.')

    df = df.rename(columns={'queryName': ID_COL1, 'hitName': ID_COL2})

    id_col1, id_col2 = ID_COL1, ID_COL2
    if backward:
        pdb_mirror2, pdb_mirror1 = pdb_mirror1, pdb_mirror2
        id_col1, id_col2 = id_col2, id_col1

    df[PATH_COL1] = df[id_col1].apply(
        lambda val: helper.utils.get_pdb_file_path(val, mirror=pdb_mirror1,
                                                   allow_obsolete=True,
                                                   ignore_missing=True))
    df[PATH_COL2] = df[id_col2].apply(
        lambda val: helper.utils.get_pdb_file_path(val, mirror=pdb_mirror2,
                                                   allow_obsolete=True,
                                                   ignore_missing=True))

    df_out = df[[id_col1, PATH_COL1, id_col2, PATH_COL2]]

    # remove duplicate ROWS, i.e. when ID and PATH are identical. Would lead to redundant
    # computations.
    df_out = df_out.drop_duplicates()
    if df_out.isna().sum(axis=1).sum() > 0:
        logger.info(f'Dropping {df_out.isna().sum(axis=1).sum()} missing values'
                    f' (of {df_out.shape[0]} total)')
        df_out = df_out.dropna()  # remove missing values, e.g. when only CIF file is available

    df_out.to_csv(outfile, sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser(
        description=
        """
        Creates a "dataset" CSV file of protein structure pair files. 
        
        Input can be one or multiple CSV output file(s) of MicroMiner (resultStatistic.csv)
        using --csv or 
        a .
        
        Can be used for TMAlign.
        """)

    parser.add_argument('--csv', '-c', required=True, type=str, nargs='+',
                        help='resultStatistic.csv file or dir in which resultStatistic.csv files'
                             ' will be searched recursively')
    parser.add_argument('--outfile', '-o', type=str, required=True,
                        help='Path to output file.')
    parser.add_argument('--backward', '-b', default=False, action='store_true',
                        help='')
    parser.add_argument('--pdb_mirror1', required=True, type=str.lower,
                        choices=DATASETS_WITH_CUSTOM_STRUCTUREFILES + ['standard'],
                        default='standard',
                        help='First PDB file mirror to look up protein structure files by IDs.'
                             ' Used for first column in result CSV (query column).')
    parser.add_argument('--pdb_mirror2', required=True, type=str.lower,
                        choices=DATASETS_WITH_CUSTOM_STRUCTUREFILES + ['standard'],
                        default='standard',
                        help='Second PDB file mirror to look up protein structure files by IDs.'
                             ' User for second column in result CSV (hit column).')

    args = parser.parse_args()

    csv_input = args.csv
    outfile = Path(args.outfile)
    backward = args.backward
    pdb_mirror1 = args.pdb_mirror1
    pdb_mirror2 = args.pdb_mirror2

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    logger.info(f'Starting scripts: {" ".join(sys.argv)}')

    csv_input = [Path(_) for _ in csv_input]
    make_from_resultStatistics(csv_input, backward, outfile, pdb_mirror1, pdb_mirror2)


if __name__ == "__main__":
    main()
