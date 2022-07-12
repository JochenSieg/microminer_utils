"""
Evaluates MicroMiner results obtained on the SCOPe data set.
"""

import os
from pathlib import Path
import argparse
import sys
import logging
from typing import List, TextIO
import pandas as pd

import helper
from helper import constants

logger = logging.getLogger(__name__)


def generate_plots(df: pd.DataFrame):
    pass  # TODO erst im jupyter notebook und dann hier automatisieren.


def add_tmalign_results(tmalign_input):
    tmalign_input = [Path(_) for _ in tmalign_input]
    tmalign_files = collect_files(tmalign_input, suffix='_tmalign.tsv')
    logger.info(f'Gathered {len(tmalign_files)} input TMAlign result files.')
    col_dtypes = {}
    with helper.utils.timer('Reading TMAlign result CSVs'):
        df_tm = pd.concat([pd.read_csv(p, sep='\t', header=0, dtype=col_dtypes)
                           for p in tmalign_files]).drop_duplicates()
    if df_tm.shape[0] == 0:
        print('WARNING: TMAlign input files are empty.')
        sys.exit(1)
    logger.info(f'Gathered {df_tm.shape[0]} TMAlign result rows.')
    return df_tm


def print_general_scope_statistic(df_cla, f: TextIO):
    df_cla['scope_class'] = df_cla['sccs'].str.split('.').str[0]
    df_cla['scope_fold'] = df_cla['sccs'].str.split('.').str[1]
    df_cla['scope_superfamily'] = df_cla['sccs'].str.split('.').str[2]
    df_cla['scope_family'] = df_cla['sccs'].str.split('.').str[3]
    f.write(f'{df_cla["sid"].unique().shape[0]} different domains\n')
    f.write(f'{df_cla["pdbid"].unique().shape[0]} PDB entries\n')
    # unique categories in class, folds, ... etc.
    f.write(f'{df_cla["sccs"].unique().shape[0]} sccs ids\n')
    f.write(f'{df_cla["scope_class"].unique().shape[0]} classes\n')
    f.write(f'{df_cla["scope_fold"].unique().shape[0]} folds\n')
    f.write(f'{df_cla["scope_superfamily"].unique().shape[0]} superfamilies\n')
    f.write(f'{df_cla["scope_family"].unique().shape[0]} families\n')


def write_scope_stats(df_mm: pd.DataFrame, outdir: Path):
    df_mm['scope_class1'] = df_mm['sccs1'].str.split('.').str[0]
    df_mm['scope_fold1'] = df_mm['sccs1'].str.split('.').str[1]
    df_mm['scope_superfamily1'] = df_mm['sccs1'].str.split('.').str[2]
    df_mm['scope_family1'] = df_mm['sccs1'].str.split('.').str[3]
    df_mm['scope_class2'] = df_mm['sccs2'].str.split('.').str[0]
    df_mm['scope_fold2'] = df_mm['sccs2'].str.split('.').str[1]
    df_mm['scope_superfamily2'] = df_mm['sccs2'].str.split('.').str[2]
    df_mm['scope_family2'] = df_mm['sccs2'].str.split('.').str[3]
    df_mm['same_sccs'] = df_mm['sccs1'] == df_mm['sccs2']
    df_mm['same_class'] = df_mm['scope_class1'] == df_mm['scope_class2']
    df_mm['same_fold'] = df_mm['scope_fold1'] == df_mm['scope_fold2']
    df_mm['same_superfamily'] = df_mm['scope_superfamily1'] == df_mm['scope_superfamily2']
    df_mm['same_family'] = df_mm['scope_family1'] == df_mm['scope_family2']

    with open(outdir / 'microminer_statistics.txt', 'w') as f:
        # read in the whole SCOPe classification file
        df_cla = helper.datasets.scope.read_scope_classification_file()
        f.write('SCOPe general statistic:\n')
        print_general_scope_statistic(df_cla, f)
        # We use the 40% identity cluster version of SCOPe
        df_cla40 = helper.datasets.scope.read_scope()
        f.write('--------------------\n')

        f.write('SCOPe 40%Ident statistic:\n')
        print_general_scope_statistic(df_cla40, f)
        f.write('--------------------\n')

        # write MicroMiner performance
        f.write(f'Same sccs {df_mm["same_sccs"].sum()/df_mm.shape[0]*100:.0f}%\n')
        f.write(f'Same class {df_mm["same_class"].sum()/df_mm.shape[0]*100:.0f}%\n')
        f.write(f'Same fold {df_mm["same_fold"].sum()/df_mm.shape[0]*100:.0f}%\n')
        f.write(f'Same superfamily {df_mm["same_superfamily"].sum()/df_mm.shape[0]*100:.0f}%\n')
        f.write(f'Same family {df_mm["same_family"].sum()/df_mm.shape[0]*100:.0f}%\n')

        # TODO ich glaube wie ich das mit der CM mache ist quatsch.
        #      Die Grundwahrheit ist gerade nicht das gesamte scope40 set!
        f.write('--------------------\n')
        from sklearn.metrics import confusion_matrix
        cm = confusion_matrix(df_mm['sccs2'].values, df_mm['sccs1'].values)
        f.write(f'SCCS: TN={cm[0,0]} FP={cm[0,1]} FN={cm[1,0]} TP={cm[1,1]}')
        cm = confusion_matrix(df_mm['scope_class2'].values, df_mm['scope_class1'].values)
        f.write(f'CLASS: TN={cm[0,0]} FP={cm[0,1]} FN={cm[1,0]} TP={cm[1,1]}')


def collect_files(paths: List[Path], suffix: str) -> List[Path]:
    files = []
    for path in paths:
        files.extend([p for p in helper.utils.scantree(path) if p.name.endswith(suffix)])
    if len(files) == 0:
        print(f'WARNING: No suffix "{suffix}" in input (and not in subdirs of any input dir).')
        # sys.exit(1)
    return files


def main():
    parser = argparse.ArgumentParser(
        description=
        """
        Evaluates MicroMiner results obtained on the SCOPe data set.
        """)

    parser.add_argument('--microminer', required=True, type=str, nargs='+',
                        help='Dir or CSV of all MicroMiner results')
    parser.add_argument('--tmalign', required=False, type=str, nargs='+',
                        help='Dir or CSV of all TMAlign results')
    # parser.add_argument('--scope_cla', required=True, type=str, nargs='+',
    #                     help='dir.cla.scope CSV file. Can be downloaded from SCOPe services.'
    #                          ' Contains full classification for each domain in SCOPe.')
    parser.add_argument('--outdir', '-o', default=os.getcwd(), type=str,
                        help='Path to output dir')

    args = parser.parse_args()

    microminer_input = args.microminer
    tmalign_input = args.tmalign
    outdir = Path(args.outdir)

    if not outdir.is_dir():
        print('Error: Specified output directory does not exist or is not a directory.')
        sys.exit(1)

    logging.basicConfig(filename=str((outdir / 'log.log').absolute()),
                        level=logging.INFO,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    logger.info(f'Starting scripts: {" ".join(sys.argv)}')

    microminer_input = [Path(_) for _ in microminer_input]
    microminer_files = collect_files(microminer_input, suffix='resultStatistic.csv')
    logger.info(f'Gathered {len(microminer_files)} input MicroMiner result files.')

    col_dtypes = {
        # PDB IDs are sometimes interpreted as a float given in scientific notation.
        # And residue position can contain iCode which are represented as string.
        constants.MICROMINER_QUERY_NAME: str,
        constants.MICROMINER_HIT_NAME: str,
        constants.MICROMINER_QUERY_POS: str,
        constants.MICROMINER_HIT_POS: str,
    }
    with helper.utils.timer('Reading MicroMiner result CSVs'):
        # collect all csv files in a single dataframe.
        df_mm = pd.concat([pd.read_csv(p, sep='\t', header=0, dtype=col_dtypes)
                           for p in microminer_files]).drop_duplicates()
    if df_mm.shape[0] == 0:
        print('Error: MicroMiner input files are empty.')
        sys.exit(1)
    logger.info(f'Gathered {df_mm.shape[0]} MicroMiner result rows.')

    df_cla = helper.datasets.scope.read_scope_classification_file()
    if df_cla.shape[0] == 0:
        print('Error: SCOPe classification file empty.')
        sys.exit(1)

    b4 = df_mm.shape[0]
    df_mm.query('queryName != hitName', inplace=True)
    logger.info(f'Removed self hits in MicroMiner. Remaining {df_mm.shape[0]} from {b4}')
    b4 = df_mm.shape[0]
    with helper.utils.timer('Merging with SCOPe classifications'):
        df = df_mm.merge(df_cla, left_on=['queryName'], right_on=['sid'], how='left',
                      suffixes=('1', '2'))
        df = df.merge(df_cla, left_on=['hitName'], right_on=['sid'], how='left',
                      suffixes=('1', '2'))
    assert b4 == df.shape[0]

    if tmalign_input is not None:
        df_tm = add_tmalign_results(tmalign_input)
        if df_tm.shape[0] > 0:
            with helper.utils.timer('Merging MicroMiner and TMAlign dataframes'):
                df = df_mm.merge(df_tm, left_on=['queryName', 'hitName'], right_on=['id1', 'id2'],
                                 how='left')
    assert b4 == df.shape[0]

    df.to_csv(outdir / f'scope_microminer_info.tsv', sep='\t', header=True, index=False)
    write_scope_stats(df, outdir)


if __name__ == "__main__":
    main()
