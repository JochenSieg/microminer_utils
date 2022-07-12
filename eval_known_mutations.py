"""
Compares mutation pairs in MicroMiner results with known mutations of a mutation dataset
and reports the number of correctly retrieved mutation structure pairs.
"""

import os
from typing import List

import pandas as pd
from pathlib import Path
import argparse
import sys
import logging

import helper
from helper import MUTATION_DATASETS_WITH_STRUCTURE_PAIRS, DATASETS_WITH_CUSTOM_STRUCTUREFILES, \
                   BAD_PDBIDS
from helper.datasets import read_mutation_dataset_single
from helper import constants

logger = logging.getLogger(__name__)


def generate_report(df_merged: pd.DataFrame, dataset_name: str, outdir: Path, backward: bool):
    pass


def run_mutation_checking(csv_input: List[Path], dataset_names: List[str], outdir: Path,
                          backward: bool):

    # gather all 'resultStatistic.csv' in the csv_input list (including recursive read of dirs)
    files = []
    for path in csv_input:
        files.extend([p for p in helper.utils.scantree(path) if p.name == 'resultStatistic.csv'])
    if len(files) == 0:
        print('Error: No resultStatistic.csv in input (and not in subdirs of any input dir).')
        sys.exit(1)
    logger.info(f'Gathered {len(files)} input resultStatistic.csv files.')

    # collect all csv files in a single dataframe.
    col_dtypes = {
        # PDB IDs are sometimes interpreted as a float given in scientific notation.
        constants.MICROMINER_QUERY_NAME: str,
        constants.MICROMINER_HIT_NAME: str,
        constants.MICROMINER_QUERY_POS: str,
        constants.MICROMINER_HIT_POS: str,
    }
    df_res = pd.concat([pd.read_csv(p, sep='\t', header=0, dtype=col_dtypes)
                        for p in files]).drop_duplicates()
    if df_res.shape[0] == 0:
        print('Error: resultStatistic.csv input files are empty.')
        sys.exit(1)
    logger.info(f'Gathered {df_res.shape[0]} input structure pairs.')

    right_on = [constants.MICROMINER_QUERY_NAME, constants.MICROMINER_HIT_NAME,
                constants.MICROMINER_QUERY_AA, constants.MICROMINER_HIT_AA,
                constants.MICROMINER_QUERY_POS]
    merged_file_suffix = '_eval.tsv'
    report_file = outdir / 'eval_report.txt'
    if backward:
        # right_on = [constants.MICROMINER_HIT_NAME, constants.MICROMINER_QUERY_NAME,
        #             constants.MICROMINER_HIT_AA, constants.MICROMINER_QUERY_AA,
        #             constants.MICROMINER_HIT_POS]
        right_on = [constants.MICROMINER_QUERY_NAME, constants.MICROMINER_HIT_NAME,
                    constants.MICROMINER_QUERY_AA, constants.MICROMINER_HIT_AA,
                    constants.MICROMINER_HIT_POS]
        merged_file_suffix = '_eval_backward.tsv'
        report_file = outdir / 'eval_report_backwards.txt'

    # Unfortunately, most mutation dataset hold no or inconsistent chain ID annotations.
    # For this reason we ignore the chain ID for matching retrieved mutations and
    # we have to drop duplicates in this sense, i.e. matches from homo-mers.
    df_res.drop_duplicates(right_on, inplace=True)

    for dataset_name in dataset_names:

        logger.info(f'Evaluating known mutations of {dataset_name}')
        df_ref = read_mutation_dataset_single(dataset_name=dataset_name, pdb_mutant_only=True)

        df_ref = df_ref[~df_ref[constants.WILD_COL].isin(BAD_PDBIDS)]
        df_ref = df_ref[~df_ref[constants.MUTANT_COL].isin(BAD_PDBIDS)]

        # drop duplicates. There are often duplicates, e.g. because of multiple ddG measurments
        df_ref.drop_duplicates(
            [helper.WILD_COL, helper.WILD_AA, helper.WILD_SEQ_NUM, helper.MUT_AA, helper.MUTANT_COL],
            inplace=True)

        from_one_to_three = lambda aa: constants.one_2_three_dict[aa] if len(aa) != 3 else aa
        df_ref['wild_aa3'] = df_ref[constants.WILD_AA].apply(from_one_to_three)
        df_ref['mutant_aa3'] = df_ref[constants.MUT_AA].apply(from_one_to_three)

        left_on = [constants.WILD_COL, constants.MUTANT_COL,
                   'wild_aa3', 'mutant_aa3', constants.WILD_SEQ_NUM]
        if backward:
            left_on = [constants.MUTANT_COL, constants.WILD_COL, 'mutant_aa3', 'wild_aa3',
                       constants.WILD_SEQ_NUM]

        df_merged = df_ref.merge(df_res,
                                 left_on=left_on,
                                 right_on=right_on,
                                 how='left')
        df_anno = df_merged.dropna(subset=right_on)
        logger.info(f'MicroMiner found {df_anno.shape[0]} of {df_merged.shape[0]} mutations in'
                    f' {dataset_name}' + (' (backward)' if backward else ''))

        df_not_found = df_ref.merge(df_res,
                                    left_on=left_on,
                                    right_on=right_on,
                                    how='left', indicator=True)
        df_not_found.query('_merge == "left_only"', inplace=True)
        df_not_found.drop('_merge', axis=1, inplace=True)

        df_merged.to_csv(outdir / f'{dataset_name}{merged_file_suffix}', sep='\t', header=True,
                         index=False)
        df_not_found.to_csv(outdir / f'{dataset_name}_not_found{merged_file_suffix}', sep='\t',
                            header=True, index=False)
        with open(report_file, 'a+') as f:
            f.write(f'{dataset_name}{merged_file_suffix}'
                    f'\t{df_anno.shape[0]}'
                    f'\t{df_merged.shape[0]}'
                    f'\n')

        # generate_report(df_merged, dataset_name, outdir, backward)


def main():

    parser = argparse.ArgumentParser(
        description=
        """
        Compares mutation pairs in MicroMiner results with known mutations of a mutation dataset. 
        """)

    parser.add_argument('--csv', '-c', required=True, type=str, nargs='+',
                        help='resultStatistic.csv file or dir in which resultStatistic.csv files'
                             ' will be searched recursively.')
    parser.add_argument('--dataset', '-d', required=True, type=str.lower,
                        choices=MUTATION_DATASETS_WITH_STRUCTURE_PAIRS, nargs='+',
                        help='Data set name. Data location is inferred from config.ini. Will be'
                             ' used as ground truth of mutation structure pairs.')
    parser.add_argument('--outdir', '-o', default=os.getcwd(), type=str,
                        help='Path to output directory')
    parser.add_argument('--backward', '-b', default=False, action='store_true',
                        help='Invert wild-type and mutant')

    args = parser.parse_args()

    csv_input = args.csv
    dataset_names = args.dataset
    outdir = Path(args.outdir)
    backward = args.backward

    if not outdir.is_dir():
        print('Error: Specified output directory does not exist or is not a directory.')
        sys.exit(1)

    logging.basicConfig(filename=str((outdir / 'log.log').absolute()),
                        level=logging.INFO,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    logger.info(f'Starting scripts: {" ".join(sys.argv)}')

    csv_input = [Path(_) for _ in csv_input]
    run_mutation_checking(csv_input, dataset_names, outdir, backward)


if __name__ == "__main__":
    main()
