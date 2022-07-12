import os
from pathlib import Path
import argparse
import sys
import logging

import helper
from helper import DATASETS_WITH_CUSTOM_STRUCTUREFILES, MUTATION_DATASETS_WITH_STRUCTURE_PAIRS, \
    BAD_PDBIDS
from helper.datasets import read_mutation_dataset_single

logger = logging.getLogger(__name__)

PATH_COL1 = 'structure_path1'
PATH_COL2 = 'structure_path2'
ID_COL1 = 'id1'
ID_COL2 = 'id2'


def make_from_predefined_dataset(dataset_names, backward, outdir):
    for dataset_name in dataset_names:

        logger.info(f'Creating pair datasets for: {dataset_name}')

        pdb_mirror1 = 'standard'
        pdb_mirror2 = 'standard'
        if dataset_name in DATASETS_WITH_CUSTOM_STRUCTUREFILES:
            pdb_mirror1 = dataset_name

        id_col1 = helper.WILD_COL
        id_col2 = helper.MUTANT_COL
        outfile_name = f'{dataset_name}_pair.tsv'
        df = read_mutation_dataset_single(dataset_name=dataset_name, pdb_mutant_only=True)

        df = df[~df[helper.WILD_COL].isin(BAD_PDBIDS)]
        df = df[~df[helper.MUTANT_COL].isin(BAD_PDBIDS)]

        if backward:
            if helper.MUTANT_COL not in df.columns:
                print(f'ERROR: backward mode war specified, but selected dataset has no mutant'
                      f' structure annotation: {dataset_name}')
                sys.exit(1)
            # backward mode always with standard PDB files.
            pdb_mirror2, pdb_mirror1 = pdb_mirror1, pdb_mirror2
            outfile_name = f'{dataset_name}_pair_backward.tsv'
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
        df_out = df_out.rename(columns={id_col1: ID_COL1, id_col2: ID_COL2})
        df_out.to_csv(outdir / outfile_name, sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser(
        description=
        """
        Creates a "dataset" CSV file of protein mutation structure pairs from pre-defined datasets. 
        """)

    parser.add_argument('--dataset', '-d', required=True, type=str.lower,
                        choices=MUTATION_DATASETS_WITH_STRUCTURE_PAIRS, nargs='+',
                        help='Data set name. Data location is inferred from config.ini')
    parser.add_argument('--outdir', '-o', default=os.getcwd(), type=str,
                        help='Path to output directory')
    parser.add_argument('--backward', '-b', default=False, action='store_true',
                        help='Invert wild-type and mutant')

    args = parser.parse_args()

    dataset_names = args.dataset
    outdir = Path(args.outdir)
    backward = args.backward

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    logger.info(f'Starting scripts: {" ".join(sys.argv)}')

    if not outdir.is_dir():
        print('Error: Specified output directory does not exist or is not a directory.')
        sys.exit(1)

    make_from_predefined_dataset(dataset_names, backward, outdir)


if __name__ == "__main__":
    main()
