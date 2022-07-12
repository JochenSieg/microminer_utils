import os
from pathlib import Path
import argparse
import sys
import logging

import helper
from helper import BAD_PDBIDS, SUPPORTED_DATASETS, MUTATION_DATASETS, \
    DATASETS_WITH_CUSTOM_STRUCTUREFILES
from helper.datasets import read_mutation_dataset_single
from helper.datasets.scope import read_scope
from helper.datasets.pisces import read_pisces

logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description=
        """
        Creates a "dataset" CSV file for a pre-defined data set for querying with MicroMiner.
        A "dataset" CSV saves protein IDs and path to structure files (i.e. PDB files).
        """)

    parser.add_argument('--dataset', '-d', required=False, type=str.lower,
                        choices=SUPPORTED_DATASETS, nargs='+', default=SUPPORTED_DATASETS,
                        help='Data set name. Data location is inferred from config.ini')
    parser.add_argument('--outdir', '-o', default=os.getcwd(), type=str,
                        help='Path to output directory')
    parser.add_argument('--backward', '-b', default=False, action='store_true',
                        help='Search the wild type with mutant (for mutation datasets only)')

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

    for dataset_name in dataset_names:

        logger.info(f'Creating datasets for: {dataset_name}')

        pdb_mirror = 'standard'
        if dataset_name in DATASETS_WITH_CUSTOM_STRUCTUREFILES:
            pdb_mirror = dataset_name

        df = None
        id_col = None
        path_col = 'structure_path'
        outfile_name = f'{dataset_name}.tsv'
        if dataset_name in MUTATION_DATASETS:
            df = read_mutation_dataset_single(dataset_name=dataset_name, pdb_mutant_only=True)

            # drop PDB Ids that are Calpha only or obsolete. We do not need those for mutation
            # evaluation experiments
            df = df[~df[helper.WILD_COL].isin(BAD_PDBIDS)]

            if backward:
                if helper.MUTANT_COL not in df.columns:
                    print(f'ERROR: backward mode war specified, but selected dataset has no mutant'
                          f' structure annotation: {dataset_name}')
                    sys.exit(1)
                pdb_mirror = 'standard'  # backward mode always with standard PDB files.
                df = df[~df[helper.MUTANT_COL].isin(BAD_PDBIDS)]

                df[path_col] = df[helper.MUTANT_COL].apply(
                    lambda val: helper.utils.get_pdb_file_path(val, mirror=pdb_mirror,
                                                               allow_obsolete=True,
                                                               ignore_missing=True))
                id_col = helper.MUTANT_COL
                outfile_name = f'{dataset_name}_backward.tsv'
            else:
                df[path_col] = df[helper.WILD_COL].apply(
                    lambda val: helper.utils.get_pdb_file_path(val, mirror=pdb_mirror,
                                                               allow_obsolete=True,
                                                               ignore_missing=True))
                id_col = helper.WILD_COL
        else:
            if dataset_name == 'scope':
                df = read_scope()
                id_col = 'sid'
            elif dataset_name == 'pisces':
                df = read_pisces()
                id_col = 'pdb'
            df[path_col] = df[id_col].apply(
                lambda val: helper.utils.get_pdb_file_path(val, mirror=pdb_mirror,
                                                           allow_obsolete=True,
                                                           ignore_missing=True))

        df_out = df[[id_col, path_col]]

        # remove duplicate ROWS, i.e. when ID and PATH are identical. Would lead to redundant
        # computations.
        df_out = df_out.drop_duplicates()

        if df_out.isna().sum(axis=1).sum() > 0:
            logger.info(f'Dropping {df_out.isna().sum(axis=1).sum()} missing values'
                        f' (of {df_out.shape[0]} total)')
            df_out = df_out.dropna()  # remove missing values, e.g. when only CIF file is available

        df_out = df_out.rename(columns={id_col: 'id'})

        df_out.to_csv(outdir / outfile_name, sep='\t', index=False)


if __name__ == "__main__":
    main()
