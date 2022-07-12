"""
Executes the SIENA version of the tool. This is the database search version of the tool.
All sites in a query protein are aligned heuristically to the structures from a database.
"""
import os
from pathlib import Path
import argparse
import sys

import helper
from helper import BAD_PDBIDS, SUPPORTED_DATASETS
from helper.eval import eval_screening_tool, eval_screening_tool2
from helper.datasets import read_mutation_dataset_single


def main():

    parser = argparse.ArgumentParser(
        description=
        """
        Executes the heuristic database search with MutScreen on a data set.
        """)

    parser.add_argument('--dataset', '-d', required=True, type=str.lower, choices=SUPPORTED_DATASETS,
                        help='Data set name. Data location is inferred from config.ini')
    parser.add_argument('--outdir', '-o', default=os.getcwd(), type=str,
                        help='Path to output directory')
    parser.add_argument('--backward', '-b', default=False, action='store_true',
                        help='Search the wild type with mutant')
    parser.add_argument('--cpus', '-c', default=1, type=int,
                        help='Number of processes to use for parallel execution')
    parser.add_argument('--cache_dirs', '-e', default=[], type=list, nargs='+',
                        help='List of directories for reading cached results.')

    args = parser.parse_args()

    dataset_name = args.dataset
    outdir = Path(args.outdir)
    backward = args.backward
    cpus = args.cpus

    if not outdir.is_dir():
        print('Error: Specified output directory does not exist or is not a directory.')
        sys.exit(1)

    # read in table of dataset
    df = read_mutation_dataset_single(dataset_name=dataset_name, pdb_mutant_only=True)

    # for the pair experiment we are only interested in the unique pairs of WILD PDB ID, WILD AA, SEQ POS, MUTANT AA,
    # MUTANT PDB ID and potentially the CHAIN ID of the mutation in the wild-type if the data set provides this.
    # Datasets with chain ids:
    #      ProTherm: has almost no chain ids -> can't do anything there
    #      ThermoMutDB: has almost only A chain ids (~700) and 7 chain ids called "unsigned" -> problematic
    #      Platinum: seems to have valid chain ids
    df.drop_duplicates([helper.WILD_COL, helper.WILD_AA, helper.WILD_SEQ_NUM, helper.MUT_AA, helper.MUTANT_COL],
                       inplace=True)

    # prepare outdir
    outdir = outdir / 'screen_experiment_{}'.format(dataset_name)
    outdir.mkdir(parents=False, exist_ok=True)

    if args.backward:
        outdir = outdir / 'backward'
    else:
        outdir = outdir / 'forward'
    outdir.mkdir(parents=False, exist_ok=False)

    df = df[~df[helper.WILD_COL].isin(BAD_PDBIDS)]
    df = df[~df[helper.MUTANT_COL].isin(BAD_PDBIDS)]

    datasets_with_custom_pdb = ['platinum', 'skempi2']
    pdb_mirror = 'standard' if backward and dataset_name in datasets_with_custom_pdb else dataset_name

    # call the tool and annotate which mutations were found.
    df['mutscreen_screen_found'] = eval_screening_tool2(df, outdir,
                                                        cpus=cpus, backward=backward, pdb_mirror=pdb_mirror)

    df.to_csv(outdir / '{}_mutscreen_screening.csv'.format(dataset_name), sep='\t', index=False)

    s = df['mutscreen_screen_found'].sum()
    print('RESULT: {}_mutscreen_screen_experiment: {}/{} ( {} ) mutations found'.format(
        dataset_name, s, df.shape[0], 100 * s / df.shape[0]))
    # TODO wie schreibe ich ne gescheite Statistik raus?


if __name__ == "__main__":
    main()
