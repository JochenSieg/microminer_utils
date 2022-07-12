"""
Executes the ASCONA version of the tool on the data set configured in config.ini
All sites in a query protein are aligned to the whole structure of a second
protein. No data base search only pairwise analysis.

"""
import os
from pathlib import Path
import argparse
import sys
import logging

import helper
from helper import BAD_PDBIDS
from helper.eval import eval_pair_tool, eval_pair_tool2
from helper.datasets import read_mutation_dataset_single

logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description=
        """
        Executes the pair alignment with MutScreen on a data set.
        """)

    parser.add_argument('--dataset', '-d', required=True, type=str.lower,
                        choices=['protherm', 'thermomutdb', 'platinum'],
                        help='Data set name. Data location is inferred from config.ini')
    parser.add_argument('--outdir', '-o', default=os.getcwd(), type=str,
                        help='Path to output directory')
    parser.add_argument('--backward', '-b', default=False, action='store_true',
                        help='Search the wild type with mutant')
    parser.add_argument('--cpus', '-c', default=1, type=int,
                        help='Number of processes to use for parallel execution')

    args = parser.parse_args()

    dataset_name = args.dataset
    outdir = Path(args.outdir)
    backward = args.backward
    cpus = args.cpus

    if not outdir.is_dir():
        print('Error: Specified output directory does not exist or is not a directory.')
        sys.exit(1)

    # prepare outdir
    outdir = outdir / 'pair_experiment_{}'.format(dataset_name)
    outdir.mkdir(parents=False, exist_ok=True)
    #
    # # Exchange query and target row if backward experiment.
    if args.backward:
        outdir = outdir / 'backward'
    else:
        outdir = outdir / 'forward'
    outdir.mkdir(parents=False, exist_ok=False)

    logging.basicConfig(filename=str((outdir / 'log.log').absolute()),
                        level=logging.INFO,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    logger.info(f'Starting scripts: {" ".join(sys.argv)}')

    df = read_mutation_dataset_single(dataset_name=dataset_name, pdb_mutant_only=True)

    # for the pair experiment we are only interested in the unique pairs of WILD PDB ID, WILD AA, SEQ POS, MUTANT AA,
    # MUTANT PDB ID and potentially the CHAIN ID of the mutation in the wild-type if the data set provides this.
    # Datasets with chain ids:
    #      ProTherm: has almost no chain ids -> can't do anything there
    #      ThermoMutDB: has almost only A chain ids (~700) and 7 chain ids called "unsigned" -> problematic
    #      Platinum: seems to have valid chain ids
    df.drop_duplicates([helper.WILD_COL, helper.WILD_AA, helper.WILD_SEQ_NUM, helper.MUT_AA, helper.MUTANT_COL],
                       inplace=True)

    df = df[~df[helper.WILD_COL].isin(BAD_PDBIDS)]
    df = df[~df[helper.MUTANT_COL].isin(BAD_PDBIDS)]

    # avoid PDB file look up in platinums prepared PDB files dictionary when running in backward mode
    datasets_with_custom_pdb = ['platinum', 'skempi2']
    pdb_mirror = 'standard' if backward and dataset_name in datasets_with_custom_pdb else dataset_name

    # call the tool and annotate which mutations were found.
    df['mutscreen_pair_found'] = eval_pair_tool2(df, outdir, cpus=cpus, backward=backward, pdb_mirror=pdb_mirror)

    #TODO TMalign und problematic annotieren?

    df.to_csv(outdir / '{}_mutscreen_pair.csv'.format(dataset_name), sep='\t', index=False)

    s = df['mutscreen_pair_found'].sum()
    logger.info('RESULT: {}_mutscreen_pair_experiment: {}/{} ( {} ) mutations found'.format(
        dataset_name, s, df.shape[0], 100 * s / df.shape[0]))
    # # TODO wie schreibe ich ne gescheite Statistik raus?


if __name__ == "__main__":
    main()
