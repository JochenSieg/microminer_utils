"""
Annotate all mutations in the data set that have a wild-type PDB with PDB structures for the mutant.
"""
import os
from pathlib import Path
import argparse
import sys

import helper
from helper import BAD_PDBIDS, SUPPORTED_DATASETS
from helper.annotate_mutations import annotate_mutations, annotate_mutations2
from helper.datasets import read_mutation_dataset_single


def main():

    parser = argparse.ArgumentParser(
        description=
        """
        Executes the heuristic database search with MutScreen on a data set and annotates
        the datasets mutations with PDB structures for the mutant.
        """)

    parser.add_argument('--dataset', '-d', required=True, type=str.lower,
                        choices=SUPPORTED_DATASETS,
                        help='Data set name. Data location is inferred from config.ini')
    parser.add_argument('--outdir', '-o', default=os.getcwd(), type=str,
                        help='Path to output directory')
    parser.add_argument('--cpus', '-c', default=1, type=int,
                        help='Number of processes to use for parallel execution')

    # TODO in der annotate und screening Anwendung kann ich ein Cache-PATH list einrichten wo fuer jede
    #      wild-type PDB aus dem Standard mirror gecheckt wird, ob die schon berechnet wurde. Dann sollte
    #      die gesamte Benchmarking platform viiiiel schneller durch laufen (weil protherm, thermomutdb
    #      und prothermdb und sicherlich auch irgendwann die fireprotdb einen 80% ueberlap haben oder so)

    args = parser.parse_args()

    dataset_name = args.dataset
    outdir = Path(args.outdir)
    cpus = args.cpus

    if not outdir.is_dir():
        print('Error: Specified output directory does not exist or is not a directory.')
        sys.exit(1)

    # read in table of dataset
    df = read_mutation_dataset_single(dataset_name=dataset_name, pdb_mutant_only=False)

    # prepare outdir
    outdir = outdir / 'mutscreen_annotate_{}'.format(dataset_name)
    outdir.mkdir(parents=False, exist_ok=True)

    df = df[~df[helper.WILD_COL].isin(BAD_PDBIDS)]

    # TODO das mit den Chains ist noch ggf ein Problem. Aktuell nutze ich die nur bei skempi
    #      platinum, prothermdb sehen auch noch aus als haben sie valide chain ids
    use_chainid = dataset_name == 'skempi2'
    df_anno = annotate_mutations2(df, outdir, use_chainid=use_chainid, cpus=cpus, pdb_mirror=dataset_name)

    df_anno.to_csv(outdir / '{}_mutscreen_annotated.csv'.format(dataset_name), sep='\t', index=False)

    u_cols = [helper.WILD_COL, helper.WILD_AA, helper.WILD_SEQ_NUM, helper.MUT_AA]
    if use_chainid:
        u_cols.append(helper.WILD_CHAIN)
    _df = df_anno.drop_duplicates(subset=u_cols)
    print('{} of {} ({}%) unique mutations could be annotated.'.format(
        _df.query('wildName == wildName').shape[0], _df.shape[0],
        100*_df.query('wildName == wildName').shape[0] / _df.shape[0]))



if __name__ == "__main__":
    main()
