"""
Runs TMAlign (Zhang group) on a list of protein structure pairs.
"""

import os
from pathlib import Path
import argparse
import sys
import logging

from helper.hpc import distribute_csv
from helper.runners import TMalign

logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description=
        """
        Structure Pair comparison with TMAlign. 
        """)

    parser.add_argument('--dataset', '-d', required=True, type=str,
                        help='Path to dataset TSV file.')
    parser.add_argument('--outdir', '-o', default=os.getcwd(), type=str,
                        help='Path to output directory')
    parser.add_argument('--cpus', '-c', default=1, type=int,
                        help='Number of processes to use for parallel execution')
    # parser.add_argument('--hpc', default=False, action='store_true',
    #                     help='Run calculation on HPC (ZBH in-house cluster).')
    parser.add_argument('--skip_same', default=False, action='store_true',
                        help='Skip same pairs of same ID and same filepath.')

    args = parser.parse_args()

    dataset_file = Path(args.dataset)
    outdir = Path(args.outdir)
    cpus = args.cpus
    # is_hpc = args.hpc
    skip_same = args.skip_same

    if not dataset_file.is_file():
        print('Error: Dataset file does not exist.')
        sys.exit(1)
    if not outdir.is_dir():
        print('Error: Specified output directory does not exist or is not a directory.')
        sys.exit(1)

    logging.basicConfig(filename=str((outdir / 'log.log').absolute()),
                        level=logging.INFO,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    logger.info(f'Starting scripts: {" ".join(sys.argv)}')

    # prepare outdir
    outdir.mkdir(parents=False, exist_ok=True)

    is_hpc = False
    if is_hpc:
        distribute_csv(dataset_file,
                       runner=TMalign(cpus=1, skip_same=skip_same, raise_error=False),
                       outdir=outdir,
                       job_name='tmalign',
                       cpus=cpus)
    else:
        runner = TMalign(cpus=cpus, skip_same=skip_same, raise_error=False)
        runner.run(dataset_file, outdir=outdir)


if __name__ == "__main__":
    main()
