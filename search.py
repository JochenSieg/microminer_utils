"""
Runs MicroMiner on a list of input protein structures in search mode.
"""
import argparse
import logging
import os
import sys
from pathlib import Path

import pandas as pd

from helper.hpc import distribute_csv
from helper.runners import MicroMinerSearch

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="""
        Search with MicroMiner. 
        """
    )

    parser.add_argument(
        "--dataset", "-d", required=True, type=str, help="Path to dataset TSV file."
    )
    parser.add_argument(
        "--mode",
        "-m",
        required=False,
        type=str,
        default="single_mutation",
        choices=["standard", "single_mutation"],
        help="Search mode to run MicroMiner.",
    )
    parser.add_argument(
        "--representation",
        "-r",
        required=False,
        type=str,
        default="monomer",
        choices=["full_complex", "monomer", "ppi"],
        help="How input structures should be represented for MicroMiner search.",
    )
    parser.add_argument(
        "--outdir", "-o", default=os.getcwd(), type=str, help="Path to output directory"
    )
    parser.add_argument(
        "--cpus",
        "-c",
        default=1,
        type=int,
        help="Number of processes to use for parallel execution",
    )
    parser.add_argument(
        "--hpc",
        default=False,
        action="store_true",
        help="Run calculation on HPC (ZBH in-house cluster).",
    )

    args = parser.parse_args()

    dataset_file = Path(args.dataset)
    outdir = Path(args.outdir)
    cpus = args.cpus
    is_hpc = args.hpc
    mm_mode = args.mode
    mm_repr = args.representation

    if not dataset_file.is_file():
        print("Error: Dataset file does not exist.")
        sys.exit(1)
    if not outdir.is_dir():
        print("Error: Specified output directory does not exist or is not a directory.")
        sys.exit(1)

    logging.basicConfig(
        filename=str((outdir / "log.log").absolute()),
        level=logging.INFO,
        format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
    )
    logger.info(f'Starting scripts: {" ".join(sys.argv)}')

    # prepare outdir
    outdir.mkdir(parents=False, exist_ok=True)

    if is_hpc:
        distribute_csv(
            dataset_file,
            runner=MicroMinerSearch(
                cpus=1, mm_mode=mm_mode, mm_repr=mm_repr, raise_error=False
            ),
            outdir=outdir,
            job_name="search",
            cpus=cpus,
        )
    else:
        runner = MicroMinerSearch(
            cpus=cpus, mm_mode=mm_mode, mm_repr=mm_repr, raise_error=False
        )
        perf_dict_list = runner.run(dataset_file, outdir=outdir)
        # df_perf = pd.DataFrame(perf_dict_list).drop(['stdout', 'stderr', 'exit_code'], axis=1)
        df_perf = pd.DataFrame(perf_dict_list)
        df_perf.to_csv(outdir / "perf.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
