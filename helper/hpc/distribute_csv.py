import getpass
import json
import logging
import pickle
import subprocess
import tempfile
import time
from collections import namedtuple
from pathlib import Path

import pandas as pd

from .check import sanity_check_microminer_result_dir
from .sge import SGEJobRunner
from helper.constants import CONFIG
from helper.runners import MicroMinerPair, MicroMinerSearch

logger = logging.getLogger(__name__)

SshSource = namedtuple("SshSource", "user, path, host, port")

ALL_QUEUES = json.loads(CONFIG.get("HPC", "QUEUES"))


def get_chunksize(runner, cpus, input_rows):
    def default_chunksize():
        """
        This is more or less the same heuristic that the
        multiprocessing module uses.
        :return: chunksize
        """
        chunksize, extra = divmod(input_rows, cpus * 4)
        if extra:
            chunksize += 1
        return chunksize

    if type(runner) == MicroMinerPair:
        return min(50, input_rows)
    return default_chunksize()


def generate_microminer_preparation_script():
    """
    Generates a Bash script that prepares necessary files for executing
    MicroMiner (mainly the k-mer look-up index).
    :return: The Bash script as string.
    """
    kmer_index = CONFIG["DATABASES"]["SITE_SEARCH_DB"]
    user = getpass.getuser()
    return f"""#! /bin/bash
    
# This script should be executed from the job script at the HPC node.

# create a cache dir for the necessary files on the local node.
cache_dir_name="microminer_distributed_cache"
node_index_dir_suffix="{user}/${{cache_dir_name}}"
node_index_dir="/local/${{node_index_dir_suffix}}"
# if [ -d "/ssd_local" ]; then
#     node_index_dir="/ssd_local/${{node_index_dir_suffix}}"
# elif [ -d "/local" ]; then
#     node_index_dir="/local/${{node_index_dir_suffix}}"
# else
#     echo "Node has no expected drive! Abort"
#     exit 1
# fi
mkdir -p ${{node_index_dir}}

# copy the index (might take a while if not cached)
# rsync -ra {kmer_index}* ${{node_index_dir}}

# create an environment variable for the helper config to specify this nodes
#   location of the lookup index.
# NOTE: that this script must be sourced by the calling script.
export SITE_SEARCH_DB="${{node_index_dir}}/$(basename {kmer_index})"
"""


def generate_runner_script(
    runner,
):
    """
    Generates a simple Python script which executes a serialized version
    of the provided runner instance for an arbitrary (correctly formatted) input file
    and writes the runners results to a result dir.

    This script is intended as the interface to this helper module on the HPC cluster.
    :param runner: A runner instance.
    :return: The Python script as string.
    """
    return f"""
# mini python runner script
from helper.runners import {type(runner).__name__}
import sys
import pickle
from pathlib import Path
import logging
assert len(sys.argv) > 2, 'Missing command line arguments!'
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    handlers=[logging.StreamHandler(sys.stdout)])
runner = pickle.loads({pickle.dumps(runner)})
runner.run(param_tsv=Path(str(sys.argv[1])), outdir=Path(str(sys.argv[2])))
"""


def generate_hpc_script(
    job_name: str,
    global_working_dir: Path,
    local_working_dir: str,
    input_file: Path,
    stdout_file: Path,
    stderr_file: Path,
    nof_jobs: int,
    cpus: int,
    conda_bin_path: Path,
    conda_env_name: str,
    queues: list,
    add_to_pythonpath: list,
    runner_script_path: Path,
    prepare_script_path: Path,
    chunksize: int,
    copy_ssh: list = [],
):
    newline = "\n"
    return f"""#! /bin/bash
#$ -N {job_name}
#$ -wd {global_working_dir.resolve()}
#$ -q {','.join(queues)}
#$ -o {stdout_file.resolve()}
#$ -e {stderr_file.resolve()}
#$ -t 1-{nof_jobs}
#$ -tc {cpus}
# -S /bin/bash
# -l os=42.2
# -l mem_free=1G
# exclude nodes which disc storage is too small to hold the PDB locally. TODO: they might change
# we need to use PDB locally because the number of reads per sec is too much for the NFS at /data
# effectively killing the HPC nodes. 
#$ -l h='!node113'


# avoid core dumps. 
ulimit -c 0

GLOBAL_WORK_DIR="{global_working_dir.resolve()}"
LOCAL_WORK_DIR="{local_working_dir}"
INPUT_FILE="{input_file.resolve()}"

echo "LOG: Task ${{SGE_TASK_ID}} running on $(hostname)"

mkdir -p ${{LOCAL_WORK_DIR}}
if [ $? -ne 0 ]; then
  echo "$0: Can't create local working dir!"
  exit 1
fi
THIS_TMPDIR=$(mktemp -d --tmpdir=${{LOCAL_WORK_DIR}})
if [ $? -ne 0 ]; then
  echo "$0: Can't create local temp dir!"
  exit 1
fi

# source the conda env
export PATH="{str(conda_bin_path)}:${{PATH}}"
source activate {conda_env_name}

# use this_tmpdir as a source for python packages/scripts
export PYTHONPATH="${{THIS_TMPDIR}}:${{PYTHONPATH}}"
{f'export PYTHONPATH="' + ':'.join([f'{a}' for a in add_to_pythonpath]) + f':${{PYTHONPATH}}"'}

{'# copy further files' if len(copy_ssh) > 0 else ''}
{f'{newline}'.join([f'scp -r -P {s.port} {s.user}@{s.host}:{s.path} ${{THIS_TMPDIR}}'
                    for s in copy_ssh])}

# run/source preparation script
source {prepare_script_path.resolve()}

THIS_INPUT_FILE="${{THIS_TMPDIR}}/input.tsv"
startline=$(( ($SGE_TASK_ID - 1) * {chunksize} ))
head -n1 ${{INPUT_FILE}} > ${{THIS_INPUT_FILE}}
tail -n "+$((startline+2))" ${{INPUT_FILE}} | head -n {chunksize} >> ${{THIS_INPUT_FILE}}

THIS_RESULTS_DIR="${{THIS_TMPDIR}}/results"
mkdir ${{THIS_RESULTS_DIR}}

echo "Calling: python {runner_script_path.resolve()} ${{THIS_INPUT_FILE}} ${{THIS_TMPDIR}}"
python {runner_script_path.resolve()} "${{THIS_INPUT_FILE}}" "${{THIS_RESULTS_DIR}}"

# collect results back in global working dir
rsync -ra "${{THIS_RESULTS_DIR}}/" "${{GLOBAL_WORK_DIR}}/results"

# clean up
rm -r ${{THIS_TMPDIR}}

echo "FINISHED das Script"
"""


def distribute_csv(
    dataset_file: Path, runner, outdir: Path, job_name: str, cpus: int = 1
) -> None:
    """Distributes the computation across the SGE cluster. Each row in the input CSV
    corresponds to a single computation. This function splits the input rows
    in chunks and submits the computation to the SGE cluster.

    :param dataset_file: The input dataset file (the input data).
    :param runner: The runner class (defines what will be computed)
    :param outdir: Directory to write results.
    :param job_name: A name for the SGE job.
    :param cpus: Number of cores to use.
    :return: None
    """
    df = pd.read_csv(dataset_file, sep="\t", header=0)

    if df.shape[0] == 0:
        raise ValueError("Input data empty")

    tmpdir_name = time.strftime("%Y%m%d_%H%M%S")
    cpus = max(1, cpus)

    chunksize = get_chunksize(runner, cpus, df.shape[0])
    nof_jobs, extra = divmod(df.shape[0], chunksize)
    if extra:
        nof_jobs += 1
    assert nof_jobs <= df.shape[0]
    assert nof_jobs * chunksize >= df.shape[0]

    with tempfile.TemporaryDirectory(
        dir=CONFIG["HPC"]["HPC_WORKING_DIR"], prefix=tmpdir_name
    ) as t:
        tmpdir = Path(t)

        # write input parameter file to disc
        input_tsv_path = tmpdir / "input.tsv"
        df.to_csv(input_tsv_path, sep="\t", header=True, index=False)

        (tmpdir / "cluster_out").mkdir()

        # preparation script to set up environment (copy additional data etc.)
        prepare_script_str = ""
        if type(runner) == MicroMinerSearch:
            prepare_script_str = generate_microminer_preparation_script()
        prepare_script_path = tmpdir / "prepare.sh"
        with open(prepare_script_path, "w") as f:
            f.write(prepare_script_str)

        # python script that runs the inner calculation using this module
        runner_script_str = generate_runner_script(runner)
        # print(runner_script_str)
        runner_script_path = tmpdir / "runner_script.py"
        with open(runner_script_path, "w") as f:
            f.write(runner_script_str)

        (tmpdir / "results").mkdir()

        # the SGE cluster job script
        job_script_str = generate_hpc_script(
            job_name=job_name,
            global_working_dir=tmpdir,
            local_working_dir=CONFIG["HPC"]["HPC_LOCAL_WORKING_DIR"],
            input_file=tmpdir / "input.tsv",
            stdout_file=tmpdir / "cluster_out",
            stderr_file=tmpdir / "cluster_out",
            nof_jobs=nof_jobs,
            cpus=cpus,
            queues=ALL_QUEUES,
            conda_bin_path=Path(CONFIG["HPC"]["CONDA_BIN_PATH"]),
            conda_env_name=CONFIG["HPC"]["CONDA_ENV_NAME"],
            # copy_ssh=[this_py_module],
            add_to_pythonpath=[Path(CONFIG["HPC"]["PYPATH_PATHS"])],
            runner_script_path=runner_script_path,
            prepare_script_path=prepare_script_path,
            chunksize=chunksize,
        )
        # print(job_script_str)
        with open(tmpdir / "job_script", "w") as f:
            f.write(job_script_str)

        SGEJobRunner.submit_and_wait(tmpdir / "job_script")

        # rsync the results to the outdir
        cmd_call = [
            "rsync",
            "-ra",
            "--info=progress2",
            (tmpdir / "results").resolve(),
            outdir.resolve(),
        ]
        res = subprocess.run(cmd_call, stdout=subprocess.PIPE)
        if res.returncode != 0:
            logger.error(res)
            raise ValueError("Encountered problem {}".format(res))

        # rsync stdout/err of HPC jobs to a log dir
        cmd_call = [
            "rsync",
            "-ra",
            "--info=progress2",
            (tmpdir / "cluster_out").resolve(),
            outdir.resolve(),
        ]
        res = subprocess.run(cmd_call, stdout=subprocess.PIPE)
        if res.returncode != 0:
            logger.error(res)
            raise ValueError("Encountered problem {}".format(res))

        if type(runner) == MicroMinerSearch or type(runner) == MicroMinerPair:
            is_sane = sanity_check_microminer_result_dir(outdir / "results", df)
