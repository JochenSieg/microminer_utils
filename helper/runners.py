import logging
import multiprocessing
from pathlib import Path
from typing import List, Dict

import pandas as pd

from .cmdl_calls import call_microminer_search, call_microminer_pair
from .utils import parse_microminer_search_stdout

logger = logging.getLogger(__name__)


def _run_parallel(func, parameter_set, cpus: int):
    """Helper function to run a function with sets of parameters in parallel.

    :param func: Compute function.
    :param parameter_set: Parameters.
    :param cpus: Number of CPU cores.
    :return:
    """

    # parallel execution with each parameter set separately
    with multiprocessing.Pool(cpus) as pool:
        return pool.starmap(
            func,
            parameter_set,
        )


class MicroMinerSearch:
    """Manages execution of MicroMiner searches in parallel."""

    MANDATORY_TSV_COLUMNS = ["id", "structure_path"]

    def __init__(
        self, mm_mode: str, mm_repr: str, cpus: int = 1, raise_error: bool = True
    ):
        """Create a new runner.

        :param mm_mode: The search mode of MicroMiner.
        :param mm_repr: The structure represention mode for MicroMiner.
        :param cpus: Number CPU cores to use.
        :param raise_error: Whether to raise an error when a MicroMiner call fails.
        """
        self.cpus = cpus
        self.raise_error = raise_error
        self.mm_mode = mm_mode
        self.mm_repr = mm_repr

    def run(self, param_tsv: Path, outdir: Path) -> List[Dict]:
        """Run MicroMiner searches.

        :param param_tsv: The parameter file with input to MicroMiner.
        :param outdir: Directory for writing results.
        :return: List of detailed MicroMiner output (return code, standard out, standard error)
        """
        df = pd.read_csv(param_tsv, sep="\t", header=0)
        logger.info(f"Read {df.shape[0]} parameter records for computation")

        if not all(c in df.columns for c in MicroMinerSearch.MANDATORY_TSV_COLUMNS):
            raise ValueError(f"Missing mandatory fields in TSV: {param_tsv}")

        # generate a list of parameter tuple containing the PDB pairs to be aligned.
        parameter_set = [
            (
                Path(getattr(row, MicroMinerSearch.MANDATORY_TSV_COLUMNS[1])),
                outdir
                / "{}".format(getattr(row, MicroMinerSearch.MANDATORY_TSV_COLUMNS[0])),
                self.mm_mode,
                self.mm_repr,
                self.raise_error,
            )
            for row in df.drop_duplicates().itertuples(index=False)
        ]

        out_list = []
        if self.cpus > 1:
            out_list = _run_parallel(call_microminer_search, parameter_set, self.cpus)
        else:
            for param_set in parameter_set:
                out = call_microminer_search(*param_set)
                out_list.append(out)
        out_parsed_list = []
        for out in out_list:
            parsed_stdout = parse_microminer_search_stdout(
                out["stdout"].decode("utf-8")
            )
            out_parsed_list.append(parsed_stdout)
        return out_parsed_list


class MicroMinerPair:
    """Manages execution of MicroMiner pair alignment in parallel."""

    MANDATORY_TSV_COLUMNS = ["id1", "structure_path1", "id2", "structure_path2"]

    def __init__(self, cpus: int = 1, raise_error: bool = True):
        """Construct a new runner.

        :param cpus: Number CPU cores to use.
        :param raise_error: Whether to raise an error when a MicroMiner call fails.
        """
        self.cpus = cpus
        self.raise_error = raise_error

    def run(self, param_tsv: Path, outdir: Path):
        """Run MicroMiner pair alignment.

        :param param_tsv: The parameter file with input to MicroMiner.
        :param outdir: Directory for writing results.
        :return: None
        """
        df = pd.read_csv(param_tsv, sep="\t", header=0)
        logger.info(f"Read {df.shape[0]} parameter records for computation")

        if not all(c in df.columns for c in MicroMinerPair.MANDATORY_TSV_COLUMNS):
            raise ValueError(f"Missing mandatory fields in TSV: {param_tsv}")

        # generate a list of parameter tuple containing the PDB pairs to be aligned.
        parameter_set = [
            (
                Path(getattr(row, MicroMinerPair.MANDATORY_TSV_COLUMNS[1])),
                Path(getattr(row, MicroMinerPair.MANDATORY_TSV_COLUMNS[3])),
                outdir
                / "{}_{}".format(
                    getattr(row, MicroMinerPair.MANDATORY_TSV_COLUMNS[0]),
                    getattr(row, MicroMinerPair.MANDATORY_TSV_COLUMNS[2]),
                ),
                self.raise_error,
            )
            for row in df.drop_duplicates().itertuples(index=False)
        ]

        if self.cpus > 1:
            _run_parallel(call_microminer_pair, parameter_set, self.cpus)
        else:
            for param_set in parameter_set:
                call_microminer_pair(*param_set)
