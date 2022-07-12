import itertools
import multiprocessing
import pandas as pd
from pathlib import Path

from .constants import WILD_COL, MUTANT_COL
from .cmdl_calls import call_microminer_search


# def run_pairs(df: pd.DataFrame, outdir: Path, cpus: int = 1, pdb_mirror: str = 'standard'):
#     # generate a list of parameter tuple containing the PDB pairs to be aligned.
#     name_pairs = [(getattr(row, WILD_COL),
#                    getattr(row, MUTANT_COL),
#                    outdir / '{}_{}'.format(getattr(row, WILD_COL), getattr(row, MUTANT_COL)),
#                    pdb_mirror)
#                   for row in df[[WILD_COL, MUTANT_COL]].drop_duplicates().itertuples(index=False)]
#
#     # parallel execution with each parameter set separately
#     with multiprocessing.Pool(cpus) as pool:
#         # a single pair job is very fast. We just use the naive chunking
#         chunksize = int(len(name_pairs) / cpus) + 1
#         itertools.chain.from_iterable(pool.starmap(call_mutscreen_pairs, name_pairs, chunksize=chunksize))
#
#
# def run_screens(df: pd.DataFrame, outdir: Path, cpus: int = 1, pdb_mirror: str = 'standard', raise_error: bool=True):
#     # generate a list of parameter tuple containing the PDB pairs to be aligned.
#     parameter_set = [(getattr(row, WILD_COL),
#                       outdir / '{}'.format(getattr(row, WILD_COL)),
#                       pdb_mirror,
#                       raise_error)
#                       for row in df[[WILD_COL]].drop_duplicates().itertuples(index=False)]
#
#     # parallel execution with each parameter set separately
#     with multiprocessing.Pool(cpus) as pool:
#         # take default chunksize. Works well.
#         pool.starmap(call_mutscreen_screen, parameter_set,)


def run_searches(df: pd.DataFrame, outdir: Path, cpus: int = 1, raise_error: bool=True):
    # generate a list of parameter tuple
    parameter_set = [(Path(getattr(row, 'structure_path')),
                      outdir / '{}'.format(getattr(row, 'id')),
                      raise_error)
                      for row in df.drop_duplicates().itertuples(index=False)]

    # parallel execution with each parameter set separately
    with multiprocessing.Pool(cpus) as pool:
        # take default chunksize. Works well.
        pool.starmap(call_microminer_search, parameter_set,)


