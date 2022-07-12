import tempfile

import pandas as pd
from pathlib import Path
import multiprocessing
import itertools
import logging

from .cmdl_calls import call_microminer_search, call_microminer_pair, call_microminer_prefilter
from .tmalign import call_and_parse

logger = logging.getLogger(__name__)


def _run_parallel(func, parameter_set, cpus):
    # parallel execution with each parameter set separately
    with multiprocessing.Pool(cpus) as pool:
        pool.starmap(func, parameter_set, )


class MicroMinerSearch:
    MANDATORY_TSV_COLUMNS = ['id', 'structure_path']

    def __init__(self, cpus: int = 1, raise_error: bool = True):
        self.cpus = cpus
        self.raise_error = raise_error

    def run(self, param_tsv: Path, outdir: Path):
        df = pd.read_csv(param_tsv, sep='\t', header=0)
        logger.info(f'Read {df.shape[0]} parameter records for computation')

        if not all(c in df.columns for c in MicroMinerSearch.MANDATORY_TSV_COLUMNS):
            raise ValueError(f'Missing mandatory fields in TSV: {param_tsv}')

        # generate a list of parameter tuple containing the PDB pairs to be aligned.
        parameter_set = [(Path(getattr(row, MicroMinerSearch.MANDATORY_TSV_COLUMNS[1])),
                          outdir / '{}'.format(getattr(row,
                                                       MicroMinerSearch.MANDATORY_TSV_COLUMNS[0])),
                          self.raise_error)
                         for row in df.drop_duplicates().itertuples(index=False)]

        if self.cpus > 1:
            _run_parallel(call_microminer_search, parameter_set, self.cpus)
        else:
            for param_set in parameter_set:
                call_microminer_search(*param_set)


class MicroMinerPair:
    MANDATORY_TSV_COLUMNS = ['id1', 'structure_path1', 'id2', 'structure_path2']

    def __init__(self, cpus: int = 1, raise_error: bool = True):
        self.cpus = cpus
        self.raise_error = raise_error

    def run(self, param_tsv: Path, outdir: Path):
        df = pd.read_csv(param_tsv, sep='\t', header=0)
        logger.info(f'Read {df.shape[0]} parameter records for computation')

        if not all(c in df.columns for c in MicroMinerPair.MANDATORY_TSV_COLUMNS):
            raise ValueError(f'Missing mandatory fields in TSV: {param_tsv}')

        # generate a list of parameter tuple containing the PDB pairs to be aligned.
        parameter_set = [(Path(getattr(row, MicroMinerPair.MANDATORY_TSV_COLUMNS[1])),
                          Path(getattr(row, MicroMinerPair.MANDATORY_TSV_COLUMNS[3])),
                          outdir / '{}_{}'.format(getattr(row,
                                                          MicroMinerPair.MANDATORY_TSV_COLUMNS[0]),
                                                  getattr(row,
                                                          MicroMinerPair.MANDATORY_TSV_COLUMNS[2])
                                                  ),
                          self.raise_error)
                         for row in df.drop_duplicates().itertuples(index=False)]

        # TODO named args would be better for the parameter set to avoid position depending errors
        # parameter_set = [{'pdb_query_path': Path(getattr(row, MicroMinerPair.MANDATORY_TSV_COLUMNS[1])),
        #                  'pdb_target_path': Path(getattr(row, MicroMinerPair.MANDATORY_TSV_COLUMNS[3])),
        #                  'outdir': outdir / '{}_{}'.format(getattr(row,
        #                                                   MicroMinerPair.MANDATORY_TSV_COLUMNS[0]),
        #                                           getattr(row,
        #                                                   MicroMinerPair.MANDATORY_TSV_COLUMNS[2])
        #                                           ),
        #                  'raise_error': self.raise_error}
        #                  for row in df.drop_duplicates().itertuples(index=False)]

        if self.cpus > 1:
            _run_parallel(call_microminer_pair, parameter_set, self.cpus)
        else:
            for param_set in parameter_set:
                call_microminer_pair(*param_set)


class MicroMinerPrefilter:
    MANDATORY_TSV_COLUMNS = ['id', 'structure_path']

    def __init__(self, cpus: int = 1, raise_error: bool = True):
        self.cpus = cpus
        self.raise_error = raise_error

    def run(self, param_tsv: Path, outdir: Path):
        df = pd.read_csv(param_tsv, sep='\t', header=0)
        logger.info(f'Read {df.shape[0]} parameter records for computation')

        if not all(c in df.columns for c in MicroMinerPrefilter.MANDATORY_TSV_COLUMNS):
            raise ValueError(f'Missing mandatory fields in TSV: {param_tsv}')

        # generate a list of parameter tuple containing the PDB pairs to be aligned.
        parameter_set = [(Path(getattr(row, MicroMinerPrefilter.MANDATORY_TSV_COLUMNS[1])),
                          outdir / '{}'.format(getattr(row,
                                                       MicroMinerPrefilter.MANDATORY_TSV_COLUMNS[
                                                           0])),
                          self.raise_error)
                         for row in df.drop_duplicates().itertuples(index=False)]

        if self.cpus > 1:
            _run_parallel(call_microminer_prefilter, parameter_set, self.cpus)
        else:
            for param_set in parameter_set:
                call_microminer_prefilter(*param_set)


class TMalign:
    MANDATORY_TSV_COLUMNS = ['id1', 'structure_path1', 'id2', 'structure_path2']

    def __init__(self, cpus: int = 1, skip_same: bool = False, raise_error: bool = True):
        self.cpus = cpus
        self.skip_same = skip_same
        self.raise_error = raise_error

    def run(self, param_tsv: Path, outdir: Path = None):
        df = pd.read_csv(param_tsv, sep='\t', header=0)
        logger.info(f'Read {df.shape[0]} parameter records for computation')

        if not all(c in df.columns for c in self.MANDATORY_TSV_COLUMNS):
            raise ValueError(f'Missing mandatory fields in TSV: {param_tsv}')

        if self.skip_same:
            # drop rows where query and target are identical by ID and filepath
            id1 = self.MANDATORY_TSV_COLUMNS[0]
            id2 = self.MANDATORY_TSV_COLUMNS[2]
            p1 = self.MANDATORY_TSV_COLUMNS[1]
            p2 = self.MANDATORY_TSV_COLUMNS[3]
            df.query(f'{id1} != {id2} or {p1} != {p2}', inplace=True)
            logger.info(f'Filtered out same inputs. {df.shape[0]} remaining.')

        # df = df.head(n=100)

        with tempfile.TemporaryDirectory() as t:
            tmpdir_path = Path(t)
            # generate a list of parameter tuple containing the PDB pairs to be aligned.
            parameter_set = [(Path(getattr(row, self.MANDATORY_TSV_COLUMNS[1])),
                              Path(getattr(row, self.MANDATORY_TSV_COLUMNS[3])),
                              tmpdir_path,
                              False,  # don't write rotation matrix file
                              self.raise_error)
                             for row in df.drop_duplicates().itertuples(index=False)]

            with multiprocessing.Pool(self.cpus) as pool:
                df_res = pd.DataFrame(pool.starmap(call_and_parse, parameter_set))

            # print(df_res)

            if outdir is not None:
                df_res.to_csv(outdir / f'{param_tsv.stem}_tmalign.tsv', sep='\t', header=True,
                              index=False)

            return df_res

        # df_res['name1'] = df_res['name1'].str[:4].str.upper() + df_res['name1'].str[
        #                                                         4:]  # first 4 chars to upper
        # df_res['name2'] = df_res['name2'].str[:4].str.upper() + df_res['name2'].str[
        #                                                         4:]  # first 4 chars to upper
        # df_res = df.merge(df_res, left_on=[self.id_col1, self.id_col2], right_on=['name2', 'name1'],
        #                   how='left')
        # TODO diese assert fliegt, wenn ich alle (chain) pairs in den mutscreen results auf der ThermoMutDB als input nehme?
        #      sind da doppelte im df drin, dass er mehr als eine Zeile fuer jede Zeile im Input generiert? Das sollte eigentlich nicht sein,
        #      weil ich die Berechnungen ja nur mit einem echten subset des df_in starte. Da ist ja alles unique... mach ich irgendwann mal upper
        #      oder lower wodurch dann zwei eigentlich unterschiedliche namen gleich heissen? ODer laueft generell was anderes schief??
        #   BEMERKE: zur Zeit dieses Kommentar hies df -> df_in und df_res -> df
        # if df.shape[0] != df_res.shape[0]:
        #     print('DAS ASSERT FLIEGT IMMER NOCH!! Was ist da los?')
        # assert (df.shape[0] == df_res.shape[0])
        # df_res.drop('name1', axis=1, inplace=True)
        # df_res.drop('name2', axis=1, inplace=True)

        return df_res
