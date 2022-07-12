import contextlib
import logging
from pathlib import Path
import os
import time

import pandas as pd
from .constants import CONFIG, WILD_COL, MUTANT_COL
from . import cmdl_calls
from .datasets import skempi2, platinum, scope

logger = logging.getLogger(__name__)


def get_divided_dir_name(pdbid: str) -> str:
    if len(pdbid) == 4 or len(pdbid) == 6 and pdbid[4] == '_':
        return pdbid[1:3].lower()
    raise ValueError('Do not recognize format of pdbid: {}'.format(pdbid))


def get_pdb_file_path(pdbid: str, allow_obsolete: bool = True, ignore_missing=True,
                      mirror: str = 'standard') -> Path:
    """
    Tries different PDB file sources to provide a path to the requested PDB file.
    :param pdbid:
    :param allow_obsolete:
    :param ignore_missing: Ignore if PDB files are missing. Can happen if PDB is invalid or
                           requested structure file is only available as CIF which NAOMI
                           does not support fully, currently. So we ignore these cases.
    :param mirror: PDB mirror to use. Default is 'standard'. If non applicable default is used.
    :return:
    """
    if mirror == 'skempi2':
        return skempi2.read_skempi_cleaned_pdb(pdbid)
    elif mirror == 'platinum':
        return platinum.read_platinum_cleaned_pdb(pdbid)
    elif mirror == 'scope':
        return scope.read_scope_pdbstyle(pdbid)
    return get_pdb_file_path_from_standard(pdbid, allow_obsolete, ignore_missing)


def get_pdb_file_path_from_standard(pdbid: str, allow_obsolete: bool = True,
                                    ignore_missing = False) -> Path:
    """
    Gets a file path to provided PDB Id if a corresponding file exists.
    :param pdbid: the PDB Id
    :param allow_obsolete: allow obsolete structures.
    :param ignore_missing: Do not raise if PDB file can not be found. Return None instead.
    :return: path to PDB file
    """
    pdb_dir = Path(CONFIG['DATA']['PDB_DIR'])
    filename = '{}{}{}'.format(CONFIG['PDB_FILE_INFO']['PREFIX'],
                               pdbid.lower() if CONFIG['PDB_FILE_INFO']['CASE'] == 'LOWER'
                               else pdbid.upper(),
                               CONFIG['PDB_FILE_INFO']['SUFFIX'])
    path = pdb_dir / filename
    if not path.is_file():
        if allow_obsolete and CONFIG['DATA']['PDB_OBSOLETE_DIR']:
            path = Path(CONFIG['DATA']['PDB_OBSOLETE_DIR']) / get_divided_dir_name(pdbid) / filename
            if path.is_file():
                return path
        if not ignore_missing:
            raise FileNotFoundError('Could not find PDB file: {}'.format(path))
        else:
            return None
    return path


def get_monomer_pdb_file_path(pdbid_chain: str) -> Path:
    """
    Gets a file path to provided PDB Id chain id if a corresponding file exists.
    :param pdbid_chain: the PDB Id and chain id. Example: 100D_A
    :return: path to PDB monomer file
    """
    pdb_dir = Path(CONFIG['DATA']['MONOMER_DIR'])
    path = pdb_dir / '{}.pdb'.format(pdbid_chain)
    if not path.is_file():
        raise FileNotFoundError('Could not found PDB monomer file: {}'.format(path))
    return path


def get_renamed_column_names(cols, rename_dict):
    """
    This is a workaround because there seems to be a problem with pandas rename function
    which raises a KeyError even when assert(key in df.columns) is True. Bug?

    This does raises the keyerror work:

        assert(THERMOMUTDB_PDB_WILD_JSON in df_single_mutation.columns)
        assert(THERMOMUTDB_PDB_MUTANT_JSON in df_single_mutation.columns)
        df_single_mutation.rename({THERMOMUTDB_PDB_WILD_JSON: helper.WILD_COL,
                                   THERMOMUTDB_PDB_MUTANT_JSON: helper.MUTANT_COL},
                                   errors='raise', inplace=True)

    :return: New columns list with renamed entries (if present).
    """
    new_cols = []
    for c in cols:
        if c in rename_dict.keys():
            new_cols.append(rename_dict[c])
        else:
            new_cols.append(c)
    return new_cols


def scantree(path: Path):
    """
    Recursively yield DirEntry objects for given directory.
    :param path: Directory path.
    :return: yields all file system entries (dirs and files) as Path object.
    """
    if path.is_file():
        yield path
    else:
        yield from _scantree(path)


def _scantree(path: Path):
    """
    Private function for yielding DirEntry objects for given directory recursively.
    :param path: Directory path.
    :return: yields all file system entries (dirs and files) as Path object.
    """
    for entry in os.scandir(path):
        if entry.is_dir(follow_symlinks=False):
            yield from _scantree(Path(entry.path))
        else:
            yield Path(entry.path)


def rename_mutscreen_wild_mutant_col(df: pd.DataFrame, wild: str = WILD_COL, mutant: str = MUTANT_COL) -> pd.DataFrame:
    """
    Rename mutscreen output columns for wild-type and mutant PDB id to specified names.
    :param df: The input table.
    :param wild: New name for wild column.
    :param mutant: New name for mutant column.
    :return: Table with renamed columns.
    :raise KeyError if mutscreen column names are not present
    """
    return df.rename({'wildName': wild, 'mutantName': MUTANT_COL}, axis=1, errors='raise')


def rename_to_lib_standard(df: pd.DataFrame, wild: str, mutant: str) -> pd.DataFrame:
    """
    Rename mutscreen output columns for wild-type and mutant PDB id to specified names.
    :param df: The input table.
    :param wild: New name for wild column.
    :param mutant: New name for mutant column.
    :return: Table with renamed columns.
    :raise KeyError if mutscreen column names are not present
    """
    return df.rename({wild: WILD_COL, mutant: MUTANT_COL}, axis=1, errors='raise')


def generate_and_get_pdb_monomer(pdbid: str, chainid: str, outdir: Path) -> Path:
    """
    Splits PDB file into single chains and returns the path to the monomer of the specified chain.
    :param pdbid:
    :param chainid:
    :return:
    """
    cmdl_calls.call_pdb2monomers(pdbid, outdir)
    path = outdir / '{}_{}.pdb'.format(pdbid, chainid)
    if not path.is_file():
        logging.warning('Warning: Monomer file not found: {}'.format(path.absolute()))
    return path


def get_switched_column_names(df: pd.DataFrame, col1: str, col2: str) -> list:
    # could also use pandas rename function.
    col_list = list(df.columns)
    a, b = col_list.index(col2), col_list.index(col1)
    col_list[b], col_list[a] = col_list[a], col_list[b]
    return col_list


def build_cache_dict(cache_dirs: list) -> dict:
    """
    Searches recursively in directory tree for 'resultStatistic.csv' files. It is assumed that the parent dir of the
    'resultStatistic.csv' is named with the PDB id that was used to generate the 'resultStatistic.csv'. A dict
    mapping upper case PDB IDs to file system paths is returned.

    Note that the 'resultStatistic.csv' files in the cache_dirs should not mix different query types, like screen,
    pair or single site queries only. Otherwise the cached results might not be equivalent to what is expected to be
    computed.

    Only use when the standard PDB mirror is used for wild and mutant! Custom prepared PDB files are not supported.
    :param cache_dirs: Directory to index 'resultStatistic.csv' files in.
    :return: dict mapping upper case PDB ids to paths of already computed 'resultStatistic.csv'.
    """
    # exclude these data sets because they are based on custom prepared PDB files for which we do not need/want caching.
    exclude_dataset = ['platinum', 'skempi2']

    def my_scantree(path):
        for entry in os.scandir(path):
            if entry.is_dir(follow_symlinks=False) and entry.name in exclude_dataset:
                yield from scantree(entry.path)
            else:
                yield entry.path

    cache_dict = {}
    for cache_dir in cache_dirs:
        cache_dir = Path(cache_dir)
        if cache_dir.is_dir():
            for entry in my_scantree(cache_dir):
                if entry.name == 'resultStatistic.csv':
                    path = Path(entry)
                    cache_dict[path.parent.name.upper()] = path.absolute()

    return cache_dict


@contextlib.contextmanager
def timer(message: str):
    logger.info(f'Starting: {message}')
    tic = time.time()
    yield
    toc = time.time()
    logger.info(f'Finished: {message} Took: {(toc - tic):.3f} seconds')
