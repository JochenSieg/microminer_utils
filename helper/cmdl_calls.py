from pathlib import Path
import subprocess

from typing import List

from .constants import CONFIG
from . import utils
import logging

logger = logging.getLogger(__name__)


def exe_cmdl_call(cmd_call: List[str], log_msg: str, raise_error: bool) -> dict:
    logger.info(f'Calling {" ".join(cmd_call)}')
    process = subprocess.Popen(cmd_call, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with utils.timer(f'{log_msg} | ' + " ".join(cmd_call)):
        stdout, stderr = process.communicate()
        exit_code = process.wait()
    if exit_code != 0:
        msg = f'Failed call {" ".join(cmd_call)}\nstdout={stdout.decode()}\nstderr={stderr.decode()}'
        logging.error(msg)
        if raise_error:
            raise ValueError(msg)
    return {'exit_code': exit_code, 'stdout': stdout, 'stderr': stderr}


def call_microminer_search(pdb_query_path: Path, outdir: Path, raise_error: bool = True):
    """
    Calls the MicroMiner executable in search mode
    :param pdb_query_path: Path to query PDB file
    :param outdir: Result dir path.
    :param raise_error: Whether to raise exception if command line call return != 0
    :return:
    """
    outdir.mkdir(parents=True, exist_ok=True)

    exe = Path(CONFIG['EXECUTABLES']['MICROMINER'])
    sitesearchdb = Path(CONFIG['DATABASES']['SITE_SEARCH_DB'])

    complexdb = None
    if 'COMPLEX_DB' in CONFIG['DATABASES']:
        complexdb = Path(CONFIG['DATABASES']['COMPLEX_DB'])

    cmd_call = [str(exe.resolve()),
                'search',
                '-q', str(pdb_query_path.resolve()),
                '-s', str(sitesearchdb.resolve()),
                '-o', str(outdir.resolve()),
                # '-d',  #write ensembles to disc
                '--cpus', str(CONFIG['MICROMINER_ALGO']['CPUS']),
                '--site_radius', str(CONFIG['MICROMINER_ALGO']['SITE_RADIUS']),
                '--identity', str(CONFIG['MICROMINER_ALGO']['IDENTITY']),
                '--fragment_length', str(CONFIG['MICROMINER_ALGO']['FRAGMENT_LENGTH']),
                '--flexibility_sensitivity',
                str(CONFIG['MICROMINER_ALGO']['FLEXIBILITY_SENSITIVITY']),
                '--kmer_matching_rate', str(CONFIG['MICROMINER_ALGO']['KMER_MATCHING_RATE']),
                ]

    # TODO es gibt wohl ein paar PDBs die vom file ohne Probleme gelesen werden, aber wenn
    #      sie aus der ComplexLib gelesen werde failed der chemie check...
    #      Bsp: 1AAR als query
    # if complexdb:
    #     cmd_call.extend(['-b', str(complexdb.resolve())])
    exe_cmdl_call(cmd_call, 'MicroMiner Search', raise_error)


def call_microminer_pair(pdb_query_path: Path, pdb_target_path: Path, outdir: Path,
                         raise_error: bool = True):
    """
    Calls the MicroMiner executable in pair mode.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    exe = Path(CONFIG['EXECUTABLES']['MICROMINER'])

    cmd_call = [str(exe.resolve()),
                'site_align',
                '-q', str(pdb_query_path.resolve()),
                '-t', str(pdb_target_path.resolve()),
                '-o', str(outdir.resolve()),
                '--cpus', str(CONFIG['MICROMINER_ALGO']['CPUS']),
                '--site_radius', str(CONFIG['MICROMINER_ALGO']['SITE_RADIUS']),
                '--identity', str(CONFIG['MICROMINER_ALGO']['IDENTITY']),
                '--fragment_length', str(CONFIG['MICROMINER_ALGO']['FRAGMENT_LENGTH']),
                '--flexibility_sensitivity',
                str(CONFIG['MICROMINER_ALGO']['FLEXIBILITY_SENSITIVITY']),
                ]
    exe_cmdl_call(cmd_call, 'MicroMiner site_align', raise_error)


def call_microminer_prefilter(pdb_query_path: Path, outdir: Path, raise_error: bool = True):
    """
    Calls the MicroMiner executable in search mode
    :param pdb_query_path: Path to query PDB file
    :param outdir: Result dir path.
    :param raise_error: Whether to raise exception if command line call return != 0
    :return:
    """
    outdir.mkdir(parents=True, exist_ok=True)

    exe = Path(CONFIG['EXECUTABLES']['MICROMINER'])
    sitesearchdb = Path(CONFIG['DATABASES']['SITE_SEARCH_DB'])

    cmd_call = [str(exe.resolve()),
                'prefilter',
                '-q', str(pdb_query_path.resolve()),
                '-s', str(sitesearchdb.resolve()),
                '-o', str(outdir.resolve()),
                ]

    exe_cmdl_call(cmd_call, 'MicroMiner Prefilter', raise_error)


def call_tmalign(pdb_query_path: Path, pdb_target_path: Path, outdir: Path, out_prefix: str,
                 write_rotation: bool, split_flag: int=0,
                 ter_flag: int=None,
                 raise_error: bool = True):
    """

    :param pdb_query_path:
    :param pdb_target_path:
    :param outdir:
    :param out_prefix:
    :param write_rotation:
    :param raise_error:
    :return:
    """
    outdir.mkdir(parents=True, exist_ok=True)

    exe = Path(CONFIG['EXECUTABLES']['TMALIGN'])

    cmd_call = [str(exe.resolve()),
                str(pdb_query_path.resolve()),
                str(pdb_target_path.resolve()),
                '-o', str(outdir.resolve()),
                '-split', str(split_flag),
                ]
    if write_rotation:
        cmd_call.extend(['-m', str((outdir / f'{out_prefix}_rotation').resolve())])
    if ter_flag is not None:
        cmd_call.extend(['-ter', str(ter_flag)])

    proc_dict = exe_cmdl_call(cmd_call, 'TMAlign', raise_error)
    if proc_dict['exit_code'] == 0:
        return proc_dict['stdout']
    else:
        return None

#
# def ALTcall_tmalign(pdbid1: str, pdbid2: str, outdir: Path, use_monomers=False):
#
#     exe = Path(CONFIG['EXECUTABLES']['TMALIGN_SCRIPT'])
#
#     if use_monomers:
#         # is use_monomer the pdbids must look like 1g9v_A with chain id
#         p1 = utils.generate_and_get_pdb_monomer(*pdbid1.split('_'), outdir)
#         p2 = utils.generate_and_get_pdb_monomer(*pdbid2.split('_'), outdir)
#     else:
#         p1 = utils.get_pdb_file_path(pdbid1)
#         p2 = utils.get_pdb_file_path(pdbid2)
#
#     cmd_call = ['sh',
#                 str(exe.resolve()),
#                 str(p1.resolve()),
#                 str(p2.resolve()),
#                 str(outdir.resolve())
#                 ]
#     # logging.info('Calling CMD: {}'.format(' '.join(cmd_call)))
#
#     res = subprocess.run(cmd_call, stdout=subprocess.PIPE, cwd=outdir.resolve())
#     if res.returncode != 0:
#         logging.error(res)
#         raise ValueError('Encountered problem {}'.format(res))

#
# def call_mutscreen_pairs(pdbid1: str, pdbid2: str, outdir: Path, pdb_mirror1='standard', pdb_mirror2='standard'):
#
#     outdir.mkdir(parents=True, exist_ok=True)
#
#     exe = Path(CONFIG['EXECUTABLES']['MICROMINER'])
#
#     p1 = utils.get_pdb_file_path(pdbid1, mirror=pdb_mirror1)
#     p2 = utils.get_pdb_file_path(pdbid2, mirror=pdb_mirror2)
#
#     cmd_call = ['time', str(exe.resolve()),
#                 'pair',
#                 '-q', str(p1.resolve()),
#                 '-t', str(p2.resolve()),
#                 '-o', str(outdir.resolve()),
#                 '--cpus', str(CONFIG['MICROMINER_ALGO']['CPUS']),
#                 '--site_radius', str(CONFIG['MICROMINER_ALGO']['SITE_RADIUS']),
#                 '--identity', str(CONFIG['MICROMINER_ALGO']['IDENTITY']),
#                 '--fragment_length', str(CONFIG['MICROMINER_ALGO']['FRAGMENT_LENGTH']),
#                 '--flexibility_sensitivity',
#                 str(CONFIG['MICROMINER_ALGO']['FLEXIBILITY_SENSITIVITY']),
#                 '--kmer_matching_rate', str(CONFIG['MICROMINER_ALGO']['KMER_MATCHING_RATE']),
#                 ]
#
#     logging.info('Calling CMD: {}'.format(' '.join(cmd_call)))
#
#     res = subprocess.run(cmd_call, stdout=subprocess.PIPE)
#     if res.returncode != 0:
#         logging.error(res)
#         raise ValueError('Encountered problem {}'.format(res))
#
#
# def call_mutscreen_screen(query_pdbid: str, outdir: Path, pdb_mirror: str, raise_error: bool = True):
#
#     outdir.mkdir(parents=True, exist_ok=True)
#
#     exe = Path(CONFIG['EXECUTABLES']['MICROMINER'])
#     sitesearchdb = Path(CONFIG['DATABASES']['SITE_SEARCH_DB'])
#
#     complexdb = None
#     if 'COMPLEX_DB' in CONFIG['DATABASES']:
#         complexdb = Path(CONFIG['DATABASES']['COMPLEX_DB'])
#
#     pdb_query_path = utils.get_pdb_file_path(query_pdbid, mirror=pdb_mirror)
#
#     cmd_call = ['time', str(exe.resolve()),
#                 'search',
#                 '-q', str(pdb_query_path.resolve()),
#                 '-s', str(sitesearchdb.resolve()),
#                 '-o', str(outdir.resolve()),
#                 # '-d',  #write ensembles to disc
#                 '--cpus', str(CONFIG['MICROMINER_ALGO']['CPUS']),
#                 '--site_radius', str(CONFIG['MICROMINER_ALGO']['SITE_RADIUS']),
#                 '--identity', str(CONFIG['MICROMINER_ALGO']['IDENTITY']),
#                 '--fragment_length', str(CONFIG['MICROMINER_ALGO']['FRAGMENT_LENGTH']),
#                 '--flexibility_sensitivity',
#                 str(CONFIG['MICROMINER_ALGO']['FLEXIBILITY_SENSITIVITY']),
#                 '--kmer_matching_rate', str(CONFIG['MICROMINER_ALGO']['KMER_MATCHING_RATE']),
#                 ]
#
#     # TODO es gibt wohl ein paar PDBs die vom file ohne Probleme gelesen werden, aber wenn
#     #      sie aus der ComplexLib gelesen werde failed der chemie check...
#     #      Bsp: 1AAR als query
#     # if complexdb:
#     #     cmd_call.extend(['-b', str(complexdb.resolve())])
#
#     logging.info('Calling MutScreen Screen CMD: {}'.format(' '.join(cmd_call)))
#
#     res = subprocess.run(cmd_call, stdout=subprocess.PIPE)
#     if res.returncode != 0:
#         logging.error(res)
#         if raise_error:
#             raise ValueError('Encountered problem {}'.format(res))
#
#
# def call_pdb2monomers(pdbid: str, outdir: Path):
#
#     exe = Path(CONFIG['EXECUTABLES']['PDB_2_MONOMER'])
#
#     p = utils.get_pdb_file_path(pdbid)
#
#     cmd_call = [str(exe.resolve()),
#                 '-c', str(p.resolve()),
#                 '-o', str(outdir.resolve()),
#                 ]
#
#     # logging.info('Calling CMD: {}'.format(' '.join(cmd_call)))
#
#     res = subprocess.run(cmd_call, stdout=subprocess.PIPE)
#     if res.returncode != 0:
#         logging.error(res)
#         raise ValueError('Encountered problem {}'.format(res))
