import logging
import subprocess
from pathlib import Path
from typing import List

from . import utils
from .constants import CONFIG

logger = logging.getLogger(__name__)


def exe_cmdl_call(cmd_call: List[str], log_msg: str, raise_error: bool) -> dict:
    """Execute a command line call.

    :param cmd_call: The command line call as list.
    :param log_msg: A message to log.
    :param raise_error: Whether to raise an error when the cmdl call returns != 0.
    :return: Dict containing the exit_code, standard out and standard error of the call.
    """
    logger.info(f'Calling {" ".join(cmd_call)}')
    process = subprocess.Popen(cmd_call, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with utils.timer(f"{log_msg} | " + " ".join(cmd_call)):
        stdout, stderr = process.communicate()
        exit_code = process.wait()
    if exit_code != 0:
        msg = (
            f'Failed call {" ".join(cmd_call)}\nstdout={stdout.decode()}\n'
            f"stderr={stderr.decode()}"
        )
        logging.error(msg)
        if raise_error:
            raise ValueError(msg)
    return {"exit_code": exit_code, "stdout": stdout, "stderr": stderr}


def call_microminer_search(
    pdb_query_path: Path,
    outdir: Path,
    mode: str = "single_mutation",
    mm_repr: str = "monomer",
    raise_error: bool = True,
) -> dict:
    """Calls the MicroMiner executable in search mode.

    :param pdb_query_path: Path to query PDB file
    :param outdir: Result dir path.
    :param mode: Search mode.
    :param mm_repr: Structure representation mode.
    :param raise_error: Whether to raise exception if command line call return != 0
    :return: Dict with details on the tool execution.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    exe = Path(CONFIG["EXECUTABLES"]["MICROMINER"])
    sitesearchdb = Path(CONFIG["DATABASES"]["SITE_SEARCH_DB"])

    cmd_call = [
        str(exe.resolve()),
        "search",
        "-q",
        str(pdb_query_path.resolve()),
        "-s",
        str(sitesearchdb.resolve()),
        "-o",
        str(outdir.resolve()),
        # '-d',  #write ensembles to disc
        "--cpus",
        str(CONFIG["MICROMINER_ALGO"]["CPUS"]),
        "--site_radius",
        str(CONFIG["MICROMINER_ALGO"]["SITE_RADIUS"]),
        "--identity",
        str(CONFIG["MICROMINER_ALGO"]["IDENTITY"]),
        "--fragment_length",
        str(CONFIG["MICROMINER_ALGO"]["FRAGMENT_LENGTH"]),
        "--fragment_distance",
        str(CONFIG["MICROMINER_ALGO"]["FRAGMENT_DISTANCE"]),
        "--score_threshold",
        str(CONFIG["MICROMINER_ALGO"]["SCORE_THRESH"]),
        "--flexibility_sensitivity",
        str(CONFIG["MICROMINER_ALGO"]["FLEXIBILITY_SENSITIVITY"]),
        "--kmer_matching_rate",
        str(CONFIG["MICROMINER_ALGO"]["KMER_MATCHING_RATE"]),
        "-m",
        mode,
        "-r",
        mm_repr,
    ]

    response = exe_cmdl_call(cmd_call, "MicroMiner Search", raise_error)
    response["params"] = " ".join(cmd_call)
    return response


def call_microminer_pair(
    pdb_query_path: Path, pdb_target_path: Path, outdir: Path, raise_error: bool = True
) -> dict:
    """Calls the MicroMiner executable in pair mode.

    :param pdb_query_path: Path to query PDB file.
    :param pdb_target_path: Path to target PDB file.
    :param outdir: Directory for writting results.
    :param raise_error: Whether to raise an exception on failure of the command line tool.
    :return: Dict with details on the tool execution.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    exe = Path(CONFIG["EXECUTABLES"]["MICROMINER"])

    cmd_call = [
        str(exe.resolve()),
        "site_align",
        "-q",
        str(pdb_query_path.resolve()),
        "-t",
        str(pdb_target_path.resolve()),
        "-o",
        str(outdir.resolve()),
        "--cpus",
        str(CONFIG["MICROMINER_ALGO"]["CPUS"]),
        "--site_radius",
        str(CONFIG["MICROMINER_ALGO"]["SITE_RADIUS"]),
        "--identity",
        str(CONFIG["MICROMINER_ALGO"]["IDENTITY"]),
        "--fragment_length",
        str(CONFIG["MICROMINER_ALGO"]["FRAGMENT_LENGTH"]),
        "--flexibility_sensitivity",
        str(CONFIG["MICROMINER_ALGO"]["FLEXIBILITY_SENSITIVITY"]),
    ]
    return exe_cmdl_call(cmd_call, "MicroMiner site_align", raise_error)


def call_tmalign(
    pdb_query_path: Path,
    pdb_target_path: Path,
    outdir: Path,
    out_prefix: str,
    write_rotation: bool,
    split_flag: int = 0,
    ter_flag: int = None,
    raise_error: bool = True,
    create_outdir: bool = True,
) -> dict:
    """Call the TM-Align tool.

    Y. Zhang, J. Skolnick, TM-align: A protein structure alignment algorithm based on TM-score,
    Nucleic Acids Research, 33: 2302-2309 (2005)

    See the TM-Align executable and help message for details on the parameters.

    :param pdb_query_path: Path to query PDF file.
    :param pdb_target_path: Path to target PDB file.
    :param outdir: Directory for result files.
    :param out_prefix: Prefix for result files.
    :param write_rotation: Whether to write rotation matrix file of the alignments' superposition.
    :param split_flag: Split strategy of TM-Align (see TM-Align help for details)
    :param ter_flag: Ter-flag used by TM-Align (see TM-Align help for details)
    :param raise_error: Whether to raise an exception if command line call fails.
    :param create_outdir: Whether to create the outdir if it not exists.
    :return: Dict with details on the tool execution.
    """

    if create_outdir:
        outdir.mkdir(parents=True, exist_ok=True)

    exe = Path(CONFIG["EXECUTABLES"]["TMALIGN"])

    cmd_call = [
        str(exe.resolve()),
        str(pdb_query_path.resolve()),
        str(pdb_target_path.resolve()),
        "-o",
        str(outdir.resolve()),
        "-split",
        str(split_flag),
    ]
    if write_rotation:
        cmd_call.extend(["-m", str((outdir / f"{out_prefix}_rotation").resolve())])
    if ter_flag is not None:
        cmd_call.extend(["-ter", str(ter_flag)])

    return exe_cmdl_call(cmd_call, "TMAlign", raise_error)
