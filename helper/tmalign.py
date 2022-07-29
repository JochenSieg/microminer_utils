"""Utils functionality for TM-align like reading and parsing results"""
import logging
import os
import tempfile
from pathlib import Path
from typing import List, Dict

from .cmdl_calls import call_tmalign
from .utils import unpack_gz

logger = logging.getLogger(__name__)


def parse_stdout(stdout_str: str, read_alignment: bool = False) -> List[Dict]:
    """Parse the standard out output of TM-Align to a list of dictionaries.

    :param stdout_str: The standard out string to parse.
    :param read_alignment: Whether to also parse the alignment.
    :return: List of dicts. One dict for each pair-wise alignment in the output.
    """
    info_dicts = [{}]
    next_line_is_alignment = False
    align_counter = 0
    for i, line in enumerate(stdout_str.split(os.linesep)):
        if i == 2 and not line.startswith(" * TM-align"):
            return []
        if i > 2 and line.startswith(" * TM-align"):
            # new entry
            next_line_is_alignment = False
            align_counter = 0
            info_dicts.append({})

        if next_line_is_alignment:
            if align_counter == 0:
                info_dicts[-1]["refseq"] = line.strip("\n")
            elif align_counter == 1:
                info_dicts[-1]["matching"] = line.strip("\n")
            elif align_counter == 2:
                info_dicts[-1]["targetseq"] = line.strip("\n")
            else:
                next_line_is_alignment = False
            align_counter += 1

        elif line.startswith("Name of Chain_1"):
            p = Path(line.split(" ")[3])
            info_dicts[-1]["id1"] = p.stem
            name_split = p.name.split(":")
            if len(name_split) > 1 and len(name_split[-1]) < 3:
                info_dicts[-1]["chain_id1"] = name_split[-1]
        elif line.startswith("Name of Chain_2"):
            p = Path(line.split(" ")[3])
            info_dicts[-1]["id2"] = p.stem
            name_split = p.name.split(":")
            if len(name_split) > 1 and len(name_split[-1]) < 3:
                info_dicts[-1]["chain_id2"] = name_split[-1]
        elif line.startswith("Length of Chain_1"):
            info_dicts[-1]["chain1_len"] = int(line.split(" ")[3])
        elif line.startswith("Length of Chain_2"):
            info_dicts[-1]["chain2_len"] = int(line.split(" ")[3])
        elif line.startswith("Aligned length="):
            splitted = line.split("=")
            info_dicts[-1]["align_len"] = int(splitted[1].strip().split(",")[0])
            info_dicts[-1]["rmsd"] = float(splitted[2].strip().split(",")[0])
            info_dicts[-1]["seqid"] = float(splitted[4].strip())
        elif line.startswith("TM-score="):
            tm_score = line.split(" ")[1]
            if "Chain_1" in line:
                info_dicts[-1]["tm_score1"] = float(tm_score)
            else:
                info_dicts[-1]["tm_score2"] = float(tm_score)
        elif read_alignment and line.startswith('(":" denotes residue pairs of'):
            next_line_is_alignment = True

    return info_dicts


def call_and_parse(
    pdb_query_path: Path,
    pdb_target_path: Path,
    id1: str,
    id2: str,
    read_alignment: bool = False,
    raise_error: bool = True,
) -> List[Dict]:
    """Run a TM-Align alignment.

    :param pdb_query_path: Path to query PDB file.
    :param pdb_target_path: Path to target PDB file.
    :param id1: Identifier for query.
    :param id2: Identifier for target.
    :param read_alignment: Whether to parse the alignment from the results.
    :param raise_error: Whether to raise exception if TM-Align fails.
    :return: List of dicts. One dict for each pair-wise alignment in the output.
    """
    with tempfile.TemporaryDirectory() as t:
        outdir = Path(t)

        # tmalign can not read gz files
        if pdb_query_path.suffixes[-1] == ".gz":
            unpack_gz(pdb_query_path, outdir / pdb_query_path.stem)
            pdb_query_path = outdir / pdb_query_path.stem
        if pdb_target_path.suffixes[-1] == ".gz":
            unpack_gz(pdb_target_path, outdir / pdb_target_path.stem)
            pdb_target_path = outdir / pdb_target_path.stem

        out_prefix = f"tmalign_stdout_{pdb_query_path.name}_{pdb_target_path.name}"

        # split_flag=2 and ter_flag=0 means that all n^2 monomer chain combinations are tried out
        res_dict = call_tmalign(
            pdb_query_path,
            pdb_target_path,
            outdir / out_prefix,
            out_prefix,
            # split_flag=2,
            # ter_flag=0,
            write_rotation=False,
            raise_error=raise_error,
        )

        if not res_dict or not "stdout" in res_dict or not res_dict["stdout"]:
            raise ValueError("Failed TMalign calculation")

        info_dicts = parse_stdout(
            res_dict["stdout"].decode(), read_alignment=read_alignment
        )

        # TMalign uses the file name as id. The filename might differ from our intern identifiers.
        # So we just overwrite TMaligns identifiers with our own.
        for i_dict in info_dicts:
            i_dict["id1"] = id1
            i_dict["id2"] = id2

        return info_dicts
