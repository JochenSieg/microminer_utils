import pandas as pd
from pathlib import Path
import logging

from helper import CONFIG
from helper.utils import scantree

logger = logging.getLogger(__name__)


def read_scope_classification_file() -> pd.DataFrame:
    """Read in the SCOPe mapping file with all metadata to the contained structures.

    :return: Table with SCOPe metadata.
    """

    mapping_file = Path(CONFIG["DATA"]["SCOPE_MAPPING"])

    # Mapping file HEADER description from SCOPe website https://scop.berkeley.edu:
    # sid
    # PDB ID
    # description
    # sccs
    # sunid
    # sunids of ancestor nodes in comma-delimited list, in the format "level=sunid". For level
    #          codes, see description for dir.des.scop(e).txt.

    # Description of IDs from SCOPe website:
    # sccs. SCOP(e) concise classification string. This is a dot notation used to concisely describe
    # a SCOP(e) class, fold, superfamily, and family. For example, a.39.1.1 references the
    # "Calbindin D9K" family, where "a" represents the class, "39" represents the fold, "1"
    # represents the superfamily, and the last "1" represents the family.
    #
    # sunid. SCOP(e) unique identifier. This is simply a number that may be used to reference any
    # entry in the SCOP(e) hierarchy, from root to leaves (Fold, Superfamily, Family, etc.).
    #
    # sid. Stable domain identifier. A 7-character sid consists of "d" followed by the 4-character
    # PDB ID of the file of origin, the PDB chain ID ('_' if none, '.' if multiple as is the case
    # in genetic domains), and a single character (usually an integer) if needed to specify the
    # domain uniquely ('_' if not). Sids are currently all lower case, even when the chain letter
    # is upper case. Example sids include d4akea1, d9hvpa_, and d1cph.1.

    columns = ["sid", "pdbid", "description", "sccs", "sunid", "sunids_ancestors"]
    usecols = [columns.index(_) for _ in ["sid", "pdbid", "sccs"]]

    # A line in the mapping file looks like this:
    # d1ux8a_ 1ux8    A:      a.1.1.1 113449  cl=46456,cf=46457,sf=46458,fa=46459,dm=46460,sp=116748,px=113449
    df = pd.read_csv(
        mapping_file,
        sep="\s+",
        comment="#",
        header=None,
        names=columns,
        usecols=usecols,
    )
    return df


def read_scope() -> pd.DataFrame:
    """Read SCOPe structure file paths, IDs and classification to a table.

    :return: Table of SCOPe data sets data.
    """
    log_prefix = "SCOPe: "

    data_dir = Path(CONFIG["DATA"]["SCOPE_DATA_DIR"])

    # path should be a dir from which we scan the files.
    if not data_dir.is_dir():
        raise ValueError(
            "SCOPe data dir must be a directory of PDB-style files."
            " Download from SCOPe website."
        )

    df = read_scope_classification_file()

    pdb_files = [p for p in scantree(data_dir) if p.suffixes[-1] == ".ent"]
    sid_set = set([p.stem for p in pdb_files])
    if len(pdb_files) != len(sid_set):
        raise ValueError("SCOPEe PDBstyle dir contains duplicate PDB files")
    df = df[df["sid"].isin(sid_set)]
    assert len(pdb_files) == df.shape[0]

    logger.info(
        log_prefix + f"read {df.shape[0]} domains, which are from"
        f' {len(df["pdbid"].unique())} PDB IDs.'
    )
    return df


def read_scope_pdbstyle(sid: str) -> Path:
    """Get the SCOPe structure file path from SCOPes sid.

    :param sid: The sid.
    :return: Path to structure file.
    """
    pdb_dir = Path(CONFIG["DATA"]["SCOPE_DATA_DIR"])
    # filename is sid + .ent. Looks like this: d5zzwa_.ent
    filename = "{}.ent".format(sid)
    pdbid = sid[1:5]
    path = pdb_dir / pdbid[1:3].lower() / filename
    if not path.is_file():
        raise FileNotFoundError("Could not find SCOPe PDB file: {}".format(path))
    return path
