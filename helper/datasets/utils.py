from pathlib import Path

from helper import CONFIG
from helper.datasets import SKEMPI2, Platinum, scope


def get_divided_dir_name(pdbid: str) -> str:
    """Get the directory name of this PDB-IDs PDB file in the PDBs "divided" file structure.

    :param pdbid: The PDB-ID
    :return: The directory name in the "divided" directory structure.
    """
    if len(pdbid) == 4 or len(pdbid) == 6 and pdbid[4] == "_":
        return pdbid[1:3].lower()
    raise ValueError("Do not recognize format of pdbid: {}".format(pdbid))


def get_pdb_file_path(
    pdbid: str, allow_obsolete: bool = True, mirror: str = "standard"
) -> Path:
    """Get filepath for PDB-ID.

    Tries different PDB file sources to provide a path to the requested PDB file.
    :param pdbid: The PDB-ID.
    :param allow_obsolete: Whether to fall back to obsolete structures if necessary.
    :param mirror: PDB mirror to use. Default is 'standard'. If non-applicable default is used.
    :return: The PDB file path for PDB-ID.
    """
    if mirror == "skempi2":
        return SKEMPI2.read_skempi_cleaned_pdb(pdbid)
    elif mirror == "platinum":
        return Platinum.read_platinum_cleaned_pdb(pdbid)
    elif mirror == "scope":
        return scope.read_scope_pdbstyle(pdbid)
    return get_pdb_file_path_from_standard(pdbid, allow_obsolete)


def get_pdb_file_path_from_standard(pdbid: str, allow_obsolete: bool = True) -> Path:
    """Gets a file path to provided PDB-ID if a corresponding file exists.

    :param pdbid: the PDB-ID
    :param allow_obsolete: allow obsolete structures.
    :return: path to PDB file
    """
    pdb_dir = Path(CONFIG["DATA"]["PDB_DIR"])
    filename = "{}{}{}".format(
        CONFIG["PDB_FILE_INFO"]["PREFIX"],
        pdbid.lower() if CONFIG["PDB_FILE_INFO"]["CASE"] == "LOWER" else pdbid.upper(),
        CONFIG["PDB_FILE_INFO"]["SUFFIX"],
    )
    path = pdb_dir / filename
    if not path.is_file():
        if allow_obsolete and CONFIG["DATA"]["PDB_OBSOLETE_DIR"]:
            path = (
                Path(CONFIG["DATA"]["PDB_OBSOLETE_DIR"])
                / get_divided_dir_name(pdbid)
                / filename
            )
            if path.is_file():
                return path
            raise FileNotFoundError("Could not find PDB file: {}".format(path))
    return path
