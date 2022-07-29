import gzip
import io
from pathlib import Path

import pandas as pd


def read_plddt_values(path: Path) -> pd.DataFrame:
    """Reads pLDDT values for each residue in the given PDB file.

    :param path: Path to alphafold PDB file.
    :return: dataframe of residue-wise pLDDT info.
    """
    if len(path.suffixes) > 0 and path.suffixes[-1] == ".gz":
        with gzip.open(path, "rt") as f_in:
            file_str = f_in.read()
    else:
        with open(path, "r") as f_in:
            file_str = f_in.read()
    pdb_fh = io.StringIO(file_str)

    name = path.stem.split(".")[0]

    # a ATOM line in PDB format looks like this:
    # ATOM      3  C   MET A   1     -33.337   4.131 -10.180  1.00 40.18           C
    residue_names = []
    chain_names = []
    residue_ids = []
    plddt_list = []

    # manually parse the information from the ATOM records.
    last_res_id = None
    for line in pdb_fh:
        current_res_id = line[23:26]
        if line.startswith("ATOM") and current_res_id != last_res_id:
            residue_names.append(line[17:20])
            chain_names.append(line[21:22])
            residue_ids.append(current_res_id.strip())
            plddt_list.append(line[61:66])
            last_res_id = current_res_id

    return pd.DataFrame(
        {
            "name": [name] * len(residue_names),
            "aa": residue_names,
            "chain": chain_names,
            "pos": residue_ids,
            "plddt": plddt_list,
        }
    )
