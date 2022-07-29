import os
from pathlib import Path
import configparser

# read project config from parent dir
CONFIG = configparser.ConfigParser()
config_path = Path(__file__).parent / ".." / "config.ini"
try:
    CONFIG.read(config_path)
except FileNotFoundError as e:
    print("ERROR: config.ini file in parent directory MUST exist")
    raise e
# overwrite configuration with environment variables if set
if len(CONFIG["OVERWRITABLE"]) > 0:
    for _, section in CONFIG["OVERWRITABLE"].items():
        for k, _ in CONFIG[section].items():
            # configparser automatically converts all keys to lower case. This is
            # bad to match against ENV variables. We use the convention that all
            # keys are upper case.
            k_upper = k.upper()
            if k_upper in os.environ:
                CONFIG[section][k_upper] = os.environ[k_upper]

# column names to describe structure (pairs) for single mutations from mutation data sets.
WILD_COL = "wild_pdb"
WILD_CHAIN = "wild_chain"
WILD_AA = "wild_aa"
WILD_SEQ_NUM = "wild_seq_num"
MUTANT_COL = "mut_pdb"
MUTANT_CHAIN = "mut_chain"
MUT_AA = "mut_aa"
MUT_SEQ_NUM = "mut_seq_num"

# MicroMiner result CSV file columns
MM_QUERY_NAME = "queryName"
MM_QUERY_AA = "queryAA"
MM_QUERY_CHAIN = "queryChain"
MM_QUERY_POS = "queryPos"
MM_HIT_NAME = "hitName"
MM_HIT_AA = "hitAA"
MM_HIT_CHAIN = "hitChain"
MM_HIT_POS = "hitPos"
MM_SITE_BACKBONE_RMSD = "siteBackBoneRMSD"
MM_SITE_TMSCORE = "siteTMScore"
MM_SITE_IDENTITY = "siteIdentity"
MM_SITE_GAPS = "siteGaps"
MM_SITE_MISMATCHES = "siteMismatches"
MM_SITE_RESIDUES = "nofSiteResidues"

# 20 standard amino acid mapping from 1- to 3-letter code
one_2_three_dict = {
    "G": "GLY",
    "A": "ALA",
    "L": "LEU",
    "M": "MET",
    "F": "PHE",
    "W": "TRP",
    "K": "LYS",
    "Q": "GLN",
    "E": "GLU",
    "S": "SER",
    "P": "PRO",
    "V": "VAL",
    "I": "ILE",
    "C": "CYS",
    "Y": "TYR",
    "H": "HIS",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "T": "THR",
}

# 20 standard amino acid mapping from 3- to 1-letter code
three_2_one_dict = {v: k for k, v in one_2_three_dict.items()}


# physico-chemical property grouping of amino acids
alipathic = {"ALA", "ILE", "LEU", "MET", "VAL"}
aromatic = {"PHE", "TRP", "TYR"}
polar_neutral = {"ASN", "CYS", "GLN", "SER", "THR"}
polar_acidic = {"ASP", "GLU"}
polar_basic = {"ARG", "HIS", "LYS"}
physicochemical_groups = [alipathic, aromatic, polar_neutral, polar_acidic, polar_basic]

# volume based size grouping of amino acids
volume_grouping = ["very_small", "small", "medium", "large", "very_large"]

# dict for looking up amino acids by size based on volume
volume_grouping_dict = {
    "GLY": 0,
    "ALA": 0,
    "SER": 0,
    "CYS": 1,
    "ASP": 1,
    "PRO": 1,
    "ASN": 1,
    "THR": 1,
    "GLU": 2,
    "VAL": 2,
    "GLN": 2,
    "HIS": 2,
    "MET": 3,
    "ILE": 3,
    "LEU": 3,
    "LYS": 3,
    "ARG": 3,
    "PHE": 4,
    "TYR": 4,
    "TRP": 4,
}

# side chain atom count based size grouping
atom_count_grouping = ["0_to_3_atoms", "4_to_5_atoms", "6_to_10_atoms"]

# dict for looking up amino acids by size based on side chain atom count
atom_count_grouping_dict = {
    "GLY": 0,
    "ALA": 0,
    "SER": 0,
    "CYS": 0,
    "PRO": 0,
    "THR": 0,
    "VAL": 0,
    "ASP": 1,
    "ASN": 1,
    "LEU": 1,
    "ILE": 1,
    "MET": 1,
    "GLU": 1,
    "GLN": 1,
    "LYS": 1,
    "HIS": 2,
    "PHE": 2,
    "ARG": 2,
    "TYR": 2,
    "TRP": 2,
}

# in the mutation data sets PDB-IDs are encountered that we want to skip
BAD_PDBIDS = [
    "1LZ2",
    "1GSB",
    "1LRP",
    "1XAS",
    "1MAR",  # C-alpha only structures
    "3GQ6",  # withdrawn from PDB(e) ?
    "1BGL",  # obsolete: superseded by: 4v40
    "4V40",  # Only mmCIF
    "1SEE",  # obsolete: removed as part of separation of theoretical model coordinate files
    # from the main archive
]
