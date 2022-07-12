import os
from pathlib import Path
import configparser

# read project config from parent dir
CONFIG = configparser.ConfigParser()
config_path = Path(__file__).parent / '..' / 'config.ini'
try:
    CONFIG.read(config_path)
except FileNotFoundError as e:
    print('config.ini file in parent directory MUST exist')
    raise e
# overwrite configuration with environment variables if set
if len(CONFIG['OVERWRITABLE']) > 0:
    for _, section in CONFIG['OVERWRITABLE'].items():
        for k, _ in CONFIG[section].items():
            # configparser automatically converts all keys to lower case. This is
            # nasty to match against ENV variables. We use the convention that all
            # keys are upper case.
            k_upper = k.upper()
            if k_upper in os.environ:
                CONFIG[section][k_upper] = os.environ[k_upper]

WILD_COL = 'wild_pdb'
WILD_CHAIN = 'wild_chain'
WILD_AA = 'wild_aa'
WILD_SEQ_NUM = 'wild_seq_num'
MUTANT_COL = 'mut_pdb'
MUTANT_CHAIN = 'mut_chain'
MUT_AA = 'mut_aa'
MUT_SEQ_NUM = 'mut_seq_num'

MICROMINER_QUERY_NAME = 'queryName'
MICROMINER_QUERY_AA = 'queryAA'
MICROMINER_QUERY_CHAIN = 'queryChain'
MICROMINER_QUERY_POS = 'queryPos'
MICROMINER_HIT_NAME = 'hitName'
MICROMINER_HIT_AA = 'hitAA'
MICROMINER_HIT_CHAIN = 'hitChain'
MICROMINER_HIT_POS = 'hitPos'
MICROMINER_SITE_BACKBONE_RMSD = 'siteBackBoneRMSD'
MICROMINER_SITE_TMSCORE = 'siteTMScore'
MICROMINER_SITE_IDENTITY = 'siteIdentity'
MICROMINER_SITE_GAPS = 'siteGaps'
MICROMINER_SITE_MISMATCHES = 'siteMismatches'
MICROMINER_SITE_RESIDUES = 'nofSiteResidues'

one_2_three_dict = {'G': 'GLY', 'A': 'ALA', 'L': 'LEU', 'M': 'MET', 'F': 'PHE',
                    'W': 'TRP', 'K': 'LYS', 'Q': 'GLN', 'E': 'GLU', 'S': 'SER',
                    'P': 'PRO', 'V': 'VAL', 'I': 'ILE', 'C': 'CYS', 'Y': 'TYR',
                    'H': 'HIS', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'T': 'THR'}

BAD_PDBIDS = [
    '1LZ2', '1GSB', '1LRP', '1XAS', '1MAR',  # C-alpha only structure
    '3GQ6',  # withdrawn from PDB(e) ?
    '1BGL',  # obsolete: superseded by: 4v40  #TODO das koennte wieder mit rein
    '1SEE',  # obsolete: removed as part of separation of theoretical model coordinate files from the main archive
    '4V40',  # Only mmCIF... TODO fix reading mmCIF
]
