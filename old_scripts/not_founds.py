from helper.datasets import read_thermomutdb_single_json
from helper import CONFIG, WILD_COL
from helper.utils import get_pdb_file_path

df = read_thermomutdb_single_json(CONFIG['DATA']['THERMOMUTDB'])

not_founds = []
for pdbid in df[WILD_COL].unique():
    try:
        get_pdb_file_path(pdbid)
    except FileNotFoundError as e:
        not_founds.append(pdbid)

print('not_founds', not_founds)

