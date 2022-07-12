import pandas as pd
from helper.datasets import read_protherm_single_csv, read_prothermdb_single, read_thermomutdb_single_json
from pathlib import Path

PROTHERM_PDB_WILD_COLUMN = 'pdb_wild'
PROTHERM_PDB_MUTANT_COLUMN = 'pdb_mutant'
PROTHERM_MUTATION_COLUMN = 'mutation'
DDG = 'dd_g'

df = pd.read_csv('data/protherm.csv', sep=',', header=0, encoding='latin-1')

rel_cols = [PROTHERM_PDB_WILD_COLUMN, PROTHERM_PDB_MUTANT_COLUMN, PROTHERM_MUTATION_COLUMN, DDG]

_df = df[rel_cols]

# -------------------------------------------------------------------------------
# all ddG mutations in Protherm
nof_all_mutations_ddg = _df.query('dd_g == dd_g and mutation != "wild" and pdb_wild == pdb_wild and mutation == mutation').drop_duplicates(['pdb_wild', 'mutation']).shape[0]
nof_all_mutations_ddg_pdbmutant = _df.query('dd_g == dd_g and mutation != "wild" and pdb_wild == pdb_wild and mutation == mutation and pdb_mutant == pdb_mutant').drop_duplicates(['pdb_wild', 'mutation']).shape[0]
anteil_mit_pdbmutant_ddg = nof_all_mutations_ddg_pdbmutant / nof_all_mutations_ddg

_df2 = _df.query('dd_g == dd_g and mutation != "wild" and pdb_wild == pdb_wild and mutation == mutation').drop_duplicates('pdb_wild').fillna('-')
print()


# -------------------------------------------------------------------------------

def get_protherm_counts():
    PROTHERM_PDB_WILD_COLUMN = 'pdb_wild'
    PROTHERM_PDB_MUTANT_COLUMN = 'pdb_mutant'
    PROTHERM_MUTATION_COLUMN = 'mutation'
    df = pd.read_csv('data/protherm.csv', sep=',', header=0, encoding='latin-1')
    orig_rows = df.shape[0]

    for c in ['pdb_wild', 'mutation', 'pdb_mutant']:
        df[c] = df[c].str.strip()

    # df[df.columns] = df.apply(lambda x: x.str.strip())

    df = df.query('pdb_wild == pdb_wild and mutation == mutation')
    nona_rows = df.shape[0]

    df_single_mutation = df[df[PROTHERM_MUTATION_COLUMN].str.split(' ').str.len() == 3]
    nof_single_mutations = df_single_mutation.shape[0]

    # note that it makes a difference of ~40 rows if pdb_mutant is included for dropping
    # df_u = df_single_mutation.drop_duplicates(subset=['pdb_wild', 'mutation', 'pdb_mutant'])
    df_u = df_single_mutation.drop_duplicates(subset=['pdb_wild', 'mutation'])
    unique_rows = df_u.shape[0]

    df_u_mutant = df_u.query('pdb_mutant == pdb_mutant')
    mutant_rows = df_u_mutant.shape[0]


    # _ = read_protherm_single_csv(Path('data/protherm.csv'))

    # def rm_second_from_first(df1, df2, cols):
    #     return df1.merge(df2, on=cols, how='left', indicator=True).query('_merge == "left_only"').drop(['_merge'],
    #
    #                                                                                                    axis=1)

    # _2 = rm_second_from_first(_, df_u_mutant, ['pdb_wild', 'mutation', 'pdb_mutant'])

    return orig_rows, nona_rows, nof_single_mutations, unique_rows, mutant_rows


def get_thermomutdb_counts():
    THERMOMUTDB_PDB_WILD_JSON = 'PDB_wild'
    THERMOMUTDB_PDB_MUTANT_JSON = 'pdb_mutant'
    THERMOMUTDB_MUTATION_JSON = 'mutation_code'

    import json
    with open('data/thermomutdb.json', 'r') as f:
        df = pd.DataFrame(json.load(f))
    orig_rows = df.shape[0]

    for c in ['PDB_wild', 'mutation_code', 'pdb_mutant']:
        df[c] = df[c].str.strip()

    df = df.query('PDB_wild == PDB_wild and mutation_code == mutation_code')
    nona_rows = df.shape[0]

    SINGLE_MUTATION_REGEX = r'^[A-Z]\d+[A-Z]$'
    df_single_mutation = df[df[THERMOMUTDB_MUTATION_JSON].str.upper().str.contains(SINGLE_MUTATION_REGEX)]
    nof_single_mutations = df_single_mutation.shape[0]

    # note that it makes a difference if pdb_mutant is included for dropping
    # df_u = df_single_mutation.drop_duplicates(subset=['pdb_wild', 'mutation', 'pdb_mutant'])
    df_u = df_single_mutation.drop_duplicates(subset=['PDB_wild', 'mutation_code'])
    unique_rows = df_u.shape[0]

    df_u_mutant = df_u.query('pdb_mutant == pdb_mutant')
    mutant_rows = df_u_mutant.shape[0]

    return orig_rows, nona_rows, nof_single_mutations, unique_rows, mutant_rows


def get_prothermdb_counts():
    PROTHERMDB_PDB_WILD_COL = 'PDB_wild'
    PROTHERMDB_PDB_MUTATION = 'PDB_Chain_Mutation'

    df = pd.read_csv('data/ProThermDB_29_march_2021.tsv', sep='\t')
    orig_rows = df.shape[0]

    for c in ['PDB_wild', 'PDB_Chain_Mutation']:
        df[c] = df[c].str.strip()

    df = df.query('PDB_wild == PDB_wild and PDB_Chain_Mutation == PDB_Chain_Mutation')
    nona_rows = df.shape[0]

    MUTATION_REGEX = r'^(?:[0-9][A-Z0-9]{3}_[A-Z|0-9]:[A-Z]\d+[A-Z]\s?)+$'  # example: 1IRO_A:V24I 1IRO_A:I33L
    df_single_mutation = df[df[PROTHERMDB_PDB_MUTATION].str.upper().str.contains(MUTATION_REGEX)]
    nof_single_mutations = df_single_mutation.shape[0]

    # note that ProThermDB has NO pdb_mutant
    df_u = df_single_mutation.drop_duplicates(subset=['PDB_wild', 'PDB_Chain_Mutation'])
    unique_rows = df_u.shape[0]

    mutant_rows = 0

    return orig_rows, nona_rows, nof_single_mutations, unique_rows, mutant_rows


# count mutations in each data set
print('#rows\tPDB wild rows\tsingle mutations\tunique\tPDB mutant rows')
print('\t'.join(str(_) for _ in get_protherm_counts()))
print('\t'.join(str(_) for _ in get_thermomutdb_counts()))
print('\t'.join(str(_) for _ in get_prothermdb_counts()))




# _ = read_protherm_single_csv(Path('data/protherm.csv'))
# _ = read_prothermdb_single(Path('data/ProThermDB_29_march_2021.tsv'))
_ = read_thermomutdb_single_json(Path('data/thermomutdb.json'))
import helper
one_2_three_dict = {'G': 'GLY', 'A': 'ALA', 'L': 'LEU', 'M': 'MET', 'F': 'PHE',
                    'W': 'TRP', 'K': 'LYS', 'Q': 'GLN', 'E': 'GLU', 'S': 'SER',
                    'P': 'PRO', 'V': 'VAL', 'I': 'ILE', 'C': 'CYS', 'Y': 'TYR',
                    'H': 'HIS', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'T': 'THR'}
_2 = _.drop_duplicates([helper.WILD_COL, helper.WILD_AA, helper.WILD_SEQ_NUM, helper.MUT_AA])
_3 = _2[_2[helper.WILD_AA].isin(one_2_three_dict)]
_4 = _3[_3[helper.MUT_AA].isin(one_2_three_dict)]
_5 = _4[~_4[helper.WILD_COL].isin(helper.BAD_PDBIDS)]

print()

from helper.utils import get_pdb_file_path
get_pdb_file_path('7kgl', True)
