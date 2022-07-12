import pandas as pd
from pathlib import Path
from helper import CONFIG, WILD_COL, MUTANT_COL, MUT_AA, WILD_SEQ_NUM, WILD_AA, BAD_PDBIDS, calc_tmalign
from helper.datasets import read_protherm_single_csv, read_thermomutdb_single_json


def rm_second_from_first(df1, df2, cols):
    return df1.merge(df2, on=cols, how='left', indicator=True).query('_merge == "left_only"').drop(['_merge'], axis=1)


def get_intersec(df1: pd.DataFrame, df2: pd.DataFrame, cols):
    return df1.merge(df2, on=cols, how='left', indicator=True).query('_merge == "both"').drop(['_merge'], axis=1)


def print_wild_and_mutant(df: pd.DataFrame, prefix: str):
    print('{} wild={} mutant={} different PDB ids'.format(
        prefix, df[WILD_COL].unique().shape[0], df[MUTANT_COL].unique().shape[0]))





relevant_cols = [WILD_COL, MUTANT_COL, MUT_AA, WILD_SEQ_NUM, WILD_AA]

df_protherm = read_protherm_single_csv(Path(CONFIG['DATA']['PROTHERM']))
df_thermomutdb = read_thermomutdb_single_json(Path(CONFIG['DATA']['THERMOMUTDB']))

df_problem = pd.read_csv(Path(CONFIG['DATA']['PROBLEMATIC_MUTATIONS']), sep='\t')

# todo das type umwandlen sollte direkt in den read funktionen passieren.
df_protherm[WILD_SEQ_NUM] = df_protherm[WILD_SEQ_NUM].astype(int)
df_thermomutdb[WILD_SEQ_NUM] = df_thermomutdb[WILD_SEQ_NUM].astype(int)
df_problem[WILD_SEQ_NUM] = df_problem[WILD_SEQ_NUM].astype(int)

# drop duplicates to avoid problems
df_protherm.drop_duplicates(subset=relevant_cols, inplace=True)
df_thermomutdb.drop_duplicates(subset=relevant_cols, inplace=True)
df_problem.drop_duplicates(subset=relevant_cols, inplace=True)

df_protherm = df_protherm[~df_protherm[WILD_COL].isin(BAD_PDBIDS)]
df_protherm = df_protherm[~df_protherm[MUTANT_COL].isin(BAD_PDBIDS)]
df_thermomutdb = df_thermomutdb[~df_thermomutdb[WILD_COL].isin(BAD_PDBIDS)]
df_thermomutdb = df_thermomutdb[~df_thermomutdb[MUTANT_COL].isin(BAD_PDBIDS)]

# CASE 1: Just the intersection of the parsed data sets
print('------------- Naively parsed ------------ ')

print_wild_and_mutant(df_protherm, 'Naively parsed: ProTherm')
print_wild_and_mutant(df_thermomutdb, 'Naively parsed: ThermoMutDB')
df_intersection = get_intersec(df_protherm, df_thermomutdb, relevant_cols)
print('Naively parsed: ProTherm and ThermoMutDB have an overlap of {}'.format(df_intersection.shape[0]))

# CASE 2: Intersection of the parsed data sets TM-Score > 0.5
print('------------- TM-Score >= 0.5 ------------ ')

protherm_b4 = df_protherm.shape[0]
thermomutdb_b4 = df_thermomutdb.shape[0]
print('START COMPUTING TM-SCORES. THIS TAKES A MINUTE...')
df_protherm = calc_tmalign(df_protherm)
df_thermomutdb = calc_tmalign(df_thermomutdb)
df_protherm = df_protherm.query('tm_score1 >= 0.5 and tm_score2 >= 0.5')
df_thermomutdb = df_thermomutdb.query('tm_score1 >= 0.5 and tm_score2 >= 0.5')
print('TM-score>=0.5: {} / {} mutations in ProTherm are with TM-score >= 0.5'.format(
    df_protherm.shape[0], protherm_b4))
print('TM-score>=0.5: {} / {} mutations in ThermoMutDB are with TM-score >= 0.5'.format(
    df_thermomutdb.shape[0], thermomutdb_b4))
print_wild_and_mutant(df_protherm, 'TM-score>=0.5: ProTherm')
print_wild_and_mutant(df_thermomutdb, 'TM-score>=0.5: ThermoMutDB')
df_intersection = get_intersec(df_protherm, df_thermomutdb, relevant_cols)
print('TM-score>=0.5: ProTherm and ThermoMutDB have an overlap of {}'.format(
    df_intersection.shape[0]))

# CASE 3: Intersection of the parsed data sets with TM-Score > 0.5 AND manually reviewed problematic cases are removed.
print('------------- NO BAD MUTATIONS ------------ ')
print('   (starting from data of TM-Score >= 0.5)')
protherm_b4 = df_protherm.shape[0]
thermomutdb_b4 = df_thermomutdb.shape[0]
df_protherm = rm_second_from_first(df_protherm, df_problem, relevant_cols)
df_thermomutdb = rm_second_from_first(df_thermomutdb, df_problem, relevant_cols)
print('NO BAD MUTATIONS: {} / {} mutations in ProTherm are NOT problematic (EXCLUDING TM-Align)'.format(
    df_protherm.shape[0], protherm_b4))
print('NO BAD MUTATIONS: {} / {} mutations in ThermoMutDB are NOT problematic (EXCLUDING TM-Align)'.format(
    df_thermomutdb.shape[0], thermomutdb_b4))
print_wild_and_mutant(df_protherm, 'NO BAD MUTATIONS: ProTherm')
print_wild_and_mutant(df_thermomutdb, 'NO BAD MUTATIONS: ThermoMutDB')
df_intersection = get_intersec(df_protherm, df_thermomutdb, relevant_cols)
print('NO BAD MUTATIONS: Without TM-score<0.5 AND problematic cases ProTherm and ThermoMutDB have an overlap of {}'.format(
    df_intersection.shape[0]))


print()
print()
print('=========== UNION ===========')
# pd.concat([df_protherm[relevant_cols], df_thermomutdb[relevant_cols]])
df_union = pd.concat([df_protherm, df_thermomutdb]).drop_duplicates(subset=relevant_cols)

print('UNION: ProTherm and ThermoMutDB together contain {} unique mutations'.format(df_union.shape[0]))
print_wild_and_mutant(df_union, 'UNION: ProTherm and ThermoMutDB')


print()
