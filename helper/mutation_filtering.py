import pandas as pd
from .constants import WILD_AA, MUT_AA, one_2_three_dict

alipathic = {'ALA', 'ILE', 'LEU', 'MET', 'VAL'}
aromatic = {'PHE', 'TRP', 'TYR'}
polar_neutral = {'ASN', 'CYS', 'GLN', 'SER', 'THR'}
polar_acidic = {'ASP', 'GLU'}
polar_basic = {'ARG', 'HIS', 'LYS'}
physicochemical_groups = [alipathic, aromatic, polar_neutral, polar_acidic, polar_basic]

# G 	A 	S
# C 	D 	P 	N 	T
# E 	V 	Q 	H
# M 	I 	L 	K 	R
# F 	Y 	W
# Die Einteilung ist irgendeine nach Volumen. Im Detail nicht ganz nachvollziehbar.
very_small = {'GLY', 'ALA', 'SER'}
small = {'CYS', 'ASP', 'PRO', 'ASN', 'THR'}
medium = {'GLU', 'VAL', 'GLN', 'HIS'}
large = {'MET', 'ILE', 'LEU', 'LYS', 'ARG'}
very_large = {'PHE', 'TYR', 'TRP'}

# side chain atom count based size grouping. Thx an Patrick an der Stelle.
few_atoms = ['GLY', 'ALA', 'SER', 'CYS', 'PRO', 'THR', 'VAL']   # 0-3 atoms
medium_atoms = ['ASP', 'ASN', 'LEU', 'ILE', 'MET', 'GLU', 'GLN', 'LYS']  # 4-5 atoms
many_atoms = ['HIS', 'PHE', 'ARG', 'TYR', 'TRP']  # 6-10 atoms
atom_count_size_groups = [few_atoms, medium_atoms, many_atoms]


def are_in_same_group(aa1, aa2, groups):
    for g in groups:
        if aa1 in g:
            return aa2 in g
    return False


def mutations_are_silent(df: pd.DataFrame) -> list:
    """
    Returns a mask for the table rows. True values indicating mutation is silent
    a.k.a. mutant side chain has same properties as wild-type. False values mean
    a change there is a change of property.
    :param df: The table.
    :return: boolean mask indicating silent mutations.
    """

    is_silent_list = []
    for row in df.itertuples(index=False):
        is_silent_list.append(
            are_in_same_group(one_2_three_dict[getattr(row, WILD_AA)],
                              one_2_three_dict[getattr(row, MUT_AA)],
                              physicochemical_groups))
    return is_silent_list


def have_same_size(df: pd.DataFrame) -> list:

    is_diff = []
    for row in df.itertuples(index=False):
        is_diff.append(are_in_same_group(one_2_three_dict[getattr(row, WILD_AA)],
                                         one_2_three_dict[getattr(row, MUT_AA)],
                                         atom_count_size_groups))
        # print('Compared: ', one_2_three_dict[getattr(row, WILD_AA)], one_2_three_dict[getattr(row, MUT_AA)], is_diff[-1])
    return is_diff


def is_proline_mutation(df: pd.DataFrame) -> list:
    """
    Return boolean mask indicating whether row contains a proline involving mutation.
    :param df: The table.
    :return: boolean mask indicating whether wild AA or mutant AA is proline.
    """
    return (df['wildAA'] == 'PRO' or df['mutantAA'] == 'PRO').tolist()
