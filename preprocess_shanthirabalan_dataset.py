"""Preprocesses the raw data from
Shanthirabalan S, Chomilier J, Carpentier M. Structural effects of point mutations in proteins.
Proteins. 2018 Aug;86(8):853-867. doi: 10.1002/prot.25499.
Thanks to Mathilde Carpentier for providing the raw data to us.

This script reads in pairs of PDB-IDs with chain identifiers. These pairs should only differ
at a single position by a point mutation. However, also some differences at the termini are
possible.
This script finds the positions of these point mutations by:
  1. aligning the PDB-files chain-wise with TMalign and extracting the mutation position
     of the alignment.
  2. mapping the position of the mutation in the alignment to the PDB-file numbering
"""
import gzip
import logging
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

import helper
from helper.constants import one_2_three_dict
from helper.tmalign import call_tmalign, parse_stdout
from helper.utils import get_pdb_file_path

logger = logging.getLogger(__name__)


def make_mutation_pair_df():
    # protein families with mutation. Taken from "LesFamilles.xlsx"
    df_families = pd.DataFrame([
        {'wild_name': '1lw9', 'cluster_id': 31255, 'wild_chain': 'A'},
        {'wild_name': '2nwd', 'cluster_id': 37522, 'wild_chain': 'X'},
        {'wild_name': '2dek', 'cluster_id': 18272, 'wild_chain': 'A'},
        {'wild_name': '2ili', 'cluster_id': 18267, 'wild_chain': 'A'},
        {'wild_name': '1ey0', 'cluster_id': 34381, 'wild_chain': 'A'},
        {'wild_name': '4bfl', 'cluster_id': 796, 'wild_chain': 'A'},
        {'wild_name': '2e3w', 'cluster_id': 38031, 'wild_chain': 'A'},
        {'wild_name': '2vb1', 'cluster_id': 37731, 'wild_chain': 'A'},
        {'wild_name': '4fi8', 'cluster_id': 37628, 'wild_chain': 'A'},
        {'wild_name': '2j8c', 'cluster_id': 13574, 'wild_chain': 'M'},
        {'wild_name': '5dei', 'cluster_id': 2739, 'wild_chain': 'A'},
    ])

    df_mutants = pd.read_csv('data/RMS_WholeProt_Mutantes.txt',
                             sep='\s+', names=['cluster_id', 'pdbid_chainid', 'rms'])
    df_mutants[['mutant_name', 'mutant_chain']] \
        = df_mutants['pdbid_chainid'].str.split('_', expand=True)
    df_mutants.drop('pdbid_chainid', axis=1, inplace=True)

    return df_families.merge(df_mutants, on=['cluster_id'], how='right')


def unpack_gz(gz_file, out_file):
    with gzip.open(gz_file, 'rb') as f_in:
        with open(out_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def get_atom_based_sequence(pdb_file):
    seq = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                seqnum = line[22:27]  # last char is iCode
                if len(seq) == 0 or seqnum != seq[-1]['seqnum']:
                    seq.append(
                        {'residue': line[17:20],
                         'chain': line[21:22],
                         'seqnum': seqnum,
                         })
    return seq


def calc_point_mutation_position(info_dict, ref_seq, tar_seq):
    if len(info_dict['refseq']) != len(info_dict['targetseq']):
        raise ValueError('Not matching alignments')
    not_matching_pos_list = []
    for i in range(len(info_dict['refseq'])):
        if info_dict['refseq'][i] != info_dict['targetseq'][i]:
            not_matching_pos_list.append(i)

    # print(not_matching_pos)
    diffs = [min(i, len(info_dict['refseq']) - i) for i in not_matching_pos_list]
    not_matching_pos = not_matching_pos_list[np.argmax(diffs)]

    # get pos in sequence
    gap_count_ref = 0
    gap_count_tar = 0
    for i in range(not_matching_pos):
        if info_dict['refseq'][i] == '-':
            gap_count_ref += 1
        if info_dict['targetseq'][i] == '-':
            gap_count_tar += 1

    ref_align_pos = not_matching_pos - gap_count_ref
    ref_aa = info_dict['refseq'][not_matching_pos]
    tar_align_pos = not_matching_pos - gap_count_tar
    tar_aa = info_dict['targetseq'][not_matching_pos]

    ref_infile_id = None
    tar_infile_id = None
    if ref_aa != '-' and tar_aa != '-':
        ref_infile_aa = ref_seq[ref_align_pos]['residue']
        tar_infile_aa = tar_seq[tar_align_pos]['residue']

        ref_aa3 = one_2_three_dict[ref_aa]
        tar_aa3 = one_2_three_dict[tar_aa]

        # print(ref_infile_aa,ref_aa3,  tar_infile_aa, tar_aa3)
        if ref_infile_aa != ref_aa3 or tar_infile_aa != tar_aa3:
            # es geht nicht in dem beispiel von 4fi8 mit ein paar Strukturen:
            #   sieht so aus als gibt es mehr als nur eine Mutation und ein 3mer SGS in 4fi8 wird
            #   als deletion algniert S-S. Dadurch stimmt die Rechnung mit den Indices oben nicht
            #   mehr.
            # An sich sollte es eigentlich korrekt sein, diese indel Faelle zu ignorieren.
            # raise ValueError('geht leider nicht')
            print('Warning: Ignoring alignment with indel')
        ref_infile_id = ref_seq[ref_align_pos]['seqnum']
        tar_infile_id = tar_seq[tar_align_pos]['seqnum']

    return ref_align_pos, ref_infile_id, ref_aa, tar_align_pos, tar_infile_id, tar_aa


def main():
    # get the data from
    df = make_mutation_pair_df()

    wild_aas = []
    wild_align_pos = []
    wild_infile_pos = []
    mutant_aas = []
    mutant_align_pos = []
    mutant_infile_pos = []
    tmalign_rmsd_list = []
    tmalign_seqid_list = []
    tmscore1_list = []
    tmscore2_list = []
    tmalign_seqlen1_list = []
    tmalign_seqlen2_list = []
    tmalign_alignlen_list = []
    for wild_name, wild_chain, mutant_name, mutant_chain in zip(df['wild_name'],
                                                                df['wild_chain'],
                                                                df['mutant_name'],
                                                                df['mutant_chain']):
        with tempfile.TemporaryDirectory() as tmpdirname:
            tmpdir_path = Path(tmpdirname)
            pdb_query_path = get_pdb_file_path(wild_name)
            pdb_target_path = get_pdb_file_path(mutant_name)

            # tmalign can not read gz files
            unpack_gz(pdb_query_path, tmpdir_path / pdb_query_path.stem)
            unpack_gz(pdb_target_path, tmpdir_path / pdb_target_path.stem)
            pdb_query_path = tmpdir_path / pdb_query_path.stem
            pdb_target_path = tmpdir_path / pdb_target_path.stem

            out_prefix = f'tmalign_stdout_{pdb_query_path.name}_{pdb_target_path.name}'
            stdout = call_tmalign(pdb_query_path, pdb_target_path, tmpdir_path / out_prefix,
                                  out_prefix,
                                  split_flag=2,
                                  ter_flag=0,
                                  write_rotation=False, raise_error=True)

            if not stdout:
                raise ValueError('Failed TMalign calculation')

            info_dicts = parse_stdout(stdout.decode(), read_alignment=True)

            relevant_info_dict = None
            for i_dict in info_dicts:
                if 'chain_id1' not in i_dict or 'chain_id2' not in i_dict:
                    continue
                if i_dict['chain_id1'] == wild_chain and i_dict['chain_id2'] == mutant_chain:
                    relevant_info_dict = i_dict
                    break
            if relevant_info_dict is None:
                raise ValueError('No alignment for chain pair')

            # print(out_prefix, relevant_info_dict)

            # get sequences from atom entries in PDB files
            ref_seq = [res_entry for res_entry in get_atom_based_sequence(pdb_query_path)
                       if res_entry['chain'] == wild_chain]
            tar_seq = [res_entry for res_entry in get_atom_based_sequence(pdb_target_path)
                       if res_entry['chain'] == mutant_chain]

            ref_align_pos, ref_infile_id, ref_aa, tar_align_pos, tar_infile_id, tar_aa = \
                calc_point_mutation_position(
                    relevant_info_dict, ref_seq, tar_seq)

            wild_aas.append(ref_aa)
            wild_align_pos.append(ref_align_pos)
            wild_infile_pos.append(ref_infile_id)
            mutant_aas.append(tar_aa)
            mutant_align_pos.append(tar_align_pos)
            mutant_infile_pos.append(tar_infile_id)

            tmalign_rmsd_list.append(i_dict['rmsd'])
            tmalign_seqid_list.append(i_dict['seqid'])
            tmscore1_list.append(i_dict['tm_score1'])
            tmscore2_list.append(i_dict['tm_score2'])
            tmalign_seqlen1_list.append(i_dict['chain1_len'])
            tmalign_seqlen2_list.append(i_dict['chain2_len'])
            tmalign_alignlen_list.append(i_dict['align_len'])

            # print(ref_align_pos, ref_infile_id, ref_aa, tar_align_pos, tar_infile_id, tar_aa)

    df[helper.WILD_AA] = wild_aas
    df[helper.WILD_SEQ_NUM] = wild_infile_pos
    df[helper.MUT_AA] = mutant_aas
    df[helper.MUT_SEQ_NUM] = mutant_infile_pos

    df['tmalign_rmsd'] = tmalign_rmsd_list
    df['tmalign_seqid'] = tmalign_seqid_list
    df['tmalign_tmscore1'] = tmscore1_list
    df['tmalign_tmscore2'] = tmscore2_list
    df['tmalign_seqlen1'] = tmalign_seqlen1_list
    df['tmalign_seqlen2'] = tmalign_seqlen2_list
    df['tmalign_alignlen'] = tmalign_alignlen_list

    df.drop('cluster_id', axis=1, inplace=True)

    df.rename(columns={'wild_name': helper.WILD_COL,
                       'wild_chain': helper.WILD_CHAIN,
                       'mutant_name': helper.MUTANT_COL,
                       'mutant_chain': helper.MUTANT_CHAIN,
                       }
              , inplace=True)

    df.to_csv('data/shanathirabalan.csv', sep='\t', index=False)


if __name__ == "__main__":
    main()
