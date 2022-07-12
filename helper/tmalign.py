"""Utils functionality for TM-align like reading and parsing results"""
from pathlib import Path
import os
import logging
from .cmdl_calls import call_tmalign
from .utils import scantree

logger = logging.getLogger(__name__)


def read_stdout_gen(dir_path: Path, prefix: str, suffix: str = ''):
    """
    Search for TMalign output files in directory structure. I called the files just 'stdout'.
    """
    for p in scantree(dir_path):
        if p.is_file():
            if p.name.startswith(prefix) and p.name.endswith(suffix):
                yield p


def parse_stdout_file(path_obj, read_alignment=False):
    with open(path_obj, 'r') as f:
        return parse_stdout(f.read())


def parse_stdout(stdout_str, read_alignment=False):
    info_dicts = [{}]
    next_line_is_alignment = False
    align_counter = 0
    for i, line in enumerate(stdout_str.split(os.linesep)):

        if i == 2 and not line.startswith(' * TM-align'):
            return None
        if i > 2 and line.startswith(' * TM-align'):
            # new entry
            next_line_is_alignment = False
            align_counter = 0
            info_dicts.append({})

        if next_line_is_alignment:
            if align_counter == 0:
                info_dicts[-1]['refseq'] = line.strip('\n')
            elif align_counter == 1:
                info_dicts[-1]['matching'] = line.strip('\n')
            elif align_counter == 2:
                info_dicts[-1]['targetseq'] = line.strip('\n')
            else:
                next_line_is_alignment = False
            align_counter += 1

        elif line.startswith("Name of Chain_1"):
            p = Path(line.split(' ')[3])
            info_dicts[-1]['id1'] = p.stem
            name_split = p.name.split(':')
            if len(name_split) > 1  and len(name_split[-1]) < 3:
                info_dicts[-1]['chain_id1'] = name_split[-1]
        elif line.startswith("Name of Chain_2"):
            p = Path(line.split(' ')[3])
            info_dicts[-1]['id2'] = p.stem
            name_split = p.name.split(':')
            if len(name_split) > 1 and len(name_split[-1]) < 3:
                info_dicts[-1]['chain_id2'] = name_split[-1]
        elif line.startswith("Length of Chain_1"):
            info_dicts[-1]['chain1_len'] = int(line.split(' ')[3])
        elif line.startswith("Length of Chain_2"):
            info_dicts[-1]['chain2_len'] = int(line.split(' ')[3])
        elif line.startswith('Aligned length='):
            splitted = line.split('=')
            info_dicts[-1]['align_len'] = int(splitted[1].strip().split(',')[0])
            info_dicts[-1]['rmsd'] = float(splitted[2].strip().split(',')[0])
            info_dicts[-1]['seqid'] = float(splitted[4].strip())
        elif line.startswith('TM-score='):
            tm_score = line.split(' ')[1]
            if 'Chain_1' in line:
                info_dicts[-1]['tm_score1'] = float(tm_score)
            else:
                info_dicts[-1]['tm_score2'] = float(tm_score)
        elif read_alignment and line.startswith('(":" denotes residue pairs of'):
            next_line_is_alignment = True

    return info_dicts


def call_and_parse(pdb_query_path: Path, pdb_target_path: Path, outdir: Path, out_prefix: str,
                 raise_error: bool = True):
    out_prefix = f'tmalign_stdout_{pdb_query_path.name}_{pdb_target_path.name}'
    stdout = call_tmalign(pdb_query_path, pdb_target_path, outdir, out_prefix, raise_error)

    info_dicts = None
    if stdout:
        info_dicts = parse_stdout(stdout.decode())
    # print('info_dict', info_dict)
    return info_dicts

# def calc_and_collect(pdb_query_path: Path, pdb_target_path: Path):
#     out_prefix = f'tmalign_stdout_{pdb_query_path.name}_{pdb_target_path.name}'
#     res_list = []
#     with tempfile.TemporaryDirectory() as tmpdirname:
#
#         tmalign_stdout = call_tmalign(pdb_query_path, pdb_target_path, Path(tmpdirname), out_prefix)
#
#         if tmalign_stdout is None:
#             return []
#
#         res_data = parse_stdout(tmalign_stdout, False)
#
#         for p in read_stdout_gen(tmpdirname, out_prefix):
#             res_data = parse_stdout_file(p, False)
#             if not res_data:
#                 logging.warning(
#                     'WARNING: Not the TM-align expected result for: {} and {}'.format(name1, name2))
#                 # raise Exception('Not the TM-align expected result')
#             res_list.append(res_data)
#     return res_list


# def calc_tmalign(df_in: pd.DataFrame, cpus=1, use_monomers=False):
#     # generate a list of parameter tuple containing the PDB pairs to be aligned.
#     name_pairs = [(getattr(row, WILD_COL), getattr(row, MUTANT_COL), use_monomers)
#                   for row in
#                   df_in[[WILD_COL, MUTANT_COL]].drop_duplicates().itertuples(index=False)]
#
#     # parallel execution with each parameter set separately
#     with multiprocessing.Pool(cpus) as pool:
#         res_list = list(itertools.chain.from_iterable(pool.starmap(calc_func, name_pairs)))
#
#     # res_list = []
#     # # make protein pairs unique. Duplicates not needed for TMalign calculation.
#     # for idx, row in df_in.drop_duplicates(subset=[WILD_COL, MUTANT_COL]).iterrows():
#     #
#     #     with tempfile.TemporaryDirectory() as tmpdirname:
#     #         call_tmalign(row[WILD_COL], row[MUTANT_COL], Path(tmpdirname))
#     #
#     #         for p in read_stdout_gen(tmpdirname, 'tmalign_stdout_', ''):
#     #             res_data = parse_stdout_file(p, False)
#     #             if not res_data:
#     #                 print('WARNING: Not the TM-align expected result for: {} and {}'.format(
#     #                     row[WILD_COL], row[MUTANT_COL]))
#     #                 # raise Exception('Not the TM-align expected result')
#     #             res_list.append(res_data)
#
#     df = pd.DataFrame(res_list)
#     df['name1'] = df['name1'].str[:4].str.upper() + df['name1'].str[4:]  # first 4 chars to upper
#     df['name2'] = df['name2'].str[:4].str.upper() + df['name2'].str[4:]  # first 4 chars to upper
#     df = df_in.merge(df, left_on=[WILD_COL, MUTANT_COL], right_on=['name2', 'name1'], how='left')
#     # TODO diese assert fliegt, wenn ich alle (chain) pairs in den mutscreen results auf der ThermoMutDB als input nehme?
#     #      sind da doppelte im df drin, dass er mehr als eine Zeile fuer jede Zeile im Input generiert? Das sollte eigentlich nicht sein,
#     #      weil ich die Berechnungen ja nur mit einem echten subset des df_in starte. Da ist ja alles unique... mach ich irgendwann mal upper
#     #      oder lower wodurch dann zwei eigentlich unterschiedliche namen gleich heissen? ODer laueft generell was anderes schief??
#     if df.shape[0] != df_in.shape[0]:
#         print('DAS ASSERT FLIEGT IMMER NOCH!! Was ist da los?')
#     assert (df.shape[0] == df_in.shape[0])
#     df.drop('name1', axis=1, inplace=True)
#     df.drop('name2', axis=1, inplace=True)
#
#     return df
