import pandas as pd
import json
from pathlib import Path
from helper.cmdl_calls import exe_cmdl_call
from helper.utils import timer
import logging
logger = logging.getLogger(__name__)

cols = ['query', 'target', 'id', 'align_len', 'mismatches', 'gap_openings', 'query_start',
        'query_end', 'target_start', 'target_end', 'evalue', 'bit_score']
df_pairs = pd.read_csv(Path('/local/sieg/data/scope95/mmseqs_clustering/scope95/align.m8'),
                       sep='\t', header=None, names=cols)
df_pairs = df_pairs.query('query != target')  # drop self hits
print(f'{df_pairs.shape[0]} pairs')

scope95_prefix = Path('/local/sieg/data/scope95/pdbstyle-2.08/')
queries = df_pairs['query']
targets = df_pairs['target']

logging.basicConfig(filename=f'{__name__}.log',
                    level=logging.INFO,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s')

#queries = queries[:2500]
#targets = targets[:2500]

res_dir = Path('test')

res_list = []
for i, (q, t) in enumerate(zip(queries, targets)):
    if i % 10 == 0:
        print(i, q, t)

    this_dir = res_dir / f'{q}_vs_{t}'
    this_dir.mkdir(exist_ok=True)

    args = [
        '/local/sieg/projekte/microminer/bin/MicroMiner_release',
        'site_align',
        '-r', str((scope95_prefix / q[2:4] / f'{q}.ent').resolve()),
        '-s', str((scope95_prefix / t[2:4] / f'{t}.ent').resolve()),
        '-o', str(this_dir.resolve())
    ]

    response = exe_cmdl_call(args, log_msg=f'{q}_vs_{t}', raise_error=True)
    # print(response)

    if response['exit_code'] == 0:
        stdout = response['stdout'].decode('utf8')
        data_dicts = []
        ascona_time = None
        mm_time = None

        for line in stdout.splitlines():
            # print(line)
            if line.startswith('{'):
                # print(line)
                try:
                    data_dicts.append(json.loads(line))
                except Exception as e:
                    print(line)
                    raise e
            elif line.strip().startswith('ascona activeSite Alignment'):
                ascona_time = float(line.split(' ')[-2])
                # print(line)
            elif line.strip().startswith('MicroMiner activeSite Alignment'):
                mm_time = float(line.split(' ')[-2])
                # print(line)
        for d in data_dicts:
            if d['method'] == "ascona":
                d['time'] = ascona_time
            elif d['method'] == 'microminer':
                d['time'] = mm_time
        assert len(data_dicts) <= 2
        res_list.extend(data_dicts)
    else:
        print('FAILED', response)


df_res = pd.DataFrame(res_list)
df_res.to_csv('ascona_vs_microminer.csv', sep='\t')
# print(res_list)
