import pandas as pd
import gzip
from pathlib import Path
from .constants import WILD_COL, WILD_AA, WILD_SEQ_NUM, one_2_three_dict
from .utils import get_pdb_file_path


def wild_type_exists(df: pd.DataFrame) -> list:

    was_found = [False for _ in range(df.shape[0])]

    for name, group_df in df.groupby([WILD_COL]):

        data = {}
        for row in df.itertuples(index=True):
            data[str(getattr(row, WILD_SEQ_NUM))] = {'AA': one_2_three_dict[getattr(row, WILD_AA)], 'idx': row.Index}

        with gzip.open(get_pdb_file_path(name), 'rt') as f:
            for line in f:
                if not line.startswith('ATOM'):
                    continue

                pos = line[23:27].strip()
                if pos in data and data[pos]['AA'] == line[18:21]:
                    row_pos_in_table = df.index.get_loc(data[pos]['idx'])
                    assert (type(row_pos_in_table) == int)  # there is something wrong if this is not int
                    was_found[row_pos_in_table] = True

