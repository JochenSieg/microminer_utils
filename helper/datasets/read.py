import pandas as pd
from . import SUPPORTED_DATASETS, DATABASE_PATH_MAP
import importlib


def read_mutation_dataset_single(dataset_name: str, pdb_mutant_only: bool) -> pd.DataFrame:
    dataset_name_lower = dataset_name.lower()

    if dataset_name_lower not in SUPPORTED_DATASETS:
        raise Exception('Unsupported dataset: {}'.format(dataset_name))

    for dataset in SUPPORTED_DATASETS:
        if dataset_name_lower == dataset:
            read_func_name = 'read_{}_single'.format(dataset)
            module = importlib.import_module('.', package=__package__)
            read_func = getattr(module, read_func_name)
            return read_func(DATABASE_PATH_MAP[dataset], pdb_mutant_only)

