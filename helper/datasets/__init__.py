"""
The :mod: dataset module for reading mutation data sets.
"""

import helper

SUPPORTED_DATASETS = [
    # mutation data sets
    'protherm', 'prothermdb', 'thermomutdb', 'skempi2', 'platinum', 'shanthirabalan',
    # other datasets
    'pisces',
    'scope',]

MUTATION_DATASETS = [
    'protherm', 'prothermdb', 'thermomutdb', 'skempi2', 'platinum', 'shanthirabalan'
]
MUTATION_DATASETS_WITH_STRUCTURE_PAIRS = ['protherm', 'thermomutdb', 'platinum', 'shanthirabalan']

DATASETS_WITH_CUSTOM_STRUCTUREFILES = ['platinum', 'skempi2', 'scope']

DATABASE_PATH_MAP = {
    'protherm': helper.CONFIG['DATA']['PROTHERM'],
    'prothermdb': helper.CONFIG['DATA']['PROTHERMDB'],
    'thermomutdb': helper.CONFIG['DATA']['THERMOMUTDB'],
    'skempi2': helper.CONFIG['DATA']['SKEMPI2'],
    'platinum': helper.CONFIG['DATA']['PLATINUM'],
    'shanthirabalan': helper.CONFIG['DATA']['SHANTHIRABALAN'],
    'pisces': helper.CONFIG['DATA']['PISCES'],
    'scope': None,  # scope IDs are inferred from the structure files in file system.
}

if len(DATABASE_PATH_MAP) != len(SUPPORTED_DATASETS) or \
        any((d not in DATABASE_PATH_MAP for d in SUPPORTED_DATASETS)):
    raise Exception('Inconsistent data structures! Keep them in sync!')

# import functions for reading single mutations from each data set
from .protherm.read import read_protherm_single
from .thermomutdb.read import read_thermomutdb_single
from .prothermdb.read import read_prothermdb_single
from .skempi2.read import read_skempi2_single
from .platinum.read import read_platinum_single
from .shanthirabalan.read import read_shanthirabalan_single

from . import skempi2
from . import platinum
from . import scope

from .read import read_mutation_dataset_single
