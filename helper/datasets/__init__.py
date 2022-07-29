"""
The :mod: dataset module for reading structure and mutation data sets.
"""
from pathlib import Path

import helper

from .afdb import AFDB
from .dataset import DatasetCollection
from .fireprotdb import FireProtDB
from .pdb import PDB
from .platinum import Platinum
from .protherm import ProTherm
from .prothermdb import ProThermDB
from .shanthirabalan import Shanthirabalan
from .skempi2 import SKEMPI2
from .thermomutdb import ThermoMutDB

# register new datasets here
dataset_collection = DatasetCollection()
dataset_collection.register_dataset(ProTherm(Path(helper.CONFIG["DATA"]["PROTHERM"])))
dataset_collection.register_dataset(Platinum(Path(helper.CONFIG["DATA"]["PLATINUM"])))
dataset_collection.register_dataset(
    ThermoMutDB(Path(helper.CONFIG["DATA"]["THERMOMUTDB"]))
)
dataset_collection.register_dataset(
    ProThermDB(Path(helper.CONFIG["DATA"]["PROTHERMDB"]))
)
dataset_collection.register_dataset(SKEMPI2(Path(helper.CONFIG["DATA"]["SKEMPI2"])))
dataset_collection.register_dataset(
    Shanthirabalan(Path(helper.CONFIG["DATA"]["SHANTHIRABALAN"]))
)
dataset_collection.register_dataset(
    FireProtDB(Path(helper.CONFIG["DATA"]["FIREPROTDB"]))
)
dataset_collection.register_dataset(PDB(Path(helper.CONFIG["DATA"]["PDB_DIR"])))
dataset_collection.register_dataset(AFDB(Path(helper.CONFIG["DATA"]["AFDB_DIR"])))


def get_dataset_collection():
    """Get the module DatasetCollection instance providing all datasets available.

    :return: The dataset collection.
    """
    return dataset_collection


DATASETS_WITH_CUSTOM_STRUCTUREFILES = ["skempi2", "platinum", "scope"]
