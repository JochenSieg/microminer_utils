from abc import ABC, abstractmethod
from typing import List

import pandas as pd


class Dataset(ABC):
    """Abstract class to work with data sets of protein structures from different sources."""

    @property
    @abstractmethod
    def name(self) -> str:
        """Abstract property holding the data sets name.

        :return: The data sets name.
        """
        pass

    @property
    @abstractmethod
    def has_structure_pairs(self) -> bool:
        """Abstract property indicating if the data set contains structure pairs.
        Structure pairs are matched structures indicating a known relationship that
        MicroMiner should find. For example, if MicroMiner searches with a
        wild-type structure it should find the corresponding mutant structure.

        :return: True if the data set contains structure pairs, false otherwise.
        """
        pass

    @property
    @abstractmethod
    def has_custom_structure_files(self) -> bool:
        """Abstract property that indicates whether the data set comes with customly
        prepared structure files (not the standard PDB files).

        :return: True if the data set comes with custom structure files, false otherwise.
        """
        pass

    @abstractmethod
    def read(self) -> pd.DataFrame:
        """Abstract function to read the raw data set into standardized format.

        :return: A table containing the data set.
        """
        pass


class MutationDataset(Dataset):
    """Abstract class to represent a mutation dataset."""

    @abstractmethod
    def read_single_mutations(self, pdb_mutant_only: bool) -> pd.DataFrame:
        """Abstract function to extract all single mutations from the mutation data set.

        :param pdb_mutant_only: Whether to read only mutations with PDB structure
               for the mutant.
        :return: A table containing the single mutations.
        """
        pass


class DatasetCollection:
    """Class to hold all datasets available through this module."""

    def __init__(self):
        """Create a new instance."""
        self.dataset_map = {}

    def register_dataset(self, dataset_obj: Dataset) -> None:
        """Register a new data set in the collection.

        :param dataset_obj: The data set object to register.
        :return: None
        """
        self.dataset_map[dataset_obj.name] = dataset_obj

    def contains(self, database_name: str) -> bool:
        """Check if ceratin data set is contained in collection.

        :param database_name: The data sets name.
        :return: True if data set is contained, false otherwise.
        """
        return database_name in self.dataset_map

    @staticmethod
    def is_mutation_dataset(database_obj: Dataset) -> bool:
        """Checks if a data set instance is a mutation data set.

        :param database_obj: The data set to check.
        :return: True if the data set is a mutation data set, false otherwise.
        """
        return isinstance(database_obj, MutationDataset)

    def get_dataset(self, database_name: str) -> Dataset:
        """Get data set instance by name.

        :param database_name: The data sets name.
        :return: The data set instance.
        """
        return self.dataset_map[database_name]

    def get_dataset_names(self) -> List[str]:
        """Get names of all data set in this collection.

        :return: List of data set names.
        """
        return list(self.dataset_map.keys())

    def get_mutation_datasets_with_structure_pairs(self) -> List[Dataset]:
        """Get all mutation data sets in the collection.

        :return: List of mutation data set instances.
        """
        return [
            dataset
            for dataset in self.dataset_map.values()
            if dataset.has_structure_pairs
            and DatasetCollection.is_mutation_dataset(dataset)
        ]

    def __iter__(self):
        """iter function to iterate contained data sets.

        :return: Iterator of the contained data sets.
        """
        return self.dataset_map.values().__iter__()


class MockDataset(Dataset):
    """Mock dataset class for testing"""

    name = "mock_dataset"
    has_custom_structure_files = False
    has_structure_pairs = True

    def read(self) -> pd.DataFrame:
        """Mock read function. Will be overwritten in tests.

        :return: Nothing here.
        """
        pass


class MockMutationDataset(MutationDataset):
    """Mock dataset class for testing"""

    name = "mock_mutation_dataset"
    has_custom_structure_files = False
    has_structure_pairs = True

    def read(self) -> pd.DataFrame:
        """Mock read function. Will be overwritten in tests.

        :return: Nothing here.
        """
        pass

    def read_single_mutations(self, pdb_mutant_only: bool) -> pd.DataFrame:
        """Mock read function for single mutations. Will be overwritten in tests.

        :param pdb_mutant_only: Whether to read only mutations with PDB structure
               for the mutant.
        :return: Nothing here.
        """
        pass
