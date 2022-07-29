"""
The :mod:`hpc` module contains functionality for distribute work on the ZBHs SGE HPC.
"""
from . import distribute_csv
from .distribute_csv import distribute_csv

__all__ = ["distribute_csv"]
