import contextlib
import gzip
import logging
import os
import shutil
import time
from pathlib import Path
from typing import Iterator

logger = logging.getLogger(__name__)


def scantree(path: Path) -> Iterator[Path]:
    """Recursively yield Path objects for given directory.

    :param path: Directory or file path.
    :return: Yields all file system entries (dirs and files) as Path object.
    """
    if path.is_file():
        yield path
    else:
        yield from _scantree(path)


def _scantree(path: Path):
    """Private function for yielding Path objects for given directory recursively.

    :param path: Directory path.
    :return: Yields all file system entries (dirs and files) as Path object.
    """
    for entry in os.scandir(path):
        if entry.is_dir(follow_symlinks=False):
            yield from _scantree(Path(entry.path))
        else:
            yield Path(entry.path)


@contextlib.contextmanager
def timer(message: str):
    """Simple context manager based timer for logging timings.

    from https://github.com/deepmind/alphafold/blob/b8accc68e93500d3d3e811e8b84675d33ad6c25c/alphafold/data/tools/utils.py
    :param message: Message to log with the measure time.
    :return:
    """
    logger.info(f"Starting: {message}")
    tic = time.time()
    yield
    toc = time.time()
    logger.info(f"Finished: {message} Took: {(toc - tic):.3f} seconds")


def unpack_gz(gz_file, out_file) -> None:
    """Unpacks a gunzipped file to a file.

    from https://stackoverflow.com/a/44712152
    :param gz_file: The gunzipped file.
    :param out_file: The new unpacked file.
    :return: None
    """
    with gzip.open(gz_file, "rb") as f_in:
        with open(out_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def _count_lines_make_gen(reader):
    """Helper generator for counting lines in file.

    from https://stackoverflow.com/a/27518377
    :param reader: The file reader.
    :return: yields the next bytes of reader.
    """
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024 * 1024)


def count_lines(path: Path):
    """Counts lines in file

    from https://stackoverflow.com/a/27518377
    :param filename: The file path to count lines for.
    :return: Number of lines in filename
    """
    with open(path, "rb") as f:
        f_gen = _count_lines_make_gen(f.raw.read)
        return sum(buf.count(b"\n") for buf in f_gen)


def parse_microminer_search_stdout(stdout: str) -> dict:
    info_dict = {
        "input_file": None,
        "index_read_time": None,
        "query_site_build_time": None,
        "kmer_search_int_time": None,
        "nof_candidates": None,
        "complex_read_time": None,
        "align_time": None,
        "full_align_time": None,  # activeSite-alignment + complex reading + overhead
        "aligned_candidates": None,
        "query_sites_with_hits": None,
        "unique_complex_hits": None,
        "all_site_alignment_hits": None,
        "search_time": None,
    }
    for line in stdout.split("\n"):
        line = line.strip()
        if line.startswith("Working on "):
            info_dict["input_file"] = line[11:]
        elif line.startswith("Reading Kmer Index Files from disc | Took:"):
            info_dict["index_read_time"] = float(line.split(" ")[-2])
        elif line.startswith("Calculate query sites | Took:"):
            info_dict["query_site_build_time"] = float(line.split(" ")[-2])
        elif line.startswith("Kmer search | Took:"):
            info_dict["kmer_search_int_time"] = float(line.split(" ")[-2])
        elif line.startswith("Kmer search generated") and line.endswith("candidates"):
            info_dict["nof_candidates"] = int(line.split(" ")[-2])
        elif line.startswith("complexReadTimeSum:"):
            info_dict["complex_read_time"] = float(line.split(" ")[-1])
        elif line.startswith("alignmentTimeSum:"):
            info_dict["align_time"] = float(line.split(" ")[-1])
        elif line.startswith("ActiveSite alignments (INT) | Took:"):
            info_dict["full_align_time"] = float(line.split(" ")[-2])
        elif line.startswith("Successfully aligned candidates:"):
            info_dict["aligned_candidates"] = int(line.split(" ")[3])
        elif line.startswith("Query sites with hits:"):
            info_dict["query_sites_with_hits"] = int(line.split(" ")[-1])
        elif line.startswith("Unique complex hits:"):
            info_dict["unique_complex_hits"] = int(line.split(" ")[-1])
        elif line.startswith("All site alignment hits:"):
            info_dict["all_site_alignment_hits"] = int(line.split(" ")[-1])
        elif line.startswith("Run search | Took:"):
            info_dict["search_time"] = float(line.split(" ")[-2])
    return info_dict
