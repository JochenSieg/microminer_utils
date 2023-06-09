import logging
from pathlib import Path

import pandas as pd

import helper

logger = logging.getLogger(__name__)


def sanity_check_microminer_result_dir(
    result_dir: Path, input_df: pd.DataFrame
) -> bool:
    """Performs a sanity check on the results from the cluster run.

    :param result_dir: Directory with results.
    :param input_df: Dataframe with query data.
    :return: True if result seem sane, false otherwise.
    """
    files = [
        p for p in helper.utils.scantree(result_dir) if p.name == "resultStatistic.csv"
    ]
    is_sane = True
    if len(files) != input_df.shape[0]:
        logger.warning(
            f"Number of resultStatistic.csv differs from input rows: file={len(files)}"
            f" rows={input_df.shape[0]} in dir={result_dir.resolve()}"
        )
        is_sane = False
    else:
        logger.info(
            f"Check successful: Number of resultStatistic.csv equals input rows:"
            f" file={len(files)}"
            f" rows={input_df.shape[0]} in dir={result_dir.resolve()}"
        )
    return is_sane
