from pathlib import Path
import argparse
import sys
import logging

from helper import CONFIG
from helper.mol_utils import read_plddt_values
from helper.utils import scantree

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="""
        Extract pLDDT values from a list of Alphafold PDB files.
        """
    )

    parser.add_argument(
        "--dir",
        "-d",
        required=False,
        type=str,
        default=CONFIG["DATA"]["AFDB_DIR"],
        help="Directory to read PDB files from.",
    )
    parser.add_argument(
        "--outfile", "-o", type=str, required=True, help="Path to output file"
    )

    args = parser.parse_args()

    dir = Path(args.dir)
    outfile = Path(args.outfile)

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s"
    )
    logger.info(f'Starting scripts: {" ".join(sys.argv)}')

    if not dir.is_dir():
        print("Error: Directory to read structure files from does not exist.")
        sys.exit(1)

    header = True
    for path in scantree(dir):
        if path.name.endswith(".pdb.gz"):
            df_lddt = read_plddt_values(path)

            df_lddt.to_csv(
                outfile,
                index=False,
                header=header,
                sep="\t",
                mode="w" if header else "a",
            )
            header = False


if __name__ == "__main__":
    main()
