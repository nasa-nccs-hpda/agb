#!/usr/bin/python

import argparse
import logging
import pathlib
import sys

from agb.model.ViAirHyper import ViAirHyper
from core.model.GeospatialImageFile import GeospatialImageFile


# -----------------------------------------------------------------------------
# main
#
# This prepares airborne hyperspectral data.
#
# agb/view/viAirHyperCLV.py -i /explore/nobackup/projects/ilab/projects/AGB/testMLBS/airHyper/clipped/2017_NEON_D07_MLBS_DP3_537000_4132000_reflectance.tif -o /explore/nobackup/projects/ilab/projects/AGB/testMLBS/airHyper/vi
# -----------------------------------------------------------------------------
def main():

    # Process command-line args.
    desc = 'Use this application to compute vegetation indicies for ' + \
           'airborne hyperspectral data.'

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('-i',
                        type=pathlib.Path,
                        required=True,
                        help='Path to the image file')

    parser.add_argument('-o',
                        type=pathlib.Path,
                        default='.',
                        help='Path to output directory')

    args = parser.parse_args()

    # ---
    # Validate input
    # ---
    if not args.i.is_file() or not args.i.exists():
        raise FileNotFoundError(args.i)

    if not args.o.is_dir():
        raise NotADirectoryError(args.o)

    # ---
    # Instantiate a logger.
    # ---
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)

    # ---
    # Instantiate the GeospatialImageFile
    # ---
    gif = GeospatialImageFile(str(args.i), logger=logger)

    # ---
    # Instantiate the VI class and run.
    # ---
    vi = ViAirHyper(gif, args.o, logger)
    vi.computeAllAndWrite()


# -----------------------------------------------------------------------------
# Invoke the main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    sys.exit(main())
