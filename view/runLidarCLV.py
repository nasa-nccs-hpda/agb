#!/usr/bin/python

import argparse
import logging
import pathlib
import sys

from agb.model.Lidar import Lidar
from agb.model.NeonSite import NeonSite


# -----------------------------------------------------------------------------
# main
#
# This prepares airborne hyperspectral data.
#
# agb/view/runLidar.py -o /explore/nobackup/people/rlgill/SystemTesting/testAGB -c /explore/nobackup/people/rlgill/innovation-lab-repositories/agb/model/tests/NeonSites.csv -s MLBS -y 2017 --subsite D07
# -----------------------------------------------------------------------------
def main():

    # Process command-line args.
    desc = 'Use this application to prepare data for the AGB project.'
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('-c',
                        type=pathlib.Path,
                        required=True,
                        help='Path to CSV site file')

    parser.add_argument('-o',
                        type=pathlib.Path,
                        default='.',
                        help='Path to output directory')

    parser.add_argument('-s',
                        type=str,
                        required=True,
                        help='Site abbreviation')

    parser.add_argument('--subsite',
                        type=str,
                        required=True,
                        help='Subsite code')

    parser.add_argument('-y',
                        type=int,
                        required=True,
                        help='Year to process')

    args = parser.parse_args()

    # ---
    # Validate input
    # ---
    if args.c and not args.c.is_file():
        raise FileNotFoundError(args.c)

    if args.o and not args.o.is_dir():
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
    # Instantiate the site.
    # ---
    site = NeonSite()
    site.fromCSV(args.c, args.s)

    # ---
    # Run it.
    # ---
    Lidar(site, args.o, logger, args.y, args.subsite).run()


# -----------------------------------------------------------------------------
# Invoke the main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    sys.exit(main())
