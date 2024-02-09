#!/usr/bin/python

import argparse
import logging
import pathlib
import sys

from agb.model.Hyper import Hyper
from agb.model.NeonSite import NeonSite


# -----------------------------------------------------------------------------
# main

# This prepares airborne hyperspectral data.
#
# agb/view/runHyperCLV.py -o /explore/nobackup/projects/ilab/projects/AGB -c /explore/nobackup/people/rlgill/innovation-lab-repositories/agb/model/tests/NeonSites.csv -s MLBS -y 2017 --subtiles 537000_4131000
#
# agb/view/runHyperCLV.py -o /explore/nobackup/projects/ilab/projects/AGB -c /explore/nobackup/people/rlgill/innovation-lab-repositories/agb/model/tests/NeonSites.csv -s MLBS --tileFile /explore/nobackup/people/rlgill/SystemTesting/testAGB3/MLBS_subtile_list.txt -y 2017 
#
# agb/view/runHyperCLV.py -o /explore/nobackup/people/rlgill/SystemTesting/testAGB3 -c /explore/nobackup/people/rlgill/innovation-lab-repositories/agb/model/tests/NeonSites.csv -s MLBS -y 2015 --subtiles 547000_4134000
# -----------------------------------------------------------------------------
def main():

    # Process command-line args.
    desc = 'Use this application to prepare hyperspectral data for the ' + \
           'AGB project.'
           
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

    parser.add_argument('-y',
                        type=int,
                        required=True,
                        help='Year to process')

    group = parser.add_mutually_exclusive_group()

    group.add_argument('--subtiles',
                       nargs='+',
                       help='537000_4131000 537000_4132000 ...')

    group.add_argument('--tileFile',
                       help='Path to file of tile IDs')

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

    subtiles = args.subtiles
    
    if args.tileFile:
        
        with open(args.tileFile) as f:
            lines = f.readlines()
        
        subtiles = [l.strip() for l in lines]
        
    # ---
    # Run it.
    # ---
    Hyper(site, args.o, logger, args.y).run(subtiles)


# -----------------------------------------------------------------------------
# Invoke the main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    sys.exit(main())
