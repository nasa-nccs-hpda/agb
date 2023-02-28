#!/usr/bin/python

import argparse
import logging
import pathlib
import sys

from agb.model.DataTypeLidar import DataTypeLidar
from agb.model.NeonSite import NeonSite


# -----------------------------------------------------------------------------
# main
#
# This prepares LIDAR data.
#
# agb/view/prepLidarCLV.py -i /explore/nobackup/projects/ilab/data/AGB/Airborne_Lidar/DP3.30024.001/neon-aop-products/2017/FullSite/D07/2017_MLBS_2/L3/DiscreteLidar/CHMGtif -o /explore/nobackup/projects/ilab/projects/AGB/testMLBS/lidar -c /explore/nobackup/people/rlgill/innovation-lab-repositories/agb/model/tests/NeonSites.csv -s MLBS --siteLength 10000
# -----------------------------------------------------------------------------
def main():

    # Process command-line args.
    desc = 'Use this application to prepare data for the AGB project.'
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('-c',
                        type=pathlib.Path,
                        required=True,
                        help='Path to CSV site file')

    parser.add_argument('-i',
                        type=pathlib.Path,
                        default='.',
                        help='Path to the input data directory')

    parser.add_argument('-o',
                        type=pathlib.Path,
                        default='.',
                        help='Path to output directory')

    parser.add_argument('-s',
                        type=str,
                        help='Site abbreviation')

    parser.add_argument('--siteLength',
                        type=int,
                        required=False,
                        help='Length in meters of box centered on site '
                             'center point')

    args = parser.parse_args()

    # ---
    # Validate input
    # ---
    if args.c and not args.c.is_file():
        raise FileNotFoundError(args.c)

    if not args.i.is_dir():
        raise NotADirectoryError(args.i)

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

    if args.siteLength:
        site.resize(args.siteLength)

    # ---
    # Instantiate the data type and clipper.
    # ---
    dataType = DataTypeLidar()
    clipper = dataType.clipper(site, args.o, logger)

    # ---
    # Find the files in the input directory.
    # ---
    logger.info('Searching for ' + dataType.name() + ' files.')
    files = dataType.findFiles(args.i, site)

    logger.info('Found ' + str(len(files)) + ' ' + dataType.name() + \
        ' files.')

    # ---
    # Prepare the files.
    # ---
    clipper.clip(files)
    # logger.warn('Only processing the first 3 files while testing.')
    # clipper.clip(hyperFiles[:3])


# -----------------------------------------------------------------------------
# Invoke the main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    sys.exit(main())
