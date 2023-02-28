#!/usr/bin/python

import argparse
import logging
import pathlib
import sys

from agb.model.DataTypeAirHyper import DataTypeAirHyper
from agb.model.NeonSite import NeonSite


# -----------------------------------------------------------------------------
# main
#
# This prepares airborne hyperspectral data.
#
# agb/view/airHyperMatchLidarCLV.py -i /explore/nobackup/projects/ilab/data/AGB/Airborne_Hyperspectral -o /explore/nobackup/projects/ilab/projects/AGB/testMLBS/airHyper -c /explore/nobackup/people/rlgill/innovation-lab-repositories/agb/model/tests/NeonSites.csv -l /explore/nobackup/projects/ilab/data/AGB/Airborne_Lidar/DP3.30024.001/neon-aop-products/2017/FullSite/D07/2017_MLBS_2/L3/DiscreteLidar/CHMGtif
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

    parser.add_argument('-l',
                        type=pathlib.Path,
                        default='.',
                        help='Path to the lidar file to match')

    parser.add_argument('-o',
                        type=pathlib.Path,
                        default='.',
                        help='Path to output directory')

    args = parser.parse_args()

    # ---
    # Validate input
    # ---
    if args.c and not args.c.is_file():
        raise FileNotFoundError(args.c)

    if not args.l.is_dir():
        raise NotADirectoryError(args.l)

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
    # Instantiate the data type.
    # ---
    dataType = DataTypeAirHyper()

    # ---
    # Instantiate the site and its clipper.  We will not know which site until
    # the first lidar image is read.  This expects all lidar files to be for
    # the same site.
    # ---
    clipper = None
    site = None

    # ---
    # Once the site is detected, find all the hyper files for it.  Later,
    # filter this list based on the lidar files passed.
    # ---
    allHyperFiles = None

    # ---
    # Find the lidar files, then find matching ones in the input directory.
    # ---
    lidarFiles = args.l.glob('*.tif')
    matchedStrings = []

    for fl in lidarFiles:

        # ---
        # lidar: 2017_NEON_D07_MLBS_DP3_537000_4132000_RUG_01101623
        # hyper: 2015_NEON_D07_MLBS_DP3_541000_4136000_reflectance.h5
        # ---
        lidarMatchString = fl.stem.split('_')[:7]

        if lidarMatchString in matchedStrings:
            continue

        logger.info('Finding mate for ' + '_'.join(lidarMatchString))

        # Instantiate the site.
        if not site:

            siteAbbreviation = lidarMatchString[3]
            site = NeonSite()
            site.fromCSV(args.c, siteAbbreviation)
            clipper = dataType.clipper(site, args.o, logger)
            allHyperFiles = dataType.findFiles(args.i, site)
            logger.info('Found ' + str(len(allHyperFiles)) + ' files.')

        # Filter that list for year, east and north.
        for fa in allHyperFiles:

            if fa.stem.split('_')[:7] == lidarMatchString:

                clipper.clip([fa])
                matchedStrings.append(lidarMatchString)
                break


# -----------------------------------------------------------------------------
# Invoke the main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    sys.exit(main())
