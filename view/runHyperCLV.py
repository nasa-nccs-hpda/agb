#!/usr/bin/python

import argparse
import logging
import pathlib
import sys

from osgeo import osr
from osgeo.osr import SpatialReference

from agb.model.Hyper import Hyper
from agb.model.NeonSite import NeonSite
from core.model.Envelope import Envelope


# -----------------------------------------------------------------------------
# main
#
# This prepares airborne hyperspectral data.
#
# agb/view/runHyperCLV.py -o /explore/nobackup/people/rlgill/SystemTesting/testAGB2 -c /explore/nobackup/people/rlgill/innovation-lab-repositories/agb/model/tests/NeonSites.csv -s MLBS -y 2017
#
# ------------------ 13km Box ------------------
# MLBS center:  (542067.6, 4136943) (easting, northing)
# ul = (542067.6 - 6500, 4136943 + 6500) = (535567.6, 4143443)
# lr = (542067.6 + 6500, 4136943 - 6500) = (548567.6, 4130443)
#
# agb/view/runHyperCLV.py -o /explore/nobackup/people/rlgill/SystemTesting/testAGB2 -c /explore/nobackup/people/rlgill/innovation-lab-repositories/agb/model/tests/NeonSites.csv -s MLBS -y 2017 -e 535567.6 4143443 548567.6 4130443 32617
# -----------------------------------------------------------------------------
def main():

    # Process command-line args.
    desc = 'Use this application to prepare data for the AGB project.'
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('-c',
                        type=pathlib.Path,
                        required=True,
                        help='Path to CSV site file')

    parser.add_argument('-e',
                        nargs=5,
                        required=False,
                        help='ulx uly lrx lry epsg-code')

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

    args = parser.parse_args()

    # ---
    # Validate input
    # ---
    if args.c and not args.c.is_file():
        raise FileNotFoundError(args.c)

    if args.o and not args.o.is_dir():
        raise NotADirectoryError(args.o)

    # ---
    # Instantiate an envelope.
    # ---
    env = None

    if args.e:

        epsgCode = int(args.e[4])
        srs = SpatialReference()
        srs.ImportFromEPSG(epsgCode)

        srs4326 = SpatialReference()
        srs4326.ImportFromEPSG(4326)

        if srs.IsSame(srs4326):
            srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

        env = Envelope()
        env.addPoint(float(args.e[0]), float(args.e[1]), 0, srs)
        env.addPoint(float(args.e[2]), float(args.e[3]), 0, srs)

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
    Hyper(site, args.o, logger, args.y).run(env)


# -----------------------------------------------------------------------------
# Invoke the main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    sys.exit(main())
