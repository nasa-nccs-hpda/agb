#!/usr/bin/python

import argparse
from pathlib import Path
import sys


# -----------------------------------------------------------------------------
# main
#
# agb/view/completionStatus.py -o /explore/nobackup/projects/ilab/projects/AGB --tileFile /explore/nobackup/people/rlgill/SystemTesting/testAGB3/MLBS_subtile_list.txt -s MLBS -y 2017 
# -----------------------------------------------------------------------------
def main():

    # Process command-line args.
    desc = 'Use this application to report the completion status of a ' + \
           'set of years and tiles for airborne hyperspectral processing.'
           
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('-o',
                        type=Path,
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

    # Validate output directory.
    if args.o and not args.o.is_dir():
        raise NotADirectoryError(args.o)

    # Collect the tile IDs to check.
    subtiles = args.subtiles
    
    if args.tileFile:
        
        with open(args.tileFile) as f:
            lines = f.readlines()
        
        subtiles = [l.strip() for l in lines]
        
    # Prepare the glob search.
    site = args.s
    year = str(args.y)
    yearSiteDir = args.o / (site + '-' + year)
    
    if not yearSiteDir.exists():
        
        raise RuntimeError('Output directory, ' + 
                           str(yearSiteDir) + 
                           ' does not exist.')

    # Check by tile ID.
    completeTids = []
    incompleteTids = []
    
    for tid in subtiles:
        
        tidGlob = year + '_NEON_D07_' + site + '_DP3_' + tid + '_reflectance'
        ref = yearSiteDir / (tidGlob + '.tif')
        anc = yearSiteDir / (tidGlob + '_ancillary.tif')
        ind = yearSiteDir / (tidGlob + '_indices.tif')
        misspelled = yearSiteDir / (tidGlob + '_indicies.tif')
        sav = yearSiteDir / (tidGlob + '_savgol.tif')
        
        if ref.exists() and \
           anc.exists() and \
           (ind.exists() or misspelled.exists()) and \
           sav.exists():
           
            completeTids.append(tid)
           
        else:
            incompleteTids.append(tid)
            
    print('Total:', str(len(subtiles)))
    print('Complete:', len(completeTids))
    print('Incomplete:', len(incompleteTids))
    print('Incomplete:', incompleteTids)
            

# -----------------------------------------------------------------------------
# Invoke the main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    sys.exit(main())        
