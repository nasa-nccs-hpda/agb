# agb

The AGB application's original purpose was to collect and prepare imagery from a number of sources.  While there exists code to process Landsat, Lidar and Ssurgo, the only data processed in production is Hyperspectral.  The hyperspectral output includes several vegetation indices, reflectances, reflectances filtered by Savgol, and certain metadata.

## Running

- The primary command-line application is runHyperCLV.py.  Use its "-h" option to for help information.  

- The input files are in h5 format.  Given the site, year and tile IDs, the application will find the correct ones and process them.  Note, this software was not built to run outside its development network.  It refers to absolute paths.

- To run, it requires a site definition file in CSV format.  There is one included with the software.  The site abbreviation, as seen in the site definition file, is used to specify which site to run.  

- You must also specify the year associated with the input files to run. 

- To specify which tiles to run, it accepts a series of subtile coordinates or the same IDs in a file with one ID per line.  Refer to the command's help information.  

## Tile File Example

537000_4131000
537000_4132000
537000_4133000
537000_4134000
537000_4135000

## Sample Command

agb/view/runHyperCLV.py -o /path/to/AGB/output -c /path/to/NeonSites.csv -s MLBS -y 2017 --subtiles 537000_4131000 537000_4132000 537000_4133000 537000_4134000 537000_4135000