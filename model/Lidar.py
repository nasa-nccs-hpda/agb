
import logging
from pathlib import Path
import tempfile

from osgeo import gdal

from agb.model.BaseProcess import BaseProcess
from agb.model.NeonSite import NeonSite
from core.model.Envelope import Envelope
from core.model.GeospatialImageFile import GeospatialImageFile


# -----------------------------------------------------------------------------
# class Lidar
# -----------------------------------------------------------------------------
class Lidar(BaseProcess):

    # -------------------------------------------------------------------------
    # __init__
    # -------------------------------------------------------------------------
    def __init__(self,
                 site: NeonSite,
                 outDir: Path,
                 logger: logging.RootLogger,
                 year: int,
                 subsite: str):

        super().__init__(site, outDir, logger)

        self._subsite = subsite
        self._year = year

    # -------------------------------------------------------------------------
    # myclip
    # -------------------------------------------------------------------------
    def _myClip(self,
                inFiles: list,
                envelope: Envelope,
                outFile: Path) -> None:

        # The target SRS is that of the site.
        tempGif = GeospatialImageFile(str(inFiles[0]))
        cEnv = self._reprojEnv(tempGif.getDataset().GetSpatialRef(), envelope)

        chmTmp = self._createBandVRT(inFiles, 'CHM')
        rugTmp = self._createBandVRT(inFiles, 'RUG')
        options = gdal.BuildVRTOptions(separate=True, bandList=[1, 2])
        gdal.BuildVRT(str(outFile), [chmTmp, rugTmp], options=options)

        # ---
        # We should be able to clip by giving the vrt the output bounds.
        # However, the resultant file is still a vrt, and not a "realized"
        # raster, despite setting the vrt to None.
        # ---
        gis = GeospatialImageFile(str(outFile), outputFormat='GTiff')
        gis.clipReproject(envelope=cEnv)
        gis.getDataset().GetRasterBand(1).SetMetadata({'Name': 'CHM'})
        gis.getDataset().GetRasterBand(2).SetMetadata({'Name': 'RUG'})
        gis = None

    # -------------------------------------------------------------------------
    # createBandVRT
    # -------------------------------------------------------------------------
    def _createBandVRT(self, inFiles: list, fileType: str) -> str:

        files = [str(f) for f in inFiles if f.stem.split('_')[-1] == fileType]
        tmp = tempfile.mkstemp(suffix='.tif')[1]
        gdal.BuildVRT(tmp, files)
        return tmp

    # -------------------------------------------------------------------------
    # defaultInputDir
    # -------------------------------------------------------------------------
    @staticmethod
    def defaultInputDir() -> Path:

        return Path('/explore/nobackup/projects/ilab/data/AGB/' +
                    'Airborne_Lidar/FSD')

    # -------------------------------------------------------------------------
    # findFiles
    # -------------------------------------------------------------------------
    def findFiles(self, searchDir: Path = None) -> list:

        searchDir = searchDir or Lidar.defaultInputDir()
        files = list(searchDir.rglob('*' + self._site.abbreviation + '*.tif'))
        filtered = [f for f in files if self._subsite in str(f)]
        filtered = [f for f in filtered if str(self._year) in str(f)]

        filtered = [f for f in filtered if
            str(f).rsplit('_', maxsplit=1)[-1] not in ['DSM.tif', 'DTM.tif']]

        return filtered

    # -------------------------------------------------------------------------
    # getOutFileName
    # -------------------------------------------------------------------------
    def _getOutFileName(self) -> str:
        return 'lidar'
