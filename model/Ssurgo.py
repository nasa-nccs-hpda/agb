
import logging
from pathlib import Path
import tempfile

from osgeo import gdal

from agb.model.BaseProcess import BaseProcess
from agb.model.NeonSite import NeonSite
from core.model.Envelope import Envelope
from core.model.GeospatialImageFile import GeospatialImageFile


# -----------------------------------------------------------------------------
# class Ssurgo
# -----------------------------------------------------------------------------
class Ssurgo(BaseProcess):

    # -------------------------------------------------------------------------
    # __init__
    # -------------------------------------------------------------------------
    # def __init__(self,
    #              site: NeonSite,
    #              outDir: Path,
    #              logger: logging.RootLogger):
    #
    #     super(Ssurgo, self).__init__(site, outDir, logger)

    # -------------------------------------------------------------------------
    # defaultInputDir
    # -------------------------------------------------------------------------
    @staticmethod
    def defaultInputDir():
        return Path('/explore/nobackup/projects/ilab/data/AGB/SSURGO')

    # -------------------------------------------------------------------------
    # findFiles
    # -------------------------------------------------------------------------
    def findFiles(self, searchDir: Path = None) -> list:

        if not searchDir:

            searchDir = self.defaultInputDir() / self._site.abbreviation

        files = list(searchDir.rglob('*.tif'))
        return files

    # -------------------------------------------------------------------------
    # getOutFileName
    # -------------------------------------------------------------------------
    def _getOutFileName(self) -> str:
        return 'ssurgo'

    # -------------------------------------------------------------------------
    # myClip
    # -------------------------------------------------------------------------
    # def _myClip(self,
    #             inFiles: list,
    #             envelope: Envelope,
    #             outFile: Path) -> None:
    #
    #     # The target SRS is that of the site.
    #     tempGif = GeospatialImageFile(str(inFiles[0]))
    #     cEnv = self._reprojEnv(tempGif.getDataset().GetSpatialRef(), envelope)
    #
    #
    #
    #
    #
    #
    #     cmd = 'ogr2ogr -f "ESRI Shapefile" ' + \
    #           outFile + ' ' + \
    #           inFile
    #
    #     SystemCommand(cmd, logger, raiseException=True)
