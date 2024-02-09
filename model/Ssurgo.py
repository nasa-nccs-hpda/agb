
from pathlib import Path

from agb.model.BaseProcess import BaseProcess


# -----------------------------------------------------------------------------
# class Ssurgo
# -----------------------------------------------------------------------------
class Ssurgo(BaseProcess):

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
