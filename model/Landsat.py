
from pathlib import Path
import shutil

from agb.model.BaseProcess import BaseProcess
from core.model.Envelope import Envelope
from core.model.GeospatialImageFile import GeospatialImageFile


# -----------------------------------------------------------------------------
# class Landsat
# -----------------------------------------------------------------------------
class Landsat(BaseProcess):

    # -------------------------------------------------------------------------
    # __init__
    # -------------------------------------------------------------------------
    # def __init__(self,
    #              site: NeonSite,
    #              outDir: Path,
    #              logger: logging.RootLogger):
    #
    #     super(Landsat, self).__init__(site, outDir, logger)

    # -------------------------------------------------------------------------
    # clip
    # -------------------------------------------------------------------------
    def _myClip(self,
                inFiles: list,
                envelope: Envelope,
                outFile: Path) -> None:

        # There should be only one file for Landsat.
        if len(inFiles) > 1:
            raise RuntimeError('Multiple Landsat files found.')

        # This is repeated in Lidar.py
        clipFile = inFiles[0]
        tempGif = GeospatialImageFile(str(clipFile))

        clipEnv = self._reprojEnv(tempGif.getDataset().GetSpatialRef(),
                                  envelope)

        shutil.copyfile(clipFile, outFile)
        gif = GeospatialImageFile(str(outFile), logger=self._logger)
        gif.clipReproject(clipEnv)

    # -------------------------------------------------------------------------
    # defaultInputDir
    # -------------------------------------------------------------------------
    @staticmethod
    def defaultInputDir():
        return Path('/explore/nobackup/projects/ilab/data/AGB/LCMS')

    # -------------------------------------------------------------------------
    # findFiles
    # -------------------------------------------------------------------------
    def findFiles(self, searchDir: Path = None) -> list:

        searchDir = searchDir or self.defaultInputDir()
        files = list(searchDir.rglob('*.tif'))

        filtered = [f for f in files if
                f.stem.split('-')[1] == '7_Most_Recent_Year_of_Fast_Loss']

        return filtered

    # -------------------------------------------------------------------------
    # processName
    # -------------------------------------------------------------------------
    def _processName(self) -> str:
        return 'landsat'
