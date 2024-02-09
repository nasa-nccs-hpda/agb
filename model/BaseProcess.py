import logging
from pathlib import Path

from osgeo.osr import SpatialReference

from agb.model.NeonSite import NeonSite
from core.model.Envelope import Envelope


# -----------------------------------------------------------------------------
# class BaseProcess
#
# Originally, this set of AGB code was supposed to clip and reproject various
# types of input data.  Only hyperspectral ended up being used.  Hyperspectral
# has its own stack of four types of images.  In light of all this, the
# generalization of this base class is unnecessary.
# -----------------------------------------------------------------------------
class BaseProcess(object):

    # -------------------------------------------------------------------------
    # __init__
    # -------------------------------------------------------------------------
    def __init__(self,
                 site: NeonSite,
                 outDir: Path,
                 logger: logging.RootLogger):

        if not outDir or not outDir.exists() or not outDir.is_dir():
            raise RuntimeError('Invalid output directory: ' + str(outDir))

        self._logger = logger
        self._outDir = outDir
        self._site = site

    # -------------------------------------------------------------------------
    # clip
    # -------------------------------------------------------------------------
    def clip(self, inFiles: list, envelope: Envelope) -> None:

        outFile = self._outDir / (self._getOutFileName() + '.tif')
        self._logger.info('Clipping to ' + str(outFile))

        if outFile.exists():

            self._logger.info('Clipped file already exists.  Not clipping.')
            return

        env = envelope or self._site.envelope
        self._myClip(inFiles, env, outFile)

    # -------------------------------------------------------------------------
    # myClip
    # -------------------------------------------------------------------------
    def _myClip(self,
                inFiles: list,
                envelope: Envelope,
                outFile: Path) -> None:

        raise NotImplementedError()

    # -------------------------------------------------------------------------
    # defaultInputDir
    # -------------------------------------------------------------------------
    @staticmethod
    def defaultInputDir():
        raise NotImplementedError()

    # -------------------------------------------------------------------------
    # findFiles
    # -------------------------------------------------------------------------
    def findFiles(self, searchDir: Path = None) -> list:
        raise NotImplementedError()

    # -------------------------------------------------------------------------
    # getOutFileName
    # -------------------------------------------------------------------------
    def _getOutFileName(self) -> str:
        raise NotImplementedError()

    # -------------------------------------------------------------------------
    # reprojEnv
    # -------------------------------------------------------------------------
    def _reprojEnv(self,
                   srs: SpatialReference,
                   envelope: Envelope) -> Envelope:

        env = envelope.Clone()

        if not env.GetSpatialReference().IsSame(srs):

            self._logger.info('Clipping env. to ' + self._site.name)
            env.TransformTo(srs)

        return env

    # -------------------------------------------------------------------------
    # run
    # -------------------------------------------------------------------------
    def run(self, envelope: Envelope = None) -> None:

        inFiles = self.findFiles()
        numFiles = len(inFiles)
        self._logger.info('Found ' + str(numFiles) + ' files.')

        if numFiles == 0:
            self._logger.info('Quitting')

        else:
            self.clip(inFiles, envelope)
