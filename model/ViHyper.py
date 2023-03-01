
import logging
from pathlib import Path

from osgeo import gdal

import numpy

from core.model.GeospatialImageFile import GeospatialImageFile


# -----------------------------------------------------------------------------
# class ViHyper
# -----------------------------------------------------------------------------
class ViHyper(object):

    # -------------------------------------------------------------------------
    # __init__
    # -------------------------------------------------------------------------
    def __init__(self,
                 airHyperImage: GeospatialImageFile,
                 outFileName: Path,
                 logger: logging) -> None:

        self._outFileName = outFileName
        self._logger = logger
        self._image = airHyperImage

        # ---
        # Index the bands to their wavelengths.  This is how to find bands
        # later.
        # ---
        numBands = self._image.getDataset().RasterCount
        self._bandIndicies = {}

        for i in range(1, numBands):

            wl = int(self._image.getDataset().GetRasterBand(i). \
                     GetMetadata()['Matched wavelength'])

            self._bandIndicies[wl] = i

        # ---
        # As we read bands, save them here.  Some of them are used with
        # multiple VIs.
        # ---
        self._bands = {}

    # -------------------------------------------------------------------------
    # computeAllAndWrite
    # -------------------------------------------------------------------------
    def computeAllAndWrite(self) -> None:

        # ---
        # Each VI is a band in the output image.  We write the VIs as we go, so
        # we must create the GDAL data set and specify the number of bands.
        # ---
        curBand = 0
        numVIs = 21  # Be sure to update this as VIs are added below.
        ds = self._createDs(numVIs)

        self._logger.info('Computing ARI')
        vi = self.computeARI()
        curBand = self._writeVi(vi, 'ARI', ds, curBand)

        self._logger.info('Computing CAI')
        vi = self.computeCAI()
        curBand = self._writeVi(vi, 'CAI', ds, curBand)

        self._logger.info('Computing CRI550')
        vi = self.computeCRI550()
        curBand = self._writeVi(vi, 'CRI550', ds, curBand)

        self._logger.info('Computing CRI700')
        vi = self.computeCRI700()
        curBand = self._writeVi(vi, 'CRI700', ds, curBand)

        self._logger.info('Computing EVI')
        vi = self.computeEVI()
        curBand = self._writeVi(vi, 'EVI', ds, curBand)

        self._logger.info('Computing EVI2')
        vi = self.computeEVI2()
        curBand = self._writeVi(vi, 'EVI2', ds, curBand)

        self._logger.info('Computing fPAR')
        vi = self.computeFPAR()
        curBand = self._writeVi(vi, 'fPAR', ds, curBand)

        self._logger.info('Computing LAI')
        vi = self.computeLAI()
        curBand = self._writeVi(vi, 'LAI', ds, curBand)

        self._logger.info('Computing MCTI')
        vi = self.computeMCTI()
        curBand = self._writeVi(vi, 'MCTI', ds, curBand)

        self._logger.info('Computing MSI')
        vi = self.computeMSI()
        curBand = self._writeVi(vi, 'MSI', ds, curBand)

        self._logger.info('Computing NDII')
        vi = self.computeNDII()
        curBand = self._writeVi(vi, 'NDII', ds, curBand)

        self._logger.info('Computing NDLI')
        vi = self.computeNDLI()
        curBand = self._writeVi(vi, 'NDLI', ds, curBand)

        self._logger.info('Computing NDNI')
        vi = self.computeNDNI()
        curBand = self._writeVi(vi, 'NDNI', ds, curBand)

        self._logger.info('Computing NDVI')
        vi = self.computeNDVI()
        curBand = self._writeVi(vi, 'NDVI', ds, curBand)

        self._logger.info('Computing NDWI')
        vi = self.computeNDWI()
        curBand = self._writeVi(vi, 'NDWI', ds, curBand)

        self._logger.info('Computing NIRv')
        vi = self.computeNIRv()
        curBand = self._writeVi(vi, 'NIRv', ds, curBand)

        self._logger.info('Computing PRIn')
        vi = self.computePRIn()
        curBand = self._writeVi(vi, 'PRIn', ds, curBand)

        self._logger.info('Computing PRIw')
        vi = self.computePRIw()
        curBand = self._writeVi(vi, 'PRIw', ds, curBand)

        self._logger.info('Computing SAVI')
        vi = self.computeSAVI()
        curBand = self._writeVi(vi, 'SAVI', ds, curBand)

        self._logger.info('Computing WBI')
        vi = self.computeWBI()
        curBand = self._writeVi(vi, 'WBI', ds, curBand)

        self._logger.info('Computing Albedo')
        vi = self.computeAlbedo()
        curBand = self._writeVi(vi, 'Albedo', ds, curBand)

        del ds

    # -------------------------------------------------------------------------
    # computeAlbedo
    # -------------------------------------------------------------------------
    def computeAlbedo(self):

        band = self._image.getDataset().ReadAsArray(band_list=[27]). \
               astype('float')

        band[band==-9999] = numpy.nan
        return band

    # -------------------------------------------------------------------------
    # computeARI
    # -------------------------------------------------------------------------
    def computeARI(self) -> numpy.ndarray:

        b550 = self._getBand(550)
        b700 = self._getBand(700)
        vi = b550 - b700
        return vi

    # -------------------------------------------------------------------------
    # computeCAI
    # -------------------------------------------------------------------------
    def computeCAI(self) -> numpy.ndarray:

        b2019 = self._getBand(2019)
        b2109 = self._getBand(2109)
        b2206 = self._getBand(2206)
        vi = 0.5 * (b2019 + b2206) - b2109
        return vi

    # -------------------------------------------------------------------------
    # computeCRI700
    # -------------------------------------------------------------------------
    def computeCRI700(self) -> numpy.ndarray:

        b510 = self._getBand(510)
        b700 = self._getBand(700)
        vi = b510 - b700
        return vi

    # -------------------------------------------------------------------------
    # computeCRI550
    # -------------------------------------------------------------------------
    def computeCRI550(self) -> numpy.ndarray:

        b510 = self._getBand(510)
        b550 = self._getBand(550)
        vi = b510 - b550
        return vi

    # -------------------------------------------------------------------------
    # computeEVI
    #
    # 2.5*(835 - 660)/(835 + 6*660 - 7.5*485 + 1)
    #
    # b660 contains 463 zeros, causing division by zero.
    # -------------------------------------------------------------------------
    def computeEVI(self) -> numpy.ndarray:

        b485 = self._getBand(485)
        b660 = self._getBand(660)
        b835 = self._getBand(835)
        denom = b835 + 6 * b660 - 7.5 * b485 + 1
        denom[denom==0] = self._getMin(denom)
        vi = 2.5 * (b835 - b660) / denom
        return vi

    # -------------------------------------------------------------------------
    # computeEVI2
    # -------------------------------------------------------------------------
    def computeEVI2(self) -> numpy.ndarray:

        b660 = self._getBand(660)
        b835 = self._getBand(835)
        denom = b835 + b660 + 1
        denom[denom==0] = self._getMin(denom)
        vi = 2.5 * (b835 - b660) / denom
        return vi

    # -------------------------------------------------------------------------
    # computeFPAR
    #
    # c*[1-a*exp(-b*LAI)] where a, b and c are set to 1, 0.4 and 1, respectively
    # -------------------------------------------------------------------------
    def computeFPAR(self) -> numpy.ndarray:

        a = 1.0
        b = 0.4
        c = 1.0
        lai = self.computeLAI()
        vi = c * (1 - a * numpy.exp(-b * lai))
        return vi

    # -------------------------------------------------------------------------
    # computeLAI
    #
    # -(1/c)ln[(a - SAVI)/b] where a, b and c are set to 0.82, 0.78 and 0.6,
    # respectively
    # -------------------------------------------------------------------------
    def computeLAI(self) -> numpy.ndarray:

        aa = 0.82
        bb = 0.78
        cc = 0.6
        savi = self.computeSAVI()
        arg = (aa - savi) / bb
        arg[arg<=0] = abs(numpy.finfo(numpy.float64).min)
        vi = -(1 / cc) * numpy.log(arg)
        return vi

    # -------------------------------------------------------------------------
    # computeMCTI
    # -------------------------------------------------------------------------
    def computeMCTI(self) -> numpy.ndarray:

        b681 = self._getBand(681)
        b709 = self._getBand(709)
        b754 = self._getBand(754)

        denom = b709 - b681
        denom[denom==0] = self._getMin(denom)
        vi = (b754 - b709) / denom
        return vi

    # -------------------------------------------------------------------------
    # computeMSI
    # -------------------------------------------------------------------------
    def computeMSI(self) -> numpy.ndarray:

        b819 = self._getBand(819)
        b1599 = self._getBand(1599)
        b819[b819==0] = self._getMin(b819)
        vi = b1599 / b819
        return vi

    # -------------------------------------------------------------------------
    # computeNDII
    # -------------------------------------------------------------------------
    def computeNDII(self) -> numpy.ndarray:

        b819 = self._getBand(819)
        b1649 = self._getBand(1649)
        denom = b819 + b1649
        denom[denom==0] = self._getMin(denom)
        vi = (b819 - b1649) / denom
        return vi

    # -------------------------------------------------------------------------
    # _computeLog
    # -------------------------------------------------------------------------
    def _computeLog(self,
                    numerator: float,
                    array: numpy.ndarray) -> numpy.ndarray:

        array[array==0] = self._getMin(array)  # 0 to very small for division
        fraction = numerator / array
        fraction[array<=0] = abs(self._getMin(fraction))  # >0 for log
        return numpy.log10(fraction)

    # -------------------------------------------------------------------------
    # computeNDLI
    #
    # [log(1/1754) - log(1/1680)]/[log(1/1754) + log(1/1680)]
    # -------------------------------------------------------------------------
    def computeNDLI(self) -> numpy.ndarray:

        b1680 = self._getBand(1680)
        b1754 = self._getBand(1754)

        log1 = self._computeLog(1, b1754)
        log2 = self._computeLog(1, b1680)

        denom = log1 + log2
        denom[denom==0] = self._getMin(denom)

        vi = (log1 - log2) / denom

        return vi

    # -------------------------------------------------------------------------
    # computeNDNI
    #
    # [log(1/1510) - log(1/1680)]/[log(1/1510) + log(1/1680)]
    # -------------------------------------------------------------------------
    def computeNDNI(self) -> numpy.ndarray:

        b1510 = self._getBand(1510)
        b1680 = self._getBand(1680)

        log1 = self._computeLog(1, b1510)
        log2 = self._computeLog(1, b1680)

        denom = log1 + log2
        denom[denom==0] = self._getMin(denom)

        vi = (log1 - log2) / denom

        return vi

    # -------------------------------------------------------------------------
    # computeNDVI
    # -------------------------------------------------------------------------
    def computeNDVI(self) -> numpy.ndarray:

        b660 = self._getBand(660)
        b835 = self._getBand(835)
        denom = b835 + b660
        denom[denom==0] = self._getMin(denom)
        vi = (b835 - b660) / denom
        return vi

    # -------------------------------------------------------------------------
    # computeNDWI
    # -------------------------------------------------------------------------
    def computeNDWI(self) -> numpy.ndarray:

        b858 = self._getBand(858)
        b1241 = self._getBand(1241)
        denom = b858 + b1241
        denom[denom==0] = self._getMin(denom)
        vi = (b858 - b1241) / denom
        return vi

    # -------------------------------------------------------------------------
    # computeNIRv
    # -------------------------------------------------------------------------
    def computeNIRv(self) -> numpy.ndarray:

        b645 = self._getBand(645)
        b858 = self._getBand(858)
        denom = (b858 + b645) * b858
        denom[denom==0] = self._getMin(denom)
        vi = (b858 - b645) / denom
        return vi

    # -------------------------------------------------------------------------
    # computePRIn
    # -------------------------------------------------------------------------
    def computePRIn(self) -> numpy.ndarray:

        b531 = self._getBand(531)
        b550 = self._getBand(550)
        denom = b531 + b550
        denom[denom==0] = self._getMin(denom)
        vi = (b531 - b550) / denom
        return vi

    # -------------------------------------------------------------------------
    # computePRIw
    # -------------------------------------------------------------------------
    def computePRIw(self) -> numpy.ndarray:

        b531 = self._getBand(531)
        b570 = self._getBand(570)
        denom = b531 + b570
        denom[denom==0] = self._getMin(denom)
        vi = (b531 - b570) / denom
        return vi

    # -------------------------------------------------------------------------
    # computeSAVI
    #
    # [(835 - 655)/(835 + 655 + L)]* (1 + L) where L represents the fractional
    # amount of green vegetation cover and is normally set to 0.5 as an
    # approximation
    # -------------------------------------------------------------------------
    def computeSAVI(self) -> numpy.ndarray:

        L = 0.5
        b655 = self._getBand(655)
        b835 = self._getBand(835)
        denom = b835 + b655 + L
        denom[denom==0] = self._getMin(denom)
        vi = ((b835 - b655) / denom) * (1 + L)
        return vi

    # -------------------------------------------------------------------------
    # computeWBI
    # -------------------------------------------------------------------------
    def computeWBI(self) -> numpy.ndarray:

        b900 = self._getBand(900)
        b970 = self._getBand(970)
        denom = b970
        denom[denom==0] = self._getMin(denom)
        vi = b900 / denom
        return vi

    # -------------------------------------------------------------------------
    # _createDs
    # -------------------------------------------------------------------------
    def _createDs(self, numBands: int) -> gdal.Dataset:

        outDs = gdal.GetDriverByName('GTiff').Create( \
            str(self._outFileName),
            self._image.getDataset().RasterXSize,
            self._image.getDataset().RasterXSize,
            numBands,
            gdal.GDT_Float32,
            options=['COMPRESS=LZW', 'BIGTIFF=YES'])

        outDs.SetSpatialRef(self._image.srs())
        outDs.SetGeoTransform(self._image.getDataset().GetGeoTransform())
        return outDs

    # -------------------------------------------------------------------------
    # _getBand
    # -------------------------------------------------------------------------
    def _getBand(self, bandNum: int) -> numpy.ndarray:

        if not bandNum in self._bands:

            bandIndex = self._bandIndicies[bandNum]

            band = self._image.getDataset(). \
                   ReadAsArray(band_list=[bandIndex]).astype('float')

            bandDs = self._image.getDataset().GetRasterBand(bandIndex)
            band[band==bandDs.GetNoDataValue()] = numpy.nan
            scaleFactor = float(bandDs.GetMetadataItem('Scale factor'))
            self._bands[bandNum] = band / scaleFactor

        return self._bands[bandNum]

    # -------------------------------------------------------------------------
    # _getMin
    # -------------------------------------------------------------------------
    def _getMin(self, array: numpy.ndarray):

        if array.dtype == numpy.int16:
            return numpy.iinfo(numpy.int16).min

        elif array.dtype == numpy.float64:
            return numpy.finfo(numpy.float64).min

        else:
            raise RuntimeError('Unexpected data type.')

    # -------------------------------------------------------------------------
    # _writeVi
    # -------------------------------------------------------------------------
    def _writeVi(self,
                 vi: numpy.ndarray,
                 suffix: str,
                 ds: gdal.Dataset,
                 bandNum: int) -> int:

        bandNum += 1
        band = ds.GetRasterBand(bandNum)

        if not band:

            raise RuntimeError('Ensure ViHyper.computeAllAndWrite ' +
                               'created enough bands.')

        vi[vi==numpy.nan] = -10001
        band.WriteArray(vi)
        band.SetMetadata({'Index Name': suffix,})
        return bandNum
