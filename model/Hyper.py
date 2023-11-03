
from collections import namedtuple
from pathlib import Path

import logging
import os
import h5py

from osgeo import gdal
from osgeo import gdalconst
from osgeo.osr import SpatialReference

from shapely.geometry import box

from agb.model.BaseProcess import BaseProcess
from agb.model.NeonSite import NeonSite
from agb.model.ViHyper import ViHyper
from core.model.Envelope import Envelope
from core.model.GeospatialImageFile import GeospatialImageFile


# -----------------------------------------------------------------------------
# class Hyper
#
# TODO:  Can some of the GeospatialInfo attributes be removed?
# -----------------------------------------------------------------------------
class Hyper(BaseProcess):

    GeospatialInfo = namedtuple('GeospatialInfo',
                                'rows, \
                                 cols, \
                                 xMin, \
                                 yMin, \
                                 xMax, \
                                 yMax, \
                                 xScale, \
                                 yScale, \
                                 srs, \
                                 xform, \
                                 wavelengths, \
                                 scaleFactor, \
                                 noDataValue, \
                                 refl')

    # -------------------------------------------------------------------------
    # __init__
    # -------------------------------------------------------------------------
    def __init__(self,
                 site: NeonSite,
                 outDir: Path,
                 logger: logging.RootLogger,
                 year: int):

        # super(Hyper, self).__init__(site, outDir, logger)
        super().__init__(site, outDir, logger)
        self._year = year

        # ---
        # Create a working directory for intermediate files, so they can be
        # reused.
        # ---
        self._workDir = self._outDir / 'hyperTempFiles'

        if not self._workDir.exists():
            os.mkdir(self._workDir)
            
    # -------------------------------------------------------------------------
    # clipOne
    # -------------------------------------------------------------------------
    def _clipOne(self,
                 mosaic: gdal.Dataset,
                 envelope: Envelope,
                 clippedName: Path) -> None:

        # The target SRS is that of the site.
        cEnv = self._reprojEnv(mosaic.GetSpatialRef(), envelope)
        pWin = [cEnv.ulx(), cEnv.uly(), cEnv.lrx(), cEnv.lry()]
        tOpts = gdal.TranslateOptions(projWin=pWin)
        gdal.Translate(str(clippedName), mosaic, options=tOpts)

    # -------------------------------------------------------------------------
    # combine
    # -------------------------------------------------------------------------
    def _combine(self,
                 h5Name: Path,
                 albName: Path,
                 geoInfo: GeospatialInfo) -> Path:

        # ---
        # Add the albedos to the h5 file.  Cannot add a band to the Geotiff
        # driver, so we must create a new file with an extra band for Albedo.
        # ---
        combinedName = self._workDir / '05-hyper-albedo-clipped.tif'
        h5Ds = gdal.Open(str(h5Name), gdalconst.GA_Update)

        combinedDs = gdal.GetDriverByName('GTiff').Create(
            str(combinedName),
            h5Ds.RasterXSize,
            h5Ds.RasterYSize,
            h5Ds.RasterCount+1,
            gdal.GDT_Float32)

        combinedDs.SetSpatialRef(geoInfo.srs)
        combinedDs.SetGeoTransform(h5Ds.GetGeoTransform())

        # Add the h5 bands.
        bandIndicies = self._getBandIndicies(geoInfo)
        outBandIndex = 0

        for bi in bandIndicies:

            outBandIndex += 1
            band = combinedDs.GetRasterBand(outBandIndex)
            npPixels = h5Ds.GetRasterBand(outBandIndex).ReadAsArray()
            band.WriteArray(npPixels)
            band.SetNoDataValue(geoInfo.noDataValue)

            band.SetMetadata({'Wavelength': geoInfo.wavelengths[bi],
                              'Matched wavelength': bandIndicies[bi][0],
                              'Scale factor': geoInfo.scaleFactor})

            band.FlushCache()
            band = None

        # Add the albedo band.
        albData = gdal.Open(str(albName)).ReadAsArray()
        outBandIndex += 1
        band = combinedDs.GetRasterBand(outBandIndex)
        band.WriteArray(albData)
        band.SetMetadata({'Name': 'Albedo'})
        del combinedDs

        return combinedName

    # -------------------------------------------------------------------------
    # defaultInputDir
    # -------------------------------------------------------------------------
    @staticmethod
    def defaultInputDir():

        return Path('/explore/nobackup/projects/ilab/data/AGB/' \
                    'Airborne_Hyperspectral')

    # -------------------------------------------------------------------------
    # findFiles
    #
    # 2015_NEON_D07_MLBS_DP3_541000_4136000_reflectance.h5
    # year      sub site
    # -------------------------------------------------------------------------
    def findFiles(self, searchDir: Path = None, ext: str = 'h5') -> list:

        if not self._year:
            raise RuntimeError('Please specify a year to process.')

        searchDir = searchDir or self.defaultInputDir() / 'DP3.30006.001'

        files = list(searchDir.rglob('*' +
                                     self._site.abbreviation +
                                     '*.' +
                                     ext))

        filtered = \
            filter(lambda f: f.stem.split('_')[0] == str(self._year),
                   files)

        return list(filtered)

    # -------------------------------------------------------------------------
    # getAlbedos
    # -------------------------------------------------------------------------
    def _getAlbedos(self) -> gdal.Dataset:

        ALBEDO_DIR = self.defaultInputDir() / 'DP3.30011.001'
        albedos = self.findFiles(ALBEDO_DIR, ext='tif')

        # Remove files like 2017_NEON_D07_MLBS_DP3_Full_Albedo.
        filtered = list(filter(lambda f: f.stem.split('_')[5] != 'Full',
                               albedos))

        self._logger.info('Found ' + str(len(filtered)) + ' Albedos.')
        vrtName = self._workDir / '02-albedo.vrt'
        filtered = [str(f) for f in filtered]
        albedos = gdal.BuildVRT(str(vrtName), filtered)

        return albedos

    # -------------------------------------------------------------------------
    # getBandIndicies
    # -------------------------------------------------------------------------
    def _getBandIndicies(self, geoInfo: GeospatialInfo) -> dict:

        # ---
        # This is the list of wavelengths we need to calculate the vegetation
        # indicies.  Find each wavelength in geoInfo, and save its band index.
        # Only read these bands below.
        # ---
        soughtWLs = [485, 510, 515, 531, 550, 565, 570, 645, 655, 660, 681,
                     700, 709, 754, 790, 819, 835, 858, 900, 970, 1241, 1510,
                     1599, 1649, 1680, 1754, 2019, 2109, 2206]

        bandIndicies = {}

        for swl in soughtWLs:

            # ---
            # There might not be an exact match from sought wavelength and the
            # wavelengths we find.  Use the closest match.
            # ---
            bi = min(range(len(geoInfo.wavelengths)),
                     key=lambda i: abs(geoInfo.wavelengths[i] - swl))

            if bi not in bandIndicies:
                bandIndicies[bi] = []

            bandIndicies[bi].append(swl)

        # Ensure two sought wavelengths did not map to the same band.
        moreThanOne = list(filter(lambda i: len(bandIndicies[i]) > 1,
                           bandIndicies))

        if moreThanOne:

            raise RuntimeError('Multiple actual wavelengths, ' +
                               str(moreThanOne) +
                               ', were mapped to one equation wavelength, ' +
                               str(geoInfo.wavelengths[bi]) +
                               '.')

        return bandIndicies

    # -------------------------------------------------------------------------
    # getGeoInfo
    # -------------------------------------------------------------------------
    def _getGeoInfo(self, inFile: Path) -> GeospatialInfo:

        refl = self._getReflectanceDs(inFile)
        mapInfoRaw = refl['Metadata']['Coordinate_System']['Map_Info'][()]
        mapInfo = str(mapInfoRaw).split(',')
        reflShape = refl['Reflectance_Data'].shape

        xScale = float(mapInfo[5])
        yScale = float(mapInfo[6])

        xMin = float(mapInfo[3])
        yMax = float(mapInfo[4])
        xMax = xMin + (reflShape[1] * xScale)
        yMin = yMax - (reflShape[0] * yScale)

        rows = int(yMax - yMin)
        cols = int(xMax - xMin)

        xform = [xMin, xScale, 0, yMax, 0, -(yScale)]

        epsg = refl['Metadata']['Coordinate_System']['EPSG Code'][()].decode()
        srs = SpatialReference()
        srs.ImportFromEPSG(int(epsg))

        wlDs = refl['Metadata']['Spectral_Data']['Wavelength']
        wavelengths = [wl for wl in wlDs] 

        scaleFactor = refl['Reflectance_Data'].attrs['Scale_Factor']
        noDataValue = refl['Reflectance_Data'].attrs['Data_Ignore_Value']

        geoInfo = Hyper.GeospatialInfo(rows,
                                       cols,
                                       xMin,
                                       yMin,
                                       xMax,
                                       yMax,
                                       xScale,
                                       yScale,
                                       srs,
                                       xform,
                                       # wavelengths,
                                       wlDs,
                                       scaleFactor,
                                       noDataValue,
                                       refl)

        return geoInfo

    # -------------------------------------------------------------------------
    # getReflectanceDs
    # -------------------------------------------------------------------------
    def _getReflectanceDs(self, inFile: Path) -> h5py._hl.group.Group:

        h5File = h5py.File(inFile)
        refl = h5File[self._site.abbreviation]['Reflectance']
        return refl

    # -------------------------------------------------------------------------
    # _h5ToTif
    # -------------------------------------------------------------------------
    def _h5ToTif(self, inFile: Path) -> str:

        # ---
        # Extract the geospatial information from the H5 file all at once,
        # instead of redundantly navigating the structures.
        # ---
        geoInfo = self._getGeoInfo(inFile)

        # These are the "subdatasets" to extract.
        bandIndicies = self._getBandIndicies(geoInfo)
        subDataSets = ['Reflectance_Data']  # Always 1?
        tempTif = None

        for subDs in subDataSets:

            ds = geoInfo.refl[subDs]  # (1000, 1000, 426)

            tempTif = self._workDir / (inFile.stem + '.tif')

            if tempTif.exists():
                return str(tempTif)

            # No compression because it caused larger files.
            outDs = gdal.GetDriverByName('GTiff').Create(
                str(tempTif),
                int(geoInfo.xMax-geoInfo.xMin),
                int(geoInfo.yMax-geoInfo.yMin),
                # len(bandIndicies),
                ds.shape[2],  # Replaces len(bandIndicies)
                gdal.GDT_Int16,
                options=['BIGTIFF=YES'])

            outDs.SetSpatialRef(geoInfo.srs)
            outDs.SetGeoTransform(geoInfo.xform)

            outBandIndex = 0

            # for bi in bandIndicies:
            for i in range(ds.shape[2]):  # Replaces for bi in bandIndicies
            
                outBandIndex += 1
                band = outDs.GetRasterBand(outBandIndex)
                # npPixels = ds[:, :, bi]
                npPixels = ds[:, :, i]
                band.WriteArray(npPixels)
                band.SetNoDataValue(geoInfo.noDataValue)

                wl = geoInfo.wavelengths[i]
                matchWl = bandIndicies[wl][0] if wl in bandIndicies else 'n/a'
                
                band.SetMetadata({'Wavelength': wl,
                                  'Matched wavelength': matchWl,
                                  'Scale factor': geoInfo.scaleFactor})

                band.FlushCache()
                band = None

            outDs = None

        return str(tempTif)

    # -------------------------------------------------------------------------
    # intersectH5s
    # -------------------------------------------------------------------------
    def _intersectH5s(self, inFiles: list, env: Envelope) -> list:

        filtered = []
        envBox = box(env.ulx(), env.lry(), env.lrx(), env.uly())

        for f in inFiles:

            minX, minY = [int(n) for n in f.stem.split('_')[5:7]]
            inBox = box(minX, minY, minX+1000, minY+1000)

            if envBox.intersects(inBox):
                filtered.append(f)

        return filtered

    # -------------------------------------------------------------------------
    # myClip
    # -------------------------------------------------------------------------
    def _myClip(self,
                inFiles: list,
                envelope: Envelope,
                outFile: Path) -> None:

        # ---
        # Filter inFiles by the envelope to reduce H5s to read.
        # ---
        self._logger.info('Intersecting H5 files.')
        filteredFiles = self._intersectH5s(inFiles, envelope)
        self._logger.info(str(len(filteredFiles)) + ' intersecting H5s.')

        # ---
        # VRTs cannot be created from these H5 files.  We must open them using
        # H5 methods.
        # ---
        self._logger.info('Reading H5 files.')
        h5tifs = [self._h5ToTif(f) for f in filteredFiles]

        # Mosaic the H5 files using a VRT.
        self._logger.info('Mosaicking H5s.')
        h5MosName = self._workDir / '01-hyper-mosaic.vrt'
        h5Vrt = gdal.BuildVRT(str(h5MosName), h5tifs)

        # Find and mosaic the Albedo files using a VRT.
        self._logger.info('Mosaicking Albedos.')
        albedoDs = self._getAlbedos()

        # Clip H5s
        self._logger.info('Clipping H5s.')
        h5Name = self._workDir / '03-hyper-clipped.tif'
        self._clipOne(h5Vrt, envelope, h5Name)

        # Clip Albedo
        self._logger.info('Clipping Albedos.')
        albName = self._workDir / '04-albedo-clipped.tif'
        self._clipOne(albedoDs, envelope, albName)

        # ---
        # Albedo are independent images, and completely overlap the h5s.  If
        # they were added to the VRT, they would overwrite or be overwritten
        # by the h5s.  Instead, create a vrt from the albedos, clip it, then
        # add it as a band to the clipped h5s.
        # ---
        self._logger.info('Combining H5 bands and Albedo.')
        geoInfo = self._getGeoInfo(filteredFiles[0])
        combined = self._combine(h5Name, albName, geoInfo)

        # Compute the VIs.
        self._logger.info('Computing vegetation indicies.')
        gis = GeospatialImageFile(str(combined))
        vi = ViHyper(gis, outFile, self._logger)
        vi.computeAllAndWrite()

    # -------------------------------------------------------------------------
    # getOutFileName
    # -------------------------------------------------------------------------
    def _getOutFileName(self) -> str:

        return self._site.abbreviation + \
               '_' + str(self._year) + \
               '_Hyperspectral_Indicies'
               