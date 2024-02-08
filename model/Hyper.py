
from collections import namedtuple
import h5py
import logging
import os
from pathlib import Path
import statistics

from osgeo import gdal
from osgeo import gdalconst
from osgeo.osr import SpatialReference

import numpy as np
from scipy.signal import savgol_filter
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
class Hyper(object):

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

        # super().__init__(site, outDir, logger)

        if not outDir or not outDir.exists() or not outDir.is_dir():
            raise RuntimeError('Invalid output directory: ' + str(outDir))

        self._logger = logger
        self._outDir = outDir
        self._site = site

        self._year = year

        # ---
        # Create a working directory for intermediate files, so they can be
        # reused.
        # ---
        self._workDir = self._outDir / (site.abbreviation + '-' + str(year))

        if not self._workDir.exists():
            os.mkdir(self._workDir)
            
        # ---
        # GeoInfo objects are used in two places: h5ToTif and
        # createMetaDataImage.  Store them here, so they need not be created
        # twice.  Key = image file name, value = GeoInfo object
        # ---
        self._geoInfos = {}
        
    # -------------------------------------------------------------------------
    # clipOne
    # -------------------------------------------------------------------------
    # def _clipOne(self,
    #              mosaic: gdal.Dataset,
    #              envelope: Envelope,
    #              clippedName: Path) -> None:
    #
    #     # The target SRS is that of the site.
    #     cEnv = self._reprojEnv(mosaic.GetSpatialRef(), envelope)
    #     pWin = [cEnv.ulx(), cEnv.uly(), cEnv.lrx(), cEnv.lry()]
    #     tOpts = gdal.TranslateOptions(projWin=pWin)
    #     gdal.Translate(str(clippedName), mosaic, options=tOpts)

    # -------------------------------------------------------------------------
    # combine
    # -------------------------------------------------------------------------
    # def _combine(self,
    #              h5Name: Path,
    #              albName: Path,
    #              geoInfo: GeospatialInfo) -> Path:
    #
    #     # ---
    #     # Add the albedos to the h5 file.  Cannot add a band to the Geotiff
    #     # driver, so we must create a new file with an extra band for Albedo.
    #     # ---
    #     combinedName = self._workDir / '05-hyper-albedo-clipped.tif'
    #     h5Ds = gdal.Open(str(h5Name), gdalconst.GA_Update)
    #
    #     combinedDs = gdal.GetDriverByName('GTiff').Create(
    #         str(combinedName),
    #         h5Ds.RasterXSize,
    #         h5Ds.RasterYSize,
    #         h5Ds.RasterCount+1,
    #         gdal.GDT_Float32)
    #
    #     combinedDs.SetSpatialRef(geoInfo.srs)
    #     combinedDs.SetGeoTransform(h5Ds.GetGeoTransform())
    #
    #     # Add the h5 bands.
    #     bandIndices = self._getBandIndices(geoInfo)
    #     outBandIndex = 0
    #
    #     for bi in bandIndices:
    #
    #         outBandIndex += 1
    #         band = combinedDs.GetRasterBand(outBandIndex)
    #         npPixels = h5Ds.GetRasterBand(outBandIndex).ReadAsArray()
    #         band.WriteArray(npPixels)
    #         band.SetNoDataValue(geoInfo.noDataValue)
    #
    #         band.SetMetadata({'Wavelength': geoInfo.wavelengths[bi],
    #                           'Matched wavelength': bandIndices[bi][0],
    #                           'Scale factor': geoInfo.scaleFactor})
    #
    #         band.FlushCache()
    #         band = None
    #
    #     # Add the albedo band.
    #     albData = gdal.Open(str(albName)).ReadAsArray()
    #     outBandIndex += 1
    #     band = combinedDs.GetRasterBand(outBandIndex)
    #     band.WriteArray(albData)
    #     band.SetMetadata({'Name': 'Albedo'})
    #     del combinedDs
    #
    #     return combinedName

    # -------------------------------------------------------------------------
    # computeBRDF
    # -------------------------------------------------------------------------
    def _computeBRDF(self, geoInfo: GeospatialInfo) -> np.ndarray:

        # key1 = list(geoInfo.refl['Metadata/Logs'].keys())[0]
        # az1Ds = geoInfo.refl['Metadata/Logs/' + key1 + '/Solar_Azimuth_Angle']
        # solZnDs = geoInfo.refl['Metadata/Logs/' + key1 + '/Solar_Zenith_Angle']

        az1s = []
        solZns = []
        
        for key in geoInfo.refl['Metadata/Logs'].keys():
            
            az1s.append(geoInfo.refl['Metadata/Logs/' + 
                                     key + 
                                     '/Solar_Azimuth_Angle'][()])
        
            # Solar Zenith
            atcorLogText = str(geoInfo.refl['Metadata/Logs/' +
                                            key +
                                            '/ATCOR_Processing_Log'][()])
                                        
            atcors = atcorLogText.split('\\n')

            # ---
            # The log data was found to have extra spaces in the "Solar
            # Zenith angle" string.  This searches considering things like
            # that.
            # ---
            for atcor in atcors:
                
                noSpaces = atcor.replace(' ', '').lower()
                
                if 'solarzenithangle' in noSpaces:
                    solZns.append(float(noSpaces.split('=')[1]))
        
        az1 = statistics.mean(az1s)
        solZn = statistics.mean(solZns)
        az2 = geoInfo.refl['Metadata']['to-sensor_azimuth_angle'][:, :]
        senZn = geoInfo.refl['Metadata']['to-sensor_zenith_angle'][:, :]
        el1 = 90 - solZn 
        el2 = 90 - senZn
        
        # brdf = arccos[cosEl1 * cosEl2 * cos(Az1−Az2) + sinEl1 * sinEl2]
        cosEl1 = np.cos(el1)
        cosEl2 = np.cos(el2)
        cos12 = np.cos(az1 - az1)
        sinEl1 = np.sin(el1)
        sinEl2 = np.sin(el2)
        brdf = np.arccos(cosEl1 * cosEl2 * cos12 * sinEl1 * sinEl2)
        
        return brdf
        
    # -------------------------------------------------------------------------
    # createMetadataImage
    # -------------------------------------------------------------------------
    def _createMetadataImage(self, inFile: Path) -> None:

        metaFileName = self._workDir / (inFile.stem + '_ancillary.tif')
        
        if metaFileName.exists():
            
            self._logger.info('Metadata image, ' + 
                              str(metaFileName) + 
                              ', already exists.')
                              
            return

        geoInfo = self._getGeoInfo(inFile)
        anc = self._getReflectanceDs(inFile)['Metadata/Ancillary_Imagery']
        
        metaKeys = ['Aerosol_Optical_Depth',
                    'Aspect',
                    'Cast_Shadow',
                    'Dark_Dense_Vegetation_Classification',
                    'Data_Selection_Index',
                    'Haze_Cloud_Water_Map',
                    'Illumination_Factor',
                    'Path_Length',
                    'Sky_View_Factor','Slope',
                    'Smooth_Surface_Elevation',
                    'Visibility_Index_Map',
                    'Water_Vapor_Column',
                    'Weather_Quality_Indicator']
                     
        # ---
        # We cannot determine how many output bands there will be because
        # some of the meta keys have more than one band.  GDAL.Create needs
        # the number of output bands.  Collect all the output arrays first,
        # then write.
        # ---
        outBands = {}  # {band name: raster}
        
        for key in metaKeys:
            
            arr = anc[key]
            
            # ---
            # Most of these have a single band, but some can have more.  We
            # Must write all of them.  Paul adds a dimension to single-band
            # arrays, so he can loop through the third dimension in every
            # case.
            # ---
            if len(arr.shape) == 2:
                arr = np.expand_dims(arr, 2)
                
            for i in range(0, arr.shape[2]):
            
                saveKey = key + '_' + str(i + 1) if i > 0 else key
                outBands[saveKey] = arr[:, :, i]

        # Create the BRDF layer.
        outBands['Ang_diff_hotspot'] = self._computeBRDF(geoInfo)
        
        # Create the output dataset.
        outDs = gdal.GetDriverByName('GTiff').Create(
            str(metaFileName),
            int(geoInfo.xMax - geoInfo.xMin),
            int(geoInfo.yMax - geoInfo.yMin),
            len(outBands),  
            gdal.GDT_Float32,
            options=['BIGTIFF=YES', 'COMPRESS=LZW'])

        outDs.SetSpatialRef(geoInfo.srs)
        outDs.SetGeoTransform(geoInfo.xform)

        outBandIndex = 0

        for key in outBands:
            
            outBandIndex += 1
            outBand = outBands[key]
            band = outDs.GetRasterBand(outBandIndex)
            band.WriteArray(outBand)
            band.SetNoDataValue(geoInfo.noDataValue)
            band.SetMetadata({'Name': key})
            band.FlushCache()
            band = None

        outDs = None
        
    # -------------------------------------------------------------------------
    # defaultInputDir
    # -------------------------------------------------------------------------
    @staticmethod
    def defaultInputDir():

        return Path('/explore/nobackup/projects/ilab/data/AGB/' \
                    'Airborne_Hyperspectral')

    # -------------------------------------------------------------------------
    # filterTids
    # -------------------------------------------------------------------------
    def _filterTids(self, tids: list, inFiles: list):

        filtered = []

        for tid in tids:

            for f in inFiles:

                if tid in str(f):

                    filtered.append(f)
                    break

        return filtered
        
    # -------------------------------------------------------------------------
    # findFiles
    #
    # 2015_NEON_D07_MLBS_DP3_541000_4136000_reflectance.h5
    # year      sub site     subtile
    # -------------------------------------------------------------------------
    def findFiles(self, 
                  tids: list,
                  searchDir: Path = None, 
                  ext: str = 'h5') -> list:

        if not self._year:
            raise RuntimeError('Please specify a year to process.')

        searchDir = searchDir or self.defaultInputDir() / 'DP3.30006.001'

        files = list(searchDir.rglob('*' +
                                     self._site.abbreviation +
                                     '*.' +
                                     ext))

        files = \
            list(filter(lambda f: f.stem.split('_')[0] == str(self._year),
                   files))

        filtered = self._filterTids(tids, files)
        
        return filtered

    # -------------------------------------------------------------------------
    # getAlbedos
    # -------------------------------------------------------------------------
    # def _getAlbedos(self) -> gdal.Dataset:
    #
    #     ALBEDO_DIR = self.defaultInputDir() / 'DP3.30011.001'
    #     albedos = self.findFiles(ALBEDO_DIR, ext='tif')
    #
    #     # Remove files like 2017_NEON_D07_MLBS_DP3_Full_Albedo.
    #     filtered = list(filter(lambda f: f.stem.split('_')[5] != 'Full',
    #                            albedos))
    #
    #     self._logger.info('Found ' + str(len(filtered)) + ' Albedos.')
    #     vrtName = self._workDir / '02-albedo.vrt'
    #     filtered = [str(f) for f in filtered]
    #     albedos = gdal.BuildVRT(str(vrtName), filtered)
    #
    #     return albedos

    # -------------------------------------------------------------------------
    # getBandIndices
    # -------------------------------------------------------------------------
    def _getBandIndices(self, geoInfo: GeospatialInfo) -> dict:
            
        # ---
        # This is the list of wavelengths we need to calculate the vegetation
        # indices.  Find each wavelength in geoInfo, and save its band index.
        # Only read these bands below.
        # --- 
        soughtWLs = [485, 510, 515, 531, 550, 565, 570, 645, 655, 660, 681,
                     700, 709, 754, 790, 819, 835, 858, 900, 970, 1241, 1510,
                     1599, 1649, 1680, 1754, 2019, 2109, 2206]
        
        bandIndices = {}
            
        for swl in soughtWLs:
            
            # ---
            # There might not be an exact match from sought wavelength and the
            # wavelengths we find.  Use the closest match.
            # ---
            bi = min(range(len(geoInfo.wavelengths)),
                     key=lambda i: abs(geoInfo.wavelengths[i] - swl))

            if bi not in bandIndices:
                bandIndices[bi] = []

            bandIndices[bi].append(swl)

        # Ensure two sought wavelengths did not map to the same band.
        moreThanOne = list(filter(lambda i: len(bandIndices[i]) > 1,
                           bandIndices))

        if moreThanOne:

            raise RuntimeError('Multiple actual wavelengths, ' +
                               str(moreThanOne) +
                               ', were mapped to one equation wavelength, ' +
                               str(geoInfo.wavelengths[bi]) +
                               '.')

        if len(soughtWLs) > len(bandIndices):
            raise RuntimeError('Not all sought wavelengths were found.')
            
        return bandIndices

    # -------------------------------------------------------------------------
    # getGeoInfo
    # -------------------------------------------------------------------------
    def _getGeoInfo(self, inFile: Path) -> GeospatialInfo:

        if inFile in self._geoInfos:
            return self._geoInfos[inFile]

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
        import pdb
        pdb.set_trace()
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
                                       
        self._geoInfos[inFile] = geoInfo

        return geoInfo

    # -------------------------------------------------------------------------
    # getOutFileName
    # -------------------------------------------------------------------------
    def _getOutFileName(self) -> str:

        return self._site.abbreviation + \
               '_' + str(self._year) + \
               '_Hyperspectral_Indices'

    # -------------------------------------------------------------------------
    # getReflectanceDs
    # -------------------------------------------------------------------------
    def _getReflectanceDs(self, inFile: Path) -> h5py._hl.group.Group:

        h5File = h5py.File(inFile)
        refl = h5File[self._site.abbreviation]['Reflectance']
        return refl

    # -------------------------------------------------------------------------
    # _h5ToTif
    #
    # https://github.com/pahbs/geoscitools/blob/master/hyper_savgol_deriv_prod.py
    # -------------------------------------------------------------------------
    def _h5ToTif(self, inFile: Path) -> Path:

        reflName = self._workDir / (inFile.stem + '.tif')
        sgName = self._workDir / (inFile.stem + '_savgol.tif')
        
        if reflName.exists() and sgName.exists():
            return reflName

        # ---
        # Extract the geospatial information from the H5 file all at once,
        # instead of redundantly navigating the structures.
        # ---
        geoInfo = self._getGeoInfo(inFile)

        # These are the "subdatasets" to extract.
        bandIndices = self._getBandIndices(geoInfo)
        subDataSets = ['Reflectance_Data']  # Always 1?

        for subDs in subDataSets:

            ds = geoInfo.refl[subDs]  # (1000, 1000, 426)

            # ---
            # Create tif with hyperspectral indices.  Use no compression
            # because it caused larger files.
            # ---
            outDs = gdal.GetDriverByName('GTiff').Create(
                str(reflName),
                int(geoInfo.xMax-geoInfo.xMin),
                int(geoInfo.yMax-geoInfo.yMin),
                ds.shape[2],  
                gdal.GDT_Int16,
                options=['BIGTIFF=YES'])

            outDs.SetSpatialRef(geoInfo.srs)
            outDs.SetGeoTransform(geoInfo.xform)
            
            # ---
            # The savgol image is the same as the hyperspectral indices,
            # but run through the savgol function.  Both of these images
            # operate on the same set of input pixels, so write them
            # simultaneously.
            # ---
            sgDs = gdal.GetDriverByName('GTiff').CreateCopy(str(sgName), outDs)

            # Loop through the input.
            outBandIndex = 0

            for i in range(ds.shape[2]): 
            
                outBandIndex += 1
                
                # Write the hyperspectral index.
                band = outDs.GetRasterBand(outBandIndex)
                npPixels = ds[:, :, i]
                band.WriteArray(npPixels)
                band.SetNoDataValue(geoInfo.noDataValue)

                wl = bandIndices[i][0] if i in bandIndices else 'n/a'

                band.SetMetadata({'Wavelength': geoInfo.wavelengths[i],
                                  'Matched wavelength': wl,
                                  'Scale factor': geoInfo.scaleFactor})

                band.FlushCache()
                band = None
                
                # Write the savgol image.  Write NaNs for certain bands.
                if i <= 4 or \
                    (i >= 189 and i <= 216) or \
                    (i >= 280 and i <= 322) or \
                    i >= 416:
                
                    nanPixels = np.full_like(npPixels, np.nan)
                    npPixels = nanPixels
                    
                savgolFiltered = savgol_filter(x=npPixels,
                                               axis=0,
                                               window_length=7,
                                               polyorder=2,
                                               deriv=1,
                                               mode='nearest')
                
                sgBand = sgDs.GetRasterBand(outBandIndex)
                sgBand.WriteArray(savgolFiltered)
                sgBand.SetMetadata({'Wavelength': wl})
                sgBand.FlushCache()
                sgBand = None

            outDs = None
            sgDs = None

        return reflName

    # -------------------------------------------------------------------------
    # intersectH5s
    # -------------------------------------------------------------------------
    # def _intersectH5s(self, inFiles: list, env: Envelope) -> list:
    #
    #     filtered = []
    #     envBox = box(env.ulx(), env.lry(), env.lrx(), env.uly())
    #
    #     for f in inFiles:
    #
    #         minX, minY = [int(n) for n in f.stem.split('_')[5:7]]
    #         inBox = box(minX, minY, minX+1000, minY+1000)
    #
    #         if envBox.intersects(inBox):
    #             filtered.append(f)
    #
    #     return filtered

    # -------------------------------------------------------------------------
    # myClip
    #-------------------------------------------------------------------------
    # def _myClip(self,
    #             inFiles: list,
    #             envelope: Envelope,
    #             outFile: Path) -> None:
    #
    #     self._logger.info('Reading H5 files.')
    #     h5Tifs = []
    #
    #     for inFile in inFiles:
    #
    #         # *_reflectance.tif, *_savgol.tif
    #         h5Tifs.append(self._h5ToTif(inFile))
    #         self._createMetadataImage(inFile)
    #
    #     # ---
    #     # VRTs cannot be created from these H5 files.  We must open them using
    #     # H5 methods.
    #     # ---
    #     # self._logger.info('Reading H5 files.')
    #     # h5Tifs = [self._h5ToTif(f) for f in filteredFiles]
    #
    #     # Mosaic the H5 files using a VRT.
    #     self._logger.info('Mosaicking H5s.')
    #     h5MosName = self._workDir / '01-hyper-mosaic.vrt'
    #     h5TifsStr = [str(f) for f in filteredFiles]
    #     h5Vrt = gdal.BuildVRT(str(h5MosName), h5TifsStr)
    #
    #     # Find and mosaic the Albedo files using a VRT.
    #     self._logger.info('Mosaicking Albedos.')
    #     albedoDs = self._getAlbedos()
    #
    #     # Clip H5s
    #     self._logger.info('Clipping H5s.')
    #     h5Name = self._workDir / '03-hyper-clipped.tif'
    #     self._clipOne(h5Vrt, envelope, h5Name)
    #
    #     # Clip Albedo
    #     self._logger.info('Clipping Albedos.')
    #     albName = self._workDir / '04-albedo-clipped.tif'
    #     self._clipOne(albedoDs, envelope, albName)
    #
    #     # ---
    #     # Albedo are independent images, and completely overlap the h5s.  If
    #     # they were added to the VRT, they would overwrite or be overwritten
    #     # by the h5s.  Instead, create a vrt from the albedos, clip it, then
    #     # add it as a band to the clipped h5s.
    #     # ---
    #     self._logger.info('Combining H5 bands and Albedo.')
    #     geoInfo = self._getGeoInfo(filteredFiles[0])
    #     combined = self._combine(h5Name, albName, geoInfo)
    #
    #     # Compute the VIs.
    #     self._logger.info('Computing vegetation indices.')
    #     gis = GeospatialImageFile(str(combined))
    #     vi = ViHyper(gis, outFile, self._logger)
    #     vi.computeAllAndWrite()

    # -------------------------------------------------------------------------
    # myClip
    # -------------------------------------------------------------------------
    def _myClip(self, inFiles: list) -> None:

        for inFile in inFiles:

            self._logger.info('Processing ' + str(inFile))
            
            # Get the albedo.
            strFile = str(inFile)
            
            albName = Path(strFile.replace('30006', '30011'). \
                                   replace('Reflectance', 'Albedo'). \
                                   replace('reflectance', 'albedo'). \
                                   replace('h5', 'tif'))
                                  
            if not albName.exists():
                
                raise RuntimeError('Albedo file, ' + \
                                    str(albName) + \
                                    ' does not exist.')
                                    
            # Albedo goes in with the VIs.  Pass it as a Numpy array.
            self._logger.info('Reading ' + str(albName))
            albedo = gdal.Open(str(albName)).ReadAsArray().astype('float')

            # *_reflectance.tif, *_savgol.tif
            self._logger.info('Writing reflectance and savgol.')
            h5Tif = self._h5ToTif(inFile)

            # *_ancillary.tif
            self._logger.info('Writing ancillary. ')
            self._createMetadataImage(inFile)
            
            # *_indices.tif
            self._logger.info('Computing vegetation indices.')
            gis = GeospatialImageFile(str(h5Tif))
            outFile = self._workDir / (inFile.stem + '_indices.tif')
            vi = ViHyper(gis, albedo, outFile, self._logger)
            vi.computeAllAndWrite()

    # -------------------------------------------------------------------------
    # run
    # -------------------------------------------------------------------------
    def run(self, tids: list = None) -> None:

        inFiles = self.findFiles(tids)
        numFiles = len(inFiles)
        self._logger.info('Found ' + str(numFiles) + ' files.')

        if numFiles == 0:
            self._logger.info('Quitting')

        self._myClip(inFiles)
        