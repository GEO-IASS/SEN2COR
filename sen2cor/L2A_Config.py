#!/usr/bin/env python

from numpy import *
import fnmatch
import sys, os, logging, inspect
import ConfigParser
from L2A_XmlParser import L2A_XmlParser
from L2A_Library import stdoutWrite, stderrWrite
from lxml import etree, objectify
from time import strftime
from datetime import datetime
from multiprocessing import Lock
from shutil import copyfile, copytree
from psutil import cpu_count
l = Lock()


def getScriptDir(follow_symlinks=True):
    if getattr(sys, 'frozen', False): # py2exe, PyInstaller, cx_Freeze
        path = os.path.abspath(sys.executable)
    else:
        path = inspect.getabsfile(getScriptDir)
    if follow_symlinks:
        path = os.path.realpath(path)
    return os.path.dirname(path)


class L2A_Config(object):
    
    def __init__(self, logger, sourceDir = False):
        self._processorName = 'Sentinel-2 Level 2A Prototype Processor (Sen2Cor)'
        self._processorVersion = '2.2.0'
        self._processorDate = '2016.03.31'
        self._productVersion = '13.1'
        self._logger = logger
        self._logLevel = 'INFO'
        self._sc_lp_blu = 1.0
        self._tEstimation = 0.0

        if(sourceDir):
            try:
                self._home = os.environ['SEN2COR_HOME']
            except:
                self._home = os.path.dirname(getScriptDir())               
            try:
                scriptDir = os.environ['SEN2COR_BIN']
            except:
                scriptDir = getScriptDir()    

            self._sourceDir = sourceDir
            self._configDir = os.path.join(scriptDir, 'cfg')
            self._logDir = os.path.join(self._home, 'log')
            self._configFn = os.path.join(self._home, 'cfg', 'L2A_GIPP.xml')
            self.configSC = os.path.join(self._configDir, 'L2A_CAL_SC_GIPP.xml')
            self.configAC = os.path.join(self._configDir, 'L2A_CAL_AC_GIPP.xml')
            
            self._atmDataFn = None
            self._tEst60 = 150.0
            self._tEst20 = self._tEst60 * 8.0
            self._tEst10 = self._tEst60 * 8.0
            self._tStart = None
            self._nrTiles = None
            self._nrProcs = None
            self._processingStatusFn = os.path.join(self._logDir, '.progress')
            self._processingEstimationFn = os.path.join(self._configDir, '.estimation')
            
            if os.path.isfile(self._processingEstimationFn) == False:
            # init processing estimation file:
                config = ConfigParser.RawConfigParser()
                config.add_section('time estimation')
                config.set('time estimation','t_est_60', self._tEst60)
                config.set('time estimation','t_est_20', self._tEst20)
                config.set('time estimation','t_est_10', self._tEst10)
                configFile = open(self._processingEstimationFn, 'w')
                config.write(configFile)
                configFile.close()
            self._ncols = -1
            self._nrows = -1
            self._nbnds = -1
            self._sza = -1
            self._saa = -1
            self._GIPP = ''
            self._ECMWF = ''
            self._DEM = ''
            self._L2A_BOA_QUANTIFICATION_VALUE = None
            self._L2A_WVP_QUANTIFICATION_VALUE = 1000.0
            self._L2A_AOT_QUANTIFICATION_VALUE = 1000.0
            self._dnScale = 1000.0
            self._adj_km = 1.0
            self._ch940 = array([8,8,9,9,0,0])
            self._cellsize = 0 # pixelsize (m), cellsize (km)
            self._dem_unit = 0 # [meter] is default DEM heigh unit unit
            self._ibrdf = 0         # brdf correction
            self._thr_g = 0.25      # lower bound for brdf correction
            self._ibrdf_dark = 0    # flag for cast shadow dark pixel processing
            self._icl_shadow = 0    # no cloud shadow calculation
            self._iclshad_mask = 3  # 1=small, 2=medium, 3=large mask with default values of thr_shad
            self._ihaze = 0
            self._ihcw = 1
            self._ihot_dynr = 2
            self._ihot_mask = 2
            self._intpol760 = 1         # always 1 for Sentinel2
            self._intpol725_825 = 1     # always 1 for Sentinel2
            self._intpol1400 = 1        # always 1 for Sentinel2
            self._intpol940_1130 = 1    # always 1 for Sentinel2
            self._istretch_type = 1     # linear
            self._iwat_shd = 0
            self._iwaterwv = 1
            self._iwv_ndvi = 0
            self._iwv_watermask = 1 # (1=average wv is used for water pixels)
            self._ksolflux = 0
            self._altit = 0.1
            self._npref = 0
            self._phi_scl_min = 0.05 # lower bound of shadow fraction, default = 0.05
            self._phi_unscl_max = -999.0 # initialize in Common/atcor_main
            self._pixelsize = 0
            self._resolution = 0
            self._rel_saturation = 1.00
            self._thr_shad = -999.0
            self._thv = 0.0   # tilt view angle, tilt azimuth angle
            self._phiv = 90.0 # sensor view zenith and azimuth angle
            self._smooth_wvmap = 100.0
            self._entityId = ''
            self._acquisitionDate = ''
            self._orbitPath = None
            self._orbitRow = None
            self._targetPath = None
            self._targetRow = None
            self._stationSgs = ''
            self._sceneStartTime = None
            self._sceneStopTime = None
            self._solaz = None
            self._solaz_arr = None
            self._solze = None
            self._solze_arr = None
            self._vaa_arr = None
            self._vza_arr = None
            self._visibility = 23.0
            self._wl940a = array([0.895, 1.000])     # range of moderate wv absorption region around  940 nm
            self._wl1130a = array([1.079, 1.180])    # range of moderate wv absorption region around 1130 nm
            self._wl1400a = array([1.330, 1.490])    # range for interpolation
            self._wl1900a = array([1.780, 1.970])    # range for interpolation
            self._wv_thr_cirrus = 0.60
            self._d2 = 1.0 # this is a constant, see explanation in line 2061
            self._timestamp = None
            self._c0 = None
            self._c1 = None
            self._e0 = None
            self._wvlsen = None
            self._fwhm = None
            self._acOnly = False
            self._L2A_INSPIRE_XML = None
            self._L2A_MANIFEST_SAFE = None
            self._L1C_UP_MTD_XML = None
            self._L1C_DS_MTD_XML = None
            self._L1C_TILE_MTD_XML = None
            self._L1C_UP_ID = None
            self._L1C_DS_ID = None
            self._L1C_TILE_ID = None
            self._L2A_UP_MTD_XML = None
            self._L2A_DS_MTD_XML = None
            self._L2A_DS_MTD_LST = []            
            self._L2A_TILE_MTD_XML = None
            self._L2A_UP_ID = None
            self._L2A_DS_LST = []
            self._L2A_DS_ID = None
            self._L2A_TILE_ID = None
            self._creationDate = None
            self._targetDirectory = None
            self._scOnly = False
            self._workDir = None
            self._processed60 = False
            self._processed20 = False
            self._processed10 = False
            self._refresh = False
            self._selectedTile = None
            self._aerosolType = None
            self._midLatitude = None
            self._ozoneContent = None
            self._waterVapour = None    
            self._ozoneMean = None        
            return

    def get_ozone_mean(self):
        return self._ozoneMean


    def set_ozone_mean(self, value):
        self._ozoneMean = value


    def del_ozone_mean(self):
        del self._ozoneMean


    def get_aerosol_type(self):
        return self._aerosolType


    def get_mid_latitude(self):
        return self._midLatitude


    def get_ozone_content(self):
        return self._ozoneContent


    def get_water_vapour(self):
        return self._waterVapour


    def set_aerosol_type(self, value):
        self._aerosolType = value


    def set_mid_latitude(self, value):
        self._midLatitude = value


    def set_ozone_content(self, value):
        self._ozoneContent = value


    def set_water_vapour(self, value):
        self._waterVapour = value


    def del_aerosol_type(self):
        del self._aerosolType


    def del_mid_latitude(self):
        del self._midLatitude


    def del_ozone_content(self):
        del self._ozoneContent


    def del_water_vapour(self):
        del self._waterVapour


    def get_processing_status_fn(self):
        return self._processingStatusFn


    def get_processing_estimation_fn(self):
        return self._processingEstimationFn


    def set_processing_status_fn(self, value):
        self._processingStatusFn = value


    def set_processing_estimation_fn(self, value):
        self._processingEstimationFn = value


    def del_processing_status_fn(self):
        del self._processingStatusFn


    def del_processing_estimation_fn(self):
        del self._processingEstimationFn



    def get_nr_tiles(self):
        if self.selectedTile != None:
            self._nrTiles = 1
            return self._nrTiles
        # else:    
        self._nrTiles = 0
        filemask = 'S2A_*_L1C_*'
        GRANULE = 'GRANULE'
        granuleDir = os.path.join(self.workDir, GRANULE)
        filelist = sorted(os.listdir(granuleDir))
        for tile in filelist:      
            if fnmatch.fnmatch(tile, filemask) == False:
                continue
            self._nrTiles += 1
        return self._nrTiles


    def get_nr_procs(self):
        nrProcs = self._nrProcs
        if nrProcs == 'AUTO':
            self._nrProcs = cpu_count()
        else:
            self._nrProcs = int(nrProcs)
        return self._nrProcs        


    def set_nr_tiles(self, value):
        self._nrTiles = value


    def set_nr_procs(self, value):
        self._nrProcs = value


    def del_nr_tiles(self):
        del self._nrTiles


    def del_nr_procs(self):
        del self._nrProcs


    def get_t_start(self):
        return self._tStart


    def set_t_start(self, value):
        self._tStart = value


    def del_t_start(self):
        del self._tStart


    def get_t_estimation(self):
        return self._tEstimation


    def set_t_estimation(self, value):
        self._tEstimation = value


    def del_t_estimation(self):
        del self._tEstimation


    def get_selected_tile(self):
        return self._selectedTile


    def set_selected_tile(self, value):
        self._selectedTile = value


    def del_selected_tile(self):
        del self._selectedTile


    def get_log_level(self):
        return self._logLevel


    def set_log_level(self, value):
        self._logLevel = value


    def del_log_level(self):
        del self._logLevel


    def get_target_directory(self):
        return self._targetDirectory


    def set_target_directory(self, value):
        self._targetDirectory = value


    def del_target_directory(self):
        del self._targetDirectory


    def get_logger(self):
        return self._logger
 
 
    def set_logger(self, value):
        self._logger = value
 
 
    def del_logger(self):
        del self._logger


    def get_sc_lp_blu(self):
        return self._sc_lp_blu


    def set_sc_lp_blu(self, value):
        self._sc_lp_blu = value


    def del_sc_lp_blu(self):
        del self._sc_lp_blu


    def get_refresh(self):
        return self._refresh


    def set_refresh(self, value):
        self._refresh = value


    def del_refresh(self):
        del self._refresh


    def get_processed_60(self):
        return self._processed60


    def get_processed_20(self):
        return self._processed20


    def get_processed_10(self):
        return self._processed10


    def set_processed_60(self, value):
        self._processed60 = value


    def set_processed_20(self, value):
        self._processed20 = value


    def set_processed_10(self, value):
        self._processed10 = value


    def del_processed_60(self):
        del self._processed60


    def del_processed_20(self):
        del self._processed20


    def del_processed_10(self):
        del self._processed10


    def get_work_dir(self):
        return self._workDir


    def set_work_dir(self, value):
        self._workDir = value


    def del_work_dir(self):
        del self._workDir


    def get_sc_only(self):
        return self._scOnly


    def set_sc_only(self, value):
        self._scOnly = value


    def del_sc_only(self):
        del self._scOnly


    def get_config_sc(self):
        return self._configSC


    def get_config_ac(self):
        return self._configAC


    def set_config_sc(self, value):
        self._configSC = value


    def set_config_ac(self, value):
        self._configAC = value


    def del_config_sc(self):
        del self._configSC


    def del_config_ac(self):
        del self._configAC


    def get_dem_directory(self):
        return self._demDirectory


    def get_dem_reference(self):
        return self._demReference


    def get_snow_map_reference(self):
        return self._snowMapReference


    def set_dem_directory(self, value):
        self._demDirectory = value


    def set_dem_reference(self, value):
        self._demReference = value


    def set_snow_map_reference(self, value):
        self._snowMapReference = value


    def del_dem_directory(self):
        del self._demDirectory


    def del_dem_reference(self):
        del self._demReference


    def del_snow_map_reference(self):
        del self._snowMapReference


    def get_median_filter(self):
        return self._medianFilter


    def set_median_filter(self, value):
        self._medianFilter = value


    def del_median_filter(self):
        del self._medianFilter


    def get_t_sdw(self):
        return self._T_SDW


    def get_t_b_02_b_12(self):
        return self._T_B02_B12


    def set_t_sdw(self, value):
        self._T_SDW = value


    def set_t_b_02_b_12(self, value):
        self._T_B02_B12 = value


    def del_t_sdw(self):
        del self._T_SDW


    def del_t_b_02_b_12(self):
        del self._T_B02_B12


    def get_no_data(self):
        return self._noData


    def get_saturated_defective(self):
        return self._saturatedDefective


    def get_dark_features(self):
        return self._darkFeatures


    def get_cloud_shadows(self):
        return self._cloudShadows


    def get_vegetation(self):
        return self._vegetation


    def get_bare_soils(self):
        return self._bareSoils


    def get_water(self):
        return self._water


    def get_low_proba_clouds(self):
        return self._lowProbaClouds


    def get_med_proba_clouds(self):
        return self._medProbaClouds


    def get_high_proba_clouds(self):
        return self._highProbaClouds


    def get_thin_cirrus(self):
        return self._thinCirrus


    def get_snow_ice(self):
        return self._snowIce


    def set_no_data(self, value):
        self._noData = value


    def set_saturated_defective(self, value):
        self._saturatedDefective = value


    def set_dark_features(self, value):
        self._darkFeatures = value


    def set_cloud_shadows(self, value):
        self._cloudShadows = value


    def set_vegetation(self, value):
        self._vegetation = value


    def set_bare_soils(self, value):
        self._bareSoils = value


    def set_water(self, value):
        self._water = value


    def set_low_proba_clouds(self, value):
        self._lowProbaClouds = value


    def set_med_proba_clouds(self, value):
        self._medProbaClouds = value


    def set_high_proba_clouds(self, value):
        self._highProbaClouds = value


    def set_thin_cirrus(self, value):
        self._thinCirrus = value


    def set_snow_ice(self, value):
        self._snowIce = value


    def del_no_data(self):
        del self._noData


    def del_saturated_defective(self):
        del self._saturatedDefective


    def del_dark_features(self):
        del self._darkFeatures


    def del_cloud_shadows(self):
        del self._cloudShadows


    def del_vegetation(self):
        del self._vegetation


    def del_bare_soils(self):
        del self._bareSoils


    def del_water(self):
        del self._water


    def del_low_proba_clouds(self):
        del self._lowProbaClouds


    def del_med_proba_clouds(self):
        del self._medProbaClouds


    def del_high_proba_clouds(self):
        del self._highProbaClouds


    def del_thin_cirrus(self):
        del self._thinCirrus


    def del_snow_ice(self):
        del self._snowIce


    def get_t_21_b_12(self):
        return self._T21_B12


    def get_t_22_b_12(self):
        return self._T22_B12


    def get_t_21_r_b_02_b_11(self):
        return self._T21_R_B02_B11


    def get_t_22_r_b_02_b_11(self):
        return self._T22_R_B02_B11


    def set_t_21_b_12(self, value):
        self._T21_B12 = value


    def set_t_22_b_12(self, value):
        self._T22_B12 = value


    def set_t_21_r_b_02_b_11(self, value):
        self._T21_R_B02_B11 = value


    def set_t_22_r_b_02_b_11(self, value):
        self._T22_R_B02_B11 = value


    def del_t_21_b_12(self):
        del self._T21_B12


    def del_t_22_b_12(self):
        del self._T22_B12


    def del_t_21_r_b_02_b_11(self):
        del self._T21_R_B02_B11


    def del_t_22_r_b_02_b_11(self):
        del self._T22_R_B02_B11


    def get_t_11_b_02(self):
        return self._T11_B02


    def get_t_12_b_02(self):
        return self._T12_B02


    def get_t_11_r_b_02_b_11(self):
        return self._T11_R_B02_B11


    def get_t_12_r_b_02_b_11(self):
        return self._T12_R_B02_B11


    def set_t_11_b_02(self, value):
        self._T11_B02 = value


    def set_t_12_b_02(self, value):
        self._T12_B02 = value


    def set_t_11_r_b_02_b_11(self, value):
        self._T11_R_B02_B11 = value


    def set_t_12_r_b_02_b_11(self, value):
        self._T12_R_B02_B11 = value


    def del_t_11_b_02(self):
        del self._T11_B02


    def del_t_12_b_02(self):
        del self._T12_B02


    def del_t_11_r_b_02_b_11(self):
        del self._T11_R_B02_B11


    def del_t_12_r_b_02_b_11(self):
        del self._T12_R_B02_B11


    def get_t_1_ndvi(self):
        return self._T1_NDVI


    def get_t_2_ndvi(self):
        return self._T2_NDVI


    def set_t_1_ndvi(self, value):
        self._T1_NDVI = value


    def set_t_2_ndvi(self, value):
        self._T2_NDVI = value


    def del_t_1_ndvi(self):
        del self._T1_NDVI


    def del_t_2_ndvi(self):
        del self._T2_NDVI


    def get_t_1_snow(self):
        return self._T1_SNOW


    def get_t_2_snow(self):
        return self._T2_SNOW


    def set_t_1_snow(self, value):
        self._T1_SNOW = value


    def set_t_2_snow(self, value):
        self._T2_SNOW = value


    def del_t_1_snow(self):
        del self._T1_SNOW


    def del_t_2_snow(self):
        del self._T2_SNOW


    def get_t_1_r_b_02_b_04(self):
        return self._T1_R_B02_B04


    def get_t_2_r_b_02_b_04(self):
        return self._T2_R_B02_B04


    def get_t_1_r_b_8_a_b_03(self):
        return self._T1_R_B8A_B03


    def get_t_2_r_b_8_a_b_03(self):
        return self._T2_R_B8A_B03


    def get_t_1_r_b_8_a_b_11(self):
        return self._T1_R_B8A_B11


    def get_t_2_r_b_8_a_b_11(self):
        return self._T2_R_B8A_B11


    def set_t_1_r_b_02_b_04(self, value):
        self._T1_R_B02_B04 = value


    def set_t_2_r_b_02_b_04(self, value):
        self._T2_R_B02_B04 = value


    def set_t_1_r_b_8_a_b_03(self, value):
        self._T1_R_B8A_B03 = value


    def set_t_2_r_b_8_a_b_03(self, value):
        self._T2_R_B8A_B03 = value


    def set_t_1_r_b_8_a_b_11(self, value):
        self._T1_R_B8A_B11 = value


    def set_t_2_r_b_8_a_b_11(self, value):
        self._T2_R_B8A_B11 = value


    def del_t_1_r_b_02_b_04(self):
        del self._T1_R_B02_B04


    def del_t_2_r_b_02_b_04(self):
        del self._T2_R_B02_B04


    def del_t_1_r_b_8_a_b_03(self):
        del self._T1_R_B8A_B03


    def del_t_2_r_b_8_a_b_03(self):
        del self._T2_R_B8A_B03


    def del_t_1_r_b_8_a_b_11(self):
        del self._T1_R_B8A_B11


    def del_t_2_r_b_8_a_b_11(self):
        del self._T2_R_B8A_B11


    def get_t_1_b_02(self):
        return self._T1_B02


    def get_t_2_b_02(self):
        return self._T2_B02


    def get_t_1_b_8_a(self):
        return self._T1_B8A


    def get_t_2_b_8_a(self):
        return self._T2_B8A


    def get_t_1_b_10(self):
        return self._T1_B10


    def get_t_2_b_10(self):
        return self._T2_B10


    def get_t_1_b_12(self):
        return self._T1_B12


    def get_t_2_b_12(self):
        return self._T2_B12


    def set_t_1_b_02(self, value):
        self._T1_B02 = value


    def set_t_2_b_02(self, value):
        self._T2_B02 = value


    def set_t_1_b_8_a(self, value):
        self._T1_B8A = value


    def set_t_2_b_8_a(self, value):
        self._T2_B8A = value


    def set_t_1_b_10(self, value):
        self._T1_B10 = value


    def set_t_2_b_10(self, value):
        self._T2_B10 = value


    def set_t_1_b_12(self, value):
        self._T1_B12 = value


    def set_t_2_b_12(self, value):
        self._T2_B12 = value


    def del_t_1_b_02(self):
        del self._T1_B02


    def del_t_2_b_02(self):
        del self._T2_B02


    def del_t_1_b_8_a(self):
        del self._T1_B8A


    def del_t_2_b_8_a(self):
        del self._T2_B8A


    def del_t_1_b_10(self):
        del self._T1_B10


    def del_t_2_b_10(self):
        del self._T2_B10


    def del_t_1_b_12(self):
        del self._T1_B12


    def del_t_2_b_12(self):
        del self._T2_B12


    def get_t_1_ndsi_snw(self):
        return self._T1_NDSI_SNW


    def get_t_2_ndsi_snw(self):
        return self._T2_NDSI_SNW


    def set_t_1_ndsi_snw(self, value):
        self._T1_NDSI_SNW = value


    def set_t_2_ndsi_snw(self, value):
        self._T2_NDSI_SNW = value


    def del_t_1_ndsi_snw(self):
        del self._T1_NDSI_SNW


    def del_t_2_ndsi_snw(self):
        del self._T2_NDSI_SNW


    def get_t_1_ndsi_cld(self):
        return self._T1_NDSI_CLD


    def get_t_2_ndsi_cld(self):
        return self._T2_NDSI_CLD


    def set_t_1_ndsi_cld(self, value):
        self._T1_NDSI_CLD = value


    def set_t_2_ndsi_cld(self, value):
        self._T2_NDSI_CLD = value


    def del_t_1_ndsi_cld(self):
        del self._T1_NDSI_CLD


    def del_t_2_ndsi_cld(self):
        del self._T2_NDSI_CLD


    def get_t_1_b_04(self):
        return self._T1_B04


    def get_t_2_b_04(self):
        return self._T2_B04


    def set_t_1_b_04(self, value):
        self._T1_B04 = value


    def set_t_2_b_04(self, value):
        self._T2_B04 = value


    def del_t_1_b_04(self):
        del self._T1_B04


    def del_t_2_b_04(self):
        del self._T2_B04


    def get_l_2_a_wvp_quantification_value(self):
        return self._L2A_WVP_QUANTIFICATION_VALUE


    def set_l_2_a_wvp_quantification_value(self, value):
        self._L2A_WVP_QUANTIFICATION_VALUE = value


    def del_l_2_a_wvp_quantification_value(self):
        del self._L2A_WVP_QUANTIFICATION_VALUE


    def get_product_version(self):
        return self._productVersion


    def set_product_version(self, value):
        self._productVersion = value


    def del_product_version(self):
        del self._productVersion


    def get_l1c_up_mtd_xml(self):
        return self._L1C_UP_MTD_XML


    def get_l1c_ds_mtd_xml(self):
        return self._L1C_DS_MTD_XML


    def get_l1c_up_id(self):
        return self._L1C_UP_ID


    def get_l1c_ds_id(self):
        return self._L1C_DS_ID


    def get_l1c_tile_id(self):
        return self._L1C_TILE_ID


    def set_l1c_up_mtd_xml(self, value):
        self._L1C_UP_MTD_XML = value


    def set_l1c_ds_mtd_xml(self, value):
        self._L1C_DS_MTD_XML = value


    def set_l1c_up_id(self, value):
        self._L1C_UP_ID = value


    def set_l1c_ds_id(self, value):
        self._L1C_DS_ID = value


    def set_l1c_tile_id(self, value):
        self._L1C_TILE_ID = value


    def del_l1c_up_mtd_xml(self):
        del self._L1C_UP_MTD_XML


    def del_l1c_ds_mtd_xml(self):
        del self._L1C_DS_MTD_XML


    def del_l1c_up_id(self):
        del self._L1C_UP_ID


    def del_l1c_ds_id(self):
        del self._L1C_DS_ID


    def del_l1c_tile_id(self):
        del self._L1C_TILE_ID



    def get_creation_date(self):
        return self._creationDate


    def set_creation_date(self, value):
        self._creationDate = value


    def del_creation_date(self):
        del self._creationDate


    def get_wvlsen(self):
        return self._wvlsen


    def get_fwhm(self):
        return self._fwhm


    def set_wvlsen(self, value):
        self._wvlsen = value


    def set_fwhm(self, value):
        self._fwhm = value


    def del_wvlsen(self):
        del self._wvlsen


    def del_fwhm(self):
        del self._fwhm


    def get_l2a_boa_quantification_value(self):
        return self._L2A_BOA_QUANTIFICATION_VALUE


    def get_l2a_wvp_quantification_value(self):
        return self._L2A_WVP_QUANTIFICATION_VALUE


    def get_l2a_aot_quantification_value(self):
        return self._L2A_AOT_QUANTIFICATION_VALUE


    def set_l2a_boa_quantification_value(self, value):
        self._L2A_BOA_QUANTIFICATION_VALUE = value


    def set_l2a_wvp_quantification_value(self, value):
        self._L2A_WVP_QUANTIFICATION_VALUE = value


    def set_l2a_aot_quantification_value(self, value):
        self._L2A_AOT_QUANTIFICATION_VALUE = value


    def del_l2a_boa_quantification_value(self):
        del self._L2A_BOA_QUANTIFICATION_VALUE


    def del_l2a_wvp_quantification_value(self):
        del self._L2A_WVP_QUANTIFICATION_VALUE


    def del_l2a_aot_quantification_value(self):
        del self._L2A_AOT_QUANTIFICATION_VALUE


    def get_l2a_up_dir(self):
        return self._L2A_UP_DIR


    def set_l2a_up_dir(self, value):
        self._L2A_UP_DIR = value


    def del_l2a_up_dir(self):
        del self._L2A_UP_DIR


    def get_l2a_up_id(self):
        return self._L2A_UP_ID


    def set_l2a_up_id(self, value):
        self._L2A_UP_ID = value


    def del_l2a_up_id(self):
        del self._L2A_UP_ID


    def get_l2a_tile_id(self):
        return self._L2A_TILE_ID


    def get_l1c_tile_mtd_xml(self):
        return self._L1C_TILE_MTD_XML


    def get_l2a_inspire_xml(self):
        return self._L2A_INSPIRE_XML


    def get_l2a_manifest_safe(self):
        return self._L2A_MANIFEST_SAFE


    def get_l2a_up_mtd_xml(self):
        return self._L2A_UP_MTD_XML


    def get_l2a_tile_mtd_xml(self):
        return self._L2A_TILE_MTD_XML


    def set_l1c_tile_mtd_xml(self, value):
        self._L1C_TILE_MTD_XML = value


    def set_l2a_tile_id(self, value):
        self._L2A_TILE_ID = value


    def set_l2a_inspire_xml(self, value):
        self._L2A_INSPIRE_XML = value


    def set_l2a_manifest_safe(self, value):
        self._L2A_MANIFEST_SAFE = value


    def set_l2a_up_mtd_xml(self, value):
        self._L2A_UP_MTD_XML = value


    def set_l2a_tile_mtd_xml(self, value):
        self._L2A_TILE_MTD_XML = value


    def del_l2a_ds_id(self):
        del self._L2A_DS_ID


    def del_l2a_tile_id(self):
        del self._L2A_TILE_ID


    def del_l1c_tile_mtd_xml(self):
        del self._L1C_TILE_MTD_XML


    def del_l2a_inspire_xml(self):
        del self._L2A_INSPIRE_XML


    def del_l2a_manifest_safe(self):
        del self._L2A_MANIFEST_SAFE


    def del_l2a_up_mtd_xml(self):
        del self._L2A_UP_MTD_XML


    def del_l2a_ds_mtd_xml(self):
        del self._L2A_DS_MTD_XML


    def del_l2a_tile_mtd_xml(self):
        del self._L2A_TILE_MTD_XML


    def get_entity_id(self):
        return self._entityId


    def get_acquisition_date(self):
        return self._acquisitionDate


    def get_orbit_path(self):
        return self._orbitPath


    def get_orbit_row(self):
        return self._orbitRow


    def get_target_path(self):
        return self._targetPath


    def get_target_row(self):
        return self._targetRow


    def get_station_sgs(self):
        return self._stationSgs


    def get_scene_start_time(self):
        return self._sceneStartTime


    def get_scene_stop_time(self):
        return self._sceneStopTime


    def set_entity_id(self, value):
        self._entityId = value


    def set_acquisition_date(self, value):
        self._acquisitionDate = value


    def set_orbit_path(self, value):
        self._orbitPath = value


    def set_orbit_row(self, value):
        self._orbitRow = value


    def set_target_path(self, value):
        self._targetPath = value


    def set_target_row(self, value):
        self._targetRow = value


    def set_station_sgs(self, value):
        self._stationSgs = value


    def set_scene_start_time(self, value):
        self._sceneStartTime = value


    def set_scene_stop_time(self, value):
        self._sceneStopTime = value


    def del_entity_id(self):
        del self._entityId


    def del_acquisition_date(self):
        del self._acquisitionDate


    def del_orbit_path(self):
        del self._orbitPath


    def del_orbit_row(self):
        del self._orbitRow


    def del_target_path(self):
        del self._targetPath


    def del_target_row(self):
        del self._targetRow


    def del_station_sgs(self):
        del self._stationSgs


    def del_scene_start_time(self):
        del self._sceneStartTime


    def del_scene_stop_time(self):
        del self._sceneStopTime


    def get_dn_scale(self):
        return self._dnScale


    def set_dn_scale(self, value):
        self._dnScale = value


    def del_dn_scale(self):
        del self._dnScale


    def get_d_2(self):
        return self._d2


    def get_c_0(self):
        return self._c0


    def get_c_1(self):
        return self._c1


    def get_e_0(self):
        return self._e0


    def set_d_2(self, value):
        self._d2 = value


    def set_c_0(self, value):
        self._c0 = value


    def set_c_1(self, value):
        self._c1 = value


    def set_e_0(self, value):
        self._e0 = value


    def del_d_2(self):
        del self._d2


    def del_c_0(self):
        del self._c0


    def del_c_1(self):
        del self._c1


    def del_e_0(self):
        del self._e0


    def get_processor_name(self):
        return self._processorName


    def get_processor_version(self):
        return self._processorVersion


    def get_processor_date(self):
        return self._processorDate


    def set_processor_name(self, value):
        self._processorName = value


    def set_processor_version(self, value):
        self._processorVersion = value


    def set_processor_date(self, value):
        self._processorDate = value


    def del_processor_name(self):
        del self._processorName


    def del_processor_version(self):
        del self._processorVersion


    def del_processor_date(self):
        del self._processorDate


    def get_ncols(self):
        return self._ncols


    def get_nrows(self):
        return self._nrows


    def get_nbnds(self):
        return self._nbnds


    def get_zenith_angle(self):
        return self._sza


    def get_azimuth_angle(self):
        return self._saa


    def get_gipp(self):
        return self._GIPP


    def get_ecmwf(self):
        return self._ECMWF


    def set_ncols(self, value):
        self._ncols = value


    def set_nrows(self, value):
        self._nrows = value


    def set_nbnds(self, value):
        self._nbnds = value


    def set_zenith_angle(self, value):
        self._sza = value


    def set_azimuth_angle(self, value):
        self._saa = value


    def set_gipp(self, value):
        self._GIPP = value


    def set_ecmwf(self, value):
        self._ECMWF = value


    def del_ncols(self):
        del self._ncols


    def del_nrows(self):
        del self._nrows


    def del_nbnds(self):
        del self._nbnds


    def del_zenith_angle(self):
        del self._sza


    def del_azimuth_angle(self):
        del self._saa


    def del_gipp(self):
        del self._GIPP


    def del_ecmwf(self):
        del self._ECMWF


    def get_home(self):
        return self._home


    def get_source_dir(self):
        return self._sourceDir


    def get_data_dir(self):
        return self._dataDir


    def get_config_dir(self):
        return self._configDir


    def get_bin_dir(self):
        return self._binDir


    def get_lib_dir(self):
        return self._libDir


    def get_log_dir(self):
        return self._logDir


    def get_config_fn(self):
        return self._configFn


    def get_input_fn(self):
        return self._inputFn


    def get_aot_fn(self):
        return self._aotFn


    def get_aspect_fn(self):
        return self._aspectFn


    def get_atm_data_fn(self):
        return self._atmDataFn


    def get_class_map_fn(self):
        return self._classMapFn


    def get_cloud_qi_map_fn(self):
        return self._cloudQiMapFn


    def get_ddv_fn(self):
        return self._ddvFn


    def get_hcw_fn(self):
        return self._hcwFn


    def get_ilumination_fn(self):
        return self._iluminationFn


    def get_sky_view_fn(self):
        return self._skyViewFn


    def get_slope_fn(self):
        return self._slopeFn


    def get_snow_qi_map_fn(self):
        return self._snowQiMapFn


    def get_vis_index_fn(self):
        return self._visIndexFn


    def get_water_vapor_fn(self):
        return self._waterVaporFn


    def get_adj_km(self):
        return self._adj_km


    def get_beta_thr(self):
        return self._beta_thr


    def get_ch_940(self):
        return self._ch940


    def get_cellsize(self):
        return self._cellsize


    def get_cloud_refl_thr_blu(self):
        return self._cloud_refl_thr_blu


    def get_dem_unit(self):
        return self._dem_unit


    def get_ibrdf(self):
        return self._ibrdf


    def get_ibrdf_dark(self):
        return self._ibrdf_dark


    def get_icl_shadow(self):
        return self._icl_shadow


    def get_iclshad_mask(self):
        return self._iclshad_mask


    def get_ihaze(self):
        return self._ihaze


    def get_ihcw(self):
        return self._ihcw


    def get_ihot_dynr(self):
        return self._ihot_dynr


    def get_ihot_mask(self):
        return self._ihot_mask


    def get_intpol_1400(self):
        return self._intpol1400


    def get_intpol_725_825(self):
        return self._intpol725_825


    def get_intpol_760(self):
        return self._intpol760


    def get_intpol_940_1130(self):
        return self._intpol940_1130


    def get_istretch_type(self):
        return self._istretch_type


    def get_iwat_shd(self):
        return self._iwat_shd


    def get_iwaterwv(self):
        return self._iwaterwv


    def get_iwv_watermask(self):
        return self._iwv_watermask


    def get_ksolflux(self):
        return self._ksolflux


    def get_altit(self):
        return self._altit


    def get_npref(self):
        return self._npref


    def get_phi_scl_min(self):
        return self._phi_scl_min


    def get_phi_unscl_max(self):
        return self._phi_unscl_max


    def get_pixelsize(self):
        return self._pixelsize


    def get_resolution(self):
        return self._resolution


    def get_ratio_blu_red(self):
        return self._ratio_blu_red


    def get_ratio_red_swir(self):
        return self._ratio_red_swir


    def get_refl_cutoff(self):
        return self._refl_cutoff


    def get_rel_saturation(self):
        return self._rel_saturation


    def get_thr_shad(self):
        return self._thr_shad


    def get_thv(self):
        return self._thv


    def get_phiv(self):
        return self._phiv


    def get_smooth_wvmap(self):
        return self._smooth_wvmap


    def get_solaz(self):
        return self._solaz


    def get_solaz_arr(self):
        return self._solaz_arr


    def get_solze(self):
        return self._solze


    def get_solze_arr(self):
        return self._solze_arr


    def get_thr_g(self):
        return self._thr_g


    def get_vaa_arr(self):
        return self._vaa_arr


    def get_visibility(self):
        return self._visibility


    def get_vza_arr(self):
        return self._vza_arr


    def get_water_refl_thr_nir(self):
        return self._water_refl_thr_nir


    def get_water_refl_thr_swir_1(self):
        return self._water_refl_thr_swir1


    def get_wl_1130a(self):
        return self._wl1130a


    def get_wl_1400a(self):
        return self._wl1400a


    def get_wl_1900a(self):
        return self._wl1900a


    def get_wl_940a(self):
        return self._wl940a


    def get_wv_thr_cirrus(self):
        return self._wv_thr_cirrus


    def set_home(self, value):
        self._home = value


    def set_source_dir(self, value):
        self._sourceDir = value


    def set_data_dir(self, value):
        self._dataDir = value

    def set_config_dir(self, value):
        self._configDir = value


    def set_bin_dir(self, value):
        self._binDir = value


    def set_lib_dir(self, value):
        self._libDir = value


    def set_log_dir(self, value):
        self._logDir = value


    def set_config_fn(self, value):
        self._configFn = value 

    def set_input_fn(self, value):
        self._inputFn = os.path.join(self._dataDir, value)


    def set_aot_fn(self, value):
        self._aotFn =  os.path.join(self._dataDir, value)


    def set_aspect_fn(self, value):
        self._aspectFn =  os.path.join(self._dataDir, value)


    def set_atm_data_fn(self, value):
        self._atmDataFn = value


    def set_class_map_fn(self, value):
        self._classMapFn =  os.path.join(self._dataDir, value)


    def set_cloud_qi_map_fn(self, value):
        self._cloudQiMapFn =  os.path.join(self._dataDir, value)


    def set_ddv_fn(self, value):
        self._ddvFn =  os.path.join(self._dataDir, value)


    def set_hcw_fn(self, value):
        self._hcwFn =  os.path.join(self._dataDir, value)


    def set_ilumination_fn(self, value):
        self._iluminationFn =  os.path.join(self._dataDir, value)


    def set_sky_view_fn(self, value):
        self._skyViewFn =  os.path.join(self._dataDir, value)


    def set_slope_fn(self, value):
        self._slopeFn =  os.path.join(self._dataDir, value)


    def set_snow_qi_map_fn(self, value):
        self._snowQiMapFn =  os.path.join(self._dataDir, value)


    def set_vis_index_fn(self, value):
        self._visIndexFn =  os.path.join(self._dataDir, value)


    def set_water_vapor_fn(self, value):
        self._waterVaporFn =  os.path.join(self._dataDir, value)


    def set_adj_km(self, value):
        self._adj_km = value


    def set_beta_thr(self, value):
        self._beta_thr = value


    def set_ch_940(self, value):
        self._ch940 = value


    def set_cellsize(self, value):
        self._cellsize = value


    def set_cloud_refl_thr_blu(self, value):
        self._cloud_refl_thr_blu = value


    def set_date(self, value):
        self._date = value


    def set_dem_unit(self, value):
        self._dem_unit = value


    def set_ibrdf(self, value):
        self._ibrdf = value


    def set_ibrdf_dark(self, value):
        self._ibrdf_dark = value


    def set_icl_shadow(self, value):
        self._icl_shadow = value


    def set_iclshad_mask(self, value):
        self._iclshad_mask = value


    def set_ihaze(self, value):
        self._ihaze = value


    def set_ihcw(self, value):
        self._ihcw = value


    def set_ihot_dynr(self, value):
        self._ihot_dynr = value


    def set_ihot_mask(self, value):
        self._ihot_mask = value


    def set_intpol_1400(self, value):
        self._intpol1400 = value


    def set_intpol_725_825(self, value):
        self._intpol725_825 = value


    def set_intpol_760(self, value):
        self._intpol760 = value


    def set_intpol_940_1130(self, value):
        self._intpol940_1130 = value


    def set_istretch_type(self, value):
        self._istretch_type = value


    def set_iwat_shd(self, value):
        self._iwat_shd = value


    def set_iwaterwv(self, value):
        self._iwaterwv = value


    def set_iwv_watermask(self, value):
        self._iwv_watermask = value


    def set_ksolflux(self, value):
        self._ksolflux = value


    def set_altit(self, value):
        self._altit = value


    def set_npref(self, value):
        self._npref = value


    def set_phi_scl_min(self, value):
        self._phi_scl_min = value


    def set_phi_unscl_max(self, value):
        self._phi_unscl_max = value


    def set_pixelsize(self, value):
        self._pixelsize = value


    def set_resolution(self, value):
        self._resolution = value
        self._pixelsize = value
        self._cellsize = value


    def set_ratio_blu_red(self, value):
        self._ratio_blu_red = value


    def set_ratio_red_swir(self, value):
        self._ratio_red_swir = value


    def set_refl_cutoff(self, value):
        self._refl_cutoff = value


    def set_rel_saturation(self, value):
        self._rel_saturation = value


    def set_thr_shad(self, value):
        self._thr_shad = value


    def set_thv(self, value):
        self._thv = value


    def set_phiv(self, value):
        self._phiv = value


    def set_smooth_wvmap(self, value):
        self._smooth_wvmap = value


    def set_solaz(self, value):
        self._solaz = value


    def set_solaz_arr(self, value):
        self._solaz_arr = value


    def set_solze(self, value):
        self._solze = value


    def set_solze_arr(self, value):
        self._solze_arr = value


    def set_thr_clear_water(self, value):
        self._thr_clear_water = value


    def set_thr_haze_water(self, value):
        self._thr_haze_water = value


    def set_thr_g(self, value):
        self._thr_g = value


    def set_vaa_arr(self, value):
        self._vaa_arr = value


    def set_visibility(self, value):
        self._visibility = value


    def set_vza_arr(self, value):
        self._vza_arr = value


    def set_water_refl_thr_nir(self, value):
        self._water_refl_thr_nir = value


    def set_water_refl_thr_swir_1(self, value):
        self._water_refl_thr_swir1 = value


    def set_wl_1130a(self, value):
        self._wl1130a = value


    def set_wl_1400a(self, value):
        self._wl1400a = value


    def set_wl_1900a(self, value):
        self._wl1900a = value


    def set_wl_940a(self, value):
        self._wl940a = value


    def set_wv_thr_cirrus(self, value):
        self._wv_thr_cirrus = value


    def del_home(self):
        del self._home


    def del_source_dir(self):
        del self._sourceDir


    def del_data_dir(self):
        del self._dataDir


    def del_config_dir(self):
        del self._configDir


    def del_bin_dir(self):
        del self._binDir


    def del_lib_dir(self):
        del self._libDir


    def del_log_dir(self):
        del self._logDir


    def del_config_fn(self):
        del self._configFn


    def del_input_fn(self):
        del self._inputFn


    def del_aot_fn(self):
        del self._aotFn


    def del_aspect_fn(self):
        del self._aspectFn


    def del_atm_data_fn(self):
        del self._atmDataFn


    def del_class_map_fn(self):
        del self._classMapFn


    def del_cloud_qi_map_fn(self):
        del self._cloudQiMapFn


    def del_ddv_fn(self):
        del self._ddvFn


    def del_hcw_fn(self):
        del self._hcwFn


    def del_ilumination_fn(self):
        del self._iluminationFn


    def del_sky_view_fn(self):
        del self._skyViewFn


    def del_slope_fn(self):
        del self._slopeFn


    def del_snow_qi_map_fn(self):
        del self._snowQiMapFn


    def del_vis_index_fn(self):
        del self._visIndexFn


    def del_water_vapor_fn(self):
        del self._waterVaporFn


    def del_adj_km(self):
        del self._adj_km


    def del_beta_thr(self):
        del self._beta_thr


    def del_ch_940(self):
        del self._ch940


    def del_cellsize(self):
        del self._cellsize


    def del_cloud_refl_thr_blu(self):
        del self._cloud_refl_thr_blu


    def del_date(self):
        del self._date


    def del_dem_unit(self):
        del self._dem_unit


    def del_ibrdf(self):
        del self._ibrdf


    def del_ibrdf_dark(self):
        del self._ibrdf_dark


    def del_icl_shadow(self):
        del self._icl_shadow


    def del_iclshad_mask(self):
        del self._iclshad_mask


    def del_ihaze(self):
        del self._ihaze


    def del_ihcw(self):
        del self._ihcw


    def del_ihot_dynr(self):
        del self._ihot_dynr


    def del_ihot_mask(self):
        del self._ihot_mask


    def del_intpol_1400(self):
        del self._intpol1400


    def del_intpol_725_825(self):
        del self._intpol725_825


    def del_intpol_760(self):
        del self._intpol760


    def del_intpol_940_1130(self):
        del self._intpol940_1130


    def del_istretch_type(self):
        del self._istretch_type


    def del_iwat_shd(self):
        del self._iwat_shd


    def del_iwaterwv(self):
        del self._iwaterwv


    def del_iwv_watermask(self):
        del self._iwv_watermask


    def del_ksolflux(self):
        del self._ksolflux


    def del_altit(self):
        del self._altit


    def del_npref(self):
        del self._npref


    def del_phi_scl_min(self):
        del self._phi_scl_min


    def del_phi_unscl_max(self):
        del self._phi_unscl_max


    def del_pixelsize(self):
        del self._pixelsize


    def del_resolution(self):
        del self._resolution


    def del_ratio_blu_red(self):
        del self._ratio_blu_red


    def del_ratio_red_swir(self):
        del self._ratio_red_swir


    def del_refl_cutoff(self):
        del self._refl_cutoff


    def del_rel_saturation(self):
        del self._rel_saturation


    def del_thr_shad(self):
        del self._thr_shad


    def del_thv(self):
        del self._thv


    def del_phiv(self):
        del self._phiv


    def del_smooth_wvmap(self):
        del self._smooth_wvmap


    def del_solaz(self):
        del self._solaz


    def del_solaz_arr(self):
        del self._solaz_arr


    def del_solze(self):
        del self._solze


    def del_solze_arr(self):
        del self._solze_arr


    def del_thr_clear_water(self):
        del self._thr_clear_water


    def del_thr_haze_water(self):
        del self._thr_haze_water


    def del_thr_g(self):
        del self._thr_g


    def del_vaa_arr(self):
        del self._vaa_arr


    def del_visibility(self):
        del self._visibility


    def del_vza_arr(self):
        del self._vza_arr


    def del_water_refl_thr_nir(self):
        del self._water_refl_thr_nir


    def del_water_refl_thr_swir_1(self):
        del self._water_refl_thr_swir1


    def del_wl_1130a(self):
        del self._wl1130a


    def del_wl_1400a(self):
        del self._wl1400a


    def del_wl_1900a(self):
        del self._wl1900a


    def del_wl_940a(self):
        del self._wl940a


    def del_wv_thr_cirrus(self):
        del self._wv_thr_cirrus


    def set_l2a_ds_id(self, value):
        self._L2A_DS_ID = value
        if value in self._L2A_DS_LST:
            return
        self._L2A_DS_LST.append(value)


    def get_l2a_ds_id(self):
        try:
            tileId = self._L2A_TILE_ID[25:40]
            for ref in self._L2A_DS_LST:
                if tileId in ref:
                    return ref
        except:
            return self._L2A_DS_ID


    def set_l2a_ds_mtd_xml(self, value):
        self._L2A_DS_MTD_XML = value
        if value in self._L2A_DS_MTD_LST:
            return
        self._L2A_DS_MTD_LST.append(value)


    def get_l2a_ds_mtd_xml(self):
        try:
            tileId = self._L2A_TILE_ID[25:40]        
            for ref in self._L2A_DS_MTD_LST:
                if tileId in ref:
                    return ref
        except:
            return self._L2A_DS_MTD_XML


    # Properties:
    processorName = property(get_processor_name, set_processor_name, del_processor_name, "processorName's docstring")
    processorVersion = property(get_processor_version, set_processor_version, del_processor_version, "processorVersion's docstring")
    processorDate = property(get_processor_date, set_processor_date, del_processor_date, "processorDate's docstring")
    home = property(get_home, set_home, del_home, "home's docstring")
    sourceDir = property(get_source_dir, set_source_dir, del_source_dir, "sourceDir's docstring")
    dataDir = property(get_data_dir, set_data_dir, del_data_dir, "dataDir's docstring")
    configDir = property(get_config_dir, set_config_dir, del_config_dir, "configDir's docstring")
    binDir = property(get_bin_dir, set_bin_dir, del_bin_dir, "binDir's docstring")
    libDir = property(get_lib_dir, set_lib_dir, del_lib_dir, "libDir's docstring")
    logDir = property(get_log_dir, set_log_dir, del_log_dir, "logDir's docstring")
    configFn = property(get_config_fn, set_config_fn, del_config_fn, "configFn's docstring")
    configSC = property(get_config_sc, set_config_sc, del_config_sc, "configSC's docstring")
    configAC = property(get_config_ac, set_config_ac, del_config_ac, "configAC's docstring")
    atmDataFn = property(get_atm_data_fn, set_atm_data_fn, del_atm_data_fn, "atmDataFn's docstring")
    adj_km = property(get_adj_km, set_adj_km, del_adj_km, "adj_km's docstring")
    beta_thr = property(get_beta_thr, set_beta_thr, del_beta_thr, "beta_thr's docstring")
    ch940 = property(get_ch_940, set_ch_940, del_ch_940, "ch940's docstring")
    cellsize = property(get_cellsize, set_cellsize, del_cellsize, "cellsize's docstring")
    cloud_refl_thr_blu = property(get_cloud_refl_thr_blu, set_cloud_refl_thr_blu, del_cloud_refl_thr_blu, "cloud_refl_thr_blu's docstring")
    dem_unit = property(get_dem_unit, set_dem_unit, del_dem_unit, "dem_unit's docstring")
    ibrdf = property(get_ibrdf, set_ibrdf, del_ibrdf, "ibrdf's docstring")
    thr_g = property(get_thr_g, set_thr_g, del_thr_g, "thr_g's docstring")
    ibrdf_dark = property(get_ibrdf_dark, set_ibrdf_dark, del_ibrdf_dark, "ibrdf_dark's docstring")
    icl_shadow = property(get_icl_shadow, set_icl_shadow, del_icl_shadow, "icl_shadow's docstring")
    iclshad_mask = property(get_iclshad_mask, set_iclshad_mask, del_iclshad_mask, "iclshad_mask's docstring")
    ihaze = property(get_ihaze, set_ihaze, del_ihaze, "ihaze's docstring")
    ihcw = property(get_ihcw, set_ihcw, del_ihcw, "ihcw's docstring")
    ihot_dynr = property(get_ihot_dynr, set_ihot_dynr, del_ihot_dynr, "ihot_dynr's docstring")
    ihot_mask = property(get_ihot_mask, set_ihot_mask, del_ihot_mask, "ihot_mask's docstring")
    intpol1400 = property(get_intpol_1400, set_intpol_1400, del_intpol_1400, "intpol1400's docstring")
    intpol725_825 = property(get_intpol_725_825, set_intpol_725_825, del_intpol_725_825, "intpol725_825's docstring")
    intpol760 = property(get_intpol_760, set_intpol_760, del_intpol_760, "intpol760's docstring")
    intpol940_1130 = property(get_intpol_940_1130, set_intpol_940_1130, del_intpol_940_1130, "intpol940_1130's docstring")
    istretch_type = property(get_istretch_type, set_istretch_type, del_istretch_type, "istretch_type's docstring")
    iwat_shd = property(get_iwat_shd, set_iwat_shd, del_iwat_shd, "iwat_shd's docstring")
    iwaterwv = property(get_iwaterwv, set_iwaterwv, del_iwaterwv, "iwaterwv's docstring")
    iwv_watermask = property(get_iwv_watermask, set_iwv_watermask, del_iwv_watermask, "iwv_watermask's docstring")
    ksolflux = property(get_ksolflux, set_ksolflux, del_ksolflux, "ksolflux's docstring")
    altit = property(get_altit, set_altit, del_altit, "altit's docstring")
    npref = property(get_npref, set_npref, del_npref, "npref's docstring")
    phi_scl_min = property(get_phi_scl_min, set_phi_scl_min, del_phi_scl_min, "phi_scl_min's docstring")
    phi_unscl_max = property(get_phi_unscl_max, set_phi_unscl_max, del_phi_unscl_max, "phi_unscl_max's docstring")
    pixelsize = property(get_pixelsize, set_pixelsize, del_pixelsize, "pixelsize's docstring")
    resolution = property(get_resolution, set_resolution, del_resolution, "resolution's docstring")
    ratio_blu_red = property(get_ratio_blu_red, set_ratio_blu_red, del_ratio_blu_red, "ratio_blu_red's docstring")
    ratio_red_swir = property(get_ratio_red_swir, set_ratio_red_swir, del_ratio_red_swir, "ratio_red_swir's docstring")
    refl_cutoff = property(get_refl_cutoff, set_refl_cutoff, del_refl_cutoff, "refl_cutoff's docstring")
    rel_saturation = property(get_rel_saturation, set_rel_saturation, del_rel_saturation, "rel_saturation's docstring")
    thr_shad = property(get_thr_shad, set_thr_shad, del_thr_shad, "thr_shad's docstring")
    thv = property(get_thv, set_thv, del_thv, "thv's docstring")
    phiv = property(get_phiv, set_phiv, del_phiv, "phiv's docstring")
    smooth_wvmap = property(get_smooth_wvmap, set_smooth_wvmap, del_smooth_wvmap, "smooth_wvmap's docstring")
    solaz = property(get_solaz, set_solaz, del_solaz, "solaz's docstring")
    solaz_arr = property(get_solaz_arr, set_solaz_arr, del_solaz_arr, "solaz_arr's docstring")
    solze = property(get_solze, set_solze, del_solze, "solze's docstring")
    solze_arr = property(get_solze_arr, set_solze_arr, del_solze_arr, "solze_arr's docstring")
    vaa_arr = property(get_vaa_arr, set_vaa_arr, del_vaa_arr, "vaa_arr's docstring")
    visibility = property(get_visibility, set_visibility, del_visibility, "visibility's docstring")
    vza_arr = property(get_vza_arr, set_vza_arr, del_vza_arr, "vza_arr's docstring")
    water_refl_thr_nir = property(get_water_refl_thr_nir, set_water_refl_thr_nir, del_water_refl_thr_nir, "water_refl_thr_nir's docstring")
    water_refl_thr_swir1 = property(get_water_refl_thr_swir_1, set_water_refl_thr_swir_1, del_water_refl_thr_swir_1, "water_refl_thr_swir1's docstring")
    wl1130a = property(get_wl_1130a, set_wl_1130a, del_wl_1130a, "wl1130a's docstring")
    wl1400a = property(get_wl_1400a, set_wl_1400a, del_wl_1400a, "wl1400a's docstring")
    wl1900a = property(get_wl_1900a, set_wl_1900a, del_wl_1900a, "wl1900a's docstring")
    wl940a = property(get_wl_940a, set_wl_940a, del_wl_940a, "wl940a's docstring")
    wv_thr_cirrus = property(get_wv_thr_cirrus, set_wv_thr_cirrus, del_wv_thr_cirrus, "wv_thr_cirrus's docstring")
    ncols = property(get_ncols, set_ncols, del_ncols, "ncols's docstring")
    nrows = property(get_nrows, set_nrows, del_nrows, "nrows's docstring")
    nbnds = property(get_nbnds, set_nbnds, del_nbnds, "nbnds's docstring")
    zenith_angle = property(get_zenith_angle, set_zenith_angle, del_zenith_angle, "zenith_angle's docstring")
    azimuth_angle = property(get_azimuth_angle, set_azimuth_angle, del_azimuth_angle, "azimuth_angle's docstring")
    GIPP = property(get_gipp, set_gipp, del_gipp, "GIPP's docstring")
    ECMWF = property(get_ecmwf, set_ecmwf, del_ecmwf, "ECMWF's docstring")
    d2 = property(get_d_2, set_d_2, del_d_2, "d2's docstring")
    c0 = property(get_c_0, set_c_0, del_c_0, "c0's docstring")
    c1 = property(get_c_1, set_c_1, del_c_1, "c1's docstring")
    e0 = property(get_e_0, set_e_0, del_e_0, "e0's docstring")
    wvlsen = property(get_wvlsen, set_wvlsen, del_wvlsen, "wvlsen's docstring")
    fwhm = property(get_fwhm, set_fwhm, del_fwhm, "fwhm's docstring")
    dnScale = property(get_dn_scale, set_dn_scale, del_dn_scale, "dnScale's docstring")
    entityId = property(get_entity_id, set_entity_id, del_entity_id, "entityId's docstring")
    acquisitionDate = property(get_acquisition_date, set_acquisition_date, del_acquisition_date, "acquisitionDate's docstring")
    orbitPath = property(get_orbit_path, set_orbit_path, del_orbit_path, "orbitPath's docstring")
    orbitRow = property(get_orbit_row, set_orbit_row, del_orbit_row, "orbitRow's docstring")
    targetPath = property(get_target_path, set_target_path, del_target_path, "targetPath's docstring")
    targetRow = property(get_target_row, set_target_row, del_target_row, "targetRow's docstring")
    stationSgs = property(get_station_sgs, set_station_sgs, del_station_sgs, "stationSgs's docstring")
    sceneStartTime = property(get_scene_start_time, set_scene_start_time, del_scene_start_time, "sceneStartTime's docstring")
    sceneStopTime = property(get_scene_stop_time, set_scene_stop_time, del_scene_stop_time, "sceneStopTime's docstring")
    L2A_INSPIRE_XML = property(get_l2a_inspire_xml, set_l2a_inspire_xml, del_l2a_inspire_xml, "L2A_INSPIRE_XML's docstring")
    L2A_MANIFEST_SAFE = property(get_l2a_manifest_safe, set_l2a_manifest_safe, del_l2a_manifest_safe, "L2A_MANIFEST_SAFE's docstring")
    L1C_UP_MTD_XML = property(get_l1c_up_mtd_xml, set_l1c_up_mtd_xml, del_l1c_up_mtd_xml, "L1C_UP_MTD_XML's docstring")
    L1C_DS_MTD_XML = property(get_l1c_ds_mtd_xml, set_l1c_ds_mtd_xml, del_l1c_ds_mtd_xml, "L1C_DS_MTD_XML's docstring")
    L1C_TILE_MTD_XML = property(get_l1c_tile_mtd_xml, set_l1c_tile_mtd_xml, del_l1c_tile_mtd_xml, "L1C_TILE_MTD_XML's docstring")
    L1C_UP_ID = property(get_l1c_up_id, set_l1c_up_id, del_l1c_up_id, "L1C_UP_ID's docstring")
    L1C_DS_ID = property(get_l1c_ds_id, set_l1c_ds_id, del_l1c_ds_id, "L1C_DS_ID's docstring")
    L1C_TILE_ID = property(get_l1c_tile_id, set_l1c_tile_id, del_l1c_tile_id, "L1C_TILE_ID's docstring")
    L2A_UP_MTD_XML = property(get_l2a_up_mtd_xml, set_l2a_up_mtd_xml, del_l2a_up_mtd_xml, "L2A_USER_PRODUCT_MTD_XML's docstring")
    L2A_DS_MTD_XML = property(get_l2a_ds_mtd_xml, set_l2a_ds_mtd_xml, del_l2a_ds_mtd_xml, "L2A_DS_MTD_XML's docstring")
    L2A_TILE_MTD_XML = property(get_l2a_tile_mtd_xml, set_l2a_tile_mtd_xml, del_l2a_tile_mtd_xml, "L2A_TILE_MTD_XML's docstring")
    L2A_TILE_ID = property(get_l2a_tile_id, set_l2a_tile_id, del_l2a_tile_id, "L2A_TILE_ID's docstring")
    L2A_DS_ID = property(get_l2a_ds_id, set_l2a_ds_id, del_l2a_ds_id, "L2A_TILE_ID's docstring")
    L2A_UP_ID = property(get_l2a_up_id, set_l2a_up_id, del_l2a_up_id, "L2A_UP_ID's docstring")
    L2A_UP_DIR = property(get_l2a_up_dir, set_l2a_up_dir, del_l2a_up_dir, "L2A_UP_DIR's docstring")
    L2A_BOA_QUANTIFICATION_VALUE = property(get_l2a_boa_quantification_value, set_l2a_boa_quantification_value, del_l2a_boa_quantification_value, "L2A_BOA_QUANTIFICATION_VALUE's docstring")
    L2A_WVP_QUANTIFICATION_VALUE = property(get_l2a_wvp_quantification_value, set_l2a_wvp_quantification_value, del_l2a_wvp_quantification_value, "L2A_WVP_QUANTIFICATION_VALUE's docstring")
    L2A_AOT_QUANTIFICATION_VALUE = property(get_l2a_aot_quantification_value, set_l2a_aot_quantification_value, del_l2a_aot_quantification_value, "L2A_AOT_QUANTIFICATION_VALUE's docstring")
    creationDate = property(get_creation_date, set_creation_date, del_creation_date, "creationDate's docstring")
    productVersion = property(get_product_version, set_product_version, del_product_version, "productVersion's docstring")
#     tTotal = property(get_t_total, set_t_total, del_t_total, "tTotal's docstring")
    T1_B04 = property(get_t_1_b_04, set_t_1_b_04, del_t_1_b_04, "T1_B04's docstring")
    T2_B04 = property(get_t_2_b_04, set_t_2_b_04, del_t_2_b_04, "T2_B04's docstring")
    T1_NDSI_CLD = property(get_t_1_ndsi_cld, set_t_1_ndsi_cld, del_t_1_ndsi_cld, "T1_NDSI_CLD's docstring")
    T2_NDSI_CLD = property(get_t_2_ndsi_cld, set_t_2_ndsi_cld, del_t_2_ndsi_cld, "T2_NDSI_CLD's docstring")
    T1_NDSI_SNW = property(get_t_1_ndsi_snw, set_t_1_ndsi_snw, del_t_1_ndsi_snw, "T1_NDSI_SNW's docstring")
    T2_NDSI_SNW = property(get_t_2_ndsi_snw, set_t_2_ndsi_snw, del_t_2_ndsi_snw, "T2_NDSI_SNW's docstring")
    T1_B02 = property(get_t_1_b_02, set_t_1_b_02, del_t_1_b_02, "T1_B02's docstring")
    T2_B02 = property(get_t_2_b_02, set_t_2_b_02, del_t_2_b_02, "T2_B02's docstring")
    T1_B8A = property(get_t_1_b_8_a, set_t_1_b_8_a, del_t_1_b_8_a, "T1_B8A's docstring")
    T2_B8A = property(get_t_2_b_8_a, set_t_2_b_8_a, del_t_2_b_8_a, "T2_B8A's docstring")
    T1_B10 = property(get_t_1_b_10, set_t_1_b_10, del_t_1_b_10, "T1_B10's docstring")
    T2_B10 = property(get_t_2_b_10, set_t_2_b_10, del_t_2_b_10, "T2_B10's docstring")
    T1_B12 = property(get_t_1_b_12, set_t_1_b_12, del_t_1_b_12, "T1_B12's docstring")
    T2_B12 = property(get_t_2_b_12, set_t_2_b_12, del_t_2_b_12, "T2_B12's docstring")
    T1_R_B02_B04 = property(get_t_1_r_b_02_b_04, set_t_1_r_b_02_b_04, del_t_1_r_b_02_b_04, "T1_R_B02_B04's docstring")
    T2_R_B02_B04 = property(get_t_2_r_b_02_b_04, set_t_2_r_b_02_b_04, del_t_2_r_b_02_b_04, "T2_R_B02_B04's docstring")
    T1_R_B8A_B03 = property(get_t_1_r_b_8_a_b_03, set_t_1_r_b_8_a_b_03, del_t_1_r_b_8_a_b_03, "T1_R_B8A_B03's docstring")
    T2_R_B8A_B03 = property(get_t_2_r_b_8_a_b_03, set_t_2_r_b_8_a_b_03, del_t_2_r_b_8_a_b_03, "T2_R_B8A_B03's docstring")
    T1_R_B8A_B11 = property(get_t_1_r_b_8_a_b_11, set_t_1_r_b_8_a_b_11, del_t_1_r_b_8_a_b_11, "T1_R_B8A_B11's docstring")
    T2_R_B8A_B11 = property(get_t_2_r_b_8_a_b_11, set_t_2_r_b_8_a_b_11, del_t_2_r_b_8_a_b_11, "T2_R_B8A_B11's docstring")
    T1_SNOW = property(get_t_1_snow, set_t_1_snow, del_t_1_snow, "T1_SNOW's docstring")
    T2_SNOW = property(get_t_2_snow, set_t_2_snow, del_t_2_snow, "T2_SNOW's docstring")
    T1_NDVI = property(get_t_1_ndvi, set_t_1_ndvi, del_t_1_ndvi, "T1_NDVI's docstring")
    T2_NDVI = property(get_t_2_ndvi, set_t_2_ndvi, del_t_2_ndvi, "T2_NDVI's docstring")
    T11_B02 = property(get_t_11_b_02, set_t_11_b_02, del_t_11_b_02, "T11_B02's docstring")
    T12_B02 = property(get_t_12_b_02, set_t_12_b_02, del_t_12_b_02, "T12_B02's docstring")
    T11_R_B02_B11 = property(get_t_11_r_b_02_b_11, set_t_11_r_b_02_b_11, del_t_11_r_b_02_b_11, "T11_R_B02_B11's docstring")
    T12_R_B02_B11 = property(get_t_12_r_b_02_b_11, set_t_12_r_b_02_b_11, del_t_12_r_b_02_b_11, "T12_R_B02_B11's docstring")
    T21_B12 = property(get_t_21_b_12, set_t_21_b_12, del_t_21_b_12, "T21_B12's docstring")
    T22_B12 = property(get_t_22_b_12, set_t_22_b_12, del_t_22_b_12, "T22_B12's docstring")
    T21_R_B02_B11 = property(get_t_21_r_b_02_b_11, set_t_21_r_b_02_b_11, del_t_21_r_b_02_b_11, "T21_R_B02_B11's docstring")
    T22_R_B02_B11 = property(get_t_22_r_b_02_b_11, set_t_22_r_b_02_b_11, del_t_22_r_b_02_b_11, "T22_R_B02_B11's docstring")
    noData = property(get_no_data, set_no_data, del_no_data, "noData's docstring")
    saturatedDefective = property(get_saturated_defective, set_saturated_defective, del_saturated_defective, "saturatedDefective's docstring")
    darkFeatures = property(get_dark_features, set_dark_features, del_dark_features, "darkFeatures's docstring")
    cloudShadows = property(get_cloud_shadows, set_cloud_shadows, del_cloud_shadows, "cloudShadows's docstring")
    vegetation = property(get_vegetation, set_vegetation, del_vegetation, "vegetation's docstring")
    bareSoils = property(get_bare_soils, set_bare_soils, del_bare_soils, "bareSoils's docstring")
    water = property(get_water, set_water, del_water, "water's docstring")
    lowProbaClouds = property(get_low_proba_clouds, set_low_proba_clouds, del_low_proba_clouds, "lowProbaClouds's docstring")
    medProbaClouds = property(get_med_proba_clouds, set_med_proba_clouds, del_med_proba_clouds, "medProbaClouds's docstring")
    highProbaClouds = property(get_high_proba_clouds, set_high_proba_clouds, del_high_proba_clouds, "highProbaClouds's docstring")
    thinCirrus = property(get_thin_cirrus, set_thin_cirrus, del_thin_cirrus, "thinCirrus's docstring")
    snowIce = property(get_snow_ice, set_snow_ice, del_snow_ice, "snowIce's docstring")
    T_SDW = property(get_t_sdw, set_t_sdw, del_t_sdw, "T_SDW's docstring")
    T_B02_B12 = property(get_t_b_02_b_12, set_t_b_02_b_12, del_t_b_02_b_12, "T_B02_B12's docstring")
    medianFilter = property(get_median_filter, set_median_filter, del_median_filter, "medianFilter's docstring")
    demDirectory = property(get_dem_directory, set_dem_directory, del_dem_directory, "demDirectory's docstring")
    demReference = property(get_dem_reference, set_dem_reference, del_dem_reference, "demReference's docstring")
    snowMapReference = property(get_snow_map_reference, set_snow_map_reference, del_snow_map_reference, "snowMapReference's docstring")
    scOnly = property(get_sc_only, set_sc_only, del_sc_only, "scOnly's docstring")
    workDir = property(get_work_dir, set_work_dir, del_work_dir, "workDir's docstring")
    processed60 = property(get_processed_60, set_processed_60, del_processed_60, "processed60's docstring")
    processed20 = property(get_processed_20, set_processed_20, del_processed_20, "processed20's docstring")
    processed10 = property(get_processed_10, set_processed_10, del_processed_10, "processed10's docstring")
    refresh = property(get_refresh, set_refresh, del_refresh, "refresh's docstring")
    logger = property(get_logger, set_logger, del_logger, "log's docstring")
    targetDirectory = property(get_target_directory, set_target_directory, del_target_directory, "targetDirectory's docstring")
    logLevel = property(get_log_level, set_log_level, del_log_level, "logLevel's docstring")
    selectedTile = property(get_selected_tile, set_selected_tile, del_selected_tile, "selectedTile's docstring")
    tEstimation = property(get_t_estimation, set_t_estimation, del_t_estimation, "tEstimation's docstring")
    tStart = property(get_t_start, set_t_start, del_t_start, "tStart's docstring")
    nrTiles = property(get_nr_tiles, set_nr_tiles, del_nr_tiles, "nrTiles's docstring")
    nrProcs = property(get_nr_procs, set_nr_procs, del_nr_procs, "nrProcs's docstring")
    processingStatusFn = property(get_processing_status_fn, set_processing_status_fn, del_processing_status_fn, "processingStatusFn's docstring")
    processingEstimationFn = property(get_processing_estimation_fn, set_processing_estimation_fn, del_processing_estimation_fn, "processingEstimationFn's docstring")
    aerosolType = property(get_aerosol_type, set_aerosol_type, del_aerosol_type, "aerosolType's docstring")
    midLatitude = property(get_mid_latitude, set_mid_latitude, del_mid_latitude, "midLatitude's docstring")
    ozoneContent = property(get_ozone_content, set_ozone_content, del_ozone_content, "ozoneContent's docstring")
    waterVapour = property(get_water_vapour, set_water_vapour, del_water_vapour, "waterVapour's docstring")
    ozoneMean = property(get_ozone_mean, set_ozone_mean, del_ozone_mean, "ozoneMean's docstring")
    
    def createOrUpdateL2A_UserProduct(self):
        firstInit = False
        L1C_UP_MASK = '*1C_*'
        L1C_UP_DIR = self.sourceDir
        if os.path.exists(L1C_UP_DIR) == False:
            self.logger.fatal('directory "%s" does not exist.' % L1C_UP_DIR)
            return False

        # detect the filename for the datastrip metadata:
        L1C_DS_DIR = os.path.join(L1C_UP_DIR, 'DATASTRIP')
        if os.path.exists(L1C_DS_DIR) == False:
            self.logger.fatal('directory "%s" does not exist.' % L1C_DS_DIR)
            return False

        L1C_DS_MASK = '*_L1C_*'
        dirlist = sorted(os.listdir(L1C_DS_DIR))
        found = False
        
        for dirname in dirlist:
            if(fnmatch.fnmatch(dirname, L1C_DS_MASK) == True):
                found = True
                break
        
        if found == False:
            self.logger.fatal('No metadata in datastrip.')
            return False

        L1C_DS_DIR = os.path.join(L1C_DS_DIR, dirname)
        L1C_DS_MTD_XML = (dirname[:-7]+'.xml').replace('_MSI_', '_MTD_')
        self.L1C_DS_MTD_XML = os.path.join(L1C_DS_DIR, L1C_DS_MTD_XML)

        dirname, basename = os.path.split(L1C_UP_DIR)
        if(fnmatch.fnmatch(basename, L1C_UP_MASK) == False):
            self.logger.fatal('%s: identifier "*1C_*" is missing.' % basename)
            return False

        GRANULE =  os.path.join(L1C_UP_DIR, 'GRANULE')
        if os.path.exists(GRANULE) == False:
            self.logger.fatal('directory "%s" does not exist.' % GRANULE)
            return False

        #
        # the product (directory) structure:
        #-------------------------------------------------------
        L2A_UP_ID = basename[:4] + 'USER' + basename[8:]
        L2A_UP_ID = L2A_UP_ID.replace('1C_', '2A_')
        # SIITBX-55: alternative output directory for PDGS:
        targetDir = self._targetDirectory
        if targetDir != 'DEFAULT':
            dirname = targetDir
            if(os.path.exists(dirname) == False):
                os.mkdir(dirname)            
            
        L2A_UP_DIR = os.path.join(dirname, L2A_UP_ID)
        self.L2A_UP_DIR = L2A_UP_DIR
        self.L2A_UP_ID = L2A_UP_ID
        L1C_INSPIRE_XML = os.path.join(L1C_UP_DIR, 'INSPIRE.xml')
        L1C_MANIFEST_SAFE = os.path.join(L1C_UP_DIR, 'manifest.safe')
        L2A_INSPIRE_XML = os.path.join(L2A_UP_DIR, 'INSPIRE.xml')
        L2A_MANIFEST_SAFE = os.path.join(L2A_UP_DIR, 'manifest.safe')
        self.L2A_MANIFEST_SAFE = L2A_MANIFEST_SAFE

        AUX_DATA = 'AUX_DATA'
        DATASTRIP = 'DATASTRIP'
        GRANULE = 'GRANULE'
        HTML = 'HTML'
        REP_INFO = 'rep_info'
        QI_DATA = 'QI_DATA'

        if os.name == 'nt':
            #special treatment for windows for long pathnames:
            L1C_UP_DIR_        = u'\\'.join([u'\\\\?', L1C_UP_DIR])
            L1C_INSPIRE_XML_   = u'\\'.join([u'\\\\?', L1C_INSPIRE_XML])
            L1C_MANIFEST_SAFE_ = u'\\'.join([u'\\\\?', L1C_MANIFEST_SAFE])
            L2A_UP_DIR_        = u'\\'.join([u'\\\\?', L2A_UP_DIR])
            L2A_INSPIRE_XML_   = u'\\'.join([u'\\\\?', L2A_INSPIRE_XML])
            L2A_MANIFEST_SAFE_ = u'\\'.join([u'\\\\?', L2A_MANIFEST_SAFE])
        else:
            L1C_UP_DIR_ = L1C_UP_DIR
            L1C_INSPIRE_XML_ = L1C_INSPIRE_XML
            L1C_MANIFEST_SAFE_ = L1C_MANIFEST_SAFE
            L2A_UP_DIR_ = L2A_UP_DIR
            L2A_INSPIRE_XML_ = L2A_INSPIRE_XML
            L2A_MANIFEST_SAFE_ = L2A_MANIFEST_SAFE
            
        if(os.path.exists(os.path.join(L2A_UP_DIR, GRANULE)) == False):
            firstInit = True
            l.acquire()
            try:
                copytree(os.path.join(L1C_UP_DIR_, AUX_DATA),  os.path.join(L2A_UP_DIR_, AUX_DATA))
                copytree(os.path.join(L1C_UP_DIR_, DATASTRIP), os.path.join(L2A_UP_DIR_, DATASTRIP))
                copytree(os.path.join(L1C_UP_DIR_, HTML), os.path.join(L2A_UP_DIR_, HTML))
                copytree(os.path.join(L1C_UP_DIR_, REP_INFO), os.path.join(L2A_UP_DIR_, REP_INFO))    
                # remove the L1C xsds:
                S2_mask = 'S2_*.xsd'
                repdir = os.path.join(L2A_UP_DIR, REP_INFO)
                filelist = os.listdir(repdir)
                
                for filename in filelist:
                    if(fnmatch.fnmatch(filename, S2_mask) == True):
                        os.remove(os.path.join(repdir, filename))
                copyfile(L1C_MANIFEST_SAFE_, L2A_MANIFEST_SAFE_)
                copyfile(L1C_INSPIRE_XML_, L2A_INSPIRE_XML_)
                os.mkdir(os.path.join(L2A_UP_DIR_, GRANULE))
            finally:
                l.release()

        #create user product:
        S2A_mask = 'S2A_*.xml'
        try:
            filelist = sorted(os.listdir(L1C_UP_DIR))
            found = False
        except:
            self.logger.fatal('L1C input is not accessible')
            return False
        
        for filename in filelist:
            if(fnmatch.fnmatch(filename, S2A_mask) == True):
                found = True
                break
        if found == False:
            self.logger.fatal('No metadata for user product')    
            return False
    
        # prepare L2A User Product metadata file
        fn_L1C = os.path.join(L1C_UP_DIR, filename)
        fn_L2A = filename[:4] + 'USER' + filename[8:]
        fn_L2A = fn_L2A.replace('L1C_', 'L2A_')
        fn_L2A = os.path.join(L2A_UP_DIR, fn_L2A)
        self.L1C_UP_MTD_XML = fn_L1C
        self.L2A_UP_MTD_XML = fn_L2A
         
        # get the bandIndex    
        if self._resolution == 10:
            bandIndex = [1,2,3,7]
            self._ch940 = [0,0,0,0,0,0]
        else:
            bandIndex = [0,1,2,3,4,5,6,8,9,10,11,12]
            self._ch940 = [8,8,9,9,0,0]
        # this must always be initialized:
        xp = L2A_XmlParser(self, 'UP1C') 
        pic = xp.getTree('General_Info', 'Product_Image_Characteristics')        
        self._dnScale = float32(pic.QUANTIFICATION_VALUE.text)
        rc = pic.Reflectance_Conversion
        # The earth sun distance correction factor, already squared:
        self._d2 = float32(rc.U.text)     
        node = pic.Reflectance_Conversion.Solar_Irradiance_List
        i=0
        self._e0 = zeros(size(bandIndex), float32)
        for index in bandIndex:
            par = node.SOLAR_IRRADIANCE[index]
            if par is None: self.parNotFound(par)
            self._e0[i]= float32(par.text)/10
            i+=1

        if firstInit == True:
            # copy L2A schemes from config_dir into rep_info:    
            xp = L2A_XmlParser(self, 'GIPP')
            cs = xp.getRoot('Common_Section')
    
            dirname, basename = os.path.split(cs.UP_Scheme_2A.text)
            copyfile(os.path.join(self.configDir, cs.UP_Scheme_2A.text), os.path.join(L2A_UP_DIR, REP_INFO, basename))
            dirname, basename = os.path.split(cs.Tile_Scheme_2A.text)            
            copyfile(os.path.join(self.configDir, cs.Tile_Scheme_2A.text), os.path.join(L2A_UP_DIR, REP_INFO, basename))          
            dirname, basename = os.path.split(cs.DS_Scheme_2A.text)
            copyfile(os.path.join(self.configDir, cs.DS_Scheme_2A.text), os.path.join(L2A_UP_DIR, REP_INFO, basename))
            # copy L2A User Product metadata file:
            copyfile(fn_L1C, fn_L2A)
            # remove old L1C entries from L1C_UP_MTD_XML:
            xp = L2A_XmlParser(self, 'UP2A')
            if(xp.convert() == False):
                self.logger.fatal('error in converting user product metadata to level 2A')
    
            xp = L2A_XmlParser(self, 'UP2A')  
            pi = xp.getTree('General_Info', 'L2A_Product_Info')
            del pi.L2A_Product_Organisation.Granule_List[:]            
            # update L2A entries from L1C_UP_MTD_XML:
            pi.PRODUCT_URI = 'http://www.telespazio-vega.de'
            pi.PROCESSING_LEVEL = 'Level-2Ap'
            pi.PRODUCT_TYPE = 'S2MSI2Ap'
            dt = datetime.utcnow()
            pi.GENERATION_TIME = strftime('%Y-%m-%dT%H:%M:%SZ', dt.timetuple())
            pi.PREVIEW_IMAGE_URL = 'http://www.telespazio-vega.de'
            qo = pi.Query_Options
            qo.Aux_List.attrib['productLevel'] = 'Level-2Ap'
            pic = xp.getTree('General_Info', 'L2A_Product_Image_Characteristics')        
            L1C_TOA_QUANTIFICATION_VALUE =pic.L1C_L2A_Quantification_Values_List
            qvl = objectify.Element('L1C_L2A_Quantification_Values_List')
            qvl.L1C_TOA_QUANTIFICATION_VALUE = L1C_TOA_QUANTIFICATION_VALUE
            rc = pic.Reflectance_Conversion
            # The earth sun distance correction factor,
            # already squared:
            self._d2 = float32(rc.U.text)
            qvl.L2A_BOA_QUANTIFICATION_VALUE = str(int(self._dnScale))
            qvl.L2A_BOA_QUANTIFICATION_VALUE.attrib['unit'] = 'none'
            qvl.L2A_AOT_QUANTIFICATION_VALUE = str(self._L2A_AOT_QUANTIFICATION_VALUE)
            qvl.L2A_AOT_QUANTIFICATION_VALUE.attrib['unit'] = 'none'
            qvl.L2A_WVP_QUANTIFICATION_VALUE = str(self._L2A_WVP_QUANTIFICATION_VALUE)
            qvl.L2A_WVP_QUANTIFICATION_VALUE.attrib['unit'] = 'cm'
            pic.L1C_L2A_Quantification_Values_List = qvl
            
            scl = objectify.Element('L2A_Scene_Classification_List')
            scid = objectify.Element('L2A_Scene_Classification_ID')    
            scid.L2A_SCENE_CLASSIFICATION_TEXT = 'SC_NODATA'
            scid.L2A_SCENE_CLASSIFICATION_INDEX = '0'
            scl.append(scid)
            
            scid = objectify.Element('L2A_Scene_Classification_ID')              
            scid.L2A_SCENE_CLASSIFICATION_TEXT = 'SC_SATURATED_DEFECTIVE'
            scid.L2A_SCENE_CLASSIFICATION_INDEX = '1'
            scl.append(scid)
    
            scid = objectify.Element('L2A_Scene_Classification_ID')                  
            scid.L2A_SCENE_CLASSIFICATION_TEXT = 'SC_DARK_FEATURE_SHADOW'
            scid.L2A_SCENE_CLASSIFICATION_INDEX = '2'
            scl.append(scid)
    
            scid = objectify.Element('L2A_Scene_Classification_ID')                  
            scid.L2A_SCENE_CLASSIFICATION_TEXT = 'SC_CLOUD_SHADOW'
            scid.L2A_SCENE_CLASSIFICATION_INDEX = '3'
            scl.append(scid)        
            
            scid = objectify.Element('L2A_Scene_Classification_ID')      
            scid.L2A_SCENE_CLASSIFICATION_TEXT = 'SC_VEGETATION'
            scid.L2A_SCENE_CLASSIFICATION_INDEX = '4'
            scl.append(scid)
    
            scid = objectify.Element('L2A_Scene_Classification_ID')                      
            scid.L2A_SCENE_CLASSIFICATION_TEXT = 'SC_BARE_SOIL_DESERT'
            scid.L2A_SCENE_CLASSIFICATION_INDEX = '5'
            scl.append(scid)
                            
            scid = objectify.Element('L2A_Scene_Classification_ID')                  
            scid.L2A_SCENE_CLASSIFICATION_TEXT = 'SC_WATER'
            scid.L2A_SCENE_CLASSIFICATION_INDEX = '6'
            scl.append(scid)
                    
            scid = objectify.Element('L2A_Scene_Classification_ID')                  
            scid.L2A_SCENE_CLASSIFICATION_TEXT = 'SC_CLOUD_LOW_PROBA'
            scid.L2A_SCENE_CLASSIFICATION_INDEX = '7'
            scl.append(scid)
            
            scid = objectify.Element('L2A_Scene_Classification_ID')                              
            scid.L2A_SCENE_CLASSIFICATION_TEXT = 'SC_CLOUD_MEDIUM_PROBA'
            scid.L2A_SCENE_CLASSIFICATION_INDEX = '8'
            scl.append(scid)
            
            scid = objectify.Element('L2A_Scene_Classification_ID')                              
            scid.L2A_SCENE_CLASSIFICATION_TEXT = 'SC_CLOUD_HIGH_PROBA'
            scid.L2A_SCENE_CLASSIFICATION_INDEX = '9'
            scl.append(scid)
            
            scid = objectify.Element('L2A_Scene_Classification_ID')                  
            scid.L2A_SCENE_CLASSIFICATION_TEXT = 'SC_THIN_CIRRUS'
            scid.L2A_SCENE_CLASSIFICATION_INDEX = '10'
            scl.append(scid)
            
            scid = objectify.Element('L2A_Scene_Classification_ID')                            
            scid.L2A_SCENE_CLASSIFICATION_TEXT = 'SC_SNOW_ICE'
            scid.L2A_SCENE_CLASSIFICATION_INDEX = '11'
            scl.append(scid) 
            pic.append(scl)
    
            auxinfo = xp.getRoot('L2A_Auxiliary_Data_Info')
            auxdata = objectify.Element('Aux_Data')   
            gipp = objectify.Element('L2A_GIPP_List')
            auxdata.append(gipp)
            auxdata.L2A_PRODUCTION_DEM_TYPE = self.demReference
            auxdata.L2A_LIBRADTRAN_LUTS = self.atmDataFn
            auxdata.L2A_SNOW_CLIMATOLOGY = self.snowMapReference
            auxinfo.append(auxdata)
            xp.export()
    
        #create datastrip ID:
        S2A_mask = 'S2A_*'
        found = False
        L2A_DS_DIR = os.path.join(self._L2A_UP_DIR, DATASTRIP)
        dirlist = sorted(os.listdir(L2A_DS_DIR))
        for dirname in dirlist:
            if(fnmatch.fnmatch(dirname, S2A_mask) == True):
                
                L1C_DS_ID = dirname
                L2A_DS_ID = L1C_DS_ID[:4] + 'USER' + L1C_DS_ID[8:]
                L2A_DS_ID = L2A_DS_ID.replace('L1C_', 'L2A_')
                self.L2A_DS_ID = L2A_DS_ID
                
                olddir = os.path.join(L2A_DS_DIR, L1C_DS_ID)
                newdir = os.path.join(L2A_DS_DIR, L2A_DS_ID)
                    
                if firstInit == True:                    
                    os.rename(olddir, newdir)
            
                #find datastrip metadada, rename and change it:
                L2A_DS_SUBDIR = newdir
                filelist = sorted(os.listdir(L2A_DS_SUBDIR))
                found = False
                for filename in filelist:
                    if(fnmatch.fnmatch(filename, S2A_mask) == True):
                        L1C_DS_MTD_XML = filename
                        found = True
                        break
                if found == False:
                    self.logger.fatal('No metadata in datastrip.')

                L2A_DS_MTD_XML = L1C_DS_MTD_XML[:4] + 'USER' + L1C_DS_MTD_XML[8:]
                L2A_DS_MTD_XML = L2A_DS_MTD_XML.replace('L1C_', 'L2A_')
                oldfile = os.path.join(L2A_DS_SUBDIR, L1C_DS_MTD_XML)
                newfile = os.path.join(L2A_DS_SUBDIR, L2A_DS_MTD_XML)
                self.L2A_DS_MTD_XML = newfile
                
                if firstInit == True:
                    os.rename(oldfile, newfile)

                    xp = L2A_XmlParser(self, 'DS2A')
                    if(xp.convert() == False):
                        self.logger.fatal('error in converting datastrip metadata to level 2A')
                    ti = xp.getTree('Image_Data_Info', 'Tiles_Information')
                    try:
                        del ti.Tile_List.Tile[:]
                    except:
                        pass
                    xp.export()
                    found = True
    
        if found == False:
            self.logger.fatal('no subdirectory in datastrip')
        if self._selectedTile != None:
            self.nrProcs = 1
            return [self._selectedTile]
        else:
            return sorted(os.listdir(os.path.join(L1C_UP_DIR, GRANULE)))


    def create_L2A_Tile(self, tile):
        # generate new Tile ID:
        GRANULE = 'GRANULE' 
        L1C_TILE_ID = tile
        L2A_TILE_ID = L1C_TILE_ID[:4] + 'USER' + L1C_TILE_ID[8:]
        L2A_TILE_ID = L2A_TILE_ID.replace('L1C_', 'L2A_')
        self.L1C_TILE_ID = L1C_TILE_ID
        self.L2A_TILE_ID = L2A_TILE_ID
        L1C_TILE_ID = os.path.join(self.sourceDir, GRANULE, L1C_TILE_ID)
        L2A_TILE_ID = os.path.join(self.L2A_UP_DIR, GRANULE, L2A_TILE_ID)
        
        if(os.path.exists(L2A_TILE_ID) == False):
            os.mkdir(L2A_TILE_ID)
        else:
            # find the tile metadata:
            filelist = sorted(os.listdir(L2A_TILE_ID))
            filemask = 'S2A_USER_MTD_*.xml'
            for filename in filelist:
                if(fnmatch.fnmatch(filename, filemask) == True):
                    break
            
            L2A_TILE_MTD_XML = filename
            self.L2A_TILE_MTD_XML = os.path.join(L2A_TILE_ID, L2A_TILE_MTD_XML)
            L1C_TILE_MTD_XML = L2A_TILE_MTD_XML[:4] + 'OPER' + L2A_TILE_MTD_XML[8:]
            L1C_TILE_MTD_XML = L1C_TILE_MTD_XML.replace('L2A_', 'L1C_')
            self.L1C_TILE_MTD_XML = os.path.join(L1C_TILE_ID, L1C_TILE_MTD_XML)
        return L2A_TILE_ID


    def postprocess(self):
        if self.L2A_TILE_MTD_XML != None:
            try:
                l2aTileId = self.L2A_TILE_ID
                xp = L2A_XmlParser(self, 'UP2A')
                auxdata = xp.getTree('L2A_Auxiliary_Data_Info', 'Aux_Data')
                gipp = auxdata.L2A_GIPP_List
                gippFn = etree.Element('GIPP_FILENAME', type='GIP_Level-2Ap', version=self.processorVersion)
                gippFn.text = l2aTileId
                gipp.append(gippFn)
                xp.export()
                return True
            except:
                return False
        return True
    

    def setTimeEstimation(self, resolution):
        if self.selectedTile != None:
            nrTiles = 1
        else:
            nrTiles = self.nrTiles
        
        factor = float32(nrTiles)
        
        config = ConfigParser.RawConfigParser(allow_no_value=True)
        l.acquire()
        try:        
            config.read(self._processingEstimationFn)
            tEst60 = config.getfloat('time estimation','t_est_60') * factor
            tEst20 = config.getfloat('time estimation','t_est_20') * factor
            tEst10 = config.getfloat('time estimation','t_est_10') * factor
        finally:
            l.release()
        
        if(resolution == 60):
            self._tEstimation = tEst60
        elif(resolution == 20):
            self._tEstimation = tEst20
        elif(resolution == 10):
            self._tEstimation = (tEst20 + tEst10) 
        else:
            self._tEstimation = (tEst60 + tEst20 + tEst10)
        return


    def writeTimeEstimation(self, tMeasure):     
        l.acquire()
        try:
            config = ConfigParser.RawConfigParser()
            config.read(self._processingEstimationFn)
              
            if(self.resolution == 60):
                tEst = config.getfloat('time estimation','t_est_60')             
                tMeasureAsString = str((tEst + tMeasure) / 2.0 )
                config.set('time estimation','t_est_60', tMeasureAsString)
                
            elif(self.resolution == 20):
                tEst = config.getfloat('time estimation','t_est_20')             
                tMeasureAsString = str((tEst + tMeasure) / 2.0 )       
                config.set('time estimation','t_est_20', tMeasureAsString)
         
            elif(self.resolution == 10):
                tEst = config.getfloat('time estimation','t_est_10')             
                tMeasureAsString = str((tEst + tMeasure) / 2.0 )
                config.set('time estimation','t_est_10', tMeasureAsString)
    
            configFile = open(self._processingEstimationFn, 'w')
            config.write(configFile)
            configFile.close()
        finally:
            l.release()

    def timestamp(self, procedure):
        import multiprocessing
        p = multiprocessing.current_process()
        l.acquire()
        try:
            tNow = datetime.now()
            tDelta = tNow - self._timestamp
            self._timestamp = tNow
            if(self.logger.getEffectiveLevel()  != logging.NOTSET):
                self.logger.info('Procedure: ' + procedure + ', elapsed time[s]: %0.3f' % tDelta.total_seconds())

            f = open(self._processingStatusFn, 'r')
            tTotal = float(f.readline()) * 0.01
            f.close()
            increment = tDelta.total_seconds() / self._tEstimation
            tTotal += increment
            if tTotal > 1.0:
                tWeighted = 100.0 - exp(-tTotal)
            elif tTotal > 0.98:
                tWeighted = tTotal * 100.0 - exp(-tTotal)
            else:
                tWeighted = tTotal * 100.0

            stdoutWrite('Progress[%%]: %03.2f : ' % tWeighted)
            stdoutWrite('PID-%d, %s, elapsed time[s]: %0.3f\n' % (p.pid, procedure, tDelta.total_seconds()))
            
            f = open(self._processingStatusFn, 'w')
            f.write(str(tWeighted) + '\n')
            f.close()                
        except:
            f = open(self._processingStatusFn, 'w')
            f.write('0.0\n')
            f.close()                
        finally:
            l.release()
        return


    def parNotFound(self, parameter):
        basename = os.path.basename(self._configFn)
        self.logger.fatal('Configuration parameter <%s> not found in %s' % (parameter, basename))
        return False


    def readPreferences(self):
        ### Classificators
    ####READING FROM L2A_CAL_SC_GIPP.xml
        xp=L2A_XmlParser(self,'SC_GIPP')
        xp.export()
        xp.validate()

    ### Snow_map_reference
        node = xp.getTree('Scene_Classification', 'References')
        
        par = node.Snow_Map
        if par is None: self.parNotFound(node)
        self.snowMapReference = par.text          
        
        node = xp.getTree('Scene_Classification', 'Classificators')

        par = node.NO_DATA
        if par is None: self.parNotFound(node)
        self.noData = int32(par.text)  
        
        par = node.SATURATED_DEFECTIVE
        if par is None: self.parNotFound(node)
        self.saturatedDefective = int32(par.text)  

        par = node.DARK_FEATURES
        if par is None: self.parNotFound(node)
        self.darkFeatures = int32(par.text)  
        
        par = node.CLOUD_SHADOWS
        if par is None: self.parNotFound(node)
        self.cloudShadows = int32(par.text)  

        par = node.VEGETATION
        if par is None: self.parNotFound(node)
        self.vegetation = int32(par.text)  

        par = node.BARE_SOILS
        if par is None: self.parNotFound(node)
        self.bareSoils = int32(par.text)  

        par = node.WATER
        if par is None: self.parNotFound(node)
        self.water = int32(par.text)  
        
        par = node.LOW_PROBA_CLOUDS
        if par is None: self.parNotFound(node)
        self.lowProbaClouds = int32(par.text)  
        
        par = node.MEDIUM_PROBA_CLOUDS
        if par is None: self.parNotFound(node)
        self.medProbaClouds = int32(par.text)
        
        par = node.HIGH_PROBA_CLOUDS
        if par is None: self.parNotFound(node)
        self.highProbaClouds = int32(par.text)
        
        par = node.THIN_CIRRUS
        if par is None: self.parNotFound(node)
        self.thinCirrus = int32(par.text)
        
        par = node.SNOW_ICE
        if par is None: self.parNotFound(node)
        self.snowIce = int32(par.text)
        
    ### Thresholds
        node = xp.getTree ('Scene_Classification', 'Thresholds')
        
        par = node.T1_B02
        if par is None: self.parNotFound(node)
        self.T1_B02 = float32(par.text)
        
        par = node.T2_B02
        if par is None: self.parNotFound(node)
        self.T2_B02 = float32(par.text)
        
        par = node.T1_B04
        if par is None: self.parNotFound(node)
        self.T1_B04 = float32(par.text)
        
        par = node.T2_B04
        if par is None: self.parNotFound(node)
        self.T2_B04 = float32(par.text)
        
        par = node.T1_B8A
        if par is None: self.parNotFound(node)
        self.T1_B8A = float32(par.text)
        
        par = node.T2_B8A
        if par is None: self.parNotFound(node)
        self.T2_B8A = float32(par.text)
        
        par = node.T1_B10
        if par is None: self.parNotFound(node)
        self.T1_B10 = float32(par.text)
        
        par = node.T2_B10
        if par is None: self.parNotFound(node)
        self.T2_B10 = float32(par.text)
        
        par = node.T1_B12
        if par is None: self.parNotFound(node)
        self.T1_B12 = float32(par.text)
        
        par = node.T2_B12
        if par is None: self.parNotFound(node)
        self.T2_B12 = float32(par.text)
        
        par = node.T1_NDSI_CLD
        if par is None: self.parNotFound(node)
        self.T1_NDSI_CLD = float32(par.text)
        
        par = node.T2_NDSI_CLD
        if par is None: self.parNotFound(node)
        self.T2_NDSI_CLD = float32(par.text)
        
        par = node.T1_NDSI_SNW
        if par is None: self.parNotFound(node)
        self.T1_NDSI_SNW = float32(par.text)
        
        par = node.T2_NDSI_SNW
        if par is None: self.parNotFound(node)
        self.T2_NDSI_SNW = float32(par.text)
        
        par = node.T1_R_B02_B04
        if par is None: self.parNotFound(node)
        self.T1_R_B02_B04 = float32(par.text)
        
        par = node.T2_R_B02_B04
        if par is None: self.parNotFound(node)
        self.T2_R_B02_B04 = float32(par.text)
        
        par = node.T1_R_B8A_B03
        if par is None: self.parNotFound(node)
        self.T1_R_B8A_B03 = float32(par.text)
        
        par = node.T2_R_B8A_B03
        if par is None: self.parNotFound(node)
        self.T2_R_B8A_B03 = float32(par.text)
        
        par = node.T1_R_B8A_B11
        if par is None: self.parNotFound(node)
        self.T1_R_B8A_B11 = float32(par.text)
        
        par = node.T2_R_B8A_B11
        if par is None: self.parNotFound(node)
        self.T2_R_B8A_B11 = float32(par.text)
        
        par = node.T1_SNOW
        if par is None: self.parNotFound(node)
        self.T1_SNOW = float32(par.text)
        
        par = node.T2_SNOW
        if par is None: self.parNotFound(node)
        self.T2_SNOW = float32(par.text)
        
        par = node.T1_NDVI
        if par is None: self.parNotFound(node)
        self.T1_NDVI = float32(par.text)
        
        par = node.T2_NDVI
        if par is None: self.parNotFound(node)
        self.T2_NDVI = float32(par.text)
        
        par = node.T1_R_B8A_B03
        if par is None: self.parNotFound(node)
        self.T1_R_B8A_B03 = float32(par.text)
        
        par = node.T2_R_B8A_B03
        if par is None: self.parNotFound(node)
        self.T2_R_B8A_B03 = float32(par.text)
        
        par = node.T11_B02
        if par is None: self.parNotFound(node)
        self.T11_B02 = float32(par.text)
        
        par = node.T12_B02
        if par is None: self.parNotFound(node)
        self.T12_B02 = float32(par.text)
        
        par = node.T11_R_B02_B11
        if par is None: self.parNotFound(node)
        self.T11_R_B02_B11 = float32(par.text)
        
        par = node.T12_R_B02_B11
        if par is None: self.parNotFound(node)
        self.T12_R_B02_B11 = float32(par.text)
        
        par = node.T21_B12
        if par is None: self.parNotFound(node)
        self.T21_B12 = float32(par.text)
        
        par = node.T22_B12
        if par is None: self.parNotFound(node)
        self.T22_B12 = float32(par.text)
        
        par = node.T21_R_B02_B11
        if par is None: self.parNotFound(node)
        self.T21_R_B02_B11 = float32(par.text)
        
        par = node.T22_R_B02_B11
        if par is None: self.parNotFound(node)
        self.T22_R_B02_B11 = float32(par.text)

        par = node.T_CLOUD_LP
        if par is None: self.parNotFound(node)
        self.T_CLOUD_LP = float32(par.text)
        
        par = node.T_CLOUD_MP
        if par is None: self.parNotFound(node)
        self.T_CLOUD_MP = float32(par.text)
        
        par = node.T_CLOUD_HP
        if par is None: self.parNotFound(node)
        self.T_CLOUD_HP = float32(par.text)
        
        par = node.T1_B10
        if par is None: self.parNotFound(node)
        self.T1_B10 = float32(par.text)
        
        par = node.T2_B10
        if par is None: self.parNotFound(node)
        self.T2_B10 = float32(par.text)

        par = node.T_SDW
        if par is None: self.parNotFound(node)
        self.T_SDW = float32(par.text)
        
        par = node.T_B02_B12
        if par is None: self.parNotFound(node)
        self.T_B02_B12 = float32(par.text)
          
    ###READING FROM L2A_CAL_AC_GIPP.xml
    ### Scaling:
        xp = L2A_XmlParser(self, 'AC_GIPP')
        xp.export()
        xp.validate()
        
        node = xp.getTree('Atmospheric_Correction', 'References')
        par = node.Lib_Dir
        if par is None: self.parNotFound(par)
        try:
            scriptDir = os.environ['SEN2COR_BIN']
        except:
            scriptDir = getScriptDir() 
        self.libDir = os.path.join(scriptDir,'lib')
        if self._resolution == 10:
            libDir = os.path.join(self.libDir,'10')
        else:
            libDir = os.path.join(self.libDir,'20_60')
        
        if self._resolution == 10:
            bandIndex = [1,2,3,7]
            self._ch940 = [0,0,0,0,0,0]
        else:
            bandIndex = [0,1,2,3,4,5,6,8,9,10,11,12]
            self._ch940 = [8,8,9,9,0,0]        
        
        sensor = xp.getTree('Atmospheric_Correction', 'Sensor')
        wavelength = sensor.Calibration.Band_List.wavelength
        i = 0
        self._c0 = zeros(size(bandIndex), float32)
        self._c1 = zeros(size(bandIndex), float32)
        self._wvlsen = zeros(size(bandIndex), float32)
        for index in bandIndex:
            self._c0[i] = float32(wavelength[index].attrib['c0'])
            self._c1[i] = float32(wavelength[index].attrib['c1'])
            self._wvlsen[i] = float32(wavelength[index].text)
            i+=1

        i = 0
        self._fwhm = zeros(size(bandIndex), float32)
        for index in bandIndex:
            par = sensor.Resolution.Band_List.fwhm[index]
            if par is None: self.parNotFound(par)
            self._fwhm[i] = float32(par.text)
            i+=1
        
    #####READING FROM L2A_GIPP.xml
        xp = L2A_XmlParser(self, 'GIPP')
        xp.export()
        xp.validate()

    ### Common_Section:
        node = xp.getRoot('Common_Section')
        if node is None: self.parNotFound(node)

        par = node.Log_Level
        if par is None: self.parNotFound(par)
        self._logLevel = par.text
        
        par = node.Nr_Processes
        if par is None: self.parNotFound(par)
        self.nrProcs = par.text
        
        # SIITBX-55: alternative output directory for PDGS:
        par = node.Target_Directory
        if par is None: self.parNotFound(par)
        self._targetDirectory = par.text  
        
        par = node.DEM_Directory
        if par is None: self.parNotFound(par)
        self._demDirectory = par.text
        
        par = node.DEM_Reference
        if par is None: self.parNotFound(par)
        self._demReference = par.text
        
    ### Scene Classification:    
    ### Filters:
        node = xp.getTree('Scene_Classification','Filters')
        if node is None: self.parNotFound(node)
        
        par = node.Median_Filter
        if par is None: self.parNotFound(node)
        self.medianFilter = int(par.text)        
        
    ### Atmospheric Correction:
    ### References:
        node = xp.getTree('Atmospheric_Correction', 'Look_Up_Tables')
        if node is None: self.parNotFound(node)
        
        if self.aerosolType == None:
            par = node.Aerosol_Type
            if par is None: self.parNotFound(par)
            self.aerosolType = par.text
        
        if self.midLatitude == None:
            par = node.Mid_Latitude
            if par is None: self.parNotFound(par)
            self.midLatitude = par.text
        
        if self.ozoneContent == None:
            par = node.Ozone_Content
            if par is None: self.parNotFound(par)
            self.ozoneContent = par.text
        
    ### Flags:
        node = xp.getTree('Atmospheric_Correction', 'Flags')
        if node is None: self.parNotFound(node)

        par = node.BRDF_Correction
        if par is None: self.parNotFound(par)
        self.ibrdf = int(par.text)

        par = node.BRDF_Lower_Bound
        if par is None: self.parNotFound(par)
        self.thr_g = float32(par.text)

        par = node.WV_Correction
        if par is None: self.parNotFound(par)
        self.iwaterwv = int(par.text)

        par = node.VIS_Update_Mode
        if par is None: self.parNotFound(par)
        self.npref = int(par.text)

        par = node.WV_Watermask
        if par is None: self.parNotFound(par)
        self.iwv_watermask = int(par.text)

        par = node.Cirrus_Correction
        if par is None: self.parNotFound(par)
        self.icirrus = int(par.text)

    ### Calibration:
        node = xp.getTree('Atmospheric_Correction', 'Calibration')
        if node is None: self.parNotFound(node)
        
        par = node.Adj_Km
        if par is None: self.parNotFound(par)
        self.adj_km = float32(par.text)

        par = node.Visibility
        if par is None: self.parNotFound(par)
        self.visibility = float32(par.text)

        par = node.Altitude
        if par is None: self.parNotFound(par)
        self.altit = float32(par.text)

        par = node.Smooth_WV_Map
        if par is None: self.parNotFound(par)
        self.smooth_wvmap = float32(par.text)
        if (self.smooth_wvmap < 0.0): self.smooth_wvmap = 0.0

        par = node.WV_Threshold_Cirrus
        if par is None: self.parNotFound(par)
        self.wv_thr_cirrus = clip(float32(par.text), 0.1, 1.0)
        return
   
   
    def readTileMetadata(self):
        xp = L2A_XmlParser(self, 'T2A')
        ang = xp.getTree('Geometric_Info', 'Tile_Angles')
        try:
            azimuthAnglesList = ang.Sun_Angles_Grid.Azimuth.Values_List.VALUES
            solaz_arr = self.getFloatArray(azimuthAnglesList)
        except:
            self.logger.warning('No azimuth angular values in tile metadata available, will be set to 0')
            solaz_arr = 0
        try:
            zenithAnglesList = ang.Sun_Angles_Grid.Zenith.Values_List.VALUES
            solze_arr = self.getFloatArray(zenithAnglesList)
        except:
            self.logger.warning('No zenith angular values in user metadata available, will be set to 0')
            solze_arr = 0
        # images may be not squared - this is the case for the current testdata used
        # angle arrays have to be adapted, otherwise the bilinear interpolation is misaligned.
        imgSizeList = xp.getTree('Geometric_Info', 'Tile_Geocoding')
        size = imgSizeList.Size
        sizelen = len(size)
        nrows = None
        ncols = None
        for i in range(sizelen):
            if int(size[i].attrib['resolution']) == self._resolution:
                nrows = int(size[i].NROWS)
                ncols = int(size[i].NCOLS)
                break

        if(nrows == None or ncols == None):
            self.logger.fatal('no image dimension in metadata specified, please correct')
        if(nrows < ncols):
            last_row = int(solaz_arr[0].size * float(nrows)/float(ncols) + 0.5)
            saa = solaz_arr[0:last_row,:]
            sza = solze_arr[0:last_row,:]
        elif(ncols < nrows):
            last_col = int(solaz_arr[1].size * float(ncols)/float(nrows) + 0.5)
            saa = solaz_arr[:,0:last_col]
            sza = solze_arr[:,0:last_col]
        else:
            saa = solaz_arr
            sza = solze_arr

        if(saa.max() < 0):
            saa *= -1
        self.solaz_arr = clip(saa, 0, 360.0)

        sza = absolute(sza)
        self.solze_arr = clip(sza, 0, 70.0)

        self.nrows = nrows
        self.ncols = ncols
        try:
            solze = float32(ang.Mean_Sun_Angle.ZENITH_ANGLE.text)
        except:
            self.logger.warning('No mean zenith angular values in tile metadata available, will be set to 0')
            solze = 0
        try:
            solaz = float32(ang.Mean_Sun_Angle.AZIMUTH_ANGLE.text)
        except:
            self.logger.warning('No mean azimuth angular values in tile metadata available, will be set to 0')
            solaz = 0

        self._solze = absolute(solze)
        if self.solze > 70.0:
            self.solze = 70.0

        if solaz < 0:
            solaz *= -1
        if solaz > 360.0:
            solaz = 360.0
        self.solaz = solaz

        #
        # ATCOR employs the Lamberts reflectance law and assumes a constant viewing angle per tile (sub-scene)
        # as this is not given, this is a workaround, which have to be improved in a future version
        #
        try:
            viewAnglesList = ang.Mean_Viewing_Incidence_Angle_List.Mean_Viewing_Incidence_Angle
            arrlen = len(viewAnglesList)
            vaa = zeros(arrlen, float32)
            vza = zeros(arrlen, float32)
            for i in range(arrlen):
                vaa[i] = float32(viewAnglesList[i].AZIMUTH_ANGLE.text)
                vza[i] = float32(viewAnglesList[i].ZENITH_ANGLE.text)
        except:
            self.logger.warning('No Mean_Viewing_Incidence_Angle values in tile metadata available, will be set to default values')
            viewAnglesList = 0
            vaa = zeros([2,2], float32)
            vza = zeros([2,2], float32)

        _min = vaa.min()
        _max = vaa.max()
        if _min < 0: _min += 360
        if _max < 0: _max += 360
        vaa_arr = array([_min,_min,_max,_max])
        self.vaa_arr = vaa_arr.reshape(2,2)

        _min = absolute(vza.min())
        _max = absolute(vza.max())
        if _min > 40.0: _min = 40.0
        if _max > 40.0: _max = 40.0
        vza_arr = array([_min,_min,_max,_max])
        self.vza_arr = vza_arr.reshape(2,2)
        return


    def _get_subNodes(self, node, valtype):
        count = int(node.attrib['count'])
        if(valtype == 'int'):
            arr = zeros([count], int)
        elif(valtype == 'float'):
            arr = zeros([count], float32)
        else:
            self.logger.error('wrong type declatarion: ' + type)
            self.parNotFound('wrong type declatarion: ' + type)

        i = 0
        for sub in node:
            if(valtype == 'int'):
                arr[i] = int(sub.text)
            else:
                arr[i] = float32(sub.text)
            i += 1
        return arr

    
    def getIntArray(self, node):
        nrows = len(node)
        if nrows < 0:
            return False

        ncols = len(node[0].split())
        a = zeros([nrows,ncols],dtype=int)

        for i in range(nrows):
            a[i,:] = array(node[i].split(),dtype(int))

        return a


    def getUintArray(self, node):
        nrows = len(node)
        if nrows < 0:
            return False

        ncols = len(node[0].split())
        a = zeros([nrows,ncols],dtype=uint)

        for i in range(nrows):
            a[i,:] = array(node[i].split(),dtype(uint))

        return a


    def getFloatArray(self, node):
        nrows = len(node)
        if nrows < 0:
            return False

        ncols = len(node[0].text.split())
        a = zeros([nrows,ncols],dtype=float32)

        for i in range(nrows):
            a[i,:] = nan_to_num(array(node[i].text.split(),dtype(float32)))

        return a


    def putArrayAsStr(self, a, node):
        set_printoptions(precision=6)
        if a.ndim == 1:
            nrows = a.shape[0]
            for i in range(nrows):
                node[i] = a[i],dtype=str
                        
        elif a.ndim == 2:
            nrows = a.shape[0]
            for i in range(nrows):
                aStr = array_str(a[i,:]).strip('[]')
                node[i] = aStr
            return True
        else:
            return False


    def getStringArray(self, node):
        nrows = len(node)
        if nrows < 0:
            return False
        
        ncols = len(node[0].text.split())
        a = zeros([nrows,ncols],dtype=str)
        
        for i in range(nrows):
            a[i,:] = array(node[i].text.split(),dtype(str))

        return a


    def _adapt(self, default, setpoint, theRange):
        setpoint *= 0.001 # convert to micron
        # check if valid range, allow a broader interval than default
        if (setpoint[0] < default[0] - theRange):
            default[0] -= theRange
            self.logger.info('Adaptation of band interval. Setpoint: ' + str(setpoint[0] * 1000.0) + ', new value: ' + str(default[0] * 1000.0))
        elif (setpoint[0] > default[1]):
            self.logger.info('Setpoint > upper limit, will be ignored! Setpoint: ' + str(setpoint[0] * 1000.0) + ', new value: ' + str(default[0] * 1000.0))
            pass
        else: default[0] = setpoint[0]

        if (setpoint[1] > default[1] + theRange):
            default[1] += theRange
            self.logger.info('Adaptation of band interval. Setpoint: ' + str(setpoint[1] * 1000.0) + ', new value: ' + str(default[1] * 1000.0))
        elif (setpoint[1] < default[0]):
            self.logger.info('Setpoint < lower limit, will be ignored! Setpoint: ' + str(setpoint[1] * 1000.0) + ', new value: ' + str(default[1] * 1000.0))
            pass
        else: default[1] = setpoint[1]
        return default


    def _getDoc(self):
        from xml.etree import ElementTree as ET
        try:
            tree = ET.parse(self.configFn)
        except Exception, inst:
            self.logger.exception("Unexpected error opening %s: %s", self.configFn, inst)
            self.logger.fatal('Error in XML document')
        doc = tree.getroot()
        return doc


    def getInt(self, label, key):
        doc = self._getDoc()
        parameter = label + '/' + key
        par = doc.find(parameter)
        if par is None: self.parNotFound(parameter)
        return int(par.text)
 
 
    def getFloat(self, label, key):
        doc = self._getDoc()
        parameter = label + '/' + key
        par = doc.find(parameter)
        if par is None: self.parNotFound(parameter)
        return float32(par.text)
 
 
    def getStr(self, label, key):
        doc = self._getDoc()
        parameter = label + '/' + key
        par = doc.find(parameter)
        if par is None: self.parNotFound(parameter)
        return par.text    


    def getNrTiles(self):
        nrTiles = 0
        filemask = 'S2A_*_L1C_*'
        GRANULE = 'GRANULE'
        granuleDir = os.path.join(self.workDir, GRANULE)
        filelist = sorted(os.listdir(granuleDir))
        for tile in filelist:      
            if fnmatch.fnmatch(tile, filemask) == False:
                continue
            nrTiles += 1
        return nrTiles

    
    def createAtmDataFilename(self):
        extension = '.atm'
        height = '99000_'
        if self.midLatitude == 'SUMMER':
            waterVapour = 'wv20_'
        elif self.midLatitude == 'WINTER':
            waterVapour = 'wv04_'
        elif self.midLatitude == 'AUTO':
            waterVapour = 'wv00_'
        else:
            waterVapour = 'wv20_' # default
                      
        par = self.aerosolType
        if par == 'RURAL':
            aerosolType = 'rura'
        elif par == 'MARITIME':
            aerosolType = 'mari'
        elif par == 'AUTO':
            aerosolType = 'auto'
        else:
            aerosolType = 'rura' # default
            
        if self.resolution == 10:
            libDir = os.path.join(self.libDir,'10')
        else:
            libDir = os.path.join(self.libDir,'20_60')            

        par = self.ozoneContent
        if par == '0':
            ozoneContent = self.assignOzoneContent()
        else:
            ozoneContent = par
        
        atmDataFn = ozoneContent + height + waterVapour + aerosolType + extension
        self.logger.info('generated file name for look up tables is: ' + atmDataFn)
        self.atmDataFn = os.path.join(libDir, atmDataFn)
        return True


    def assignOzoneContent(self):
        # get the ozone value from metadata:
        columns = None
        ozoneMeasured = self.ozoneMean
        if self.midLatitude == 'SUMMER':
            columns = {
                       "f": abs(250 - ozoneMeasured),
                       "g": abs(290 - ozoneMeasured),
                       "h": abs(331 - ozoneMeasured),
                       "i": abs(370 - ozoneMeasured),
                       "j": abs(410 - ozoneMeasured),
                       "k": abs(450 - ozoneMeasured)
                       }    
        elif self.midLatitude == 'WINTER':
            columns = {
                       "t": abs(250 - ozoneMeasured),
                       "u": abs(290 - ozoneMeasured),
                       "v": abs(330 - ozoneMeasured),
                       "w": abs(377 - ozoneMeasured),
                       "x": abs(420 - ozoneMeasured),
                       "y": abs(460 - ozoneMeasured)                       
                       }

        ozoneContent = min(columns, key=columns.get)
        return ozoneContent
