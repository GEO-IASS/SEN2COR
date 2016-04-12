#!/usr/bin/env python
'''
Created on Feb 24, 2012
@author: umuellerwilm
'''
import fnmatch
import warnings
import subprocess
import sys, os, tempfile, logging, shutil
import glob
from L2A_Library import *
from time import sleep
import glymur
from skimage.measure import block_reduce
from skimage.transform import resize as skit_resize
from tables import *
from numpy import *
from tables.description import *
from shutil import copyfile, copytree
from scipy.ndimage.interpolation import zoom
from lxml import etree, objectify
from L2A_XmlParser import L2A_XmlParser

from subprocess import CalledProcessError
from osgeo.gdal_array import BandReadAsArray
from multiprocessing import Lock
l = Lock()

try:
    from osgeo import gdal,osr
    from osgeo.gdalconst import *
    gdal.TermProgress = gdal.TermProgress_nocb
except ImportError:
    import gdal,osr
    from gdalconst import *

    
# SIITBX-47: to suppress user warning due to the fact that 
# http://trac.osgeo.org/gdal/ticket/5480 is not implemented
# in the current openJPEG driver for windows used by ANACONDA:
gdal.PushErrorHandler('CPLQuietErrorHandler')
gdal.UseExceptions()

class Particle(IsDescription):
    bandName = StringCol(8)
    projectionRef = StringCol(512)
    geoTransformation = Int32Col(shape=6)
    rasterXSize = UInt16Col()
    rasterYSize = UInt16Col()
    rasterCount = UInt8Col()


class L2A_Tables(object):
    def __init__(self, config):
        self._config = config
        self._logger = config.logger
        self._DEV0 = None
        self._SHELL = None
        L1C_UP_MASK = '*1C_*'

        AUX_DATA = 'AUX_DATA'
        IMG_DATA = 'IMG_DATA'
        QI_DATA = 'QI_DATA'
        GRANULE = 'GRANULE'
        
        self.aux_src = None

        if self._logger.level == logging.DEBUG:
            self._DEV0 = ''
            if os.name == 'posix':
                self._SHELL = True
            else:
                self._SHELL = False
        else:
            if os.name == 'posix':
                self._DEV0 = ' &>/dev/null'
                self._SHELL = True
            else:
                self._DEV0 = ''
                self._SHELL = False
        # Resolution:
        self._resolution = int(self.config.resolution)
        if(self._resolution == 10):
            self._bandIndex = [1,2,3,7]
            self._nBands = 4
            bandDir = 'R10m'
        elif(self._resolution == 20):
            self._bandIndex = [0,1,2,3,4,5,6,8,9,10,11,12]
            self._nBands = 9
            bandDir = 'R20m'
        elif(self._resolution == 60):
            self._bandIndex = [0,1,2,3,4,5,6,8,9,10,11,12]
            self._nBands = 11
            bandDir = 'R60m'

        BANDS = bandDir
        # generate new Tile ID:
        L1C_TILE_ID = os.path.join(config.sourceDir, GRANULE, config.L1C_TILE_ID)
        L2A_TILE_ID = os.path.join(config.L2A_UP_DIR, GRANULE, config.L2A_TILE_ID)
        self.logger.info('tile directory is: ' + L2A_TILE_ID)

        if os.name == 'nt':
            #special treatment for windows for long pathnames:
            L1C_TILE_ID_ = u'\\'.join([u'\\\\?', L1C_TILE_ID])
            L2A_TILE_ID_ = u'\\'.join([u'\\\\?', L2A_TILE_ID])
            configFn     = u'\\'.join([u'\\\\?', config.configFn])
        else:
            L1C_TILE_ID_ = L1C_TILE_ID
            L2A_TILE_ID_ = L2A_TILE_ID
            configFn = config.configFn
        
        qiData = os.path.join(L2A_TILE_ID_, QI_DATA)
        if(os.path.exists(qiData) == False):
            copytree(os.path.join(L1C_TILE_ID_, QI_DATA),  os.path.join(L2A_TILE_ID_, QI_DATA))
            filelist = sorted(os.listdir(L1C_TILE_ID_))
            found = False
            for filename in filelist:
                if(fnmatch.fnmatch(filename, L1C_UP_MASK) == True):
                    found = True
                    break
            if found == False:
                self.logger.fatal('No metadata in tile')

            L1C_TILE_MTD_XML = os.path.join(L1C_TILE_ID_, filename)
            L2A_TILE_MTD_XML = filename
            L2A_TILE_MTD_XML = L2A_TILE_MTD_XML[:4] + 'USER' + L2A_TILE_MTD_XML[8:]
            L2A_TILE_MTD_XML = L2A_TILE_MTD_XML.replace('L1C_', 'L2A_')
            L2A_TILE_MTD_XML = os.path.join(L2A_TILE_ID_, L2A_TILE_MTD_XML)
            self.config.L1C_TILE_MTD_XML = L1C_TILE_MTD_XML
            self.config.L2A_TILE_MTD_XML = L2A_TILE_MTD_XML        

            xp = L2A_XmlParser(self.config, 'T1C')
            copyfile(L1C_TILE_MTD_XML, L2A_TILE_MTD_XML)
            xp = L2A_XmlParser(self.config, 'T2A')
            if(xp.convert() == False):
                self.logger.fatal('error in converting tile metadata to level 2A')

        tileId = self.config.L2A_TILE_ID
        if('L2A_CO_' in tileId):
            self.logger.info('consolidated tile ' + tileId + ': no entry in datastrip metadata generated')
        else:
            xp = L2A_XmlParser(self.config, 'DS2A')
            ti = xp.getTree('Image_Data_Info', 'Tiles_Information')
            try:
                Tile = objectify.Element('Tile', tileId = self.config.L2A_TILE_ID)
                ti.Tile_List.append(Tile)
                xp.export()
                xp.validate()
            except:
                self.logger.fatal('no Tile_List in datastrip metadata found')

        L1C_ImgDataDir = os.path.join(L1C_TILE_ID_, IMG_DATA)
        self._L2A_ImgDataDir = os.path.join(L2A_TILE_ID_, IMG_DATA)        

        self._L1C_bandDir = L1C_ImgDataDir
        self._L2A_bandDir = os.path.join(self._L2A_ImgDataDir, BANDS)

        if(os.path.exists(self._L2A_bandDir) == False):
            os.makedirs(self._L2A_bandDir)

        self._L1C_QualityMasksDir = os.path.join(L1C_TILE_ID_, QI_DATA)
        self._L2A_QualityDataDir = os.path.join(L2A_TILE_ID_, QI_DATA)
        self._L2A_AuxDataDir = os.path.join(L2A_TILE_ID_, AUX_DATA)

        if(os.path.exists(self._L2A_AuxDataDir) == False):
            copytree(os.path.join(L1C_TILE_ID_, AUX_DATA),  os.path.join(L2A_TILE_ID_, AUX_DATA))
            filelist = sorted(os.listdir(self._L2A_AuxDataDir))
            filemask = 'S2A_*'
            for filenameAux in filelist:
                if fnmatch.fnmatch(filenameAux, filemask) == True:
                    self.aux_src=filenameAux
                    break       
            # copy configuration to AUX dir:
            basename = os.path.basename(config._L2A_TILE_MTD_XML)
            fnAux = basename.replace('_MTD_', '_GIP_')
            target = os.path.join(self._L2A_AuxDataDir,fnAux)
            if(os.path.exists(target) == False):   
                try:
                    copyfile(configFn, target)
                except:
                    self.logger.error('cannot copy configuration to %s' % target)
                    
        filelist = sorted(os.listdir(self._L2A_AuxDataDir))
        filemask = 'S2A_*'
        for filenameAux in filelist:
            if fnmatch.fnmatch(filenameAux, filemask) == True:
                self.aux_src=filenameAux
                break   

        if(os.path.exists(self._L2A_QualityDataDir) == False):
            os.makedirs(self._L2A_QualityDataDir)  
        #
        # the File structure:
        #-------------------------------------------------------
        L2A_TILE_ID_SHORT = config.L2A_TILE_ID[:55]
        pre = L2A_TILE_ID_SHORT[:8]
        post = L2A_TILE_ID_SHORT[13:]
        self._L2A_Tile_BND_File = os.path.join(self._L2A_bandDir        , L2A_TILE_ID_SHORT + '_BXX_' + str(self._resolution) + 'm.jp2')
        self._L2A_Tile_VIS_File = os.path.join(self._L2A_bandDir        , pre + '_VIS_' + post + '_' + str(self._resolution) + 'm.jp2')
        self._L2A_Tile_AOT_File = os.path.join(self._L2A_bandDir        , pre + '_AOT_' + post + '_' + str(self._resolution) + 'm.jp2')
        self._L2A_Tile_WVP_File = os.path.join(self._L2A_bandDir        , pre + '_WVP_' + post + '_' + str(self._resolution) + 'm.jp2')
        self._L2A_Tile_SDW_File = os.path.join(self._L2A_bandDir        , pre + '_SDW_' + post + '_' + str(self._resolution) + 'm.jp2')
        self._L2A_Tile_SLP_File = os.path.join(self._L2A_bandDir        , pre + '_SLP_' + post + '_' + str(self._resolution) + 'm.jp2')
        self._L2A_Tile_ASP_File = os.path.join(self._L2A_bandDir        , pre + '_ASP_' + post + '_' + str(self._resolution) + 'm.jp2')
        self._L2A_Tile_DEM_File = os.path.join(self._L2A_AuxDataDir     , pre + '_DEM_' + post + '_' + str(self._resolution) + 'm.jp2')
        self._L2A_Tile_CLD_File = os.path.join(self._L2A_QualityDataDir , pre + '_CLD_' + post + '_' + str(self._resolution) + 'm.jp2')
        self._L2A_Tile_SNW_File = os.path.join(self._L2A_QualityDataDir , pre + '_SNW_' + post + '_' + str(self._resolution) + 'm.jp2')
        self._L2A_Tile_SCL_File = os.path.join(self._L2A_ImgDataDir     , pre + '_SCL_' + post + '_' + str(self._resolution) + 'm.jp2')
        self._L2A_Tile_PVI_File = os.path.join(self._L2A_QualityDataDir , pre + '_PVI_' + post + '_' + str(self._resolution) + 'm.jp2')
        self._L2A_Tile_DDV_File = os.path.join(self._L2A_QualityDataDir , pre + '_DDV_' + post + '_' + str(self._resolution) + 'm.jp2')

        self._ImageDataBase = os.path.join(self._L2A_bandDir, '.database.h5')
        self._TmpFile       = os.path.join(self._L2A_bandDir, '.tmpfile.tif')
        self._TmpDemFile    = os.path.join(self._L2A_bandDir, '.tmpdem.tif')
        self._TmpPwcFile    = os.path.join(self._L2A_bandDir, '.tmppwc.tif')
        self._TmpMslFile    = os.path.join(self._L2A_bandDir, '.tmpmsl.tif')
        self._TmpOzoFile    = os.path.join(self._L2A_bandDir, '.tmpozo.tif')
        self._acMode = False # default setting for scene classification

        # Geodata from image metadata:
        self._cornerCoordinates = None
        self._geoTransformation = None
        self._geoExtent = None
        self._projectionRef = None

        # Band Names:
        self._bandNames = ['B01','B02','B03','B04','B05','B06','B07','B08','B8A',\
                        'B09','B10','B11','B12','DEM','SCL','SNW','CLD','AOT',\
                        'WVP','VIS','SCM','PRV','ILU','SLP','ASP','HAZ','SDW',\
                        'DDV','HCW','ELE', 'PWC', 'MSL', 'OZO']

        # the mapping of the channels and bands
        self._B01 = 0
        self._B02 = 1
        self._B03 = 2
        self._B04 = 3
        self._B05 = 4
        self._B06 = 5
        self._B07 = 6
        self._B08 = 7
        self._B8A = 8
        self._B09 = 9
        self._B10 = 10
        self._B11 = 11
        self._B12 = 12
        self._DEM = 13
        self._SCL = 14
        self._SNW = 15
        self._CLD = 16
        self._AOT = 17
        self._WVP = 18
        self._VIS = 19
        self._SCM = 20
        self._PRV = 21
        self._ILU = 22
        self._SLP = 23
        self._ASP = 24
        self._HAZ = 25
        self._SDW = 26
        self._DDV = 27
        self._HCW = 28
        self._ELE = 29
        self._PWC = 30
        self._MSL = 31
        self._OZO = 32

        self.logger.debug('Module L2A_Tables initialized with resolution %d' % self._resolution)

    def get_ac_mode(self):
        return self._acMode


    def set_ac_mode(self, value):
        self._acMode = value


    def del_ac_mode(self):
        del self._acMode


    def get_logger(self):
        return self._logger


    def set_logger(self, value):
        self._logger = value


    def del_logger(self):
        del self._logger


    def get_corner_coordinates(self):
        return self._cornerCoordinates


    def get_geo_extent(self):
        return self._geoExtent


    def get_projection(self):
        return self._projection


    def set_corner_coordinates(self, value):
        self._cornerCoordinates = value


    def set_geo_extent(self, value):
        self._geoExtent = value


    def set_projection(self, value):
        self._projection = value


    def del_corner_coordinates(self):
        del self._cornerCoordinates


    def del_geo_extent(self):
        del self._geoExtent


    def del_projection(self):
        del self._projection


    def getBandNameFromIndex(self, index):
        return self._bandNames[index]


    def get_band_index(self):
        return self._bandIndex


    def get_n_bands(self):
        return self._nBands


    def get_db_name(self):
        return self._dbName


    def set_band_index(self, value):
        self._bandIndex = value


    def set_n_bands(self, value):
        self._nBands = value


    def set_db_name(self, value):
        self._dbName = value


    def del_band_index(self):
        del self._bandIndex


    def del_n_bands(self):
        del self._nBands


    def del_db_name(self):
        del self._dbName


        # end mapping of channels and bands

    def __exit__(self):
        sys.exit(-1)


    def __del__(self):
        self.logger.debug('Module L2A_Tables deleted')


    def get_config(self):
        return self._config


    def set_config(self, value):
        self._config = value


    def del_config(self):
        del self._config


    def get_b01(self):
        return self._B01


    def get_b02(self):
        return self._B02


    def get_b03(self):
        return self._B03


    def get_b04(self):
        return self._B04


    def get_b05(self):
        return self._B05


    def get_b06(self):
        return self._B06


    def get_b07(self):
        return self._B07


    def get_b08(self):
        return self._B08


    def get_b8a(self):
        return self._B8A


    def get_b09(self):
        return self._B09


    def get_b10(self):
        return self._B10


    def get_b11(self):
        return self._B11


    def get_b12(self):
        return self._B12


    def get_dem(self):
        return self._DEM


    def get_scl(self):
        return self._SCL


    def get_qsn(self):
        return self._SNW


    def get_qcl(self):
        return self._CLD


    def get_aot(self):
        return self._AOT


    def get_wvp(self):
        return self._WVP


    def get_vis(self):
        return self._VIS


    def get_scm(self):
        return self._SCM


    def get_prv(self):
        return self._PRV
    
   
    def get_pwc(self):
        return self._PWC


    def get_msl(self):
        return self._MSL


    def get_ozo(self):
        return self._OZO


    def set_b01(self, value):
        self._B01 = value


    def set_b02(self, value):
        self._B02 = value


    def set_b03(self, value):
        self._B03 = value


    def set_b04(self, value):
        self._B04 = value


    def set_b05(self, value):
        self._B05 = value


    def set_b06(self, value):
        self._B06 = value


    def set_b07(self, value):
        self._B07 = value


    def set_b08(self, value):
        self._B08 = value


    def set_b8a(self, value):
        self._B8A = value


    def set_b09(self, value):
        self._B09 = value


    def set_b10(self, value):
        self._B10 = value


    def set_b11(self, value):
        self._B11 = value


    def set_b12(self, value):
        self._B12 = value


    def set_dem(self, value):
        self._DEM = value


    def set_scl(self, value):
        self._SCL = value


    def set_qsn(self, value):
        self._SNW = value


    def set_qcl(self, value):
        self._CLD = value


    def set_aot(self, value):
        self._AOT = value


    def set_wvp(self, value):
        self._WVP = value


    def set_vis(self, value):
        self._VIS = value


    def set_scm(self, value):
        self._SCM = value


    def set_prv(self, value):
        self._PRV = value
        
        
    def set_pwc(self, value):
        self._PWC = value


    def set_msl(self, value):
        self._MSL = value


    def set_ozo(self, value):
        self._OZO = value


    def del_b01(self):
        del self._B01


    def del_b02(self):
        del self._B02


    def del_b03(self):
        del self._B03


    def del_b04(self):
        del self._B04


    def del_b05(self):
        del self._B05


    def del_b06(self):
        del self._B06


    def del_b07(self):
        del self._B07


    def del_b08(self):
        del self._B08


    def del_b8a(self):
        del self._B8A


    def del_b09(self):
        del self._B09


    def del_b10(self):
        del self._B10


    def del_b11(self):
        del self._B11


    def del_b12(self):
        del self._B12


    def del_dem(self):
        del self._DEM


    def del_scl(self):
        del self._SCL


    def del_qsn(self):
        del self._SNW


    def del_qcl(self):
        del self._CLD


    def del_aot(self):
        del self._AOT


    def del_wvp(self):
        del self._WV


    def del_vis(self):
        del self._VIS


    def del_scm(self):
        del self._SCM


    def del_prv(self):
        del self._PRV
        
        
    def del_pwc(self):
        del self._PWC


    def del_msl(self):
        del self._MSL


    def del_ozo(self):
        del self._OZO


    def get_ilu(self):
        return self._ILU


    def get_slp(self):
        return self._SLP


    def get_asp(self):
        return self._ASP


    def set_ilu(self, value):
        self._ILU = value


    def set_slp(self, value):
        self._SLP = value


    def set_asp(self, value):
        self._ASP = value


    def del_ilu(self):
        del self._ILU


    def del_slp(self):
        del self._SLP


    def del_asp(self):
        del self._ASP


    def get_sdw(self):
        return self._SDW


    def set_sdw(self, value):
        self._SDW = value


    def del_sdw(self):
        del self._SDW


    def get_ddv(self):
        return self._DDV


    def set_ddv(self, value):
        self._DDV = value


    def del_ddv(self):
        del self._DDV

    def get_hcw(self):
        return self._HCW


    def get_ele(self):
        return self._ELE


    def set_hcw(self, value):
        self._HCW = value


    def set_ele(self, value):
        self._ELE = value


    def del_hcw(self):
        del self._HCW


    def del_ele(self):
        del self._ELE


    B01 = property(get_b01, set_b01, del_b01, "B01's docstring")
    B02 = property(get_b02, set_b02, del_b02, "B02's docstring")
    B03 = property(get_b03, set_b03, del_b03, "B03's docstring")
    B04 = property(get_b04, set_b04, del_b04, "B04's docstring")
    B05 = property(get_b05, set_b05, del_b05, "B05's docstring")
    B06 = property(get_b06, set_b06, del_b06, "B06's docstring")
    B07 = property(get_b07, set_b07, del_b07, "B07's docstring")
    B08 = property(get_b08, set_b08, del_b08, "B08's docstring")
    B8A = property(get_b8a, set_b8a, del_b8a, "B8A's docstring")
    B09 = property(get_b09, set_b09, del_b09, "B09's docstring")
    B10 = property(get_b10, set_b10, del_b10, "B10's docstring")
    B11 = property(get_b11, set_b11, del_b11, "B11's docstring")
    B12 = property(get_b12, set_b12, del_b12, "B12's docstring")
    DEM = property(get_dem, set_dem, del_dem, "DEM's docstring")
    SCL = property(get_scl, set_scl, del_scl, "SCL's docstring")
    SNW = property(get_qsn, set_qsn, del_qsn, "SNW's docstring")
    CLD = property(get_qcl, set_qcl, del_qcl, "CLD's docstring")
    AOT = property(get_aot, set_aot, del_aot, "AOT's docstring")
    WVP = property(get_wvp, set_wvp, del_wvp, "WVP's docstring")
    VIS = property(get_vis, set_vis, del_vis, "VIS's docstring")
    SCM = property(get_scm, set_scm, del_scm, "SCM's docstring")
    PRV = property(get_prv, set_prv, del_prv, "PRV's docstring")
    PWC = property(get_pwc, set_pwc, del_pwc, "PWC's docstring")
    MSL = property(get_msl, set_msl, del_msl, "MSL's docstring")
    OZO = property(get_ozo, set_ozo, del_ozo, "OZO's docstring")
    ILU = property(get_ilu, set_ilu, del_ilu, "ILU's docstring")
    SLP = property(get_slp, set_slp, del_slp, "SLP's docstring")
    SDW = property(get_sdw, set_sdw, del_sdw, "SDW's docstring")
    ASP = property(get_asp, set_asp, del_asp, "ASP's docstring")
    DDV = property(get_ddv, set_ddv, del_ddv, "DDV's docstring")
    HCW = property(get_hcw, set_hcw, del_hcw, "HCW's docstring")
    ELE = property(get_ele, set_ele, del_ele, "ELE's docstring")
    config = property(get_config, set_config, del_config, "config's docstring")
    logger = property(get_logger, set_logger, del_logger, "logger's docstring")
    bandIndex = property(get_band_index, set_band_index, del_band_index, "bandIndex's docstring")
    nBands = property(get_n_bands, set_n_bands, del_n_bands, "nBands's docstring")
    dbName = property(get_db_name, set_db_name, del_db_name, "dbName's docstring")
    cornerCoordinates = property(get_corner_coordinates, set_corner_coordinates, del_corner_coordinates, "cornerCoordinates's docstring")
    geoExtent = property(get_geo_extent, set_geo_extent, del_geo_extent, "geoExtent's docstring")
    projection = property(get_projection, set_projection, del_projection, "projection's docstring")
    acMode = property(get_ac_mode, set_ac_mode, del_ac_mode, "acMode's docstring")

    def ReprojectCoords(self,coords,src_srs,tgt_srs):
        ''' Reproject a list of x,y coordinates.
            @type geom:     C{tuple/list}
            @param geom:    List of [[x,y],...[x,y]] coordinates
            @type src_srs:  C{osr.SpatialReference}
            @param src_srs: OSR SpatialReference object
            @type tgt_srs:  C{osr.SpatialReference}
            @param tgt_srs: OSR SpatialReference object
            @rtype:         C{tuple/list}
            @return:        List of transformed [[x,y],...[x,y]] coordinates
        '''
        trans_coords=[]
        transform = osr.CoordinateTransformation( src_srs, tgt_srs)
        for x,y in coords:
            x,y,z = transform.TransformPoint(x,y)
            trans_coords.append([x,y])
        return trans_coords


    def GetExtent(self, gt,cols,rows):
        ''' Return list of corner coordinates from a geotransform

            @type gt:   C{tuple/list}
            @param gt: geotransform
            @type cols:   C{int}
            @param cols: number of columns in the dataset
            @type rows:   C{int}
            @param rows: number of rows in the dataset
            @rtype:    C{[float,...,float]}
            @return:   coordinates of each corner
        '''
        ext=[]
        xarr=[0,cols]
        yarr=[0,rows]

        for px in xarr:
            for py in yarr:
                x=gt[0]+(px*gt[1])+(py*gt[2])
                y=gt[3]+(px*gt[4])+(py*gt[5])
                ext.append([x,y])
            yarr.reverse()
        return ext


    def getResolutionIndex(self):
        res = self._config.resolution
        if res == 10:
            return 0
        elif res == 20:
            return 1
        elif res == 60:
            return 2
        else:
            return False


    def checkAotMapIsPresent(self, resolution):
        sourceDir = os.path.join(self._L2A_ImgDataDir, 'R' + str(resolution) + 'm')
        try:
            dirs = sorted(os.listdir(sourceDir))
            filemask = '*_AOT_L2A_*' + '.jp2'
            for filename in dirs:
                if fnmatch.fnmatch(filename, filemask) == True:
                    return True
            return False
        except:
            return False
        
    
    def checkBandCount(self):    
        sourceDir = self._L1C_bandDir
        dirs = sorted(os.listdir(sourceDir))
        bandIndex = self.bandIndex
        bandCount = 0
        for i in bandIndex:
            for filename in dirs:
                bandName = self.getBandNameFromIndex(i)
                filemask = '*_L1C_*_%3s.jp2' % bandName
                if fnmatch.fnmatch(filename, filemask) == False:
                    continue
                bandCount += 1
                break
        if len(bandIndex) > bandCount:
            return False
        return True


    def importBandList(self):
        # convert JPEG-2000 input files to H5 file format
        # initialize H5 database for usage:
        sourceDir = self._L1C_bandDir
        database = os.path.join(sourceDir, self._ImageDataBase)
        rasterX = False
        if(os.path.isfile(database)):
            os.remove(database)
            self.logger.info('Old database removed')
        self.config.timestamp('L2A_Tables: start import')

        dirs = sorted(os.listdir(sourceDir))
        bandIndex = self.bandIndex
        for i in bandIndex:
            for filename in dirs:
                bandName = self.getBandNameFromIndex(i)
                filemask = '*_L1C_*_%3s.jp2' % bandName
                if fnmatch.fnmatch(filename, filemask) == False:
                    continue
                if(rasterX == False):
                    self.setCornerCoordinates()
                    rasterX = True
                res = self.importBand(i, os.path.join(sourceDir, filename))
                if res == False:
                    return False
                break

        upsampling = False
        # 10m bands only: perform an up sampling of SCL, AOT, and VIS from 20 m channels to 10
        if(self._resolution == 10):
            self.logger.info('perform up sampling of SCL, AOT and VIS from 20m channels to 10m')
            srcResolution = '_20m'
            channels = [17,19]
            sourceDir = self._L2A_bandDir.replace('R10m', 'R20m')
            upsampling = True
        
        if upsampling == True:
            dirs = sorted(os.listdir(sourceDir))
            for i in channels:
                for filename in dirs:
                    bandName = self.getBandNameFromIndex(i)
                    if (bandName == 'VIS') or (bandName == 'AOT') or (bandName == 'WVP'):
                        filemask = '*_' + bandName + '_L2A_*' + srcResolution + '.jp2'
                    if fnmatch.fnmatch(filename, filemask) == False:
                        continue
                    res = self.importBand(i, os.path.join(sourceDir, filename))
                    if res == False:
                        return False
                    break
            # scene class is in different directory:
            sourceDir = self._L2A_ImgDataDir
            dirs = sorted(os.listdir(sourceDir))
            bandName = 'SCL'
            channel = 14
            filemask = '*_' + bandName + '_L2A_*' + srcResolution + '.jp2'
            for filename in dirs:                        
                if fnmatch.fnmatch(filename, filemask) == False:
                    continue
                res = self.importBand(channel, os.path.join(sourceDir, filename))
                if res == False:
                    return False
                break

        self.dem = False        
        demDir =  self.config.demDirectory
        if demDir == 'NONE':
            self.logger.info('DEM directory not specified, flat surface is used')
            return True
        
        demDir = os.path.join(self.config.home, self.config.demDirectory)
        # check if DEM is a DTED type, these files must exist in the given directory:
        DTED = False
        filemask = 'e*.dt1'
        if os.path.exists(demDir):
            files = sorted(os.listdir(demDir))
            for filename in files:
                if fnmatch.fnmatch(filename, filemask) == True:
                    DTED = True
                    break
            
        # yes it is, run dem preparation for DTED:
        if DTED == True:
            if(self.gdalDEM_dted() == False):
                self.config.demDirectory = 'NONE'
                # continue with flat surface ...
                return True
                      
        else: # run DEM preparation for SRTM:
            if(self.gdalDEM_srtm() == False):
                self.config.demDirectory = 'NONE'
                # continue with flat surface ...
                return True

        # generate hill shadow, slope and aspect using DEM:
        if(self.gdalDEM_Shade() == False):
            self.logger.fatal('shell execution error generating DEM shadow')
            return False
        
        if(self.gdalDEM_Slope() == False):
            self.logger.fatal('shell execution error generating DEM slope')
            return False
        
        if(self.gdalDEM_Aspect() == False):
            self.logger.fatal('shell execution error generating DEM aspect')
            return False

        try:
            os.remove(self._TmpDemFile)
        except:
            pass
        
        self.config.timestamp('L2A_Tables: stop import')

        return True


    def setCornerCoordinates(self):
        # get the target resolution and metadata for the resampled bands below:
        xp = L2A_XmlParser(self.config, 'T2A')           
        tg = xp.getTree('Geometric_Info', 'Tile_Geocoding')
        nrows = self.config.nrows
        ncols = self.config.ncols
        idx = self.getResolutionIndex()
        ulx = tg.Geoposition[idx].ULX
        uly = tg.Geoposition[idx].ULY        
        res = float32(self.config.resolution)
        geoTransformation = [ulx,res,0.0,uly,0.0,-res]
        extent = self.GetExtent(geoTransformation, ncols, nrows)
        self._cornerCoordinates = asarray(extent)
        return


    def get_utm_zone(self, longitude):
        return (int(1+(longitude+180.0)/6.0))


    def is_northern(self, latitude): # Determines if given latitude is a northern for UTM
            if (latitude < 0.0):
                return 0
            else:
                return 1


    def transform_utm_to_wgs84(self, easting, northing, zone1, zone2):
        utm_coordinate_system = osr.SpatialReference() #Create a new spatial reference object using a named parameter
        utm_coordinate_system.SetWellKnownGeogCS("WGS84") # Set geographic coordinate system to handle lat/lon
        zone = zone1
        hemi = zone2
        # SIITBX-48:
        if(hemi == 'N'): # N is Northern Hemisphere
            utm_coordinate_system.SetUTM(zone, 1) # call sets detailed projection transformation parameters
        else:
            utm_coordinate_system.SetUTM(zone, 0)
        wgs84_coordinate_system = utm_coordinate_system.CloneGeogCS() # Clone ONLY the geographic coordinate system
        # create transform component
        utm_to_wgs84_geo_transform = osr.CoordinateTransformation(utm_coordinate_system, wgs84_coordinate_system)
        return utm_to_wgs84_geo_transform.TransformPoint(easting, northing, 0) # returns lon, lat, altitude


    def transform_wgs84_to_utm(self, lon, lat):
        utm_coordinate_system = osr.SpatialReference()
        utm_coordinate_system.SetWellKnownGeogCS("WGS84") # Set geographic coordinate system to handle lat/lon
        utm_coordinate_system.SetUTM(self.get_utm_zone(lon), self.is_northern(lat))
        wgs84_coordinate_system = utm_coordinate_system.CloneGeogCS() # Clone ONLY the geographic coordinate system
        # create transform component
        wgs84_to_utm_geo_transform = osr.CoordinateTransformation(wgs84_coordinate_system, utm_coordinate_system)
        return wgs84_to_utm_geo_transform.TransformPoint(lon, lat, 0) # returns easting, northing, altitude


    def getAuxData(self, bandIndex):
        '''
        PWC (Precipitable Water Content), Grib Unit [kg/m^2]
        MSL (Mean Sea Level pressure),    Grib Unit [Pa]
        OZO (Ozone),                      Grib Unit [kg/m^2]
        
        calculation for Ozone according to R. Richter (20/1/2016):
        ----------------------------------------------------------
        GRIB_UNIT = [kg/m^2]
        standard ozone column is 300 DU (Dobson Units),
        equals to an air column of 3 mm at STP (standard temperature (0 degree C) and pressure of 1013 mbar).
        
        Thus, molecular weight of O3 (M = 48): 2.24 g (equals to 22.4 liter at STP)
        
        300 DU = 3 mm  (equals to (0.3*48 / 2.24) [g/m^2])
         = 6.428 [g/m^2] = 6.428 E-3 [kg/m^2]        

        Example:
        
        ozone (GRIB) = 0.005738 (equals to DU = 300 * 0.005738/6.428 E-3)
        ozone (DU)   = 267.4 DU
        
        Thus, ozone GRIB will be weighted with factor 155.5694 (equals to 1/6.428 E-3)
        in order to receive ozone in DU
        '''
        auxBands = [self.PWC, self.MSL, self.OZO]
        if bandIndex in auxBands == False:
            self.logger.error('wrong band index for aux data')
            return False
        
        bandIndex -= 29 # bandIndex starts at 30
        ozoneFactor = 155.5694 # 1/6.428 E-3
        standardOzoneColumn = 300.0
        
        straux_src = os.path.join(self._L2A_AuxDataDir, self.aux_src)
        curdir = os.path.curdir
        head, tail = os.path.split(straux_src)
        l.acquire()
        os.chdir(head)
        arr = False
        while True:
            try:
                dataSet = gdal.Open(tail, GA_ReadOnly)
                band = dataSet.GetRasterBand(bandIndex)
                arr = BandReadAsArray(band)
                if bandIndex == 3: # recalculate to 300 DU:
                    arr = arr * standardOzoneColumn * ozoneFactor
                break
            except:
                self.logger.error('error in reading aux data from band %s' % self.getBandNameFromIndex(bandIndex))
            finally:
                os.chdir(curdir)
                l.release()
                return arr

    
    def gdalDEM_dted(self):
        demDir = self.config.demDirectory
        if demDir == 'NONE':
            self.logger.info('DEM directory not specified, flat surface is used')
            return False

        self.logger.info('Start DEM alignment for tile')
        sourceDir = os.path.join(self.config.home, demDir)

        xy = self.cornerCoordinates
        xp = L2A_XmlParser(self.config, 'T2A')
        tg = xp.getTree('Geometric_Info', 'Tile_Geocoding')
        hcsName = tg.HORIZONTAL_CS_NAME.text
        zone = hcsName.split()[4]
        zone1 = int(zone[:-1])
        zone2 = zone[-1:].upper()
        lonMin, latMin, dummy = self.transform_utm_to_wgs84(xy[1,0], xy[1,1], zone1, zone2)
        lonMax, latMax, dummy = self.transform_utm_to_wgs84(xy[3,0], xy[3,1], zone1, zone2)
        
        lat_cen = int((latMax+latMin)/2)
        lon_cen = int((lonMax+lonMin)/2)
        lonMin = int(lonMin)
        lonMax = int(lonMax)
        latMin = int(latMin)
        latMax = int(latMax)
        
        if lon_cen<0:
            lonMask='w'
            lon_cen=-lon_cen
        else:
            lonMask='e'
            
        if lat_cen<0:
            latMask='s'
            lat_cen=-lat_cen
        else:
            latMask='n'
            
        if lon_cen < 100:
            if lon_cen < 10:
                lonMask = lonMask + '00'
            else:
                lonMask = lonMask + '0'
        
        srtmf_src = ''
        
        filelist = sorted(os.listdir(sourceDir))
        found = False
        for i in range(lonMin-1,lonMax+1):
            for j in range(latMin-1,latMax+1):
                file_mask = lonMask + str(i) + '_' + latMask + str(j) + '.dt1'
                    
                for filename in filelist:
                    if(fnmatch.fnmatch(filename, file_mask) == True):
                        found = True
                        srtmf_src += ' ' + os.path.join(sourceDir, filename)  
                        break                      
                            
        if found == False:  
            self.logger.info('DEM not found, flat surface is used')
            return False

        hcsCode = tg.HORIZONTAL_CS_CODE.text
        t_srs = '-t_srs ' + hcsCode

        te = ' -te %f %f %f %f' % (xy[0,0], xy[2,1], xy[2,0], xy[0,1])
        tr = ' -tr %d %d' % (self.config.resolution, self.config.resolution)
        t_warp = te + tr + ' -r cubicspline '

        if(os.path.isfile(self._TmpFile) == True):
            os.remove(self._TmpFile)

        tmpFile = os.path.join(sourceDir, 'srtm_' + self.config.L2A_TILE_ID + '.tif')
        command = 'gdalwarp '
        arguments = '-ot Int16 '
        callstr = command + arguments + t_srs + t_warp + srtmf_src + ' ' + tmpFile + self._DEV0
        try:
            if(os.path.isfile(tmpFile) == True):
                os.remove(tmpFile)
            if(os.path.isfile(self._TmpFile) == True):
                os.remove(self._TmpFile)
            subprocess.call(callstr, shell=self._SHELL)
        except:
            stderrWrite('Error reading DEM, flat surface will be used.\n')
            self.logger.fatal('Error reading DEM, flat surface will be used')
            os.remove(srtmf_src)
            return False
        finally:
            l.release()
            
        while True:
            if(os.path.isfile(tmpFile) == True):
                break
            else:
                sleep(1)
                continue            

        copyfile(tmpFile, self._TmpFile)
        self.importBand(self.DEM)
        os.rename(self._TmpFile, self._TmpDemFile)
        os.remove(srtmf_src)
        self.logger.info('DEM received and prepared')     
        return True
    
    
    def gdalDEM_srtm(self):
        import urllib, urllib2
        import zipfile
        import multiprocessing
        p = multiprocessing.current_process()
 
        demDir = self.config.demDirectory
        if demDir == 'NONE':
            self.logger.info('DEM directory not specified, flat surface is used')
            return False
 
        self.logger.info('Start DEM alignment for tile')
        sourceDir = os.path.join(self.config.home, demDir)
        l.acquire()
        try:
            if(os.path.exists(sourceDir) == False):
                os.makedirs(sourceDir)
        finally:
            l.release()
 
        xy = self.cornerCoordinates
        xp = L2A_XmlParser(self.config, 'T2A')
        tg = xp.getTree('Geometric_Info', 'Tile_Geocoding')
        hcsName = tg.HORIZONTAL_CS_NAME.text
        zone = hcsName.split()[4]
        zone1 = int(zone[:-1])
        zone2 = zone[-1:].upper()
        lonMin, latMin, dummy = self.transform_utm_to_wgs84(xy[1,0], xy[1,1], zone1, zone2)
        lonMax, latMax, dummy = self.transform_utm_to_wgs84(xy[3,0], xy[3,1], zone1, zone2)
        lonMinId = int((-180-lonMin)/(-360)*72)+1
        lonMaxId = int((-180-lonMax)/(-360)*72)+1
        latMinId = int((60-latMax)/(120)*24)+1 # this is inverted by intention
        latMaxId = int((60-latMin)/(120)*24)+1 # this is inverted by intention

        if(lonMinId < 1) or (lonMaxId > 72) or (latMinId < 1) or (latMaxId > 24):
            self.logger.error('no SRTM dataset available for this tile, flat surface will be used')
            stdoutWrite('no SRTM dataset available for this tile, flat surface will be used')
            return False

        for i in range(lonMinId, lonMaxId+1):
            for j in range(latMinId, latMaxId+1):        
                tifFn = 'srtm_{:0>2d}_{:0>2d}.tif'.format(i,j)
                zipFn = 'srtm_{:0>2d}_{:0>2d}.zip'.format(i,j)
                tmpDir = 'srtm_{:0>2d}_{:0>2d}_tmp'.format(i,j)
                tifFn_path = os.path.join(sourceDir,tifFn)
                zipFn_path = os.path.join(sourceDir,zipFn)
                tmpDir = os.path.join(sourceDir,tmpDir)
                     
                while True:           
                    try:
                        # does the tiff file already exist?
                        with open (tifFn_path) as fp:
                            if self.logger.level == logging.DEBUG: 
                                print p.pid, 'DEM exist: %s', tifFn_path
                            break
                        
                    except IOError:
                        # block the zip download for the other processes:
                        if not os.path.isdir(tmpDir):
                            os.makedirs(tmpDir)
                            tifFn_src = os.path.join(tmpDir, tifFn)
                            try:
                                # zipfile needs to be downloaded ...
                                if self.logger.level == logging.DEBUG: 
                                    print p.pid, 'read zipfile ...'
                                prefix = self.config._demReference
                                stdoutWrite('Trying to retrieve DEM from URL %s this may take some time ...\n' % prefix)
                                self.logger.info('Trying to retrieve DEM from URL: %s', prefix)
                                url = prefix + zipFn
                                webFile = urllib.urlopen(url)
                                localFile = open(os.path.join(sourceDir, url.split('/')[-1]), 'wb')
                                localFile.write(webFile.read())
                                webFile.close()
                                localFile.close()
                                if self.logger.level == logging.DEBUG: 
                                    print p.pid, 'zip downloaded.'
                            except urllib2.URLError as e:
                                stderrWrite('Download error %s, flat surface will be used.\n' % e.reason)
                                self.logger.error('Download error %s, flat surface will be used' % e.reason)
                                return True

                            zipf = zipfile.ZipFile(zipFn_path, mode='r')
                            if(zipf.testzip() == None):
                                zipf.extractall(tmpDir)
                                zipf.close()
                                os.rename(tifFn_src, tifFn_path)
                                self.logger.info('DEM unpacked and moved: %s', tifFn_path)
                                if self.logger.level == logging.DEBUG: 
                                    print p.pid, 'DEM unpacked and moved: %s', tifFn_path
                                shutil.rmtree(tmpDir)
                                os.remove(zipFn_path)                    
                                if self.logger.level == logging.DEBUG: 
                                    print p.pid, 'zipfile removed: %s', zipf
                                break
                        else:
                            if self.logger.level == logging.DEBUG: 
                                print p.pid, 'wait for termination of download ...'
                            sleep(1)
 
        command = 'gdalwarp '
        arguments = '-ot Int16 '

        if(lonMinId == lonMaxId) & (latMinId == latMaxId):
            srtmf_src = os.path.join(sourceDir,'srtm_{:0>2d}_{:0>2d}.tif'.format(i,j))
        else:
            # more than 1 DEM needs to be concatenated:
            for i in range(lonMinId, lonMaxId+1):
                for j in range(latMinId, latMaxId+1):
                    tifFn = os.path.join(sourceDir,'srtm_{:0>2d}_{:0>2d}.tif'.format(i,j))
                    arguments += tifFn + ' '
 
            srtmf_src = os.path.join(sourceDir, 'srtm_' + self.config.L2A_TILE_ID + '_src.tif')
            if(os.path.isfile(srtmf_src) == True):
                os.remove(srtmf_src)
            callstr = command + arguments + srtmf_src  + self._DEV0
            l.acquire()
            try:
                subprocess.call(callstr, shell=self._SHELL)
            except:
                stderrWrite('shell execution error using gdalwarp.\n')
                self.logger.fatal('shell execution error using gdalwarp')
                return False
            finally:
                l.release()

        while True:
            try:
                with open (srtmf_src) as fp:
                    break
            except:
                if self.logger.level == logging.DEBUG:
                    print(p.pid, 'waiting for DEM extraction ...')
                sleep(1)
                                
        hcsCode = tg.HORIZONTAL_CS_CODE.text
        t_srs = '-t_srs ' + hcsCode
 
        te = ' -te %f %f %f %f' % (xy[0,0], xy[2,1], xy[2,0], xy[0,1])
        tr = ' -tr %d %d' % (self.config.resolution, self.config.resolution)
        t_warp = te + tr + ' -r cubicspline '
        arguments = '-ot Int16 '
        tmpFile = os.path.join(sourceDir, 'srtm_' + self.config.L2A_TILE_ID + '.tif')
        callstr = command + arguments + t_srs + t_warp + srtmf_src + ' ' + tmpFile + self._DEV0
        l.acquire()
        try:
            subprocess.call(callstr, shell=self._SHELL)
        except:
            stderrWrite('Error reading DEM, flat surface will be used.\n')
            self.logger.fatal('Error reading DEM, flat surface will be used')
            if((os.path.isfile(srtmf_src) == True) and ('_src.tif' in srtmf_src)):
                os.remove(srtmf_src)
            return False
        finally:
            l.release()

        while True:
            try:
                with open (tmpFile) as fp:
                    break
            except:
                sleep(1)

        if(os.path.isfile(self._TmpFile) == True):
            os.remove(self._TmpFile)
        copyfile(tmpFile, self._TmpFile)
        self.importBand(self.DEM)
        if(os.path.isfile(self._TmpDemFile) == True):
            os.remove(self._TmpDemFile)
        if((os.path.isfile(srtmf_src) == True) and ('_src.tif' in srtmf_src)):
            os.remove(srtmf_src)
            
        os.rename(self._TmpFile, self._TmpDemFile)
        
        if(os.path.isfile(tmpFile) == True):
            os.remove(tmpFile)
        
        self.logger.info('DEM received and prepared')      
        return True


    def gdalDEM_Shade(self):
        head, tmpFile = os.path.split(self._TmpFile)
        head, tmpDemFile = os.path.split(self._TmpDemFile)
        curdir = os.path.curdir
        
        altitude = 90.0 - float32(mean(self.config.solze_arr))
        azimuth = float32(mean(self.config.solaz_arr))
        command = 'gdaldem hillshade '
        options = '-compute_edges -az ' + str(azimuth) + ' -alt ' + str(altitude)
        callstr = command + options + ' ' + tmpDemFile + ' ' + tmpFile + self._DEV0

        l.acquire()        
        try:
            os.chdir(head)
            subprocess.call(callstr, shell=self._SHELL)  
            os.chdir(curdir)
        except:
            self.logger.fatal('shell execution error using gdaldem shade')
            return False
        finally:
            l.release()
        
        self.importBand(self.SDW)
        return True


    def gdalDEM_Slope(self):
        head, tmpFile = os.path.split(self._TmpFile)
        head, tmpDemFile = os.path.split(self._TmpDemFile)
        curdir = os.path.curdir
                
        command = 'gdaldem slope '
        options = '-compute_edges'
        callstr = command + options + ' ' + tmpDemFile + ' ' + tmpFile + self._DEV0

        l.acquire()
        try:
            os.chdir(head)
            subprocess.call(callstr, shell=self._SHELL)
            os.chdir(curdir)
        except:
            self.logger.fatal('shell execution error using gdaldem slope')
            return False
        finally:
            l.release()

        self.importBand(self.SLP)    
        return True


    def gdalDEM_Aspect(self):
        head, tmpFile = os.path.split(self._TmpFile)
        head, tmpDemFile = os.path.split(self._TmpDemFile)
        curdir = os.path.curdir

        command = 'gdaldem aspect '
        options = '-compute_edges'
        callstr = command + options + ' ' + tmpDemFile + ' ' + tmpFile + self._DEV0 

        l.acquire()
        try:
            os.chdir(head)
            subprocess.call(callstr, shell=self._SHELL)   
            os.chdir(curdir)
        except:
            self.logger.fatal('shell execution error using gdaldem aspect')
            return False
        finally:
            l.release()

        self.importBand(self.ASP)
        return True


    def importBand(self, index, filename=None):
        bandName = self.getBandNameFromIndex(index)
        
        if filename == None: # input via GDAL (GeoTiff data, DEM):
            tmpFile = self._TmpFile
            curdir = os.path.curdir
            head, tail = os.path.split(tmpFile)
            l.acquire()
            os.chdir(head)
            while True:
                try:
                    indataArr = gdal.Open(tail, GA_ReadOnly)
                    break
                except:
                    sleep(1)
                    continue

            os.chdir(curdir)
            l.release()           
            tgt_nrows = indataArr.RasterYSize
            tgt_ncols = indataArr.RasterXSize
        else: # the new input for JP2 data:
            warnings.filterwarnings("ignore")
            kwargs = {"tilesize": (2048, 2048), "prog": "RPCL"}
            indataset = glymur.Jp2k(filename, **kwargs)  
            # now the resamling:
            src_nrows = indataset.shape[0]
            src_ncols = indataset.shape[1]
            tgt_nrows = self.config.nrows
            tgt_ncols = self.config.ncols
            fullRes = indataset[:]            
            if src_nrows == tgt_nrows:
            # no resamling required:
                indataArr = fullRes
            elif (src_nrows / tgt_nrows) == 2:
                # mean per 2x2 block of pixels of 10m band for 20m res
                indataArr = uint16(block_reduce(fullRes,block_size=(2,2), func=mean)+0.5)
            elif (src_nrows / tgt_nrows) == 3:
                # mean per 3x3 block of pixels of 20m band for 60m res
                indataArr = uint16(block_reduce(fullRes,block_size=(3,3), func=mean)+0.5)
            elif (src_nrows / tgt_nrows) == 6:
                # mean per 6x6 block of pixels of 10m band for 60m res
                indataArr = uint16(block_reduce(fullRes,block_size=(6,6), func=mean)+0.5)
            elif tgt_nrows > src_nrows:
                # upsampling is required:
                indataArr = indataset[:]
                sizeUp = tgt_nrows
                # order 3 is for cubic spline:
                indataArr = (skit_resize(fullRes.astype(uint16),([sizeUp,sizeUp]),order=3)*65535.).round().astype(uint16)
            indataset = None
                
        if (index < 10) and (indataArr.max() == 0):
            self.logger.fatal('Band ' + bandName + ' does not contain any data')
            return False
            
        h5file = openFile(self._ImageDataBase, mode='a', title =  str(self._resolution) + 'm bands')
        if(h5file.__contains__('/tmp') == False):
            h5file.createGroup('/', 'tmp', 'temporary data')

        if(h5file.__contains__('/metadata')):
            group = h5file.root.metadata
            table = h5file.root.metadata.META
        else:
            group = h5file.createGroup('/', 'metadata', 'metadata information')
            table = h5file.createTable(group, 'META', Particle, "Meta Data")

        if(h5file.__contains__('/arrays')):
            garrays = h5file.root.arrays
        else:
            garrays = h5file.createGroup('/', 'arrays', 'band arrays')

        if(h5file.__contains__('/arrays/' + bandName)):
            self.logger.error( bandName + ' already exists in database ' + self._ImageDataBase)
            h5file.close()
            indataArr = None
            return True
        else:
            filters = Filters(complib="zlib", complevel=1)
            if filename == None:
                inband = indataArr.GetRasterBand(1)
                dtOut = self.setDataType(inband.DataType)
                eArray = h5file.createEArray(garrays, bandName, dtOut, (0, tgt_ncols), bandName, filters=filters)
                for i in range(inband.YSize):
                    scanline = inband.ReadAsArray(0, i, inband.XSize, 1, inband.XSize, 1)
                    scanline = choose( equal( scanline, None), (scanline, None) )
                    eArray.append(scanline)
            else:
                dtOut = self.setDataType(indataArr.dtype)
                eArray = h5file.createEArray(garrays, bandName, dtOut, (0, tgt_ncols), bandName, filters=filters)
                eArray.append(indataArr)

        particle = table.row
        particle['bandName'] = bandName
        #particle['geoTransformation'] = self._geoTransformation
        #particle['projectionRef'] = projectionRef
        particle['rasterYSize'] = tgt_nrows
        particle['rasterXSize'] = tgt_ncols
        particle['rasterCount'] = 1
        particle.append()
        table.flush()
        h5file.close()
        indataArr = None
        self.config.timestamp('L2A_Tables: band ' + bandName + ' imported')
        return True

    def exportBandList(self):
        sourceDir = self._L2A_bandDir
        if(os.path.exists(sourceDir) == False):
            self.logger.fatal('missing directory %s:' % sourceDir)
            return False

        database = self._ImageDataBase

        self.config.timestamp('L2A_Tables: start export')
        if(self._resolution == 10):
            if(self.acMode == True):
                bandIndex = [1,2,3,7,13,17,18]
            else:
                bandIndex = [1,2,3,7,13]
        elif(self._resolution == 20):
            if(self.acMode == True):
                bandIndex = [1,2,3,4,5,6,8,11,12,13,14,15,16,17,18,19,27]
            else:
                bandIndex = [1,2,3,4,5,6,8,11,12,13,14,15,16]
        elif(self._resolution == 60):
            if(self.acMode == True):
                bandIndex = [0,1,2,3,4,5,6,8,9,11,12,13,14,15,16,17,18,19,27]
            else:
                bandIndex = [0,1,2,3,4,5,6,8,9,11,12,13,14,15,16]
        #prepare the xml export
        tileId = self.config.L2A_TILE_ID
        consolidatedTile = False
        if ('L2A_CO_' in tileId):
            consolidatedTile = True
            self.logger.info('consolidated tile ' + tileId + ': no entry in user product metadata generated')
        else:
            Granules = objectify.Element('Granules')
            Granules.attrib['granuleIdentifier'] = self.config.L2A_TILE_ID
            Granules.attrib['datastripIdentifier'] = self.config.L2A_DS_ID
            Granules.attrib['imageFormat'] = 'JPEG2000'
        h5file = openFile(self._ImageDataBase, mode='r')
        for index in bandIndex:
            bandName = self.getBandNameFromIndex(index)
            filename = self._L2A_Tile_BND_File
            filename = filename.replace('BXX', bandName)
            if (bandName == 'VIS'):
                filename = self._L2A_Tile_VIS_File                
            elif (bandName == 'SNW') :
                filename = self._L2A_Tile_SNW_File
            elif(bandName == 'CLD'):
                filename = self._L2A_Tile_CLD_File
            elif(bandName == 'SCL'):
                filename = self._L2A_Tile_SCL_File
            elif(bandName == 'AOT'):
                filename = self._L2A_Tile_AOT_File
            elif(bandName == 'WVP'):
                filename = self._L2A_Tile_WVP_File
            elif(bandName == 'DEM'):
                filename = self._L2A_Tile_DEM_File
                demDir = self.config._demDirectory
                if demDir == 'NONE':
                    continue
            elif(bandName == 'DDV'):
                filename = self._L2A_Tile_DDV_File
                if self.logger.level != logging.DEBUG:
                    continue
            try:
                node = h5file.getNode('/arrays', bandName)
                band = node.read()
                kwargs = {"tilesize": (2048, 2048), "prog": "RPCL"}
                glymur.Jp2k(filename, band.astype(uint16), **kwargs)
                self.config.timestamp('L2A_Tables: band ' + bandName + ' exported')
                filename = os.path.basename(filename.strip('.jp2'))
                if (bandName != 'VIS' and consolidatedTile == False):
                    imageId = etree.Element('IMAGE_ID_2A')
                    imageId.text = filename
                    Granules.append(imageId)
            except:
                continue
        
        h5file.close()
        # update on UP level:
        xp = L2A_XmlParser(self.config, 'UP2A')
        # SIITBX-64: remove unsupported bands 8 and 10:
        try:
            pi = xp.getTree('General_Info', 'L2A_Product_Info')
            bl = pi.Query_Options.Band_List.BAND_NAME
            if(self._resolution == 60):
                for i in range(len(bl)):
                    if bl[i].text == 'B8':
                        if self.checkAotMapIsPresent('10') == False:
                            del bl[i]
                        continue
                    if bl[i].text == 'B10':
                        del bl[i]
                        break
        except:
            self.logger.info('Unsupported band entries already removed')

        # SIITBX-62: set to current resolution:
        try:
            pic = xp.getTree('General_Info', 'L2A_Product_Image_Characteristics')
            si = pic.Spectral_Information_List.Spectral_Information
            if self._resolution == 60:
                for i in range(len(si)):        
                    si[i].RESOLUTION = 60
            elif self._resolution == 20:
                for i in [1,2,3,4,5,6,8,11,12]:        
                    si[i].RESOLUTION = 20           
            elif self._resolution == 10:
                for i in [1,2,3,7]:
                    si[i].RESOLUTION = 10
                # SIITBX-64: add info for Band 8:  
                bn = etree.Element('BAND_NAME')            
                bn.text = 'B8'
                bl = pi.Query_Options.Band_List
                bl.insert(7,bn)
        except:
            self.logger.warning('No Spectral_Information in user metadata available')
            
        if consolidatedTile == False:
            gl = objectify.Element('Granule_List')     
            gl.append(Granules)
            po = pi.L2A_Product_Organisation
            po.append(gl)

        xp.export()

        # update on tile level:
        if(self._resolution > 10):
            xp = L2A_XmlParser(self.config, 'T2A')
            gi = xp.getRoot('General_Info')
            tiOld = xp.getTree('General_Info', 'TILE_ID_2A')
            tiNew = etree.Element('TILE_ID_2A')
            tiNew.text = self.config.L2A_TILE_ID
            gi.replace(tiOld, tiNew)
            dsOld = xp.getTree('General_Info', 'DATASTRIP_ID_2A')
            dsNew = etree.Element('DATASTRIP_ID_2A')
            dsNew.text = self.config.L2A_DS_ID
            gi.replace(dsOld, dsNew)
            pxlqi2a = xp.getTree('Quality_Indicators_Info', 'L2A_Pixel_Level_QI')
            if pxlqi2a == False:
                pxlqi2a = objectify.Element('L2A_Pixel_Level_QI')
                fn = os.path.basename(self._L2A_Tile_CLD_File)
                fn = fn.replace('_60m.jp2', '')
                fn = fn.replace('_20m.jp2', '')
                pxlqi2a.CLOUD_CONFIDENCE_MASK = fn
                fn = os.path.basename(self._L2A_Tile_SNW_File) 
                fn = fn.replace('_60m.jp2', '')            
                fn = fn.replace('_20m.jp2', '')
                pxlqi2a.SNOW_ICE_CONFIDENCE_MASK = fn		
                qii = xp.getRoot('Quality_Indicators_Info')
                qii.insert(3, pxlqi2a)
                pviOld = xp.getTree('Quality_Indicators_Info', 'PVI_FILENAME')
                pviNew = etree.Element('PVI_FILENAME')
                fn = os.path.basename(self._L2A_Tile_PVI_File)
                fn = fn.replace('.jp2', '')  
                pviNew.text = fn
                qii.replace(pviOld, pviNew)
                
            xp.export()

            
        self.createPreviewImage()
        self.config.timestamp('L2A_Tables: preview image exported')   
        # cleanup:
        if(os.path.isfile(self._TmpFile)):
            os.remove(self._TmpFile)
        if(os.path.isfile(database)):
            os.remove(database)
            
        self.config.timestamp('L2A_Tables: stop export')
        return True


    def createPreviewImage(self):
        from PIL import Image
        filename = self._L2A_Tile_PVI_File
        if os.path.exists(filename):
            self.logger.info('Preview Image already exists')
            return True
        else:
            self.logger.info('Creating Preview Image')
               
        acMode = self.acMode
        self.acMode = False
        b = self.getBand(self.B02)
        g = self.getBand(self.B03)
        r = self.getBand(self.B04)
        self.acMode = acMode
        
        b = self.scaleImgArray(b)
        g = self.scaleImgArray(g)
        r = self.scaleImgArray(r)

        b = Image.fromarray(b)
        g = Image.fromarray(g)
        r = Image.fromarray(r)

        try:
            out = Image.merge('RGB', (r,g,b))
            a = array(out)
            kwargs = {"tilesize": (2048, 2048), "prog": "RPCL"}
            glymur.Jp2k(filename, a.astype(uint8), **kwargs)   
            self.logger.debug('preview image exported')
            return True
        except:
            self.logger.fatal('preview image export failed')
            return False


    def scaleImgArray(self, arr):
        if(arr.ndim) != 2:
            self.logger.fatal('must be a 2 dimensional array')
            return False
        if self.logger.level != logging.DEBUG:
            src_ncols = self.config.ncols
            src_nrows = self.config.nrows
            tgt_ncols = 343.0
            tgt_nrows = 343.0
            zoomX = float64(tgt_ncols)/float64(src_ncols)
            zoomY = float64(tgt_nrows)/float64(src_nrows)
            arr = zoom(arr, ([zoomX,zoomY]), order=0)

        arrclip = arr.copy()
        min_ = 0.0
        max_ = 0.250
        scale = 255.0
        arr = clip(arrclip, min_, max_)
        #SIITBX-50: wrong scale was used: 
        scaledArr = uint8(arr*scale/max_)
        return scaledArr


    def testDb(self):
        result = False
        try:
            h5file = openFile(self._ImageDataBase, mode='r')
            h5file.getNode('/arrays', 'B02')
            status = 'Database ' + self._ImageDataBase + ' exists and can be used'
            result = True
        except:
            status = 'Database  ' + self._ImageDataBase + ' will be removed due to corruption'
            os.remove(self._ImageDataBase)
            result = False

        h5file.close()
        self.logger.info(status)
        return result


    def hasBand(self, index):
        result = False
        bandName = self.getBandNameFromIndex(index)
        try:
            h5file = openFile(self._ImageDataBase, mode='r')
            h5file.getNode('/arrays', bandName)
            self.logger.debug('Channel %s is present', self.getBandNameFromIndex(index))
            result = True
        except:
            self.logger.debug('Channel %s is not available', self.getBandNameFromIndex(index))
            result = False

        h5file.close()
        return result


    def getBandSize(self, index):
        h5file = openFile(self._ImageDataBase, mode='r')
        bandName = self.getBandNameFromIndex(index)
        table = h5file.root.metadata.META
        for x in table.iterrows():
            if(x['bandName'] == bandName):
                nrows = x['rasterYSize']
                ncols = x['rasterXSize']
                count = x['rasterCount']
        table.flush()
        h5file.close()
        return(nrows, ncols, count)


    def getBand(self, index, dataType=uint16):
        bandName = self.getBandNameFromIndex(index)
        # the output is context sensitive
        # it will return TOA_reflectance (0:1) if self.acMode = False (this is the scene classification mode)
        # it will return the unmodified value for cirrus and all channels > 12, these are all generated products
        # it will return the radiance if self.acMode = True (this is the atmospheric correction mode)

        h5file = openFile(self._ImageDataBase, mode='r')
        node = h5file.getNode('/arrays', bandName)
        nrows, ncols, count = self.getBandSize(index)
        if (count < 1):
            self.logger.fatal('Insufficient band size: ' + count)
            h5file.close()

        array = node.read()
        h5file.close()

        if(index > 12):
            return array # no further modification
        elif(self.acMode == True):
            L = self.TOA_refl2rad(index, array) # return radiance
            return L
        else: # return reflectance value:
            return (array / float32(self.config.dnScale)) # scaling from 0:1

       
    def TOA_refl2rad(self, index, indataArr):
        ''' Converts the TOA reflectance to radiance.

            :param indataArray: the digital numbers representing TOA reflectance.
            :type indataArray: a 2 dimensional numpy array (row x column) of type unsigned int 16.
            :return: the pixel data converted to radiance.
            :rtype: a 2 dimensional numpy array (row x column) of type unsigned int 16, representing radiance.
            
            Additional inputs from L1 user Product_Image_Characteristics metadata:
            * QUANTIFICATION_VALUE: the scaling factor for converting DN to reflectance.
            * SOLAR_IRRADIANCE: the mean solar exoatmospheric irradiances for each band.
            * d2: the sun-earth distance, calculated from metadata as d2 = 1.0 / U

            Additional inputs from L1 tile Geometric_Info metadata:
            * Sun_Angles_Grid.Zenith.Values: the interpolated zenith angles grid.

        '''
        nrows = self.config.nrows
        ncols = self.config.ncols

        if(self.config.resolution == 10):
            validBand = [None,0,1,2,None,None,None,3,None,None,None,None,None]
        else:
            validBand = [0,1,2,3,4,5,6,None,7,8,9,10,11]
        bandIndex = validBand[index]
        if(bandIndex == None):
            self.logger.fatal('Wrong band index %02d for selected resolution %02d', index, self.config.resolution)

        # The solar zenith array:
        x = arange(nrows, dtype=float32) / (nrows-1) * self.config.solze_arr.shape[0]
        y = arange(ncols, dtype=float32) / (ncols-1) * self.config.solze_arr.shape[1]
        sza = float32(rectBivariateSpline(x,y,self.config.solze_arr))
 
        # Reflectance to ratiance modification according to Rolf, 27/10/2015:
        rho = indataArr.astype(float32)
        c1  = float32(self.config.c1[bandIndex])
        Es  = float32(self.config.e0[bandIndex])
        sc  = float32(1 / (c1 * self.config.dnScale))
        
        # adding again the earth sun correction d2 = 1.0/U:
        pi32_d2 = float32(pi) * self.config.d2

        # The final formaula is:
        # rad = rho * cos(radians(sza)) * Es * sc / (pi * d2)
        # where: d2 = 1.0 / U
        # scale: 1 / (0.001 * 1000) = 1 (default)

        # To reduce the memory consumption for 10m images:
        rad = zeros_like(rho)
        for i in range(nrows):
            rad_sza = float32(radians(sza[i,:]))
            cos_sza = float32(cos(rad_sza))
            rho_cos_sza_Es_sc = float32(rho[i,:] * cos_sza * Es * sc)        
            rad[i,:] = rho_cos_sza_Es_sc / pi32_d2

        return rad


    def getDataType(self, index):
        h5file = openFile(self._ImageDataBase, mode='r')
        bandName = self.getBandNameFromIndex(index)
        node = h5file.getNode('/arrays', bandName)
        dt = node.dtype
        h5file.close()
        return(dt)


    def setBand(self, index, array):
        h5file = openFile(self._ImageDataBase, mode='a')
        bandName = self.getBandNameFromIndex(index)

        if(h5file.__contains__('/arrays/' + bandName)):
            node = h5file.getNode('/arrays', bandName)
            node.remove()

        arr = h5file.root.arrays
        dtIn = self.setDataType(array.dtype)
        filters = Filters(complib="zlib", complevel=1)
        node = h5file.createEArray(arr, bandName, dtIn, (0,array.shape[1]), bandName, filters=filters)
        self.logger.debug('Channel %02d %s added to table', index, self.getBandNameFromIndex(index))
        node.append(array)

        table = h5file.root.metadata.META
        update = False
        # if row exists, change it:
        for row in table.iterrows():
            if(row['bandName'] == bandName):
                row['rasterYSize'] = array.shape[0]
                row['rasterXSize'] = array.shape[1]
                row['rasterCount'] = 1
                row.update()
                update = True
        # else append it:
        if(update == False):
            row = table.row
            row['bandName'] = bandName
            row['rasterYSize'] = array.shape[0]
            row['rasterXSize'] = array.shape[1]
            row['rasterCount'] = 1
            row.append()

        table.flush()
        h5file.close()
        return


    def removeBand(self, index):
        h5file = openFile(self._ImageDataBase, mode='a')
        bandName = self.getBandNameFromIndex(index)
        table = h5file.root.metadata.META
        if(h5file.__contains__('/arrays/' + bandName)):
            node = h5file.getNode('/arrays', bandName)
            node.remove()
            table.flush()
            self.logger.debug('Channel %02d %s removed from table', index, self.getBandNameFromIndex(index))


        h5file.close()
        return


    def removeAllBands(self):
        h5file = openFile(self._ImageDataBase, mode='a')
        if(h5file.__contains__('/arrays/')):
            node = h5file.getNode('/', 'arrays')
            del node
            self.logger.debug('All channels removed from table')

        h5file.close()
        self.removeAllTmpBands()
        return


    def getTmpBand(self, index, dataType=int16):
        h5file = openFile(self._ImageDataBase, mode='r')
        bandName = self.getBandNameFromIndex(index)
        node = h5file.getNode('/tmp', bandName)

        if (node.dtype != dataType):
            h5file.close()
            self.logger.fatal('wrong data type, must be: ' + str(node.dtype))
            
        nrows, ncols, count = self.getBandSize(index)
        if (count < 1):
            h5file.close()
            self.logger.fatal('insufficient band size: ' + count)

        array = node.read()
        h5file.close()
        return array


    def setTmpBand(self, index, array):
        h5file = openFile(self._ImageDataBase, mode='a')
        bandName = self.getBandNameFromIndex(index)
        table = h5file.root.metadata.META

        if(h5file.__contains__('/tmp/' + bandName)):
            node = h5file.getNode('/tmp', bandName)
            node.remove()
            table.flush()
        tmp = h5file.root.tmp
        dtIn = self.setDataType(array.dtype)
        filters = Filters(complib="zlib", complevel=1)
        node = h5file.createEArray(tmp, bandName, dtIn, (0,array.shape[1]), bandName, filters=filters)
        self.logger.debug('Temporary channel ' + str(index) + ' added to table')
        node.append(array)
        table.flush()
        h5file.close()
        return


    def removeTmpBand(self, index):
        h5file = openFile(self._ImageDataBase, mode='a')
        bandName = self.getBandNameFromIndex(index)
        table = h5file.root.metadata.META

        if(h5file.__contains__('/tmp/' + bandName)):
            node = h5file.getNode('/tmp', bandName)
            node.remove()
            table.flush()
            self.logger.debug('Temporary channel ' + str(index) + ' removed from table')

        h5file.close()
        return


    def removeAllTmpBands(self):
        h5file = openFile(self._ImageDataBase, mode='a')
        table = h5file.root.metadata.META
        if(h5file.__contains__('/tmp/' +'')):
            node = h5file.getNode('/tmp')
            del node
            table.flush()
            self.logger.debug('All temporary channels removed from table')
        h5file.close()
        return


    def setDataType(self, dtIn):
        if(dtIn == uint8):
            dtOut = UInt8Atom()
        elif(dtIn == uint16):
            dtOut = UInt16Atom()
        elif(dtIn == int16):
            dtOut = Int16Atom()
        elif(dtIn == uint32):
            dtOut = UInt32Atom()
        elif(dtIn == int32):
            dtOut = Int32Atom()
        elif(dtIn == float32):
            dtOut = Float32Atom()
        elif(dtIn == float64):
            dtOut = Float64Atom()
        elif(dtIn == GDT_Byte):
            dtOut = UInt8Atom()
        elif(dtIn == GDT_UInt16):
            dtOut = UInt16Atom()
        elif(dtIn == GDT_Int16):
            dtOut = Int16Atom()
        elif(dtIn == GDT_UInt32):
            dtOut = UInt32Atom()
        elif(dtIn == GDT_Int32):
            dtOut = Int32Atom()
        elif(dtIn == GDT_Float32):
            dtOut = Float32Atom()
        elif(dtIn == GDT_Float64):
            dtOut = Float64Atom()

        return dtOut


    def getArray(self, filename):
        filename = self._testdir + filename + '.npy'
        if((os.path.isfile(filename)) == False):
            self.logger.critical('File ' + filename + ' not present')
            return False

        return load(filename)



    def getDiffFromArrays(self, filename1, filename2):
        filename1 = self._testdir + filename1 + '.npy'
        filename2 = self._testdir + filename2 + '.npy'
        if((os.path.isfile(filename1)) == False):
            self.logger.critical('File ' + filename1 + ' not present')
            return False

        if((os.path.isfile(filename2)) == False):
            self.logger.critical('File ' + filename2 + ' not present')
            return False

        arr1 = load(filename1)
        arr2 = load(filename2)
        return (arr1-arr2)


    def saveArray(self, filename, arr):
        filename = self._testdir + filename + '.npy'
        save(filename, arr)

        if(os.path.exists(self._L2A_bandDir) == False):
            os.makedirs(self._L2A_bandDir)
            self.logger.info('File ' + filename + ' saved to disk')
        return


    def sceneCouldHaveSnow(self):
        from PIL import Image
        globalSnowMapFn = self.config.snowMapReference
        globalSnowMapFn = os.path.join(self.config.libDir, globalSnowMapFn)
        if((os.path.isfile(globalSnowMapFn)) == False):
            self.logger.error('global snow map not present, snow detection will be performed')
            return True
        
        l.acquire()
        try:
            img = Image.open(globalSnowMapFn)
        finally:
            l.release()
            
        globalSnowMap = array(img)
        xy = self.cornerCoordinates
        xp = L2A_XmlParser(self.config, 'T2A')
        tg = xp.getTree('Geometric_Info', 'Tile_Geocoding')
        hcsName = tg.HORIZONTAL_CS_NAME.text
        zone = hcsName.split()[4]
        zone1 = int(zone[:-1])
        zone2 = zone[-1:].upper()
        lonMin, latMin, dummy = self.transform_utm_to_wgs84(xy[1,0], xy[1,1], zone1, zone2)
        lonMax, latMax, dummy = self.transform_utm_to_wgs84(xy[3,0], xy[3,1], zone1, zone2)

        # Snow map should have a dimension of 7200 x 3600, 20 pixels per degree:
        xMin = int((lonMin + 180.0) * 20.0 + 0.5)
        xMax = int((lonMax + 180.0) * 20.0 + 0.5)
        yMin = 3600 - int((latMax + 90.0) * 20.0 + 0.5) # Inverted by intention
        yMax = 3600 - int((latMin + 90.0) * 20.0 + 0.5) # Inverted by intention
        aoi = globalSnowMap[yMin:yMax,xMin:xMax]
        if(aoi.max() > 0):
            return True
        
        return False
