#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

from numpy import *
from tables import *
import sys
import os, logging
from time import time
from datetime import datetime

from L2A_Manifest import L2A_Manifest
from L2A_Tables import L2A_Tables
from L2A_SceneClass import L2A_SceneClass
from L2A_AtmCorr import L2A_AtmCorr
from L2A_XmlParser import L2A_XmlParser
from L2A_Library import stdoutWrite, stderrWrite
from L2A_Logger import SubProcessLogHandler, getLevel

import multiprocessing
import cPickle as pickle
SUCCESS = 0
FAILURE = 1
l = multiprocessing.Lock()
formatter = logging.Formatter('<check>\n<inspection execution=\"%(asctime)s\" level=\"%(levelname)s\" module=\"%(module)s\" function=\"%(funcName)s\" line=\"%(lineno)d\"/>\n<message contentType=\"Text\">%(message)s</message>\n</check>')
DEFAULT_LEVEL = logging.INFO


class L2A_ProcessTile(multiprocessing.Process):
    def __init__(self, tileId, res, msgQueue, queue):
        multiprocessing.Process.__init__(self)
        self.queue = queue
        self._msgQueue = msgQueue
        self._tileId = tileId
        self._finalRes = res
        self._scOnly = False
        self._configFn = None
        self._logger = None
        self._localTimestamp = time()
        self._res = res
        picFn = os.path.join(tileId,'configPic.p')
        try:
            f = open(picFn, 'rb')
            self.config = pickle.load(f)
            f.close()
        except:
            stderrWrite('Cannot load configuration\n.' % picFn)
 
        self.config._timestamp = datetime.now()
        self.scOnly = self.config.scOnly
        self.config.resolution = self._res


    def get_logger(self):
        return self._logger


    def set_logger(self, value):
        self._logger = value


    def del_logger(self):
        del self._logger

          
    def get_sc_only(self):
        return self._scOnly
 
 
    def set_sc_only(self, value):
        self._scOnly = value
 
 
    def del_sc_only(self):
        del self._scOnly
 
 
    def get_tables(self):
        return self._tables
 
 
    def set_tables(self, value):
        self._tables = value
 
 
    def del_tables(self):
        del self._tables


    def get_config(self):
        return self._config
 
 
    def set_config(self, value):
        self._config = value
 
 
    def del_config(self):
        del self._config
 
 
    def __exit__(self):
            sys.exit(-1)
 
    config = property(get_config, set_config, del_config, "config's docstring")
    tables = property(get_tables, set_tables, del_tables, "tables's docstring")
    scOnly = property(get_sc_only, set_sc_only, del_sc_only, "scOnly's docstring")
    logger = property(get_logger, set_logger, del_logger, "logger's docstring")

    def run(self):
        if self.config.resolution == 0:
            stdoutWrite('no resolution specified, will process all resolutions.\n')            
        else:
            stdoutWrite('selected resolution is %s m.\n' % self.config.resolution)
        
        self.setupLogger()
        logger = self.logger
        logger.level = getLevel(self.config.logLevel)
        self.config.logger = logger

        if(self.config.resolution == 0):
            if self.process_60() == False:
                return False      
            if self.process_20() == False:
                return False
            if self.process_10() == False:
                return False      
            return True
        elif(self.config.resolution == 60):
            if self.process_60() == False:
                return False
        elif(self.config.resolution == 20):
            if self.process_20() == False:
                return False
        elif(self.config.resolution == 10):
            if (self.scOnly == True):
                stdoutWrite('No scene classification processing for 10 m resolution in sc_only mode.\n')
                self.logger.info('no scene classification processing for 10 m resolution in sc_only mode')   
                tMeasure = time() - self.config.tStart
                self.config.writeTimeEstimation(tMeasure)           
                self._msgQueue.put(FAILURE)
                return False
            if self.process_20() == False:
                return False
            elif self.process_10() == False:
                return False
        return True
        

    def process_60(self):
        self.config.resolution = 60
        self.config.readPreferences()
        self.config.createOrUpdateL2A_UserProduct()
        self.tables = L2A_Tables(self.config)
        return self.process()


    def process_20(self):
        self.config.resolution = 20   
        self.config.readPreferences()
        self.config.createOrUpdateL2A_UserProduct()
        self.tables = L2A_Tables(self.config)                  
        return self.process()


    def process_10(self):
        if (self.scOnly == True):
            stdoutWrite('No scene classification processing for 10 m resolution in sc_only mode.\n')
            self.logger.info('no scene classification processing for 10 m resolution in sc_only mode')   
            tMeasure = time() - self.config.tStart
            self.config.writeTimeEstimation(tMeasure)           
            self._msgQueue.put(FAILURE)
            return False        
        
        res20 = 20
        if self.tables.checkAotMapIsPresent(res20) == False:
            stdoutWrite('20 m resolution must be processed first.\n')
            self.logger.info('20 m resolution must be processed first')
            if self.process_20() == False:
                return False
        
        self.config.resolution = 10  
        self.config.readPreferences()  
        self.config.createOrUpdateL2A_UserProduct()
        self.tables = L2A_Tables(self.config)            
        return self.process()


    def process(self):     
        p = multiprocessing.current_process()
        self.config.readTileMetadata()
        if self.tables.checkAotMapIsPresent(self.config.resolution):
            self.config.timestamp('L2A_ProcessTile: resolution '+ str(self.config.resolution) + 'm already processed')
            self._msgQueue.put(SUCCESS)
            return True
        
        astr = 'L2A_ProcessTile: processing with resolution ' + str(self.config.resolution) + ' m'
        self.config.timestamp(astr)
        self.config.timestamp('L2A_ProcessTile: start of pre processing')
        if(self.preprocess() == False):
            self.logger.fatal('Module %s - %s failed' %(p.name, self.config.L2A_TILE_ID))
            self._msgQueue.put(FAILURE)
            return False
     
        if(self.config.resolution > 10):
            self.config.timestamp('L2A_ProcessTile: start of Scene Classification')
            sc = L2A_SceneClass(self.config, self.tables)
            self.logger.info('performing Scene Classification with resolution %d m' % self.config.resolution)
            if(sc.process() == False):
                self.logger.fatal('Module %s - %s failed' %(p.name, self.config.L2A_TILE_ID))
                self._msgQueue.put(FAILURE)
                return False
           
        if(self.scOnly == False):
            ac = L2A_AtmCorr(self.config, self.tables)
            if(self.config.resolution > 10):
                if ((self.config.midLatitude == 'AUTO') or (self.config.aerosolType == 'AUTO')):
                    self.config.timestamp('L2A_ProcessTile: start of Automatic Aerosol Type Detection')
                    self.logger.info('performing aerosol type detection with resolution %d m' % self.config.resolution)
                    if(ac.automaticAerosolDetection() == False):
                        self.logger.fatal('Module %s - %s failed' %(p.name, self.config.L2A_TILE_ID))
                        self._msgQueue.put(FAILURE)   
                        return False                
 
            self.config.timestamp('L2A_ProcessTile: start of Atmospheric Correction')
            self.logger.info('performing Atmospheric Correction with resolution %d m' % self.config.resolution)
            ac.aerosolDetection = False
            if(ac.process() == False):
                self.logger.fatal('Module %s - %s failed' %(p.name, self.config.L2A_TILE_ID))
                self._msgQueue.put(FAILURE)
                return False
          
        self.config.timestamp('L2A_ProcessTile: start of post processing')
        if(self.postprocess() == False):
            self.logger.fatal('Module %s - %s failed' %(p.name, self.config.L2A_TILE_ID))
            self._msgQueue.put(FAILURE)
            return False
                 
        self._msgQueue.put(SUCCESS)
        return True
 
 
    def setupLogger(self):
        # create the logger to use.
        logger = logging.getLogger('sen2cor.subprocess')


        for handler in logger.handlers:
            assert not isinstance(handler, SubProcessLogHandler)
            logger.removeHandler(handler)
        # add the handler
        handler = SubProcessLogHandler(self.queue)
        handler.setFormatter(formatter)
        logger.addHandler(handler)

        logger.setLevel(DEFAULT_LEVEL)
        self._logger = logger
 
         
    def preprocess(self):
        self.logger.info('pre-processing with resolution %d m', self.config.resolution)
        # this is to check the config for the L2A_AtmCorr in ahead.
        # This has historical reasons due to ATCOR porting.
        # Should be moved to the L2A_Config for better design:
        dummy = L2A_AtmCorr(self.config, self.logger)
        dummy.checkConfiguration()
 
        # validate the meta data:
        xp = L2A_XmlParser(self.config, 'T1C')
        xp.export()
        xp.validate()
 
        if(self.tables.checkBandCount() == False):
            self.logger.fatal('insufficient nr. of bands in tile: ' + self.config.L2A_TILE_ID)
            return False

        if(self.config.ozoneContent == '0'):
            try:
                ozoneMeasured = self.tables.getAuxData(self.tables.OZO)
                self.config.ozoneMean = ozoneMeasured.mean()
                self.logger.info('ozone mean value is: ' + str(self.config.ozoneMean))
            except:
                if self.config.midLatitude == 'SUMMER':
                    self.logger.info('no ozone data present, standard mid summer will be used')
                    self.config.ozoneContent = 'h'
                elif self.config.midLatitude == 'WINTER':
                    self.logger.info('no ozone data present, standard mid winter will be used')
                    self.config.ozoneContent = 'w'
                else:
                    self.logger.info('no ozone data present and no mid latitude configured, standard mid summer will be used')
                    self.config.ozoneContent = 'h'

        if(self.tables.importBandList() == False):
            self.logger.fatal('import of band list failed')
            return False

        if((self.config.aerosolType != 'AUTO') and (self.config.midLatitude != 'AUTO')):
            self.config.createAtmDataFilename()

        return True

 
    def postprocess(self):
        self.logger.info('post-processing with resolution %d m', self.config.resolution)
         
        res = True 
        if self.tables.exportBandList() == False:
            res = False
        # validate the meta data:
        xp = L2A_XmlParser(self.config, 'T2A')
        xp.export()
        xp.validate()
        if self.config.postprocess() == False:
            res = False
        
        if multiprocessing.active_children() == False:
            #Create the manifest.safe (L2A)
            mn = L2A_Manifest(self.config)
            mn.generate(self.config._L2A_UP_DIR, self.config.L2A_MANIFEST_SAFE)
            xp = L2A_XmlParser(self.config, 'Manifest')
            xp.export()
            xp.validate()
 
        tMeasure = time() - self._localTimestamp
        self.config.writeTimeEstimation(tMeasure)
        return res

