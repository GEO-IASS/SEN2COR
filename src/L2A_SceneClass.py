#!/usr/bin/env pythonimport osimport sysfrom numpy import *from scipy.ndimage.morphology import *from scipy.ndimage.interpolation import *from scipy.ndimage.filters import median_filterfrom time import timefrom lxml import etree, objectifyfrom L2A_Borg import Borgfrom L2A_Config import L2A_Configfrom L2A_Library import *from L2A_Tables import L2A_Tablesfrom L2A_XmlParser import L2A_XmlParserfrom PIL import Imageset_printoptions(precision = 7, suppress = True)class L2A_SceneClass(Borg):    def __init__(self, config, tables):        self._notClassified = 100        self._notSnow = 50        self._config = config        self._tables = tables        self.tables.acMode = False        x,y,n = tables.getBandSize(self.tables.B02)        self.classificationMask = ones([x,y], uint8) * self._notClassified        self.confidenceMaskSnow = zeros_like(tables.getBand(self.tables.B02))        self.confidenceMaskCloud = zeros_like(tables.getBand(self.tables.B02))        self._meanShadowDistance = 0        self.filter =  None        self.LOWEST = 0.000001        self._noData = self.config.getInt('Scene_Classification/Classificators', 'NO_DATA')        self._saturatedDefective = self.config.getInt('Scene_Classification/Classificators', 'SATURATED_DEFECTIVE')        self._darkFeatures = self.config.getInt('Scene_Classification/Classificators', 'DARK_FEATURES')        self._cloudShadows = self.config.getInt('Scene_Classification/Classificators', 'CLOUD_SHADOWS')        self._vegetation = self.config.getInt('Scene_Classification/Classificators', 'VEGETATION')        self._bareSoils = self.config.getInt('Scene_Classification/Classificators', 'BARE_SOILS')        self._water = self.config.getInt('Scene_Classification/Classificators', 'WATER')        self._lowProbaClouds = self.config.getInt('Scene_Classification/Classificators', 'LOW_PROBA_CLOUDS')        self._medProbaClouds = self.config.getInt('Scene_Classification/Classificators', 'MEDIUM_PROBA_CLOUDS')        self._highProbaClouds = self.config.getInt('Scene_Classification/Classificators', 'HIGH_PROBA_CLOUDS')        self._thinCirrus = self.config.getInt('Scene_Classification/Classificators', 'THIN_CIRRUS')        self._snowIce = self.config.getInt('Scene_Classification/Classificators', 'SNOW_ICE')        self.config.logger.debug('Module L2A_SceneClass initialized')        self._processingStatus = True        self._sumPercentage = 0.0    def assignClassifcation(self, arr, treshold, classification):        cm = self.classificationMask        cm[(arr == treshold) & (cm == self._notClassified)] = classification        self.confidenceMaskCloud[(cm == classification)] = 0        return    def get_config(self):        return self._config    def get_tables(self):        return self._tables    def set_config(self, value):        self._config = value    def set_tables(self, value):        self._tables = value    def del_config(self):        del self._config    def del_tables(self):        del self._tables    config = property(get_config, set_config, del_config, "config's docstring")    tables = property(get_tables, set_tables, del_tables, "tables's docstring")    def preprocess(self):        B03 = self.tables.getBand(self.tables.B03)        B8A = self.tables.getBand(self.tables.B8A)        self.classificationMask[(B03==0) & (B8A==0)] = self._noData        return    def postprocess(self):        if(self._processingStatus == False):            return False        CM = self.classificationMask        CM[(CM == self._notClassified)] = self._saturatedDefective             value = self.config.getInt('Scene_Classification/Calibration', 'Median_Filter')        if(value > 0):            CM = median_filter(CM, value)            self.config.logger.info('Filtering output with level: ' + str(value))                self.config.logger.info('Storing final Classification Mask')        self.tables.setBand(self.tables.SCL,(CM).astype(uint8))        self.config.logger.info('Storing final Snow Confidence Mask')        self.tables.setBand(self.tables.SNW,(self.confidenceMaskSnow*100+0.5).astype(uint8))        self.config.logger.info('Storing final Cloud Confidence Mask')        self.tables.setBand(self.tables.CLD,(self.confidenceMaskCloud*100+0.5).astype(uint8))                # add L2A quality info on tile level:        self.updateQualityIndicators(1, 'T2A')        # add L2A quality info on user level:        xp = L2A_XmlParser(self.config, 'DS2A')        ti = xp.getTree('Image_Data_Info', 'Tiles_Information')        nrTilesProcessed = len(ti.Tile_List.Tile)        self.updateQualityIndicators(nrTilesProcessed, 'UP2A')    def __exit__(self):        sys.exit(-1)    def __del__(self):        self.config.logger.info('Module L2A_SceneClass deleted')    def L2A_CSND_1_1(self):        # Step 1a: Brightness threshold on red (Band 4)        T1_B04 = self.config.getFloat('Scene_Classification/Thresholds', 'T1_B04')        T2_B04 = self.config.getFloat('Scene_Classification/Thresholds', 'T2_B04')        T1_B08 = 0.04        T2_B08 = 0.15        B04 = self.tables.getBand(self.tables.B04)        B08 = self.tables.getBand(self.tables.B8A)        self.confidenceMaskCloud = clip(B04, T1_B04, T2_B04)        #self.confidenceMaskCloud = ((self.confidenceMaskCloud - T1_B04)/(T2_B04-T1_B04))**2        self.confidenceMaskCloud = ((self.confidenceMaskCloud - T1_B04)/(T2_B04-T1_B04))        CM = self.classificationMask        CM[(B04<T1_B04) & (B08>T1_B08) & (B08<T2_B08) & (CM==self._notClassified)] = self._darkFeatures        self.confidenceMaskCloud[(CM == self._darkFeatures)] = 0        self.config.tracer.debug(statistics(self.confidenceMaskCloud, 'CM Cloud step 1.1'))        return    def L2A_CSND_1_2(self):        # Step 1b: Normalized Difference Snow Index (NDSI)        T1_NDSI_CLD = self.config.getFloat('Scene_Classification/Thresholds', 'T1_NDSI_CLD')        T2_NDSI_CLD = self.config.getFloat('Scene_Classification/Thresholds', 'T2_NDSI_CLD')        f1 = self.confidenceMaskCloud > 0        B03 = self.tables.getBand(self.tables.B03)        B11 = self.tables.getBand(self.tables.B11)        NDSI = (B03 - B11) / maximum((B03 + B11), self.LOWEST)        CMC = clip(NDSI, T1_NDSI_CLD, T2_NDSI_CLD)        CMC = ((CMC - T1_NDSI_CLD)/(T2_NDSI_CLD-T1_NDSI_CLD))        CM = self.classificationMask        CM[(CMC==0)] = self._notClassified        self.confidenceMaskCloud *= CMC        self.config.tracer.debug(statistics(self.confidenceMaskCloud, 'CM Cloud step 1.2'))        return    def L2A_CSND_2_0(self):        return    def L2A_CSND_2_1(self):        # Snow filter 1: Normalized Difference Snow Index (NDSI)        T1_NDSI_SNW = self.config.getFloat('Scene_Classification/Thresholds', 'T1_NDSI_SNW')        T2_NDSI_SNW = self.config.getFloat('Scene_Classification/Thresholds', 'T2_NDSI_SNW')        B03 = self.tables.getBand(self.tables.B03)        B11 = self.tables.getBand(self.tables.B11)        NDSI = (B03 - B11) / maximum((B03 + B11), self.LOWEST)        CMS = clip(NDSI, T1_NDSI_SNW, T2_NDSI_SNW)        CMS = ((CMS - T1_NDSI_SNW)/(T2_NDSI_SNW-T1_NDSI_SNW))        CM = self.classificationMask        CM[(CMS == 0) & (CM == self._notClassified)] = self._notSnow        self.confidenceMaskSnow = CMS        return    def L2A_CSND_2_2(self):        # Snow filter 2: Band 8 thresholds        T1_B8A = self.config.getFloat('Scene_Classification/Thresholds', 'T1_B8A')        T2_B8A = self.config.getFloat('Scene_Classification/Thresholds', 'T2_B8A')        B8A = self.tables.getBand(self.tables.B8A)        CMS = clip(B8A, T1_B8A, T2_B8A)        CMS = ((CMS - T1_B8A) / (T2_B8A - T1_B8A))        CM = self.classificationMask        CM[(CMS == 0) & (CM == self._notClassified)] = self._notSnow        self.confidenceMaskSnow *= CMS        return    def L2A_CSND_2_3(self):        # Snow filter 3: Band 2 thresholds        T1_B02 = self.config.getFloat('Scene_Classification/Thresholds', 'T1_B02')        T2_B02 = self.config.getFloat('Scene_Classification/Thresholds', 'T2_B02')        B02 = self.tables.getBand(self.tables.B02)        CMS = clip(B02, T1_B02, T2_B02)        CMS = ((CMS - T1_B02) / (T2_B02 - T1_B02))        CM = self.classificationMask        CM[(CMS == 0) & (CM == self._notClassified)] = self._notSnow        self.confidenceMaskSnow *= CMS        return    def L2A_CSND_2_4(self):        # Snow filter 4: Ratio Band 2 / Band 4        T1_R_B02_B04 = self.config.getFloat('Scene_Classification/Thresholds', 'T1_R_B02_B04')        T2_R_B02_B04 = self.config.getFloat('Scene_Classification/Thresholds', 'T2_R_B02_B04')        B02 = self.tables.getBand(self.tables.B02)        B04 = self.tables.getBand(self.tables.B04)        RB02_B04 = B02 / maximum(B04,self.LOWEST)        CMS = clip(RB02_B04, T1_R_B02_B04, T2_R_B02_B04)        CMS = ((CMS - T1_R_B02_B04) / (T2_R_B02_B04 - T1_R_B02_B04))        CM = self.classificationMask        CM[(CMS == 0) & (CM == self._notClassified)] = self._notSnow        self.confidenceMaskSnow *= CMS        CM = self.classificationMask        return    def L2A_CSND_2_5(self):        # Snow filter 5: snow boundary zones        T1_SNOW = self.config.getFloat('Scene_Classification/Thresholds', 'T1_SNOW')        T2_SNOW = self.config.getFloat('Scene_Classification/Thresholds', 'T2_SNOW')        B12 = self.tables.getBand(self.tables.B12)        CM = self.classificationMask        CMS = self.confidenceMaskSnow        CMS[B12 > T2_SNOW] = 0        CM[(B12 > T2_SNOW) & (CM == self._notClassified)] = self._notSnow        CM[CM == self._notClassified] = self._snowIce        # important, if classified as snow, this should not become cloud:        self.confidenceMaskCloud[CM == self._snowIce] = 0        # release the lock for the non snow classification        CM[CM == self._notSnow] = self._notClassified        return    def L2A_CSND_3(self):        # Step 3: Normalized Difference Vegetation Index (NDVI)        T1_NDVI = self.config.getFloat('Scene_Classification/Thresholds', 'T1_NDVI')        T2_NDVI = self.config.getFloat('Scene_Classification/Thresholds', 'T2_NDVI')        T1_B2T = 0.15        B02 = self.tables.getBand(self.tables.B02)        B04 = self.tables.getBand(self.tables.B04)        B8A = self.tables.getBand(self.tables.B8A)        NDVI = (B8A - B04) / maximum((B8A + B04), self.LOWEST)        CMC = clip(NDVI, T1_NDVI, T2_NDVI)        CMC = ((CMC - T1_NDVI)/(T2_NDVI-T1_NDVI))        CM = self.classificationMask        CM[(CMC==1) & (CM == self._notClassified) & (B02 < T1_B2T)] = self._vegetation        CMC[(CM== self._vegetation)] = 0        FLT = [(CMC>0) & (CMC < 1.0)]        CMC[FLT] = CMC[FLT] * -1 + 1        self.confidenceMaskCloud[FLT] *= CMC[FLT]        self.config.tracer.debug(statistics(self.confidenceMaskCloud, 'CM Cloud step 3'))        return    def L2A_CSND_4(self):        # Step 4: Ratio Band 8 / Band 3 for senescing vegetation        T1_R_B8A_B03 = self.config.getFloat('Scene_Classification/Thresholds', 'T1_R_B8A_B03')        T2_R_B8A_B03 = self.config.getFloat('Scene_Classification/Thresholds', 'T2_R_B8A_B03')        B03 = self.tables.getBand(self.tables.B03)        B8A = self.tables.getBand(self.tables.B8A)        rb8b3 = B8A/maximum(B03,self.LOWEST)        CMC = clip(rb8b3, T1_R_B8A_B03, T2_R_B8A_B03)        CMC = (CMC - T1_R_B8A_B03) / (T2_R_B8A_B03 - T1_R_B8A_B03)        CM = self.classificationMask        CM[(CMC==1) & (CM == self._notClassified)] = self._vegetation        CMC[(CM== self._vegetation)] = 0        FLT = [(CMC>0) & (CMC < 1.0)]        CMC[FLT] = CMC[FLT] * -1 + 1        self.confidenceMaskCloud[FLT] *= CMC[FLT]        self.config.tracer.debug(statistics(self.confidenceMaskCloud, 'CM Cloud step 4'))        return    def L2A_CSND_5_1(self):        # Step 5.1: Ratio Band 2 / Band 11 for soils        T11_B02 = self.config.getFloat('Scene_Classification/Thresholds', 'T11_B02') # -0.40        T12_B02 = self.config.getFloat('Scene_Classification/Thresholds', 'T12_B02') #  0.46        T11_R_B02_B11 = self.config.getFloat('Scene_Classification/Thresholds', 'T11_R_B02_B11') # 0.55        T12_R_B02_B11 = self.config.getFloat('Scene_Classification/Thresholds', 'T12_R_B02_B11') # 0.80        B02 = self.tables.getBand(self.tables.B02)        B11 = self.tables.getBand(self.tables.B11)        R_B02_B11 = clip((B02/maximum(B11,self.LOWEST)),0,100)        B02_FT = clip(R_B02_B11*T11_B02+T12_B02, 0.15, 0.32)        R_B02_B11_GT_T12_R_B02_B11 = where((R_B02_B11 > T12_R_B02_B11) | (B02 > B02_FT), True, False)        CM = self.classificationMask        CM[(R_B02_B11_GT_T12_R_B02_B11 == False) & (CM == self._notClassified)] = self._bareSoils        self.confidenceMaskCloud[CM == self._bareSoils] = 0        R_B02_B11_GT_T11_R_B02_B11_LE_T12_R_B02_B11 = where((R_B02_B11 > T11_R_B02_B11) & (R_B02_B11 < T12_R_B02_B11), True, False)        a = 1 / (T12_R_B02_B11 - T11_R_B02_B11)        b = -T11_R_B02_B11 * a        CMC = self.confidenceMaskCloud        FLT = (R_B02_B11_GT_T11_R_B02_B11_LE_T12_R_B02_B11 == True) & (R_B02_B11_GT_T12_R_B02_B11 == False) & (CM == self._notClassified)        CMC[FLT] = a * R_B02_B11[FLT] + b        self.confidenceMaskCloud[FLT] *= CMC[FLT]        self.config.tracer.debug(statistics(self.confidenceMaskCloud, 'CM Cloud step 5.1'))        return    def L2A_CSND_5_2(self):        # Step 5.2: Ratio Band 2 / Band 11 for water bodies, dependent on Band 12        T21_B12 = self.config.getFloat('Scene_Classification/Thresholds', 'T21_B12') # 0.1        T22_B12 = self.config.getFloat('Scene_Classification/Thresholds', 'T22_B12') # -0.09        T21_R_B02_B11 = self.config.getFloat('Scene_Classification/Thresholds', 'T21_R_B02_B11') # 2.0        T22_R_B02_B11 = self.config.getFloat('Scene_Classification/Thresholds', 'T22_R_B02_B11') # 4.0        B02 = self.tables.getBand(self.tables.B02) # for Istanbul, add 0.5, else no water discrimination !!!        B11 = self.tables.getBand(self.tables.B11)        B12 = self.tables.getBand(self.tables.B12)        R_B02_B11 = B02 / maximum(B11,self.LOWEST)        B12_FT = clip(R_B02_B11*T21_B12+T22_B12, 0.07, 0.21)        R_B02_B11_GT_T22_R_B02_B11 = where((R_B02_B11 > T22_R_B02_B11) & (B12 < B12_FT), True, False)        CM = self.classificationMask # this is a reference, no need to reassign        CM[(R_B02_B11_GT_T22_R_B02_B11 == True) & (CM == self._notClassified)] = self._water        self.confidenceMaskCloud[CM == self._water] = 0        R15_AMB = (R_B02_B11 < T22_R_B02_B11) & (R_B02_B11 >= T21_R_B02_B11) & (B12 < B12_FT)        if(R15_AMB.size > 0):            a = -1 / (T22_R_B02_B11 - T21_R_B02_B11)            b = -T21_R_B02_B11 * a + 1            CMC = a * R_B02_B11[R15_AMB] + b            self.confidenceMaskCloud[R15_AMB] *= CMC        # second part, modification for improvement of water classification:        T_24 = 0.034        B04 = self.tables.getBand(self.tables.B04)        DIFF24_AMB = B02-B04        #CM = self.classificationMask        F1 = DIFF24_AMB > T_24        CM[F1 & (CM == self._notClassified)] = self._water        self.confidenceMaskCloud[F1 & (CM == self._water)] = 0        self.config.tracer.debug(statistics(self.confidenceMaskCloud, 'CM Cloud step 5.2'))        return    def L2A_CSND_6(self):        # Step 6: Ratio Band 8 / Band 11 for rocks and sands in deserts        T1_R_B8A_B11 = self.config.getFloat('Scene_Classification/Thresholds', 'T1_R_B8A_B11')        T2_R_B8A_B11 = self.config.getFloat('Scene_Classification/Thresholds', 'T2_R_B8A_B11')        B8A = self.tables.getBand(self.tables.B8A)        B11 = self.tables.getBand(self.tables.B11)        R_B8A_B11 = B8A/maximum(B11,self.LOWEST)        CMC = clip(R_B8A_B11, T1_R_B8A_B11, T2_R_B8A_B11)        CMC = (CMC - T1_R_B8A_B11) / (T2_R_B8A_B11 - T1_R_B8A_B11)        self.assignClassifcation(CMC, 0, self._bareSoils)        self.confidenceMaskCloud *= CMC        self.config.tracer.debug(statistics(self.confidenceMaskCloud, 'CM Cloud step 6'))        return    def L2A_CSND_7(self):        T_CLOUD_LP = self.config.getFloat('Scene_Classification/Thresholds', 'T_CLOUD_LP')        T_CLOUD_MP = self.config.getFloat('Scene_Classification/Thresholds', 'T_CLOUD_MP')        T_CLOUD_HP = self.config.getFloat('Scene_Classification/Thresholds', 'T_CLOUD_HP')        T1_B10 = self.config.getFloat('Scene_Classification/Thresholds', 'T1_B10')        T2_B10 = self.config.getFloat('Scene_Classification/Thresholds', 'T2_B10')        B02 = self.tables.getBand(self.tables.B02)        B10 = self.tables.getBand(self.tables.B10)        LPC = self._lowProbaClouds        MPC = self._medProbaClouds        HPC = self._highProbaClouds        CIR = self._thinCirrus        CM = self.classificationMask        CMC = self.confidenceMaskCloud        CM[(CMC > T_CLOUD_LP) & (CMC < T_CLOUD_MP) & (CM == self._notClassified)] = LPC        self.config.tracer.debug(statistics(CMC[(CM == LPC)], 'CM LOW_PROBA_CLOUDS'))        CM[(CMC >= T_CLOUD_MP) & (CMC < T_CLOUD_HP) & (CM == self._notClassified)] = MPC        self.config.tracer.debug(statistics(CMC[(CM == MPC)], 'CM MEDIUM_PROBA_CLOUDS'))        CM[(CMC >= T_CLOUD_HP) & (CM == self._notClassified)] = HPC        self.config.tracer.debug(statistics(CMC[(CM == HPC)], 'CM HIGH_PROBA_CLOUDS'))        CM[(B10 > T1_B10) & (B10 < T2_B10) & (CMC < 0.8) & (B02 < 0.50) & (CM != HPC)] = CIR        self.config.tracer.debug(statistics(CMC[(CM == CIR)], 'CM THIN_CIRRUS'))        CM[(B10 >= T2_B10) & (CM == self._notClassified)]= MPC        self.config.tracer.debug(statistics(CMC[(CM == MPC)], 'CM MEDIUM_PROBA_CLOUDS, step2'))        return    def L2A_SHD(self):        csd1 = self.L2A_CSHD_1()        csd2 = self.L2A_CSHD_2()        CSP = (csd1 * csd2 > 0)        CM = self.classificationMask        if(self.tables.hasBand(self.tables.SDW)):            T_SDW = self.config.getFloat('Scene_Classification/Thresholds', 'T_SDW')            shadow = self.tables.getBand(self.tables.SDW, uint8)            tShadow = array(shadow, float32) / 255.0            CM[(CM == self._darkFeatures) & (tShadow > T_SDW) & (CSP == True)] = self._cloudShadows            CM[(CM == self._water) & (CSP == True)] = self._cloudShadows        else:            CM[(CM == self._darkFeatures) & (CSP == True)] = self._cloudShadows            CM[(CM == self._water) & (CSP == True)] = self._cloudShadows        return    def L2A_CSHD_2(self):        # Part2: radiometric input:        x,y,n = self.tables.getBandSize(2)        BX = zeros((6,x,y), float32)        BX[0,:,:] = self.tables.getBand(self.tables.B02)        BX[1,:,:] = self.tables.getBand(self.tables.B03)        BX[2,:,:] = self.tables.getBand(self.tables.B04)        BX[3,:,:] = self.tables.getBand(self.tables.B8A)        BX[4,:,:] = self.tables.getBand(self.tables.B11)        BX[5,:,:] = self.tables.getBand(self.tables.B12)        RV_MEAN = array([0.0696000, 0.0526667, 0.0537708, 0.0752000, 0.0545000, 0.0255000], dtype=float32)        distance = zeros((6,x,y), float32)        for i in range(0,6):            distance[i,:,:] = abs(BX[i,:,:] - RV_MEAN[i])        T_B02_B12 = self.config.getFloat('Scene_Classification/Thresholds', 'T_B02_B12')        msd = mean(distance, axis=0)        msd = median_filter(msd, 3)        msd = 1.0 - msd        T0 = 1.0 - T_B02_B12        msd[msd < T0] = 0.0        return msd    def L2A_CSHD_1(self):        #Part1 geometric input:        y = self.confidenceMaskCloud.shape[0]        x = self.confidenceMaskCloud.shape[1]        cloud_mask = self.confidenceMaskCloud        filt_b = zeros([y,x], float32)        mask_shadow = zeros([y,x], float32)        # Read azimuth and elevation solar angles        solar_azimuth = -int(self.config.solaz + 0.5)        solar_elevation = int(90.0 - self.config.solze + 0.5)                # Median Filter 7x7        cloud_mask = median_filter(cloud_mask, (7,7))        # Dilatation cross-shape operator        shape = generate_binary_structure(2,1)        cloud_mask = binary_dilation(cloud_mask > 0.33, shape).astype(cloud_mask.dtype)        # Create cloud height distribution (for 30m pixel resolution)        distr_clouds = concatenate([reverse(1. / (1.0 + (arange(51) / 30.0) ** (2 * 5))), 1 / (1.0 + (arange(150) / 90.0) ** (2 * 5))])        # Create projected cloud shadow distribution        npts_shad = distr_clouds.size / tan(solar_elevation * pi / 180.)        factor = npts_shad/distr_clouds.size        distr_shad = zoom(distr_clouds, factor)        # Create filter for convolution (depends on azimuth solar angle)        ys = float(y/2.0)        xs = float(x/2.0)        ds = float(distr_shad.size/2.0)        filt_b[0:distr_shad.size,0] = distr_shad                # keep the original value of first distr_shad pixel        # and mark first distr_shad pixel with -1:        filt_b0 = distr_shad[0]        filt_b[0,0]= -1.0        filt_b[0,1]= -1.0        # Place into center for rotation, subtract 90 degree:        filt_b = roll(filt_b, int(ys-ds), axis=0)        filt_b = roll(filt_b, int(xs), axis=1)        filt_b = reverse(rotate(filt_b, solar_azimuth, reshape=False, order=0))        # identify first distr_shad pixel after rotation and keep these values        # for retranslation. Multiple entries can occurr due to rotation,        # so take the first one:        y0, x0 = where(filt_b < 0.0)        if(y0.size) > 1: y0 = min(y0)        if(x0.size) > 1: x0 = min(x0)        filt_b[filt_b < 0.0] = filt_b0        y0 = y0.astype(int)        x0 = x0.astype(int)        # Now perform the convolution:        fft1 = fft.rfft2(cloud_mask)        fft2 = fft.rfft2(filt_b)        shadow_prob = fft.irfft2(fft1 * fft2)        # Move back to corners:        shadow_prob = roll(shadow_prob, y0, axis=0)        shadow_prob = roll(shadow_prob, x0, axis=1)        CM = self.classificationMask                # Shadow_prob can be smaller as CM, correct this:        dy = CM.shape[0] - shadow_prob.shape[0]        dx = CM.shape[1] - shadow_prob.shape[1]        if(dy > 0):            b = zeros([x,dy],float32)            shadow_prob = concatenate((shadow_prob, b), axis=0)        elif(dx > 0):            b = zeros([y,dx],float32)            shadow_prob = concatenate((shadow_prob, b), axis=1)        # Remove data outside of interest:        shadow_prob[CM == self._noData] = 0        # Normalisation:        shadow_prob = shadow_prob * (1.0 / maximum(shadow_prob.max(), 1.0))        # Remove cloud_mask from Shadow probability:        shadow_prob = maximum((shadow_prob - cloud_mask), 0)        return shadow_prob    def L2A_DarkVegetationRecovery(self):        B04 = self.tables.getBand(self.tables.B04)        B8A = self.tables.getBand(self.tables.B8A)        NDVI = (B8A - B04) / maximum((B8A + B04), self.LOWEST)        T2_NDVI = self.config.getFloat('Scene_Classification/Thresholds', 'T2_NDVI')        F1 = NDVI > T2_NDVI        CM = self.classificationMask        CM[F1 & (CM == self._darkFeatures)] = self._vegetation        CM[F1 & (CM == self._notClassified)] = self._vegetation        T2_R_B8A_B03 = self.config.getFloat('Scene_Classification/Thresholds', 'T2_R_B8A_B03')        B03 = self.tables.getBand(self.tables.B03)        rb8b3 = B8A/maximum(B03,self.LOWEST)        F2 = rb8b3 > T2_R_B8A_B03        CM[F2 & (CM == self._darkFeatures)] = self._vegetation        CM[F2 & (CM == self._notClassified)] = self._vegetation    def L2A_WaterPixelRecovery(self):        B02 = self.tables.getBand(self.tables.B02)        B11 = self.tables.getBand(self.tables.B11)        R_B02_B11 = B02/maximum(B11,self.LOWEST)        T3 = 4.0        F3 = R_B02_B11 > T3        CM = self.classificationMask        CM[F3 & (CM == self._darkFeatures)] = self._water        CM[F3 & (CM == self._notClassified)] = self._water    def L2A_SoilRecovery(self):        T4 = 0.65        B02 = self.tables.getBand(self.tables.B02)        B11 = self.tables.getBand(self.tables.B11)        R_B02_B11 = B02/maximum(B11,self.LOWEST)        F4 = R_B02_B11 < T4        CM = self.classificationMask        CM[F4 & (CM == self._darkFeatures)] = self._bareSoils    def average(self, oldVal, classifier, count):        newVal = self.getClassificationPercentage(classifier)        result = (float32(oldVal) * float32(count) + float32(newVal)) / float32(count + 1.0)           return format('%f' % result)    def getClassificationPercentage(self, classificator):        cm = self.classificationMask        if(classificator == self._noData):            # count all for no data pixels:            nrEntriesTotal = float32(size(cm))            nrEntriesClassified = float32(size(cm[cm == self._noData]))            self._sumPercentage = 0.0        else:            # count percentage of classified pixels:            nrEntriesTotal = float32(size(cm[cm != self._noData]))            nrEntriesClassified = float32(size(cm[cm == classificator]))        fraction = nrEntriesClassified / nrEntriesTotal        percentage = fraction * 100        self._sumPercentage += percentage        self.config.logger.info('Classificator: %d' % classificator)        self.config.logger.info('Percentage: %f' % percentage)        self.config.logger.info('Sum Percentage: %f' % self._sumPercentage)         if(classificator == self._noData):            self._sumPercentage = 0.0                       percentageStr = format('%f' % percentage)        return percentageStr    def updateQualityIndicators(self, nrTilesProcessed, metadata):        xp = L2A_XmlParser(self.config, metadata)        if(nrTilesProcessed == 1) and (self.config.resolution == 60):            # Node must be created:            if(metadata == 'T2A'):                icqi = objectify.Element('L2A_Image_Content_QI')            else:                icqi = objectify.Element('Image_Content_QI')            icqi.NODATA_PIXEL_PERCENTAGE = self.getClassificationPercentage(self._noData)            icqi.SATURATED_DEFECTIVE_PIXEL_PERCENTAGE = self.getClassificationPercentage(self._saturatedDefective)            icqi.DARK_FEATURES_PERCENTAGE = self.getClassificationPercentage(self._darkFeatures)            icqi.CLOUD_SHADOW_PERCENTAGE = self.getClassificationPercentage(self._cloudShadows)            icqi.VEGETATION_PERCENTAGE = self.getClassificationPercentage(self._vegetation)            icqi.BARE_SOILS_PERCENTAGE = self.getClassificationPercentage(self._bareSoils)            icqi.WATER_PERCENTAGE = self.getClassificationPercentage(self._water)            icqi.LOW_PROBA_CLOUDS_PERCENTAGE = self.getClassificationPercentage(self._lowProbaClouds)            icqi.MEDIUM_PROBA_CLOUDS_PERCENTAGE = self.getClassificationPercentage(self._medProbaClouds)            icqi.HIGH_PROBA_CLOUDS_PERCENTAGE = self.getClassificationPercentage(self._highProbaClouds)            icqi.THIN_CIRRUS_PERCENTAGE = self.getClassificationPercentage(self._thinCirrus)            icqi.SNOW_ICE_PERCENTAGE = self.getClassificationPercentage(self._snowIce)            icqi.RADIATIVE_TRANSFER_ACCURAY = 0.0            icqi.WATER_VAPOUR_RETRIEVAL_ACCURACY = 0.0            icqi.AOT_RETRIEVAL_ACCURACY = 0.0            if(metadata == 'T2A'):                qii = xp.getRoot('Quality_Indicators_Info')                qii.insert(1, icqi)            else:                qii = xp.getRoot('L2A_Quality_Indicators_Info')                qii.append(icqi)        else: # Node exists and has to be updated:            if(metadata == 'T2A'):                icqi = xp.getTree('Quality_Indicators_Info', 'L2A_Image_Content_QI')            else:                icqi = xp.getTree('L2A_Quality_Indicators_Info', 'Image_Content_QI')            icqi.NODATA_PIXEL_PERCENTAGE = self.average(icqi.NODATA_PIXEL_PERCENTAGE, self._noData, nrTilesProcessed)            icqi.SATURATED_DEFECTIVE_PIXEL_PERCENTAGE = self.average(icqi.SATURATED_DEFECTIVE_PIXEL_PERCENTAGE, self._saturatedDefective, nrTilesProcessed)            icqi.DARK_FEATURES_PERCENTAGE = self.average(icqi.DARK_FEATURES_PERCENTAGE, self._darkFeatures, nrTilesProcessed)            icqi.CLOUD_SHADOW_PERCENTAGE = self.average(icqi.CLOUD_SHADOW_PERCENTAGE, self._cloudShadows, nrTilesProcessed)            icqi.VEGETATION_PERCENTAGE = self.average(icqi.VEGETATION_PERCENTAGE, self._vegetation, nrTilesProcessed)            icqi.BARE_SOILS_PERCENTAGE = self.average(icqi.BARE_SOILS_PERCENTAGE, self._bareSoils, nrTilesProcessed)            icqi.WATER_PERCENTAGE = self.average(icqi.WATER_PERCENTAGE, self._water, nrTilesProcessed)            icqi.LOW_PROBA_CLOUDS_PERCENTAGE = self.average(icqi.LOW_PROBA_CLOUDS_PERCENTAGE, self._lowProbaClouds, nrTilesProcessed)            icqi.MEDIUM_PROBA_CLOUDS_PERCENTAGE = self.average(icqi.MEDIUM_PROBA_CLOUDS_PERCENTAGE, self._medProbaClouds, nrTilesProcessed)            icqi.HIGH_PROBA_CLOUDS_PERCENTAGE = self.average(icqi.HIGH_PROBA_CLOUDS_PERCENTAGE, self._highProbaClouds, nrTilesProcessed)            icqi.THIN_CIRRUS_PERCENTAGE = self.average(icqi.THIN_CIRRUS_PERCENTAGE, self._thinCirrus, nrTilesProcessed)            icqi.SNOW_ICE_PERCENTAGE = self.average(icqi.SNOW_ICE_PERCENTAGE, self._snowIce, nrTilesProcessed)            icqi.RADIATIVE_TRANSFER_ACCURAY = 0.0            icqi.WATER_VAPOUR_RETRIEVAL_ACCURACY = 0.0            icqi.AOT_RETRIEVAL_ACCURACY = 0.0        xp.export()    def process(self):        ts = time.time()        self.config.timestamp('Pre process  ')        self.preprocess()        self.config.timestamp('L2A_SC init  ')        self.L2A_CSND_1_1()        self.config.timestamp('L2A_CSND_1_1 ')        self.L2A_CSND_1_2()        self.config.timestamp('L2A_CSND_1_2 ')        if(self.tables.sceneCouldHaveSnow() == True):            self.config.logger.info('Snow probality, detection will be performed')            self.L2A_CSND_2_0()            self.config.timestamp('L2A_CSND_2_0 ')            self.L2A_CSND_2_1()            self.config.timestamp('L2A_CSND_2_1 ')            self.L2A_CSND_2_2()            self.config.timestamp('L2A_CSND_2_2 ')            self.L2A_CSND_2_3()            self.config.timestamp('L2A_CSND_2_3 ')            self.L2A_CSND_2_4()            self.config.timestamp('L2A_CSND_2_4 ')            self.L2A_CSND_2_5()            self.config.timestamp('L2A_CSND_2_5 ')        else:            self.config.logger.info('Now snow probality, detection will be ignored')        self.L2A_CSND_3()        self.config.timestamp('L2A_CSND_3   ')        self.L2A_CSND_4()        self.config.timestamp('L2A_CSND_4   ')        self.L2A_CSND_5_1()        self.config.timestamp('L2A_CSHD_5_1 ')        self.L2A_CSND_5_2()        self.config.timestamp('L2A_CSND_5_2 ')        self.L2A_CSND_6()        self.config.timestamp('L2A_CSND_6   ')        self.L2A_CSND_7()        self.config.timestamp('L2A_CSND_7   ')        self.L2A_SHD()        self.config.timestamp('L2A_SHD      ')        self.L2A_SoilRecovery()        self.config.timestamp('Soil recovery')        self.L2A_DarkVegetationRecovery()        self.config.timestamp('DV recovery  ')        self.L2A_WaterPixelRecovery()        self.config.timestamp('WP recovery  ')        self.postprocess()        self.config.timestamp('Post process ')        tDelta = time.time() - ts        self.config.logger.info('Procedure L2A_SceneClass overall time [s]: %0.3f' % tDelta)        if(self.config.traceLevel == 'DEBUG'):            stdoutWrite('Procedure L2A_SceneClass, overall time[s]: %0.3f.\n' % tDelta)        return True