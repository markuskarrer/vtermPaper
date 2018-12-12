'''
here are some functions to postprocess the radar data (mainly taken from Jose)
'''
import numpy as np
import pandas as pd

def getVar( var, tempDataSet,fileList):
    
    for filePath in joyrad10FileList:
    
        tempDS = xr.open_dataset(filePath)
        tempDS.time.attrs['units']='seconds since 1970-01-01 00:00:00 UTC'
        tempDS = xr.decode_cf(tempDS)
        tempDSVar = tempDS[var]#.where(tempDS['Zg']>0)
        #tempDSZe = 10*np.log10(tempDSZe) 
    
        tempDataSet = xr.merge([tempDataSet, tempDSVar])
        
      
    #joyrad10[var].attrs = tempDSVar.attrs    
    return tempDataSet

def calcRadarDeltaGrid(refGrid, radarGrid):
    
    radGrid2d = np.ones((len(refGrid),
                         len(radarGrid)))*radarGrid
    
    deltaGrid = radGrid2d - np.reshape(refGrid,(len(refGrid),1))
    
    return deltaGrid

def getNearestIndexM2(deltaGrid, tolerance):

    gridIndex = np.argmin(abs(deltaGrid), axis=1)
    deltaGridMin = np.min(abs(deltaGrid), axis=1)
    gridIndex = np.array(gridIndex, np.float)
    gridIndex[deltaGridMin>tolerance] = np.nan
    
    return gridIndex

def getResampledVar(var, xrDataset, timeIndexArray, rangeIndexArray):

    resampledArr = np.ones((timeIndexArray.shape[0],rangeIndexArray.shape[0]))*np.nan
    resampledTimeArr = np.ones((timeIndexArray.shape[0], xrDataset.range.values.shape[0]))*np.nan

    for t, timeIndex in enumerate(timeIndexArray):

        try:
            resampledTimeArr[t]=  xrDataset[var].values[int(timeIndex)]

        except:
            pass


    resampledTimeArrTra = resampledTimeArr.T
    resampledArr = resampledArr.T

    for r, rangeIndex in enumerate(rangeIndexArray):

            try:
                resampledArr[r] = resampledTimeArrTra[int(rangeIndex)]
            except:
                pass

    return resampledArr.T

def getTimeRef(startDay, endDay):

    
    #start = pd.datetime(date.year,
    #                   date.month,
    #                   date.day, 
    #                    0, 0, 0)

    #end = pd.datetime(date.year,
    #                  date.month,
    #                  date.day,
    #                  23, 59, 59)

    timeRef = pd.date_range(startDay, endDay, freq='2s')
    
    return timeRef

def getCutSpec(dataSet):

    newDoppler=np.ones(524)
    newDoppler[0:262]=dataSet.doppler.values[3834:]
    newDoppler[262:]=dataSet.doppler.values[0:262]
    
    newSpect = np.ones((dataSet.time.shape[0],
                        dataSet.range.shape[0],
                        524)) * np.nan 

    newSpect[:,:,:262]=dataSet.SPCco.values[:,:,3834:]
    newSpect[:,:,262:]=dataSet.SPCco.values[:,:,:262]
    
    return newDoppler, newSpect

def getCalSpect(dataSet, specArr):
    
    rangeLen = dataSet.range.shape[0]
    timeLen = dataSet.time.shape[0]
    
    rangeData = dataSet.range.values.reshape(1, rangeLen, 1)
    radConst = dataSet.RadarConst.values.reshape(timeLen,1,1)
    SNRCorFaCo = dataSet.SNRCorFaCo.values.reshape(timeLen,rangeLen,1)
    SPCco = specArr
    npw = dataSet.npw1.values.reshape(timeLen,1,1)
    
    calSPC = radConst*SNRCorFaCo*(rangeData**2/5000.**2)*SPCco/npw
    
    return calSPC

def getCalNoise(dataSet):
    
    rangeLen = dataSet.range.shape[0]
    timeLen = dataSet.time.shape[0]
    
    rangeData = dataSet.range.values.reshape(1, rangeLen, 1)
    radConst = dataSet.RadarConst.values.reshape(timeLen,1,1)
    SNRCorFaCo = dataSet.SNRCorFaCo.values.reshape(timeLen,rangeLen,1)
    HSDco = dataSet.HSDco.values.reshape(timeLen, rangeLen, 1)
    npw = dataSet.npw1.values.reshape(timeLen,1,1)
    
    calNoise = radConst*SNRCorFaCo*(rangeData**2/5000.**2)*HSDco/npw
    
    return calNoise

###
def read_tripex_pol_radar(obsData,year,month,day,hourBegin,hourEnd,i_freq=0,minBegin='00',minEnd='00',secBegin='00',secEnd='00',radar_name='Ka',onlymoments=False):
    '''
    this function reads the data from the tripex_pol campaign for the three vertically pointing radars
    obsData is the dictionary to write in the data (can contain already some data
    radar_name: specifies which radar to take 
    onlymoments: read only the moments to save time if the spectra is not relevant
    '''
    import glob
    from netCDF4 import Dataset
    import xarray as xr

    ###
    #get indices for the time array
    ###
    #define start and end as pandas time
    start = pd.datetime(int(year), int(month), int(day),
                        int(hourBegin), int(minBegin),
                        int(secBegin))


    end = pd.datetime(int(year), int(month), int(day),
                        int(hourEnd), int(minEnd), int(secEnd))
    
    print "reading data from: " + radar_name

    if radar_name in ('joyrad10','joyrad35'):
        Z_in_linear_units = True
        if radar_name=='joyrad10':
            joyradPath = '/data/obs/site/jue/joyrad10'
            filePath=''
            while(filePath==''): #find file which contains this timestep
                filePath = ('/').join([joyradPath, year, month, day, year+month+day +'_' + hourBegin + minBegin +'??.znc']) #for Ka-Band just 1200.znc #for X-Band 1200??.znc
                print filePath
                try:
                    filePath = glob.glob(filePath)[0]
                except:
                    hourBegin = str(int(hourBegin)-1).zfill(2)
                    print filePath + ' not found; try with one hour earlier'
                    filePath=''
                    print filePath

        elif radar_name=='joyrad35':
            Z_in_linear_units = True
            joyradPath = '/data/obs/site/jue/joyrad35' 
            filePath=''
            while(filePath==''): #find file which contains this timestep
                filePath = ('/').join([joyradPath, year, month, day, year+month+day +'_' + hourBegin + minBegin +'.znc']) #for Ka-Band just 1200.znc #for X-Band 1200??.znc
                try:
                    filePath = glob.glob(filePath)[0]
                except:
                    print filePath + ' not found; try with one hour earlier'
                    filePath=''
                    print filePath
                    hourBegin = str(int(hourBegin)-1).zfill(2)
        joyrad = xr.open_dataset(filePath)
        joyrad.time.attrs['units']='seconds since 1970-01-01 00:00:00 UTC'
        joyrad = xr.decode_cf(joyrad)
        joyrad = joyrad.sel(time=slice(start, end)).copy()


        #read times from radar-data
        times = joyrad.time
        epoch = pd.DatetimeIndex([pd.datetime(2001,1,1,0,0,0)]).astype(np.int)/10**9
        radarTimes = times+epoch.values
        refTime_startaverage = pd.DatetimeIndex([start]).astype(np.int)/10**9
        refTime_endaverage = pd.DatetimeIndex([end]).astype(np.int)/10**9

        index_start = (np.abs(np.array(radarTimes.values,float)/10**9 - refTime_startaverage.values)).argmin()
        index_end = (np.abs(np.array(radarTimes.values,float)/10**9 - refTime_endaverage.values)).argmin()    
        if radar_name=='joyrad10':
            obsData["height"] = joyrad.range+20.0
        elif radar_name=='joyrad35':
            obsData["height"] = joyrad.range
        
        '''
        ##copy and average (median) moments
        obsData["Ze"] = np.nanmedian(joyrad.Zg[index_start:index_end],axis=0); obsData["Ze"] = obsData["Ze"][:,None] #TODO: include dimension for frequency
        obsData["Radar_MeanDopplerVel"] = np.nanmedian(joyrad.VELg[index_start:index_end],axis=0) ; obsData["Radar_MeanDopplerVel"] = -obsData["Radar_MeanDopplerVel"][:,None]
        obsData["Radar_SpectrumWidth"] = np.nanmedian(joyrad.RMSg[index_start:index_end],axis=0); obsData["Radar_SpectrumWidth"] = obsData["Radar_SpectrumWidth"][:,None]
        #obsData["Radar_Skewness"] = np.nanmedian(-joyrad[index_start:index_end],axis=0) ; obsData["Radar_Skewness"] = obsData["Radar_Skewness"][:,None] #no skewness directly available
        obsData["Radar_Skewness"] = np.ones_like(obsData["Radar_SpectrumWidth"][:,None])*np.nan
        '''
        if i_freq==0:
            obsData["Ze"] = np.nanmedian(joyrad.Zg[index_start:index_end],axis=0); obsData["Ze"] = obsData["Ze"][:,None] #TODO: include dimension for frequency
            obsData["Radar_MeanDopplerVel"] = np.nanmedian(joyrad.VELg[index_start:index_end],axis=0) ; obsData["Radar_MeanDopplerVel"] = -obsData["Radar_MeanDopplerVel"][:,None]
            obsData["Radar_SpectrumWidth"] = np.nanmedian(joyrad.RMSg[index_start:index_end],axis=0); obsData["Radar_SpectrumWidth"] = obsData["Radar_SpectrumWidth"][:,None]
            #obsData["Radar_Skewness"] = np.nanmedian(-joyrad[index_start:index_end],axis=0) ; obsData["Radar_Skewness"] = obsData["Radar_Skewness"][:,None] #no skewness directly available
            obsData["Radar_Skewness"] = np.ones_like(obsData["Radar_SpectrumWidth"][:,None])*np.nan
        else: #this is not used (!?)
            Ze_now = np.nanmedian(joyrad.Zg[index_start:index_end],axis=0); obsData["Ze"][:,i_freq] = Ze_now #TODO: include dimension for frequency
            MeanDopplerVel_now = np.nanmedian(joyrad.VELg[index_start:index_end],axis=0) ; obsData["Radar_MeanDopplerVel"][:,i_freq] = MeanDopplerVel_now
            MeanSpectrumWidth_now = np.nanmedian(joyrad.RMSg[index_start:index_end],axis=0); obsData["Radar_SpectrumWidth"][:,i_freq] = MeanSpectrumWidth_now
            MeanSkewness_now = np.ones_like(MeanSpectrumWidth_now) ; obsData["Radar_Skewness"][:,i_freq] = MeanSkewness_now

        #convert from SI to dBz if necessary
        if Z_in_linear_units:
            obsData["Ze"] = 10*np.log10(obsData["Ze"])
        if radar_name=='joyrad35':
            offset = 3.0
            offset_inlinunits = 10**(offset/10)

            print "applied offset correction of " + str(offset) + "dB for " + radar_name
            obsData["Ze"] += offset
        elif radar_name=='joyrad10':
            print 'offset not defined for ' + radar_name
            offset= 0.0
            offset_inlinunits = 10**(0.0/10)

        if not onlymoments: #TODO: does this work here?
            #read in the spectra
            calSpect = getCalSpect(joyrad, joyrad.SPCco.values)
            calNoise = getCalNoise(joyrad)
            specCalNoise = calSpect - (np.ones_like(calSpect)*calNoise)
            #from IPython.core.debugger import Tracer ; Tracer()()
            obsData["Radar_Velocity"] = joyrad.doppler; obsData["Radar_Velocity"] = np.expand_dims(obsData["Radar_Velocity"],axis=0)
            obsData["Radar_Spectrum"] = specCalNoise[index_start,:]; obsData["Radar_Spectrum"] = np.expand_dims(obsData["Radar_Spectrum"],axis=1)
            if Z_in_linear_units:
                obsData["Radar_Spectrum"] = 10*np.log10(obsData["Radar_Spectrum"])+offset
    elif radar_name in ('grarad94','Wbandtripex'):
        if radar_name=='grarad94':
            dataPath = ('/').join(['/data/obs/campaigns/tripex-pol/wband_gra/l1',year,month,day])
            fileName = ('_').join([year[2:]+month+day+'_'+hourBegin+'*ZEN.nc'])

            #merge the strings and find the files
            filePath = ('/').join([dataPath, fileName])
            newFilePath = glob.glob(filePath)[0]
            Z_in_linear_units=True
            epoch_start='2001'
            chirpSeqmissing= True
        elif radar_name=='Wbandtripex':
            dataPath = ('/').join(['/data/hatpro/jue/data/joyrad94/l0',year+month,day])
            fileName = ('_').join(['joyrad94_joyce', year+month+day+hourBegin+'.nc'])
            #merge the strings and find the files
            filePath = ('/').join([dataPath, fileName])
            FilePath = filePath #glob.glob(filePath)[0]
            
            Z_in_linear_units=False
            epoch_start='2001'
            chirpSeqmissing= True

        #load the data
        rootgrp = Dataset(newFilePath)
        

        
        times = rootgrp.variables['time'][:]
        epoch = pd.DatetimeIndex([pd.datetime(2001,1,1,0,0,0)]).astype(np.int)/10**9
        radarTimes = times+epoch.values
        refTime_startaverage = pd.DatetimeIndex([start]).astype(np.int)/10**9
        refTime_endaverage = pd.DatetimeIndex([end]).astype(np.int)/10**9

        index_start = (np.abs(radarTimes - refTime_startaverage.values)).argmin()
        index_end = (np.abs(radarTimes - refTime_endaverage.values)).argmin()
            
        #read in the moments
        Ze = 10*np.log10(np.array(rootgrp.variables['Ze']))
        vDoppler =  np.array(rootgrp.variables['vm']) #ms-1
        swidth = np.array(rootgrp.variables['sigma']) #ms-1
        skew = np.array(rootgrp.variables['skew']) #1
        if i_freq==0:
            ##copy and average (median) moments
            obsData["Ze"] = np.nanmedian(Ze[index_start:index_end],axis=0); obsData["Ze"] = obsData["Ze"][:,None] #TODO: include dimension for frequency
            ###
            #add offset
            ###
            if radar_name=='grarad94':
                offset = -2.
                print "applied offset correction of " + str(offset) + "dB for " + radar_name
                obsData["Ze"] += offset
            elif radar_name=='Wbandtripex':
                print 'offset not defined'
                offset = 0.0
            obsData["Radar_MeanDopplerVel"] = np.nanmedian(vDoppler[index_start:index_end],axis=0) ; obsData["Radar_MeanDopplerVel"] = -obsData["Radar_MeanDopplerVel"][:,None]
            obsData["Radar_SpectrumWidth"] = np.nanmedian(swidth[index_start:index_end],axis=0); obsData["Radar_SpectrumWidth"] = obsData["Radar_SpectrumWidth"][:,None]
            obsData["Radar_Skewness"] = np.nanmedian(-skew[index_start:index_end],axis=0) ; obsData["Radar_Skewness"] = obsData["Radar_Skewness"][:,None]
        else:
            Ze_now = np.nanmedian(Ze[index_start:index_end],axis=0); obsData["Ze"][:,i_freq] = Ze_now #TODO: include dimension for frequency
            MeanDopplerVel_now = np.nanmedian(vDoppler[index_start:index_end],axis=0) ; obsData["Radar_MeanDopplerVel"][:,i_freq] = MeanDopplerVel_now
            MeanSpectrumWidth_now = np.nanmedian(swidth[index_start:index_end],axis=0); obsData["Radar_SpectrumWidth"][:,i_freq] = MeanSpectrumWidth_now
            MeanSkewness_now = np.nanmedian(-skew[index_start:index_end],axis=0) ; obsData["Radar_Skewness"][:,i_freq] = MeanSkewness_now
        if not onlymoments:
            #read in the spectra
        ###
            #read in the variables
            ###
            ranges = rootgrp.variables['range'][:] #heights
            if not chirpSeqmissing:
                chirpSeq = rootgrp.variables['chirp_sequences'][:]
            elif chirpSeqmissing: #so far this is not properly saved in the tripex-pol files
                chirpSeq = np.arange(1,5)
            rangeOff = rootgrp.variables['range_offsets'][:]
            rangeIndex = rangeOff-1

            vel = rootgrp.variables['velocity'][:] #velocity array
                
            #something with chirp
            try: #sometimes there are only 3 chirp sequences
                firstRange = np.ones([rangeIndex[1]-rangeIndex[0],vel.shape[1]])*vel[0]
                secondRange = np.ones([rangeIndex[2]-rangeIndex[1],vel.shape[1]])*vel[1]
                thirdRange = np.ones([rangeIndex[3]-rangeIndex[2],vel.shape[1]])*vel[2]
                fourthRange = np.ones([ranges.shape[0]-rangeIndex[3],vel.shape[1]])*vel[3]
                #combining the chirp-sequences?
                vel2Darr = np.concatenate([firstRange, secondRange, thirdRange, fourthRange])
            except:
                firstRange = np.ones([rangeIndex[1]-rangeIndex[0],vel.shape[1]])*vel[0]
                secondRange = np.ones([rangeIndex[2]-rangeIndex[1],vel.shape[1]])*vel[1]
                thirdRange = np.ones([ranges.shape[0]-rangeIndex[2],vel.shape[1]])*vel[2]
                #combining the chirp-sequences?
                vel2Darr = np.concatenate([firstRange, secondRange, thirdRange])
            ranges2D = (np.ones_like(vel2Darr).T*ranges).T/10**3
            ranges2D= np.ma.masked_invalid(ranges2D)
            vel2Darr = np.ma.masked_invalid(vel2Darr)

            #read the spectrum and mask low values
            spec = rootgrp.variables['spec'][:]
                #select spectra
            newSpec = spec[index_start]        
            #convert from SI to dBz if necessary
            if Z_in_linear_units:
                newSpec = 10*np.log10(newSpec)
                
            #copy the Spectra to the obsData dictionary
            obsData["Radar_Velocity"] = vel2Darr; obsData["Radar_Velocity"] = -obsData["Radar_Velocity"][None,:]  #TODO: if there are different frequencies None should be replaced by some index
            obsData["Radar_Spectrum"] = newSpec; obsData["Radar_Spectrum"] = obsData["Radar_Spectrum"][:,None,:]  #TODO: if there are different frequencies None should be replaced by some index
            #correct offset
            obsData["Radar_Spectrum"] += offset
            #obsData["frequency"] = np.array([94.0]) #this should fit to the next line
            obsData["height"] = ranges2D[:,0]*1000;
    
    return obsData


def convert_linear_to_dBz(Ze_lin):
    return 10*np.log10(Ze_lin)