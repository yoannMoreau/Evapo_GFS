#-*- coding: utf-8 -*-
'''
Created on 16 déc. 2013

@author: yoann Moreau

All controls operations :
return true if control ok 
'''
import os
import errno
from datetime import date,datetime,timedelta
import ogr,osr
import re
import gdal
import numpy as np
import subprocess
import shutil
import math
import scipy.ndimage as ndimage
import pygrib
import urllib2
import pyproj as pp
    
def checkForFile(pathToFile):
    if os.path.isfile(pathToFile):
        return True
    else:
        return False

def createParamFile(pathFile,user,key):
    
    f = open(pathFile, 'w+')
    f.write("{\n")
    f.write(' "url"   : "https://api.ecmwf.int/v1",\n')
    f.write('"key"   : "'+key+'",\n')
    f.write('"email" : "'+user+'"\n')
    f.write("}")
    f.close()
    

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def checkForFolder(pathToFolder):
    try:
        os.makedirs(pathToFolder)
    except OSError as exception:
        if exception.errno != errno.EEXIST:            
            exit('Path for downloaded Era Interim could not be create. Check your right on the parent folder...')
            
def checkForDate(dateC):
    #convert string to date from YYYY-MM-DD
    if len(dateC)==10:
        YYYY=dateC[0:4]
        MM=dateC[5:7]
        DD=dateC[8:10]
        if (YYYY.isdigit() and MM.isdigit()  and DD.isdigit()):
            try:
                date(int(YYYY),int(MM),int(DD))
            except ValueError:
                exit('Error on Date Format... please give a date in YYYY-MM-DD format')
            
            return date(int(YYYY),int(MM),int(DD))

        else:
            exit('Error on Date Format... please give a date in YYYY-MM-DD format')
    else: 
        exit('Error on Date Format... please give a date in YYYY-MM-DD format')


def convertShpToExtend(pathToShp):
    """
    reprojette en WGS84 et recupere l'extend
    """ 
    
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataset = driver.Open(pathToShp)
    if dataset is not None:
        # from Layer
        layer = dataset.GetLayer()
        spatialRef = layer.GetSpatialRef()
        # from Geometry
        feature = layer.GetNextFeature()
        geom = feature.GetGeometryRef()
        spatialRef = geom.GetSpatialReference()
        
        #WGS84
        outSpatialRef = osr.SpatialReference()
        outSpatialRef.ImportFromEPSG(4326)

        coordTrans = osr.CoordinateTransformation(spatialRef, outSpatialRef)

        env = geom.GetEnvelope()

        pointMAX = ogr.Geometry(ogr.wkbPoint)
        pointMAX.AddPoint(env[1], env[3])
        pointMAX.Transform(coordTrans)
        
        pointMIN = ogr.Geometry(ogr.wkbPoint)
        pointMIN.AddPoint(env[0], env[2])
        pointMIN.Transform(coordTrans)


        return [pointMAX.GetPoint()[1],pointMIN.GetPoint()[0],pointMIN.GetPoint()[1],pointMAX.GetPoint()[0]]
    else:
        exit(" shapefile not found. Please verify your path to the shapefile")


def is_float_re(element):
    _float_regexp = re.compile(r"^[-+]?(?:\b[0-9]+(?:\.[0-9]*)?|\.[0-9]+\b)(?:[eE][-+]?[0-9]+\b)?$").match
    return True if _float_regexp(element) else False


def checkForExtendValidity(extendList):
    
    if len(extendList)==4 and all([is_float_re(str(x)) for x in extendList]) and extendList[0]>extendList[2] and extendList[1]<extendList[3]:
        if float(extendList[0]) > -180 and float(extendList[2]) <180 and float(extendList[1]) <90 and  float(extendList[3]) > -90:
            extendArea=[str(x) for x in extendList]
            return extendArea
        else:
            exit('Projection given is not in WGS84. Please verify your -t parameter')
    else:
        exit('Area specified is not conform to a  ymax xmin ymin xmax  extend. please verify your declaration')

def checkForStepValidity(listStep,typeData):
    
    validParameters=(0,6,12,18)
    if typeData=="forecast":
        if len(listStep)>0 and isinstance(listStep, list) and all([int(x) in validParameters for x in listStep]):
            listStep=[int(x) for x in listStep]
            return listStep
        else: 
            exit('step parameters not conform to eraInterim posibility : '+ ",".join([str(x) for x in validParameters]))
    else:
        if len(listStep)>0:
            exit('step parameters not conform to eraInterim posibility : '+ ",".join([str(x) for x in validParameters])+ 'for analyse')
        else:
            return listStep

def checkForGridValidity(grid):
    
    if (is_float_re(grid)):
        grid=float(grid)
        validParameters=(0.25,0.5,1,2.5)
        
        if grid in validParameters:
            return grid
        else:
            exit('grid parameters not conform to eraInterim posibility : '+ ",".join([str(x) for x in validParameters]))
    else:
        exit('grid parameters not conform to eraInterim posibility : '+ ",".join([str(x) for x in validParameters]))
    

def create_request_gfs(dateStart,dateEnd,stepList,levelList,grid,extent,paramList,typeData):
    """
        Genere la structure de requete pour le téléchargement de données GFS
        
        INPUTS:\n
        -date : au format annee-mois-jour\n
        -heure : au format heure:minute:seconde\n
        -coord : une liste des coordonnees au format [N,W,S,E]\n
        -dim_grille : taille de la grille en degree \n
    """
    
    URLlist=[]
    #Control datetype
    listForcastSurface=['GUST','HINDEX','PRES','HGT','TMP','WEASD','SNOD','CPOFP','WILT','FLDCP','SUNSD','LFTX','CAPE','CIN','4LFTX','HPBL','LAND']
    if (0 not in [int(x) for x in stepList]):
        listForcastSurfaceExt=listForcastSurface+['PEVPR','CPRAT','PRATE','APCP','ACPCP','WATR','CSNOW','CICEP','CFPER','CRAIN','LHTFL','SHTFL','SHTFL','GFLUX','UFLX','VFLX','U-GWD','V-GWD','DSWRF','DLWRF','ULWRF','USWRF','ALBDO']
    listAnalyseSurface=['HGT','PRES','LFTX','CAPE','CIN','4LFTX']
    
    if typeData == 'analyse' and all([x in listAnalyseSurface for x in paramList]):
        typeData= 'analyse'
        validChoice = None
        prbParameters =  None
    else :
        if all([x in listForcastSurface for x in paramList]) and typeData != 'cycleforecast':
            typeData= 'forecast'
            validChoice = typeData
            indexParameters=[i for i, elem in enumerate([x in listAnalyseSurface for x in paramList], 1) if not elem]
            prbParameters=[]
            for i in indexParameters:
                prbParameters.append(paramList[i-1])
        else:
            typeData= 'cycleforecast'
            validChoice = typeData
            indexParameters=[i for i, elem in enumerate([x in listAnalyseSurface for x in paramList], 1) if not elem]
            prbParameters=[]
            for i in indexParameters:
                prbParameters.append(paramList[i-1])
    
      
    #Control si date/timeList disponible
    today=date.today()
    lastData = today - timedelta(days=14)
    if dateStart < lastData or dateEnd > today : 
        exit('date are not in 14 days range from today' )
    else:
        #Pour chaque jour souhaité
        nbDays=(dateEnd -dateStart).days+1 
        for i in range(0,nbDays):
            #on crontrole pour les timeList
            if dateStart + timedelta(days=i) == today:
                maxT=datetime.now().hour-5
                timeListCorr=[ x for x in stepList if x<maxT ]
            else:
                timeListCorr=stepList
              
            for t in timeListCorr:
                URL='http://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_'
                #grid
                URL=URL+"{:.2f}".format(grid).replace('.','p')+'.pl?file=gfs.'
                #time ( attention limiter avec décalage horaire for today
                URL=URL+'t'+str(t).zfill(2)+'z.'
                if (grid==0.5):
                    URL=URL+'pgrb2full.'
                else:
                    URL=URL+'pgrb2.'
                URL=URL+"{:.2f}".format(grid).replace('.','p')+'.'
                
                if typeData=='cycleforecast':
                    URL=URL+'f006&lev_'
                elif typeData=='forecast':
                    URL=URL+'f000&lev_'
                else:
                    URL=URL+'anl&lev_'
                URL=URL+"=on&lev_".join(levelList)+"=on&var_"
                URL=URL+"=on&var_".join(paramList)+"=on&subregion=&"
                URL=URL+"leftlon="+str(round(float(extent[1])-0.05,1))+"&rightlon="+str(round(float(extent[3])+0.05,1))+"&toplat="+str(round(float(extent[0])+0.5,1))+"&bottomlat="+str(round(float(extent[2])-0.5,1))
                URL=URL+"&dir=%2Fgfs."+"{:%Y%m%d}".format(dateStart+timedelta(days=i))+str(t).zfill(2)
                URLlist.append(URL)
        
        return (URLlist,validChoice,prbParameters)

def GFSDownload(pathToFile,pathToOutputFile):

    response = urllib2.urlopen(pathToFile)
    
    try:
        html = response.read()
    except:
        exit("error while downloading file")
    
    if len(html) > 0:
        f = open(pathToOutputFile, 'wb')
        f.write(html)
        return True
    else:
        return False

def moveFile(inputImg,outputImg):
        "move a file to a directory"
        
        #TODO controls to check if exist
        #on déplace le fichier dans le bon répertoire
        shutil.move(inputImg, outputImg)

def reprojRaster(pathToImg,output,shape,pathToShape):
    
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(pathToShape, 0)
    layer = dataSource.GetLayer()
    srs = layer.GetSpatialRef()
    
    Xres=shape[1]
    Yres=shape[0]
    
    subprocess.call(["gdalwarp","-q","-s_srs","EPSG:4326","-t_srs",srs.ExportToWkt(),pathToImg,output,'-ts',str(Xres),str(Yres),'-overwrite','-dstnodata',"0"])
    return output

def writeTiffFromDicoArray(DicoArray,outputImg,shape,geoparam,proj=None,format=gdal.GDT_Float32):
    
    gdalFormat = 'GTiff'
    driver = gdal.GetDriverByName(gdalFormat)

    dst_ds = driver.Create(outputImg, shape[1], shape[0], len(DicoArray), format)
    
    j=1
    for i in DicoArray.values():
        dst_ds.GetRasterBand(j).WriteArray(i, 0)
        band = dst_ds.GetRasterBand(j)
        band.SetNoDataValue(0)
        j+=1
    
    originX =  geoparam[0]
    originY =  geoparam[1]
    pixelWidth = geoparam[2]
    pixelHeight  = geoparam[3]

    dst_ds.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    if proj==None:
        spatialRef = osr.SpatialReference()
        spatialRef.ImportFromEPSG(4326) 
        dst_ds.SetProjection(spatialRef.ExportToWkt())

def convertGribToDicoArray(listeFile,listParam,listLevel,liststep,grid,startDate,endDate):
    """ Convert GRIB to Tif"""
    
    dicoValues={}
    
    for l in listeFile:
        grbs = pygrib.open(l)
        grbs.seek(0)
        index=1
        for j in range(len(listLevel),0,-1):
            for i in range(len(listParam)-1,-1,-1):
                grb = grbs[index]
                p=grb.name.replace(' ','_')
                if grb.level != 0:
                    l=str(grb.level)+'_'+grb.typeOfLevel
                else:
                    l=grb.typeOfLevel
                if p+'_'+l not in dicoValues.keys():
                    dicoValues[p+'_'+l]=[]
                dicoValues[p+'_'+l].append(grb.values)
                shape=grb.values.shape
                lat,lon=grb.latlons()
                geoparam=(lon.min(),lat.max(),grid,grid)
                index+= 1

    nbJour=(endDate-startDate).days+1
    #on joute des arrayNan si il manque des fichiers
    for s in range(0, (len(liststep)*nbJour-len(listeFile))):
        for k in dicoValues.keys():
            dicoValues[k].append(np.full(shape, np.nan))

    #On écrit pour chacune des variables dans un fichier
    for i in range(len(dicoValues.keys())-1,-1,-1):
        dictParam=dict((k,dicoValues[dicoValues.keys()[i]][k]) for k in range(0,len(dicoValues[dicoValues.keys()[i]])))
        sorted(dictParam.items(), key=lambda x: x[0])
        
    return dictParam

def convertKToD(dicoT):
    
    Degree={}
    
    for i in dicoT:
        mask=np.logical_not(dicoT[i] > 0).astype(int)
        DegreeArray=dicoT[i]-273.15
        np.putmask(DegreeArray,mask,np.nan)
        Degree[i]=DegreeArray
    return Degree

def convertToHectoPascal(dicoP):
    
    pressure= {}
    
    for i in dicoP:
        mask=np.logical_not(dicoP[i] > 0).astype(int)
        PArray=dicoP[i]/100
        np.putmask(PArray,mask,np.nan)
        pressure[i]=PArray
    return pressure

def convertPaToKgPa(dicoP):
    
    pressure= {}
    
    for i in dicoP:
        mask=np.logical_not(dicoP[i] > 0).astype(int)
        PArray=dicoP[i]/1000
        np.putmask(PArray,mask,np.nan)
        pressure[i]=PArray
    return pressure

def convertMToMm(dicoEvap):
    evap={}
    for i in dicoEvap:
        RArray=(dicoEvap[i]*(10**3))

        evap[i]=RArray
        
    return evap

def convertWToMJ(dicoRay):

    rayonnement={}
    for i in dicoRay:
        mask=np.logical_not(dicoRay[i] > 0).astype(int)
        RArray=(dicoRay[i]*6*3600/(10**6))
        np.putmask(RArray,mask,np.nan)
        rayonnement[i]=RArray
        
    return rayonnement

def convertGeoToAlt(dicoGeo):
    
    def mean(values):
        return np.nanmean(values)
    
    Altitude={}
    cstGravit=9.80665 
    footprint = np.array([[0,1,0],
                          [1,0,1],
                          [0,1,0]])

    for i in dicoGeo:
        mask=np.logical_not(dicoGeo[i] > 0).astype(int)
        GeoArray=np.divide(dicoGeo[i],cstGravit)
        np.putmask(GeoArray,mask,np.nan)
        indices = np.where(np.isnan(GeoArray))
        results = ndimage.generic_filter(GeoArray, mean, footprint=footprint)
        for row, col in zip(*indices):
            GeoArray[row,col] = results[row,col]
        Altitude[i]=GeoArray
        
    return Altitude


def computeDailyAccumulation(dicoBand,nbBandByDay,typeData):
    
    accumulation={}
    for i in range(0,len(dicoBand.keys())/nbBandByDay):
        maxRange=nbBandByDay+i*nbBandByDay
            
        #on ne prend pas la dernière bande... correspondante à 00-->3h
        for j in range (i*nbBandByDay,maxRange):
            if "array" in locals():
                array=array+dicoBand.items()[j][1]
            else:
                array=dicoBand.items()[j][1]
        accumulation[i]=array
        del array
    
    return accumulation

def computeDailyMean(dicoBand,nbBandByDay,typeData):

    def meanCalc(values):
        return np.nanmean(values)

    mean={}
    footprint = np.array([[0,1,0],
                          [1,0,1],
                          [0,1,0]])
    
    for i in range(0,len(dicoBand.keys())/nbBandByDay):
        maxRange=nbBandByDay+i*nbBandByDay
        #on ne prend pas la dernière bande... correspondante à 00-->3h
        for j in range (i*nbBandByDay,maxRange):
            if "array" in locals():
                array=array+dicoBand.items()[j][1]
                np.putmask(dicoBand.items()[j][1], dicoBand.items()[j][1]==0, 0)
                mask=mask+(dicoBand.items()[j][1] > 0).astype(int)
            else:
                array=dicoBand.items()[j][1]
                np.putmask(dicoBand.items()[j][1], dicoBand.items()[j][1]==0, 0)
                mask=(dicoBand.items()[j][1] > 0).astype(int)

        mean[i]=array
        del array

        #utilisation de la fonction nanmean --> bcp plus simple

        mean[i]=mean[i]/mask
        indices = np.where(np.isnan(mean[i]))
        results = ndimage.generic_filter(mean[i], meanCalc, footprint=footprint)
        for row, col in zip(*indices):
            mean[i][row,col] = results[row,col]    
    
    return mean

def computeDailyMax(dicoBand,nbBandByDay,typeData=None):
    maxB={}
    for i in range(0,len(dicoBand.keys())/nbBandByDay):
        maxRange=nbBandByDay+i*nbBandByDay
        #on ne prend pas la dernière bande... correspondante à 00-->3h si 
        for j in range (i*nbBandByDay,maxRange):
            if "array" in locals():
                array=np.fmax(array,dicoBand.items()[j][1])
            else:
                array=dicoBand.items()[j][1]
        maxB[i]=array 
        del array
    
    return maxB

def computeDailyMin(dicoBand,nbBandByDay,typeData=None):
    minB={}
    for i in range(0,len(dicoBand.keys())/nbBandByDay):
        maxRange=nbBandByDay+i*nbBandByDay
        #on ne prend pas la dernière bande... correspondante à 00-->3h
        for j in range (i*nbBandByDay,maxRange):
            np.putmask(dicoBand.items()[j][1],dicoBand.items()[j][1]==0,np.nan)
            if "array" in locals():
                array=np.fmin(array,dicoBand.items()[j][1])
            else:
                array=dicoBand.items()[j][1]
        minB[i]=array 
        del array
    
    return minB
            
def fusVentFromDict(dicToFus,nbBandByDay,zmesure=80):
    """ Wind profile relationship [m.s-1]
        Estimate wind speed at 2m
        uz wind speed at height zmesure above ground surface
        
        wind is the norm of U and V direction speed
    """
    wind={}
    keys=dicToFus.keys()
    if (len(dicToFus)==2):
        for i in dicToFus[keys[0]]:
            #Math.log = ln 
            u=dicToFus[keys[0]][i]*4.87/math.log(67.8*zmesure-5.42);
            v=dicToFus[keys[1]][i]*4.87/math.log(67.8*zmesure-5.42);
            wind[i]=np.sqrt(pow(u,2)+pow(v,2))
    
    return wind

def ComputeHumidityFromPT(pressureDico,TDico,DewDico):
    """ Compute Humidity for each Band and each day based on pressure,Temperature and Dew Point"""
    Humidity={}
    
    for i in pressureDico:
        Humidity[i]=esat(pressureDico[i],DewDico[i])/esat(pressureDico[i],TDico[i])*100
        np.putmask(Humidity[i], pressureDico[i]==0, 0)
        np.putmask(Humidity[i], DewDico[i]==0, 0)
        np.putmask(Humidity[i], TDico[i]==0, 0)
    
    return Humidity
    
    
def esat(pressure,T):
    """ Compute partial presure depending on P and T
        P(T)=0.61121*exp(17.502*T/(T+240.97))*(1.00072+pression*(3.2+0.00059*temperature²)/100000.0)
    
        From Wexler and al. 1976
        Pressure en hpa --> convert to kPa
        T en °C
    """
    pressure=pressure/10
    d_es = 0.61121*np.exp(np.multiply(T,17.502)/(T+240.97))
    d_f = 1.00072+pressure*(3.2+0.00059*pow(T,2))/100000.0
                        
    return d_es*d_f

def eocalc(T):
    """ Saturation vapor pressure at the air temperature [KPa]
        T en °C
    """
    eo_calc=0.6108*np.exp(17.27*T/(T+237.3))
    return eo_calc

def delta_calc(T):
    # Slope of saturation vapour pressure curve at air temperature [kPa.°C-1]
    # T air temperature in °C
    # Equation 13 FAO
    delta=4098*(0.6108*np.exp(17.27*T/(T+237.3)))/(T+237.3)**2;
    return delta

def doy(datetoConvert,deltaDays):
    
    deltaJ=timedelta(days=deltaDays)
    datetoConvert=datetoConvert+deltaJ
    J = datetoConvert.timetuple().tm_yday
    
    return J

def getGeoTransform(pathToImg):
        
    srcImage = gdal.Open(pathToImg)
    geoTrans = srcImage.GetGeoTransform()
    
    xOrigin = geoTrans[0]
    yOrigin = geoTrans[3]
    pixelWidth = geoTrans[1]
    pixelHeight = geoTrans[5]

    return (xOrigin,yOrigin,pixelWidth,pixelHeight)
     
def getProj(pathToShape):
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(pathToShape, 0)
    layer = dataSource.GetLayer()
    srs = layer.GetSpatialRef()
    
    return srs.ExportToProj4()    

def getShape(pathToImg):
        
    raster = gdal.Open(pathToImg)
    
    transform = raster.GetGeoTransform()
    pixelWidth = transform[1]
    pixelHeight = transform[5]
    
    return (pixelWidth,pixelHeight)

def getCentroidLatFromArray(shape,geotransform,grid):

    lat = np.zeros(shape)
    lon = np.zeros(shape)

    originX =  geotransform[0]
    originY =  geotransform[1]
    
    for index in np.ndenumerate(lat):
        lat.itemset(index[0], float(originX)+float(index[0][1])*float(grid)+(float(grid)/2))
        lon.itemset(index[0], float(originY)-float(index[0][0])*float(grid)-(float(grid)/2))

    
    dicoLatLong={}
    dicoLatLong[0]=lat
    dicoLatLong[1]=lon
    
    return dicoLatLong

<<<<<<< HEAD
def WriteTxtFileForEachPixel(outputFolder,et0_0,DateList,DoyList,RayShort,Tmean,Tmax,Tmin,Hmean,Hmax,Hmin,vent,precipitation,pressure,Geo,latlon,projShape):
    """ Write a Txtfile """
    
    for i in range(0,et0_0[0].shape[0]): #row
        for j in range(0,et0_0[0].shape[1]): #col
=======
def writeTiffFromDicoArray(DicoArray,outputImg,shape,geoparam,proj=None,format=gdal.GDT_Float32):
    
    gdalFormat = 'GTiff'
    driver = gdal.GetDriverByName(gdalFormat)

    dst_ds = driver.Create(outputImg, shape[1], shape[0], len(DicoArray), format)
    
    j=1
    for i in DicoArray.values():
        dst_ds.GetRasterBand(j).WriteArray(i, 0)
        band = dst_ds.GetRasterBand(j)
        band.SetNoDataValue(0)
        j+=1
    
    originX =  geoparam[0]
    originY =  geoparam[1]
    pixelWidth = geoparam[2]
    pixelHeight  = geoparam[3]

    
    dst_ds.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))

def WriteTxtFileForEachPixel(outputFolder,et0_0,et0_1,et0_2,DateList,DoyList,Ray,RayShort,RayLong,Tmean,Tmax,Tmin,Hmean,Hmax,Hmin,vent,precipitation,pressure,Geo,latlon):
    """ Write a Txtfile """
    
    
    for i in range(0,et0_0[0].shape[0]):
        for j in range(0,et0_0[0].shape[1]):
>>>>>>> 60ca0156007a438a81c6e95150f141bf4f4d38e9
            lat=latlon[0][i][j]
            lon=latlon[1][i][j]
            p1 = pp.Proj(projShape)
            latP,lonP = p1(lat,lon)

            numero = str(round(lat,2)).replace('.','')+str(round(lon,2)).replace('.','')
            
            pathTodateFolder=outputFolder+'/POINT_'+numero+'.txt'
            f = open(pathTodateFolder,'w+')
<<<<<<< HEAD
            f.write('numero;altitude;lat/lon(WGS84);lat/lon(initial)\n')
            f.write(str(numero)+'\t ; '+str(Geo[0][i][j])+'\t ; '+str(lat)+'/'+str(lon)+';'+str(latP)+'/'+str(lonP)+'\n')
            f.write('ANNEE\tMOIS\tJOUR\tDOY\tRGSHORT\tTAMEAN\tTAMAX\tTAMIN\tRHMEAN\tRHMAX\tRHMIN\tVUMEAN\tPRECIP\tPRESSURE\tET0FAO56\n')
            f.write('[YYYY]\t[MM]\t[DD]\t[1-365]\t[MJ.m-2.jour-1]\t[Kelvin]\t[Kelvin]\t[Kelvin]\t[%]\t[%]\t[%]\t[m.s-1]\t[kPa]\t[mm.d-1]\t[mm.d-1]\n')
=======
            f.write('numero;altitude;lat/lon(WGS84)\n')
            f.write(str(numero)+'\t ; '+str(Geo[0][i][j])+'\t ; '+str(lat)+'/'+str(lon)+'\n')
            f.write('ANNEE\tMOIS\tJOUR\tDOY\tRGSURF\tRGLONG\tRGSHORT\tTAMEAN\tTAMAX\tTAMIN\tRHMEAN\tRHMAX\tRHMIN\tVUMEAN\tPRECIP\tPRESSURE\tET0FAO56\tET0SolarEra\tEvapEraInterim\n')
            f.write('[YYYY]\t[MM]\t[DD]\t[1-365]\t[MJ.m-2.jour-1]\t[MJ.m-2.jour-1]\t[MJ.m-2.jour-1]\t[Kelvin]\t[Kelvin]\t[Kelvin]\t[%]\t[%]\t[%]\t[m.s-1]\t[kPa]\t[mm.d-1]\t[mm.d-1]\t[mm.d-1]\t[mm.d-1]\n')
>>>>>>> 60ca0156007a438a81c6e95150f141bf4f4d38e9
            for d in range(0,len(DateList)):
                year=DateList[d].year
                month=DateList[d].month
                day=DateList[d].day
<<<<<<< HEAD
                f.write(str(year)+'\t'+str(month)+'\t'+str(day)+'\t'+ str(DoyList[d])+'\t'+str(RayShort[d][i][j])+'\t'+str(Tmean[d][i][j])+'\t'+str(Tmax[d][i][j])+'\t'+str(Tmin[d][i][j])+'\t'+str(Hmean[d][i][j])+'\t'+str(Hmax[d][i][j])+'\t'+str(Hmin[d][i][j])+'\t'+ str(vent[d][i][j])+'\t'+str(precipitation[d][i][j])+'\t'+str(pressure[d][i][j])+'\t'+str(et0_0[d][i][j])+'\n')
=======
                f.write(str(year)+'\t'+str(month)+'\t'+str(day)+'\t'+ str(DoyList[d])+'\t'+str(Ray[d][i][j])+'\t'+str(RayShort[d][i][j])+'\t'+str(RayLong[d][i][j])+'\t'+str(Tmean[d][i][j])+'\t'+str(Tmax[d][i][j])+'\t'+str(Tmin[d][i][j])+'\t'+str(Hmean[d][i][j])+'\t'+str(Hmax[d][i][j])+'\t'+str(Hmin[d][i][j])+'\t'+ str(vent[d][i][j])+'\t'+str(precipitation[d][i][j])+'\t'+str(pressure[d][i][j])+'\t'+str(et0_0[d][i][j])+'\t'+str(et0_1[d][i][j])+'\t'+str(et0_2[d][i][j])+'\n')
>>>>>>> 60ca0156007a438a81c6e95150f141bf4f4d38e9
            f.close()
    return pathTodateFolder
    
def WritePointList(outputFolder,latlon,projShape):
    
    pathTodateFolder=outputFolder+'/ListeStations.txt'
    f = open(pathTodateFolder,'w+')
    f.write('numero;lat/lon(WGS84)\n')
    for i in range(0,latlon[0].shape[0]):
        for j in range(0,latlon[0].shape[1]):
            lat=latlon[0][i][j]
            lon=latlon[1][i][j]
            p1 = pp.Proj(projShape)
            latP,lonP = p1(lat,lon)
            
            numero = str(round(lat,2)).replace('.','')+str(round(lon,2)).replace('.','')
            f.write(str(numero)+';'+str(lat)+';'+str(lon)+';'+str(latP)+';'+str(lonP)+'\n')
    f.close()


