#-*- coding: utf-8 -*-
'''
Created on 25 mars. 2014

Toolbox for ET0 calculation based on EraInterim product
for real time purpose

@author: yoann Moreau
@author: vincent Rivallant
'''

import sys
import getopt
import os
import numpy as np
import math

from datetime import date,datetime,timedelta
import utils as utils


def main(argv):

    try:
        opts,argv = getopt.getopt(argv,":h:i:e:s:E:o:g:P:t:f:r:",['help','[outFile]','code','[shapeFile]','start','end','[tr]'])
    except getopt.GetoptError:
        print 'error in parameter for evapo_GFS. type evapo_GFS.py -help for more detail on use '
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            print 'evapo_GFS.py  '
            print '    [mandatory] : '
            print '        --init <dateStart YYYY-MM-DD>'
            print '        --end <dateEnd YY-MM-DD>'
            print '        --shapefile <shapefile> OU -Extend < xmin,ymax,xmax,ymin>'
            print '    [optional] :'
            print '        --typeData  < analyse , forecast, cycleforecast> (default forecast)'
            print '        --grid <EraInterim Time> (default 0.75)'
            print '        --outfile <outfolder> (default /home/user/eraInterim)'
            print '        --proxy <True/False> (default False)'
            print '        --temporaryFile <True/False> (default False)'
            print '        --result < TxtFile / RasterFile> default RasterFile'
            print ''
            sys.exit() 
        elif opt in ('-o','--outFolder'):
            oFolder = arg
        elif opt in ('-i','--start'):
            startDate = arg
        elif opt in ('-e','--end'):
            endDate = arg
        elif opt in ('-s','--shapefile'):
            pathToShapefile = arg
        elif opt in ('-E','--tr'):
            extend = arg.split(',')
        elif opt in ('-g','--grid'):
            grid = arg
        elif opt in ('-P','--proxy'):
            proxy = arg
        elif opt in ('-t','--typeData'):
            typeData = arg
        elif opt in ('-f','--temporaryFile'):
            temporaryFile = arg
        elif opt in ('-r','--result'):
            typeOutput = arg
    
    if len(sys.argv) < 7:
        print 'evapo_GFS.py'
        print '    -i <dateStart YYYY-MM-DD> '
        print '    -e <dateEnd YY-MM-DD>'
        print '    -s <shapefile> '
        print '  or'
        print '    -E < xmin,ymax,xmax,ymin>]'
        print ''
        print '    [-t < analyse , forecast,cycleforecast> (default analyse)]'
        print '    [-g <size of grid in 0.25/0.5/1/2.5> (default 0.25)]'
        print '    [-o <outfolder> (default /home/user/eraInterim)]'
        print '    [-P <proxy : True/False> (default False)]'
        print '    [-f <temporaryFile : True/False> (default False)]'
        print '    [-r <resultOutput : TxtFile/RasterFile> (default RasterFile)]'
        print ''
        sys.exit(2)
        
    try:
        oFolder
    except NameError:
        oFolder = os.path.expanduser('~')
        oFolder = oFolder + '/GFS'
        print "output folder not precised : downloaded GFS images on "+oFolder
    
    # verification du folder/or creation if not exists
    utils.checkForFolder(oFolder) 
        
    try:
        startDate
    except NameError:
        exit ('init Date not precised')
    # verification si sartDate est une date
    startDate=utils.checkForDate(startDate) 
    
    try:
        endDate
    except NameError:
        exit ('end Date not specified')
    # verification si sartDate est une date
    endDate=utils.checkForDate(endDate) 
    
    if (startDate>endDate):
        exit('startDate could not be greater than endDate')

    try:
        pathToShapefile
    except NameError:
        try:
            extend
        except NameError:
            exit ('no Area of interest have been specified. please use -shp or -tr to declare it')
    
    if 'pathToShapefile' in locals():
        extendArea=utils.convertShpToExtend(pathToShapefile)
    else:
        extendArea=extend

    extendArea=utils.checkForExtendValidity(extendArea)
      
    try:
        typeData
    except NameError:
        typeData='analyse'

    try:
        grid
    except NameError:
        grid='0.25'
    grid=utils.checkForGridValidity(grid)
            
    try:
        proxy
    except NameError:
        proxy=False

    #Proxy parameteres needed
    if(proxy):
        login = raw_input('login proxy : ')
        pwd = raw_input('password proxy :  : ')
        site = raw_input('site (surf.cnes.fr) : ')
        os.environ["http_proxy"] = "http://%s:%s@%s:8050"%(login,pwd,site)
        os.environ["https_proxy"] = "http://%s:%s@%s:8050"%(login,pwd,site)
    
    try:
        temporaryFile
    except NameError:
        temporaryFile=False
    
    try:
        typeOutput
    except NameError:
        typeOutput='RasterFile'
    """----------------------------------------"""
        
    delta = endDate - startDate
    nbDays = delta.days + float(delta.seconds) / 86400 + 1

    #--------------------------On charge les rasters
    step = [0,6,12,18]
    levelList = ["surface"]
    nbBandByDay=len(step)
    
    """ altitude de la grille GFS """
    print('downloading altitude grid')
    # Only Forecast possible
    codeGeopot= ['HGT']
    struct=utils.create_request_gfs(startDate, endDate, step, levelList, grid, extendArea, codeGeopot, typeData)        
    listeGeoFile=[]
    if len(struct[0])==0:
        exit("No data founded")
    else:
        for i in struct[0]:
            try :
                GeoFile = oFolder+'/'+",".join(codeGeopot)+'_'+i.rsplit('.',1)[1]+'.grb'
                listeGeoFile.append(GeoFile)
                result=utils.GFSDownload(i,GeoFile)
            except:
                print("---")
                exit('Error in GFS server')
        
        if result:
            Geo = utils.convertGribToDicoArray(listeGeoFile, codeGeopot, levelList, step, grid, startDate, endDate)
            Geo = utils.convertGeoToAlt(Geo)
            #un peu inutile car ne change pas ... mais bon! 
            Geo = utils.computeDailyMax(Geo, nbBandByDay, typeData)
        else:
            exit("PARAM needed is not compatible with level selected")
    
    """ Vitesse du vent """
    print('downloading wind velocity')
    codeVent= ['UGRD','VGRD']
    levelListVent = ["80_m_above_ground"]
    vent={}
    ventFile=[]
    for v in codeVent:
        struct=utils.create_request_gfs(startDate, endDate, step, levelListVent, grid, extendArea, codeVent, typeData)        
        listeFile=[]
        if len(struct[0])==0:
            exit("No data founded")
        else:
            for i in struct[0]:
                try :
                    tmpVent = oFolder+'/'+v+'_'+i.rsplit('.',1)[1]+'.grb'
                    ventFile.append(tmpVent)
                    result=utils.GFSDownload(i,tmpVent)
                except:
                    print("---")
                    exit('Error in GFS server')
            
            if result:
                vent[v] = utils.convertGribToDicoArray(ventFile, codeVent, levelList, step, grid, startDate, endDate)
            else:
                exit("PARAM needed is not compatible with level selected")
    
    vent = utils.fusVentFromDict(vent,nbBandByDay)
    vent=utils.computeDailyMean(vent,nbBandByDay,typeData)
    exit
    """ Humidité relative """
    print('downloading relative humidity')
    codePressure= ['PRES']
    struct=utils.create_request_gfs(startDate, endDate, step, levelList, grid, extendArea, codePressure, typeData)        
    listePressureFile=[]
    if len(struct[0])==0:
        exit("No data founded")
    else:
        for i in struct[0]:
            try :
                pressureFile = oFolder+'/'+",".join(codePressure)+'_'+i.rsplit('.',1)[1]+'.grb'
                listePressureFile.append(pressureFile)
                result=utils.GFSDownload(i,pressureFile)
            except:
                print("---")
                exit('Error in GFS server')
        
        if result:
            pressure = utils.convertGribToDicoArray(listePressureFile, codePressure, levelList, step, grid, startDate, endDate)
            pressureMean = utils.convertPaToKgPa(pressure)
            pressure = utils.convertToHectoPascal(pressure)
            pressureMean = utils.computeDailyMean(pressureMean,nbBandByDay,typeData)
        else:
            exit("PARAM needed is not compatible with level selected")
    
    """ Temperature 2m """
    print('downloading Temperature')
    levelListTmp = ["2_m_above_ground"]
    codeT2m= ['TMP']
    struct=utils.create_request_gfs(startDate, endDate, step, levelListTmp, grid, extendArea, codeT2m, typeData)        
    listeTmpFile=[]
    if len(struct[0])==0:
        exit("No data founded")
    else:
        for i in struct[0]:
            try :
                T2mFile = oFolder+'/'+",".join(codeT2m)+'_'+i.rsplit('.',1)[1]+'.grb'
                listeTmpFile.append(T2mFile)
                result=utils.GFSDownload(i,T2mFile)
            except:
                print("---")
                exit('Error in GFS server')
        
        if result:
            T2m = utils.convertGribToDicoArray(listeTmpFile, codeT2m, levelList, step, grid, startDate, endDate)
            Tmean = utils.computeDailyMean(T2m, nbBandByDay, typeData)
            Tmax = utils.computeDailyMax(T2m,nbBandByDay)
            Tmin = utils.computeDailyMin(T2m,nbBandByDay)
            T2m = utils.convertKToD(T2m)
        else:
            exit("PARAM needed is not compatible with level selected")
    
    """ DewPoint """
    print('downloading DewPoint')
    levelListdew = ["2_m_above_ground"]
    codeDewP= ['DPT']
    struct=utils.create_request_gfs(startDate, endDate, step, levelListdew, grid, extendArea, codeDewP, typeData)        
    listeDewFile=[]
    if len(struct[0])==0:
        exit("No data founded")
    else:
        for i in struct[0]:
            try :
                DewPFile = oFolder+'/'+",".join(codeDewP)+'_'+i.rsplit('.',1)[1]+'.grb'
                listeDewFile.append(DewPFile)
                result=utils.GFSDownload(i,DewPFile)
            except:
                print("---")
                exit('Error in GFS server')
        
        if result:
            DewP = utils.convertGribToDicoArray(listeDewFile, codeDewP, levelList, step, grid, startDate, endDate)
            DewP = utils.convertKToD(DewP)
        else:
            exit("PARAM needed is not compatible with level selected")

    
    humidity = utils.ComputeHumidityFromPT(pressure,T2m,DewP)
    #humidity = utils.computeDailyMean(humidity,nbBandByDay,typeData)
    Hmax = utils.computeDailyMax(humidity,nbBandByDay)
    Hmin = utils.computeDailyMin(humidity,nbBandByDay)
    Hmean = utils.computeDailyMean(humidity,nbBandByDay,typeData)
    

    """ Rayonnement global incident journalier """
    print('downloading downscale shortwave radiation')
    # Only Forecast possible
    levelListRay = ["surface"]
    codeRay= ["DSWRF"]
    typeData="cycleforecast"
    struct=utils.create_request_gfs(startDate, endDate, step, levelListRay, grid, extendArea, codeRay, typeData)        
    listeRayFile=[]
    if len(struct[0])==0:
        exit("No data founded")
    else:
        for i in struct[0]:
            try :
                RayFileDownShort = oFolder+'/'+",".join(codeRay)+'_'+i.rsplit('.',1)[1]+'.grb'
                listeRayFile.append(RayFileDownShort)
                result=utils.GFSDownload(i,RayFileDownShort)
            except:
                print("---")
                exit('Error in GFS server')
        
        if result:
            RayDownShort = utils.convertGribToDicoArray(listeRayFile, codeRay, levelList, step, grid, startDate, endDate)
            RayDownShort = utils.computeDailyMean(RayDownShort,nbBandByDay,typeData)
            RayDownShort = utils.convertWToMJ(RayDownShort)
        else:
            exit("PARAM needed is not compatible with level selected")
    
    """ Precipitation """
    #NOT NEEDED FOR ETO BUT Exported
    print('downloading Precipitation')
    codePrecipitation= ["APCP"]
    struct=utils.create_request_gfs(startDate, endDate, step, levelList, grid, extendArea, codePrecipitation, typeData)        
    listePrecFile=[]
    if len(struct[0])==0:
        exit("No data founded")
    else:
        for i in struct[0]:
            try :
                precipitationFile = oFolder+'/'+",".join(codePrecipitation)+'_'+i.rsplit('.',1)[1]+'.grb'
                listePrecFile.append(precipitationFile)
                result=utils.GFSDownload(i,precipitationFile)
            except:
                print("---")
                exit('Error in GFS server')
        
        if result:
            precipitation = utils.convertGribToDicoArray(listePrecFile, codePrecipitation, levelList, step, grid, startDate, endDate)
            precipitation=utils.computeDailyAccumulation(precipitation,nbBandByDay,typeData)
        else:
            exit("PARAM needed is not compatible with level selected")
    
    """ Grid of latitude [0],longitude[1] in WGS84"""
    geoTransform=utils.getGeoTransform(RayFileDownShort)
    shape=RayDownShort[0].shape
    latlon = utils.getCentroidLatFromArray(shape,geoTransform,grid)

    """ --------------------- Compute ET0---------------------- """
    
    ET0_0={}
    DoyList=[]
    DateList=[]
    
    for i in range(0,int(nbDays)):
        #jour Julien
        J = utils.doy(startDate,i)
        dateEnCours=startDate+ timedelta(days=i)
        DateList.append(dateEnCours)
        DoyList.append(J)
        Hmax[i] = np.where(Hmax[i]>100,100,Hmax[i])
        
        # --- Constants ---# 
        #Solar constant
        Gsc = 0.0820 # [MJ.m-2.min-1]
        #Albedo - grass reference crop
        a = 0.23
        #Ratio of molecular weight of water vapor/dry air
        epsilon=0.622 
        #Latente heat of vaporisation
        Lv=2.45 # [MJ.kg-1]
        # Specific heat at constant pressure
        Cp= 1.013e-3; # [MJ.kg-1.°C-1]
        # Stefan-Boltzmann constant [MJ.K-4.m-2.day-1]
        StefBoltz=4.903e-9; #FAO
        
        # --- Equations ---# 
        # Psychometric constant [kPa.°C-1]
        cte_psy = (Cp*pressureMean[i])/(epsilon*Lv) # Equation 8 Chap 3 FAO
        #Mean sturation vapor presure [kPa]
        #es = (utils.esat(pressureMean[i],Tmax[i]) + utils.esat(pressureMean[i],Tmin[i]))/2;    #Equation 12 Chap 3
        es = (utils.eocalc(Tmax[i]-273.16)+utils.eocalc(Tmin[i]-273.16))/2    #Equation 12 Chap 3
        # Slope of saturation vapour pressure curve at air temperature [kPa.°C-1]
        delta = utils.delta_calc(Tmean[i]);                    # Equation 13 Chap 3
        # Actual vapour pressure derived from relative humidity [kPa]
        #ea = (utils.esat(pressureMean[i]/100,Tmax[i]-273.16)*(Hmax[i]/100) + utils.esat(pressureMean[i]/100,Tmin[i]-273.16)*(Hmin[i]/100))/2;      # Equation 17 Chap 3
        ea = (utils.eocalc(Tmax[i]-273.16)*(Hmax[i]/100)+utils.eocalc(Tmin[i]-273.16)*(Hmin[i]/100))/2
        # Conversion of latitude from degrees to radians
        phi = (np.pi/180)* latlon[1];     
        # Relative distance Earth-Sun
        dr = 1 + 0.033*math.cos(2*math.pi*J/365);         # Equation 23 Chap 3
        # Solar declination
        d = 0.4093*math.sin(2*math.pi*J/365 - 1.39);      # Equation 24 Chap 3
        # sunset hour angle
        ws = np.arccos(-np.tan(phi)*math.tan(d));                # Equation 25 Chap 3
        
        """Classical calculation FAO """
        
        # Extraterestrial radiation for daily periods
        Ra = (24.*60/np.pi)*Gsc*dr*(ws*np.sin(phi)*np.sin(d) + np.cos(phi)*np.cos(d)*np.sin(ws))    # Equation 21 Chap 3
        # Clear sky solar radiation [MJ.m-2.day-1]
        Rso = (0.75 + 2e-5*Geo[i])*Ra;                 # Equation 37 Chap 3
        # Net solar radiation [MJ.m-2.day-1]
        Rns = (1 - a)*RayDownShort[i];                          # Equation 38 Chap 3
        #
        f=(1.35*(np.fmin(RayDownShort[i]/Rso,1)) - 0.35);
        # Net longwave radiation [MJ.m-2.day-1]
        Rnl = StefBoltz*((Tmax[i]**4 + Tmin[i]**4)/2)*(0.34 - 0.14*np.sqrt(ea))*f; # Equation 39 Chap 3
        # Net Radiation [MJ.m-2.day-1]
        Rn =Rns - Rnl;                              # Equation 40 Chap 3
        G = 0;                                      # Equation 42 Chap 3
        ET0_0[i] = ( 0.408*delta*( Rn-G )+ cte_psy*( 900/(Tmean[i] + 273) )*(es - ea)*vent[i] )/( delta + cte_psy*(1 + 0.34*vent[i]) );  # Equation 6 Chap 4
    
    if typeOutput=='RasterFile':
        #On ecrit le fichier ET0 
        geoTransform=utils.getGeoTransform(RayFileDownShort)
        shape=RayDownShort[0].shape
        utils.writeTiffFromDicoArray(ET0_0,oFolder+"/tmp.tif",shape,geoTransform)
        if 'pathToShapefile' in locals():
            utils.reprojRaster(oFolder+"/tmp.tif", oFolder+"/ET0.tif",shape, pathToShapefile)
            os.remove(oFolder+"/tmp.tif")
        else : 
            utils.moveFile(oFolder+"/tmp.tif", oFolder+"/ET0.tif")

        #On écrit le fichier Precipitation
        geoTransform=utils.getGeoTransform(precipitationFile)
        shape=precipitation[0].shape
        utils.writeTiffFromDicoArray(precipitation,oFolder+"/tmp.tif",shape,geoTransform)
    
        if 'pathToShapefile' in locals():
            utils.reprojRaster(oFolder+"/tmp.tif", oFolder+"/precipitationAcc.tif",shape, pathToShapefile)
            os.remove(oFolder+"/tmp.tif")
        else : 
            utils.moveFile(oFolder+"/tmp.tif", oFolder+"/precipitationAcc.tif")
        
    else:
        #On ecrit le fichier au format Txt
        proj=utils.getProj(pathToShapefile)
        utils.WriteTxtFileForEachPixel(oFolder,ET0_0,DateList,DoyList,RayDownShort,Tmean,Tmax,Tmin,Hmean,Hmax,Hmin,vent,precipitation,pressureMean,Geo,latlon,proj)
        utils.WritePointList(oFolder,latlon,proj)
    
    """ ------------------------------------------- """
    if(temporaryFile):
        #On ecrit le fichier latlon 
        geoTransform=utils.getGeoTransform(GeoFile)
        shape=Geo[0].shape
        utils.writeTiffFromDicoArray(Geo,oFolder+"/tmp.tif",shape,geoTransform)
    
        if 'pathToShapefile' in locals():
            utils.reprojRaster(oFolder+"/tmp.tif", oFolder+"/altitude.tif",shape, pathToShapefile)
            os.remove(oFolder+"/tmp.tif")
        else : 
            utils.moveFile(oFolder+"/tmp.tif", oFolder+"/altitude.tif")
    
        #On ecrit le fichier latlon 
        geoTransform=utils.getGeoTransform(RayFileDownShort)
        shape=RayDownShort[0].shape
        utils.writeTiffFromDicoArray(latlon,oFolder+"/tmp.tif",shape,geoTransform)
    
        if 'pathToShapefile' in locals():
            utils.reprojRaster(oFolder+"/tmp.tif", oFolder+"/latLon.tif",shape, pathToShapefile)
            os.remove(oFolder+"/tmp.tif")
        else : 
            utils.moveFile(oFolder+"/tmp.tif", oFolder+"/latLon.tif")
        
        
        #On ecrit le fichier vent --> a enlever
        geoTransform=utils.getGeoTransform(ventFile[-1])
        shape=vent[0].shape
        utils.writeTiffFromDicoArray(vent,oFolder+"/tmp.tif",shape,geoTransform)
    
        if 'pathToShapefile' in locals():
            utils.reprojRaster(oFolder+"/tmp.tif", oFolder+"/ventMean.tif",shape, pathToShapefile)
            os.remove(oFolder+"/tmp.tif")
        else : 
            utils.moveFile(oFolder+"/tmp.tif", oFolder+"/ventMean.tif")
        
        #On ecrit le fichier Rhmin
        geoTransform=utils.getGeoTransform(pressureFile)
        shape=pressureMean[0].shape
        utils.writeTiffFromDicoArray(pressureMean,oFolder+"/tmp.tif",shape,geoTransform)
    
        if 'pathToShapefile' in locals():
            utils.reprojRaster(oFolder+"/tmp.tif", oFolder+"/pressureMean.tif",shape, pathToShapefile)
            os.remove(oFolder+"/tmp.tif")
        else : 
            utils.moveFile(oFolder+"/tmp.tif", oFolder+"/pressureMean.tif")
        
        #On ecrit le fichier Rhmax
        geoTransform=utils.getGeoTransform(pressureFile)
        shape=Hmax[0].shape
        utils.writeTiffFromDicoArray(Hmax,oFolder+"/tmp.tif",shape,geoTransform)
    
        if 'pathToShapefile' in locals():
            utils.reprojRaster(oFolder+"/tmp.tif", oFolder+"/humidityMax.tif",shape, pathToShapefile)
            os.remove(oFolder+"/tmp.tif")
        else : 
            utils.moveFile(oFolder+"/tmp.tif", oFolder+"/humidityMax.tif")
        
        #On ecrit le fichier Rhmin
        geoTransform=utils.getGeoTransform(pressureFile)
        shape=Hmin[0].shape
        utils.writeTiffFromDicoArray(Hmin,oFolder+"/tmp.tif",shape,geoTransform)
    
        if 'pathToShapefile' in locals():
            utils.reprojRaster(oFolder+"/tmp.tif", oFolder+"/humidityMin.tif",shape, pathToShapefile)
            os.remove(oFolder+"/tmp.tif")
        else : 
            utils.moveFile(oFolder+"/tmp.tif", oFolder+"/humidityMin.tif")
    
        #On ecrit le fichier Tmax
        geoTransform=utils.getGeoTransform(T2mFile)
        shape=Tmax[0].shape
        utils.writeTiffFromDicoArray(Tmax,oFolder+"/tmp.tif",shape,geoTransform)
    
        if 'pathToShapefile' in locals():
            utils.reprojRaster(oFolder+"/tmp.tif", oFolder+"/TemperatureMax.tif",shape, pathToShapefile)
            os.remove(oFolder+"/tmp.tif")
        else : 
            utils.moveFile(oFolder+"/tmp.tif", oFolder+"/TemperatureMax.tif")
    
        #On ecrit le fichier Tmin
        geoTransform=utils.getGeoTransform(T2mFile)
        shape=Tmin[0].shape
        utils.writeTiffFromDicoArray(Tmin,oFolder+"/tmp.tif",shape,geoTransform)
    
        if 'pathToShapefile' in locals():
            utils.reprojRaster(oFolder+"/tmp.tif", oFolder+"/TemperatureMin.tif",shape, pathToShapefile)
            os.remove(oFolder+"/tmp.tif")
        else : 
            utils.moveFile(oFolder+"/tmp.tif", oFolder+"/TemperatureMin.tif")
        
        #On ecrit le fichier Rayonnement
        geoTransform=utils.getGeoTransform(RayFileDownShort)
        shape=RayDownShort[0].shape
        utils.writeTiffFromDicoArray(RayDownShort,oFolder+"/tmp.tif",shape,geoTransform)
    
        if 'pathToShapefile' in locals():
            utils.reprojRaster(oFolder+"/tmp.tif", oFolder+"/RayonnementMean.tif",shape, pathToShapefile)
            os.remove(oFolder+"/tmp.tif")
        else : 
            utils.moveFile(oFolder+"/tmp.tif", oFolder+"/RayonnementMean.tif")
    
    #on supprime les fichier intermédiare !
    for i in listePressureFile:
        os.remove(i)
    for i in listeTmpFile:
        os.remove(i)
    for i in listeDewFile:
        os.remove(i)
    for i in listeRayFile:
        os.remove(i)
    for i in listeGeoFile:
        os.remove(i)
    for i in ventFile:
        os.remove(i)
    for i in listePrecFile:
        os.remove(i)
    
if __name__ == '__main__':
    main(sys.argv[1:])

    pass