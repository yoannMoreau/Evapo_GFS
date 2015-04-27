=======
# Evapo_EraInterim
Tool to calculate refence Evapo-transpiration on an area based on Era Interim models
The calculation is based on FAO equation. 


<h2>Installation<b></h2>

mkdir PATH/TO/INSTALL <br>
cd  PATH/TO/INSTALL/Evapo_EraInterim <br>
git clone git+git://github.com/yoannMoreau/Evapo_EraInterim.git <br>
sudo pip install https://software.ecmwf.int/wiki/download/attachments/23694554/ecmwf-api-client-python.tgz 
python python/calculateETfromERAI.py.py -help <br>

( or you cas Use pip to install Evapo_EraInterim )

<h2>Overview: What can Evapo_EraInterim do?</h2>

Evapo_EraInterim has a main function, allow download parameters for evapo-transpiration : <br>
Temperature, <br>
Wind speed, <br>
Net Solar radiation <br>
Pressure <br>
Precipitation <br>
<br>
All these data are used to compute Et0 wherever the user need it. The result could be export on tif or in txt file
for precise analysis. 
It's possible to compute Et0 over large period of time, the final result will have as band as years needed.

<u>Three paramaters are mandatory: <br><br></u>
<b>--Interval needed : </b><br>
init date <dateStart YYYY-MM-DD>' AND end date <dateEnd YY-MM-DD>'
init date and end date coul not exceed, at the time of writing 2014-31-12
<br><br> 
<b>--Area needed </b><br>
shapefile <pathToShapefile> (srs is not important because it will be reprojected in WGS84)
OR 
--Extend <xmin,ymax,xmax,ymin> in WGS84
<br><br>

<b>EXAMPLES :</b><br>
<i>EvapoTranspiration during january 2014 <br></i>
python calculateETfromERAI.py -i 2014-01-01 -e 2014-01-31 -s PATH_TO_SHAPE'<br>
<br><br>
<u>Five paramaters are optional: </u><br><br>
<b>----typeData <Type of Data : analyse or forcast >  (default analyse)'</b><br>
The type of data needed. 
Analyse is post production data. At time of writing it's possible to use until the end of 2014
Forecast are less precise but operationnal data. 
Forecast use step of 3 hour during all the day (start a 00h until midnight). 
The last part (from 21h00 to 00h) is not used due to some error on Era-Interim modelisaton
Analyse use 4 time : 00h00 06h00 12h00 and 18h00
<br><br>
python calculateETfromERAI.py -i 2013-11-08 -e 2013-12-09 -E xmin,ymax,xmax,ymin -t analyse
<br><br>
<b>--grid <grid size> (default 0.75)' </b><br>
the grid of output data. It should be from 0.125 to 3 (almost  12km to 300km pixel width) But the model 
run for an output of 80km, the data are interpolasized for finest grid.
default is 0.75
<br><br>
python calculateETfromERAI.py  -i 2013-11-08 -e 2013-12-09 -E xmin,ymax,xmax,ymin -g 0,125'
<br><br>
<b>--calculateETfromERAI <Path to downloaded EvapoTranspiration> (default /home/user/eraInterim)'</b>
<br><br>
python calculateETfromERAI.py -i 2011-10-01 -e 2011-10-02 -s PATH/TO/SHAPEFILE -o PATH/TO/FILE'
All downloaded raster are TIF 
<br><br>
<b>--proxy <proxy : True/False></b> (default False)
<br><br>
Sometimes a proxy definition is needed for downloading from external network.
When this option is activated, a proxy user/key/site could be defined to overpass it
<br><br>
<b>--temporaryFile <True/False></b> (default False)
<br><br>
It could be useful to have intermediary product to check how computation have been done.
default is False :
Temporary files are : 
altitude / humidityMax / humidityMin / latLon / precipitationAcc / pressureMean / RayonnementMean / TemperatureMax / TemperatureMin / ventMean
<br><br>
python calculateETfromERAI.py -i 2014-01-01 -e 2014-01-04 -s /home/ouaf/landsat/zone_etude.shp -f True
<br><br>
<b>--resultOutput <TxtFile/RasterFile></b> (default False)
<br><br>
For some reason it could be useful to have txtfiles for each pixel. it could be done using -r TxtFile option.
Each pixel will be export in a txt format.
Defaut is RasterFile.
<br><br>
python calculateETfromERAI.py -i 2014-01-01 -e 2014-01-04 -s /home/ouaf/landsat/zone_etude.shp -r TxtFile


<h2>Important Notes </h2>

All downloaded and processed images are stored by default in your home directory in eraInterim forlder: ~/eraInterim
<br><br>
To Do List
<br><br>
Improve console output<br>
Maintain with bug error <br>

