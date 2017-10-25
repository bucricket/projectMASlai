#!/usr/bin/env python2
#!/Applications/anaconda/envs/root3 python
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 13:50:04 2016

@author: mschull
"""
#python

#from .search import Search
import os
import numpy as np
import subprocess
import glob
import pandas as pd
import datetime
import argparse
import getpass
import keyring
from pyproj import Proj
from .utils import folders
from getlandsatdata import getlandsatdata
import pycurl
from .landsatTools import landsat_metadata
import logging
from pymodis.downmodis import downModis
logging.getLogger("urllib3").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.DEBUG)

# The current URL hosting the ESPA interfaces has reached a stable version 1.0
host = 'https://espa.cr.usgs.gov/api/v1/'
TIMEOUT=86400

base = os.getcwd()
cacheDir = os.path.abspath(os.path.join(base,os.pardir,"SATELLITE_DATA","LANDSAT"))
Folders = folders(base)   
modis_base = Folders['modis_base']
landsat_SR = Folders['landsat_SR']
landsat_LAI = Folders['landsat_LAI']
landsat_NDVI = Folders['landsat_NDVI']
landsat_Mask = Folders['landsat_Mask']
landsat_temp = os.path.join(landsat_SR,'temp')
if not os.path.exists(landsat_temp):
    os.mkdir(landsat_temp)



def get_modis_lai(tiles,product,version,start_date,end_date,auth):  
    startdd = datetime.datetime.strptime(start_date, '%Y-%m-%d')
    enddd = datetime.datetime.strptime(end_date, '%Y-%m-%d')
    numDays= (enddd-startdd).days
    
    laidates = np.array(range(1,366,4))

    if product.startswith('MCD'):
        folder = "MOTA"
    elif product.startswith('MOD'):
        folder = "MOLT"
    else:
        folder = "MOTA"
    product_path = os.path.join(modis_base,product)   
    if not os.path.exists(product_path):
        os.mkdir(product_path)
        
    modisOgg = downModis(url="https://e4ftl01.cr.usgs.gov", destinationFolder=product_path, 
                         user=auth[0], password=auth[1], tiles=tiles, path=folder, 
                         product="%s.%s" % (product,version),today=start_date,enddate=end_date)

    # connect to http or ftp
    modisOgg.connect()
    if modisOgg.nconnection <= 20:
        # download data
        for i in range(numDays+1):
            dd = startdd+datetime.timedelta(days=i)            
            year = dd.year
            doy = (dd-datetime.datetime(year,1,1,1,1)).days
            rday = laidates[laidates>=doy][0]
            dd = datetime.datetime(year,1,1,1,1)+datetime.timedelta(days=rday-1)
            year = dd.year
            month = dd.month
            day = dd.day
            dayIn = '%d.%02d.%02d' % (year,month,day)
            listNewFile = modisOgg.getFilesList(dayIn)
            listFilesDown = modisOgg.checkDataExist(listNewFile)
            modisOgg.dayDownload(dayIn, listFilesDown)
#        modisOgg.downloadsAllDay()
    else:
        print("A problem with the connection occured")

                     
def latlon_2modis_tile(lat,lon):
    # reference: https://code.env.duke.edu/projects/mget/wiki/SinusoidalMODIS
    p_modis_grid = Proj('+proj=sinu +R=6371007.181 +nadgrids=@null +wktext')
    x, y = p_modis_grid(lon, lat)
    # or the inverse, from x, y to lon, lat
    lon, lat = p_modis_grid(x, y, inverse=True)
    tileWidth = 1111950.5196666666
    ulx = -20015109.354
    uly = -10007554.677
    H = (x-ulx)/tileWidth
    V = 18-((y-uly)/tileWidth)
    return int(V),int(H)
    
#def geotiff_2envi():   
#    #geotiffConvert = 'GeoTiff2ENVI'
#    bands = ["blue","green","red","nir","swir1","swir2","cloud"]
#    l8bands = ["sr_band2","sr_band3","sr_band4","sr_band5","sr_band6","sr_band7","pixel_qa"] 
#    
#    landsat_files = glob.glob(os.path.join(landsat_temp,"*_MTL.txt"))
#    for i in range(len(landsat_files)):
#        
#        fn = landsat_files[i][:-8]
#        meta = landsat_metadata(landsat_files[i])
#        fstem = os.path.join(os.sep.join((fn.split(os.sep)[:-1])),meta.LANDSAT_SCENE_ID)
#
#        for i in range(len(bands)):
#            tifFile = fn+"_%s.tif" % l8bands[i]
#            datFile = fstem+"_%s.%s.dat" % (l8bands[i],bands[i])
#            #subprocess.call(["%s" % geotiffConvert ,"%s" % tifFile, "%s" % datFile])
#            subprocess.call(["gdal_translate","-q", "-of", "ENVI","%s" % tifFile, "%s" % datFile ])
#            os.rename("%s.hdr" % datFile[:-4],"%s.dat.hdr" % datFile[:-4])
            
def geotiff_2envi(paths,productIDs):   
    #geotiffConvert = 'GeoTiff2ENVI'
    bands = ["blue","green","red","nir","swir1","swir2","cloud"]
    l8bands = ["sr_band2","sr_band3","sr_band4","sr_band5","sr_band6","sr_band7","pixel_qa"] 
    count=0
    for productID in productIDs:
        fstem = os.path.join(paths[count],productID)
        count+=1
        for i in range(len(bands)):
            tifFile = fstem+"_%s.tif" % l8bands[i]
            datFile = fstem+"_%s.%s.dat" % (l8bands[i],bands[i])
            #subprocess.call(["%s" % geotiffConvert ,"%s" % tifFile, "%s" % datFile])
            subprocess.call(["gdal_translate","-q", "-of", "ENVI","%s" % tifFile, "%s" % datFile ])
            os.rename("%s.hdr" % datFile[:-4],"%s.dat.hdr" % datFile[:-4])

#def sample():    
#    sample = 'lndlai_sample'
#    convert = 'lndqa2fmask'
#    bands = ["blue","green","red","nir","swir1","swir2","cloud"]
#    l8bands = ["sr_band2","sr_band3","sr_band4","sr_band5","sr_band6","sr_band7","cfmask"] 
#    
#    landsat_files = glob.glob(os.path.join(landsat_temp,"*_MTL.txt"))
#    
#    for i in range(len(landsat_files)):
#        #sceneID = landsat_files[i].split(os.sep)[-1][:-4]
#        meta = landsat_metadata(landsat_files[i])
#        sceneID = meta.LANDSAT_SCENE_ID
#        
#        # extract the Landsat doy and year
#        d = meta.DATETIME_OBJ
#        year = d.year
#        ldoy = sceneID[13:16]
##        year = int(sceneID[9:13])
#        # convert to date    
##        dd = datetime.datetime(year, 1, 1) + datetime.timedelta(int(ldoy) - 1)
##        date = '%d-%02d-%02d' % (dd.year,dd.month,dd.day)
#        date = meta.DATE_ACQUIRED
#        # find the 4 day MODIS doy prior to the Landsat doy
#        mdoy = int((int((float(ldoy)-1.)/4.)*4.)+1)
#        
#        modFiles = glob.glob(os.path.join(modis_base,"MCD15A3H","MCD15A3H.A%s%s.*.hdf" % (year,mdoy)))
#
#        #fstem = landsat_files[i][:-4]
#        fn = landsat_files[i][:-8]
#        fstem = os.path.join(os.sep.join((fn.split(os.sep)[:-1])),meta.LANDSAT_SCENE_ID)
#        lai_path = landsat_LAI
#        if not os.path.exists(lai_path):
#            os.mkdir(lai_path)
#        sam_file = os.path.join(lai_path,"SR_LAI.%s.%s.MCD15A3H_A%s%s.txt" %(date,sceneID,year,mdoy))
#        #====convert the qa to cfmask=====
#        datFile_qa = fstem+"_%s.%s.dat" % ("pixel_qa",bands[6])
#        datFile_cfmask = fstem+"_%s.%s.dat" % (l8bands[6],bands[6])
#        subprocess.call(["%s" % convert, "-lndsr", "%s" % datFile_qa, "-cmask", "%s" % datFile_cfmask])
#        
#        for i in range(len(modFiles)):  
#            fn = os.path.join(lai_path,"slai%d.inp" % i)
#            file = open(fn, "w")
#            file.write("LANDSAT_BASE_BLUE = %s_%s.%s.dat\n" % (fstem,l8bands[0],bands[0]))
#            file.write("LANDSAT_BASE_GREEN = %s_%s.%s.dat\n" % (fstem,l8bands[1],bands[1]))
#            file.write("LANDSAT_BASE_RED = %s_%s.%s.dat\n" % (fstem,l8bands[2],bands[2]))
#            file.write("LANDSAT_BASE_NIR = %s_%s.%s.dat\n" % (fstem,l8bands[3],bands[3]))
#            file.write("LANDSAT_BASE_SWIR1 = %s_%s.%s.dat\n" % (fstem,l8bands[4],bands[4]))
#            file.write("LANDSAT_BASE_SWIR2 = %s_%s.%s.dat\n" % (fstem,l8bands[5],bands[5]))
#            file.write("LANDSAT_BASE_CLOUD = %s_%s.%s.dat\n" % (fstem,l8bands[6],bands[6]))
#            file.write("MODIS_BASE_FILE = %s\n" % modFiles[i])
#            file.write("SAMPLE_FILE_OUT = %s\n" % sam_file)
#            file.write("PURE_SAMPLE_TH = 0.2\n")
#            file.close()
#        
#            subprocess.call(["%s" % sample , "%s" % fn])
#            #os.remove(os.path.join(lai_path,"slai%d.inp" % i))
            
def sample(paths,productIDs,MODIS_product):    
    sample = 'lndlai_sample'
    convert = 'lndqa2fmask'
    bands = ["blue","green","red","nir","swir1","swir2","cloud"]
    l8bands = ["sr_band2","sr_band3","sr_band4","sr_band5","sr_band6","sr_band7","cfmask"] 

    count = 0
    for productID in productIDs:
        fstem = os.path.join(paths[count],productID)
        count+=1
        meta = landsat_metadata(fstem+"_MTL.txt")
        sceneID = meta.LANDSAT_SCENE_ID
        
        # extract the Landsat doy and year
        d = meta.DATETIME_OBJ
        year = d.year
        ldoy = sceneID[13:16]
        date = meta.DATE_ACQUIRED
        # find the 4 day MODIS doy prior to the Landsat doy
        mdoy = int((int((float(ldoy)-1.)/4.)*4.)+1)
        
        modFiles = glob.glob(os.path.join(modis_base,"%s","%s.A%s%s.*.hdf" % (MODIS_product,MODIS_product,year,mdoy)))
        sam_file = os.path.join(landsat_LAI,"SR_LAI.%s.%s.%s_A%s%s.txt" %(date,sceneID,MODIS_product,year,mdoy))
        #====convert the qa to cfmask=====
        datFile_qa = fstem+"_%s.%s.dat" % ("pixel_qa",bands[6])
        datFile_cfmask = fstem+"_%s.%s.dat" % (l8bands[6],bands[6])
        subprocess.call(["%s" % convert, "-lndsr", "%s" % datFile_qa, "-cmask", "%s" % datFile_cfmask])
        
        for i in range(len(modFiles)):  
            fn = os.path.join(landsat_LAI,"slai%d.inp" % i)
            file = open(fn, "w")
            file.write("LANDSAT_BASE_BLUE = %s_%s.%s.dat\n" % (fstem,l8bands[0],bands[0]))
            file.write("LANDSAT_BASE_GREEN = %s_%s.%s.dat\n" % (fstem,l8bands[1],bands[1]))
            file.write("LANDSAT_BASE_RED = %s_%s.%s.dat\n" % (fstem,l8bands[2],bands[2]))
            file.write("LANDSAT_BASE_NIR = %s_%s.%s.dat\n" % (fstem,l8bands[3],bands[3]))
            file.write("LANDSAT_BASE_SWIR1 = %s_%s.%s.dat\n" % (fstem,l8bands[4],bands[4]))
            file.write("LANDSAT_BASE_SWIR2 = %s_%s.%s.dat\n" % (fstem,l8bands[5],bands[5]))
            file.write("LANDSAT_BASE_CLOUD = %s_%s.%s.dat\n" % (fstem,l8bands[6],bands[6]))
            file.write("MODIS_BASE_FILE = %s\n" % modFiles[i])
            file.write("SAMPLE_FILE_OUT = %s\n" % sam_file)
            file.write("PURE_SAMPLE_TH = 0.2\n")
            file.close()
        
            subprocess.call(["%s" % sample , "%s" % fn])

#def train():    
#    cubist = 'cubist'
#    landsat_files = glob.glob(os.path.join(landsat_LAI,"*.txt"))
#    #======combine input data======================================
#    df = pd.DataFrame(columns=['ulx','uly','blue',
#        'green','red','nir','swir1','swir2','ndvi','ndwi','lai','weight','satFlag'])
#    for i in range(len(landsat_files)):
#        sam_file = landsat_files[i]
#    
#        df = df.append(pd.read_csv(sam_file,delim_whitespace=True,names=['ulx','uly','blue',
#        'green','red','nir','swir1','swir2','ndvi','ndwi','lai','weight','satFlag']),ignore_index=True)
#    
#    #=====create filestem.data====================================
#    df = df[(df.satFlag=='N')]
#    df = df.sort_values(by='weight')
#    start_date='200'
#    end_date = '300'
#    filestem = os.path.join(landsat_LAI,"lndsr_modlai_samples.combined_%s-%s" %(start_date,end_date))
#    df.to_csv(os.path.join(landsat_LAI,filestem+".data"), columns = ['blue','green','red',
#    'nir','swir1','swir2','ndvi','ndwi','lai','weight'],header=None, 
#    index=None, mode='w',  sep="\t", encoding='utf-8')
#    
#    #====create filestem.names====================================
#    fn = os.path.join(landsat_LAI,"%s.names" % filestem)
#    file = open(fn, "w")
#    file.write("lai.\n")
#    file.write("B1: continuous\n")
#    file.write("B2: continuous\n")
#    file.write("B3: continuous\n")
#    file.write("B4: continuous\n")
#    file.write("B5: continuous\n")
#    file.write("B7: continuous\n")
#    file.write("ndvi: continuous\n")
#    file.write("ndwi: continuous\n")
#    file.write("lai: continuous\n")
#    file.write("case weight: continuous\n")
#    file.write("attributes excluded: B1, B2, B7, ndvi, ndwi\n")
#    file.close()
#    
#    nrules = 5
#    subprocess.call(["%s" % cubist , "-f" ,"%s" % filestem, "-r", "%d" % nrules, "-u"])
            
def train(paths,productIDs,MODIS_product):    
    cubist = 'cubist'
    #======combine input data======================================
    df = pd.DataFrame(columns=['ulx','uly','blue',
        'green','red','nir','swir1','swir2','ndvi','ndwi','lai','weight','satFlag'])
    count = 0
    for productID in productIDs:
        fstem = os.path.join(paths[count],productID)
        count+=1
        meta = landsat_metadata(fstem+"_MTL.txt")
        sceneID = meta.LANDSAT_SCENE_ID
        # extract the Landsat doy and year
        d = meta.DATETIME_OBJ
        year = d.year
        ldoy = sceneID[13:16]
        date = meta.DATE_ACQUIRED
        # find the 4 day MODIS doy prior to the Landsat doy
        mdoy = int((int((float(ldoy)-1.)/4.)*4.)+1)
        sam_file = os.path.join(landsat_LAI,"SR_LAI.%s.%s.%s_A%s%s.txt" %(date,sceneID,MODIS_product,year,mdoy))
    
        df = df.append(pd.read_csv(sam_file,delim_whitespace=True,names=['ulx','uly','blue',
        'green','red','nir','swir1','swir2','ndvi','ndwi','lai','weight','satFlag']),ignore_index=True)
    
    #=====create filestem.data====================================
    df = df[(df.satFlag=='N')]
    df = df.sort_values(by='weight')
    start_date='200'
    end_date = '300'
    filestem = os.path.join(landsat_LAI,"lndsr_modlai_samples.combined_%s-%s" %(start_date,end_date))
    df.to_csv(os.path.join(landsat_LAI,filestem+".data"), columns = ['blue','green','red',
    'nir','swir1','swir2','ndvi','ndwi','lai','weight'],header=None, 
    index=None, mode='w',  sep="\t", encoding='utf-8')
    
    #====create filestem.names====================================
    fn = os.path.join(landsat_LAI,"%s.names" % filestem)
    file = open(fn, "w")
    file.write("lai.\n")
    file.write("B1: continuous\n")
    file.write("B2: continuous\n")
    file.write("B3: continuous\n")
    file.write("B4: continuous\n")
    file.write("B5: continuous\n")
    file.write("B7: continuous\n")
    file.write("ndvi: continuous\n")
    file.write("ndwi: continuous\n")
    file.write("lai: continuous\n")
    file.write("case weight: continuous\n")
    file.write("attributes excluded: B1, B2, B7, ndvi, ndwi\n")
    file.close()
    
    nrules = 5
    subprocess.call(["%s" % cubist , "-f" ,"%s" % filestem, "-r", "%d" % nrules, "-u"])
            
#def compute():    
#    lndbio ='lndlai_compute'
#    bands = ["blue","green","red","nir","swir1","swir2","cloud"]
#    l8bands = ["sr_band2","sr_band3","sr_band4","sr_band5","sr_band6","sr_band7","cfmask"] 
#    
#    landsat_files = glob.glob(os.path.join(landsat_temp,"*_MTL.txt"))
#    for i in range(len(landsat_files)):
##        sceneID = landsat_files[i].split(os.sep)[-1][:-4]  
#        meta = landsat_metadata(landsat_files[i])
#        sceneID = meta.LANDSAT_SCENE_ID
#        #fstem = landsat_files[i][:-4]       
#        fn = landsat_files[i][:-8]
#        fstem = os.path.join(os.sep.join((fn.split(os.sep)[:-1])),meta.LANDSAT_SCENE_ID)
#        # create a folder for lai if it does not exist
#        #lai_path = os.path.join(landsat_LAI,'%s' % sceneID[9:16])
#        lai_path = os.path.join(landsat_LAI,'%s' % sceneID[3:9])
#        if not os.path.exists(lai_path):
#            os.mkdir(lai_path)
#        ndvi_path = os.path.join(landsat_NDVI,'%s' % sceneID[3:9])
#        if not os.path.exists(ndvi_path):
#            os.mkdir(ndvi_path)
#        cfmask_path = os.path.join(landsat_Mask,'%s' % sceneID[3:9])
#        if not os.path.exists(cfmask_path):
#            os.mkdir(cfmask_path)
#        start_date='200'
#        end_date = '300'
#        filestem = os.path.join(landsat_LAI,"lndsr_modlai_samples.combined_%s-%s" %(start_date,end_date))
#        laiFN = os.path.join(landsat_LAI,"lndlai.%s.hdf" % sceneID)
#        fn = os.path.join(landsat_LAI,"compute_lai.inp")
#        file = open(fn, "w")
#        file.write("LANDSAT_BASE_BLUE = %s_%s.%s.dat\n" % (fstem,l8bands[0],bands[0]))
#        file.write("LANDSAT_BASE_GREEN = %s_%s.%s.dat\n" % (fstem,l8bands[1],bands[1]))
#        file.write("LANDSAT_BASE_RED = %s_%s.%s.dat\n" % (fstem,l8bands[2],bands[2]))
#        file.write("LANDSAT_BASE_NIR = %s_%s.%s.dat\n" % (fstem,l8bands[3],bands[3]))
#        file.write("LANDSAT_BASE_SWIR1 = %s_%s.%s.dat\n" % (fstem,l8bands[4],bands[4]))
#        file.write("LANDSAT_BASE_SWIR2 = %s_%s.%s.dat\n" % (fstem,l8bands[5],bands[5]))
#        file.write("LANDSAT_BASE_CLOUD = %s_%s.%s.dat\n" % (fstem,l8bands[6],bands[6]))
#        file.write("LANDSAT_ANC_FILE = %s\n" % filestem)
#        file.write("BIOPHYSICS_PARA_FILE_OUT = %s\n" % laiFN)
#        file.close()
#        
#        subprocess.call(["%s" % lndbio , "%s" % fn])
#        #====convert to geotiff=========
#        outlaifn = os.path.join(lai_path,'%s_lai.tiff' % sceneID)
#        outndvifn = os.path.join(ndvi_path,'%s_ndvi.tiff' % sceneID)
#        outcfmaskfn = os.path.join(cfmask_path,'%s_Mask.tiff' % sceneID)
##        tempfn = os.path.join(lai_path,'temp.tiff')
#        subprocess.call(["gdal_translate", 'HDF4_EOS:EOS_GRID:"%s":LANDSAT:LAI' % laiFN, "%s" % outlaifn])
##        subprocess.call(["gdal_calc.py", "-A %s" % tempfn,  "--outfile=%s" % outlaifn,
##                                 "--type=UInt16", "--overwrite", '--calc="A"'])
#        subprocess.call(["gdal_translate", 'HDF4_EOS:EOS_GRID:"%s":LANDSAT:NDVI' % laiFN, "%s" % outndvifn])
#        
#        subprocess.call(["gdal_translate", 'HDF4_EOS:EOS_GRID:"%s":LANDSAT:cfmask' % laiFN, "%s" % outcfmaskfn])
##        subprocess.call(["gdal_calc.py", "-A %s" % tempfn,  "--outfile=%s" % outndvifn, 
##                         "--type=UInt16", "--overwrite", '--calc="A"'])
##        shutil.move(laiFN,os.path.join(lai_path,"lndlai.%s.hdf" % sceneID))
#        os.remove(fn)
##        os.remove(tempfn)
#    #=====CLEANING UP========
#    filelist = [ f for f in os.listdir(landsat_LAI) if f.startswith("lndsr_modlai_samples") ]
#    for f in filelist:
#        os.remove(os.path.join(landsat_LAI,f))

def compute(paths,productIDs,MODIS_product,sat):    
    lndbio ='lndlai_compute'
    bands = ["blue","green","red","nir","swir1","swir2","cloud"]
    l8bands = ["sr_band2","sr_band3","sr_band4","sr_band5","sr_band6","sr_band7","cfmask"] 
    
    count = 0
    for productID in productIDs:
        fstem = os.path.join(paths[count],productID)
        count+=1
        meta = landsat_metadata(fstem+"_MTL.txt")
        sceneID = meta.LANDSAT_SCENE_ID
        scene = sceneID[3:9]
        folder = os.path.join(cacheDir,"L%d" % sat,scene)
        lai_path = os.path.join(folder,"LAI")
        if not os.path.exists(lai_path):
            os.mkdir(lai_path)
        ndvi_path = os.path.join(folder,"NDVI")
        if not os.path.exists(ndvi_path):
            os.mkdir(ndvi_path)
        cfmask_path = os.path.join(folder,"CF_MASK")
        if not os.path.exists(cfmask_path):
            os.mkdir(cfmask_path)
        start_date='200'
        end_date = '300'
        filestem = os.path.join(landsat_LAI,"lndsr_modlai_samples.combined_%s-%s" %(start_date,end_date))
        laiFN = os.path.join(landsat_LAI,"lndlai.%s.hdf" % sceneID)
        fn = os.path.join(landsat_LAI,"compute_lai.inp")
        file = open(fn, "w")
        file.write("LANDSAT_BASE_BLUE = %s_%s.%s.dat\n" % (fstem,l8bands[0],bands[0]))
        file.write("LANDSAT_BASE_GREEN = %s_%s.%s.dat\n" % (fstem,l8bands[1],bands[1]))
        file.write("LANDSAT_BASE_RED = %s_%s.%s.dat\n" % (fstem,l8bands[2],bands[2]))
        file.write("LANDSAT_BASE_NIR = %s_%s.%s.dat\n" % (fstem,l8bands[3],bands[3]))
        file.write("LANDSAT_BASE_SWIR1 = %s_%s.%s.dat\n" % (fstem,l8bands[4],bands[4]))
        file.write("LANDSAT_BASE_SWIR2 = %s_%s.%s.dat\n" % (fstem,l8bands[5],bands[5]))
        file.write("LANDSAT_BASE_CLOUD = %s_%s.%s.dat\n" % (fstem,l8bands[6],bands[6]))
        file.write("LANDSAT_ANC_FILE = %s\n" % filestem)
        file.write("BIOPHYSICS_PARA_FILE_OUT = %s\n" % laiFN)
        file.close()
        
        subprocess.call(["%s" % lndbio , "%s" % fn])
        #====convert to geotiff=========
        outlaifn = os.path.join(lai_path,'%s_lai.tiff' % sceneID)
        outndvifn = os.path.join(ndvi_path,'%s_ndvi.tiff' % sceneID)
        outcfmaskfn = os.path.join(cfmask_path,'%s_Mask.tiff' % sceneID)
        subprocess.call(["gdal_translate", 'HDF4_EOS:EOS_GRID:"%s":LANDSAT:LAI' % laiFN, "%s" % outlaifn])
        subprocess.call(["gdal_translate", 'HDF4_EOS:EOS_GRID:"%s":LANDSAT:NDVI' % laiFN, "%s" % outndvifn])        
        subprocess.call(["gdal_translate", 'HDF4_EOS:EOS_GRID:"%s":LANDSAT:cfmask' % laiFN, "%s" % outcfmaskfn])
        os.remove(fn)
    #=====CLEANING UP========
    filelist = [ f for f in os.listdir(landsat_LAI) if f.startswith("lndsr_modlai_samples") ]
    for f in filelist:
        os.remove(os.path.join(landsat_LAI,f))    
def get_LAI(loc,start_date,end_date,usgs_user,usgs_pass,earth_user,
            earth_pass,cloud,sat,cacheDir):    
    # find MODIS tiles that cover landsat scene
    # MODIS products   
    MODIS_product = 'MCD15A3H'
    version = '006'
    [v,h] = latlon_2modis_tile(loc[0],loc[1])
    tiles = "h%02dv%02d" %(h,v)
    #====search for available data=============================================
    available = 'Y'
    search_df = getlandsatdata.search(loc[0],loc[1],start_date,end_date,cloud,available,cacheDir,sat)
    productIDs = search_df.LANDSAT_PRODUCT_ID
    paths = search_df.local_file    
    # download MODIS LAI over the same area and time
    print("Downloading MODIS data...")
    get_modis_lai(tiles,MODIS_product,version,start_date,end_date,("%s"% earth_user,"%s"% earth_pass))
    print(paths)
    # Convert Landsat SR downloads to ENVI format
    # Note:  May be some warnings about unknown field - ignore.
    print("Converting Landsat SR to ENVI format...")
    geotiff_2envi(paths,productIDs)
    
    # Generate MODIS-Landsat samples for LAI computation
    print("Generating MODIS-Landsat samples...")
    sample(paths,productIDs,MODIS_product)    
    
    # Compute Landsat LAI
    print("Computing Landsat LAI...")
    train(paths,productIDs,MODIS_product)
    compute(paths,productIDs,MODIS_product)    

def main():
    # Get time and location from user
    parser = argparse.ArgumentParser()
    parser.add_argument("lat", type=float, help="latitude")
    parser.add_argument("lon", type=float, help="longitude")
    parser.add_argument("start_date", type=str, help="Start date yyyy-mm-dd")
    parser.add_argument("end_date", type=str, help="Start date yyyy-mm-dd")
    parser.add_argument("cloud", type=int, help="cloud coverage")
    parser.add_argument('-s','--sat', nargs='?',type=int, default=8, help='which landsat to search or download, i.e. Landsat 8 = 8')
    args = parser.parse_args()
      
    loc = [args.lat,args.lon] 
    start_date = args.start_date
    end_date = args.end_date
    cloud = args.cloud
    sat = args.sat
    
    # =====USGS credentials===============
     # need to get this from pop up
    usgs_user = str(getpass.getpass(prompt="usgs username:"))
    if keyring.get_password("usgs",usgs_user)==None:
        usgs_pass = str(getpass.getpass(prompt="usgs password:"))
        keyring.set_password("usgs",usgs_user,usgs_pass)
    else:
        usgs_pass = str(keyring.get_password("usgs",usgs_user)) 
    
    
     # =====earthData credentials===============
    earth_user = str(getpass.getpass(prompt="earth login username:"))
    if keyring.get_password("nasa",earth_user)==None:
        earth_pass = str(getpass.getpass(prompt="earth login password:"))
        keyring.set_password("nasa",earth_user,earth_pass)
    else:
        earth_pass = str(keyring.get_password("nasa",earth_user)) 
    
    get_LAI(loc,start_date,end_date,usgs_user,usgs_pass,earth_user,
            earth_pass,cloud,sat,cacheDir)

    print("All done with LAI")


if __name__ == "__main__":
    try:
        main()
    except (KeyboardInterrupt, pycurl.error):
        exit('Received Ctrl + C... Exiting! Bye.', 1)   
    