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
import sqlite3
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
cacheDir = os.path.abspath(os.path.join(base,os.pardir,"SATELLITE_DATA"))
Folders = folders(base)   
modis_base = Folders['modis_base']
landsat_SR = Folders['landsat_SR']
landsat_LAI = Folders['landsat_LAI']
landsat_NDVI = Folders['landsat_NDVI']
landsat_Mask = Folders['landsat_Mask']
landsat_temp = os.path.join(landsat_SR,'temp')
if not os.path.exists(landsat_temp):
    os.mkdir(landsat_temp)

def updateModisDB(filenames,cacheDir):
    
    db_fn = os.path.join(cacheDir,"modis_db.db")
    fn = filenames[0].split(os.sep)[-1]
    product = fn.split('.')[0]
    years = []
    doys = []
    tiles = []
    fns = []
    for filename in filenames:
        fn = filename.split(os.sep)[-1]
        fns.append(filename)
        years.append(fn.split('.')[1][1:5])
        doys.append(fn.split('.')[1][5:9])
        tiles.append(fn.split('.')[2])
    if not os.path.exists(db_fn):
        conn = sqlite3.connect( db_fn )
        modis_dict = {"TILE":tiles,"YEAR":years,"DOY": doys,"filename":fns}
        modis_df = pd.DataFrame.from_dict(modis_dict)
        modis_df.to_sql("%s" % product, conn, if_exists="replace", index=False)
        conn.close()
    else:
        conn = sqlite3.connect( db_fn )
        orig_df = pd.read_sql_query("SELECT * from %s" % product,conn)
        modis_dict = {"TILE":tiles,"YEAR":years,"DOY": doys,"filename":fns}
        modis_df = pd.DataFrame.from_dict(modis_dict)
        orig_df = orig_df.append(modis_df,ignore_index=True)
        orig_df = orig_df.drop_duplicates(keep='last')
        orig_df.to_sql("%s" % product, conn, if_exists="replace", index=False)
        conn.close()

def searchModisDB(tiles,start_date,end_date,product,cacheDir):
    db_fn = os.path.join(cacheDir,"modis_db.db")
    conn = sqlite3.connect( db_fn )
    startdd = datetime.datetime.strptime(start_date, '%Y-%m-%d')
    enddd = datetime.datetime.strptime(end_date, '%Y-%m-%d')
    numDays= (enddd-startdd).days
    laidates = np.array(range(1,366,4))
    df1 = pd.DataFrame.from_dict({"TILE":[],"YEAR":[],"DOY": [],"filename":[]})
    df2 = pd.DataFrame.from_dict({"TILE":[],"YEAR":[],"DOY": []})
    for tile in tiles:
        for i in range(numDays+1):
            dd = startdd+datetime.timedelta(days=i)            
            year = dd.year
            doy = (dd-datetime.datetime(year,1,1,0,0)).days+1
            rday = laidates[laidates>=doy][0]
            if (doy==rday):
                dd = datetime.datetime(year,1,1,0,0)+datetime.timedelta(days=rday-1)
                year = dd.year
                df = pd.read_sql_query("SELECT * from %s WHERE (TILE = '%s')"
                               "AND (YEAR =  '%d') AND (DOY = '%03d' )" % 
                               (product,tile,year,rday),conn)
                df1 = df1.append(df,ignore_index=True)
                df1 = df1[["DOY","TILE","YEAR"]]
                row = pd.Series({"TILE":str(tile),"YEAR":str(year),"DOY": str(rday)})
                df2 = df2.append(row,ignore_index=True)
    merged = df2.merge(df1, indicator=True, how='outer')
    df3 = merged[merged['_merge'] != 'both' ]
    out_df = df3[["DOY","TILE","YEAR"]]
    conn.close()
    return out_df

    
def updateLandsatProductsDB(landsatDB,filenames,cacheDir,product):
    
    db_fn = os.path.join(cacheDir,"landsat_products.db")
    
    date = landsatDB.acquisitionDate
    ullat = landsatDB.upperLeftCornerLatitude
    ullon = landsatDB.upperLeftCornerLongitude
    lllat = landsatDB.lowerRightCornerLatitude
    lllon = landsatDB.lowerRightCornerLongitude
    productIDs = landsatDB.LANDSAT_PRODUCT_ID
    
    if not os.path.exists(db_fn):
        conn = sqlite3.connect( db_fn )
        landsat_dict = {"acquisitionDate":date,"upperLeftCornerLatitude":ullat,
                      "upperLeftCornerLongitude":ullon,
                      "lowerRightCornerLatitude":lllat,
                      "lowerRightCornerLongitude":lllon,
                      "LANDSAT_PRODUCT_ID":productIDs,"filename":filenames}
        landsat_df = pd.DataFrame.from_dict(landsat_dict)
        landsat_df.to_sql("%s" % product, conn, if_exists="replace", index=False)
        conn.close()
    else:
        conn = sqlite3.connect( db_fn )
        res = conn.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = res.fetchall()[0]
        if (product in tables):
            orig_df = pd.read_sql_query("SELECT * from %s" % product,conn)
        else:
            orig_df = pd.DataFrame()
            
        landsat_dict = {"acquisitionDate":date,"upperLeftCornerLatitude":ullat,
                      "upperLeftCornerLongitude":ullon,
                      "lowerRightCornerLatitude":lllat,
                      "lowerRightCornerLongitude":lllon,
                      "LANDSAT_PRODUCT_ID":productIDs,"filename":filenames}
        landsat_df = pd.DataFrame.from_dict(landsat_dict)
        orig_df = orig_df.append(landsat_df,ignore_index=True)
        orig_df = orig_df.drop_duplicates(keep='last')
        orig_df.to_sql("%s" % product, conn, if_exists="replace", index=False)
        conn.close()

def searchLandsatProductsDB(lat,lon,start_date,end_date,product,cacheDir):
    db_fn = os.path.join(cacheDir,"landsat_products.db")
    conn = sqlite3.connect( db_fn )

    out_df = pd.read_sql_query("SELECT * from %s WHERE (acquisitionDate >= '%s')"
                   "AND (acquisitionDate < '%s') AND (upperLeftCornerLatitude > %f )"
                   "AND (upperLeftCornerLongitude < %f ) AND "
                   "(lowerRightCornerLatitude < %f) AND "
                   "(lowerRightCornerLongitude > %f)" % 
                   (product,start_date,end_date,lat,lon,lat,lon),conn)   
    conn.close()
    return out_df
        
def get_modis_lai(tiles,product,version,start_date,end_date,auth,cacheDir): 
    db_fn = os.path.join(cacheDir,"modis_db.db")
    
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
    product_path = os.path.join(cacheDir,product)   
    if not os.path.exists(product_path):
        os.makedirs(product_path)
    if os.path.exists(db_fn):
        out_df = searchModisDB(tiles,start_date,end_date,product,cacheDir) 
        print(len(out_df))
        filenames = []
        for i in range(len(out_df)):
            year = int(out_df.YEAR[i])
            doy = int(out_df.DOY[i])
            fdate = datetime.date.fromordinal(datetime.date(year, 1, 1).toordinal() + doy - 1).isoformat()
            modisOgg = downModis(url="https://e4ftl01.cr.usgs.gov", user=auth[0],
                                           password=auth[1],
                                           destinationFolder=product_path,
                                           tiles=tiles, path=folder,
                                           product="%s.%s" % (product,version), delta=1,
                                           today=fdate)
            modisOgg.connect()
            day = modisOgg.getListDays()[0]
            listAllFiles = modisOgg.getFilesList(day)
            listFilesDown = modisOgg.checkDataExist(listAllFiles)
            filenames.append(os.path.join(product_path,listFilesDown))
            modisOgg.dayDownload(day, listFilesDown)
            modisOgg.closeFilelist()
            
    else:
        
        modisOgg = downModis(url="https://e4ftl01.cr.usgs.gov", destinationFolder=product_path, 
                             user=auth[0], password=auth[1], tiles=tiles, path=folder, 
                             product="%s.%s" % (product,version),today=start_date,enddate=end_date)
                # connect to http or ftp
        modisOgg.connect()
        if modisOgg.nconnection <= 20:
            # download data
            filenames = []
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
                print(listFilesDown)
                if (len(listFilesDown) < 1):
                    continue
                filenames.append(os.path.join(product_path,listFilesDown[0]))
                modisOgg.dayDownload(dayIn, listFilesDown)
    #        modisOgg.downloadsAllDay()
        else:
            print("A problem with the connection occured")
    return filenames

                     
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
        
        modFiles = glob.glob(os.path.join(modis_base,"%s" % MODIS_product,"%s.A%d%03d.*.hdf" % (MODIS_product,year,mdoy)))
        sam_file = os.path.join(landsat_LAI,"SR_LAI.%s.%s.%s_A%d%03d.txt" %(date,sceneID,MODIS_product,year,mdoy))
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
        sam_file = os.path.join(landsat_LAI,"SR_LAI.%s.%s.%s_A%d%03d.txt" %(date,sceneID,MODIS_product,year,mdoy))
    
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

def compute(paths,productIDs,MODIS_product,sat,cacheDir):    
    lndbio ='lndlai_compute'
    bands = ["blue","green","red","nir","swir1","swir2","cloud"]
    l8bands = ["sr_band2","sr_band3","sr_band4","sr_band5","sr_band6","sr_band7","cfmask"] 
    
    count = 0
    lai_fns = []
    ndvi_fns = []
    mask_fns = []
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
        lai_fns.append(outlaifn)
        outndvifn = os.path.join(ndvi_path,'%s_ndvi.tiff' % sceneID)
        ndvi_fns.append(outndvifn)
        outcfmaskfn = os.path.join(cfmask_path,'%s_Mask.tiff' % sceneID)
        mask_fns.append(outcfmaskfn)
        subprocess.call(["gdal_translate", 'HDF4_EOS:EOS_GRID:"%s":LANDSAT:LAI' % laiFN, "%s" % outlaifn])
        subprocess.call(["gdal_translate", 'HDF4_EOS:EOS_GRID:"%s":LANDSAT:NDVI' % laiFN, "%s" % outndvifn])        
        subprocess.call(["gdal_translate", 'HDF4_EOS:EOS_GRID:"%s":LANDSAT:cfmask' % laiFN, "%s" % outcfmaskfn])
        os.remove(fn)

    #=====CLEANING UP========
    filelist = [ f for f in os.listdir(landsat_LAI) if f.startswith("lndsr_modlai_samples") ]
    for f in filelist:
        os.remove(os.path.join(landsat_LAI,f))    
    return lai_fns,ndvi_fns,mask_fns
def get_LAI(loc,start_date,end_date,earth_user,earth_pass,cloud,sat,cacheDir): 
    landsatCacheDir = os.path.join(cacheDir,"LANDSAT")
    modisCacheDir = os.path.join(cacheDir,"MODIS")
    # find MODIS tiles that cover landsat scene
    # MODIS products   
    MODIS_product = 'MCD15A3H'
    version = '006'
    [v,h] = latlon_2modis_tile(loc[0],loc[1])
    tiles = ["h%02dv%02d" %(h,v)]
    #====search for available data=============================================
    available = 'Y'
    search_df = getlandsatdata.search(loc[0],loc[1],start_date,end_date,cloud,available,landsatCacheDir,sat)
    productIDs = search_df.LANDSAT_PRODUCT_ID
    paths = search_df.local_file_path   
    # download MODIS LAI over the same area and time
    print("Downloading MODIS data...")
    modis_files = get_modis_lai(tiles,MODIS_product,version,start_date,end_date,("%s"% earth_user,"%s"% earth_pass),modisCacheDir)
    updateModisDB(modis_files,modisCacheDir)
    
    #====check what products are done against what Landsat data is available===
    processedProductIDs = searchLandsatProductsDB(loc[0],loc[1],start_date,end_date,"LAI",landsatCacheDir)
    df1 = processedProductIDs[["LANDSAT_PRODUCT_ID"]]
    merged = df1.merge(productIDs, indicator=True, how='outer')
    df3 = merged[merged['_merge'] != 'both' ]
    productIDs = df3[["LANDSAT_PRODUCT_ID"]]
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
    lai_fns,ndvi_fns,mask_fns = compute(paths,productIDs,MODIS_product,sat,landsatCacheDir)  
    
    updateLandsatProductsDB(search_df,lai_fns,landsatCacheDir,'LAI')
    updateLandsatProductsDB(search_df,ndvi_fns,landsatCacheDir,'NDVI')
    updateLandsatProductsDB(search_df,mask_fns,landsatCacheDir,'CF_MASK')

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
#    usgs_user = str(getpass.getpass(prompt="usgs username:"))
#    if keyring.get_password("usgs",usgs_user)==None:
#        usgs_pass = str(getpass.getpass(prompt="usgs password:"))
#        keyring.set_password("usgs",usgs_user,usgs_pass)
#    else:
#        usgs_pass = str(keyring.get_password("usgs",usgs_user)) 
    
    
     # =====earthData credentials===============
    earth_user = str(getpass.getpass(prompt="earth login username:"))
    if keyring.get_password("nasa",earth_user)==None:
        earth_pass = str(getpass.getpass(prompt="earth login password:"))
        keyring.set_password("nasa",earth_user,earth_pass)
    else:
        earth_pass = str(keyring.get_password("nasa",earth_user)) 
    
    get_LAI(loc,start_date,end_date,earth_user,earth_pass,cloud,sat,cacheDir)

    print("All done with LAI")


if __name__ == "__main__":
    try:
        main()
    except (KeyboardInterrupt, pycurl.error):
        exit('Received Ctrl + C... Exiting! Bye.', 1)   
    