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
import shutil
import pandas as pd
from datetime import datetime
import argparse
import getpass
import keyring
import json
from pyproj import Proj
from .utils import folders,search,check_order_cache
from .Clients import Client
from .Downloaders import BaseDownloader
import pycurl
from .landsatTools import landsat_metadata
import requests
from time import sleep
import logging
from pymodis.downmodis import downModis
logging.getLogger("urllib3").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.DEBUG)

# The current URL hosting the ESPA interfaces has reached a stable version 1.0
host = 'https://espa.cr.usgs.gov/api/v1/'
TIMEOUT=86400

base = os.getcwd()
Folders = folders(base)   
modis_base = Folders['modis_base']
landsat_SR = Folders['landsat_SR']
landsat_LAI = Folders['landsat_LAI']
landsat_temp = os.path.join(landsat_SR,'temp')
if not os.path.exists(landsat_temp):
    os.mkdir(landsat_temp)

def espa_api(endpoint, verb='get', body=None, uauth=None):
    """ Suggested simple way to interact with the ESPA JSON REST API """
#    auth_tup = uauth if uauth else print "need USGS creds!" exit()
    if uauth:
        auth_tup = uauth
    else:
        print "need USGS creds!" 
        exit()
        
    response = getattr(requests, verb)(host + endpoint, auth=auth_tup, json=body)
    print('{} {}'.format(response.status_code, response.reason))
    data = response.json()
    if isinstance(data, dict):
        messages = data.pop("messages", None)  
        if messages:
            print(json.dumps(messages, indent=4))
    try:
        response.raise_for_status()
    except Exception as e:
        print(e)
        return None
    else:
        return data
    
def download_order_gen(order_id, auth, downloader=None, sleep_time=300, timeout=86400, **dlkwargs):
    """
    This function is a generator that yields the results from the input downloader classes
    download() method. This is a generator mostly so that data pipeline functions that operate
    upon freshly downloaded files may immediately get started on them.

    :param order_id:            order name
    :param downloader:          optional downloader for tiles. child of BaseDownloader class
                                of a Downloaders.BaseDownloader or child class
    :param sleep_time:          number of seconds to wait between checking order status
    :param timeout:             maximum number of seconds to run program
    :param dlkwargs:            keyword arguments for downloader.download() method.
    :returns:                   yields values from the input downloader.download() method.
    """

    complete = False
    reached_timeout = False
    starttime = datetime.now()

    if downloader is None:
        downloader = BaseDownloader('espa_downloads')

    while not complete and not reached_timeout:
        # wait a while before the next ping and check timeout condition
        elapsed_time = (datetime.now() - starttime).seconds
        reached_timeout = elapsed_time > timeout
        print("Elapsed time is {0}m".format(elapsed_time / 60.0))

        # check order completion status, and list all items which ARE complete
#            complete_items = self._complete_items(order_id, verbose=False)
        
        filters = {"status": "complete"}  # Here, we ignore any purged orders
        complete_items = espa_api('item-status/{0}'.format(order_id), uauth=auth, body=filters)[order_id]
        for c in complete_items:
            if isinstance(c, dict):
                url = c["product_dload_url"]
            elif isinstance(c, requests.Request):
                url = c.json()["product_dload_url"]
            else:
                raise Exception("Could not interpret {0}".format(c))
            yield downloader.download(url, **dlkwargs)
        resp = espa_api('item-status/{0}'.format(order_id), uauth=auth)
        all_items = resp[order_id]


        active_items = [item for item in all_items
        if item['status'] != 'complete' and
        item['status'] != 'error' and
        item['status'] != 'unavailable' and 
        item['status'] != 'purged']

        complete = (len(active_items) < 1)
        if not complete:
            sleep(sleep_time)

def get_landsat_data(collection,loc,start_date,end_date,auth,cloud):
    username = auth[0]
    password = auth[1]
#    client = Client(auth)
    def api_request(endpoint, verb='get', json=None, uauth=None):
        """
        Here we can see how easy it is to handle calls to a REST API that uses JSON
        """
        auth_tup = uauth if uauth else (username, password)
        response = getattr(requests, verb)(host + endpoint, auth=auth_tup, json=json)
        return response.json()
    


    #=====set products=======
    l8_prods = ['sr','bt']
    #=====search for data=======
    print("Searching...")
    sceneIDs = search(collection,loc[0],loc[1],start_date,end_date,cloud)
    ordered_data = check_order_cache(auth)
    l8_tiles =[]
    orderedIDs_completed = []
    orderedIDs_not_completed = []
    sceneIDs_completed = []
    sceneIDs_not_completed =[]
    for sceneID in sceneIDs:
        if sceneID.startswith('LC'):
            if collection == 0:
                scene = sceneID[3:9]
            else:
                scene = sceneID.split('_')[2]    
            dataFN = os.path.join(landsat_SR,"%s" % scene,"%s_MTL.txt" % sceneID)
            if not os.path.exists(dataFN): # check if we have the data on our system
                if collection > 0: # can't order pre-collection data anymore
                    if not ordered_data.empty:
                        if np.sum(ordered_data.productID==sceneID)>0:
                            completed_test = (ordered_data.productID==sceneID) & (ordered_data.status=='complete')
                            not_complete_test = (ordered_data.productID==sceneID) & (ordered_data.status!='complete')
                            if len(ordered_data[completed_test])>0:
                                orderedIDs_completed.append(list(ordered_data.orderid[completed_test])[0])
                                sceneIDs_completed.append(sceneID)
                            else:
                                orderedIDs_not_completed.append(list(ordered_data.orderid[not_complete_test][0]))
                                sceneIDs_not_completed.append(sceneID)
                        else:
                            l8_tiles.append(sceneID)
                    else:
                        l8_tiles.append(sceneID)
            else:
                files = glob.glob("%s*" % dataFN[:-8])
                for file in files:
                    linked_file = os.path.join(landsat_temp,file.split(os.sep)[-1])
                    if not os.path.exists(linked_file):
                        os.symlink(file,linked_file)

    
    
    if l8_tiles:
        print("Ordering new data...")
        #========setup order=========
        order = espa_api('available-products',uauth=auth, body=dict(inputs=l8_tiles))
        for sensor in order.keys():
            if isinstance(order[sensor], dict) and order[sensor].get('inputs'):
                order[sensor]['products'] = l8_prods
        
        order['format'] = 'gtiff'
        # =======order the data============
        resp = espa_api('order', verb='post',uauth=auth, body=order)
        print(json.dumps(resp, indent=4))
        orderidNew = resp['orderid']
                    
    if orderedIDs_completed:
        
        print("downloading completed existing orders...")
        print orderedIDs_completed
        i = -1
        for orderid in orderedIDs_completed:
            i+=1
            sceneID = sceneIDs_completed[i]
            complete = False
            reached_TIMEOUT = False
            starttime = datetime.now()
            while not complete and not reached_TIMEOUT:
                resp = espa_api('item-status/{0}'.format(orderid),uauth=auth)
                for item in resp[orderid]:
                    if item.get('name')==sceneID:
                        url = item.get('product_dload_url')                      
                        elapsed_time = (datetime.now() - starttime).seconds
                        reached_TIMEOUT = elapsed_time > TIMEOUT
                        print("Elapsed time is {0}m".format(elapsed_time / 60.0))
                        if len(url)>0:
                            downloader = BaseDownloader('espa_downloads')
                            downloader.download(url)
                        #if os.path.exists(os.path.join(os.getcwd,'espa_downloads',url.split(os.sep)[-1][:-7])):
                            complete = True
                        
                        if not complete:
                            sleep(300)
                        
                        
    if orderedIDs_not_completed:
        print("waiting for cached existing orders...")
        i = -1
        for orderid in orderedIDs_not_completed:
            i+=1
            complete = False
            reached_TIMEOUT = False
            starttime = datetime.now()
            sceneID = sceneIDs_not_completed[i]
            while not complete and not reached_TIMEOUT:
                resp = api_request('item-status/{0}'.format(orderid))
                for item in resp[orderid]:
                    if item.get('name')==sceneID:
                        url = item.get('product_dload_url')
                        elapsed_time = (datetime.now() - starttime).seconds
                        reached_TIMEOUT = elapsed_time > TIMEOUT
                        print("Elapsed time is {0}m".format(elapsed_time / 60.0))
                        if len(url)>0:
                            downloader = BaseDownloader('espa_downloads')                            
                            downloader.download(url)
                        #if os.path.exists(os.path.join(os.getcwd,'espa_downloads',url.split(os.sep)[-1][:-7])):
                            complete = True
                        
                        if not complete:
                            sleep(300)

    if l8_tiles:
        print("Download new data...")       
        #======Download data=========    
#        for download in client.download_order_gen(orderidNew):
#            print(download)
        for download in download_order_gen(orderidNew,auth):
            print(download)

def get_modis_lai(tiles,product,version,start_date,end_date,auth):    

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
        modisOgg.downloadsAllDay()
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
    
def geotiff_2envi():   
    #geotiffConvert = 'GeoTiff2ENVI'
    bands = ["blue","green","red","nir","swir1","swir2","cloud"]
    l8bands = ["sr_band2","sr_band3","sr_band4","sr_band5","sr_band6","sr_band7","cfmask"] 
    
    landsat_files = glob.glob(os.path.join(landsat_temp,"*_MTL.txt"))
    for i in range(len(landsat_files)):
        
        fn = landsat_files[i][:-8]
        meta = landsat_metadata(landsat_files[i])
        fstem = os.path.join(os.sep.join((fn.split(os.sep)[:-1])),meta.LANDSAT_SCENE_ID)

        for i in range(len(bands)):
            tifFile = fn+"_%s.tif" % l8bands[i]
            datFile = fstem+"_%s.%s.dat" % (l8bands[i],bands[i])
            #subprocess.call(["%s" % geotiffConvert ,"%s" % tifFile, "%s" % datFile])
            subprocess.call(["gdal_translate","-q", "-of", "ENVI","%s" % tifFile, "%s" % datFile ])
            os.rename("%s.hdr" % datFile[:-4],"%s.dat.hdr" % datFile[:-4])

def sample():    
    sample = 'lndlai_sample'
    bands = ["blue","green","red","nir","swir1","swir2","cloud"]
    l8bands = ["sr_band2","sr_band3","sr_band4","sr_band5","sr_band6","sr_band7","cfmask"] 
    
    landsat_files = glob.glob(os.path.join(landsat_temp,"*_MTL.txt"))
    
    for i in range(len(landsat_files)):
        #sceneID = landsat_files[i].split(os.sep)[-1][:-4]
        meta = landsat_metadata(landsat_files[i])
        sceneID = meta.LANDSAT_SCENE_ID
        
        # extract the Landsat doy and year
        d = meta.DATETIME_OBJ
        year = d.year
        ldoy = sceneID[13:16]
#        year = int(sceneID[9:13])
        # convert to date    
#        dd = datetime.datetime(year, 1, 1) + datetime.timedelta(int(ldoy) - 1)
#        date = '%d-%02d-%02d' % (dd.year,dd.month,dd.day)
        date = meta.DATE_ACQUIRED
        # find the 4 day MODIS doy prior to the Landsat doy
        mdoy = int((int((float(ldoy)-1.)/4.)*4.)+1)
        
        modFiles = glob.glob(os.path.join(modis_base,"MCD15A3H","MCD15A3H.A%s%s.*.hdf" % (year,mdoy)))

        #fstem = landsat_files[i][:-4]
        fn = landsat_files[i][:-8]
        fstem = os.path.join(os.sep.join((fn.split(os.sep)[:-1])),meta.LANDSAT_SCENE_ID)
        lai_path = landsat_LAI
        if not os.path.exists(lai_path):
            os.mkdir(lai_path)
        sam_file = os.path.join(lai_path,"SR_LAI.%s.%s.MCD15A3H_A%s%s.txt" %(date,sceneID,year,mdoy))
        
        for i in range(len(modFiles)):  
            fn = os.path.join(lai_path,"slai%d.inp" % i)
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
            #os.remove(os.path.join(lai_path,"slai%d.inp" % i))
            
def train():    
    cubist = 'cubist'
    landsat_files = glob.glob(os.path.join(landsat_LAI,"*.txt"))
    #======combine input data======================================
    df = pd.DataFrame(columns=['ulx','uly','blue',
        'green','red','nir','swir1','swir2','ndvi','ndwi','lai','weight','satFlag'])
    for i in range(len(landsat_files)):
        sam_file = landsat_files[i]
    
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
    
def compute():    
    lndbio ='lndlai_compute'
    bands = ["blue","green","red","nir","swir1","swir2","cloud"]
    l8bands = ["sr_band2","sr_band3","sr_band4","sr_band5","sr_band6","sr_band7","cfmask"] 
    
    landsat_files = glob.glob(os.path.join(landsat_temp,"*_MTL.txt"))
    for i in range(len(landsat_files)):
#        sceneID = landsat_files[i].split(os.sep)[-1][:-4]  
        meta = landsat_metadata(landsat_files[i])
        sceneID = meta.LANDSAT_SCENE_ID
        #fstem = landsat_files[i][:-4]       
        fn = landsat_files[i][:-8]
        fstem = os.path.join(os.sep.join((fn.split(os.sep)[:-1])),meta.LANDSAT_SCENE_ID)
        # create a folder for lai if it does not exist
        #lai_path = os.path.join(landsat_LAI,'%s' % sceneID[9:16])
        lai_path = os.path.join(landsat_LAI,'%s' % sceneID[3:9])
        if not os.path.exists(lai_path):
            os.mkdir(lai_path)
        
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
        shutil.move(laiFN,os.path.join(lai_path,"lndlai.%s.hdf" % sceneID))
        os.remove(fn)
    #=====CLEANING UP========
    filelist = [ f for f in os.listdir(landsat_LAI) if f.startswith("lndsr_modlai_samples") ]
    for f in filelist:
        os.remove(os.path.join(landsat_LAI,f))
    
def get_LAI(collection,loc,start_date,end_date,usgs_user,usgs_pass,earth_user,
            earth_pass,cloud):    
    
    #start Landsat order process
    get_landsat_data(collection,loc,start_date,end_date,("%s"% usgs_user,"%s"% usgs_pass),cloud)
    
    # move surface relectance files and estimate get LAI
    download_folder = os.path.join(base,'espa_downloads')
    folders_2move = glob.glob(os.path.join(download_folder ,'*'))

    for i in range(len(folders_2move)):
        inputFN = folders_2move[i]
        metFN = glob.glob(os.path.join(inputFN,'*MTL.txt'))[0]
        meta = landsat_metadata(metFN)
        #sceneID = (inputFN).split(os.sep)[-1].split('-')[0]
        sceneID = meta.LANDSAT_SCENE_ID
        scene = sceneID[3:9]
        folder = os.path.join(landsat_SR,scene)
        if not os.path.exists(folder):
            os.mkdir(folder)
    
        for filename in glob.glob(os.path.join(inputFN, '*.*')):
            shutil.copy(filename, folder)  
            os.symlink(os.path.join(folder,filename.split(os.sep)[-1]),
            os.path.join(landsat_temp,filename.split(os.sep)[-1]))
 
    if len(folders_2move)>0:
            #======Clean up folder===============================
            shutil.rmtree(download_folder)
    # find MODIS tiles that cover landsat scene
    # MODIS products   
    product = 'MCD15A3H'
    version = '006'
    [v,h] = latlon_2modis_tile(loc[0],loc[1])
    tiles = "h%02dv%02d" %(h,v)
    #tiles = 'h10v04,h10v05'
    
    # download MODIS LAI over the same area and time
    print("Downloading MODIS data...")
    get_modis_lai(tiles,product,version,start_date,end_date,("%s"% earth_user,"%s"% earth_pass))
    
    # Convert Landsat SR downloads to ENVI format
    # Note:  May be some warnings about unknown field - ignore.
    print("Converting Landsat SR to ENVI format...")
    geotiff_2envi()
    
    # Generate MODIS-Landsat samples for LAI computation
    print("Generating MODIS-Landsat samples...")
    sample()    
    
    # Compute Landsat LAI
    print("Computing Landsat LAI...")
    train()
    compute()    

def main():
    # Get time and location from user
    parser = argparse.ArgumentParser()
    parser.add_argument("lat", type=float, help="latitude")
    parser.add_argument("lon", type=float, help="longitude")
    parser.add_argument("start_date", type=str, help="Start date yyyy-mm-dd")
    parser.add_argument("end_date", type=str, help="Start date yyyy-mm-dd")
    parser.add_argument("cloud", type=int, help="cloud coverage")
    args = parser.parse_args()
      
    loc = [args.lat,args.lon] 
    start_date = args.start_date
    end_date = args.end_date
    cloud = args.cloud
    collection = 1
    
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
    
    get_LAI(collection,loc,start_date,end_date,usgs_user,usgs_pass,earth_user,
            earth_pass,cloud)

    print("All done with LAI")


if __name__ == "__main__":
    try:
        main()
    except (KeyboardInterrupt, pycurl.error):
        exit('Received Ctrl + C... Exiting! Bye.', 1)   
    