#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 11:38:31 2017

@author: mschull
"""
import os
import numpy as np
import pandas as pd
import requests
import json
from datetime import datetime
import wget
from osgeo import gdal,osr

# The current URL hosting the ESPA interfaces has reached a stable version 1.0
host = 'https://espa.cr.usgs.gov/api/v1/'
#from numba import jit

def folders(base):
    database = os.path.join(base,'data')
    landsat_database = os.path.join(database,'Landsat-8')
    landsat_SR = os.path.join(landsat_database,'SR')
    if not os.path.exists(landsat_SR):
        os.makedirs(landsat_SR)
    landsat_LAI = os.path.join(landsat_database,'LAI')
    if not os.path.exists(landsat_LAI):
        os.makedirs(landsat_LAI)
    landsat_NDVI = os.path.join(landsat_database,'NDVI')
    if not os.path.exists(landsat_NDVI):
        os.makedirs(landsat_NDVI)
    landsat_Mask = os.path.join(landsat_database,'Mask')
    if not os.path.exists(landsat_Mask):
        os.makedirs(landsat_Mask)
    modis_base = os.path.join(base,'data','MODIS')
    if not os.path.exists(modis_base):
        os.mkdir(modis_base)
    out = {'landsat_SR':landsat_SR,'modis_base':modis_base,'landsat_LAI':landsat_LAI,
           'landsat_NDVI':landsat_NDVI,'landsat_Mask':landsat_Mask}
    return out


def search(collection,lat,lon,start_date,end_date,cloud):
    end = datetime.strptime(end_date, '%Y-%m-%d')
    # this is a landsat-util work around when it fails
    if collection==0:
        metadataUrl = 'https://landsat.usgs.gov/landsat/metadata_service/bulk_metadata_files/LANDSAT_8.csv'
    else:
        metadataUrl = 'https://landsat.usgs.gov/landsat/metadata_service/bulk_metadata_files/LANDSAT_8_C1.csv'
    fn  = metadataUrl.split(os.sep)[-1]
    # looking to see if metadata CSV is available and if its up to the date needed
    if os.path.exists(fn):
        d = datetime.fromtimestamp(os.path.getmtime(fn))
        if ((end.year>d.year) and (end.month>d.month) and (end.day>d.day)):
            wget.download(metadataUrl)
    else:
        wget.download(metadataUrl)
        
    metadata= pd.read_csv(fn)
    if collection==0:
        output = metadata[(metadata.acquisitionDate >= start_date) & (metadata.acquisitionDate < end_date) & 
             (metadata.upperLeftCornerLatitude > lat ) & (metadata.upperLeftCornerLongitude < lon )& 
             (metadata.lowerRightCornerLatitude < lat ) & (metadata.lowerRightCornerLongitude > lon)  & 
             (metadata.cloudCoverFull <= cloud)].sceneID
    else:
        output = metadata[(metadata.acquisitionDate >= start_date) & (metadata.acquisitionDate < end_date) & 
             (metadata.upperLeftCornerLatitude > lat ) & (metadata.upperLeftCornerLongitude < lon )& 
             (metadata.lowerRightCornerLatitude < lat ) & (metadata.lowerRightCornerLongitude > lon)  & 
             (metadata.cloudCoverFull <= cloud)].LANDSAT_PRODUCT_ID
    return output.values

def check_order_cache(auth):
    # This program checks the last month of orders from ESPA
    username = auth[0]
    password = auth[1]
    def api_request(endpoint, verb='get', json=None, uauth=None):
        """
        Here we can see how easy it is to handle calls to a REST API that uses JSON
        """
        auth_tup = uauth if uauth else (username, password)
        response = getattr(requests, verb)(host + endpoint, auth=auth_tup, json=json)
        return response.json()
    
    def espa_api(endpoint, verb='get', body=None, uauth=None):
        """ Suggested simple way to interact with the ESPA JSON REST API """
        auth_tup = uauth if uauth else (username, password)
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
    
#    usr = api_request('user')
    
#    order_list = api_request('list-orders/%s' % usr['email'])
    filters = {"status": ["complete", "ordered"]}  # Here, we ignore any purged orders
    order_list = espa_api('list-orders', body=filters)
    orderID =[]
    fName = []
    order_status=[]
#    for i in range(len(order_list['orders'])):
#        orderid = order_list['orders'][i]
    for i in range(len(order_list)):
        orderid = order_list[i]
#        date = orderid.split('-')[1]
#        if len(date)>8:
#            continue
#        year = int(date[4:])
#        day = int(date[2:4])
#        month = int(date[:2])
#        dOrder = datetime(year,month,day)
#        delt_date = (datetime.now()-dOrder).days
#        if delt_date>10:
#            continue
#        resp = api_request('item-status/{0}'.format(orderid))
        resp = espa_api('item-status/{0}'.format(orderid))
        ddd = json.loads(json.dumps(resp))
#        if not ddd['orderid']['%s' % orderid][0]['status']=='purged':
        for j in range(len(ddd['%s' % orderid])):
            fname = ddd['%s' % orderid][j]['name']
            status = ddd['%s' % orderid][j]['status']
            orderID.append(orderid)
            fName.append(fname)
            order_status.append(status)
                
    output = {'orderid':orderID,'productID':fName,'status':order_status}
    outDF = pd.DataFrame(output)  
    
    return outDF

def writeArray2Tiff(data,res,UL,inProjection,outfile,outFormat):

    xres = res[0]
    yres = res[1]

    ysize = data.shape[0]
    xsize = data.shape[1]

    ulx = UL[0] #- (xres / 2.)
    uly = UL[1]# - (yres / 2.)
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(outfile, xsize, ysize, 1, outFormat)
    #ds = driver.Create(outfile, xsize, ysize, 1, gdal.GDT_Int16)
    
    srs = osr.SpatialReference()
    
    if isinstance(inProjection, basestring):        
        srs.ImportFromProj4(inProjection)
    else:
        srs.ImportFromEPSG(inProjection)
        
    ds.SetProjection(srs.ExportToWkt())
    
    gt = [ulx, xres, 0, uly, 0, -yres ]
    ds.SetGeoTransform(gt)
    
    ds.GetRasterBand(1).WriteArray(data)
    #ds = None
    ds.FlushCache() 
# helper function
def _test_outside(testx, lower, upper):
    """
    True if testx, or any element of it is outside [lower, upper].

    Both lower bound and upper bound included
    Input: Integer or floating point scalar or Numpy array.
    """
    test = np.array(testx)
    return np.any(test < lower) or np.any(test > upper)

# custom exception
class RasterError(Exception):
    """Custom exception for errors during raster processing in Pygaarst"""
    pass