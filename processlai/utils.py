#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 11:38:31 2017

@author: mschull
"""
import os
import numpy as np
from osgeo import gdal,osr

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