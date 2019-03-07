#!/usr/bin/env python2
# !/Applications/anaconda/envs/root3 python
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 13:50:04 2016

@author: mschull
"""
# python

# from .search import Search
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
import wget
import pycurl
from .landsatTools import landsat_metadata
import logging
from pymodis.downmodis import downModis

logging.getLogger("urllib3").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.DEBUG)

# The current URL hosting the ESPA interfaces has reached a stable version 1.0
host = 'https://espa.cr.usgs.gov/api/v1/'
TIMEOUT = 86400

base = os.getcwd()
cacheDir = os.path.abspath(os.path.join(base, "SATELLITE_DATA"))
Folders = folders(base)
modis_base = Folders['modis_base']
landsat_SR = Folders['landsat_SR']
landsat_LAI = Folders['landsat_LAI']
landsat_NDVI = Folders['landsat_NDVI']
landsat_Mask = Folders['landsat_Mask']
landsat_temp = os.path.join(landsat_SR, 'temp')
if not os.path.exists(landsat_temp):
    os.mkdir(landsat_temp)

def search(lat, lon, start_date, end_date, cloud, cacheDir, sat):
    columns = ['acquisitionDate', 'acquisitionDate', 'upperLeftCornerLatitude', 'upperLeftCornerLongitude',
               'lowerRightCornerLatitude', 'lowerRightCornerLongitude', 'cloudCover', 'sensor', 'LANDSAT_PRODUCT_ID']
    end = datetime.datetime.strptime(end_date, '%Y-%m-%d')
    # this is a landsat-util work around when it fails
    if sat == 7:
        metadataUrl = 'https://landsat.usgs.gov/landsat/metadata_service/bulk_metadata_files/LANDSAT_ETM_C1.csv'
    else:
        metadataUrl = 'https://landsat.usgs.gov/landsat/metadata_service/bulk_metadata_files/LANDSAT_8_C1.csv'

    fn = os.path.join(cacheDir, metadataUrl.split(os.sep)[-1])
    # looking to see if metadata CSV is available and if its up to the date needed
    if os.path.exists(fn):
        d = datetime.datetime.fromtimestamp(os.path.getmtime(fn))
        if (end.year > d.year) and (end.month > d.month) and (end.day > d.day):
            wget.download(metadataUrl, out=fn)
            df = pd.read_csv(fn, usecols=columns)
            df.to_csv(fn)
        df = pd.read_csv(fn)
        index = ((df.acquisitionDate >= start_date) & (df.acquisitionDate < end_date) & (
                df.upperLeftCornerLatitude > lat) & (df.upperLeftCornerLongitude < lon) & (
                         df.lowerRightCornerLatitude < lat) & (df.lowerRightCornerLongitude > lon) & (
                         df.cloudCover <= cloud) & (df.sensor == 'OLI_TIRS'))
        df = df[index]

    else:
        wget.download(metadataUrl, out=fn)
        df = pd.read_csv(fn, usecols=columns)
        df.to_csv(fn)
        index = ((df.acquisitionDate >= start_date) & (df.acquisitionDate < end_date) & (
                df.upperLeftCornerLatitude > lat) & (df.upperLeftCornerLongitude < lon) & (
                         df.lowerRightCornerLatitude < lat) & (df.lowerRightCornerLongitude > lon) & (
                         df.cloudCover <= cloud) & (df.sensor == 'OLI_TIRS'))
        df = df[index]

    return df


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


def find_already_downloaded(df, cache_dir):
    usgs_available = list(df.LANDSAT_PRODUCT_ID.values)
    # find sat
    sat = usgs_available[0].split("_")[0][-1]
    # find scenes
    scenes = [x.split("_")[2] for x in usgs_available]
    scenes = list(set(scenes))
    available_list = []
    for scene in scenes:
        path_to_search = os.path.join(cache_dir, 'L%s/%s/RAW_DATA/*MTL*' % (sat, scene))
        available = [os.path.basename(x) for x in
                     glob.glob(path_to_search)]
        available = [x[:-8] for x in available]
        available_list = available_list + available
    return intersection(usgs_available, available_list)


def find_not_processed(downloaded, cache_dir):
    """finds the files that are downloaded but still need to process LAI data"""
    # find sat
    sat = downloaded[0].split("_")[0][-1]
    # find scenes
    scenes = [x.split("_")[2] for x in downloaded]
    scenes = list(set(scenes))
    available_list = []
    for scene in scenes:
        path_to_search = os.path.join(cache_dir, 'L%s/%s/LAI/*_lai.tif' % (sat, scene))
        available = [os.path.basename(x) for x in
                     glob.glob(path_to_search)]
        available = [x[:-8] for x in available]
        available_list = available_list + available
    for x in available_list:
        if x in downloaded:
            downloaded.remove(x)
    return downloaded


def get_modis_lai(tiles, product, version, start_date, end_date, auth, cacheDir):
    startdd = datetime.datetime.strptime(start_date, '%Y-%m-%d')
    enddd = datetime.datetime.strptime(end_date, '%Y-%m-%d')
    numDays = (enddd - startdd).days + 1

    laidates = np.array(range(1, 366, 4))

    if product.startswith('MCD'):
        folder = "MOTA"
    elif product.startswith('MOD'):
        folder = "MOLT"
    else:
        folder = "MOTA"
    product_path = os.path.join(cacheDir, product)
    if not os.path.exists(product_path):
        os.makedirs(product_path)


    modisOgg = downModis(url="https://e4ftl01.cr.usgs.gov", destinationFolder=product_path,
                         user=auth[0], password=auth[1], tiles=tiles, path=folder,
                         product="%s.%s" % (product, version), today=start_date, enddate=end_date)
    # connect to http or ftp
    modisOgg.connect()
    if modisOgg.nconnection <= 20:
        # download data
        filenames = []
        for i in range(numDays):
            dd = startdd + datetime.timedelta(days=i)
            year = dd.year
            doy = (dd - datetime.datetime(year, 1, 1, 1, 1)).days
            rday = laidates[laidates >= doy][0]
            if (rday == doy):
                dd = datetime.datetime(year-1 , 12, 31, 1, 1) + datetime.timedelta(days=rday)
                year = dd.year
                month = dd.month
                day = dd.day
                dayIn = '%d.%02d.%02d' % (year, month, day)
                listNewFile = modisOgg.getFilesList(dayIn)
                listFilesDown = modisOgg.checkDataExist(listNewFile)
                if (len(listFilesDown) < 1):
                    continue
                file2download = filter(lambda x: x.endswith('.hdf'), listFilesDown)
                if np.shape(file2download)[0] >0:
                    filenames.append(os.path.join(product_path, file2download[0]))
                    #                    filenames.append(os.path.join(product_path,listFilesDown[0]))
                    modisOgg.dayDownload(dayIn, listFilesDown)
    #        modisOgg.downloadsAllDay()
    else:
        print("A problem with the connection occured")
    return None


def latlon_2modis_tile(lat, lon):
    # reference: https://code.env.duke.edu/projects/mget/wiki/SinusoidalMODIS
    p_modis_grid = Proj('+proj=sinu +R=6371007.181 +nadgrids=@null +wktext')
    x, y = p_modis_grid(lon, lat)
    # or the inverse, from x, y to lon, lat
    lon, lat = p_modis_grid(x, y, inverse=True)
    tileWidth = 1111950.5196666666
    ulx = -20015109.354
    uly = -10007554.677
    H = (x - ulx) / tileWidth
    V = 18 - ((y - uly) / tileWidth)
    return int(V), int(H)

def get_modis_tiles(ulLat, ulLon, lrLat, lrLon, modisCacheDir):
    tileList = set()
    for i in range(int(lrLat), int(ulLat)):
        for j in range(int(ulLon), int(lrLon)):
            v, h = latlon_2modis_tile(i, j)
            tileList.add(v * 1000 + h)
    tileList = list(tileList)
    modis_tiles = []
    for i in range(len(tileList)):
        v = tileList[i] / 1000
        h = tileList[i] - v * 1000
        mod_tile = "h%02dv%02d" % (h, v)
        fns = glob.glob(os.path.join(modisCacheDir, '*%s*.hdf' % mod_tile))
        if np.shape(fns)[0] == 0:
            modis_tiles.append(mod_tile)
    return modis_tiles

def geotiff_2envi(paths, productIDs):
    # geotiffConvert = 'GeoTiff2ENVI'
    bands = ["blue", "green", "red", "nir", "swir1", "swir2", "cloud"]
    l8bands = ["sr_band2", "sr_band3", "sr_band4", "sr_band5", "sr_band6", "sr_band7", "pixel_qa"]
    count = 0
    for productID in productIDs:
        fstem_in = os.path.join(paths[count], productID)
        fstem = os.path.join(landsat_temp, productID)
        count += 1
        for i in range(len(bands)):
            tifFile = fstem_in + "_%s.tif" % l8bands[i]
            datFile = fstem + "_%s.%s.dat" % (l8bands[i], bands[i])
            # subprocess.call(["%s" % geotiffConvert ,"%s" % tifFile, "%s" % datFile])
            subprocess.call(["gdal_translate", "-q", "-of", "ENVI", "%s" % tifFile, "%s" % datFile])
            os.rename("%s.hdr" % datFile[:-4], "%s.dat.hdr" % datFile[:-4])


def sample(paths, productIDs, MODIS_product, modisCacheDir):
    sample = 'lndlai_sample'
    convert = 'lndqa2fmask'
    bands = ["blue", "green", "red", "nir", "swir1", "swir2", "cloud"]
    l8bands = ["sr_band2", "sr_band3", "sr_band4", "sr_band5", "sr_band6", "sr_band7", "cfmask"]

    count = 0
    for productID in productIDs:
        fstem = os.path.join(paths[count], productID)
        count += 1
        meta = landsat_metadata(fstem + "_MTL.txt")
        sceneID = meta.LANDSAT_SCENE_ID

        # extract the Landsat doy and year
        d = meta.DATETIME_OBJ
        year = d.year
        ldoy = sceneID[13:16]
        date = meta.DATE_ACQUIRED
        # find the 4 day MODIS doy prior to the Landsat doy
        mdoy = int((int((float(ldoy) - 1.) / 4.) * 4.) + 1)

        modFiles = glob.glob(
            os.path.join(modisCacheDir, "%s" % MODIS_product, "%s.A%d%03d.*.hdf" % (MODIS_product, year, mdoy)))
        sam_file = os.path.join(landsat_LAI, "SR_LAI.%s.%s.%s_A%d%03d.txt" % (date, sceneID, MODIS_product, year, mdoy))
        # ====convert the qa to cfmask=====
        fstem = os.path.join(landsat_temp, productID)
        datFile_qa = fstem + "_%s.%s.dat" % ("pixel_qa", bands[6])
        datFile_cfmask = fstem + "_%s.%s.dat" % (l8bands[6], bands[6])
        subprocess.call(["%s" % convert, "-lndsr", "%s" % datFile_qa, "-cmask", "%s" % datFile_cfmask])
        for i in range(len(modFiles)):
            fn = os.path.join(landsat_LAI, "slai%d.inp" % i)
            file = open(fn, "w")
            file.write("LANDSAT_BASE_BLUE = %s_%s.%s.dat\n" % (fstem, l8bands[0], bands[0]))
            file.write("LANDSAT_BASE_GREEN = %s_%s.%s.dat\n" % (fstem, l8bands[1], bands[1]))
            file.write("LANDSAT_BASE_RED = %s_%s.%s.dat\n" % (fstem, l8bands[2], bands[2]))
            file.write("LANDSAT_BASE_NIR = %s_%s.%s.dat\n" % (fstem, l8bands[3], bands[3]))
            file.write("LANDSAT_BASE_SWIR1 = %s_%s.%s.dat\n" % (fstem, l8bands[4], bands[4]))
            file.write("LANDSAT_BASE_SWIR2 = %s_%s.%s.dat\n" % (fstem, l8bands[5], bands[5]))
            file.write("LANDSAT_BASE_CLOUD = %s_%s.%s.dat\n" % (fstem, l8bands[6], bands[6]))
            file.write("MODIS_BASE_FILE = %s\n" % modFiles[i])
            file.write("SAMPLE_FILE_OUT = %s\n" % sam_file)
            file.write("PURE_SAMPLE_TH = 0.2\n")
            file.close()

            subprocess.call(["%s" % sample, "%s" % fn])


def train(paths, productIDs, MODIS_product):
    cubist = 'cubist'
    # ======combine input data======================================
    df = pd.DataFrame(columns=['ulx', 'uly', 'blue',
                               'green', 'red', 'nir', 'swir1', 'swir2', 'ndvi', 'ndwi', 'lai', 'weight', 'satFlag'])
    count = 0
    for productID in productIDs:
        fstem = os.path.join(paths[count], productID)
        count += 1
        meta = landsat_metadata(fstem + "_MTL.txt")
        sceneID = meta.LANDSAT_SCENE_ID
        # extract the Landsat doy and year
        d = meta.DATETIME_OBJ
        year = d.year
        ldoy = sceneID[13:16]
        date = meta.DATE_ACQUIRED
        # find the 4 day MODIS doy prior to the Landsat doy
        mdoy = int((int((float(ldoy) - 1.) / 4.) * 4.) + 1)
        sam_file = os.path.join(landsat_LAI, "SR_LAI.%s.%s.%s_A%d%03d.txt" % (date, sceneID, MODIS_product, year, mdoy))

        df = df.append(pd.read_csv(sam_file, delim_whitespace=True, names=['ulx', 'uly', 'blue',
                                                                           'green', 'red', 'nir', 'swir1', 'swir2',
                                                                           'ndvi', 'ndwi', 'lai', 'weight', 'satFlag']),
                       ignore_index=True)

    # =====create filestem.data====================================
    df = df[(df.satFlag == 'N')]
    df = df.sort_values(by='weight')
    start_date = '200'
    end_date = '300'
    filestem = os.path.join(landsat_LAI, "lndsr_modlai_samples.combined_%s-%s" % (start_date, end_date))
    df.to_csv(os.path.join(landsat_LAI, filestem + ".data"), columns=['blue', 'green', 'red',
                                                                      'nir', 'swir1', 'swir2', 'ndvi', 'ndwi', 'lai',
                                                                      'weight'], header=None,
              index=None, mode='w', sep="\t", encoding='utf-8')

    # ====create filestem.names====================================
    fn = os.path.join(landsat_LAI, "%s.names" % filestem)
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
    subprocess.call(["%s" % cubist, "-f", "%s" % filestem, "-r", "%d" % nrules, "-u"])


def compute(paths, productIDs, sat, cacheDir):
    lndbio = 'lndlai_compute'
    bands = ["blue", "green", "red", "nir", "swir1", "swir2", "cloud"]
    l8bands = ["sr_band2", "sr_band3", "sr_band4", "sr_band5", "sr_band6", "sr_band7", "cfmask"]

    count = 0
    lai_fns = []
    ndvi_fns = []
    mask_fns = []
    for productID in productIDs:
        fstem = os.path.join(paths[count], productID)
        count += 1
        meta = landsat_metadata(fstem + "_MTL.txt")
        sceneID = meta.LANDSAT_SCENE_ID
        scene = sceneID[3:9]
        folder = os.path.join(cacheDir, "L%d" % sat, scene)
        lai_path = os.path.join(folder, "LAI")
        if not os.path.exists(lai_path):
            os.mkdir(lai_path)
        ndvi_path = os.path.join(folder, "NDVI")
        if not os.path.exists(ndvi_path):
            os.mkdir(ndvi_path)
        cfmask_path = os.path.join(folder, "CF_MASK")
        if not os.path.exists(cfmask_path):
            os.mkdir(cfmask_path)
        start_date = '200'
        end_date = '300'
        filestem = os.path.join(landsat_LAI, "lndsr_modlai_samples.combined_%s-%s" % (start_date, end_date))
        laiFN = os.path.join(landsat_LAI, "lndlai.%s.hdf" % sceneID)
        fn = os.path.join(landsat_LAI, "compute_lai.inp")
        fstem = os.path.join(landsat_temp, productID)
        file = open(fn, "w")
        file.write("LANDSAT_BASE_BLUE = %s_%s.%s.dat\n" % (fstem, l8bands[0], bands[0]))
        file.write("LANDSAT_BASE_GREEN = %s_%s.%s.dat\n" % (fstem, l8bands[1], bands[1]))
        file.write("LANDSAT_BASE_RED = %s_%s.%s.dat\n" % (fstem, l8bands[2], bands[2]))
        file.write("LANDSAT_BASE_NIR = %s_%s.%s.dat\n" % (fstem, l8bands[3], bands[3]))
        file.write("LANDSAT_BASE_SWIR1 = %s_%s.%s.dat\n" % (fstem, l8bands[4], bands[4]))
        file.write("LANDSAT_BASE_SWIR2 = %s_%s.%s.dat\n" % (fstem, l8bands[5], bands[5]))
        file.write("LANDSAT_BASE_CLOUD = %s_%s.%s.dat\n" % (fstem, l8bands[6], bands[6]))
        file.write("LANDSAT_ANC_FILE = %s\n" % filestem)
        file.write("BIOPHYSICS_PARA_FILE_OUT = %s\n" % laiFN)
        file.close()

        subprocess.call(["%s" % lndbio, "%s" % fn])
        # ====convert to geotiff=========
        outlaifn = os.path.join(lai_path, '%s_lai.tif' % sceneID)
        lai_fns.append(outlaifn)
        outndvifn = os.path.join(ndvi_path, '%s_ndvi.tif' % sceneID)
        ndvi_fns.append(outndvifn)
        outcfmaskfn = os.path.join(cfmask_path, '%s_Mask.tif' % sceneID)
        mask_fns.append(outcfmaskfn)
        subprocess.call(["gdal_translate", 'HDF4_EOS:EOS_GRID:"%s":LANDSAT:LAI' % laiFN, "%s" % outlaifn])
        subprocess.call(["gdal_translate", 'HDF4_EOS:EOS_GRID:"%s":LANDSAT:NDVI' % laiFN, "%s" % outndvifn])
        subprocess.call(["gdal_translate", 'HDF4_EOS:EOS_GRID:"%s":LANDSAT:cfmask' % laiFN, "%s" % outcfmaskfn])
        os.remove(fn)

    # =====CLEANING UP========
    filelist = [f for f in os.listdir(landsat_LAI) if f.startswith("lndsr_modlai_samples")]
    for f in filelist:
        os.remove(os.path.join(landsat_LAI, f))
    return lai_fns, ndvi_fns, mask_fns


def get_LAI(loc, start_date, end_date, earth_user, earth_pass, cloud, sat, cacheDir):
    landsatCacheDir = os.path.join(cacheDir, "LANDSAT")
    modisCacheDir = os.path.join(cacheDir, "MODIS")
    # find MODIS tiles that cover landsat scene
    # MODIS products   
    MODIS_product = 'MCD15A3H'
    MODIS_path = os.path.join(modisCacheDir, MODIS_product)
    version = '006'
    # ====search for available data=============================================
    output_df = search(loc[0], loc[1], start_date, end_date, cloud, landsatCacheDir, sat)
    downloaded = find_already_downloaded(output_df, landsatCacheDir)
    productIDs = find_not_processed(downloaded, landsatCacheDir)
    paths = []
    for productID in productIDs:
        sat_str = productID.split("_")[0][-1]
        scene = productID.split("_")[2]
        paths = paths + [os.path.join(landsatCacheDir, 'L%s/%s/RAW_DATA/' % (sat_str, scene))]
    # paths = list(set(paths))

    # download MODIS LAI over the same area and time
    print("Downloading MODIS data...")
    # find the MODIS tiles for Landsat scenes
    tiles = []
    for i in range(len(output_df)):
        tiles = tiles + get_modis_tiles(output_df.iloc[i].upperLeftCornerLatitude,
                                        output_df.iloc[i].upperLeftCornerLongitude,
                                        output_df.iloc[i].lowerRightCornerLatitude,
                                        output_df.iloc[i].lowerRightCornerLongitude, modisCacheDir)
    tiles = list(set(tiles))

    get_modis_lai(tiles, MODIS_product, version, start_date, end_date,
                                ("%s" % earth_user, "%s" % earth_pass), modisCacheDir)

    # Convert Landsat SR downloads to ENVI format
    # Note:  May be some warnings about unknown field - ignore.
    if len(productIDs) > 0:
        print("Converting Landsat SR to ENVI format...")
        geotiff_2envi(paths, productIDs)

        # Generate MODIS-Landsat samples for LAI computation
        print("Generating MODIS-Landsat samples...")
        sample(paths, productIDs, MODIS_product, modisCacheDir)

        # Compute Landsat LAI
        print("Computing Landsat LAI...")
        train(paths, productIDs, MODIS_product)
        lai_fns, ndvi_fns, mask_fns = compute(paths, productIDs, sat, landsatCacheDir)
    else:
        print("Nothing to process!!")


def main():
    # Get time and location from user
    parser = argparse.ArgumentParser()
    parser.add_argument("lat", type=float, help="latitude")
    parser.add_argument("lon", type=float, help="longitude")
    parser.add_argument("start_date", type=str, help="Start date yyyy-mm-dd")
    parser.add_argument("end_date", type=str, help="Start date yyyy-mm-dd")
    parser.add_argument("cloud", type=int, help="cloud coverage")
    parser.add_argument('-s', '--sat', nargs='?', type=int, default=8,
                        help='which landsat to search or download, i.e. Landsat 8 = 8')
    args = parser.parse_args()

    loc = [args.lat, args.lon]
    start_date = args.start_date
    end_date = args.end_date
    cloud = args.cloud
    sat = args.sat

    # =====USGS credentials===============

    # =====earthData credentials===============
    earth_user = str(getpass.getpass(prompt="earth login username:"))
    if keyring.get_password("nasa", earth_user) == None:
        earth_pass = str(getpass.getpass(prompt="earth login password:"))
        keyring.set_password("nasa", earth_user, earth_pass)
    else:
        earth_pass = str(keyring.get_password("nasa", earth_user))

    get_LAI(loc, start_date, end_date, earth_user, earth_pass, cloud, sat, cacheDir)

    print("All done with LAI")


if __name__ == "__main__":
    try:
        main()
    except (KeyboardInterrupt, pycurl.error):
        exit('Received Ctrl + C... Exiting! Bye.', 1)
