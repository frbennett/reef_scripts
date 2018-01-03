"""
Zonal Statistics
Raster-Raster Analysis

Copyright 2013 Matthew Perry

Usage:
  zonal_stats.py ogr.osgeo.Layer String RASTER
  zonal_stats.py -h | --help
  zonal_stats.py --version

Options:
  -h --help     Show this screen.
  --version     Show version.
"""
from osgeo import gdal, ogr
from osgeo.gdalconst import *
import numpy as np
import sys
gdal.PushErrorHandler('CPLQuietErrorHandler')


def bbox_to_pixel_offsets(gt, bbox, colCount, rowCount):
    originX = gt[0]
    originY = gt[3]
    pixel_width = gt[1]
    pixel_height = gt[5]
    
    
    x1 = int((bbox[0] - originX) / pixel_width)
    x2 = int((bbox[1] - originX) / pixel_width) + 1

    y1 = int((bbox[3] - originY) / pixel_height)
    y2 = int((bbox[2] - originY) / pixel_height) + 1
    
    xsize = x2 - x1
    ysize = y2 - y1
    
    # Need to manipulate xsize and ysize if our polygon happens to have
    # an envelope that coincides exatcly with our template raster's extent
    if (x1 + xsize) > colCount:
    		xsize = colCount - x1

    if (y1 + ysize) > rowCount:
    		ysize = rowCount - y1    
    
    return (x1, y1, xsize, ysize)


#def zonal_stats(vector_path, raster_path, nodata_value=None, global_src_extent=False):
#def zonal_stats(polygon_layer, IDattribute, data_raster, nodata_value=None):
def zonal_stats(zone_dataset, data_raster):
		# assume that zone_dataset and data_raster are 100% compatible for noe
    #print("Working with Unique ID: " + IDattribute)
    #data_ds = gdal.Open(raster_path, GA_ReadOnly)
    zone_ds = zone_dataset
    data_ds = data_raster
    #assert(data_ds)
    data_band = data_ds.GetRasterBand(1)
    data_geotransform = data_ds.GetGeoTransform()
    
    ZoneNoDataVal = zone_ds.GetRasterBand(1).GetNoDataValue()
    templateColCount = zone_dataset.RasterXSize
    templateRowCount = zone_dataset.RasterYSize
    
    nodata_value = data_ds.GetRasterBand(1).GetNoDataValue()
    
    print("DataRaster NoDataValue: " + str(nodata_value))
    print("ZoneRaster NoDataValue: " + str(ZoneNoDataVal))

    data_array = data_band.ReadAsArray()
    
    DataInputNoDataMask = data_array == nodata_value

    mem_drv = ogr.GetDriverByName('Memory')
    driver = gdal.GetDriverByName('MEM')

    # Loop through uniq Values
    zonesRasterArray = zone_ds.ReadAsArray()
    uniq = np.unique(zonesRasterArray)
    
    print("Doing zonal stats on " + str(len(uniq) + 1) + " unique zones")
    
    stats = []
    #feat = vlyr.GetNextFeature()
    
    #while feat is not None:
    for zid in uniq:
    
        if zid == ZoneNoDataVal:
            print("Ignoring value: " + str(zid))
            continue
				        
        #print("Processing feature with UID: " + str(feat.GetField(IDattribute)))
				
        # Mask the source data array with our current feature
        # we take the logical_not to flip 0<->1 to get the correct mask effect
        # we also mask out nodata values explictly
        
        #dataMaskedForFu = data_array.copy()
        #dataMaskedForFu[zoneCellsToIgnore] = nodata_value;
                
        # This gets the right cell count per FU, but that mean/min/max are wrong... min always -9999, max always -9999
        masked = np.ma.masked_where(np.not_equal(zonesRasterArray, zid), data_array)
        
        
        
        
        #This gets the appropriate cell count for the FU, but seems to include -9999 cells, thus stuffing up the mean
        #Interesting the min & max always exactly -9999, despite slightly different means
        #masked = np.ma.MaskedArray(
        #    data_array,
        #    mask=np.logical_not(
        #        zonesMasked
        #    )
        #)
        
        #This version uses all Kfactor cells (16309852 of them), doesn't consider the FU footprint
        #masked = np.ma.MaskedArray(
        #    data_array,
        #    mask=np.logical_or(
        #        data_array == nodata_value,
        #        zonesRasterArray == zid
        #    )
        #)
        
        #This version uses all Kfactor cells (16309852 of them), doesn't consider the FU footprint
        #masked = np.ma.MaskedArray(
        #    data_array,
        #    mask=np.logical_or(
        #        data_array == nodata_value,
        #        zonesMasked
        #    )
        #)
        
        #This yields the 'converting a masked element to nan' error, but then also gives no stats (presumably all masks are made up of no cells)
        #masked = np.ma.MaskedArray(
        #    data_array,
        #    mask=np.logical_or(
        #        data_array == nodata_value,
        #        np.logical_not(zonesMasked)
        #    )
        #)
				
        feature_stats = {
            'min': float(masked.min()),
            'mean': float(masked.mean()),
            'max': float(masked.max()),
            'std': float(masked.std()),
            'sum': float(masked.sum()),
            'count': int(masked.count()),
            'zoneid': int(zid)}
				
        stats.append(feature_stats)
        
        print("Result for feat ZID: " + str(zid) + " Mean: " + str(masked.mean()) + " Count: " + str(masked.count()) + " Min: " + str(masked.min()) + " Max: " + str(masked.max()) + " Sum: " + str(masked.sum()))
				
        zonesMasked = None
        #rvds = None
        #mem_ds = None
        #feat = vlyr.GetNextFeature()

    vds = None
    data_ds = None
    
    print("Length of stats: " + str(len(stats)))
    
    return stats


if __name__ == "__main__":
    opts = {'LAYER': sys.argv[1], 'STRING': sys.argv[2], 'RASTER': sys.argv[3]}
    stats = zonal_stats(opts['LAYER'], opts['STRING'], opts['RASTER'])

    try:
        from pandas import DataFrame
        print(DataFrame(stats))
    except ImportError:
        import json
        print(json.dumps(stats, indent=2))
