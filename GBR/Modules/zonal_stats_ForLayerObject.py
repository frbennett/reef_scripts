"""
Zonal Statistics
Vector-Raster Analysis

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
def zonal_stats(polygon_layer, IDattribute, data_raster, nodata_value=None):
    #print("Working with Unique ID: " + IDattribute)
    #data_ds = gdal.Open(raster_path, GA_ReadOnly)
    data_ds = data_raster
    #assert(data_ds)
    data_band = data_ds.GetRasterBand(1)
    data_geotransform = data_ds.GetGeoTransform()
    
    templateColCount = data_raster.RasterXSize
    templateRowCount = data_raster.RasterYSize

    if nodata_value:
        nodata_value = float(nodata_value)
        data_band.SetNoDataValue(nodata_value)

    #vds = ogr.Open(vector_path, GA_ReadOnly)  # TODO maybe open update if we want to write stats
    #assert(vds)
    #vlyr = vds.GetLayer(0)
    
    vlyr = polygon_layer
    
    print("Doing zonal stats on " + str(vlyr.GetFeatureCount()) + " features")

    global_src_extent = False

    # create an in-memory numpy array of the source raster data
    # covering the whole extent of the vector layer
    if global_src_extent:
        # use global source extent
        # useful only when disk IO or raster scanning inefficiencies are your limiting factor
        # advantage: reads raster data in one pass
        # disadvantage: large vector extents may have big memory requirements
        src_offset = bbox_to_pixel_offsets(data_geotransform, vlyr.GetExtent(), templateColCount, templateRowCount)
        src_array = data_band.ReadAsArray(*src_offset)

        # calculate new geotransform of the layer subset
        new_gt = (
            (data_geotransform[0] + (src_offset[0] * data_geotransform[1])),
            data_geotransform[1],
            0.0,
            (data_geotransform[3] + (src_offset[1] * data_geotransform[5])),
            0.0,
            data_geotransform[5]
        )

    mem_drv = ogr.GetDriverByName('Memory')
    driver = gdal.GetDriverByName('MEM')

    # Loop through vectors
    stats = []
    feat = vlyr.GetNextFeature()
    while feat is not None:

        if not global_src_extent:
            # use local source extent
            # fastest option when you have fast disks and well indexed raster (ie tiled Geotiff)
            # advantage: each feature uses the smallest raster chunk
            # disadvantage: lots of reads on the source raster
            src_offset = bbox_to_pixel_offsets(data_geotransform, feat.geometry().GetEnvelope(), templateColCount, templateRowCount)
            src_array = data_band.ReadAsArray(*src_offset)

            # calculate new geotransform of the feature subset
            new_gt = (
                (data_geotransform[0] + (src_offset[0] * data_geotransform[1])),
                data_geotransform[1],
                0.0,
                (data_geotransform[3] + (src_offset[1] * data_geotransform[5])),
                0.0,
                data_geotransform[5]
            )

        # Create a temporary vector layer in memory
        mem_ds = mem_drv.CreateDataSource('out')
        mem_layer = mem_ds.CreateLayer('poly', None, ogr.wkbPolygon)
        mem_layer.CreateFeature(feat.Clone())

        # Rasterize it
        rvds = driver.Create('', src_offset[2], src_offset[3], 1, gdal.GDT_Byte)
        rvds.SetGeoTransform(new_gt)
        gdal.RasterizeLayer(rvds, [1], mem_layer, burn_values=[1])
        #print("Attempting Rasterize for feature")
        #gdal.RasterizeLayer(rvds, [1], mem_layer, options=["ATTRIBUTE=" + IDattribute])
        rv_array = rvds.ReadAsArray()
        
        #print("Processing feature with UID: " + str(feat.GetField(IDattribute)))

        # Mask the source data array with our current feature
        # we take the logical_not to flip 0<->1 to get the correct mask effect
        # we also mask out nodata values explictly
        masked = np.ma.MaskedArray(
            src_array,
            mask=np.logical_or(
                src_array == nodata_value,
                np.logical_not(rv_array)
            )
        )

        feature_stats = {
            'min': float(masked.min()),
            'mean': float(masked.mean()),
            'max': float(masked.max()),
            'std': float(masked.std()),
            'sum': float(masked.sum()),
            'count': int(masked.count()),
            IDattribute: int(feat.GetField(IDattribute))}

        stats.append(feature_stats)
        
        #print("Result for feat UID: " + str(feat.GetField(IDattribute)) + " Mean: " + str(masked.mean()) + " Count: " + str(masked.count()))

        rvds = None
        mem_ds = None
        feat = vlyr.GetNextFeature()

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
