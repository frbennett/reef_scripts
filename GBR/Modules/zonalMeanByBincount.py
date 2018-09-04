"""
Zonal Mean using Numpy Array representations of Rasters

Rasters as Arrays MUST have been made the same size
Maked Arrays don't seem to work
Use the NoData as an index that will be FIRST processed,
as the results will come back in the order of index

Usage:
  zonalMeanByBincount.py RASTER as Numpy Array, RASTER as Numpy Array
  zonalMeanByBincount.py -h | --help
  zonalMeanByBincount.py --version

Options:
  -h --help     Show this screen.
  --version     Show version.
"""

import numpy as np

def zonalMeanByBincount(zoneRasterAsNumpyArray, dataRasterAsNumpyArray):
	
	#This will flatten the input zones array, the NoData val will be included as a value!
	unique_labels, zoneLabels = np.unique(zoneRasterAsNumpyArray, return_inverse=True)
	
	#flatten the data raster
	newData = dataRasterAsNumpyArray.ravel()
	#Convert NoData to NaN
	newData = newData.filled(np.nan)
	
	#Identify just the components of the data array to retain
	keep = ~np.isnan(newData)
	
	counts = np.bincount(zoneLabels[keep])
	sums = np.bincount(zoneLabels[keep], weights=newData[keep])
	return sums/np.asanyarray(counts).astype(np.float)
	
