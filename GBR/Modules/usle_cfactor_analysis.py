"""
Rob's USLE cfactor analysis, to enable paralelisation
"""

import pandas
import numpy
from osgeo import gdal, gdalconst
from osgeo import ogr


def analyseThisCFactorRaster(theCFactorRas, theZonesRasterArray, ZonesRowCount, ZonesColCount, ZoneND, theDefaultC, theReqValidProp):
    # do pixel by pixel analysis
    # zones raster of intersected FUs will exist, and KLSFines raster will be co-located with this
    # send back a dictionary that also includes the valid/invalid cFactor numbers for each FU
    # should be able to access rasters etc without passing them into this function
    # same for the dictionary of results, it could already exist
    
    # Will need to watch the Y min/max interpretation for our Aus Albers stuff
    
    zoneIDStr = "zoneID"
    zoneCountStr = "zoneCount"
    validKLSCSum = "validKLSCSum"
    allKLSCSum = "allKLSCSum"
    validKLSCFinesSum = "validKLSCFinesSum"
    allKLSCFinesSum = "allKLSCFinesSum"
    validDataCountStr = "dataCount"
    
    KLSCvalStr = "KLSC"
    KLSCFinesValStr = "KLSCFines"
    
    cFactRas = gdal.Open(theCFactorRas)
    cFact_GeoTrans = cFactRas.GetGeoTransform()
    print(cFact_GeoTrans)
    cFactND = cFactRas.GetRasterBand(1).GetNoDataValue
    cFactNumpy = cFactRas.ReadAsArray()
    
    cFactRows = cFactRas.RasterYSize
    cFactCols = cFactRas.RasterXSize
    
    zGT = theZones.GetGeoTransform()
    
    #theZonesRasterArray = theZones.ReadAsArray()
    

    statsDict = {}
    statsDictToReturn = {}
    uniq = numpy.unique(theZonesRasterArray)
    for num in uniq:
        #using masked, so this check for NoData should be redundant
        #print("trying this UID: " + str(num))
        if num == NoData_value:
            continue
        
        statsDict[num] = {zoneCountStr: int(0), validKLSCSum: float(0), validKLSCFinesSum: float(0),
                          allKLSCSum: float(0), allKLSCFinesSum: float(0), validDataCountStr: int(0)}

    rowid = 0
    while rowid < ZonesRowCount:
    #while rowid < theZones.RasterYSize:
        #reset colid
        colid = 0
        if rowid % 500 == 0:
            print("Starting row: " + str(rowid))
        
        while colid < ZonesColCount:
        #while colid < theZones.RasterXSize:
            #ffff
            zoneVal = theZonesRasterArray[rowid, colid]
            
            if zoneVal == ZoneND:
                colid+=1
                continue
            
            #we have explicitly ensured that KLSFines matches our pixels
            #so just need to get the requisite cFactor row/col
            
            zoneCoords = centreCoordsForPixel(zGT, rowid, colid)
            cFactPix = pixelForCentreCoords(cFact_GeoTrans, zoneCoords[0], zoneCoords[1])
            
            #if cFactPix[0] >= cFactRows:
            #    print("Exceeded row count at " + str(cFactPix[0]) + " for coords: " + str(zoneCoords)
            #          + " where zone pixels were: " + str(rowid) + ", " + str(colid))
            #    SystemExit(0)
            
            #if cFactPix[1] >= cFactCols:
            #    print("Exceeded col count " + str(cFactPix[1]) + " for coords: " + str(zoneCoords)
            #          + " where zone pixels were: " + str(rowid) + ", " + str(colid))
            #    SystemExit(0)
            
            theC = theDefaultC
            
            isValid = False
            if not cFactNumpy[cFactPix[0], cFactPix[1]] == cFactND:
                theC = cFactNumpy[cFactPix[0], cFactPix[1]]
                isValid = True
            
            KLSC = KLSNumpy[rowid, colid] * theC
            KLSCFines = KLSFinesNumpy[rowid, colid] * theC
            
            statsDict[zoneVal][zoneCountStr] += 1
            statsDict[zoneVal][allKLSCSum] += KLSC
            statsDict[zoneVal][allKLSCFinesSum] += KLSCFines
            
            if(isValid):
                statsDict[zoneVal][validKLSCSum] += KLSC
                statsDict[zoneVal][validKLSCFinesSum] += KLSCFines
                statsDict[zoneVal][validDataCountStr] += 1
            
            
            
            colid += 1
        
        
        rowid += 1
        
        
    #now work out which stats to use
    #look to see if we'll use the default Cfactor derived values in lieu of valid entries
    for uid in statsDict.keys():
        
        if statsDict[uid][zoneCountStr] == 0:
            #nothing to report, add zeros
            statsDictToReturn[uid] = {KLSCvalStr: float(0), KLSCFinesvalStr: float(0)}
            continue
        
        valProp = statsDict[uid][validDataCountStr] / statsDict[uid][zoneCountStr]
        
        if valProp >= theReqValidProp:
            # use valid only stats
            statsDictToReturn[uid] = {KLSCvalStr: float(statsDict[uid][validKLSCSum] / statsDict[uid][validDataCountStr]),
                                      KLSCFinesValStr: float(statsDict[uid][validKLSCFinesSum] / statsDict[uid][validDataCountStr])}
        else:
            # Use stats that include the default where necessary
            statsDictToReturn[uid] = {KLSCvalStr: float(statsDict[uid][allKLSCSum] / statsDict[uid][zoneCountStr]),
                                      KLSCFinesValStr: float(statsDict[uid][allKLSCFinesSum] / statsDict[uid][zoneCountStr])}
        
    
    
    
    
    return statsDictToReturn