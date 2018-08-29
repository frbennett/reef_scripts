"""
Robin
"""
from osgeo import gdal, gdalconst

def resample_raster_to_match_template(template_raster, input_raster, outfile, resample_method=None):
	#inputfile = #Path to input file
	#input = gdal.Open(inputfile, gdalconst.GA_ReadOnly)
	inputProj = input_raster.GetProjection()
	inputTrans = input_raster.GetGeoTransform()
	inputRef = input_raster.GetRasterBand(1)

	#referencefile = #Path to reference file
	#reference = gdal.Open(referencefile, gdalconst.GAReadOnly)
	referenceProj = template_raster.GetProjection()
	referenceTrans = template_raster.GetGeoTransform()
	bandreference = template_raster.GetRasterBand(1)
	x = template_raster.RasterXSize 
	y = template_raster.RasterYSize


	outputfile = outfile
	#driver = gdal.GetDriverByName('GTiff')
	#output = driver.Create(outputfile,x,y,1,bandreference.DataType)
	output = gdal.GetDriverByName("GTiff").Create(outputfile, x, y, 1, inputRef.DataType)
	output.GetRasterBand(1).SetNoDataValue(-9999)
	#output = gdal.GetDriverByName("MEM").Create("", x, y, 1, inputRef.DataType)
	output.SetGeoTransform(referenceTrans)
	output.SetProjection(referenceProj)
	
	theResampleMethod = gdalconst.GRA_Bilinear

	if resample_method:
		theResampleMethod = resample_method

	gdal.ReprojectImage(input_raster,output,inputProj,referenceProj,theResampleMethod)

	return output
