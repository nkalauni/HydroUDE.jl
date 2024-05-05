## Readme created by:  Andy Newman (NCAR/RAL/HAP)
##	               20 May 2014             
##	     Updated:  24 Nov 2014
##           Updated:  09 Feb 2015
##                     
## Describes USGS shapefiles for 671 CONUS basins 
## Publication in HESS, Newman et al. (2015): Development of a large-sample watershed-scale hydrometeorological 
##                                            dataset for the contiguous USA: Dataset characteristics and assessment 
##                                            of regional variability in hydrologic model performance

R
Readme:

The shapefiles/ directory contains the shapefiles for every basin.  There are two .zip archives for each HUC level 02 region.  One: "Region_XX_nhru_simplify_100.zip" contains the shapefile
for the outer basin boundaries (attribute: GAGEID) and for the various HRUs (attribute: hru_id)  where XX is the HUC level 02 region code.  The other .zip archive:
"Region_XX_contours.zip" contains the 100-m elevation band shapefile (attribute: GAGEID_BAND).  

Finally the file huc_02.zip contains the USGS Regional HUC 2 shapefile.
