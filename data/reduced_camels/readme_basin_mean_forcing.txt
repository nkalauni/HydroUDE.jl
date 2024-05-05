## Readme created by:  Andy Newman (NCAR/RAL/HAP)
##                     17 Dec 2013 
##           Updated:  17 Apr 2014
##	     Updated:  20 May 2014
##	     Updated:  24 Nov 2014
##           Updated:  09 Feb 2015
##                     
## Describes basin mean forcing (lump) dataset for 671 CONUS basins 
## Publication in HESS, Newman et al. (2015): Development of a large-sample watershed-scale hydrometeorological 
##                                            dataset for the contiguous USA: Dataset characteristics and assessment 
##                                            of regional variability in hydrologic model performance

Readme:

This readme describes the directory and file structure for the basin mean forcing data.


Directory structure for basin mean forcing data:

  The directory: basin_mean_forcing/ contains three directories: daymet (Thornton et al. 2012, http://daymet.ornl.gov/), maurer (Maruer, et al. 2002) and nldas (Xia et al. 2012).  These are three distinct forcing datasets generated in different manners with different spatial resolutions, all at a daily time step.
  Each forcing data directory contains 18 subdirectories.
  These subdirectories are 01..18 and correspond to the USGS hydrologic unit code (HUC) level 1, 2 digit region codes
  Within each HUC subdirectory, there will be forcing files named:  gggggggg_lump_cida_forcing_leap.txt for daymet;  gggggggg_lump_maurer_forcing_leap.txt for maruer;   gggggggg_lump_nldas_forcing_leap.txt for nldas
  where gggggggg  is the 8-digit zero padded USGS stream gauge ID.  
 
File info:

  The files have a 4 line header.  The first three lines are:

  latitude of gauge
  elevation of gauge (m)
  area of basin (m^2)

  Then the fourth line is a text header line describing the rest of the file except the first four columns, which are: year month day hour
  Hour is set to zero since these are daily values.

