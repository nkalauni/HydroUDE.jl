## Readme created by:  Andy Newman (NCAR/RAL/HAP)
##                     17 Dec 2013 
##           Updated:  17 Apr 2014
##	     Updated:  20 May 2014
##	     Updated:  24 Nov 2014
##           Updated:  09 Feb 2015
##                     
## Describes USGS streamflow dataset for 671 CONUS basins 
## Publication in HESS, Newman et al. (2015): Development of a large-sample watershed-scale hydrometeorological 
##                                            dataset for the contiguous USA: Dataset characteristics and assessment 
##                                            of regional variability in hydrologic model performance

Readme:

This readme describes the directory and file structure for the USGS observed streamflow data.


Directory structure for streamflow data:

  The directory usgs_streamflow/  contains 18 subdirectories
  These subdirectories are 01..18 and correspond to the USGS hydrologic unit code (HUC) level 1, 2 digit region codes
  They contain streamflow files named:  gggggggg_streamflow_qc.txt
  where gggggggg  is the zero padded USGS stream gauge ID.  

File info:

  The files have no header with five columns of data:  GAGEID Year Month Day Streamflow(cubic feet per second) QC_flag

  Streamflow data that are missing are given the streamflow value -999.0

  The files are daily mean streamflow, so the hour is always just set to 0.
  The QC flag is a string and can have to values: "A","A:e", "M"
	A   -> streaflow value is certified by USGS as the actual daily mean flow
	A:e -> streamflow value is certified by the USGS as the actual ESTIMATED daily mean flow
        M   -> streamflow value is missing from USGS record
