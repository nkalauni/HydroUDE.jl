## Readme created by:  Andy Newman (NCAR/RAL/HAP)
##                     30 Dec 2013 
##           Updated:  17 Apr 2014
##	     Updated:  24 Nov 2014
##           Updated:  09 Feb 2015
##
## Describes elevation band forcing dataset for 671 CONUS basins
## Publication in HESS, Newman et al. (2015): Development of a large-sample watershed-scale hydrometeorological 
##                                            dataset for the contiguous USA: Dataset characteristics and assessment 
##                                            of regional variability in hydrologic model performance 


Readme:

Directory structure forcing data:

  The directory elev_bands_forcing/ contains only the daymet forcing data (see readme_basin_mean_forcing.txt) with 18 subdirectories
  These subdirectories are 01..18 and correspond to the USGS hydrologic unit code (HUC) level 1, 2 digit region codes
  Within each HUC subdirectory, there will be forcing files named:  
  

  For the elevation band forcing:
  gggggggg_elev_band_bbb_cida_forcing_leap.txt

  where gggggggg is the zero padded USGS stream gauge ID and bbb is the zero padded elevation band number.
  
  The bands are numbered according to elevation starting from sea level such that band 000 is the 0 - 100 m elevation band, band 001 is the 100-200 m elevation band and so forth.


File info:

  The files have a 4 line header.  The first three lines are:

  latitude of gauge
  base elevation of band (m)
  area of band (m^2)

  Then the fourth line is a text header line describing the rest of the file except the first four columns, which are: year month day hour
  Hour is set to zero since these are daily values.

  For the elevation band files, the elevation given is the LOWEST elevation of the band, the mid-point is 50 m higher.


  *.list text files:

  The .list files contain one header line that gives the number of elevation bands for the basin, and corresponds to the number of lines remaining in the text file.  There are two columns in the .list files, the first is a string giving the file name of each elevation band forcing file and the second column is the area of that elevation band (m^2).  The elevation band hypsometric data can be recreated from the .list file using the elevation band number where:  band_elev(m) = band_num*100 + 50 and the area given in the second column.
