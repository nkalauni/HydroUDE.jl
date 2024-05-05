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

  The directory hru_forcing/ contains only the daymet forcing data (see readme_basin_mean_forcing.txt) with 18 subdirectories
  These subdirectories are 01..18 and correspond to the USGS hydrologic unit code (HUC) level 1, 2 digit region codes
  Within each HUC subdirectory, there will be forcing files named:  
  

  For HRU forcing:
  gggggggg_hru_hhhhh_cida_forcing_leap.txt    

  where gggggggg  is the zero padded USGS stream gauge ID and hhhhh is the zero padded HRU number.
 
File info:

  The files have a 4 line header.  The first three lines are:

  latitude of gauge
  elevation of HRU (m)
  area of HRU (m^2)

  Then the fourth line is a text header line describing the rest of the file except the first four columns, which are: year month day hour
  Hour is set to zero since these are daily values.

  For the HRU, the elevation given is the mean elevation of the HRU.



  .list text files:

  The .list files contain one header line that gives the number of HRUs for the basin, and corresponds to the number of lines remaining in the text file.  There are two columns in the .list files, the first is a string giving the file name of each HRU forcing file and the second column is the area of that HRU (m^2).  The hypsometric data for the HRUs can be recreated using the header information in the individual forcing files or a combination of the individual forcing files and the .list files.
