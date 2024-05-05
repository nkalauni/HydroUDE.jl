## Readme created by:  Andy Newman (NCAR/RAL/HAP)
##               17 Dec 2013  Internal Version 0
##     Updated:  20 May 2014  Updated Internal Version 0.1
##     Updated:  24 Nov 2014  First Public Release Version 1.0
##     Updated:  09 Feb 2015  Version 1.1
##     Updated:  11 Mar 2016  Version 1.2 
##                     
## Publication in HESS, Newman et al. (2015): Development of a large-sample watershed-scale hydrometeorological 
##                                            dataset for the contiguous USA: Dataset characteristics and assessment 
##                                            of regional variability in hydrologic model performance

Readme:

There are many sub-directories that contain forcing data, model output, observed streamflow, basin shapefiles, etc.  There are several other readmes that
describe the directory structure and file structure in more detail.  All the data (except shapefiles) are in whitespace delimited ASCII files.  Basin shapefiles are in standard shapefile format and summary graphics are in Portable Network Graphics (PNG) format.

The basin_metadata directory that contains basic information about every basin included in calibration.
The file gauge_information.txt contains the Regional Hydroloic Unit Code from the USGS, the USGS stream gauge ID, gauge name, latitude, longitude and basin size for all 671 basins.  
This directory also contains basin hydrometeorological and physical characteristic information for each basin.  The hydrometeorological data are based on the following
time periods:

Daymet:  1980-2014
NLDAS:   1980-2014
Maurer:  1980-2008

Please also carefully read the dataset_summary.txt file for a great summary of
specific important characteristics of this dataset.  Thanks to Nans Addor of
NCAR for putting this together and Drew Guswa of Smith College for further input.

If you have any questions please contact me via email at: anewman@ucar.edu



Version 1.0:

First public release of dataset


Version 1.1:
Updates from Version 1.0:

Found optimal parameter files in model_output directory had erroreous values in V1.0, reprocessed to correct.
Also changed all the gauge naming conventions to 8-digits zero padded from 9-digits.  This was done because while
the USGS does have streamflow gages/basins with 9-digit identifiers, none of those are included in this dataset.

Version 1.2:
Updates from Version 1.1:

Updated Daymet and NLDAS basin mean forcing through 12/31/2014 for all basins.
Updated Daymet and NLDAS model output through 12/31/2014 (or end of observed streamflow record) for all basins.

Added dataset_summary.txt file (developed by Nans Addor of NCAR/RAL) with a summary and short descriptions of important
dataset characteristics.  This includes identification of basin size differences between the USGS geospatial
fabric and the USGS GAGES-II watershed delineations.  Thanks to Carolina Massmann at Bristol University, UK for 
originally pointing this out.  Also see paper by Bock et al. (2015), "Parameter regionalization of a monthly
water balance model for the conterminous United States", in HESS-D for a more complete description of
the GF and differences with GAGES-II.

We were able to manually fix 4 basins in version 1.1 from the USGS GF.  All forcing, model output, etc
with these four basins are now corrected in this version.  Those four basins are:

02108000
05120500
07067000
09492400

Fixed issue with missing data in observed streamflow.  Now represented with -999 flow value and M for the QC flag instead of
blanks.

Fixed issue with a missing day of model output between the calibration and validation periods.  Note that the
calibraition and validation periods are two distinct model runs with unique spin-up periods.  Therefore,
the flow on the last day of the calibration period will not be continuous into the first day of the validation period.

Fixed issue with small negative PET for some basins in winter.  P-T PET was not limited to values >= 0.  There are very
minor changes to the model output (PET, SAC ET and runoff) as compared to version 1.1 for these basins.
