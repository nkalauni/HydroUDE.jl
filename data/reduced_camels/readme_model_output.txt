## Readme created by:  Andy Newman (NCAR/RAL/HAP)
##                     17 Apr 2014
##	     Updated:  24 Nov 2014
##           Updated:  09 Feb 2015
##
## Describes model output, summary graphics and statistics for 671 CONUS basins
## Publication in HESS, Newman et al. (2015): Development of a large-sample watershed-scale hydrometeorological
##                                            dataset for the contiguous USA: Dataset characteristics and assessment
##                                            of regional variability in hydrologic model performance

Readme:

There are four subdirectories in the model_output directory

flow_timeseries/

  Directory structure for flow timeseries data:

    Within the directory flow_timeseres/   there are three subdirectories, one for each forcing data: daymet, maurer, and nldas with 18 subdirectories
    These subdirectories are 01..18 and correspond to the USGS hydrologic unit code (HUC) level 1, 2 digit region codes
    Within each HUC subdirectory, there will be 30 files per basin gggggggg_ss_model_output.txt, gggggggg_ss_model_soil_state_output.txt, and gggggggg_ss_model_parameters.txt

    where gggggggg  is the zero padded USGS stream gauge ID and ss  is the random number starting seed for the shuffeled complex evolution algorithm (SCE).  SCE was run 10 times for each basin (see paper) and all 10 optimal parameter sets and model output are given here.

    The *model_output.txt files contain the Snow-17/SAC-SMA snow water equivalent (SWE, mm), input precipitation (PRCP (mm/day), uncorrected),
    surface water input (RAIM (mm/day), corrected using parameter SCF for frozen precipitation), mean daily air temperature (TAIR (^oC)),
    potential evapotranspiration (PET, mm/day), evapotranspiration from the SAC model (ET, mm/day), model runoff (MOD_RUN (mm/day) and observed runoff (OBS_RUN, mm/day).

    The *model_soil_state_output.txt files contain the SAC-SMA soil state variables all in units of mm:
      UZTWC -> upper zone tension water storage content
      UZFWC -> upper zone free water storage content
      LZTWC -> lower zone tension water storage content
      LZFSC -> lower zone free supplemental water storage content
      LZFPC -> lower zone free primary water storage content
      ADIMC -> additional impervious area content


    The *model_parameters.txt contain the benchmark parameter sets used to generate the model output
    following the spin-up procedure noted in the paper.
    See Snow-17 and SAC-SMA documentation for a description of the parameters.
    
    The last three parameters in the parameter files are non-standard parameters for the 
    Snow-17/SAC-SMA system.  They are the calibrated Priestly-Taylor coefficient (PT_COEF) used to 
    generate daily PET from the meteorological forcing and two parameters for a unit hydrograph
    unit_shape and unit_scale.
  

summary_plots/

  Summary_plots directory contains a panel plot for each basin for the optimal seed.  Monthly flow timeseries are in the top two rows, a two year 
  subsection of daily flow around the beginning of the validation phase (left subplot) and the flow duration curve (right subplot) are in the third
  row.  The fourth row contains a monthly flow scatter plot with daily NSE values listed (left subplot) and long-term monthly observed and modeled
  flow from the calibration and validation time periods (right subplot).

basin_stats/

  The basin_stats directory contains text files with computed statistics for each basin for the calibration and validation periods for each of the three forcing datasets and the daymet second split sample (see paper).  See the paper for a more complete description of the metrics:

  fdc_slope_bias:		The bias in the slope of the flow duration curve from the 20th to 70th percentiles
  high_flow_bias:		The bias in the top 2% flow volumes
  low_flow_bias:		The bias in the lowest 30% flow volumes
  NSE:			The Nash-Sutcliffe efficiency
  NSE_decomp_alpha:	The alpha term (flow variance) in the decomposition of NSE.  Closer to 1 is better.
  NSE_decomp_beta:	The beta term (multaplicative mean flow volume bias) in the decomposition of NSE. Closer to 1 is better.
  NSE_decomp_rho:	The correlation term (daily obs-model flow correlation) in the decomposition of NSE. Closer to 1 is better.
  MNSE:			NSE using a 31-day window long-term flow average instead of the long-term mean flow.


code/

  The code for the Priestly-Taylor implementation used here and the unit hydrograph routine
  are provided in the code.tar.gz file.