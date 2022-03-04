############################################################################################
# DGVMresilience.R                                                                         #
# R SCRIPT FOR ANALYSING GLOBAL VEGETATION RESILIENCE TRENDS IN A DYNAMIC GENERAL          #
#  VEGETATION MODEL (DGVM)                                                                 #
# Code author: Dr. David I. Armstrong McKay [Georesilience Analytics; Uni. Exeter]         #
# Written for: Lenton et al., 2022: "A resilience sensing system for the biosphere",       # 
#  PhilTransB, doi:10.1098/rstb.2021.0383                                                  #
# ISIMIP2b data from: https://doi.org/10.5880/PIK.2019.012 and https://data.isimip.org/    #
############################################################################################
#
### R Initialisation ###
#
## Check for & install required libraries
#
if (!require("ncdf4")) install.packages("ncdf4")
if (!require("raster")) install.packages("raster")
if (!require("foreach")) install.packages("foreach")
if (!require("doParallel")) install.packages("doParallel")
if (!require("parallel")) install.packages("parallel")
if (!require("tseries")) install.packages("tseries")
if (!require("bfast")) install.packages("bfast")
if (!require("R.utils")) install.packages("R.utils")
if (!require("Kendall")) install.packages("Kendall")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggthemes")) install.packages("ggthemes")
#
## Load libraries
#
library(ncdf4)      #to read netcdf data
library(raster)     #to manipulate raster data
library(foreach)    #for parallel processing
library(doParallel) #for parallel processing
library(parallel)   #for parallel processing
library(tseries)    #for na.remove
library(bfast)      #for bfast function (STL decomposition)
library(R.utils)    #for isZero
library(Kendall)    #for Kendall tau calculation
library(bigmemory)  #for time-series analysis output
library(ggplot2)    #for plotting
library(ggthemes)   #for plotting
#
### Input data selection ###
#
## Data pre-processing in linux command line (*not* in R, requires cdo & nco packages)
# (This script analyses LPJmL output for 2000-2019 (driven by GFDL-ESM2M historical emissions & land-use forcing 2000-2005 
# & RCP8.5 emissions & 2005 land use 2006-2019 w/ EWEMBI bias correction) from ISIMIP2b as an example - see links above.
# After downloading appropriate DGVM model output data, adapt & use the following code to select data for the required date-range:)
#
#cdo -selyear,2000/2005 lpjml_gfdl-esm2m_ewembi_historical_histsoc_co2_npp_global_monthly_1970_2005.nc lpjml_gfdl-esm2m_ewembi_historical_histsoc_co2_npp_global_monthly_2000_2005.nc #select years
#cdo -selyear,2006/2019 lpjml_gfdl-esm2m_ewembi_rcp85_2005soc_co2_npp_global_monthly_2006_2099.nc lpjml_gfdl-esm2m_ewembi_rcp85_2005soc_co2_npp_global_monthly_2006_2019.nc #select years
#ncrcat lpjml_gfdl-esm2m_ewembi_historical_histsoc_co2_npp_global_monthly_2000_2005.nc lpjml_gfdl-esm2m_ewembi_rcp85_2005soc_co2_npp_global_monthly_2006_2019.nc lpjml_gfdl-esm2m_ewembi_historical_rcp85_co2_npp_global_monthly_2000_2019.nc #stitch netcdf files together
#ncks -4 -L 1 lpjml_gfdl-esm2m_ewembi_historical_rcp85_co2_npp_global_monthly_2000_2019.nc lpjml_gfdl-esm2m_ewembi_historical_rcp85_co2_npp_global_monthly_2000_2019.nc #quick deflate to save space
#ncdump -h lpjml_gfdl-esm2m_ewembi_historical_rcp85_co2_npp_global_monthly_2000_2019.nc #inspect file header
#
## Select netcdf filename & variable
#
setwd("/home/DGVMresilience/") #set top directory for analysis - replace with your own
filepath <- ("Input_Data/ISIMIP_2b/LPJmL/GFDL-ESM2M/") #filepath to input data folder
filename <- "lpjml_gfdl-esm2m_ewembi_historical_rcp85_co2_npp_global_monthly_2000_2019" # name of pre-processed model output file
filenamenc <- paste(filepath,filename,".nc", sep="")
dname <- "npp"  #use "npp" for NPP (adjust if otherwise - could be any model output)
#
### Prepare data for analysis ###
#
## Create stacked raster for NPP, convert to dataframe, & and trim out NAs (to slim filesize)
#
ncfile <- stack(filenamenc,varname=dname) #make a raster stack for selected variable
ncfile_df <- as.data.frame(ncfile, xy = TRUE) #convert to dataframe (df), to enable time-series analysis
ncfile_trim_df<-ncfile_df[rowSums(is.na(ncfile_df)) != (ncol(ncfile_df)-2),] #remove NAs (i.e. ocean)
rm(ncfile_df) #remove original untrimmed dataframe to save space
#
## Make raster of mean across time-range, and make a mask of low NPP areas (to avoid bias)
#
ncfile_mean <- mean(ncfile) #create mean of rasterstack
plot(ncfile_mean) #plot to inspect
ncfile_mean_df <- as.data.frame(ncfile_mean, xy = TRUE) #convert to dataframe 
ncfile_mean_trim_df <- ncfile_mean_df[rowSums(is.na(ncfile_mean_df)) != (ncol(ncfile_mean_df)-2),] #remove NAs (i.e. ocean)
rm(ncfile_mean_df) #remove untrimmed df to save space
ncfile_mean_trim_masked_df <- ncfile_mean_trim_df #initialise data to mask
ncfile_mean_trim_masked_df[ncfile_mean_trim_masked_df[,3] < 0.1e-08,3] <- NA #replace all NPP values below threshold with NAs
mask <- ncfile_mean_trim_masked_df[is.na(ncfile_mean_trim_masked_df[,3]),] #create layer showing where data is masked
#
### Time-series analysis on a rolling window ###
#
## Prepare for analysis loop 
# (here set up for autocorrelation function (acf), but standard deviation (sd) also an option)
#
n <- nrow(ncfile_trim_df) #length (i.e. no. pixels) of input data
ncl <- ncol(ncfile_trim_df) #number of columns (i.e. time-slices plus first 2 co-ordinate columns) in input data
N <- ceiling(ncl/2) #window length for 50% rolling window (50% is standard, but can be adjusted)
Q <- ncl-2-(N-1) #last column to start rolling window on
bfast_res_ni_acf <- big.matrix(nrow=n,ncol=Q,init=NA) #empty matrix for bfast acf results to write to
bfast_roll_acf <- rep(as.numeric(NA), Q) #empty vector to write acf results to
#bfast_res_ni_sd <- matrix(NA,nrow=n,ncol=Q) #uncomment these (& comment 2 above) to look at SD instead
#bfast_roll_sd <- rep(as.numeric(NA), Q)
#
## Time-series analysis to detrend data with STL & calculate ACF(/SD) - loop in series 
# (much slower - takes days at 0.5x0.5degree resolution - use if no multiple cores available)
#
for (i in 1:n) { #for each input row (i.e. each cell)
  for (q in 1:Q) { #repeat a 50% rolling window starting until ~50% through row
    x <- ts(as.numeric(ncfile_trim_df[i,(q+2):(q+2+N-1)]),frequency=12) #convert 50% segment of input to ts object for (+ sampling frequency, LPJmL=12 for 12 months per year)
    if (all(is.na(x))==FALSE) { #double-check no NA cells remaining (or bfast breaks)
      X_na <- na.remove(x) #remove NAs (an issue for acf/bfast - shouldn't be common for model output though)
      if ((all(isZero(X_na))==FALSE) & (length(X_na)>45)) { #check not all zeroes and not too many entries removed in previous step
        bfast_output <- bfast(X_na,h=1,season=c("harmonic"),max.iter=1) #bfast to detrend input (linear & seasonal via STL)
        bfast_acf <- acf(bfast_output$output[[1]]$Nt,1,plot=FALSE)$acf[2] #post-detrending ACF for once cell at a time
        #bfast_sd  <- sd(bfast_output$output[[1]]$Nt)
      } else { #if all zeroes or too few entries in cell then result is NA
        bfast_acf <- NA
        #bfast_sd <- NA
      }
    } else {bfast_acf <- NA} #if all NAs in cell then result is NA
    bfast_roll_acf[q] <- bfast_acf  #add bfast_acf to matrix  
    #bfast_roll_sd[q] <- bfast_sd
  }
  bfast_res_ni_acf[i,] <- bfast_roll_acf #add to results matrix
  #bfast_res_ni_sd[i,] <- bfast_roll_sd
  message(i) #report progress
}
save.image("./ISIMIP2b_LPJmL_GFDLESM2M_20002019_roll.RData") #save progress
#
## Time-series analysis to detrend data with STL & calculate ACF(/SD) - loop in parallel 
# (much faster, use if multiple cores available - N.B. still some bugs here, doesn't always work)
#
no_cores <- 20 #number of cores to use - use detectCores() to find out what's available
registerDoParallel(makeCluster(no_cores)) #activate cores
#
bfast_res_ni_acf <- foreach(i=1:n,.combine='rbind',.packages=c("bfast","tseries","R.utils","bigmemory")) %dopar% {
  for (q in 1:Q) { #repeat a 50% rolling window starting until ~50% through row
    x <- ts(as.numeric(ncfile_trim_df[i,(q+2):(q+2+N-1)]),frequency=12) #convert 50% segment of input to ts object for (+ sampling frequency, LPJmL=12 for 12 months per year)
    if (all(is.na(x))==FALSE) { #double-check no NA cells remaining (or bfast breaks)
      X_na <- na.remove(x) #remove NAs (an issue for acf/bfast - shouldn't be common for model output though)
      if ((all(isZero(X_na))==FALSE) & (length(X_na)>45)) { #check not all zeroes and not too many entries removed in previous step
        bfast_output <- bfast(X_na,h=1,season=c("harmonic"),max.iter=1) #bfast to detrend input (linear & seasonal via STL)
        bfast_acf <- acf(bfast_output$output[[1]]$Nt,1,plot=FALSE)$acf[2] #post-detrending ACF for once cell at a time
        #bfast_sd  <- sd(bfast_output$output[[1]]$Nt)
      } else { #if all zeroes or too few entries in cell then result is NA
        bfast_acf <- NA
        #bfast_sd <- NA
      }
    } else {bfast_acf <- NA} #if all NAs in cell then result is NA
    bfast_roll_acf[q] <- bfast_acf  #add bfast_acf to matrix  
    #bfast_roll_sd[q] <- bfast_sd
  }
  bfast_res_ni_acf[i,] <- bfast_roll_acf #add to results matrix
  #bfast_res_ni_sd[i,] <- bfast_roll_sd
  return(bfast_res_ni_acf) #(bfast_res_ni_sd) #send raster back
  #return(bfast_res_ni_sd)
}
stopImplicitCluster() #turns off parallel core (must do!)
save.image("./ISIMIP2b_LPJmL_GFDLESM2M_20002019_roll.RData") #save progress
#
### Create mapped results ###
#
## Create dataframe of NPP ACF(/SD) results with co-ordinates
#
GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf <- array(NA, dim=c(n,2)) #create empty matrix
GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf[,1] <- ncfile_trim_df[1:n,1] #insert X co-ordinate column
GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf[,2] <- ncfile_trim_df[1:n,2] #insert y co-ordinate column
GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf <- as.data.frame(cbind(GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf,bfast_res_ni_acf[1:n,1:Q])) #insert NPP ACF results
#
# Create spatial pixels dataframe (spdf) & raster of NPP ACF(/SD) results
#
plot_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf <- GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf #initialise variable
names(plot_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf) <- c('x','y','acf') #label columns
coordinates(plot_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf) = ~x+y #assign co-ordinates for spdf gridding
gridded(plot_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf) = TRUE #turn into an spdf
plot_GFDL_rcp8p5_20002019_npp_roll_bfast_outraster_ni_acf <- brick(plot_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf) #turn spdf in to a raster
#
## Create mean of NPP ACF(/SD) results & convert to spdf/df
#
plot_GFDL_rcp8p5_20002019_npp_roll_bfast_outraster_ni_acfmean <- mean(plot_GFDL_rcp8p5_20002019_npp_roll_bfast_outraster_ni_acf) #calculate mean of NPP ACF(/SD) 
plot_GFDL_rcp8p5_20002019_npp_roll_bfast_outraster_ni_acfmeanspdf <- as(plot_GFDL_rcp8p5_20002019_npp_roll_bfast_outraster_ni_acfmean,"SpatialPixelsDataFrame") #convert to spdf
plot_GFDL_rcp8p5_20002019_npp_roll_bfast_outraster_ni_acfmeandf <- as.data.frame(plot_GFDL_rcp8p5_20002019_npp_roll_bfast_outraster_ni_acfmeanspdf, xy = TRUE) #convert to df
colnames(plot_GFDL_rcp8p5_20002019_npp_roll_bfast_outraster_ni_acfmeandf) <- c("acf","x", "y") #name columns of NPP ACF(/SD) results df
#
### Kendall Tau trend calculation & plot ###
#
## Calculate Mann Kendall Tau coefficient for each pixel on a loop
#
kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50 <- rep(NA, n) #create empty vector
#kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_sd_RW50 <- rep(NA, n)
for (l in 1:n) { #for every row (i.e. pixel)
  kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50[l] <- MannKendall(GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf[l,3:Q+2])$tau #calculate kendall tau
  #kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_sd_RW50[l] <- MannKendall(GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_sd[l,3:Q+2])$tau
  message(l) #report progress
}
#
## Convert Kendall Tau results to spdf/raster
#
plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50 <- cbind(GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf[,1:2],kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50) #add co-ordinates
names(plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50) <- c('x','y','kendall_tau') #rename columns
plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50_outnidf <- plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50 #save a copy as a df for analysis below
coordinates(plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50) = ~x+y #assign co-ordinates for spdf gridding
gridded(plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50) = TRUE #turn into an spdf
plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50_outniraster <- raster(plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50) #turn in to raster for plotting
#
### Calculate mean Kendall Tau value across unmasked area (& where NPP ACF > 0)
#
plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50_outnidf <- plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50_outnidf[order(-plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50_outnidf$y,plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50_outnidf$x),] #make sure df order matches ncfile_df (otherwise misaligned)
plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50_outnidf_masked <- plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50_outnidf[!is.na(ncfile_mean_trim_masked_df[,3]),] #mask NPP ACF Kendall Tau layer with low NPP mask
KendallTau_unmaskedarea_mean <- mean(plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50_outnidf_masked[,3],na.rm=TRUE) #calculate mean of Kendall Tau in unmasked area
KendallTau_unmaskedpositiveNPPACFarea_mean <- mean(plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50_outnidf_masked[plot_GFDL_rcp8p5_20002019_npp_roll_bfast_outraster_ni_acfmeandf[,1]>0,3],na.rm=TRUE) #calculate mean of Kendall Tau in unmasked area where NPP ACF > 0
save.image("./ISIMIP2b_LPJmL_GFDLESM2M_20002019_roll.RData") #save progress

#
### Final plotting for results for both NPP ACF(/SD) & associated Kendall Tau
#
## Prep datasets for plotting
#
plot_GFDL_rcp8p5_20002019_npp_roll_bfast_outraster_ni_acfmean[is.na(plot_GFDL_rcp8p5_20002019_npp_roll_bfast_outraster_ni_acfmean)] <- -2 #set oceans to out-of-range value for plotting
plot_GFDL_rcp8p5_20002019_npp_roll_bfast_outraster_ni_acfmeanspdf <- as(plot_GFDL_rcp8p5_20002019_npp_roll_bfast_outraster_ni_acfmean,"SpatialPixelsDataFrame") #convert to spdf/df for ggplot
plot_GFDL_rcp8p5_20002019_npp_roll_bfast_outraster_ni_acfmeandf <- as.data.frame(plot_GFDL_rcp8p5_20002019_npp_roll_bfast_outraster_ni_acfmeanspdf, xy = TRUE)
colnames(plot_GFDL_rcp8p5_20002019_npp_roll_bfast_outraster_ni_acfmeandf) <- c("NPP ACF","x","y") #rename variables for plot legend
#
plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50_outniraster[is.na(plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50_outniraster)] <- -2 #set oceans to out-of-range value for plotting
plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50_outnispdf <- as(plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50_outniraster,"SpatialPixelsDataFrame") #convert to spdf/df for ggplot
plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50_outnidf <- as.data.frame(plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50_outnispdf, xy = TRUE)
colnames(plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50_outnidf) <- c("NPP ACF Kendall Tau","x", "y") #rename variables for plot legend
#
## Plot NPP ACF layer
#
ggplot() +
  geom_tile(data=plot_GFDL_rcp8p5_20002019_npp_roll_bfast_outraster_ni_acfmeandf,aes(x=x,y=y,fill=`NPP ACF`)) + 
  layer(geom="tile",stat="identity",data=mask,aes(x=x,y=y),position="identity",inherit.aes=FALSE,show.legend=FALSE) +
  coord_quickmap() + #uncomment this (& comment next line) for a quick plot
  #coord_map('rectangular',lat0=30) + #uncomment this (& comment previous line) for a full (slow) plot
  scale_fill_gradient2(low="#ffa201",mid="white",high="#42016f",midpoint=0,space="rgb",
                       na.value="#d5e8ec",guide="colourbar",limits=c(-1,1)) +
  theme_map() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(legend.title=element_text(colour="black",size=15,face="bold")) 
  theme(legend.text=element_text(colour="black",size=15)) + 
  theme(legend.position=c(0,0)) + 
  theme(legend.key.width=unit(3,"cm"))
#
## Plot NPP ACF Kendall Tau layer
#
ggplot() +
  geom_tile(data=plot_kendtau_GFDL_rcp8p5_20002019_npp_roll_bfast_out_ni_acf_RW50_outnidf,aes(x=x,y=y,fill=`NPP ACF Kendall Tau`)) + 
  layer(geom="tile",stat="identity",data=mask,aes(x=x,y=y),position="identity",inherit.aes=FALSE,show.legend=FALSE) +
  coord_quickmap() + #uncomment this (& comment next line) for a quick plot
  #coord_map('rectangular',lat0=30) + #uncomment this (& comment previous line) for a full (slow) plot
  scale_fill_gradient2(low="#2b83c3",mid="white",high="#ca0020",midpoint=0,space="rgb",
                       na.value="#d5e8ec",guide="colourbar",limits=c(-1,1)) + 
  theme_map() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(legend.title=element_text(colour="black",size=15,face="bold")) +
  theme(legend.text=element_text(colour="black",size=15)) +
  theme(legend.position=c(0,0)) +
  theme(legend.key.width=unit(3,"cm"))
#
### CODE END ###