#########################################################################
#
# PRRFECT INPUT FILE v1
#
# Prequisites:
#       - R (free download from http://www.r-project.org/)
#       - R packages: "randomForest", "xtable", "caret", "prospectr", "signal", "wrt", "KernSmooth", "Peaks", and "hyperSpec"
#               (-- to install packages from the R command prompt, type:
#                       - install.packages("randomForest,  dependencies = TRUE)
#                       - install.packages("xtable"),  dependencies = TRUE )
#			- install.packages("caret"),  dependencies = TRUE)
#			- install.packages("prospectr"),  dependencies = TRUE)
#			- etc.
#
# Before using this script, search for "#User-defined variables" below and set
# the following variables to the appropriate values:
#
# k 				       #What value of k should be used in k-fold cross-validation?  (Try k=5 or k=10)
# se_and_ci                            #Should standard errors and confidence intervals be included in text result files? 
# se_and_ci2                           #Should standard errors and confidence intervals be included in tabular result files? 
# infile_tr			       #The name of the file containing the training dataset
# infile_te			       #The name of the file containing the testing dataset
#
# Notes:
#
#########################################################################

#######################################################
#######################################################
#       INPUT USER-DEFINED VARIABLES        #
#       ###################################
#       


infile_tr<-"input_training_dataset.ssv" #Name of file containing the training dataset
infile_te<-"input_testing_dataset.ssv"  #Name of file containing the testing dataset
#######################################################
#######################################################
#       CLASSIFIER AND OUTPUT USER-DEFINED VARIABLES
#       ###################################
#       

k<-5                                 	#What value of k should be used in k-fold cross-validation?
se_and_ci<-FALSE                      	#Should standard errors and confidence intervals be included in text result files?
se_and_ci2<-FALSE                     	#Should standard errors and confidence intervals be included in tabular result files?


bin_factor <- 2       #What factor should wavenumber resolution be reduced by? 1 for none.  
smooth_choice <- 1         #Smoothing selection.  Choices are: 0 for none.  1 for SG filter.  2 for wavelet denoising.  3 for local polynomial fit with Gaussian weighting.
smooth_par <- 2            #Parameter for smoothing selection. If option 1 - filter order.  If option 2 - length of filter.  If option 3 - bandwith of Gaussian.
norm_choice <- 2           #Normalisation choice. 0 for none.   1 for min/max to 0-1.  2 for vector.  3 for normalisation to amide I band.
bg_choice <- 1             #Baseline correction choice. 0 for none.  1 for first derivative.  2 for second derivative. 3 for rubberband.  4 for polynomial.
bg_par <- 0                #ONLY for options 3 and 4 above. If option 3 - noise cut-off level.  If option 4 - polynomial degree.
RBp <- 0                 #ONLY for Rubberband - noise cutoff parameter - match this to approximate noise floor of spectra

spectra <- TRUE            # Plot all spectra graph?
average_spectra <- TRUE    # Plot average spectra graph?
imp_and_all_spectra <-TRUE # Plot importance and spectra combo graph?
avg_and_gini <- TRUE


#Select columns to classify. This will build a classification model based upon the classes in this column(s). Can set more than 1 column at a time.
yvariables<-c("country")

#Options
min_wavenumber<-800                  #Must be multiple of 100!  min_wavenumber <= selected frequencies <= max_wavenumber
max_wavenumber<-1800                 	#Must be multiple of 100!  min_wavenumber <= selected frequencies <= max_wavenumber 
patients<-"sample_number"         #Name of column containing patient IDs or sample numbers etc.

###NOTE ABOUT _par SETTINGS - These will only affect the current selected _choice.  If option 0 is selected for _choice, _par is ignored.
###   

#example1: smooth_choice of 1. smooth_par of 3 will give a SGfilter of order 3.
#example2: smooth_choice of 3. smooth_par of 24 will give a local polynomial filter with Gaussian bandwidth of 24.
#example3: bg_choice of 3. bg_par of 4 will give a rubberband function with quadratic equation of factor 4

