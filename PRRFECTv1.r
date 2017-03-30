###Version 0.0.1  Very similar to our original code.  I have disabled the d2ss derivative, and implemented >
#< prospectr savitzkyGolay command.  Not in function as yet.
###Version 0.0.2  Attempted to fix rubberband and tested various functions

###Version 0.0.3 All functions now in this program 

### Note 1 - Sticking with x1N system we already had.  Problem with whilecount.

### Note 2 - Ran into problem with binning.  It gets rid of datapoints, which doesn't work nicely with out code.

###Version 0.0.4 Added functionality for plotting spectra only, and mean spectra from each class.
###Also, non-cancer class 3 is removed from the dataset when yvar is B_Met.v_.... this allows for the goal of the current project,
###where we want to do a 2 level classification.
###Things for next version
### - add choice of which graphs are wanted, all 4 are not required in all situations
### - generalise the selective removal of classes from columns.

### Version 0.0.5 Got most things working. Fixed Binning problems, all features work. More testing required. Should move
### function switches to top in next.  For readability.

### Version 0.0.6 All functions moved to top, in right order. Tidying needed.
### 
### Version 0.0.7 To do - add chopco2 function, and get graphs to plot correctly. <- DONE, kind of.

### Version 0.0.8 More things working.  Did some testing.

### Version 0.1.0 Tidy up the code. Make input file format <- DONE.  Also, made diagnostic mode for quick graphs.  Now we also have avg_and_gini for even faster information.  Decided to remove plotting of meandecrease_accuracy.

### Version 0.1.1 Lightened the resulting RData files considerably, removed all datasets at the end, saved patients and seed.

### Version 0.1.1a Archie version where input files are deleted immediately after reading in.  Solves problem of failed runs leaving large files in place
source("user_defined_input.R")
### Version 0.1.2  Rubberband now has sensible parameters, noise level is bg.p/100

### Version 0.1.3  Wavelet now has daub coeff of sm.p*2 

### Version 0.1.4  Derivatives now more sensibly written into function

### Fruit 0.1.4a  Rubberband has bg.p/5000 based on overall signal level of spectra

### Olive oil 0.1.4a Rubberband set to bg.p/200

### Coffee 0.1.4a Rubberband set to bg.p*3
#######################################################
#
# Functions
#
#######################################################

#Libraries
#library(colorout)    # For BRS to see things easier
#setOutputColors256() # Same
library(lattice)     # randomForest auto-attaches these.  I was getting a masking problem, so have just explicitly loaded these.
library(ggplot2)     # ^
library(randomForest)	#randomForest() used to train Random Forests
library(xtable)		#xtable() used to write .tex and .html files
library(caret)		#createDataPartition() used to create training/testing split
library(sfsmisc)	#D1ss() and D2ss() used for numerical first and second derivatives
library(prospectr) #implements new derivatives with SG smoothing
library(signal)    # for sgolayfilt
library(rwt)      #for wavelet denoising
library(KernSmooth)  #for local polynomial fit
library(Peaks)     # 
library.dynam('Peaks', 'Peaks', lib.loc=NULL)
library(hyperSpec)   #for rubberband

#Function: modify D2ss() so that it returns second derivative as a numeric vector rather than a list
D2ssy<-function(x=NULL, y=NULL){
  return(as.numeric(D2ss(x=x,y=y)$y))
}

#POSSIBLE future updates, not active at the moment.  Please leave set to defaults.
preprocessing<-"none"               	
normalize<-FALSE                      	
minzero<-FALSE                       	
se_and_ci<-FALSE                      	
se_and_ci2<-FALSE     
PCNR <- 0
chop_co2 <- FALSE


#Function: Print results to screen and to file
print_results<-function(file=NULL, mystats=NULL, ftr=NULL, fte=NULL, myclass=NULL, yvar=NULL, what=NULL, se_and_ci=TRUE){

  sink(file=file, split=TRUE)
  cat("\n\n", what, " statistics for: y-variable=",yvar," and class=", myclass,"\n")
  cat("No. of items in training dataset:", dim(ftr)[1],"\n")
  cat("No. of items in testing dataset:", dim(fte)[1],"\n")
  cat("No. of items labelled ", myclass, "in training dataset:", sum(ftr[,y]==myclass),"\n")
  cat("No. of items labelled ", myclass, "in testing dataset:", sum(fte[,y]==myclass),"\n")
  cat("\nStandard errors are given in parentheses, 95% Confidence intervals are given in square brackets\n")
  cat("Prediction accuracy (te):", mystats$Q, "\n")
  cat("Matthews correlation coefficient (te):", mystats$C, "\n")

  if(se_and_ci){
    cat("Sensitivity - (+ accuracy) (te):", mystats$sensitivity, " (",mystats$se_sensitivity,")", " [",paste(mystats$ci_sensitivity$ci95,collapse=",") ,"]\n", sep="")
    cat("Specificity - (- accuracy) (te):", mystats$specificity, " (",mystats$se_specificity,")", " [",paste(mystats$ci_specificity$ci95,collapse=",") ,"]\n", sep="")
    cat("Precision - Positive (te):", mystats$precision_p, " (",mystats$se_precision_p,")", " [",paste(mystats$ci_precision_p$ci95,collapse=",") ,"]\n", sep="")
    cat("Precision - Negative (te):", mystats$precision_n, " (",mystats$se_precision_n,")", " [",paste(mystats$ci_precision_n$ci95,collapse=",") ,"]\n", sep="")
    cat("\n\n")
  } else {
    cat("Sensitivity - (+ accuracy) (te):", mystats$sensitivity, "\n")
    cat("Specificity - (- accuracy) (te):", mystats$specificity, "\n")
    cat("Precision - Positive (te):", mystats$precision_p, "\n")
    cat("Precision - Negative (te):", mystats$precision_n, "\n")
  }
  sink(NULL)

}

#Function: Split data into train and test sets so that proportions of different classes are correct and so no patient is represented in both training and test sets.
#f is the dataset
#y is the name of the column containing the classes
#patients is the name of the column containing the patient IDs
get_tr_te_split<-function(f=NULL, y=NULL, patients=patients){

  #Get two-column matrix containing classification in column 1 (i.e. cancer .v. non-cancer) and patient number in column 2
  patients_classifications<-matrix(unlist(strsplit(unique(paste(f[,y],f[,patients], sep=";")), ";")), ,2, byrow=TRUE)

  #Select 10% of patients in each class
  tr_patients<-vector()
  te_patients<-vector()
  for(i in unique(f[,y])){

    #Get row numbers of all items in class
    all_items_in_class<-which(patients_classifications[,1]==i)
 
    #Split items into training (90%) and testing (10%) sets, while checking that each class has at least 10 elements.
    tmp<-1:length(all_items_in_class)
    if(length(tmp)>=10){
      tmp2<-sample(tmp,round(length(tmp)/10), replace=FALSE)
    } else {
      cat("\nWARNING: Class ", i, " in dataset ", yvar, " has only got ", length(tmp), " items in it!\n\n")
      tmp2<-1
    }
    tr_patients<-c(tr_patients,all_items_in_class[-tmp2])
    te_patients<-c(te_patients,all_items_in_class[tmp2])
 
    #Print
    #cat("Train:", length(all_items_in_class[-tmp2]),"\n")
    #cat("Test:", length(all_items_in_class[tmp2]),"\n")
    #cat("--\n")
  }
 
  #Convert patient numbers to spectra numbers  (there are multiple spectra for each patient)
  tr<-vector()
  for(i in tr_patients){
    patient<-patients_classifications[i,2]
    tr<-c(tr,which(f[,patients]==patient))
  }
  te<-vector()
  for(i in te_patients){
    patient_no<-patients_classifications[i,2]
    te<-c(te,which(f[,patients]==patient_no))
  }

  #Check
  all<-c(tr,te)
  if(length(all)!=dim(f)[1]){
    stop("ERROR! Training and test sets have wrong number of molecules in them")
  }
  if(length(all)!=length(unique(all))){
    stop("ERROR! Training and test sets allocated incorrectly")
  }
  tmp<-1:dim(f)[1]
  tmp2<-sort(all, decreasing=FALSE)
  if(sum((1:dim(f)[1]) == sort(all, decreasing=FALSE))!=length(tmp)){
    stop("ERROR! Wrong number of training and test set molecules")
  }
 
  return(list(tr=tr,te=te))
}




#Function: compute classification statistics
cstats<-function(tp=NULL, tn=NULL, fp=NULL,fn=NULL){

  tp<-as.numeric(tp); tn<-as.numeric(tn); fp<-as.numeric(fp); fn<-as.numeric(fn);   #HACK: gets round problems caused by .Machine$integer.max

  Q<-round((tp+tn) / (tp+tn+fp+fn),3)    #Prediction Accuracy

  #Matthews Correlation Coefficient
  C_denominator<-sqrt((tp+fn)*(tp+fp)*(tn+fn)*(tn+fp))
  if(C_denominator==0){ C_denominator<-1}	
  C<-round((tp*tn-fn*fp)/C_denominator,3)       

  #Positive
  precision_p<-round(tp/(tp+fp),3)	#Precision,"Positive Predictive Value"
  sensitivity<-round(tp/(tp+fn),3)	#Recall,"True Positive Rate",Sensitivity, BBB+ accuracy (Li_2005)

  #Negative
  precision_n<-round(tn/(tn+fn),3)       #Precision,"Negative Predictive Value"
  specificity<-round(tn/(tn+fp),3)	#"True negative rate","Specificity", BBB- accuracy (Li_2005)

  ###############################################
  #
  # Useful references:
  # http://onlinestatbook.com/2/estimation/proportion_ci.html
  # http://stattrek.com/estimation/confidence-interval-proportion.aspx?Tutorial=AP 
  #
  # Note - first link suggests an additional (0.5/N) term in confidence interval calculation 
  #        to correct for the fact that a normal distribution is used to approximate a discrete distribution
  #
  ###############################################

  #Function: Compute standard error of a statistic that is defined as x/y.
  #x_y is the value of x/y (i.e. the value of the chosen statistic)
  #y is the value of the denominator
  standard_error_of_proportion<-function(x_y, y){
    return(round(sqrt((x_y*(1-x_y))/y),3))
  }

  #Standard errors
  se_precision_p<-standard_error_of_proportion(x_y=precision_p, y=(tp+fp))
  se_precision_n<-standard_error_of_proportion(x_y=precision_n, y=(tn+fn))
  se_sensitivity<-standard_error_of_proportion(x_y=sensitivity, y=(tp+fn))
  se_specificity<-standard_error_of_proportion(x_y=specificity, y=(tn+fp))

  #Confidence intervals
  #Using:
  #Z(90%)=1.645
  #Z(95%)=1.96
  #Z(99%)=2.575
  #x is the statistic for which CI is calculated, (e.g. sensitivity, specificity, etc)
  #se is the standard error of the statistic
  ci<-function(x, se){
    z90=1.645; z95=1.96; z99=2.575
    se99=se*z99; se95=se*z95; se90=se*z90
    return(list(ci99=round(c(x-se99,x+se99),3), ci95=round(c(x-se95,x+se95),3) , ci90=round(c(x-se90,x+se90),3)))
  }
  ci_precision_p<-ci(precision_p, se_precision_p)
  ci_precision_n<-ci(precision_n, se_precision_n)
  ci_sensitivity<-ci(sensitivity, se_sensitivity)
  ci_specificity<-ci(specificity, se_specificity)

  #return(list(tp=tp,tn=tn,fp=fp,fn=fn, Q=Q, C=C, precision_p=precision_p, sensitivity=sensitivity, precision_n=precision_n, specificity=specificity, info="Q is prediction accuracy and C is Matthews correlation coefficient.  See equations 11 and 12 in Li et al, JCIM, 2005. Sensitivity is called BBB+ accuracy, while Specificity is called BBB- accuracy by Li et al."))

  return(list(tp=tp,tn=tn,fp=fp,fn=fn, Q=Q, C=C, 
	precision_p=precision_p, sensitivity=sensitivity, precision_n=precision_n, specificity=specificity, 
	se_precision_p=se_precision_p, se_precision_n=se_precision_n, se_sensitivity=se_sensitivity, se_specificity=se_specificity,
        ci_precision_p=ci_precision_p, ci_precision_n=ci_precision_n, ci_sensitivity=ci_sensitivity, ci_specificity=ci_specificity)) 

}

#Function: normalize a vector
vecnorm<-function(x){
  l<-sqrt(sum(x^2))
  return(x/l)
}

#Function: centre a vector
center<-function(x){
  return(x-mean(x))
}

#Function: zero function
zero<-function(x){
  return(x-min(x))
}

#Function: convert a multi-class confusion matrix into a binary confusion matrix 
#"confusion" is an n-by-n multiclass confusion matrix (with predicted classes on the top and real classes on the left-hand-side) 
#"x" is the name of the column of interest from the multi-class confusion matrix.
#Output is a 2-by-2 confusion matrix with row and column names of c(x,"not-"x)
confusion_matrix_multiclass2binary<-function(confusion=NULL,x=NULL){
  x<-as.character(x) #To distinguish between column names (e.g. "1") and column numbers (e.g. 1)

  tmp<-which(colnames(confusion)==x)

  output<-matrix(,2,2)
  output[1,1]<-confusion[tmp,tmp]
  output[2,2]<-sum(confusion[-tmp,-tmp])
  output[1,2]<-sum(confusion[tmp,-tmp]) 
  output[2,1]<-sum(confusion[-tmp,tmp])

  colnames(output)<-c(x, paste("not-",x,sep=""))
  rownames(output)<-c(x, paste("not-",x,sep=""))

  return(output)
}

#Function: compute a multiclass confusion matrix from the real and predicted classifications
get_multiclass_confusion_matrix<-function(yexp=NULL, ypred=NULL){
  classes<-unique(yexp)
  output<-matrix(,length(classes), length(classes))
  colnames(output)<-classes
  rownames(output)<-classes
  for(pred in classes){
    for(exp in classes){
      output[pred,exp]<-sum(((ypred==pred) + (yexp==exp))==2)
    }
  }
  return(output) 
}

#########################################################################
#
# Pre-processing section
#
#########################################################################
bg.p <- bg_par
#Big wrapper function to do most preproc.
preprocess<-function(x,Mat,fs=1,bin=4,PC=10,sm=2,sm.p=2,bg=0,bg.p=4,norm=0)
{
 
  {
    ############################################# binning
    cat("binning... ")
    if (bin>=1)
    {
      Mat<-binning(x=x,Mat=Mat,by=bin)
      x<-as.numeric(colnames(Mat))
    }
    #library(signal)
    #library(rwt)
    #library(KernSmooth)
    #library(Peaks)
    #library.dynam('Peaks', 'Peaks', lib.loc=NULL)
    #library(hyperSpec)
    
    if (bg>2) 
    {der<-0}
    else
    {der<-bg}
    
    Nzeroes<-2^ceiling(log2(NCOL(Mat)))-NCOL(Mat) # for wavelet transform
    
    
    ############################################# PCA
    cat("smoothing... ")
    if (PC>0)       #PCA noise reduction
    {
      Dim<-dim(Mat)
      PCs<-1:as.integer(PC)
      PCA<-prcomp(Mat,center=TRUE) 
      Mat.new<-matrix(PCA$center,Dim[1],Dim[2],byrow = TRUE) 
      for (i in PCs)
      {
        for (j in 1:Dim[2])
        {
          Mat.new[,j]<-Mat.new[,j]+PCA$rotation[j,i]*PCA$x[,i]
        }
      }
      Mat<-Mat.new
    }
    ############################################# SMOOTHING and derivatisation
    Mat.new<-c()
    for (i in 1:NROW(Mat))
    {
      tmp<-Mat[i,]
      if (sm==1)        # SG filter
      {
        tmp<-sgolayfilt(tmp,p=sm.p,m=der)
      }
      else if (sm==2)    # Wavelet denoising filter
      {
        capture.output(suppressMessages(
          tmp <- as.vector(denoise.dwt(c(tmp,rep(0,times=Nzeroes)),daubcqf(sm.p*2)$h.0)$xd[1:NCOL(Mat),1])
        ))
        if (der==1)
        {
          tmp<-diff(tmp)
          x<-diff(x)/2+x[-length(x)]
        }
        else if (der==2)
        {
          tmp<-diff(tmp,differences = 2)
          x<-x[2:(length(x)-1)]
        } 
      }
      else if (sm==3)   # local polynomial using gaussian weighting as smoothing
      {
        TEMP<-locpoly(y=tmp,x=x,bandwidth=sm.p,gridsize=NCOL(Mat),drv=der)        
        tmp<-approx(x=TEMP$x,y=TEMP$y,xout=x)$y  
      }
      else 
      {
        if (der==1)
        {
          tmp<-diff(tmp)
          x<-diff(x)/2+x[-length(x)]
        }
        else if (der==2)
        {
          tmp<-diff(tmp,differences = 2)
          x<-x[2:(length(x)-1)]
        } 
      }
      Mat.new<-rbind(Mat.new,tmp)
      if(i%%100==0)
      {
        cat(".",sep="")
      }
    }
    #rownames(Mat.new)<-rownames(Mat)
    if (identical(length(x),NCOL(Mat.new)))
    {
      colnames(Mat.new)<-x
    }
    
    ############################################# normalisations
    cat("normalisations... ")
    amide<-which(x>1600 & x<1700)                 # for normalisation to amide I band 
    for (i in 1:NROW(Mat.new))
    {
      tmp<-Mat.new[i,]
      if(norm==1)
        tmp<-(tmp-min(tmp))/diff(range(tmp))
      else if(norm==2)
        tmp<-tmp/sqrt(sum(tmp^2))
      else if(norm==3)
        tmp<-tmp/max(Mat[i,amide])
      Mat.new[i,]<-tmp
    }
    #whilecount<-whilecount+1
    if (fs==3)   #This is not needed anymore. It used to sort out problems with gaps in the spectra (CO2)
    {
      Mat.new[,c(1,NCOL(Mat.new))]<-0
      Mat.out<-cbind(Mat.out,Mat.new)
      x<-x.raw[region2]
      Mat<-Mat.raw[,region2]
      cat("region 2 is processed ")
    }
    else
    {
      #Mat.out<-Mat.new
      Mat.new <- Mat.new
    }
    
    ############################################# background corrections
    cat("background correction... ")
    #I have set p to 3 and w to 9 as a standard after derivatives
    if(bg==1)
    {
      #Take first derivative
      
      firstderiv_output <- savitzkyGolay(Mat.new,p=3,w=5,m=1) 
     
      Mat.new <- firstderiv_output # Replacing ftr[,x1N] with 1st deriv data
      
      
      
      
      
      
      
      
    }
    if(bg==2)
    {
      #Take second derivative
      x<-as.numeric(gsub("X_", "", names(ftr)[x1N]))
      #ftr[,x1N]<-t(apply(ftr[,x1N], 1, D2ssy, x=x)) DSP derivative code
      #fte[,x1N]<-t(apply(fte[,x1N], 1, D2ssy, x=x))
      secderiv_output <- savitzkyGolay(Mat.new,p=3,w=5,m=2)
      Mat.new <- secderiv_output
      
      
    }
    
    if(bg==3)
    {
      #if (!bg.p==0)
      #{
       # Fac<-diff(range(Mat.new))/10000*bg.p
      #  Mat.new<-Mat.new+rep(((x-median(x))*Fac)^2,each=NROW(Mat.new))
      #}
      Mat.new<-Mat.new - rubberband(x=x,Mat=Mat.new,spline=TRUE,noise=bg.p*RBp)
    }
    else if(bg==4)
    {
      library(baseline)
      Mat.new<-baseline.modpolyfit(Mat.new,x,degree = bg.p)$corrected
    }
    
  
  }
  cat("done")
  return(Mat.new)
}


binning<-function (x,Mat, by = stop("reduction factor needed"), na.rm = FALSE, 
          ...) 
{
  n <- ceiling(length(x)/by)
  small <- NCOL(Mat)%%by
  #if (small != 0) 
  #  warning(paste(c("Last data point averages only ", small, 
  #                 " points.")))
  bin <- rep(seq_len(n), each = by, length.out = NCOL(Mat))
  na <- is.na(Mat)
  if ((na.rm > 0) && any(na)) {
    if (na.rm == 1) {
      na <- apply(!na, 1, tapply, bin, sum, na.rm = FALSE)
      Mat <- t(apply(Mat, 1, tapply, 
                              bin, sum, na.rm = TRUE)/na)
    }
    else {
      tmp <- t(apply(Mat, 1, tapply, bin, sum, 
                     na.rm = FALSE))
      tmp <- sweep(tmp, 2, rle(bin)$lengths, "/")
      na <- which(is.na(tmp), arr.ind = TRUE)
      bin <- split(x, bin)
      for (i in seq_len(nrow(na))) {
        tmp[na[i, 1], na[i, 2]] <- mean(Mat[na[i,1], bin[[na[i, 2]]]], na.rm = TRUE)
      }
      Mat <- tmp
    }
  }
  else {
    Mat <- t(apply(Mat, 1, tapply, bin, 
                            sum, na.rm = FALSE))
    Mat <- sweep(Mat, 2, rle(bin)$lengths, 
                          "/")
  }
  x <- as.numeric(tapply(x, bin, mean, 
                               na.rm = na.rm > 0))
  colnames(Mat)<-x
  return(Mat)
}

rubberband<-function (x, Mat, noise=noise, spline=TRUE, ...) #This spline setting may be something to change
{
  if (!is.matrix(Mat))
  {
    y<-Mat
    pts <- chull(x, y)
    imax <- which.max(pts)
    pts <- c(pts[-seq_len(imax)], pts[seq_len(imax)])
    neg <- which(diff(pts) < 0)
    pts <- c(length(y), pts[c(1, neg + 1)])
    pts <- sort(unique(pts))
    tmp <- approx(x = x[pts], y = y[pts], xout = x, method = "linear")$y
    if (spline) 
    {
      pts <- which(y <= tmp + noise)
      if (length(pts) > 3) 
        tmp <- predict(smooth.spline(x[pts], y[pts], 
                                     ...)$fit, x, 0)$y
      else tmp <- spline(x[pts], y[pts], xout = x)$y
    }
    Mat <- tmp
  }
  else
  {
    for (s in seq_len(nrow(Mat)))
    {
      pts <- chull(x, Mat[s, ])
      imax <- which.max(pts)
      pts <- c(pts[-seq_len(imax)], pts[seq_len(imax)])
      neg <- which(diff(pts) < 0)
      pts <- c(ncol(Mat), pts[c(1, neg + 1)])
      pts <- sort(unique(pts))
      tmp <- approx(x = x[pts], y = Mat[s, pts], xout = x, method = "linear")$y
      if (spline) 
      {
        pts <- which(Mat[s, ] <= tmp + noise)
        if (length(pts) > 3) 
          tmp <- predict(smooth.spline(x[pts], Mat[s, pts], 
                                       ...)$fit, x, 0)$y
        else tmp <- spline(x[pts], Mat[s, pts], xout = x)$y
      }
      Mat[s, ] <- tmp
    }
  }
  Mat
}









#########################################################################
#
# Main Program
#
#########################################################################

#Read data
ftr<-read.table(infile_tr, header=TRUE, sep=";")
fte<-read.table(infile_te, header=TRUE, sep=";")
#file.remove("input_testing_dataset.ssv")
#file.remove("input_training_dataset.ssv")
#Check
if(sum(names(ftr)==names(fte))!=dim(ftr)[2]){
  stop("ERROR!  Training and Testing datasets have different column names!")
}

#########################################################################
#
# 
#
#########################################################################

#Select columns containing frequency data
x1N<-grep("X_", names(ftr))

#Select spectral region for classification.
tmp <-as.numeric(gsub("X_","", names(ftr)[x1N])) >= min_wavenumber

tmp2<-as.numeric(gsub("X_","", names(ftr)[x1N])) <= max_wavenumber #& >= 2800
x1N<-x1N[which((tmp+tmp2)==2)]  #AND gate
if(chop_co2 ==TRUE){
tmp <-as.numeric(gsub("X_","", names(ftr)[x1N])) <= minCO2
tmp2<-as.numeric(gsub("X_","", names(ftr)[x1N])) >= maxCO2
x1N<-x1N[which((tmp+tmp2)==1)]  #OR gate
}
#x1N <- x1N[which(x1N != co2region)] ##DOESN'T SEEM TO ALWAYS WORK!
#Check
if(!isTRUE(all.equal(min_wavenumber, as.integer(min_wavenumber)))){
  stop("ERROR! min_wavenumber not a multiple of 100 (which means the x-axis in barplot will be wrong)")
}
if(!isTRUE(all.equal(max_wavenumber, as.integer(max_wavenumber)))){
  stop("ERROR! max_wavenumber not a multiple of 100 (which means the x-axis in barplot will be wrong)")
}
out_fte_preproc <-preprocess(x=(as.numeric(gsub("X_", "", names(fte)[x1N]))),Mat=fte[,x1N],bin=bin_factor,PC=PCNR,sm=smooth_choice,sm.p=smooth_par,norm=norm_choice,bg=bg_choice,bg.p=bg_par)

cat("TEST SET PREPROC DONE")
out_ftr_preproc <-preprocess(x=(as.numeric(gsub("X_", "", names(ftr)[x1N]))),Mat=ftr[,x1N],bin=bin_factor,PC=PCNR,sm=smooth_choice,sm.p=smooth_par,norm=norm_choice,bg=bg_choice,bg.p=bg_par)

cat("TRAIN SET PREPROC DONE")
#x1N <- colnames(out_ftr_preproc)
out_ftr_preproc <-data.frame(out_ftr_preproc)
out_fte_preproc <-data.frame(out_fte_preproc)
names(out_ftr_preproc) <- paste("X_",names(out_ftr_preproc), sep="")
names(out_ftr_preproc) <- gsub("X_X", "X_", names(out_ftr_preproc))
#ftr[,x1N] <- NULL
#fte[,x1N]
names(out_fte_preproc) <- paste("X_",names(out_ftr_preproc), sep="")
names(out_fte_preproc) <- gsub("X_X", "X_", names(out_ftr_preproc))
#Here I have just stuck the processed data onto the end of the starting dataset
ftr_plus <- cbind(ftr,out_ftr_preproc)
fte_plus <- cbind(fte,out_fte_preproc)
#This is so that x1N can be set as the end section of the dataset
x1N <- ((ncol(ftr)+1):ncol(ftr_plus))
#Set the names back
ftr <- ftr_plus
fte <- fte_plus

#Normalize data?
if(normalize){
  #Normalize
  ftr[,x1N]<-t(apply(ftr[,x1N], 1, vecnorm))
  fte[,x1N]<-t(apply(fte[,x1N], 1, vecnorm))

  #Centre data
  #ftr[,x1N]<-t(apply(ftr[,x1N], 1, center))
  #fte[,x1N]<-t(apply(fte[,x1N], 1, center))

}

#Zero data
if(minzero){
  #Zero data
  ftr[,x1N]<-t(apply(ftr[,x1N], 1, zero))
  fte[,x1N]<-t(apply(fte[,x1N], 1, zero))
}

#Reverse order of frequencies
x1N<-rev(x1N)
tmp<-min_wavenumber	#HACK
min_wavenumber<-max_wavenumber #HACK
max_wavenumber<-tmp #HACK

#################################################################
#
# Loop over all y variables
#
#################################################################

#Variables
for(yvar in yvariables){

  #Print
  cat("Working on y-variable ",yvar," ...\n\n")
  
  #Get column number
  y<-which(names(ftr)==yvar)
  
  #Rename classes (names are used in plot legends)
  if(yvar=="A_Cancer_v_NonCancer"){
    ftr[,y]<-gsub("1","cancer", ftr[,y])
    ftr[,y]<-gsub("2","non-cancer", ftr[,y])
    fte[,y]<-gsub("1","cancer", fte[,y])
    fte[,y]<-gsub("2","non-cancer", fte[,y])
  } else if(yvar=="B_Met.v_Brain_v_NonCancer"){
    ftr[,y]<-gsub("1","Met", ftr[,y])
    ftr[,y]<-gsub("2","GBM", ftr[,y])
    ftr[,y]<-gsub("3","Non-Cancer", ftr[,y])
    fte[,y]<-gsub("1","Met", fte[,y])
    fte[,y]<-gsub("2","GBM", fte[,y])
    fte[,y]<-gsub("3","Non-Cancer", fte[,y])
  # This gets rid of non-cancer class in met brain non.  Capital letters are 
  # there for tracking that they have actually gone
    ftr <- ftr[which(ftr[yvar]!="Non-Cancer"),]
    fte <- fte[which(fte[yvar]!="Non-Cancer"),]
  }
  #Set factor
  ftr[,y]<-as.factor(ftr[,y])
  fte[,y]<-as.factor(fte[,y])

  
  #################################################################y
  #
  # Select k-folds for cross-validation on training data (split by patient and stratified)
  #
  #################################################################
  
  #Get two-column matrix containing classification in column 1 (i.e. cancer .v. non-cancer) and patient number in column 2
  patients_classifications<-matrix(unlist(strsplit(unique(paste(ftr[,y],ftr[,patients], sep=";")), ";")), ,2, byrow=TRUE)
  
  #Select k-folds (by patient)
  kfolds_by_patients<-createFolds(patients_classifications[,1], k = k, list = TRUE, returnTrain = FALSE)
  
  #Select k-folds (by spectra)
  kfolds_by_spectra<-list()
  for(i in names(kfolds_by_patients)){
    for(j in kfolds_by_patients[[i]]){
      kfolds_by_spectra[[i]]<-c(kfolds_by_spectra[[i]], which(ftr[,patients]==patients_classifications[j,2]))
    }
  }
  
  #Check
  if(length(unlist(kfolds_by_spectra))!=length(unique(unlist(kfolds_by_spectra)))){
    stop("ERROR! Duplicate spectra found in kfolds_by_spectra")
  }
  
  #################################################################
  #
  # k-fold cross-validation 
  #
  #################################################################
  
  #k-fold cross-validation.   (Will only work properly if k-folds are independent, i.e. no duplication, no resampling)
  ycv<-vector(mode="character",length=dim(ftr)[1])
  for(i in names(kfolds_by_spectra)){
  
    #Debug
    cat("\nRunning cross-validation ...\n")
    print(i)
  
    #Define tr/te split for this k-fold
    tecv<-as.numeric(unlist(kfolds_by_spectra[[i]]))
    trcv<-(1:dim(ftr)[1])[-tecv]
  
    #Train RF on k-fold training split
    myrf<-randomForest(x=ftr[trcv,x1N], y=ftr[trcv,y], importance=FALSE)
  
    #Predict results for k-fold testing split
    ycv[tecv]<-as.character(predict(myrf,ftr[tecv,x1N]))
  }
  
  #################################################################
  #
  # Build RF model using all training data
  #
  # n.b. From help(randomForest): "confusion: (classification only) the confusion matrix of the prediction (based on OOB data)."
  #
  #################################################################
  
  #Train RF model
  myrf<-randomForest(x=ftr[,x1N], y=ftr[,y], importance=TRUE, keep.forest=TRUE) 
  
  ################
  #
  # Analyze each class
  #
  ################
  
  output<-matrix(,0,13)
  colnames(output)<-c(    "Class",
                          "pacCV", "mccCV", "sensCV","specCV", "pprCV", "nprCV",
                          "pacTS", "mccTS", "sensTS","specTS", "pprTS", "nprTS")
  
  myclasses<-colnames(myrf$confusion)
  myclasses<-myclasses[-length(myclasses)]
  
  for(myclass in myclasses){
  
    ################
    #
    # k-fold cross-validation analysis
    #
    ################
    
    #Get confusion matrix for test set predictions (and save it to file)
    myconf_cv_multiclass<-get_multiclass_confusion_matrix(yexp=ftr[,y], ypred=ycv)
    #write.table(myconf_cv_multiclass, paste(yvar,"_confusion_matrix_cv_", myclass,".txt",sep=""))
    
    #Reduce multiclass confusion matrix to binary
    myconf_cv<-confusion_matrix_multiclass2binary(confusion=myconf_cv_multiclass,x=myclass)
    tp<-myconf_cv[1,1]       #True positives
    tn<-myconf_cv[2,2]       #True negatives
    fp<-myconf_cv[2,1]       #False positives
    fn<-myconf_cv[1,2]       #False negatives
    mystats_cv<-cstats(tp=tp, tn=tn, fp=fp,fn=fn)
    
    #Print results to screen and to file
    #print_results(file=paste(yvar,"_results_rf_cv_",myclass,".txt", sep=""), mystats=mystats_cv,  ftr=ftr, fte=fte, myclass=myclass, yvar=yvar, what=paste(k,"-fold cross-validation", sep=""))
    
    ################
    #
    # Test Set analysis
    #
    ################
    
    #Get confusion matrix for test set predictions
    myconf_te_multiclass<-get_multiclass_confusion_matrix(yexp=fte[,y], ypred=predict(myrf,fte[,x1N]))
    #write.table(myconf_te_multiclass, paste(yvar,"_confusion_matrix_te_", myclass,".txt",sep=""))
    
    #Reduce multiclass confusion matrix to binary
    myconf_te<-confusion_matrix_multiclass2binary(confusion=myconf_te_multiclass,x=myclass)
    tp<-myconf_te[1,1]       #True positives
    tn<-myconf_te[2,2]       #True negatives
    fp<-myconf_te[2,1]       #False positives
    fn<-myconf_te[1,2]       #False negatives
    mystats_te<-cstats(tp=tp, tn=tn, fp=fp,fn=fn)
    
    #Print results to screen and to file
    #print_results(file=paste(yvar,"_results_rf_te_",myclass,".txt", sep=""), mystats=mystats_te,  ftr=ftr, fte=fte, myclass=myclass, yvar=yvar, what="Test set prediction")
    
    #Save results
    if(se_and_ci2){
      output<-rbind(output, c(	myclass,

				#Cross-validation
                              	mystats_cv$Q,mystats_cv$C,
				paste(mystats_cv$sensitivity," (",mystats_cv$se_sensitivity,") [",paste(mystats_cv$ci_sensitivity$ci95,collapse=","),"]", sep=""),
                                paste(mystats_cv$specificity," (",mystats_cv$se_specificity,") [",paste(mystats_cv$ci_specificity$ci95,collapse=","),"]", sep=""),
                                paste(mystats_cv$precision_p," (",mystats_cv$se_precision_p,") [",paste(mystats_cv$ci_precision_p$ci95,collapse=","),"]", sep=""),
                                paste(mystats_cv$precision_n," (",mystats_cv$se_precision_n,") [",paste(mystats_cv$ci_precision_n$ci95,collapse=","),"]", sep=""),


				#Test set                              
			    	mystats_te$Q,mystats_te$C,
                                paste(mystats_te$sensitivity," (",mystats_te$se_sensitivity,") [",paste(mystats_te$ci_sensitivity$ci95,collapse=","),"]", sep=""),
                                paste(mystats_te$specificity," (",mystats_te$se_specificity,") [",paste(mystats_te$ci_specificity$ci95,collapse=","),"]", sep=""),
                                paste(mystats_te$precision_p," (",mystats_te$se_precision_p,") [",paste(mystats_te$ci_precision_p$ci95,collapse=","),"]", sep=""),
                                paste(mystats_te$precision_n," (",mystats_te$se_precision_n,") [",paste(mystats_te$ci_precision_n$ci95,collapse=","),"]", sep="")))


    } else {
      output<-rbind(output, c(myclass,
                              mystats_cv$Q,mystats_cv$C,mystats_cv$sensitivity,mystats_cv$specificity,mystats_cv$precision_p,mystats_cv$precision_n,
                              mystats_te$Q,mystats_te$C,mystats_te$sensitivity,mystats_te$specificity,mystats_te$precision_p,mystats_te$precision_n))
    }
    
    #Print
    cat("------------\n")
    
  }
  
  #Write results and patients
  write.table(output, file=paste(yvar,"_results_ALL.txt",sep=""))
  write.table(fte$patients, file="test_patients.txt")
  write.table(ftr$pateints, file="train_patients.txt")
  #Test
  xresults<-xtable(output, caption="caption", label="label")
  xCV <- xresults[,1:7]
  xTS <- xresults[,c(1,8:13)]
  
  print(xresults, type="html", file=paste(yvar,"_results_ALL.html",sep=""), include.rownames=FALSE)
  #print(xCV, file=paste(yvar,"_results_CV.tex",sep=""), type="latex",row.names=FALSE)
  #print(xTS, type="latex", file=paste(yvar,"_results_TS.tex",sep=""), row.names=FALSE)
  #################################################################
  #
  # RF Importance   NOTE - Only looking at Gini for now
  #
  #################################################################
  suppressWarnings(
  for(importance_measure in list("MeanDecreaseGini")){
    
    #Check importance
    imp<-sort(myrf$importance[,importance_measure])
    write.table(imp, file=paste(yvar, "_random_forest_importance_", importance_measure, ".txt", sep=""))
    
    #Plot importance
    tmp<-order(as.numeric(gsub("X_","",names(imp))), decreasing=(sign(max_wavenumber-min_wavenumber)<0))
    imp2<-imp[tmp]
    x<-round(as.numeric(gsub("X_","",names(imp2))))
    pdf(file=paste(yvar, "_random_forest_importance_", importance_measure, ".pdf", sep=""), paper="special", width=16, height=8)
    par(cex=1.3,cex.axis=1.3, cex.lab=1.5,lwd=1, mar=c(4,5,1,1), oma=c(4,4.1,4,4)) 
    mybarplot2<-plot(-x,imp2,typ='h',xlab=expression(paste("Wavenumber (cm"^"-1",")")), ylab="RF Importance")#,xaxt='n')

    mylabels<-seq(min_wavenumber,max_wavenumber,100*sign(max_wavenumber-min_wavenumber))

    myticks<-seq(imp2[1], imp2[length(imp2)], length.out=length(mylabels))
    axis(1)
    axis(1,at=myticks,labels=mylabels)
    axis(1)
    graphics.off()
      
    #Plot importance and all spectra (colour spectra by classification)
    #if(min(ftr[,x1N])>0){  #HACKED. changed Dave's preproc="none" option.
    #if(preprocessing=="none"){  # bg_choice here gives different results that don't make sense for everything.  If set to 0, plots are on x axis. 
    #If set to higher, plot will be in middle of y axis (good for plain un-derivatised spectra). This depends on bg_choice of course.  Caveat emptor.
    if(bg_choice >=5){
      ylim_barplot<-c(0, max(imp2))
      ylim_plot<-c(0, max(ftr[,x1N]))
    } else {
      ylim_barplot<-c(-max(imp2), max(imp2))
      ylim_plot<-c(-max(abs(ftr[,x1N])), max(abs(ftr[,x1N])))
    }
    
    y<-which(names(ftr)==yvar)
    tmp<-as.character(ftr[,y])
    classes<-unique(tmp)
    mycolours<-c("red","blue","green","orange","purple", "grey","brown","yellow","pink","cyan")
    if(length(classes)<=10){
      for(i in 1:length(classes)){
        tmp[which(tmp==classes[i])]<-mycolours[i]
      }
    }
    if(imp_and_all_spectra ==TRUE){
    pdf(file=paste(yvar,"_random_forest_importance_", importance_measure, "_and_all_spectra_colour-coded.pdf", sep=""), paper="special", width=16, height=8)
    par(cex=1.3,cex.axis=1.3, cex.lab=1.5,lwd=1, mar=c(4,5,1,1), oma=c(4,4.1,4,4)) 
    
    mybarplot2<-plot(-x,imp2,typ='h',xlab=expression(paste("Wavenumber (cm"^"-1",")")), ylab="RF Importance")#,xaxt='n')

    mylabels<-seq(min_wavenumber,max_wavenumber,100*sign(max_wavenumber-min_wavenumber))

    myticks<-seq(imp2[1], imp2[length(imp2)], length.out=length(mylabels))
    axis(1,labels=mylabels,at=myticks)#,labels=mylabels)
    par(new=TRUE)
    plot(x,ftr[1,x1N], typ="l", axes=FALSE, ylab="", xlab="", ylim=ylim_plot,xlim=c(min_wavenumber,max_wavenumber)) #xlim=c(max(x),min(x)))
    axis(4)
    mtext("Normalised Intensity", side = 4, line = 0, outer = TRUE, padj=3.0, cex=2, adj=0.6)#, col="red") , adj=0.75)
    ys<-2:dim(ftr)[1]
    cols<-rainbow(length(ys))
    i<-1
    for(y in ys){
      points(x,na.omit(ftr[y,x1N]), typ="l", col=tmp[i])
      i=i+1
    }
    legend("topright", legend=classes, col=mycolours[1:length(classes)], lty=1, lwd=2)
    graphics.off()
    }
    
  ###Plot spectra only
  ###
    if(spectra ==TRUE){
    pdf(file=paste(yvar,"_", "spectra_only.pdf", sep=""), paper="special", width=16, height=8)
    par(cex=1.1,cex.axis=1.1, cex.lab=1.1,lwd=1, mar=c(1,2,1,1), oma=c(4,4.1,4,4)) 
    
    plot(x,ftr[1,x1N], typ="l", axes=FALSE, ylab="", xlab="")#, ylim=ylim_plot,xlim=c(min_wavenumber,max_wavenumber)) #xlim=c(max(x),min(x))
    axis(1,at=myticks,labels=mylabels)
    axis(1)
    axis(2)

    ys<-2:dim(ftr)[1]
    cols<-rainbow(length(ys))
    i<-1
    for(y in ys){
      points(x,ftr[y,x1N], typ="l", col=tmp[i])
      i=i+1
    }
    legend("topright", legend=classes, col=mycolours[1:length(classes)], lty=1, lwd=2)
    graphics.off()
    }
    if(average_spectra ==TRUE){
    #Plot average spectra
    
    #pdf(file=paste(yvar,"_", "mean_spectra.pdf", sep=""), paper="special", width=16, height=8) # Standard plotting procedure
    pdf(file=paste("bf",bin_factor,"_","pcnr",PCNR,"_","smooth",smooth_choice,"_","smoothpar",smooth_par,"_","norm",norm_choice,"_","baseline",bg_choice,"_","basepar",bg_par,"_", "mean_spectra.pdf", sep=""), paper="special", width=16, height=8) #BRS code diagnostics method
    par(cex=1.1,cex.axis=1.1, cex.lab=1.1,lwd=1, mar=c(1,2,1,1), oma=c(4,4.1,4,4)) 
    plot(x,ftr[1,x1N], typ="none", axes=FALSE, ylab="", xlab="", ylim=ylim_plot, xlim=c(min_wavenumber,max_wavenumber))#xlim=c(max(x),min(x)))
    axis(1,at=myticks,labels=mylabels)
    axis(1)
    axis(2)

    ys<-2:dim(ftr)[1]
    cols<-rainbow(length(ys))
    
     ii <-1
    for(j in myclasses){
      points(x,colMeans(ftr[which(ftr[yvar] ==j),x1N]), typ="l", col=mycolours[ii])
      ii=ii+1
      }
    legend("topright", legend=classes, col=mycolours[1:length(classes)], lty=1, lwd=2)
    graphics.off()
    }
    
    if(avg_and_gini ==TRUE){
      #pdf(file=paste(yvar,"_random_forest_importance_", importance_measure, "_and_all_spectra_colour-coded.pdf", sep=""), paper="special", width=16, height=8)
      pdf(file=paste("bf",bin_factor,"_","pcnr",PCNR,"_","smooth",smooth_choice,"_","smoothpar",smooth_par,"_","norm",norm_choice,"_","baseline",bg_choice,"_","basepar",bg_par,"_", "avg_and_gini.pdf", sep=""), paper="special", width=16, height=8) #BRS code diagnostics method
      par(cex=1.1,cex.axis=1.1, cex.lab=1.1,lwd=1, mar=c(1,2,1,1), oma=c(4,4.1,4,4)) 
      
      mybarplot2<-plot(-x,imp2,typ='h',xlab=expression(paste("Wavenumber (cm"^"-1",")")), ylab="RF Importance")#,xaxt='n')
      
      mylabels<-seq(min_wavenumber,max_wavenumber,100*sign(max_wavenumber-min_wavenumber))
      
      myticks<-seq(imp2[1], imp2[length(imp2)], length.out=length(mylabels))
      axis(1,labels=mylabels,at=myticks)#,labels=mylabels)
      par(new=TRUE)
      par(cex=1.1,cex.axis=1.1, cex.lab=1.1,lwd=1, mar=c(1,2,1,1), oma=c(4,4.1,4,4)) 
      plot(x,ftr[1,x1N], typ="none", axes=FALSE, ylab="", xlab="", ylim=ylim_plot, xlim=c(min_wavenumber,max_wavenumber))#xlim=c(max(x),min(x)))
      #axis(1,at=myticks,labels=mylabels)
      #axis(1)
      #axis(2)
      axis(4)
      mtext("Normalised Intensity", side = 4, line = 0, outer = TRUE, padj=3.0, cex=2, adj=0.6)#, col="red") , adj=0.75)
     
      ys<-2:dim(ftr)[1]
      cols<-rainbow(length(ys))
      
      ii <-1
      for(j in myclasses){
        points(x,colMeans(ftr[which(ftr[yvar] ==j),x1N]), typ="l", col=mycolours[ii])
        ii=ii+1
      }
      legend("topright", legend=classes, col=mycolours[1:length(classes)], lty=1, lwd=2)
      graphics.off()
    }
  })
  #Remove large datasets to save diskspace - CANCELLED - selective save instead
  #rm(fte,ftr,fte_plus,ftr_plus,out_fte_preproc,out_ftr_preproc)
  #Save state of current seed
  currentseed <- .Random.seed
  write.table(currentseed,file = "seed.txt")
  #Save workspace
  save.image(file = paste(yvar, "_rf.RData",sep=""))
  #save(list = c("myrf","currentseed"),file = paste(yvar, "_rf.RData",sep=""))
  
  
}
#Remove large input files generated by split script
#file.remove("input_testing_dataset.ssv")
#file.remove("input_training_dataset.ssv")
