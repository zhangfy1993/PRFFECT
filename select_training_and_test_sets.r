########################################################################
#
# Split IR Spectra dataset into training and testing dataset
# 
#
#########################################################################

#Libraries
library(caret)

#Read data - original full dataset
f<-read.table("olive_oil_BRS.csv", header=TRUE, sep=",")

#Variables
patients<-"sample_number"   #Name of column containing patient IDs, sample numbers etc. 

#yvar sets which column to split data by - aiming to have a proportional split of classes in that column between tr and te
yvar="country"
trprop <- 0.75

##############################################
x1N<-grep("X_", names(f))       #Select columns containing frequency data - i.e. those with the X_ tag
y<-which(names(f)==yvar)
 
#Get two-column matrix containing classification in column 1 (i.e. cancer .v. non-cancer) and patient number in column 2
patients_classifications<-matrix(unlist(strsplit(unique(paste(f[,y],f[,patients], sep=";")), ";")), ,2, byrow=TRUE)

#Select training and testing datasets (by patients) - CHANGE p= VALUE TO SET TRAINING SET SIZE (proportion)
tr_by_patients<-as.numeric(unlist(createDataPartition(patients_classifications[,1], times=1,p=trprop,list=TRUE)))
te_by_patients<-(1:dim(patients_classifications)[1])[-tr_by_patients]

#Select training and testing datasets (by ID)
tr<-vector()
for(i in tr_by_patients){
  patient<-patients_classifications[i,2]
  tr<-c(tr,which(f[,patients]==patient))
}
te<-vector()
for(i in te_by_patients){
  patient_no<-patients_classifications[i,2]
  te<-c(te,which(f[,patients]==patient_no))
}

#Check
if(length(c(tr,te))!=length(unique(c(tr,te)))){
  stop("ERROR!  Duplicates in indices for training and test sets")
}

#Write data
write.table(f[tr,], file="input_training_dataset.ssv", row.names=FALSE, col.names=TRUE, sep=";", quote=FALSE)
write.table(f[te,], file="input_testing_dataset.ssv", row.names=FALSE, col.names=TRUE, sep=";", quote=FALSE)


#Check that each item in each class is represented properly in both training and test sets. Generates on-screen report.
sink("select_training_and_test_sets.log", split=TRUE)
for(i in c(yvar)){
  cat("Class: ",i,"\n")
  for(j in unique(f[,i])){
  
    itr<-sum(f[tr,i]==j)
    ite<-sum(f[te,i]==j)
    cat("  Item: ", j, "\n")
    cat("    Training: ", itr, "\n")
    cat("    Testing: ", ite, "\n")
    cat("    Fraction (tr): ", round(itr/(itr+ite),3), "\n")
  }
  cat("--\n")
}
sink(NULL)

