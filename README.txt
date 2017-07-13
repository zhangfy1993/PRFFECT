# PRFFECT

Please read the manual.

All packages should be installed using install.packages("<pkg_name>", dependencies = TRUE)

The archive trial_run contains files which are intended to be a demonstration.  Please use the example provided to test your R setup.  This folder is self-contained, and the whole process of set selection, pre-processing and classification can be done in the same directory.  Please extract this to a folder of your choice on your computer.  Please read the README in the extracted folder to begin the demonstration.

Archive example_input.zip contains a full spectral dataset and a training and test set created from that dataset.

Archive example_output.zip contains all output files associated with running the main program


If you have problems which are not solved by reading the manual, please contact benjamin.r.smith@strath.ac.uk



Bug reports: please contact benjamin.r.smith@strath.ac.uk


### PREREQUISITES ###

Please install the packages: "randomForest", "xtable", "caret", "prospectr", "signal", "rwt", "KernSmooth", "Peaks", "hyperSpec".

These should be installed using the dependencies = TRUE option, e.g. install.packages("randomForest", dependencies = TRUE)

### QUICK START GUIDE/DEMONSTRATION ###

STEP 1 Unzip the archive trial_run.zip .  This will produce a folder called trial_run within the main PRFFECT-master directory containing the following files: olive_oil_BRS.csv, select_training_and_test_sets.r, user_defined_input.R, PRFFECTv1.r.

STEP 2* Open R or RStudio, and ensure the correct directory (trial_run) is selected as the working directory, and the aforementioned files are present.

STEP 3 Run the set selection script by typing source("select_training_and_test_sets.r").  This will produce two new files, input_training_dataset.ssv and input_testing_dataset.ssv.  Optionally, tr_prop can be changed, to choose proportion of data used as a training set.

STEP 4 Click user_defined_input.R in the "files" pane (RStudio), or open in a text editor of your choice (R).

STEP 5 Examine and (optionally) choose a pre-processing type in user_defined_input.R.  Be sure to save any changes.

STEP 6 Run the main program by typing: source("PRFFECTv1.r")

STEP 7 When the program has finished, new files containing graphs, statistical results and importance information will appear.

STEP 8 The output files can then be examined, the html file is useful for quickly ascertaining the statistics of the classification.  See contents of 'example_output.zip' for a full example of expected output. 


*  To change directory in R, the setwd command is used. e.g. setwd("~/Downloads/PRFFECT-master/trial_run")

### QUICK WORKFLOW OVERVIEW OF DEMONSTRATION (> indicates an R command) ###


i.   Unzip trial_run.zip

ii.  > setwd("~/Downloads/PRFFECT-master/trial_run")

iii. > source("select_training_and_test_sets.r")

iv.  Choose pre-processing in user_defined_input.R, save changes.

v.   > source("PRFFECTv1.r")


### FOR YOUR OWN DATA ###

The simplest way to progress from this demonstration to using other datasets is to create a new directory containing:

yourdataset.csv
select_training_and_test_sets.r
user_defined_input.R
PRFFECTv1.r

A workflow similar to the one above can then be carried out.
See manual for expected formatting of dataset.






