dgc_predict: This includes code, data and results for the manuscript entitled "Cell-specific prediction of drug-induced gene expression profiles" by Rachel Hodos, et al. Note that 'dgc' stands for drug, gene, and cell.

Author: Rachel Hodos
Date: September 28, 2017 
Email: rachelhodos@gmail.com.

** CONFIGURATION INSTRUCTIONS **
You need to configure your own directory pointers in config.txt.

** OVERVIEW **
This code does the following:
- Compiles LINCS L1000 Characteristic Direction signatures into a tensor (in a configurable way) with dimensions of drugs x genes x cell types.
- Completes missing entries of the tensor using a variety of techniques, including the four outlined in the paper (see matlab/complete_tensor/CompleteTensor.m for all available methods).
- Runs cross validation to evaluate prediction accuracy
- Runs a variety of performance evaluations on prediction results, as described in the paper

All of the actual predictions are computed in MATLAB, and then all of the analysis of the results is in R. However, there are wrapper functions to call the main MATLAB functions (see R/scripts/example_complete_tensor_and_cv.R) from R.

Note that the tensor is assumed to have dimensions of drugs x genes x cell types, where the gene dimension is either completely observed or completely missing for each drug/cell combination.

** WHERE TO START **

In R, change you working directory to this one, then run:
>> source('R/src/init.R')

Note that if you don't want to run all unit tests, do:
>> testAll = FALSE
>> source('R/src/init.R')

Then in R/scripts, see: 
example_complete_tensor_and_cv.R
define_data_tensors.R for the code that compiles the Characteristic Direction drug signatures and subsets into the four tensors referred to in the paper.
regenerate_all_figures_from_paper.R. This reruns some of the computations, but some has been precomputed to save time.

** FIGURES **
All figures should be automatically written to the plot directory (defined in your config file). 
