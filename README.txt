dgc_predict: This includes code, data and results for the manuscript entitled "Cell-specific prediction of drug-induced gene expression profiles" by Rachel Hodos, et al. Note that 'dgc' stands for drug, gene, and cell.

Author: Rachel Hodos
Date: February 13, 2017 
Email: hodos@cims.nyu.edu.

** CONFIGURATION INSTRUCTIONS **
You need to configure your own directory pointers in config/config.txt.

** OVERVIEW **
This code does the following:
- Compiles Characteristic Direction signatures into a tensor (in a configurable way) with dimensions of drugs x genes x cell types.
- Completes missing entries of the tensor using a variety of techniques, including the four outlined in the paper (see matlab/complete_tensor/CompleteTensor.m for all available methods).
- Runs cross validation to evaluate prediction accuracy
- Runs a variety of performance evaluations on prediction results, as described in the paper

All of the actual predictions are computed in MATLAB, and then all of the analysis of the results is in R. However, there are wrapper functions to call the main MATLAB functions (see R/scripts/example_complete_tensor_and_cv.R), and hence users should be able to stay in R if they wish.

Note that the tensor is assumed to have dimensions of drugs x genes x cell types, where the gene dimension is either completely observed or completely missing for each drug/cell combination.

** WHERE TO START **
In R/scripts, see: 
example_complete_tensor_and_cv.R
define_data_tensors.R for the code that compiles the Characteristic Direction drug signatures and subsets into the four tensors referred to in the paper.
regenerate_all_figures_from_paper.R, which as you can probably guess, generates all the figures from the paper. This reruns some of the computations, but a lot of it has been precomputed, in order to avoid extensive compute time.

** FIGURES **
All figures should be automatically written to the plot directory (defined in your config file). 

