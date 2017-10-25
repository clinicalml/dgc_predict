dgc_predict: This code contains the source for the manuscript entitled "Cell-specific prediction and application of drug-induced gene expression profiles" by Rachel Hodos, et al. (PSB 2018). It applies and evaluates a variety of methods to complete a partially-observed data tensor, e.g. comprising gene expression profiles corresponding to various drugs, applied in various cellular contexts. Note that 'dgc' stands for drug, gene, and cell.

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

*In R, change you working directory to this one, then run:
>> source('R/src/init.R')

*The following warnings are normal:
1: replacing previous import ‘igraph::union’ by ‘GSEABase::union’ when loading ‘HTSanalyzeR’ 
2: In ComputeAUC(est, labels) : Taking absolute value before computing AUC

*Note that if you don't want to run all unit tests, do:
>> testAll = FALSE
>> source('R/init.R')

*Then in R/scripts, see: 
example_complete_tensor_and_cv.R
define_data_tensors.R for the code that compiles the Characteristic Direction drug signatures and subsets into the four tensors referred to in the paper.
regenerate_all_figures_from_paper.R. This reruns some of the computations, but some have been precomputed to save time.

** FIGURES **
All figures will be automatically written to the plot directory (defined in your config file). 

** DATA **
The code interacts with lots of different datasets and files, all available at goo.gl/nTy8sH. After untarring, you should have a 'data' and 'results' directory. You can place these wherever you like, just make sure you point to these directories in the config.txt file.

The data tensors (both the small and large tensors referred to in the paper) are included, in data/tensors/[large/small].mat. 

The cross-validated tensors are available in results/<sz>/<sz>_tensor_cv_results.mat, where sz is either small or large.

The final completed tensors are available for the larger tensor, in results/large/final_pred. 
