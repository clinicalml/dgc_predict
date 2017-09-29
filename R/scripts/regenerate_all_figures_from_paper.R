
# This generates all figures except for Figure 5 (see Rmd file).

# Load all tensors (input and cross-validated)
if(!exists('tensors')){
  tensors = LoadTensors('large', print=TRUE)
}

# Generate figures
source('R/figures/FIGURE_1_L1000_cell_spec.R')
source('R/figures/FIGURE_3A_scatter.R')
source('R/figures/FIGURE_3B_best_method_per_profile.R')
source('R/figures/FIGURE_3C_obs_density.R')
source('R/figures/FIGURE_3D_DEG_ROC.R')
source('R/figures/FIGURE_4_and_S2_cell_specific_drugs.R')
source('R/figures/SUPP_FIGURE_S1_benchmarking.R')
source('R/figures/SUPP_FIGURE_S3_chem_bias.R')