
# Run all necessary setup
testAll = TRUE # If true, runs all unit tests
source('R/init.R')

# Load all tensors (input and cross-validated)
print(sprintf('Loading data (warning: this takes a few minutes)'))
sz = 'large'
fullDataset = TRUE
loadMergeAndPred = TRUE
source('R/scripts/load_all_tensors.R')

# Generate figures
source('R/figures/FIGURE_1_L1000_cell_spec.R')
source('R/figures/FIGURE_3A_scatter.R')
source('R/figures/FIGURE_3B_DEG_ROC.R')
source('R/figures/FIGURE_3C_compare_cmap.R')
source('R/figures/FIGURE_4AB_cell_specific_drugs.R')
source('R/figures/FIGURE_4C_and_S10_ABC_preservation_of_cell_specificity.R')
source('R/figures/FIGURE_5A_gene_cor_heatmaps.R')
source('R/figures/FIGURE_5B_gene_cor_bargraphs.R')
source('R/figures/FIGURE_6B-D_entity_specific_accuracy.R')
source('R/figures/FIGURE_7B_tsize_results.R')
source('R/figures/FIGURE_7C_obs_density.R')
source('R/figures/FIGURE_8AB_drug_repurposing.R')
source('R/figures/FIGURE_S4_benchmarking.R')
source('R/figures/FIGURE_S5_tensor_DEG_method.R')
