testAll = FALSE
makePlots = TRUE
source('R/init.R')
library(R.matlab)

sz = 'large'
fullDataset = TRUE
loadMergeAndPred = FALSE

if(fullDataset){
  # just load measured and CV tensors (takes a long time to run)
  source('R/scripts/load_all_tensors.R')
}else{
    tensors = list()
    tensors$meas = LoadTensorMat(DataDir('tensors/T50_1.RData'))
    
    matlab = StartMatlab()

    for(method in c('mean', 'mean2', 'knn', 'fa_lrtc')){
      if(method == 'fa_lrtc'){
        mthd = 'tensor'
      }else{mthd = method}
      tensors$cv[[mthd]] = CrossValidateTensor(matlab, tensors$meas, method, dataset, 
                                               nFolds = NumSigs(tensors$meas), 
                                               maxFolds = NumSigs(tensors$meas))
      dimnames(tensors$cv[[mthd]]) = dimnames(tensors$meas)
    }
    close(matlab)
    rm(matlab)
}

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
