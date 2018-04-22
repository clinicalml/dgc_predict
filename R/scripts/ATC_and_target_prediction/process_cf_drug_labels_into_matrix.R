GetKnownCorrectors2 = function(addNewCompounds=TRUE){
  # the integers are pubchem ids
  correctors = c(
    'VALPROIC ACID' = 3121,
    'CHLORZOXAZONE' = 2733,
    'NEOSTIGMINE BROMIDE' = 8246,
    'MIDODRINE' = 4195,
    'TRICHOSTATIN A' = 444732,
    'MS-275' = 4261,
    'GLAFENINE' = 3474,
    'IBUPROFEN' = 3672,
    'GLYCEROL' = 753,
    'PYRIDOSTIGMINE' = 4991,
    'VARDENAFIL' = 110634,
    'TADALAFIL' = 110635,
    'KM11060' = 1241327,
    'RDR01752' =        9566157,
    'MIGLUSTAT' = 51634,
    'CHLORAMPHENICOL' = 5959,
    'OUABAIN' = 439501,
    'DIGOXIN' = 2724385,
    'PIZOTIFEN' = 27400,
    'BIPERIDEN' = 2381,
    'MG-132' = 462382,
    'VORINOSTAT' = 5311,
    'SCRIPTAID' = 5186,
    'CURCUMIN' = 969516)
  if(addNewCompounds){
    correctors = c(correctors,
                   '15-DELTA PROSTAGLANDIN J2' = 5311211,
                   'AZACITIDINE' = 23760137,
                   'BRD-K94991378' = 6731789,
                   'CD1530' = 24868309,
                   'LDN-193189' = 25195294,
                   'MD-II-008-P' = 60194076,
                   'STROPHANTHIDIN' = 6185,
                   'WITHAFERIN-A' = 73707417)
  }
  return(correctors)
}

load('/Users/rhodos/Desktop/Research/LINCS/submission/data/metadata/tensor_annot.RData')

load(DataDir('expr/drug/tensor_features_for_drug_property_prediction_10cells_knn.RData'))
L = L[setdiff(names(L), c('allcell','pca200','pca978'))]

cf2 = GetKnownCorrectors2(addNewCompounds=TRUE)
idx = which(annot$pubchemIds %in% as.character(cf2)) %>%
  union(which(toupper(annot$pertName) %in% names(cf2)))
cfDrugs = annot$pertId[idx]

L_perts = rownames(L$mean$obs)
Y = data.frame(cf_drugs=ifelse(L_perts %in% cfDrugs, 1, 0))
rownames(Y) = L_perts

save(Y, file='~/Desktop/Box/CFDR/output/ML_drug_pred/cf_drug_labels.RData')
