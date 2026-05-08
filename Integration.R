

# Load packages
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Rcpp)
library(harmony)
library(generics)
library(Seurat)
library(SoupX)
library(gridExtra)
library(sctransform)
library(celldex)
library(SingleR)
library(DoubletFinder)
library(glmGamPoi)
library(generics)
library(pracma)
library(fields)

wd="/GitHub.com/SCRNA/Batch_11/QC/" 

# set WD
setwd(paste0(wd,"/basic_qc/individual_post_qc_datasets_new_v5/")) 


post_qc_datasets = lapply(files,function(x){
  message("reading in ",x)
  df = readRDS(x)
  return(df)
})

all_combo = merge(x = post_qc_datasets[[1]], y = post_qc_datasets[2:length(post_qc_datasets)])

# Integrate datasets
var.features = SelectIntegrationFeatures(post_qc_datasets, nfeatures = 5000)
VariableFeatures(all_combo) = var.features

# Run PCA
all_combo = all_combo %>% RunPCA(assay.use="SCT")
print("PCA done)")

# Plot PCA reduction
ElbowPlot(all_combo)

#Run Harmony 
gc(full=T)
all_combo = all_combo %>% RunHarmony(assay.use="SCT",group.by.vars="Pool_id")
print("Harmony done")

head(all_combo@meta.data)
# save

saveRDS(all_combo, file = paste0(wd,"/datasets/individual_post_qc_datasets_new_v5_integrated.rds"))


