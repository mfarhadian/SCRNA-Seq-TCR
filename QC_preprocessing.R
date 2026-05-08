
# Load packages

library(dplyr)
library(readr)
library(Cairo)
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

# set WD
hd="/data_storage/Cellranger_RunOutput"
setwd(hd)

# Create QC directories
dir.create(paste0(hd,"/basic_qc_v5"))
dir.create(paste0(hd,"/basic_qc_v5/individual_cluster_plots_new_v5"))
dir.create(paste0(hd,"/basic_qc_v5/individual_qc_plots_new_v5"))
dir.create(paste0(hd,"/basic_qc_v5/individual_post_qc_datasets_new_v5"))
dir.create(paste0(hd,"/datasets_v5"))


# make empty df where we'll store QC stats
qc_df = data.frame()
start_time=Sys.time()

cluster_plots = list()
vln_plots = list()

# define function to do the following:
# 1 read in sample for the batch
# 2 calculate basic QC metrics
# 3 get rid of bad cells
# 4 label with donor ids and case status
# 5 get rid of individuals with unpaired samples


do_gex_qc = function(x){
  message("Reading in GEX data for...")
  message("Batch:", x)
  message("Starting at ",start_time)
  
  # run SoupX
  print(paste0(hd,"/",x,"/outs/multi/count"))
  raw = load10X(dataDir=paste0(hd,"/",x,"/outs/count")) 
  toc = Seurat::Read10X(data.dir=paste0(hd,"/",x,"/outs/count/filtered_feature_bc_matrix"))

  # manual method
  nonExpressedGeneList = list(HB = c("HBB", "HBA2","HBA1"), IG = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC"))
  nonzeroraw  = raw$toc[rowSums(raw$toc[,-1]) != 0,]
  hbs = c("HBB", "HBA2","HBA1"); hbsc = length(hbs[hbs %in% rownames(nonzeroraw) == T])
  if(hbsc  == 0 | x %in% c("Pool4")){
    nonExpressedGeneList = list(IG = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC"))
  }
  useToEst = estimateNonExpressingCells(raw, nonExpressedGeneList = nonExpressedGeneList)
  sc = calculateContaminationFraction(raw,nonExpressedGeneList,useToEst=useToEst)
  out = adjustCounts(sc)
  df = out
  
  print("Making Seurat object")
  df = df %>% CreateSeuratObject()
  
  # stash initial cell number
  no_cells = df@meta.data %>% nrow()
  rho = sc$metaData$rho[1]

  # define % MT genes and Hb
  df[["percent.mt"]] = PercentageFeatureSet(df,pattern="^MT-")
  df[["percent.hb"]] = PercentageFeatureSet(df,pattern="^HB-")

  # Histograms of UMA and Gene counts
  dh1 <- data.frame(UMI=as.numeric(as.character(df@meta.data[['nCount_RNA']])))
  dh2 <- data.frame(Feature=as.numeric(as.character(df@meta.data[['nFeature_RNA']])))
  dh <- cbind(dh1,dh2)
  cust_theme <- theme(text = element_text(family = "Arial", face="plain",size=12,color = "black"))
  h1 <- qplot(dh$'UMI', bins=100)+ cust_theme
  h2 <- qplot(dh$'Feature', bins=100)+cust_theme
  h3 <- qplot(dh[dh$'UMI' < 3000,]$UMI, bins=50)+cust_theme
  h4 <- qplot(dh[dh$'Feature' < 2000,]$Feature, bins=50)+cust_theme
  CairoPNG(paste0("./basic_qc_v5/individual_qc_plots_new_v5/",x,"_CountHistograms.png"),height=8,width=8,units="in",res=300)
  grid.arrange(h1,h2,h3,h4,nrow=2)
  dev.off()
  
  # nCount & MT
  df = subset(df, subset = nFeature_RNA > 100 & nCount_RNA > 1000 & percent.mt < 5)

  # stash cell number after basic qc
  post_mt_nfeature_filters_cells = df@meta.data %>% nrow()

  # initial seurat workflow
  df = df %>% SCTransform()  %>% RunPCA() %>% FindNeighbors(dims=1:30) %>% FindClusters() %>% RunUMAP(dims=1:30)
  
  # Run DoubletFinder
  # Homotypic Doublet Proportion Estimate
  homotypic.prop = modelHomotypic(df@meta.data$seurat_clusters)
  nExp_poi = round(0.075*nrow(df@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj = round(nExp_poi*(1-homotypic.prop))
  
  df =  doubletFinder(df, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE) 
  colnames(df@meta.data)[ncol(df@meta.data)] = "dbfinder_call"
  no_of_doublets_dbfinder = df@meta.data %>% filter(dbfinder_call=="Doublet") %>% nrow()
  df = subset(df, subset = dbfinder_call =="Singlet")
  
  # now do QC
  final_cells = df@meta.data %>% nrow()

  # store QC metrics
  qc_metrics = data.frame(rho,no_of_doublets_dbfinder,no_cells,post_mt_nfeature_filters_cells,final_cells)
  write_csv(qc_metrics,paste0("./basic_qc_v5/qc_metrics_",x,".csv"))

  now_time=Sys.time()
  message("Started at: ",start_time,". Currently: ", now_time)

  # now normalise for the batch using SC transform
  df = df %>% SCTransform(assay="RNA",vars.to.regress="percent.mt",method = "glmGamPoi",variable.features.n = 10000,return.only.var.genes = FALSE)
  
  # plots
  CairoPNG(paste0("./basic_qc_v5/individual_cluster_plots_new_v5/",x,"individual_cellranger_cluster_plots_qc_plot.png"),height=8,width=8,units="in",res=300)
  print(DimPlot(df,label=T)+ggtitle(paste0("Batch:",x)))
  dev.off()

  output_file_name = paste0("./basic_qc_v5/individual_post_qc_datasets_new_v5/post_qc_run_",x,".rds")
  saveRDS(df,output_file_name)
  
  
}

for (pool in c("Pool1","Pool2","Pool10")){
  do_gex_qc(pool)
}

# The rest of the sample must be add

