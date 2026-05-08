
import celltypist
from celltypist import models
import scanpy as sc
import pandas as pd
from scipy import sparse, io
import numpy as np



genenames = pd.read_table("/data_storage/datasets/all_combo_celltypist_counts_names.tsv")
adata = sc.read_mtx("/data_storage/all_combo_celltypist_counts.tsv")
adata.var_names = genenames["x"]

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

model = models.Model.load(model = 'Immune_All_High.pkl')
predictions1 = celltypist.annotate(adata,model="Immune_All_High.pkl",majority_voting=True)
adata1 = predictions1.to_adata()

adata1.obs.to_csv("/data_storage/celltypist_predictions_high.csv")

model = models.Model.load(model = 'Immune_All_Low.pkl')
predictions2 = celltypist.annotate(adata,model="Immune_All_Low.pkl",majority_voting=True)
adata2 = predictions2.to_adata()

adata2.obs.to_csv("/data_storage/celltypist_predictions_low.csv")
