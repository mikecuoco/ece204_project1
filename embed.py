#!/usr/bin/env python
# Created on: Apr 10, 2023 at 9:03:04 PM
__author__ = 'Michael Cuoco'

import numpy as np
import pandas as pd
from anndata import AnnData
from sklearn.preprocessing import StandardScaler

# import embedding functions
from sklearn.manifold import TSNE, MDS
from sklearn.decomposition import PCA
from umap import UMAP

# read expression data, put into AnnData object
expr = pd.read_csv('data/TCGA.HNSC.expression.txt', sep='\t', header=0, index_col=[0,1])
meta = pd.read_csv('data/TCGA.HNSC.metadata.txt', sep='\t', header=0, index_col=0)
meta["tissue_source_site"] = [s.split("-")[1] for s in meta.index.values] 

# add sample_id to metadata, reindex
meta = (
    expr
    .reset_index()[['patient_id','sample_id']]
    .join(meta, on='patient_id')
    .set_index(['patient_id','sample_id'])
)

# put into Annotated Data object for nice organization
adata = AnnData(expr, obs=meta, var=expr.columns.to_frame(name='gene'), dtype=np.float32)
adata.obs_names = expr.index.get_level_values(1)
adata.var_names = expr.columns

# remove zeros
adata = adata[:, adata.X.sum(axis=0) > 0]

# log(x+1) transform each value in matrix
adata.obsm['X_normed'] = np.log1p(adata.X)

# scale only, don't center (PCA centers internally)
adata.obsm["X_scaled"] = StandardScaler().fit(adata.obsm['X_normed'], with_mean=False)

# scale only, don't center (PCA centers internally)
adata.obsm["X_centered"] = StandardScaler().fit(adata.obsm['X_normed'], with_std=False)

# scale and center
adata.obsm['X_standardized'] = StandardScaler().fit(adata.obsm['X_normed'])

# embed data
# PCA
print("Running PCA...")
adata.obsm["pca"] = PCA().fit_transform(adata.obsm["X_scaled"]).astype(np.float32)

# MDS
print("Running MDS...")
adata.obsm["mds"] = MDS(normalized_stress='auto').fit_transform(adata.obsm['X_standardized']).astype(np.float32)

# t-SNE
print("Running t-SNE...")
adata.obsm["tsne"] = TSNE().fit_transform(adata.obsm['X_standardized']).astype(np.float32)

# UMAP
print("Running UMAP...")
adata.obsm["umap"] = UMAP().fit_transform(adata.obsm['X_standardized']).astype(np.float32)

# save data
adata._sanitize()
adata.write_h5ad('data/TCGA.HNSC.embedded.h5ad')