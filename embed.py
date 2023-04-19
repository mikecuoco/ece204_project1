#!/usr/bin/env python
# Created on: Apr 10, 2023 at 9:03:04 PM
__author__ = 'Michael Cuoco'

import numpy as np
import pandas as pd
from anndata import AnnData
from sklearn.manifold import TSNE, MDS
from sklearn.decomposition import PCA
from umap import UMAP

def read_tcga():
    '''
    Reads TCGA HNSC expression data and metadata, 
    extracts tissue_source_site and adds as new column to metadata
    compiles all into into AnnData object
    '''

    # read expression data, put into AnnData object
    print("Reading data/TCGA.HNSC.expression.txt...")
    expr = pd.read_csv('data/TCGA.HNSC.expression.txt', sep='\t', header=0, index_col=[0,1])
    print("Reading data/TCGA.HNSC.metadata.txt...")
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
    print("Creating AnnData object...")
    adata = AnnData(expr, obs=meta, var=expr.columns.to_frame(name='gene'), dtype=np.float32)
    adata.obs_names = expr.index.get_level_values(1)
    adata.var_names = expr.columns
    adata.raw = adata.copy()

    return adata

def normalize(adata: AnnData):
    '''
    Remove unexpressed genes, log transform, scale and center
    '''
    import scanpy as sc

    # remove zeros
    print("Removing unexpressed genes...")
    sc.pp.filter_genes(adata, min_cells=1)

    # log(x+1) transform each value in matrix
    print("Log transforming...")
    sc.pp.log1p(adata)

    # scale and center
    print("Scaling and centering...")
    sc.pp.scale(adata, max_value=10, zero_center=True)

    return adata

if __name__ == '__main__':

    adata = read_tcga()
    adata = normalize(adata)

    for nc in range(1,11):
        print(f"Embedding with {nc} components...")
        adata.obsm[f"pca_{nc}"] = PCA(n_components=nc).fit_transform(adata.X).astype(np.float32)
        adata.obsm[f"mds_{nc}"] = MDS(n_components=nc, normalized_stress='auto').fit_transform(adata.X).astype(np.float32)
        if nc <= 3:
            adata.obsm[f"tsne_{nc}"] = TSNE(n_components=nc).fit_transform(adata.X).astype(np.float32)
        adata.obsm[f"umap_{nc}"] = UMAP(n_components=nc).fit_transform(adata.X).astype(np.float32)

    # do a full PCA
    print("Embedding with full PCA...")
    adata.obsm["pca"] = PCA().fit_transform(adata.X).astype(np.float32)

    # save data
    print("Saving to data/TCGA.HNSC.embedded.h5ad...")
    adata._sanitize()
    adata.write_h5ad('data/TCGA.HNSC.embedded.h5ad')