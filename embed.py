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

	print("Cleaning metadata...")
	meta["tissue_source_site"] = [s.split("-")[1] for s in meta.index.values] 
	
	# add sample_id to metadata, reindex
	meta = (
		expr
		.reset_index()[['patient_id','sample_id']]
		.join(meta, on='patient_id')
	)

	sample_type = []
	for s in meta["sample_id"].values:
		if s.split("-")[3][:2] == "01":
			sample_type.append("Primary Solid Tumor")
		elif s.split("-")[3][:2] == "11":
			sample_type.append("Solid Tissue Normal")
		elif s.split("-")[3][:2] == "06":
			sample_type.append("Metastatic")
		else:
			raise ValueError("Unknown sample type")
		
	hist_grade_num = []
	for g in meta["histological_grade"].values:
		if g == "GX":
			hist_grade_num.append(np.nan)
		elif g == "G1":
			hist_grade_num.append(1)
		elif g == "G2":
			hist_grade_num.append(2)
		elif g == "G3":
			hist_grade_num.append(3)
		elif g == "G4":
			hist_grade_num.append(4)
		else:
			hist_grade_num.append(np.nan)

	meta["race"] = [r if r in ["WHITE", "BLACK OR AFRICAN AMERICAN", "ASIAN", "AMERICAN INDIAN OR ALASKA NATIVE"] else "UNKNOWN" for r in meta["race"]]
	meta["sample_type"] = sample_type
	meta["histological_grade_num"] = hist_grade_num
	meta.set_index(["patient_id","sample_id"], inplace=True)

	for i in ["ajcc_pathologic_tumor_stage", "histological_type", "histological_grade", "histological_grade_num"]:
		meta[i] = meta.apply(lambda x: x[i] if x["sample_type"] == "Primary Solid Tumor" else np.nan, axis=1)

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
	adata.obsm["mds"] = MDS().fit_transform(adata.X).astype(np.float32)
	adata.obsm["tsne"] = TSNE().fit_transform(adata.X).astype(np.float32)
	adata.obsm["umap"] = UMAP().fit_transform(adata.X).astype(np.float32)
	# save data
	print("Saving to data/TCGA.HNSC.embedded.h5ad...")
	adata._sanitize()
	adata.write_h5ad('data/TCGA.HNSC.embedded.h5ad')