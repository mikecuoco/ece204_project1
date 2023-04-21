# BENG 285/BNFO 285/ECE 204. Statistical Learning in Bioinformatics: Project 1

## Setup

Clone this repo, create a virtual python environment, and install the required python libraries

```bash
git clone https://github.com/mikecuoco/ece204_project1.git
cd ece204_project1
python -m venv .venv # must be python 3.9
source .venv/bin/activate
pip install -r requirements.txt
```

If you need additional python libraries for the project, add them to this file

## Download the data

To download the data (143 Mb) from Dropbox, run the following command

```bash
bash download_data.sh
```

This will create the following data to the `data` directory.

```bash
data/
├── TCGA.HNSC.expression.txt
├── TCGA.HNSC.metadata.txt
└── TCGA.HNSC.mutations.txt
```

## Generate the embedded data

To generate the embedded data, run the following command

```bash
python embed.py
```

This will create the following data to the `data` directory.

```bash
data/
└── TCGA.HNSC.embedded.h5ad
```

## Analyze

Create your own jupyter notebook and select the python installation in `.venv` as the kernel to analyze the data

To read in the embedded data, use the following code.

```python
from anndata import read_h5ad
adata = read_h5ad('data/TCGA.HNSC.embedded.h5ad')
adata.X # log-normalized and standardized expression matrix
adata.obs # see metadata
adata.obsm['pca'] # PCA
adata.obsm['mds'] # MDS
adata.obsm['umap'] # UMAP
adata.obsm['tsne'] # t-SNE
```

[Anndata](https://anndata.readthedocs.io/en/latest/index.html) is a handy library for working with a data matrix and metadata. 
