# BENG 285/BNFO 285/ECE 204. Statistical Learning in Bioinformatics: Project 1

## Setup

Clone this repo, create a virtual python environment, and install the required python libraries

```bash
git clone https://github.com/mikecuoco/ece204_project1.git
python3.9 -m venv .venv # must be python >=3.7,<3.11
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

Create your own `{name}_analysis.ipynb` and select the python installation in `.venv` as the kernel to analyze the data

To read in the embedded data, use the following code.

```python
from anndata import read_h5ad
adata = read_h5ad('data/TCGA.HNSC.embedded.h5ad')
adata.obsm['pca'] # PCA embedding
```

[Anndata](https://anndata.readthedocs.io/en/latest/index.html) is a handy library for working with a data matrix and metadata. See [`./explore.ipynb`](./explore.ipynb) for how I used it.
