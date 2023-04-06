# BENG 285/BNFO 285/ECE 204. Statistical Learning in Bioinformatics - Project 1

## Setup

Create a virtual python environment for this project and install the required python libraries

```bash
python3 -m venv .venv 
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

## Analyze
