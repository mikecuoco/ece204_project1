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

## Analyze

Create your own `{name}_analysis.ipynb` to analyze the data
