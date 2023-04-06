#!/usr/bin/env bash
# Author: Mike Cuoco
# Created on: 4/6/23, 1:22 PM

# exit if any non-zero, exit if undefined var
set -euo pipefail

mkdir -p data && cd data

printf "Downloading TCGA.HNSC.expression.txt..."
curl -L -s -o TCGA.HNSC.expression.txt https://www.dropbox.com/sh/57resib0deyyb2s/AAAOvQLe1m_dSq8updEyAwt8a/Team_2_HNSC/TCGA.HNSC.expression.txt?dl=0
printf "Done\n"

printf "Downloading TCGA.HNSC.metadata.txt..."
curl -L -s -o TCGA.HNSC.metadata.txt https://www.dropbox.com/sh/57resib0deyyb2s/AACXagtvlw-ig6LQK1_dMmK0a/Team_2_HNSC/TCGA.HNSC.metadata.txt?dl=0
printf "Done\n"

printf "Downloading TCGA.HNSC.mutations.txt..."
curl -L -s -o TCGA.HNSC.mutations.txt https://www.dropbox.com/sh/57resib0deyyb2s/AABc3MZBMis2QTuq731-hBIha/Team_2_HNSC/TCGA.HNSC.mutations.txt?dl=0
printf "Done\n"

