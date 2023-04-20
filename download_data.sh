#!/usr/bin/env bash
# Author: Mike Cuoco
# Created on: 4/6/23, 1:22 PM

# exit if any non-zero, exit if undefined var

mkdir -p data && cd data

printf "Downloading TCGA.HNSC.expression.txt..."
wget https://www.dropbox.com/sh/57resib0deyyb2s/AAAOvQLe1m_dSq8updEyAwt8a/Team_2_HNSC/TCGA.HNSC.expression.txt
printf "Done\n"

printf "Downloading TCGA.HNSC.metadata.txt..."
wget https://www.dropbox.com/sh/57resib0deyyb2s/AACXagtvlw-ig6LQK1_dMmK0a/Team_2_HNSC/TCGA.HNSC.metadata.txt
printf "Done\n"

printf "Downloading TCGA.HNSC.mutations.txt..."
wget https://www.dropbox.com/sh/57resib0deyyb2s/AABc3MZBMis2QTuq731-hBIha/Team_2_HNSC/TCGA.HNSC.mutations.txt
printf "Done\n"

