#!/usr/bin/env bash

## download endomtrial cancer data from GEO (GSE111976)
cd /home/degan/Ov_Endo_MultiOmics/Inputs/ref_datasets && curl -O -J "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE111976&format=file&file=GSE111976%5Fct%5Fendo%5F10x%2Erds%2Egz"
curl -O -J "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111976/suppl/GSE111976_summary_10x_day_donor_ctype.csv.gz"

curl -O -J "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE49910&format=file"