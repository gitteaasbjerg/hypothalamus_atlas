#!/bin/bash
if [ ! -d "Data" ]; then
mkdir -p Data
cd Data
wget https://s3.amazonaws.com/mssm-seq-matrix/mouse_matrix.h5
cd ..
fi
if [ ! -d "Outputs" ]; then
mkdir -p Outputs
fi
Rscript 1_extraction.R
Rscript 2_normalisation.R
Rscript 3_outlier_removal.R
Rscript 4_batch_correction.R