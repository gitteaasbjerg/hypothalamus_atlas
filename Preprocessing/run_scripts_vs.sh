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
Rscript Preprocessing/Code/1_extraction.R
Rscript Preprocessing/Code/2.2_variance_stabilization.R
Rscript Preprocessing/Code/3_outlier_removal.R
Rscript Preprocessing/Code/4_batch_correction.R