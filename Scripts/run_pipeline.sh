#!/bin/bash
if [ ! -d "Data" ]; then
    mkdir -p Data
    cd Data
    wget https://s3.amazonaws.com/mssm-seq-matrix/mouse_matrix.h5
    cd ..
fi
if [ ! -d "$Outputs" ]; then
    mkdir -p Outputs
fi
Rscript Scripts/get_sample_ids_from_meta.R Settings/get_sample_ids_from_meta.yaml
Rscript Scripts/get_gene_idxs.R Outputs/sample_idx_from_meta.yaml 3
Rscript Scripts/WGCNA_outlier_removal.R
Rscript Scripts/batch_correction.R