#!/bin/bash
Rscript Scripts/get_sample_ids_from_meta.R Settings/get_sample_ids_from_meta.yaml 
Rscript Scripts/get_gene_idxs.R Outputs/sample_idx_from_meta.yaml 3
Rscript Scripts/WGCNA_outlier_removal.R
Rscript Scripts/batch_correction.R