library("yaml")
library("rhdf5")

h5_path <- "~/Projects/rotation_Tune/Data/mouse_matrix.h5"
internal_path <- "data/expression"
sample_idx_path <- "~/Projects/rotation_Tune/Outputs/sample_idx_from_meta.yaml"
gene_idx_path <- "~/Projects/rotation_Tune/Outputs/gene_idx.yaml"
drop_sample_idx_path <- "~/Projects/rotation_Tune/Outputs/drop_sample_idx.yaml"
output_path <- "~/Projects/rotation_Tune/Outputs/batch_corrected_expression.csv"

study_idx <- read_yaml(sample_idx_path)
drop_sample_idx <- read_yaml(drop_sample_idx_path)
gene_idx <- read_yaml(gene_idx_path)

retain_sample_idx <- list()
for (i in 1:length(study_idx)){
    s_idx <- study_idx[[i]]
    s_d_idx <- drop_sample_idx[[i]]
    n <- names(study_idx[i])
    if (typeof(s_d_idx) == typeof(list())){
        retain_s_idx <- s_idx
    }
    else{
        throwaway_idx <- which(s_idx %in% s_d_idx)
        retain_s_idx <- s_idx[-(throwaway_idx)]
    }
    retain_sample_idx[[n]] <- retain_s_idx
}

retain_sample_idx <- Reduce(union, retain_sample_idx)
retain_gene_idx <- Reduce(union, gene_idx)

data <- h5read(h5_path, internal_path, index=list(retain_gene_idx, retain_sample_idx))
H5close()
data = log2(data+1)
data = preprocessCore::normalize.quantiles(data)
series <- h5read(h5_path, "meta/Sample_series_id", index=list(retain_sample_idx))
H5close()
gene_names <- h5read(h5_path, "meta/genes", index=list(retain_gene_idx))
H5close()
gsm_accession <- h5read(h5_path, "meta/Sample_geo_accession", index=list(retain_sample_idx))
H5close()
batchid = match(series, unique(series))
data_batchcorrected <- sva::ComBat(dat=data, batch=batchid, par.prior=TRUE, prior.plots=FALSE)
rownames(data_batchcorrected) = gene_names
colnames(data_batchcorrected) = gsm_accession
write.csv(data_batchcorrected, output_path, quote=F)