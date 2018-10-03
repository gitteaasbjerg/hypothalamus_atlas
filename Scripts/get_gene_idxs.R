library("rhdf5")
library("yaml")

h5_path <- "~/Projects/rotation_Tune/Data/mouse_matrix.h5"
internal_path <- "data/expression"
output_path <- "~/Projects/rotation_Tune/Outputs/gene_idx.yaml"


cmd_args = commandArgs(trailingOnly = TRUE)
yaml_path <- cmd_args[1]
cutoff_study <- as.integer(cmd_args[2])
cutoff_all <- as.integer(cmd_args[3])

print(paste("Genes must be non-zero in ", cutoff_study, " sample(s) per single study to be retained."))
print(paste("Genes must be non-zero in ", cutoff_all, " sample(s) for all studies to be retained."))

gene_idx <- list()
study_idx <- read_yaml(yaml_path)

for (i in 1:length(study_idx)){
    s_idx <- study_idx[[i]]
    n <- names(study_idx[i])
    data <- h5read(h5_path, internal_path, index=list(NULL, s_idx))
    H5close()
    data = log2(data+1)
    data = preprocessCore::normalize.quantiles(data)
    retain_idx <- which(apply(data, 1, function(c)sum(c!=0)) >= cutoff_study)
    gene_idx[[n]] <- retain_idx
    print(paste0("Study ", n, " has ", length(retain_idx), " genes that meet the single study cutoff."))
}

all_sample_idx <- Reduce(union, study_idx)
all_gene_idx <- Reduce(union, gene_idx)
data <- h5read(h5_path, internal_path, index=list(all_gene_idx, all_sample_idx))
data = log2(data+1)
data = preprocessCore::normalize.quantiles(data)
H5close()
retain_idx <- which(apply(data, 1, function(c)sum(c!=0)) >= cutoff_all)
all_gene_idx <- all_gene_idx[retain_idx]
for (i in 1:length(gene_idx)){
    g_idx <- gene_idx[[i]]
    n <- names(gene_idx[i])
    retain_idx <- which(g_idx %in% all_gene_idx)
    g_idx <- g_idx[retain_idx]
    gene_idx[[n]] <- g_idx
    print(paste0("Study ", n, " has ", length(g_idx), " genes that meet all studies cutoff."))
}


all_gene_idx <- Reduce(union, gene_idx)
print(paste(length(all_gene_idx), " genes survive per-study and all-study filtering."))
write_yaml(gene_idx, output_path)
