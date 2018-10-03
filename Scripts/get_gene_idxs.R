library("rhdf5")
library("yaml")

h5_path <- "~/Projects/rotation_Tune/Data/mouse_matrix.h5"
internal_path <- "data/expression"
output_path <- "~/Projects/rotation_Tune/Outputs/gene_idx.yaml"


cmd_args = commandArgs(trailingOnly = TRUE)
yaml_path <- cmd_args[1]
cutoff <- cmd_args[2]

print(paste("Genes must be non-zero in ", cutoff, " samples to be retained."))

gene_idx <- list()
study_idx <- read_yaml(yaml_path)

for (i in 1:length(study_idx)){
    s_idx <- study_idx[[i]]
    n <- names(study_idx[i])
    data <- h5read(h5_path, internal_path, index=list(NULL, s_idx))
    H5close()
    data = log2(data+1)
    data = preprocessCore::normalize.quantiles(data)
    retain_idx <- which(apply(data, 1, function(c)sum(c!=0)) >= cutoff)
    gene_idx[[n]] <- retain_idx
    print(paste0("Study ", n, " has ", length(retain_idx), " genes that meet cutoff."))
}

write_yaml(gene_idx, output_path)

