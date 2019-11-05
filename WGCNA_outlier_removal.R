library("yaml")
library("rhdf5")
library("preprocessCore")
library("WGCNA")

h5_path <- "~/Downloads/mouse_matrix.h5"
internal_path <- "data/expression"
sample_idx_path <- "~/Desktop/9th_semester/RNA-seq/hypothalamus_atlas/Outputs/sample_idx_from_meta.yaml"
gene_idx_path <- "~/Desktop/9th_semester/RNA-seq/hypothalamus_atlas/Outputs/gene_idx.yaml"
output_path <- "~/Desktop/9th_semester/RNA-seq/hypothalamus_atlas/Outputs/drop_sample_idx.yaml"

study_idx <- read_yaml(sample_idx_path)
gene_idx <- read_yaml(gene_idx_path)

drop_idx <- list()
for (i in 1:length(study_idx)){
    g_idx <- gene_idx[[i]]
    s_idx <- study_idx[[i]]
    n <- names(study_idx[i])
    print('-------')
    print(paste0("WGCNA outlier removal for study: ", n))
    data <- h5read(h5_path, internal_path, index=list(g_idx, s_idx))
    H5close()
    data = log2(data+1)
    data = preprocessCore::normalize.quantiles(data)
    sdout <- 2
    normadj <- (0.5+0.5*bicor(data, use='pairwise.complete.obs'))^2
    netsummary <- fundamentalNetworkConcepts(normadj)
    K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
    C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
    outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
    temp_idx <- which(outliers == TRUE)
    outlier_idx <- s_idx[temp_idx]
    print(paste("Found ", length(outlier_idx), " outliers."))
    drop_idx[[n]] <- outlier_idx
}
print("")
print("")
write_yaml(drop_idx, output_path)
