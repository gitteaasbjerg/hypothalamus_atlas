#!/usr/bin/env Rscript

############## Preparing data ######################
library("rhdf5")
library("yaml")

cat('------------------------------------------------\nExtracting studies according to provided filters\n------------------------------------------------\n')
yaml_path <- '/projects/gitte/hypothalamus_atlas/get_sample_ids_from_meta.yaml'
y <- read_yaml(yaml_path)

h5_path <- y[['input']][['file_path']]
h5_path_meta <- y[['input']][['internal_path_meta']]
h5_path_exp <- y[['input']][['internal_path_exp']]

meta <- h5read(h5_path, h5_path_meta)

############### Getting study index ##################
studies <- y[['studies']]
study_idx <- list()
for (s in studies){
  num_filters <- length(s)
  name_filters <- names(s)
  temp_idx_list <- list()
  for (i in 1:num_filters){
    temp_idx_list[[name_filters[i]]] <- which(meta[[name_filters[i]]] %in% s[[i]])
  }
  reduced_idx <- Reduce(intersect, temp_idx_list)
  study_idx[[s[['Sample_series_id']]]] <- reduced_idx
  print(paste("Study", s[['Sample_series_id']], "has", length(reduced_idx), "sample(s)"))
}

all_sample_idx <- Reduce(union, study_idx)
print(paste0("Found ", length(all_sample_idx), " sample(s) in ", length(study_idx), " studie(s)"))

################ Creating data file ####################
data <- h5read(h5_path, h5_path_exp, index=list(NULL,all_sample_idx)) 
H5close()
gene_names <- h5read(h5_path, "meta/genes")
H5close()
gsm_accession <- h5read(h5_path, "meta/Sample_geo_accession", index=list(all_sample_idx))
H5close()

rownames(data) = gene_names
colnames(data) = gsm_accession

write.csv(data, y[['output']][['data_path']], quote=F)

############# Extracting new meta file #################

meta_extract <- meta[ c('Sample_series_id', 'Sample_geo_accession','Sample_platform_id','Sample_submission_date',
                        'Sample_last_update_date','total_reads','reads_aligned','Sample_title')]
meta_extract <- data.frame(lapply(meta_extract, function(x) x[all_sample_idx]))
write.csv(meta_extract, y[['output']][['meta_path']] , quote=F)