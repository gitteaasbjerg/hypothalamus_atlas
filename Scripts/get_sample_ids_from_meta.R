#!/usr/bin/env Rscript
library("rhdf5")
library("yaml")

cmd_args = commandArgs(trailingOnly = TRUE)
yaml_path <- cmd_args[1]

y <- read_yaml(yaml_path)

fpath <- y[['input']][['file_path']]
h5_path <- y[['input']][['internal_path']]
meta <- h5read(fpath, h5_path)


print("Getting indicies for each study according to provided filters.")
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
  print(paste("Study ", s[['Sample_series_id']], "has ", length(reduced_idx), " samples."))
}

output_path <- y[['output']][['file_path']]
print(paste("Writing sample indicies to YAML file at ", output_path))
write_yaml(study_idx, output_path)
print("")
print("")
