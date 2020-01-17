#!/usr/bin/env Rscript
library("yaml")

yaml_path <- '/projects/gitte/hypothalamus_atlas/get_sample_ids_from_meta.yaml'
y <- read_yaml(yaml_path)
cat('------------------------------------------------\n            Normalizing studie(s)\n------------------------------------------------\n')


data = read.csv(y[['output']][['data_path']],row.names = 1, header = TRUE)
meta = read.csv(y[['output']][['meta_path']],row.names = 3, header = TRUE)


all_data <- NULL

studies <- y[['studies']]
for (i in 1:length(studies)){
  sample_data <- data[,which(meta$Sample_series_id==studies[[i]]$Sample_series_id)]
  
  cname <- colnames(sample_data)
  rname <- rownames(sample_data)
  
  sample_data <- log2(sample_data+1)
  sample_data <- preprocessCore::normalize.quantiles(as.matrix(sapply(sample_data, as.numeric)))

  colnames(sample_data) <- cname
  rownames(sample_data) <- rname
  
  all_data <- cbind(all_data, sample_data)
  print(paste0('Successfully normalized ',i,' out of ', length(studies), ' studies'))
}

write.csv(as.data.frame(all_data),y[['output']][['data_path']], quote=F)
print(paste0('The ', length(studies),' studie(s) has been successfully normalized'))

