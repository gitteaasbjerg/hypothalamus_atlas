#!/usr/bin/env Rscript
library("yaml")
library('ggplot2')
library('ggrepel')

yaml_path <- '/projects/gitte/hypothalamus_atlas/get_sample_ids_from_meta.yaml'
y <- read_yaml(yaml_path)
cat('------------------------------------------------\n             Removing outliers \n------------------------------------------------\n')


data = read.csv(y[['output']][['data_path']],row.names = 1, header = TRUE)
meta = read.csv(y[['output']][['meta_path']],row.names = 1, header = TRUE)

outlier_idx = list()
studies <- y[['studies']]

for (i in 1:length(studies)){
  sample_meta <- data[,meta$Sample_series_id == studies[[i]]$Sample_series_id]
  pca <- as.data.frame(prcomp(t(sample_meta),rank=2)$x)
  pca$outlier <- scale(pca$PC1) > 3 | scale(pca$PC2) > 3
  
  outlier = row.names(pca[which(pca$outlier),])
  outlier_idx = c(outlier_idx,which(meta$Sample_geo_accession %in% outlier))
  
  p <- ggplot(pca,aes(x=PC1,y=PC2,color=outlier)) + scale_color_manual(values=c("black", "red")) + geom_point() +ggtitle(studies[[i]]$Sample_series_id) + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  p <- p + geom_label_repel(aes(label=ifelse(pca$outlier==TRUE,as.character(row.names(pca)),'')),
                            box.padding   = 0.35, 
                            point.padding = 0.5)
  
  print(paste0("Removed ", length(outlier), " outlier(s) from study ", studies[[i]]$Sample_series_id))
  ggsave(paste0(y[['output']][['output_path']],'/PCA_',studies[[i]]$Sample_series_id,'.pdf'), p, width = 10, height=10)
  
}

outlier_idx <- unlist(outlier_idx, recursive=FALSE)

meta <- meta[-outlier_idx,]
write.csv(meta, y[['output']][['meta_path']], quote=F)

data <- data[, -outlier_idx]
write.csv(data, y[['output']][['data_path']], quote=F)

print(paste0("Removed a total of ", length(outlier_idx)," outlier(s) across ", length(studies), ' studie(s). ',ncol(data),' studie(s) remains'))
