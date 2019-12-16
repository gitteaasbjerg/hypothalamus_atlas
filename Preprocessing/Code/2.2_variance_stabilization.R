#!/usr/bin/env Rscript
library("yaml")
library('DESeq2')
yaml_path <- '/projects/gitte/hypothalamus_atlas/get_sample_ids_from_meta.yaml'
y <- read_yaml(yaml_path)
cat('------------------------------------------------\n            Variance stabilization\n------------------------------------------------\n')

make_expr <- function(vsd, meta){
  expr<-t(assay(vsd))
  expr<-expr[(match(rownames(meta), rownames(expr))),]
  datExpr<-as.matrix(expr)
  datExpr
}

data = read.csv(y[['output']][['data_path']],row.names = 1, header = TRUE)
meta = read.csv(y[['output']][['meta_path']],row.names = 3)


batches <- y[['batches']]
for (i in 1:length(batches)){
  meta[[batches[i]]] <- factor(meta[[batches[1]]])
}

all_data <- NULL
studies <- y[['studies']]
for (i in 1:length(studies)){
  sample_meta <- meta[which(meta$Sample_series_id==studies[[i]]$Sample_series_id),]
  sample_data <- data[,which(meta$Sample_series_id==studies[[i]]$Sample_series_id)]
  
  dds <- DESeqDataSetFromMatrix(as.matrix(sample_data), colData = sample_meta, design = ~1)
  dds<-DESeq(dds)
  
  vsd <- vst(dds, blind=TRUE)
  datExpr <- make_expr(vsd, sample_meta)
  
  all_data <- rbind(all_data, datExpr)
  print(paste0('Successfully stablilized variation in ',i,' out of ', length(studies), ' studies'))
}

write.csv(as.data.frame(all_data),y[['output']][['data_path']], quote=F)
print(paste0('The variation in ', length(studies),' studie(s) has been successfully stabilized'))

