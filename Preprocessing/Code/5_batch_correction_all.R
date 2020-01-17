#!/usr/bin/env Rscript
library('limma')
library("yaml")
library('ggplot2')
library('ggrepel')
#warning=FALSE
yaml_path <- '/projects/gitte/hypothalamus_atlas/get_sample_ids_from_meta.yaml'
y <- read_yaml(yaml_path)
studies <- y[['studies']]
cat('------------------------------------------------\n  Removing batch effects from expression data \n------------------------------------------------\n')


data = read.csv(y[['output']][['data_path']],row.names = 1, header = TRUE)
meta = read.csv(y[['output']][['meta_path']],row.names = 1, header = TRUE)

#Initiate batch correction
batches <- y[['batches']]
covariets <- y[['covariates']]


#If no batches are given use 'Sample_series_id' as default
if (is.null(batches)){
  batches <- 'Sample_series_id'
}

#Create PCA before batch correction
pca <- as.data.frame(prcomp(t(data),rank=2)$x)
pca$study <- as.character(meta[batches][,1])
p <- ggplot(pca,aes(x=PC1,y=PC2,color=study)) + geom_point() 
ggsave(paste0(y[['output']][['output_path']],'/PCA_pre_batch_correction','.pdf'), p, width = 10, height=10)


#Extract columns & rows relevant for batch correction
sample_covariets <- meta[covariets]
sample_batches <- meta[batches]

#Remove batches with only one distinct value
sample_batches <- sample_batches[sapply(sample_batches, function(col) length(unique(col))) > 1]

batch <- list()
bb    <- list()
model_matrix_b_input <- ''

#If batch effect is observed prepare for correction
if (length(sample_batches) > 0){
  for (b in 1:length(sample_batches)){
    bb            <- droplevels(sample_batches[[b]])
    contrasts(bb) <- contr.sum(levels(bb))
    batch[[b]]    <- bb
    
    #Create batch input-string for model matrix
    model_matrix_b_input <- paste0(model_matrix_b_input,'batch[[',b,']]+')
  }
}

#Create covariet input-string for model matrix if covariate is given
model_matrix_c_input <- ''
if (length(sample_covariets) > 0){
  for (i in 1:length(sample_covariets)){
    model_matrix_c_input <- paste0(model_matrix_c_input,'meta[,covariets[',i,']]+')
  }
}

#Combining input-string for model matrix
model_matrix_input <- paste0('~',model_matrix_b_input,model_matrix_c_input)
model_matrix_input <- substr(model_matrix_input,1,nchar(model_matrix_input)-1)

#Only run batch correction if needed
if (model_matrix_input == ''){
  print('No batch corrections has been made')
} else{ 
  cova            <- model.matrix(as.formula(model_matrix_input))
  design_matrix   <- matrix(1,ncol(data),1)
  covariates      <- cova[,-1]
  data_corrected  <- removeBatchEffect(data, covariates=cova, design=design_matrix)
}

#Create PCA after batch correction
pca1 <- as.data.frame(prcomp(t(data_corrected),rank=2)$x)
pca1$study <- as.character(meta['Sample_series_id'][,1])
p <- ggplot(pca1,aes(x=PC1,y=PC2,color=study)) + geom_point() 
ggsave(paste0(y[['output']][['output_path']],'/PCA_batch_corrected','.pdf'), p, width = 10, height=10)


write.csv(as.data.frame(data_corrected), y[['output']][['data_path']], quote=F)
print(paste0("Batch correction of ", length(sample_batches)," batches and ", length(sample_covariets), ' covariates has been successfully conducted on all ', length(studies), ' studies'))
