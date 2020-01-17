#!/usr/bin/env Rscript
library('limma')
library("yaml")

yaml_path <- '/projects/gitte/hypothalamus_atlas/get_sample_ids_from_meta.yaml'
y <- read_yaml(yaml_path)
cat('------------------------------------------------\n  Removing batch effects from expression data \n------------------------------------------------\n')

data = read.csv(y[['output']][['data_path']],row.names = 1, header = TRUE)
meta = read.csv(y[['output']][['meta_path']],row.names = 1, header = TRUE)


batches <- y[['batches']]
covariets <- y[['covariates']]

all_data <- NULL

studies <- y[['studies']]
for (i in 1:length(studies)){
  study <- studies[[i]]$Sample_series_id
  
  #Extract relevant samples for study i
  sample_meta <- meta[which(meta$Sample_series_id==study),]
  sample_data <- data[,which(meta$Sample_series_id==study)]
  
  #Extract columns & rows relevant for batch correction
  sample_covariets <- sample_meta[covariets]
  sample_batches <- sample_meta[batches]
  
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
      model_matrix_c_input <- paste0(model_matrix_c_input,'sample_meta[,covariets[',i,']]+')
    }
  }

  #Combining input-string for model matrix
  model_matrix_input <- paste0('~',model_matrix_b_input,model_matrix_c_input)
  model_matrix_input <- substr(model_matrix_input,1,nchar(model_matrix_input)-1)

  if (model_matrix_input == ''){
    print('No batch corrections have been made for study',studies[[i]]$Sample_series_id)
  } else {
    formula         <- as.formula(model_matrix_input)
    cova            <- model.matrix(formula)
    design_matrix   <- matrix(1,ncol(sample_data),1)
    covariates      <- cova[,-1]
    data_corrected  <- removeBatchEffect(sample_data, covariates=cova, design=design_matrix)
    print(paste0("Batch correction of ", length(sample_batches)," batches and ", length(sample_covariets), ' covariates has been successfully conducted on study ',study))
    
  }
  all_data <- cbind(all_data, data_corrected)
}

write.csv(as.data.frame(all_data),y[['output']][['data_path']], quote=F)
