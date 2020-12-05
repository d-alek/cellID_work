library(tidyverse)
library(ica)

# Generate 10 permuted data matrices, run an ica on each for each n_comp... for two versions of input data
# Save vafs (variance accounted for) - this needed for scree-like plots to compare with vafs from real ica runs 

data.version <- "spikeF" ### CHANGE ACCORDINGLY ! ###  ("sumF" or other...)

#-----------------------------------------------------------------------------------------------------------#
if(data.version == "spikeF") {
  
  ### 1. TM data normalised in a way that the variability in total RNA across cell types was RETAINED:
  exp_mat <- readRDS(file = "/data/cephfs/punim0613/icafast/data/tm.smartseq2.spikeF.mat.hvg.123.rds")
  exp_mat <- as.matrix(exp_mat)
  
  facs_exp_null <- replicate(n = 10, expr = t(apply(exp_mat, 1, sample)), simplify = F) # list of 10 permuted mats
  rm(exp_mat)
  
  # do one ica run on each permuted matrix, across all ncomp values (total 10*5 runs)
  # note this time Rmat is set to default diagonal matrix so different ncomp runs are comparable
  # only vafs are saved - one matrix of 10 vafs (columns) for each ncomp run
  for(ncomp in c(50,100,150,200,250)) {
    
    ica.null <- map(.x = facs_exp_null, .f = icafast, nc = ncomp, maxit = 750, tol = 1e-6, Rmat = diag(ncomp))
    
    ica.null.vafs <- list()
    for(i in 1:10) {
      ica.null.vafs[[i]] <- ica.null[[i]]$vafs
    } 
    mat <- do.call(cbind, ica.null.vafs)
    rm(ica.null, ica.null.vafs)
    
    fname <- paste0("/data/cephfs/punim0613/icafast/output/null_ica_spikeF_", ncomp, ".rds")
    saveRDS(mat, file = fname)
    
    print(paste(ncomp, "done...")); print(Sys.time()); cat("\n")
    
  }
  
}

#-----------------------------------------------------------------------------------------------------------#
if(data.version == "sumF") {
  
  ### 2. TM data normalised in a way that the variability in total RNA across cell types was REMOVED:
  exp_mat <- readRDS(file = "/data/cephfs/punim0613/icafast/data/tm.smartseq2.sumF.mat.hvg.123.rds")
  exp_mat <- as.matrix(exp_mat)
  
  facs_exp_null <- replicate(n = 10, expr = t(apply(exp_mat, 1, sample)), simplify = F) # list of 10 permuted mats
  rm(exp_mat)
  
  # do one ica run on each permuted matrix, across all ncomp values (total 10*5 runs)
  # note this time Rmat is set to default diagonal matrix so different ncomp runs are comparable
  # only vafs are saved - one matrix of 10 vafs (columns) for each ncomp run
  for(ncomp in c(50,100,150,200,250)) {
    
    ica.null <- map(.x = facs_exp_null, .f = icafast, nc = ncomp, maxit = 750, tol = 1e-6, Rmat = diag(ncomp))
    
    ica.null.vafs <- list()
    for(i in 1:10) {
      ica.null.vafs[[i]] <- ica.null[[i]]$vafs
    } 
    mat <- do.call(cbind, ica.null.vafs)
    rm(ica.null, ica.null.vafs)
    
    fname <- paste0("/data/cephfs/punim0613/icafast/output/null_ica_sumF_", ncomp, ".rds")
    saveRDS(mat, file = fname)
    
    print(paste(ncomp, "done...")); print(Sys.time()); cat("\n")
    
  }
  
}
  
  