library(ica)
# lib.loc = "/home/adakic/R/x86_64-pc-linux-gnu-library/3.5"
# but .libPaths() contains it already                         

# This script is accompanied with the sbatch script with the same name (with suffix _job)
# I created a separate script for each tested n_comp... so this could be simplified further 

# Function for 100 ICA reps, where the first rep is with the default diagonal starting rotation matrix,
# others start with a random rotation matrix. If faster run required increase the tolerance to 1e-5.
L_icafast_250 <- function(datmat, n_comp = 250, n_rep = 100, max_it = 900, tol_it = 1e-6){
  x <- list()
  x[[1]] <- icafast(datmat, nc = n_comp, maxit = max_it, tol = tol_it, Rmat = diag(n_comp))
  
  for(i in 2:n_rep) {
    x[[i]] <- icafast(datmat, nc = n_comp, maxit = max_it, tol = tol_it, Rmat = matrix(rnorm(n_comp^2), n_comp, n_comp))
    print(paste("rep", i))
  }
  return(x)
}


# I just re-use the the function above conditionally for running on various versions of data normalisation
# the best to run just one at the time and serialise on the cluster, as it takes 15-150 hours depending on the size of data

data.version <- "spikeF" # Change accordingly - "sumF" or other


if(data.version == "spikeF") {
  
  ### 1. ICA on TM data normalised in a way that the variability in total RNA across cell types was RETAINED:
  #
  # data on punim - change dir accordingly -----------------#
  exp_mat <- readRDS(file = "/data/cephfs/punim0613/icafast/data/tm.smartseq2.spikeF.mat.hvg.123.rds")
  exp_mat <- as.matrix(exp_mat)
  
  system.time(out <- L_icafast_250(datmat = exp_mat))
  
  # save on punim - change dir accordingly -----------------#
  saveRDS(out, file = "/data/cephfs/punim0613/icafast/output/ica_spikeF_list_250.rds")
  rm(exp_mat)
  
}


if(data.version == "sumF") {
  
  ### 2. ICA runs on TM data normalised in a way that the variability in total RNA across cell types was REMOVED:
  #
  # data on punim - change dir accordingly -----------------#
  exp_mat <- readRDS(file = "/data/cephfs/punim0613/icafast/data/tm.smartseq2.sumF.mat.hvg.123.rds")
  exp_mat <- as.matrix(exp_mat)
  
  system.time(out <- L_icafast_250(datmat = exp_mat))
  
  # save on punim - change dir accordingly -----------------#
  saveRDS(out, file = "/data/cephfs/punim0613/icafast/output/ica_sumF_list_250.rds")
  rm(exp_mat)
  
}
