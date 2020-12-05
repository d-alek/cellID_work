# general ICA extractor code (combination of cluster and R session work)

library(tidyverse)
library(dplyr)
library(magrittr)
library(stringr)
library(ggplot2)
library(RColorBrewer)

library(Seurat)

# FACS dataset - SCT
facs.exp <- readRDS("data/facs.exp.707.rds") # 31,438 cells; 22,921 genes (exclused: 0 exp, ERCC, transgene)
facs.exp.sct <- SCTransform(facs.exp, variable.features.n = 7500, return.only.var.genes = TRUE)
saveRDS(facs.exp.sct, "/Users/adakic/projects/icatest/data/facs.exp.sct.707.rds")

# RESID
facs.exp.sct.mat.RESID <- as.matrix(GetAssayData(facs.exp.sct, assay = "SCT", slot = "scale.data")) # 1.8 gb
saveRDS(facs.exp.sct.mat.RESID, file = "/Users/adakic/projects/icatest/data/facs.exp.sct.7500hvg.mat.RESID.rds") # 7500 31438

# likewse for COUNT matrix - but now developing code for RESID option

# then run L_icafast_xxx.r and ica_loop_job_xxx_r100.sh functions on the cluster
# this results in ica_list_RESID_xxx.rds files containing outputs of 100 ica runs (each is a list of 12)

# from there extract S and A matrices and vafs vectors - first develop the code for 50 comps, then functionalise it for 50-400

#--------------------------------------------ON CLUSTER---#

library(stringr)
library(cluster)

# Based on the code below I made script for the cluster that extracts S and A matrices, and vafs from ica runs
# then reorients them consistently, derives Spearman distances for clustering, does pam clustering of components 
# and saves everything for ncomp = 50, 100, 150, 200, 250, 300, 350, 400
# in: ica_extractor_func.r and ica_extractor_job.sh ~ 24 job for 7.5k genes, 30k cells (can be parallelised if needed)
# just change names of the in/out files (depending on which data normalisatin protocol is taken) 

facs_resid_ica_list_50 <- readRDS("/data/cephfs/punim0613/icafast/output/ica_list_RESID_50.rds") 

# create names for S, A and distance matrix
reps <- rep(1:50, each = 100)
comps <- rep(1:50, times = 100)
names_50 <- str_c(reps, comps, sep = "_")

# S--------------------#
# list of 100 matrices
facs_resid_S_50 <- list()
for(i in 1:100) {
  facs_resid_S_50[[i]] <- facs_resid_ica_list_50[[i]]$S
} 
# put them all in one matrix
facs_resid_S_50_mat <- do.call(cbind, facs_resid_S_50)

# orient them consistently
col.medians <- apply(facs_resid_S_50_mat, 2, median)
sign <- ifelse(col.medians > 0, -1, 1)
facs_resid_S_50_mat <- t(t(facs_resid_S_50_mat)*sign)
colnames(facs_resid_S_50_mat) <- names_50
saveRDS(facs_resid_S_50_mat, "/data/cephfs/punim0613/icafast/output/facs_resid_S_50_mat.rds")

# get spearman distances for clustering
facs_resid_Spear_50 <- cor(facs_resid_S_50_mat, method = "spearman")
facs_resid_Spear_50 <- 1 - abs(facs_resid_Spear_50)
colnames(facs_resid_Spear_50) <- names_50
rownames(facs_resid_Spear_50) <- names_50
saveRDS(facs_resid_Spear_50, "/data/cephfs/punim0613/icafast/output/facs_resid_Spear_50.rds")

# PAM clustering
facs_resid_50_pam <- pam(x = facs_resid_Spear_50, k = 50, diss = TRUE, medoids = NULL, stand = F, cluster.only = F, 
                         do.swap = TRUE, keep.diss = F, keep.data = F, pamonce = 5, trace.lev = 1)
saveRDS(facs_resid_50_pam, "/data/cephfs/punim0613/icafast/output/facs_resid_50_pam.rds")

rm(facs_resid_S_50, facs_resid_S_50_mat, facs_resid_Spear_50)

# A--------------------#
# list of 100 matrices
facs_resid_A_50 <- list()
for(i in 1:100) {
  facs_resid_A_50[[i]] <- facs_resid_ica_list_50[[i]]$M
} 
# put them all in one matrix
facs_resid_A_50_mat <- do.call(cbind, facs_resid_A_50)

# orient them same way as S
facs_resid_A_50_mat <- t(t(facs_resid_A_50_mat)*sign)
colnames(facs_resid_A_50_mat) <- names_50
saveRDS(facs_resid_A_50_mat, "/data/cephfs/punim0613/icafast/output/facs_resid_A_50_mat.rds")

rm(facs_resid_A_50, facs_resid_A_50_mat)

# vaf------------------#
# list of 100 vectors
facs_resid_vaf_50 <- list()
for(i in 1:100) {
  facs_resid_vaf_50[[i]] <- facs_resid_ica_list_50[[i]]$vafs
} 
# put them all in one matrix
facs_resid_vaf_50_vec <- do.call(cbind, facs_resid_vaf_50)
saveRDS(facs_resid_vaf_50_vec, "/data/cephfs/punim0613/icafast/output/facs_resid_vaf_50.rds")

rm(facs_resid_vaf_50, facs_resid_vaf_50_vec)
rm(facs_resid_ica_list_50)

print("50 done ..."); Sys.time()


#---------------------------------------------R SESSION---#

# use pam clustering outputs to extract average and per-cluster silhouette stats for plotting, and medoids
# also derive reconstitution errors ON CLUSTER (based on one out of 100 runs or on medoids from pam clustering) - DONE ?
# plot and based on this choose optimal ncomp 
# for the optimal ncomp choose: a) most consitent ica run out of 100 and 
# b) the consensus set of components - maybe after removing low silhouette components and pam re-clustering (two approaches)
# for case b) order the S and A components according to the vafs or the "relative data power" method 


# Set the dataset to work on ...
folder <- c("/Users/adakic/projects/icatest/output/facs/",
            "/Users/adakic/projects/icatest/output/drop/")
folder <- folder[2] #--------------------------------------------CHANGE ACC.!
name <- c("facs_vst_",
          "facs_resid_",
          "facs_count_",
          "drop_vst_",
          "drop_resid_",
          "drop_count_")
name <- name[5] #------------------------------------------------CHANGE ACC.!


# Extract average, cluster-average and individual silhouette stats for all ncomps for plotting
# also medoids names and indices for extraction (computed on cluster)
sil.ave <- rep(0,10)
sil.clust.ave <- list()
sil.ind <- list()
medoid.names <- list()
medoid.idx <- list()


for(ncomp in c(50,100,150,200,250,300,350,400,450,500)) {
  connect <- paste0(folder, name, "PAM_", ncomp, ".rds") # when these are rsynced locally 
  pam_out <- readRDS(file = connect)
  
  i <- ncomp / 50
  
  sil.ave[i] <- pam_out$silinfo$avg.width
  sil.clust.ave[[i]] <- pam_out$silinfo$clus.avg.widths
  sil.ind[[i]] <- pam_out$silinfo$widths[,3]
  medoid.names[[i]] <- pam_out$medoids
  medoid.idx[[i]] <- pam_out$id.med
  
  print(paste(ncomp, "done..."))
}

# Import reconstitution errors (computed on cluster)
connect <- paste0(folder, name, "reconst_errors.rds")
reconst_errors <- readRDS(file = connect)

# Calculate summed vafs
vafs <- rep(0,10)
for(ncomp in c(50,100,150,200,250,300,350,400,450,500)) {
  connect <- paste0(folder, name, "vaf_", ncomp, ".rds")
  vaf <- readRDS(file = connect)
  vaf <- sum(vaf[,1])
  
  i <- ncomp / 50
  
  vafs[i] <- vaf
  
  print(paste(ncomp, "done..."))
}

# Import correlation scores (computed on cluster)
connect <- paste0(folder, name, "corr_scores.rds")
corr_scores <- readRDS(file = connect)


#####################   PLOT IT ALL   #####################   

nComps <- c(50,100,150,200,250,300,350,400,450,500)
silhStat <- sil.ave
fnStat <- reconst_errors$norm.F.error.m
fnStat.r <- reconst_errors$norm.F.error.r
vafs <- vafs

# StabilityVsError PLOT:
fname <- paste0(folder, "figures/", name, "StabilityVsError_10", ".png")
png(filename = fname, width = 600, height = 500, units = "px", pointsize = 13)
## add extra space to right margin of plot within frame
par(mar=c(5, 4, 4, 6) + 0.1)
## Plot first set of data and draw its axis
plot(nComps, silhStat, pch=16, axes=FALSE, ylim=c(0,1), xlab="", ylab="", 
     type="b",col="dodgerblue4", main="Stability vs. Reconstitution error")
mtext("Component stability  (Average silhouette width)", side=2, col="dodgerblue4", line=2.9)
axis(2, ylim=c(0,1), col="dodgerblue4", col.axis="dodgerblue4", las=1)  ## las=1 makes horizontal labels
box()
## Allow a second plot on the same graph
par(new=TRUE)
## Plot the second plot and put axis scale on right
plot(nComps, fnStat, pch=21,  xlab="", ylab="", ylim=c(0.3,0.65), 
     axes=FALSE, type="b", col="darkred")
## a little farther out (line=4) to make room for labels
mtext("Reconstitution error  (normalised)", side=4, col="darkred", line=3) 
axis(4, ylim=c(0.3,0.65), col="darkred",col.axis="darkred", las=1)
## Draw the time axis
axis(1, at = nComps)
mtext("Number of components",side=1,col="black",line=2.5)  
## Add Legend
legend("topright",legend=c("Component stability","Reconstitution error"),
       text.col=c("dodgerblue4","darkred"),pch=c(16,21),col=c("dodgerblue4","darkred"))
dev.off()
#-----------------------------------------------------------------------------#

# StabilityVsError_Dg PLOT:
fname <- paste0(folder, "figures/", name, "StabilityVsError_10_Dg", ".png")
png(filename = fname, width = 600, height = 500, units = "px", pointsize = 13)
## add extra space to right margin of plot within frame
par(mar=c(5, 4, 4, 6) + 0.1)
## Plot first set of data and draw its axis
plot(nComps, silhStat, pch=16, axes=FALSE, ylim=c(0,1), xlab="", ylab="", 
     type="b",col="dodgerblue4", main="Stability vs. Reconstitution error")
mtext("Component stability  (Average silhouette width)", side=2, col="dodgerblue4", line=2.9)
axis(2, ylim=c(0,1), col="dodgerblue4", col.axis="dodgerblue4", las=1)  ## las=1 makes horizontal labels
box()
## Allow a second plot on the same graph
par(new=TRUE)
## Plot the second plot and put axis scale on right
plot(nComps, fnStat.r, pch=21,  xlab="", ylab="", ylim=c(0.3,0.65), 
     axes=FALSE, type="b", col="darkred")
## a little farther out (line=4) to make room for labels
mtext("Reconstitution error  (normalised)", side=4, col="darkred", line=3) 
axis(4, ylim=c(0.3,0.65), col="darkred",col.axis="darkred", las=1)
## Draw the time axis
axis(1, at = nComps)
mtext("Number of components",side=1,col="black",line=2.5)  
## Add Legend
legend("topright",legend=c("Component stability","Reconstitution error"),
       text.col=c("dodgerblue4","darkred"),pch=c(16,21),col=c("dodgerblue4","darkred"))
dev.off()
#-----------------------------------------------------------------------------#

# StabilityVsVaf PLOT:
fname <- paste0(folder, "figures/", name, "StabilityVsVaf_10", ".png")
png(filename = fname, width = 600, height = 500, units = "px", pointsize = 13)
## add extra space to right margin of plot within frame
par(mar=c(5, 4, 4, 6) + 0.1)
## Plot first set of data and draw its axis
plot(nComps, silhStat, pch=16, axes=FALSE, ylim=c(0,1), xlab="", ylab="", 
     type="b",col="dodgerblue4", main="Stability vs. Variance accounted for")
axis(2, ylim=c(0,1),col="dodgerblue4", col.axis="dodgerblue4" , las=1)  ## las=1 makes horizontal labels
mtext("Component stability  (Average silhouette width)",side=2,col="dodgerblue4",line=2.9) # "Average silhouette width"
box()
## Allow a second plot on the same graph
par(new=TRUE)
## Plot the second plot and put axis scale on right
plot(nComps, vafs, pch=21,  xlab="", ylab="", ylim=c(0.5,0.75), 
     axes=FALSE, type="b", col="deeppink2")
## a little farther out (line=4) to make room for labels
mtext("Variance accounted for",side=4,col="deeppink2",line=3) 
axis(4, ylim=c(0.5,0.7), col="deeppink2",col.axis="deeppink2",las=1)
## Draw the time axis
axis(1, at = nComps)
mtext("Number of components",side=1,col="black",line=2.5)  
## Add Legend
legend("topright",legend=c("Component stability","Variance accounted for"),
       text.col=c("dodgerblue4","deeppink2"),pch=c(16,21),col=c("dodgerblue4","deeppink2"))
dev.off()
#-----------------------------------------------------------------------------#

# ScreeVsStability PLOT:
# here correlation stability - can add silhouette stability as well, or plot it separately as below
# png
fname <- paste0(folder, "figures/", name, "ScreeVsStability_10", ".png")
png(filename = fname, width = 850, height = 1400, units = "px")
par(mfrow = c(5,2))

for(i in c(50,100,150,200,250,300,350,400,450,500)) {
  
  comps <- 1:i
  j <- i/50
  
  connect <- paste0(folder, name, "vaf_", i, ".rds")
  vaf <- readRDS(file = connect) # i x 100
  vaf <- apply(vaf, 1, mean)
  print(paste("last", i, "vaf:", vaf[i]))
  
  connect <- paste0(folder, name, "null_vaf_", i, ".rds")
  vaf_null <- readRDS(file = connect) # i x 10
  vaf_null <- apply(vaf_null, 1, mean)
  print(paste("null crosses", i, "vaf at ncomp:" , min(which(vaf_null >= vaf))))
  
  connect <- paste0(folder, name, "corr_scores", ".rds")
  corr_scores <- readRDS(file = connect) # list of 10
  corr_scores <- corr_scores[[j]]
  print(paste("min", i, "corr score", min(corr_scores)))
  
  title <- paste(i, "components")
  par(mar = c(5, 5, 3, 5)) # default c(5.1, 4.1, 4.1, 2.1)
  plot(x=comps, y=vaf, type = "h", lwd = 2, xlab="", ylab="", axes=FALSE, main = title) # "Variance accounted for vs. Stability"
  lines(x=comps, y=vaf_null, col = "deeppink3", lwd = 2)
  axis(2, col="deeppink4", col.axis="deeppink4") # las=1
  mtext("Variance accounted for (data and null)", side=2, col="deeppink4", line=2.5) 
  axis(1)
  mtext("Component number", side=1, line = 2.5) 
  box()
  par(new=TRUE)
  plot(x=comps, y=corr_scores, type = "p", pch=21, col="dodgerblue3", xlab="", ylab="", axes=FALSE,ylim=c(0,1))
  axis(4, ylim=c(0,1), col="dodgerblue4",col.axis="dodgerblue4") # las=1
  mtext("Component stability (Correlation score)",side=4, col="dodgerblue4", line=2.5)
  
}
dev.off()


# pdf
fname <- paste0(folder, "figures/", name, "ScreeVsStability_10", ".pdf")
pdf(file = fname, width = 7, height = 13, pointsize = 8)
par(mfrow = c(5,2))

for(i in c(50,100,150,200,250,300,350,400,450,500)) {
  
  comps <- 1:i
  j <- i/50
  
  connect <- paste0(folder, name, "vaf_", i, ".rds")
  vaf <- readRDS(file = connect) # i x 100
  vaf <- apply(vaf, 1, mean)
  print(paste("last", i, "vaf:", vaf[i]))
  
  connect <- paste0(folder, name, "null_vaf_", i, ".rds")
  vaf_null <- readRDS(file = connect) # i x 10
  vaf_null <- apply(vaf_null, 1, mean)
  print(paste("null crosses", i, "vaf at ncomp:" , min(which(vaf_null >= vaf))))
  
  connect <- paste0(folder, name, "corr_scores", ".rds")
  corr_scores <- readRDS(file = connect) # list of 10
  corr_scores <- corr_scores[[j]]
  print(paste("min", i, "corr score", min(corr_scores)))
  
  title <- paste(i, "components") 
  par(mar = c(5, 5, 3, 5)) # default c(5.1, 4.1, 4.1, 2.1)
  plot(x=comps, y=vaf, type = "h", lwd = 2, xlab="", ylab="", axes=FALSE, main = title) # "Variance accounted for vs. Stability"
  lines(x=comps, y=vaf_null, col = "deeppink3", lwd = 2)
  axis(2, col="deeppink4", col.axis="deeppink4") # las=1
  mtext("Variance accounted for (data and null)", side=2, col="deeppink4", line=2.5) 
  axis(1)
  mtext("Component number", side=1, line = 2.5) 
  box()
  par(new=TRUE)
  plot(x=comps, y=corr_scores, type = "p", pch=21, col="dodgerblue3", xlab="", ylab="", axes=FALSE,ylim=c(0,1))
  axis(4, ylim=c(0,1), col="dodgerblue4",col.axis="dodgerblue4") # las=1
  mtext("Component stability (Correlation score)",side=4, col="dodgerblue4", line=2.5)
  
}
dev.off()


# MAYBE better incorporate info from these figures in the above figure ...

# analyse & plot cluster silh. stats across ncomps
#---------------------------------------------------------#
x.min <- map(.x = sil.clust.ave, .f = min) %>% unlist() %>% min()
x.min <- min(0, x.min)

fname <- paste0(folder, "figures/", name, "Clust_Silhuettes_10", ".png")
png(filename = fname, width = 700, height = 1200, units = "px", pointsize = 13)
par(mfrow = c(5,2))

for(i in c(50,100,150,200,250,300,350,400,450,500)) {
  
  j <- i/50
  
  title <- paste(i, "components")
  hist(sil.clust.ave[[j]], breaks = 50, col = "dodgerblue", xlim = c(x.min, 1), xlab = "Cluster Silhouette scores", main = title)
}
dev.off()



# analyse & plot individual silh. stats across ncomps
#---------------------------------------------------------#
x.min <- map(.x = sil.ind, .f = min) %>% unlist() %>% min()
x.min <- min(0, x.min)

fname <- paste0(folder, "figures/", name, "Ind_Silhuettes_10", ".png")
png(filename = fname, width = 700, height = 1200, units = "px", pointsize = 13)
par(mfrow = c(5,2))

for(i in c(50,100,150,200,250,300,350,400,450,500)) {
  
  j <- i/50
  
  title <- paste(i, "components")
  hist(sil.ind[[j]], breaks = 80, col = "forestgreen", xlim = c(x.min, 1), xlab = "Individual Silhouette scores", main = title)
}
dev.off()


#-----------------------------------------------------------------------------#
# Continue with 100-400 components ... (on the cluster and rsync them):
# - using medoids extract consensus S and A matrices (_c), ordered them according to vafs
# - note which consensus components have average cluster silh. stats below ~ 0.5 !
# - also extract non-consensus S and A matrices (_r)

#---ON CLUSTER (CONSENSUS_EXTRACTOR.R):

#-----------------------------------------------------------------------------#


S.100 <- readRDS(paste0(folder, name, "S_100_c.rds"))
S.100.Dg <- readRDS(paste0(folder, name, "S_100_r.rds"))

S.200 <- readRDS(paste0(folder, name, "S_200_c.rds"))
S.200.Dg <- readRDS(paste0(folder, name, "S_200_r.rds"))

S.300 <- readRDS(paste0(folder, name, "S_300_c.rds"))
S.300.Dg <- readRDS(paste0(folder, name, "S_300_r.rds"))



# To get entrez gene neames etc. for pathway analyses.. a handy way is to use scater's biomart connection
# but for this data needs to be in SingleCellExperiment format - so get one from BioC TabulaMurisData pack
library(ExperimentHub)
library(SingleCellExperiment)
library(TabulaMurisData)

eh <- ExperimentHub()
query(eh, "TabulaMurisData")

#-------------------------------------#
# TM.BioC.facs <- eh[["EH1618"]]
# TM.BioC.facs <- TabulaMurisSmartSeq2()
# tm.hivar <- TM.BioC.facs[rownames(S100), ]

#-------------------------------------#
TM.BioC.drop <- eh[["EH1617"]]
TM.BioC.drop <- TabulaMurisDroplet()
tm.hivar <- TM.BioC.drop[rownames(S.100), ]

# get "mgi_symbol", "ensembl_gene_id", "entrezgene" on 7.5k hv genes
library("biomaRt")
library("scater")
ensembl <- useEnsembl(biomart = "ensembl", dataset = 'mmusculus_gene_ensembl')
searchAttributes(mart = ensembl, pattern = "entrez")
tm.hivar <- getBMFeatureAnnos(tm.hivar, 
                              filters = c("mgi_symbol"), 
                              attributes = c("mgi_symbol", "ensembl_gene_id", "entrezgene_id", "entrezgene_accession"),  
                              biomart = "ENSEMBL_MART_ENSEMBL", 
                              dataset = "mmusculus_gene_ensembl", 
                              host = "www.ensembl.org")
#---------------------------------------------------------#
names(rowData(tm.hivar))

length(unique(rowData(tm.hivar)$entrezgene_id)) # 6749 (out of 7500)
#---------------------------------------------------------#

all.equal(rowData(tm.hivar)$Symbol, rownames(tm.hivar))
all.equal(rowData(tm.hivar)$Symbol, rownames(S.100)) # !!! same order

# Create ordered list of three types of gene names to attach to S matrix before defining modules
genenames <- tibble(Symbol = rowData(tm.hivar)$Symbol, 
                    ensembl_gene_id = rowData(tm.hivar)$ensembl_gene_id, 
                    entrezgene = rowData(tm.hivar)$entrezgene_id)

# attach
# 100
S.100.lst <- list()
for(i in 1:100) {
  S.100.lst[[i]] <- genenames
  S.100.lst[[i]]$S <- S.100[ ,i]
}
S.100.Dg.lst <- list()
for(i in 1:100) {
  S.100.Dg.lst[[i]] <- genenames
  S.100.Dg.lst[[i]]$S <- S.100.Dg[ ,i]
}
# 200
S.200.lst <- list()
for(i in 1:200) {
  S.200.lst[[i]] <- genenames
  S.200.lst[[i]]$S <- S.200[ ,i]
}
S.200.Dg.lst <- list()
for(i in 1:200) {
  S.200.Dg.lst[[i]] <- genenames
  S.200.Dg.lst[[i]]$S <- S.200.Dg[ ,i]
}
# 300
S.300.lst <- list()
for(i in 1:300) {
  S.300.lst[[i]] <- genenames
  S.300.lst[[i]]$S <- S.300[ ,i]
}
S.300.Dg.lst <- list()
for(i in 1:300) {
  S.300.Dg.lst[[i]] <- genenames
  S.300.Dg.lst[[i]]$S <- S.300.Dg[ ,i]
}

# some checks
sum(abs(S.100.lst[[6]]$S) >= 4)
sum(abs(S.100.Dg.lst[[6]]$S) >= 4)

sum(abs(S.200.lst[[6]]$S) >= 4)
sum(abs(S.200.Dg.lst[[6]]$S) >= 4)

sum(abs(S.300.lst[[6]]$S) >= 4)
sum(abs(S.300.Dg.lst[[6]]$S) >= 4)

# filter on 2, 3 and 4 sd
# 100
S.100.2sd.lst <- purrr::map(S.100.lst, .f = dplyr::filter, abs(S) >= 2)
S.100.3sd.lst <- purrr::map(S.100.lst, .f = dplyr::filter, abs(S) >= 3)
S.100.4sd.lst <- purrr::map(S.100.lst, .f = dplyr::filter, abs(S) >= 4)

S.100.Dg.2sd.lst <- purrr::map(S.100.Dg.lst, .f = dplyr::filter, abs(S) >= 2)
S.100.Dg.3sd.lst <- purrr::map(S.100.Dg.lst, .f = dplyr::filter, abs(S) >= 3)
S.100.Dg.4sd.lst <- purrr::map(S.100.Dg.lst, .f = dplyr::filter, abs(S) >= 4)

# 200
S.200.2sd.lst <- purrr::map(S.200.lst, .f = dplyr::filter, abs(S) >= 2)
S.200.3sd.lst <- purrr::map(S.200.lst, .f = dplyr::filter, abs(S) >= 3)
S.200.4sd.lst <- purrr::map(S.200.lst, .f = dplyr::filter, abs(S) >= 4)

S.200.Dg.2sd.lst <- purrr::map(S.200.Dg.lst, .f = dplyr::filter, abs(S) >= 2)
S.200.Dg.3sd.lst <- purrr::map(S.200.Dg.lst, .f = dplyr::filter, abs(S) >= 3)
S.200.Dg.4sd.lst <- purrr::map(S.200.Dg.lst, .f = dplyr::filter, abs(S) >= 4)

# 300
S.300.2sd.lst <- purrr::map(S.300.lst, .f = dplyr::filter, abs(S) >= 2)
S.300.3sd.lst <- purrr::map(S.300.lst, .f = dplyr::filter, abs(S) >= 3)
S.300.4sd.lst <- purrr::map(S.300.lst, .f = dplyr::filter, abs(S) >= 4)

S.300.Dg.2sd.lst <- purrr::map(S.300.Dg.lst, .f = dplyr::filter, abs(S) >= 2)
S.300.Dg.3sd.lst <- purrr::map(S.300.Dg.lst, .f = dplyr::filter, abs(S) >= 3)
S.300.Dg.4sd.lst <- purrr::map(S.300.Dg.lst, .f = dplyr::filter, abs(S) >= 4)

# some module stats:
# module size
range(sapply(S.100.2sd.lst, nrow))
range(sapply(S.100.3sd.lst, nrow))
range(sapply(S.100.4sd.lst, nrow))

range(sapply(S.200.2sd.lst, nrow)) # max = 359
range(sapply(S.200.3sd.lst, nrow))
range(sapply(S.200.4sd.lst, nrow))

range(sapply(S.300.2sd.lst, nrow))
range(sapply(S.300.3sd.lst, nrow))
range(sapply(S.300.4sd.lst, nrow)) # min = 1

# number of unique genes contributing to this (out of 7500)
uniq.100.2 <- purrr::map(S.100.2sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length() 
uniq.100.3 <- purrr::map(S.100.3sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length() 
uniq.100.4 <- purrr::map(S.100.4sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length() 

uniq.200.2 <- purrr::map(S.200.2sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length() 
uniq.200.3 <- purrr::map(S.200.3sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length() 
uniq.200.4 <- purrr::map(S.200.4sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length() 

uniq.300.2 <- purrr::map(S.300.2sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length() 
uniq.300.3 <- purrr::map(S.300.3sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length() 
uniq.300.4 <- purrr::map(S.300.4sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length() 

# average gene reuse factor
reus.100.2 <- round(sum(sapply(S.100.2sd.lst, nrow)) / purrr::map(S.100.2sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length(), 1)
reus.100.3 <- round(sum(sapply(S.100.3sd.lst, nrow)) / purrr::map(S.100.3sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length(), 1)
reus.100.4 <- round(sum(sapply(S.100.4sd.lst, nrow)) / purrr::map(S.100.4sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length(), 1) 

reus.200.2 <- round(sum(sapply(S.200.2sd.lst, nrow)) / purrr::map(S.200.2sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length(), 1)
reus.200.3 <- round(sum(sapply(S.200.3sd.lst, nrow)) / purrr::map(S.200.3sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length(), 1)
reus.200.4 <- round(sum(sapply(S.200.4sd.lst, nrow)) / purrr::map(S.200.4sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length(), 1) 

reus.300.2 <- round(sum(sapply(S.300.2sd.lst, nrow)) / purrr::map(S.300.2sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length(), 1)
reus.300.3 <- round(sum(sapply(S.300.3sd.lst, nrow)) / purrr::map(S.300.3sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length(), 1)
reus.300.4 <- round(sum(sapply(S.300.4sd.lst, nrow)) / purrr::map(S.300.4sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length(), 1) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# extent of overlap between consensus and "diagonal" ICA runs in terms of module gene composition (proportion of mudule genes shared)
symbol.100.3sd <- S.100.3sd.lst %>% purrr::map(.f = dplyr::select, Symbol) %>% flatten()
names(symbol.100.3sd) <- as.character(1:100)
symbol.100.Dg.3sd <- S.100.Dg.3sd.lst %>% purrr::map(.f = dplyr::select, Symbol) %>% flatten()
names(symbol.100.Dg.3sd) <- as.character(1:100)

olap1 <- list(); olap2 <- list()
for(i in 1:100) {
  olap1[[i]] <- purrr::map(symbol.100.Dg.3sd, .f = function(a,b) length(intersect(a,b)) / length(b), b = symbol.100.3sd[[i]]) %>% unlist() %>% max()
  olap1 <- unlist(olap1)
  
  olap2[[i]] <- purrr::map(symbol.100.Dg.3sd, .f = function(a,b) length(intersect(a,b)) / length(b), b = symbol.100.3sd[[i]]) %>% unlist() %>% sort(decreasing=T) %>% dplyr::nth(n=2)
  olap2 <- unlist(olap2)
}
min(olap1); max(olap2)

# compare this with overlap between different consensus modules
olap0 <- list()
for(i in 1:100) {
  olap0[[i]] <- purrr::map(symbol.100.3sd, .f = function(a,b) length(intersect(a,b)) / length(b), b = symbol.100.3sd[[i]]) %>% unlist() %>% sort(decreasing=T) %>% dplyr::nth(n=2)
  olap0 <- unlist(olap0)
}
max(olap0)

# a quick plot of it just for 100 comps & 3sd
plot(x=1:100, y=olap0, type = "h", lwd = 1.25, pch = 3, col = "slategrey", ylim = c(-0.125,1), xlab="module", ylab="Overlap (proportion of consensus module)", main = "Overlap between consensus and diagonal ICA run modules") # "Variance accounted for vs. Stability"
points(x=1:100, y=olap1, type = "p", pch = 16, col = "red")
points(x=1:100, y=olap2, type = "p", pch = 21, col = "red")
legend("bottom",legend=c("Best module overlap (consensus vs. diagonal)", "Second best module overlap (consensus vs. diagonal)", "Second best overlap ( between consensus modules)"), 
       col=c("red","red","slategrey"), 
       pch=c(16,21,3), lty=c(0,0,1), ncol = 1, inset=0.01, cex = 0.8)

##
# similarly, extent of overlap between consensus modules of facs and drop datasets in terms of gene composition (proportion of mudule genes shared)
facs.S.100 <- readRDS(paste0("/Users/adakic/projects/icatest/output/facs/", "facs_resid_", "S_100_c.rds"))

facs.S.100.lst <- list()
for(i in 1:100) {
  facs.S.100.lst[[i]] <- tibble::enframe(facs.S.100[ ,i], name = "Symbol", value = "S")
}

facs.S.100.3sd.lst <- purrr::map(facs.S.100.lst, .f = dplyr::filter, abs(S) >= 3)

facs.symbol.100.3sd <- facs.S.100.3sd.lst %>% purrr::map(.f = dplyr::select, Symbol) %>% flatten()
names(facs.symbol.100.3sd) <- as.character(1:100)

olap1 <- list(); olap2 <- list()
for(i in 1:100) {
  olap1[[i]] <- purrr::map(facs.symbol.100.3sd, .f = function(a,b) length(intersect(a,b)) / length(b), b = symbol.100.3sd[[i]]) %>% unlist() %>% max()
  olap1 <- unlist(olap1)
  
  olap2[[i]] <- purrr::map(facs.symbol.100.3sd, .f = function(a,b) length(intersect(a,b)) / length(b), b = symbol.100.3sd[[i]]) %>% unlist() %>% sort(decreasing=T) %>% dplyr::nth(n=2)
  olap2 <- unlist(olap2)
}
min(olap1); max(olap2)

plot(x=1:100, y=olap1, type = "p", pch = 16, col = "red", ylim = c(-0.1,1), xlab="module", ylab="Overlap (proportion of consensus 10X module)", main = "Overlap between SmartSeq2 and 10X consensus ICA modules")
points(x=1:100, y=olap2, type = "p", pch = 21, col = "red")
legend("bottom",legend=c("Best module overlap (SmartSeq2 vs. 10X)", "Second best module overlap (SmartSeq2 vs. 10X)"), 
       col=c("red","red"), 
       pch=c(16,21), ncol = 1, inset=0.01, cex = 0.9)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



########## Histograms of module sizes (across ncomp and sd parameters) 
# including inset indicating the number of unigue genes and average gene reuse facror [ uniq (reus x)]
fname <- paste0(folder, "figures/", name, "module_sizes.png")
png(filename = fname, width = 800, height = 500, units = "px")
par(mfrow=c(3,3)) # breaks = seq(0, 500, by = 20)

hist(sapply(S.100.2sd.lst, nrow), breaks = seq(0, 360, by = 10), col = "thistle" , xlab = "Module size", main = "100 modules @ 2 sd")
text(300, 15, paste0(uniq.100.2, "  (", reus.100.2, " x", ")"), cex = 1.5, col = "blue")
hist(sapply(S.100.3sd.lst, nrow), breaks = seq(0, 360, by = 10), col = "thistle" , xlab = "Module size", main = "100 modules @ 3 sd")
text(300, 15, paste0(uniq.100.3, "  (", reus.100.3, " x", ")"), cex = 1.5, col = "blue")
hist(sapply(S.100.4sd.lst, nrow), breaks = seq(0, 360, by = 10), col = "thistle" , xlab = "Module size", main = "100 modules @ 4 sd")
text(300, 15, paste0(uniq.100.4, "  (", reus.100.4, " x", ")"), cex = 1.5, col = "blue")
  
hist(sapply(S.200.2sd.lst, nrow), breaks = seq(0, 360, by = 10), col = "thistle" , xlab = "Module size", main = "200 modules @ 2 sd")
text(300, 15, paste0(uniq.200.2, "  (", reus.200.2, " x", ")"), cex = 1.5, col = "blue")
hist(sapply(S.200.3sd.lst, nrow), breaks = seq(0, 360, by = 10), col = "thistle" , xlab = "Module size", main = "200 modules @ 3 sd")
text(300, 15, paste0(uniq.200.3, "  (", reus.200.3, " x", ")"), cex = 1.5, col = "blue")
hist(sapply(S.200.4sd.lst, nrow), breaks = seq(0, 360, by = 10), col = "thistle" , xlab = "Module size", main = "200 modules @ 4 sd")
text(300, 15, paste0(uniq.200.4, "  (", reus.200.4, " x", ")"), cex = 1.5, col = "blue")

hist(sapply(S.300.2sd.lst, nrow), breaks = seq(0, 360, by = 10), col = "thistle" , xlab = "Module size", main = "300 modules @ 2 sd")
text(300, 15, paste0(uniq.300.2, "  (", reus.300.2, " x", ")"), cex = 1.5, col = "blue")
hist(sapply(S.300.3sd.lst, nrow), breaks = seq(0, 360, by = 10), col = "thistle" , xlab = "Module size", main = "300 modules @ 3 sd")
text(300, 15, paste0(uniq.300.3, "  (", reus.300.3, " x", ")"), cex = 1.5, col = "blue")
hist(sapply(S.300.4sd.lst, nrow), breaks = seq(0, 360, by = 10), col = "thistle" , xlab = "Module size", main = "300 modules @ 4 sd")
text(300, 15, paste0(uniq.300.4, "  (", reus.300.4, " x", ")"), cex = 1.5, col = "blue")

dev.off()


########## Barplot of re-use of genes by modules (across ncomp and sd parameters) 
# x = in how many modules genes contribute (y = how many genes like that)
fname <- paste0(folder, "figures/", name, "module_gene_reuse.png")
png(filename = fname, width = 900, height = 800, units = "px")
par(mfrow=c(3,3))

purrr::map(S.100.2sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% table() %>% sort(decreasing = T) %>% table() %>% 
  barplot(xlab = "number of modules", ylab = "number of unique genes", main = "100 modules @ 2 sd", ylim = c(0,2000), xlim = c(1,80), col = "tomato")
text(20, 1500, paste0(uniq.100.2), cex = 1.5, col = "blue")
purrr::map(S.100.3sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% table() %>% sort(decreasing = T) %>% table() %>% 
  barplot(xlab = "number of modules", ylab = "number of unique genes", main = "100 modules @ 3 sd", ylim = c(0,2000), xlim = c(1,80), col = "tomato") 
text(20, 1500, paste0(uniq.100.3), cex = 1.5, col = "blue")
purrr::map(S.100.4sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% table() %>% sort(decreasing = T) %>% table() %>% 
  barplot(xlab = "number of modules", ylab = "number of unique genes", main = "100 modules @ 4 sd", ylim = c(0,2000), xlim = c(1,80), col = "tomato") 
text(20, 1500, paste0(uniq.100.4), cex = 1.5, col = "blue")

purrr::map(S.200.2sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% table() %>% sort(decreasing = T) %>% table() %>% 
  barplot(xlab = "number of modules", ylab = "number of unique genes", main = "200 modules @ 2 sd", ylim = c(0,2000), xlim = c(1,80), col = "tomato") 
text(20, 1500, paste0(uniq.200.2), cex = 1.5, col = "blue")
purrr::map(S.200.3sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% table() %>% sort(decreasing = T) %>% table() %>% 
  barplot(xlab = "number of modules", ylab = "number of unique genes", main = "200 modules @ 3 sd", ylim = c(0,2000), xlim = c(1,80), col = "tomato") 
text(20, 1500, paste0(uniq.200.3), cex = 1.5, col = "blue")
purrr::map(S.200.4sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% table() %>% sort(decreasing = T) %>% table() %>% 
  barplot(xlab = "number of modules", ylab = "number of unique genes", main = "200 modules @ 4 sd", ylim = c(0,2000), xlim = c(1,80), col = "tomato") 
text(20, 1500, paste0(uniq.200.4), cex = 1.5, col = "blue")

purrr::map(S.300.2sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% table() %>% sort(decreasing = T) %>% table() %>% 
  barplot(xlab = "number of modules", ylab = "number of unique genes", main = "300 modules @ 2 sd", ylim = c(0,2000), xlim = c(1,80), col = "tomato") 
text(20, 1500, paste0(uniq.300.2), cex = 1.5, col = "blue")
purrr::map(S.300.3sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% table() %>% sort(decreasing = T) %>% table() %>% 
  barplot(xlab = "number of modules", ylab = "number of unique genes", main = "300 modules @ 3 sd", ylim = c(0,2000), xlim = c(1,80), col = "tomato") 
text(20, 1500, paste0(uniq.300.3), cex = 1.5, col = "blue")
purrr::map(S.300.4sd.lst, .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% table() %>% sort(decreasing = T) %>% table() %>% 
  barplot(xlab = "number of modules", ylab = "number of unique genes", main = "300 modules @ 4 sd", ylim = c(0,2000), xlim = c(1,80), col = "tomato") 
text(20, 1500, paste0(uniq.300.4), cex = 1.5, col = "blue")

dev.off()


########## Histogram of module weights (across ncomp and sd parameters)
# i.e. what proportion of component's total weight is being captured by the module genes 

# e.g. how to calculate proportion of component's total "weight" i.e. variance captured by the module genes 
# sum(S.100[,1]^2) # 7500 - comonent's weigh is always 7500 = number of genes 
# sum(S.100.2sd.lst[[1]]$S^2)/7500 # 91%

fname <- paste0(folder, "figures/", name, "module_weight.png")
png(filename = fname, width = 800, height = 600, units = "px")
par(mfrow=c(3,3))

purrr::map(S.100.2sd.lst, .f = dplyr::select, S) %>% purrr::map(.f = function(x) x^2) %>% purrr::map(.f = sum) %>% unlist() %>% (function(x) x / 7500) %>% 
  hist(breaks = seq(0, 1, by = 0.025), col = "linen" , xlab = "Module weight (proportion of component total weight)", main = "100 modules @ 2 sd")
purrr::map(S.100.3sd.lst, .f = dplyr::select, S) %>% purrr::map(.f = function(x) x^2) %>% purrr::map(.f = sum) %>% unlist() %>% (function(x) x / 7500) %>% 
  hist(breaks = seq(0, 1, by = 0.025), col = "linen" , xlab = "Module weight (proportion of component total weight)", main = "100 modules @ 3 sd")
purrr::map(S.100.4sd.lst, .f = dplyr::select, S) %>% purrr::map(.f = function(x) x^2) %>% purrr::map(.f = sum) %>% unlist() %>% (function(x) x / 7500) %>% 
  hist(breaks = seq(0, 1, by = 0.025), col = "linen" , xlab = "Module weight (proportion of component total weight)", main = "100 modules @ 4 sd")

purrr::map(S.200.2sd.lst, .f = dplyr::select, S) %>% purrr::map(.f = function(x) x^2) %>% purrr::map(.f = sum) %>% unlist() %>% (function(x) x / 7500) %>% 
  hist(breaks = seq(0, 1, by = 0.025), col = "linen" , xlab = "Module weight (proportion of component total weight)", main = "200 modules @ 2 sd")
purrr::map(S.200.3sd.lst, .f = dplyr::select, S) %>% purrr::map(.f = function(x) x^2) %>% purrr::map(.f = sum) %>% unlist() %>% (function(x) x / 7500) %>% 
  hist(breaks = seq(0, 1, by = 0.025), col = "linen" , xlab = "Module weight (proportion of component total weight)", main = "200 modules @ 3 sd")
purrr::map(S.200.4sd.lst, .f = dplyr::select, S) %>% purrr::map(.f = function(x) x^2) %>% purrr::map(.f = sum) %>% unlist() %>% (function(x) x / 7500) %>% 
  hist(breaks = seq(0, 1, by = 0.025), col = "linen" , xlab = "Module weight (proportion of component total weight)", main = "200 modules @ 4 sd")

purrr::map(S.300.2sd.lst, .f = dplyr::select, S) %>% purrr::map(.f = function(x) x^2) %>% purrr::map(.f = sum) %>% unlist() %>% (function(x) x / 7500) %>% 
  hist(breaks = seq(0, 1, by = 0.025), col = "linen" , xlab = "Module weight (proportion of component total weight)", main = "300 modules @ 2 sd")
purrr::map(S.300.3sd.lst, .f = dplyr::select, S) %>% purrr::map(.f = function(x) x^2) %>% purrr::map(.f = sum) %>% unlist() %>% (function(x) x / 7500) %>% 
  hist(breaks = seq(0, 1, by = 0.025), col = "linen" , xlab = "Module weight (proportion of component total weight)", main = "300 modules @ 3 sd")
purrr::map(S.300.4sd.lst, .f = dplyr::select, S) %>% purrr::map(.f = function(x) x^2) %>% purrr::map(.f = sum) %>% unlist() %>% (function(x) x / 7500) %>% 
  hist(breaks = seq(0, 1, by = 0.025), col = "linen" , xlab = "Module weight (proportion of component total weight)", main = "300 modules @ 4 sd")

dev.off()

## entrezgenes ##
...

## Gene Symbols ##
...

## Diagnostic on medoid vs. diagonal modules - HOW MUCH THEY OVERLAP IN TERMS OF MODULE GENE CONTENT ##   !!!
#-         Essentially you need to also do all the above plots with diagonal modules version         -#
...

## Doing the gene set analyses ##   !!! LOTS WORK HERE !!!
library("ReactomePA")
library("clusterProfiler")
...


# Focus on the 100 3sd: S.100.3sd.lst
  
## entrezgenes ##
entrez.3sd <- S.100.3sd.lst %>% purrr::map(.f = dplyr::select, entrezgene) %>% flatten()
names(entrez.3sd) <- as.character(1:100)

s.3sd <- S.100.3sd.lst %>% purrr::map(.f = dplyr::select, S) %>% flatten()
names(s.3sd) <- as.character(1:100)

# Trying REactomePA:
enrichPathway() # takes just a geneList, can limit universe; readable = T gives gene symbols / names listed !!!
gsePathway() # takes order ranked geneList
viewPathway("M phase", organism = "mouse")

# Some selected modules at Myeloid level: 
# Monocytes: +: 89! 90 91 38, -: 45! 60 69 80 36
# Macrophages: +: 12! 69 52, -: 95! 46 59 75

entrez.myelo.subset <- entrez.3sd[c(89, 90, 91, 38, 45, 60, 69, 80, 36, 12, 52, 95, 46, 59, 75)]
s.myelo.subset <- s.3sd[c(89, 90, 91, 38, 45, 60, 69, 80, 36, 12, 52, 95, 46, 59, 75)]

###
# compareCluster(geneClusters, fun = "enrichGO", data = "", 
#                ...) # "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway"
reacto.myelo <- compareCluster(entrez.myelo.subset, fun = "enrichPathway", organism = "mouse", readable = T)
dim(reacto.myelo) # 150  10
###
# dotplot(object, x = "GeneRatio",
#         color = "p.adjust", showCategory = 10, size = NULL, split = NULL,
#         font.size = 12, title = "", ...)
dotplot(reacto.myelo, showCategory = 4, font.size = 10, title = "Main Myeloid Cell Identity Modules")


### ENRICHMENTS ###

###
# enrichPathway(gene, organism = "human", pvalueCutoff = 0.05,
#               pAdjustMethod = "BH", qvalueCutoff = 0.2, universe, minGSSize = 10,
#               maxGSSize = 500, readable = FALSE)
enrich.myelo.89 <- enrichPathway(entrez.myelo.subset$`89`, organism = "mouse", pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE)
enrich.myelo.45 <- enrichPathway(entrez.myelo.subset$`45`, organism = "mouse", pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE)
enrich.myelo.12 <- enrichPathway(entrez.myelo.subset$`12`, organism = "mouse", pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE)
enrich.myelo.38 <- enrichPathway(entrez.myelo.subset$`38`, organism = "mouse", pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE)
enrich.myelo.69 <- enrichPathway(entrez.myelo.subset$`69`, organism = "mouse", pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE)
enrich.myelo.95 <- enrichPathway(entrez.myelo.subset$`95`, organism = "mouse", pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE)



### PLOTS ###

###
# emapplot(x, showCategory = 30, color = "p.adjust", layout = "kk",
#          ...)
emapplot(enrich.myelo.89)
emapplot(enrich.myelo.45)
emapplot(enrich.myelo.12)
emapplot(enrich.myelo.38)
emapplot(enrich.myelo.69)
emapplot(enrich.myelo.95)


###
# cnetplot(x, showCategory = 5, foldChange = NULL, layout = "kk", colorEdge = FALSE, circular = FALSE, node_label = TRUE, 
#          ...)
names(s.myelo.subset$`89`) <- entrez.myelo.subset$`89` # need to give it named vector...
cnetplot(enrich.myelo.89, showCategory = 6, foldChange = s.myelo.subset$`89`, colorEdge = T) + theme(legend.position="none")
names(s.myelo.subset$`45`) <- entrez.myelo.subset$`45` # need to give it named vector...
cnetplot(enrich.myelo.45, showCategory = 6, foldChange = s.myelo.subset$`45`, colorEdge = T) + theme(legend.position="none")
names(s.myelo.subset$`12`) <- entrez.myelo.subset$`12` # need to give it named vector...
cnetplot(enrich.myelo.12, showCategory = 6, foldChange = s.myelo.subset$`12`, colorEdge = T) + theme(legend.position="none")
names(s.myelo.subset$`38`) <- entrez.myelo.subset$`38` # need to give it named vector...
cnetplot(enrich.myelo.38, showCategory = 7, foldChange = s.myelo.subset$`38`, colorEdge = T) + theme(legend.position="none")
names(s.myelo.subset$`69`) <- entrez.myelo.subset$`69` # need to give it named vector...
cnetplot(enrich.myelo.69, showCategory = 7, foldChange = s.myelo.subset$`69`, colorEdge = T) + theme(legend.position="none")
names(s.myelo.subset$`95`) <- entrez.myelo.subset$`95` # need to give it named vector...
cnetplot(enrich.myelo.95, showCategory = 7, foldChange = s.myelo.subset$`95`, colorEdge = T) + theme(legend.position="none")

###
# heatplot(x, showCategory = 30, foldChange = NULL)
heatplot(enrich.myelo.89, showCategory = 10, foldChange = s.myelo.subset$`89`)
heatplot(enrich.myelo.45, showCategory = 10, foldChange = s.myelo.subset$`45`)
heatplot(enrich.myelo.12, showCategory = 10, foldChange = s.myelo.subset$`12`)
heatplot(enrich.myelo.38, showCategory = 10, foldChange = s.myelo.subset$`38`)
heatplot(enrich.myelo.69, showCategory = 10, foldChange = s.myelo.subset$`69`)
heatplot(enrich.myelo.95, showCategory = 10, foldChange = s.myelo.subset$`95`)



library("enrichplot")
library("DOSE")
###
# gseNCG(geneList, exponent = 1, nPerm = 1000, minGSSize = 10,
#        maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH",
#        verbose = TRUE, seed = FALSE, by = "fgsea")
gsea.myelo.89 <- gseNCG(sort(s.myelo.subset$`89`, decreasing = T), exponent = 1, nPerm = 1000)

gsea.myelo.89.abs <- gseNCG(sort(abs(s.myelo.subset$`89`), decreasing = T), exponent = 1, nPerm = 1000)

# maybe this gsea works on human names only ...?

###
# gseaplot(x, geneSetID, by = "all", title = "", color = "black", color.line = "green", color.vline = "#FA5860",
#          ...)


#-----------------------------------------------------------------------------#
# To explore meaning of + and - in S matrix, A matrix, Lasso coefficents matrix and relate it to gene expression matrix
# focus on say modules 89, 45, 12, 95 / 75, 51, 88

#- select those components (columns) from a.100 mat - summarise by cell type (all or 7 myelo and others) / [or pick 50-100 representatives from each 7 cell types] - heatplot
#- select 3sd genes from those components (columns) from S.100 mat - heatplot
#- select the same genes and summarise their expression by cell type (all or 7 myelo and others) / [or pick 50-100 representatives from each 7 cell types] - heatplot
# (but for gene exp you need to use log values ??? those "residual" values are good too as they are "relative to the population average")


subset894512 <- S.100.3sd.lst[c(89,45,12)] %>% purrr::map(.f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unname() %>% unique()
#
names.100.3sd <- S.100.3sd.lst %>% purrr::map(.f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unname() %>% unique() # 3994

exp.subset894512 <- exp.df[, c("use_celltype", subset894512)]

exp.summ.subset894512 <- exp.subset894512 %>% dplyr::group_by(use_celltype) %>% dplyr::summarise_all(.funs = mean)

tst <- exp.summ.subset894512 %>% as.data.frame() %>% column_to_rownames(var = "use_celltype") %>% as.matrix() %>% t()
library("limma")
coolmap(x = tst, cluster.by="de pattern", col=NULL, linkage.row="ward", linkage.col="ward", show.dendrogram="both")
coolmap(x = tst, cluster.by="expression level", col=NULL, linkage.row="ward", linkage.col="ward", show.dendrogram="both")


### 89 alone
# Exp per celltype
subset89 <- S.100.3sd.lst[c(89)] %>% purrr::map(.f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unname() %>% unique()
exp.subset89 <- exp.df[, c("use_celltype", subset89)] 
exp.summ.subset89 <- exp.subset89 %>% dplyr::group_by(use_celltype) %>% dplyr::summarise_all(.funs = mean)

tst <- exp.summ.subset89 %>% as.data.frame() %>% column_to_rownames(var = "use_celltype") %>% as.matrix() %>% t()
library("limma")
coolmap(x = tst, cluster.by="de pattern", col=NULL, linkage.row="average", linkage.col="average", show.dendrogram="both")
coolmap(x = tst, cluster.by="expression level", col=NULL, linkage.row="average", linkage.col="average", show.dendrogram="both")

# Exp
...
# better expore this in seurat plots ... !!!

# S
subset89.S <- S.100.3sd.lst[c(89)] %>% purrr::map(.f = dplyr::select, Symbol, S) %>% as.data.frame() %>% column_to_rownames(var = "Symbol") %>% as.matrix()
# cannot coolmap with just one column
subset89.S <- subset89.S[order(subset89.S[,1]), ]

Heatmap(subset89.S, show_heatmap_legend=T, name = "S.89")

# A per celtype 
subset89.A <- a.100[, c("use_celltype", 89)] 
subset89.A <- subset89.A %>% dplyr::group_by(use_celltype) %>% dplyr::summarise_all(.funs = median) %>% as.data.frame() %>% column_to_rownames(var = "use_celltype") %>% as.matrix() # better use median

Heatmap(subset89.A, show_heatmap_legend=T, name = "S.89")

# A 
...






# Some selected modules at Tissue-Macrophage level: 
# bm 75+!, spl 71+ 22+ 58-, lung 40+, kidney 6+! 54+, mammary 51+! 21+ 88-! 25- 99-

entrez.mac.subset <- entrez.3sd[c(75, 71, 22, 58, 40, 6, 54, 51, 21, 88, 25, 99)]














#-----A matrix focus - LASSO logistic regression (bi- and multi-nomial; global and tissue-specific)
###################################################################################################


# prepare data first
A.100 <- readRDS(paste0(folder, name, "A_100_c.rds"))
A.100.Dg <- readRDS(paste0(folder, name, "A_100_r.rds"))

A.200 <- readRDS(paste0(folder, name, "A_200_c.rds"))
A.200.Dg <- readRDS(paste0(folder, name, "A_200_r.rds"))

A.300 <- readRDS(paste0(folder, name, "A_300_c.rds"))
A.300.Dg <- readRDS(paste0(folder, name, "A_300_r.rds"))


# use appropriate Seurat objects to attach all potentially relevant covariates (per-cell annotation / metadata) to the A matrix 
# NOTE: here only important to pick the right dataset (facs or drop), normalisation (vst or sct) does not affect cell annotations 

dataset <- c("facs.", "drop.")
dataset <- dataset[2] #------------------------------------------------CHANGE ACC.!

transf <- c("exp.vst.7500hvg.mat.rds" , "exp.sct.7500hvg.mat.RESID.rds")
transf <- transf[2] #--------------------------------------------------CHANGE ACC.!

exp <- readRDS(paste0("/Users/adakic/projects/icatest/data/", dataset, "exp.sct.707.rds"))
exp.mat <- readRDS(paste0("/Users/adakic/projects/icatest/data/", dataset, transf)) 
exp.df <- as.data.frame(t(exp.mat)) %>% rownames_to_column(var = "cellnames")   # (need it later)!

# Make data frames first
A.100.df <- as.data.frame(A.100)
A.100.Dg.df <- as.data.frame(A.100.Dg)

A.200.df <- as.data.frame(A.200)
A.200.Dg.df <- as.data.frame(A.200.Dg)

A.300.df <- as.data.frame(A.300)
A.300.Dg.df <- as.data.frame(A.300.Dg)

# Before attaching annotations / metadata, check cellnames as id
all.equal(colnames(exp.mat), colnames(exp)) # TRUE

# Tuck-in cellnames as id (for matching below)
A.100.df$cellnames <- colnames(exp.mat)
A.100.Dg.df$cellnames <- colnames(exp.mat)

A.200.df$cellnames <- colnames(exp.mat)
A.200.Dg.df$cellnames <- colnames(exp.mat)

A.300.df$cellnames <- colnames(exp.mat)
A.300.Dg.df$cellnames <- colnames(exp.mat)


# NOTE: DATASET-SPECIFIC TASK (facs vs. drop) 
# Two versions of the below code segment depending on "facs" or "drop" - then goes back to common code:

# Prep df for lasso by removing key ambiguous cell annotations, fusing inconsistent annotation types etc.




######################## 
# FACS prep df for lasso: --- THIS FACS SECTION NEEDS TO BE DOUBLE CHECKED AND MAKE SURE ALL WORKS WELL AS IN DROP SECTION BELOW ---
if(dataset == "facs.") {
  
  print(paste("Doing", dataset)); cat("\n")
  
  # Attach full annotations to A.100.df first - then transfer them to other data frames
  
  # fix cell annotation strings
  cell_ontology_class <- str_replace_all(exp@meta.data$cell_ontology_class, " ", "_")
  A.100.df$celltype <- cell_ontology_class # do not facor it so you can change according to your needs
  A.100.df$tissue <- as.factor(exp@meta.data$tissue)
  A.100.df$mouse <- as.factor(exp@meta.data$mouse.id)
  A.100.df$sex <- as.factor(exp@meta.data$mouse.sex)
  A.100.df$plate <- as.factor(exp@meta.data$plate.barcode)
  A.100.df$nCount <- exp@meta.data$nCount_RNA
  A.100.df$nFeature <- exp@meta.data$nFeature_RNA
  A.100.df$nCount_SCT <- exp@meta.data$nCount_SCT
  A.100.df$nFeature_SCT <- exp@meta.data$nFeature_SCT
  
  print(unique(A.100.df$celltype) %>% sort()); cat("\n")
  
  # reorganise and delete ambigous cells annotations - but keep a version of intact annoations here
  a.100 <- as_tibble(A.100.df) %>% 
    dplyr::select("cellnames", "celltype", "tissue", "1":"100", "mouse", "sex", "plate", "nCount", "nFeature", "nCount_SCT", "nFeature_SCT") %>% 
    dplyr::mutate(celltype = if_else(celltype == "", "unknown", celltype)) %>% 
    dplyr::filter(!(celltype %in% c("unknown", "blood_cell", "leukocyte", "myeloid_cell", "professional_antigen_presenting_cell", "granulocytopoietic_cell")))
  
  # make a universal keep vector for all A matrices
  keep <- a.100$cellnames
  
  # organise and/or change annotations for several different downstream purposes
  # "use_celltype" will be currently-in-use cell annotation
  a.100 <- a.100 %>% 
    dplyr::mutate(use_celltype = if_else(celltype == "classical_monocyte", "monocyte", celltype)) %>% 
    dplyr::mutate(use_celltype = if_else(celltype %in% c("macrophage", "monocyte", "Kupffer_cell", "microglial_cell", "granulocyte"), celltype, "other"))
  
  a.100 %>% dplyr::filter(celltype =="classical_monocyte") %>% dplyr::select(tissue, celltype, use_celltype)
  
  # "use_celltype_tissue" will be currently-in-use cell-tissue annotation
  a.100 <- unite(a.100, use_celltype, tissue, col = "use_celltype_tissue", sep = "_", remove = F) %>% 
    dplyr::mutate(use_celltype_tissue = if_else(str_detect(use_celltype_tissue, "other_"), "other", use_celltype_tissue)) %>% 
    dplyr::mutate(use_celltype = factor(use_celltype, levels = c("granulocyte", "monocyte", "macrophage", "Kupffer_cell", "microglial_cell", "other"))) %>% 
    dplyr::mutate(use_celltype_tissue = factor(use_celltype_tissue, levels = c("granulocyte_Marrow", "monocyte_Marrow", "macrophage_Marrow", "macrophage_Spleen", "macrophage_Limb_Muscle", "macrophage_Diaphragm", "macrophage_Kidney", "macrophage_Brain_Myeloid", "Kupffer_cell_Liver", "microglial_cell_Brain_Myeloid", "other"))) %>% 
    dplyr::select("cellnames", "celltype", "tissue", "use_celltype", "use_celltype_tissue", "1":"100", "mouse", "sex", "plate", "nCount", "nFeature", "nCount_SCT", "nFeature_SCT")
  
  print(unique(a.100$use_celltype) %>% sort()); cat("\n")
  print(unique(a.100$use_celltype_tissue) %>% sort()); cat("\n")
  
  # Get pure annotations / metadata
  a.annos <- dplyr::select(a.100, "cellnames", "celltype", "tissue", "use_celltype", "use_celltype_tissue", "mouse", "sex", "plate", "nCount", "nFeature", "nCount_SCT", "nFeature_SCT")
  
  # Transfer the same annotations / metadata to all other A data frames 100-300 (trimmed to "keeps")
  # [revise clever usages of "assign", "get" (work in reverse) and "eval" in loops]
  for(i in c(200, 300)){
    i.str <- as.character(i)
    a <- as_tibble(get(paste0("A.", i, ".df"))) %>% 
      dplyr::filter(cellnames %in% keep) %>% 
      left_join(y = a.annos, by = "cellnames") %>% 
      dplyr::select("cellnames", "celltype", "tissue", "use_celltype", "use_celltype_tissue", "1":i.str, "mouse", "sex", "plate", "nCount", "nFeature", "nCount_SCT", "nFeature_SCT")
    
    assign(paste0("a.", i), a)
    rm(a)
  }
  # same for Diagonal version
  for(i in c(100, 200, 300)){
    i.str <- as.character(i)
    a <- as_tibble(get(paste0("A.", i, ".Dg.df"))) %>% 
      dplyr::filter(cellnames %in% keep) %>% 
      left_join(y = a.annos, by = "cellnames") %>% 
      dplyr::select("cellnames", "celltype", "tissue", "use_celltype", "use_celltype_tissue", "1":i.str, "mouse", "sex", "plate", "nCount", "nFeature", "nCount_SCT", "nFeature_SCT")
    
    assign(paste0("a.", i, ".d"), a)
    rm(a)
  }
  
  # remove original A mat. and df.
  for (i in c(100, 200, 300)) {
    rm(list = c(paste0("A.", i), paste0("A.", i, ".df"), paste0("A.", i, ".Dg"), paste0("A.", i, ".Dg.df")))
  }
  
  # Keep the same cells and Transfer the same annotations/metadata as above to Gene expression Data matrix
  hvg <- colnames(exp.df)[2:length(exp.df)]
  exp.df <- as_tibble(exp.df) %>% 
    dplyr::filter(cellnames %in% keep) %>% 
    left_join(y = a.annos, by = "cellnames") %>% 
    dplyr::select("cellnames", "celltype", "tissue", "use_celltype", "use_celltype_tissue", hvg, "mouse", "sex", "nCount", "nFeature", "nCount_SCT", "nFeature_SCT")
  
}


########################
# DROP prep df for lasso:
if(dataset == "drop.") {
  
  print(paste("Doing", dataset)); cat("\n")
  
  # Attach full annotations to A.100.df first - then transfer them to other data frames 
  
  # fix cell annotation strings
  cell_ontology_class <- str_replace_all(exp@meta.data$cell_ontology_class, " ", "_")
  A.100.df$celltype <- cell_ontology_class # do not facor it so you can change according to your needs
  A.100.df$tissue <- as.factor(exp@meta.data$tissue)
  A.100.df$mouse <- as.factor(exp@meta.data$mouse.id)
  A.100.df$sex <- as.factor(exp@meta.data$mouse.sex)
  #A.100.df$plate <- as.factor(exp@meta.data$plate.barcode) #--------------no plates in DROP dataset !!!
  A.100.df$nCount <- exp@meta.data$nCount_RNA
  A.100.df$nFeature <- exp@meta.data$nFeature_RNA
  A.100.df$nCount_SCT <- exp@meta.data$nCount_SCT
  A.100.df$nFeature_SCT <- exp@meta.data$nFeature_SCT
  
  print(unique(A.100.df$celltype) %>% sort()); cat("\n")
  
  # reorganise and delete ambigous cells annotations - but keep a version of intact annoations here
  a.100 <- as_tibble(A.100.df) %>% 
    dplyr::select("cellnames", "celltype", "tissue", "1":"100", "mouse", "sex", "nCount", "nFeature", "nCount_SCT", "nFeature_SCT") %>% 
    dplyr::mutate(celltype = if_else(celltype == "", "unknown", celltype)) %>% 
    dplyr::filter(!(celltype %in% c("unknown", "blood_cell", "leukocyte", "myeloid_cell", "granulocytopoietic_cell")))
  
  print(unique(a.100$celltype) %>% sort()); cat("\n")
  
  # make a universal keep vector for all A matrices
  keep <- a.100$cellnames
  
  # organise and/or change annotations for several different downstream purposes
  # "use_celltype" will be currently-in-use cell annotation
  a.100 <- dplyr::mutate(a.100, use_celltype = celltype)
  a.100 <- dplyr::mutate(a.100, use_celltype = case_when(use_celltype == "classical_monocyte" ~ "monocyte", 
                                                         use_celltype == "non-classical_monocyte" ~ "monocyte", 
                                                         use_celltype == "alveolar_macrophage" ~ "macrophage", 
                                                         TRUE ~ as.character(use_celltype))) %>% 
    dplyr::mutate(use_celltype = if_else(use_celltype %in% c("macrophage", "monocyte", "dendritic_cell", "promonocyte", "granulocyte", "mast_cell"), use_celltype, "other"))
  
  # "use_celltype_tissue" will be currently-in-use cell-tissue annotation
  a.100 <- unite(a.100, use_celltype, tissue, col = "use_celltype_tissue", sep = "_", remove = F) %>% 
    dplyr::mutate(use_celltype_tissue = if_else(str_detect(use_celltype_tissue, "other_"), "other", use_celltype_tissue))
  
  # some tests
  print(unique(a.100$use_celltype)); cat("\n"); print(unique(a.100$use_celltype_tissue))
  
  a.100 <- dplyr::mutate(a.100, use_celltype = factor(use_celltype, levels = c("granulocyte", "promonocyte", "monocyte", "macrophage", "dendritic_cell", "mast_cell", "other")))
  a.100 <- dplyr::mutate(a.100, use_celltype_tissue = factor(use_celltype_tissue, levels = c("granulocyte_Marrow", "promonocyte_Marrow", "monocyte_Marrow", "monocyte_Lung", 
                                                                                             "macrophage_Marrow", "macrophage_Spleen", "macrophage_Lung", "macrophage_Kidney", "macrophage_Mammary_Gland", "macrophage_Limb_Muscle",  
                                                                                             "dendritic_cell_Spleen", "mast_cell_Lung", "other")))
  a.100 <- dplyr::select(a.100, "cellnames", "celltype", "tissue", "use_celltype", "use_celltype_tissue", "1":"100", "mouse", "sex", "nCount", "nFeature", "nCount_SCT", "nFeature_SCT")
  # "promonocyte" are left as they are.. consider blocking them, same as "granulocytopoietic_cell"
  
  # some tests
  print(unique(a.100$use_celltype)); cat("\n"); print(unique(a.100$use_celltype_tissue))
  
  a.100 %>% dplyr::filter(is.na(use_celltype_tissue)) %>% dplyr::select("celltype", "tissue", "use_celltype", "use_celltype_tissue") %>% print(n=500)
  
  a.100 %>% dplyr::filter(celltype =="non-classical_monocyte") %>% dplyr::select(tissue, celltype, use_celltype, use_celltype_tissue)
  a.100 %>% dplyr::filter(celltype =="macrophage") %>% dplyr::select(tissue, celltype, use_celltype, use_celltype_tissue)
  a.100 %>% dplyr::filter(celltype =="fibroblast") %>% dplyr::select(tissue, celltype, use_celltype, use_celltype_tissue)
  
  # Get pure annotations / metadata for keep cells
  a.annos <- dplyr::select(a.100, "cellnames", "celltype", "tissue", "use_celltype", "use_celltype_tissue", "mouse", "sex", "nCount", "nFeature", "nCount_SCT", "nFeature_SCT")
  fname <- paste0(folder, "fits/", name, "a.annos", ".rds")
  saveRDS(a.annos, file = fname)
  
  # Transfer the same annotations / metadata to all other A data frames 100-300 (trimmed to "keeps")
  # [revise clever usages of "assign", "get" (work in reverse) and "eval" in loops]
  for(i in c(200, 300)){
    i.str <- as.character(i)
    a <- as_tibble(get(paste0("A.", i, ".df"))) %>% 
      dplyr::filter(cellnames %in% keep) %>% 
      left_join(y = a.annos, by = "cellnames") %>% 
      dplyr::select("cellnames", "celltype", "tissue", "use_celltype", "use_celltype_tissue", "1":i.str, "mouse", "sex", "nCount", "nFeature", "nCount_SCT", "nFeature_SCT")
    
    assign(paste0("a.", i), a)
    rm(a)
  }
  # same for Diagonal version
  for(i in c(100, 200, 300)){
    i.str <- as.character(i)
    a <- as_tibble(get(paste0("A.", i, ".Dg.df"))) %>% 
      dplyr::filter(cellnames %in% keep) %>% 
      left_join(y = a.annos, by = "cellnames") %>% 
      dplyr::select("cellnames", "celltype", "tissue", "use_celltype", "use_celltype_tissue", "1":i.str, "mouse", "sex", "nCount", "nFeature", "nCount_SCT", "nFeature_SCT")
    
    assign(paste0("a.", i, ".d"), a)
    rm(a)
  }
  
  # Keep the same cells and Transfer the same annotations/metadata as above to Gene expression Data matrix
  hvg <- colnames(exp.df)[2:length(exp.df)]
  exp.df <- as_tibble(exp.df) %>% 
    dplyr::filter(cellnames %in% keep) %>% 
    left_join(y = a.annos, by = "cellnames") %>% 
    dplyr::select("cellnames", "celltype", "tissue", "use_celltype", "use_celltype_tissue", hvg, "mouse", "sex", "nCount", "nFeature", "nCount_SCT", "nFeature_SCT")
  
  # some tests
  print(unique(exp.df$use_celltype)); cat("\n"); print(unique(exp.df$use_celltype_tissue))
  
  # remove original A mat. and df.
  for (i in c(100, 200, 300)) {
    rm(list = c(paste0("A.", i), paste0("A.", i, ".df"), paste0("A.", i, ".Dg"), paste0("A.", i, ".Dg.df")))
  }
  
}

rm(exp, exp.mat)

# for (i in c(100, 200, 300)) {
#   rm(list = c(paste0("a.", i), paste0("a.", i, ".d")))
# }



#---------Split to Training and Testing datasets----------#
library(caret)

# These cell indices are same for Gen.Exp. Df and all A matrices
# Best to make them based on use_celltype_tissue variable (to keep celltype-tissue representation)!
set.seed(707)
inTrain.tissue <- createDataPartition(y = exp.df$use_celltype_tissue, p = .85, list = FALSE) # SHOULD USE THIS ONE IN ALL CASSES !!!
fname <- paste0(folder, "fits/", name, "inTrain.tissue", ".rds")
saveRDS(inTrain.tissue, file = fname)
#---------------------------------------------------------#


# At this stage also good to save a.100-300 and exp.df:
fname <- paste0(folder, "fits/", name, "exp.df", ".rds")
saveRDS(exp.df, file = fname)

for (ncomp in c(100, 200, 300)) {
  
  a <- get(paste0("a.", ncomp))
  
  fname <- paste0(folder, "fits/", name, "a.", ncomp, ".rds")
  saveRDS(a, file = fname)
  
  a <- get(paste0("a.", ncomp, ".d"))
  
  fname <- paste0(folder, "fits/", name, "a.", ncomp, ".d", ".rds")
  saveRDS(a, file = fname)
  
}
  


# LASSO GLM options to consider:
#-----------------------------------------------#
# 1. bi- vs. multi-nomial
# 2. pure lasso vs. "strong" elastic net
# 3. global then tissue-specific on macs
# 4. with A matrix values vs. gene expression values !
# 5. Group-Overlap-LASSO: with module genes at different sd (this one would be apart - try "grpregOverlap" (extension of "grpreg") or "mlgl" packages)!

# Below: 

### I: multi-nomial, -pure lasso, global, .A mat values (MLGA)
### Ia: multi-nomial, -pure lasso, global to tissue-specific, .A mat values (MLG2A)
### II: multi-nomial, -"strong" elastic net, global, .A mat values (MEGA)
### IIa: multi-nomial, -"strong" elastic net, global to tissue-specific, .A mat values (MEG2A)
### III: multi-nomial, -pure lasso, global, .gene expression values (MLGExp)
### IV: multi-nomial, -"strong" elastic net, global, .gene expression values (MEGExp)
### IIIa: multi-nomial, -pure lasso, global to tissue-specific, .gene expression values (MLG2Exp)
### IVa: multi-nomial, -"strong" elastic net, global to tissue-specific, .gene expression values (MEG2Exp)
### V: multi-nomial, pure lasso, -tissue-specific, A mat values (MLTA)
### VI: multi-nomial, pure lasso, -tissue-specific, gene expression values (MLTExp)

### VII: bi-nomial, pure lasso, global-then tissue-specific, gene exp values (maybe mac vs. others) then (mac_spl vs. other macs)
### VIII (Overlap-LASSO): bi-nomial, pure lasso, global-then tissue-specific, gene exp values (modules!) (maybe mac vs. others) then (mac_spl vs. other macs)




library("glmnet")
# try parallel
library(parallel)
detectCores() # 4
detectCores(logical = FALSE) # 2 (2 physical, but can do 4...)
library(doParallel)
cores <- makeCluster(detectCores(logical = FALSE), type='PSOCK') # fine to just put 2 or 4 in braces




### I: multi-nomial, -pure lasso, global, .A mat values (MLGA)
#-----------------------------------------------------------------------------#


######################## Create loop for just saving fit.cv

cores <- makeCluster(detectCores(logical = FALSE), type='PSOCK')
registerDoParallel(cores)

for (ncomp in c(100, 200, 300)) {
  
  a <- get(paste0("a.", ncomp))
  
  # global
  training <- a[ inTrain.tissue, ]
  testing  <- a[-inTrain.tissue, ]
  
  print(dplyr::count(training, use_celltype))
  print(dplyr::count(testing, use_celltype))
  
  x <- as.matrix(training[, as.character(1:ncomp)]) 
  y <- training$use_celltype 
  
  # Fit & Cross-validation
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  # here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
  fname <- paste0(folder, "fits/", name, "I_MLGA_fitCV_", ncomp, ".rds")
  saveRDS(fit.cv, file = fname)
  plot(fit.cv, main = paste("I", ncomp))
  
}


for (ncomp in c(100, 200, 300)) {
  
  a <- get(paste0("a.", ncomp, ".d"))
  
  # global
  training <- a[ inTrain.tissue, ]
  testing  <- a[-inTrain.tissue, ]
  
  print(dplyr::count(training, use_celltype))
  print(dplyr::count(testing, use_celltype))
  
  x <- as.matrix(training[, as.character(1:ncomp)]) 
  y <- training$use_celltype 
  
  # Fit & Cross-validation
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  # here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
  fname <- paste0(folder, "fits/", name, "I_MLGA_fitCV_", ncomp, "_D.rds")
  saveRDS(fit.cv, file = fname)
  plot(fit.cv, main = paste("I D", ncomp))
}
print(paste("I done", Sys.time())) # ~ 1hr 10min

######################## 



### Ia: multi-nomial, -pure lasso, global to tissue-specific, .A mat values (MLG2A)
#-----------------------------------------------------------------------------#
# here filter on Macrophages only and do second level: tissue-specific marophage modules 

######################## Create loop for just saving fit.cv

for (ncomp in c(100, 200, 300)) {
  
  a <- get(paste0("a.", ncomp))
  
  # tissue-specific macs only 
  training <- a[ inTrain.tissue, ] %>% dplyr::filter(str_detect(use_celltype_tissue, "macrophage_")) %>% dplyr::mutate(use_celltype_tissue = fct_drop(use_celltype_tissue))
  testing  <- a[-inTrain.tissue, ] %>% dplyr::filter(str_detect(use_celltype_tissue, "macrophage_")) %>% dplyr::mutate(use_celltype_tissue = fct_drop(use_celltype_tissue))
  
  print(dplyr::count(training, use_celltype_tissue))
  print(dplyr::count(testing, use_celltype_tissue))
  
  x <- as.matrix(training[, as.character(1:ncomp)]) 
  y <- training$use_celltype_tissue 
  
  # Fit & Cross-validation
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  # here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
  fname <- paste0(folder, "fits/", name, "Ia_MLG2A_fitCV_", ncomp, ".rds")
  saveRDS(fit.cv, file = fname)
  plot(fit.cv, main = paste("Ia", ncomp))
  
}


for (ncomp in c(100, 200, 300)) {
  
  a <- get(paste0("a.", ncomp, ".d"))
  
  # tissue-specific macs only
  training <- a[ inTrain.tissue, ] %>% dplyr::filter(str_detect(use_celltype_tissue, "macrophage_")) %>% dplyr::mutate(use_celltype_tissue = fct_drop(use_celltype_tissue))
  testing  <- a[-inTrain.tissue, ] %>% dplyr::filter(str_detect(use_celltype_tissue, "macrophage_")) %>% dplyr::mutate(use_celltype_tissue = fct_drop(use_celltype_tissue))
  
  print(dplyr::count(training, use_celltype_tissue))
  print(dplyr::count(testing, use_celltype_tissue))
  
  x <- as.matrix(training[, as.character(1:ncomp)]) 
  y <- training$use_celltype_tissue 
  
  # Fit & Cross-validation
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  # here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
  fname <- paste0(folder, "fits/", name, "Ia_MLG2A_fitCV_", ncomp, "_D.rds")
  saveRDS(fit.cv, file = fname)
  plot(fit.cv, main = paste("Ia D", ncomp))
}
print(paste("Ia done", Sys.time())) # ~ 2min ~

######################## 



### II: multi-nomial, -"strong" elastic net, global, .A mat values (MEGA)
#-----------------------------------------------------------------------------#


######################## Create loop for just saving fit.cv

for (ncomp in c(100, 200, 300)) {
  
  a <- get(paste0("a.", ncomp))
  
  # global
  training <- a[ inTrain.tissue, ]
  testing  <- a[-inTrain.tissue, ]
  
  print(dplyr::count(training, use_celltype))
  print(dplyr::count(testing, use_celltype))
  
  x <- as.matrix(training[, as.character(1:ncomp)]) 
  y <- training$use_celltype 
  
  # Fit & Cross-validation
  fit.cv <- cv.glmnet(x=x, y=y, alpha=0.95, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  # here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
  fname <- paste0(folder, "fits/", name, "II_MEGA_fitCV_", ncomp, ".rds")
  saveRDS(fit.cv, file = fname)
  plot(fit.cv, main = paste("II", ncomp))
  
}


for (ncomp in c(100, 200, 300)) {
  
  a <- get(paste0("a.", ncomp, ".d"))
  
  # global
  training <- a[ inTrain.tissue, ]
  testing  <- a[-inTrain.tissue, ]
  
  print(dplyr::count(training, use_celltype))
  print(dplyr::count(testing, use_celltype))
  
  x <- as.matrix(training[, as.character(1:ncomp)]) 
  y <- training$use_celltype
  
  # Fit & Cross-validation
  fit.cv <- cv.glmnet(x=x, y=y, alpha=0.95, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  # here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
  fname <- paste0(folder, "fits/", name, "II_MEGA_fitCV_", ncomp, "_D.rds")
  saveRDS(fit.cv, file = fname)
  plot(fit.cv, main = paste("II D", ncomp))
}
print(paste("II done", Sys.time())) # ~ 1hr 10min

######################## 



### IIa: multi-nomial, -"strong" elastic net, global to tissue-specific, .A mat values (MEG2A)
#-----------------------------------------------------------------------------#
# here filter on Macrophages only and do second level: tissue-specific marophage modules 

######################## Create loop for just saving fit.cv

for (ncomp in c(100, 200, 300)) {
  
  a <- get(paste0("a.", ncomp))
  
  # tissue-specific macs only 
  training <- a[ inTrain.tissue, ] %>% dplyr::filter(str_detect(use_celltype_tissue, "macrophage_")) %>% dplyr::mutate(use_celltype_tissue = fct_drop(use_celltype_tissue))
  testing  <- a[-inTrain.tissue, ] %>% dplyr::filter(str_detect(use_celltype_tissue, "macrophage_")) %>% dplyr::mutate(use_celltype_tissue = fct_drop(use_celltype_tissue))
  
  print(dplyr::count(training, use_celltype_tissue))
  print(dplyr::count(testing, use_celltype_tissue))
  
  x <- as.matrix(training[, as.character(1:ncomp)]) 
  y <- training$use_celltype_tissue 
  
  # Fit & Cross-validation
  fit.cv <- cv.glmnet(x=x, y=y, alpha=0.95, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  # here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
  fname <- paste0(folder, "fits/", name, "IIa_MEG2A_fitCV_", ncomp, ".rds")
  saveRDS(fit.cv, file = fname)
  plot(fit.cv, main = paste("IIa", ncomp))
  
}


for (ncomp in c(100, 200, 300)) {
  
  a <- get(paste0("a.", ncomp, ".d"))
  
  # tissue-specific macs only 
  training <- a[ inTrain.tissue, ] %>% dplyr::filter(str_detect(use_celltype_tissue, "macrophage_")) %>% dplyr::mutate(use_celltype_tissue = fct_drop(use_celltype_tissue))
  testing  <- a[-inTrain.tissue, ] %>% dplyr::filter(str_detect(use_celltype_tissue, "macrophage_")) %>% dplyr::mutate(use_celltype_tissue = fct_drop(use_celltype_tissue))
  
  print(dplyr::count(training, use_celltype_tissue))
  print(dplyr::count(testing, use_celltype_tissue))
  
  x <- as.matrix(training[, as.character(1:ncomp)]) 
  y <- training$use_celltype_tissue 
  
  # Fit & Cross-validation
  fit.cv <- cv.glmnet(x=x, y=y, alpha=0.95, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  # here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
  fname <- paste0(folder, "fits/", name, "IIa_MEG2A_fitCV_", ncomp, "_D.rds")
  saveRDS(fit.cv, file = fname)
  plot(fit.cv, main = paste("IIa D", ncomp))
}
print(paste("IIa done", Sys.time())) # ~ 2min ~

######################## 



### III: multi-nomial, -pure lasso, global, .gene expression values (MLGExp)
#-----------------------------------------------------------------------------#


# global 
exp.training <- exp.df[ inTrain.tissue, ]
exp.testing  <- exp.df[-inTrain.tissue, ]

dplyr::count(exp.training, use_celltype)
dplyr::count(exp.testing, use_celltype)

x <- as.matrix(exp.training[, hvg]) # make sure hvg are for the appropriate dataset & normalisation
y <- exp.training$use_celltype

# Fit & Cross-validation
fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
# here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
fname <- paste0(folder, "fits/", name, "III_MLGExp_fitCV.rds")
saveRDS(fit.cv, file = fname)
plot(fit.cv, main = paste("III"))

print(paste("III done", Sys.time())) # ~ 1hr 15min



### IV: multi-nomial, -"strong" elastic net, global, .gene expression values (MEGExp)
#-----------------------------------------------------------------------------#

# setup is same as above

# Fit & Cross-validation
fit.cv <- cv.glmnet(x=x, y=y, alpha=0.95, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
# here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
fname <- paste0(folder, "fits/", name, "IV_MEGExp_fitCV.rds")
saveRDS(fit.cv, file = fname)
plot(fit.cv, main = paste("IV"))

print(paste("IV done", Sys.time())) # ~ 1hr 10min



### IIIa: multi-nomial, -pure lasso, global to tissue-specific, .gene expression values (MLG2Exp)
#-----------------------------------------------------------------------------#
# here filter on Macrophages only and do second level: tissue-specific marophage modules 

# tissue-specific macs only 
exp.training <- exp.df[ inTrain.tissue, ] %>% dplyr::filter(str_detect(use_celltype_tissue, "macrophage_")) %>% dplyr::mutate(use_celltype_tissue = fct_drop(use_celltype_tissue))
exp.testing  <- exp.df[-inTrain.tissue, ] %>% dplyr::filter(str_detect(use_celltype_tissue, "macrophage_")) %>% dplyr::mutate(use_celltype_tissue = fct_drop(use_celltype_tissue))

dplyr::count(exp.training, use_celltype_tissue)
dplyr::count(exp.testing, use_celltype_tissue)

x <- as.matrix(exp.training[, hvg]) # make sure hvg are for the appropriate dataset & normalisation
y <- exp.training$use_celltype_tissue

# Fit & Cross-validation
fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
# here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
fname <- paste0(folder, "fits/", name, "IIIa_MLG2Exp_fitCV.rds")
saveRDS(fit.cv, file = fname)
plot(fit.cv, main = paste("IIIa"))

print(paste("IIIa done", Sys.time())) # ~ 4min ~



### IVa: multi-nomial, -"strong" elastic net, global to tissue-specific, .gene expression values (MEG2Exp)
#-----------------------------------------------------------------------------#
# here filter on Macrophages only and do second level: tissue-specific marophage modules 

# setup is same as above

# Fit & Cross-validation
fit.cv <- cv.glmnet(x=x, y=y, alpha=0.95, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
# here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
fname <- paste0(folder, "fits/", name, "IVa_MEG2Exp_fitCV.rds")
saveRDS(fit.cv, file = fname)
plot(fit.cv, main = paste("IVa"))

print(paste("IVa done", Sys.time())) # ~ 4min ~



### V: multi-nomial, pure lasso, -tissue-specific, A mat values (MLTA)
#-----------------------------------------------------------------------------#


######################## Create loop for just saving fit.cv

for (ncomp in c(100, 200, 300)) {
  
  a <- get(paste0("a.", ncomp))
  
  # tissue-specific !
  training <- a[ inTrain.tissue, ]
  testing  <- a[-inTrain.tissue, ]
  
  print(dplyr::count(training, use_celltype_tissue))
  print(dplyr::count(testing, use_celltype_tissue))
  
  x <- as.matrix(training[, as.character(1:ncomp)]) 
  y <- training$use_celltype_tissue 
  
  # Fit & Cross-validation
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  # here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
  fname <- paste0(folder, "fits/", name, "V_MLTA_fitCV_", ncomp, ".rds")
  saveRDS(fit.cv, file = fname)
  plot(fit.cv, main = paste("V", ncomp))
  
}


for (ncomp in c(100, 200, 300)) {
  
  a <- get(paste0("a.", ncomp, ".d"))
  
  # tissue-specific !
  training <- a[ inTrain.tissue, ]
  testing  <- a[-inTrain.tissue, ]
  
  print(dplyr::count(training, use_celltype_tissue))
  print(dplyr::count(testing, use_celltype_tissue))
  
  x <- as.matrix(training[, as.character(1:ncomp)]) 
  y <- training$use_celltype_tissue 
  
  # Fit & Cross-validation
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  # here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
  fname <- paste0(folder, "fits/", name, "V_MLTA_fitCV_", ncomp, "_D.rds")
  saveRDS(fit.cv, file = fname)
  plot(fit.cv, main = paste("V D", ncomp))
}
print(paste("V done", Sys.time())) # ~ 1hr 30min

######################## 



### VI: multi-nomial, pure lasso, -tissue-specific, gene expression values (MLTExp)
#-----------------------------------------------------------------------------#

# global 
exp.training <- exp.df[ inTrain.tissue, ]
exp.testing  <- exp.df[-inTrain.tissue, ]

dplyr::count(exp.training, use_celltype_tissue)
dplyr::count(exp.testing, use_celltype_tissue)

x <- as.matrix(exp.training[, hvg]) # make sure hvg are for the appropriate dataset & normalisation
y <- exp.training$use_celltype_tissue

# Fit & Cross-validation
fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
# here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
fname <- paste0(folder, "fits/", name, "VI_MLTExp_fitCV.rds")
saveRDS(fit.cv, file = fname)
plot(fit.cv, main = paste("VI", ncomp))

print(paste("VI done", Sys.time())) # 2hr

stopCluster(cores)

# Chunks I-VI done in ~ 9hr 30min ~



### VII: bi-nomial, pure lasso, global (maybe mac vs. others) then (mac_spl vs. mac_kidney)
#-----------------------------------------------------------------------------#



### VIII (Overlap-LASSO): bi-nomial, pure lasso, global (maybe mac vs. others) then (mac_spl vs. mac_kidney)
#-----------------------------------------------------------------------------#





# Now analysis of all the saved model fits can be performed in the loop ...

# The best is to follow the order and structure of the model fitting above -
# as testing data sets need to be created on the go depending on the I-VI analysis

# Each analysis consists of two parts.. 
# in the 1st: extract coefficients, save them, calculate sparsity and npar per celltype (tuck in list), plot coefficient heatmap, save
# in the 2nd: get unseen test sample (15%), predict classes and probabilities, get confussion matrix, save, get BER Bal.Accuracy and per celltype Bal.Accuracy (tuck in list)


########################
########################
######-UNCOMENT-########



library("ComplexHeatmap")
library("circlize")
library("mixOmics")
library("caret")

# Start a list with 6 lists for I-VI analyses
A_stats <- list()
A_stats$I <- list()


### I analysis: (MLGA)
#-----------------------------------------------------------------------------#

# for I II
outcomes <- c("granulocyte", "promonocyte", "monocyte", "macrophage", "dendritic_cell", "mast_cell", "other")
#-nbetas:
A_stats$I$nbetas <- matrix(0, nrow = 3, ncol = length(outcomes)-1)
colnames(A_stats$I$nbetas) <- outcomes[1:length(outcomes)-1]
rownames(A_stats$I$nbetas) <- c(100,200,300)
#-ave.sparsity:
A_stats$I$ave.sparsity <- vector(mode = "numeric", length = 3)
names(A_stats$I$ave.sparsity) <- c(100,200,300)
#-ave.ber:
A_stats$I$ave.ber <- vector(mode = "numeric", length = 3)
names(A_stats$I$ave.ber) <- c(100,200,300)
#-ave.kappa:
A_stats$I$ave.kappa <- vector(mode = "numeric", length = 3)
names(A_stats$I$ave.kappa) <- c(100,200,300)
#-ave.bacc:
A_stats$I$ave.bacc <- matrix(0, nrow = 3, ncol = 3)
colnames(A_stats$I$ave.bacc) <- c("Accuracy", "loCI", "hiCI")
rownames(A_stats$I$ave.bacc) <- c(100,200,300)
#-bacc:
A_stats$I$bacc <- matrix(0, nrow = 3, ncol = length(outcomes))
colnames(A_stats$I$bacc) <- outcomes[1:length(outcomes)]
rownames(A_stats$I$bacc) <- c(100,200,300)
# duplicate for .d analysis:
A_stats$I_d <- A_stats$I


# Saves: coefficients table, coefficients heatmap, confussion table; and builds a list with stats 
for (ncomp in c(100, 200, 300)) {
  # import a. if not loaded
  a <- get(paste0("a.", ncomp))
  testing <- testing  <- a[-inTrain.tissue, ]
  # import cv.fit
  fname <- paste0(folder, "fits/", name, "I_MLGA_fitCV_", ncomp, ".rds")
  fit.cv <- readRDS(file = fname)
  
  # use fuctions: plot(), print(), coef(), predict()
  
  
  #---All betas and analyses are at lambda.1se (lambda.min would result in less sparsity)
  betas <- coef(fit.cv, s = "lambda.1se") # a list
  outcomes <- names(betas)
  
  betas <- purrr::map(betas, .f = as.matrix)
  betas <- do.call(cbind, betas)
  colnames(betas) <- outcomes
  
  betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble()
  betas <- betas %>% dplyr::filter_at(vars(-compID, -other), any_vars(.!=0)) # disregard "other only" != 0 coefficients
  print(betas, n=ncomp) # also kable it !
  ave.spars <- nrow(betas)-1
  
  fname <- paste0(folder, "fits/", name, "I_MLGA_betas_", ncomp, ".csv")
  write_excel_csv(betas, path = fname)
  
  # make heatplot for easier assesment of coefficients
  # library("ComplexHeatmap")
  b.mat <- as.matrix(betas[2:nrow(betas), 2:ncol(betas)])
  rownames(b.mat) <- betas$compID[2:nrow(betas)]
  colnames(b.mat) <- colnames(b.mat) %>% str_replace_all("_", " ")
  b.mat <- t(b.mat)
  
  # maybe exclude weak coefficients from visualisation?.. 0.5sd away from zero?..
  plot(density(b.mat))
  cutof <- sd(b.mat)/2
  
  # library("circlize")
  col_fun = colorRamp2(breaks = c(min(b.mat)+2, -cutof, cutof, max(b.mat)-2), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
  cols <- col_fun(seq(min(b.mat), max(b.mat), by = 2))
  ratio <- ncol(b.mat) / nrow(b.mat)
  
  #---first see how to exclude small values +/- ~ 0.5 or 0.25 sd and less !
  ht <- Heatmap(b.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
                row_names_max_width = max_text_width(rownames(b.mat)),
                column_dend_height = unit(0.7, "cm"), row_dend_width = unit(0.8, "cm"), name = "",
                clustering_distance_rows="spearman", clustering_distance_columns="spearman",
                clustering_method_rows="average", clustering_method_columns="average",
                width = unit(10, "in"), height = unit(2*10/ratio, "in"),
                show_heatmap_legend=FALSE, na_col = "white", col=cols)
  # col=cols, cluster_rows=F, cluster_columns=F, clustering_distance_rows="spearman", clustering_method_rows="ward.D", "ward.D2", "single", "complete", "average"
  ht <- draw(ht)
  w <- as.numeric(ComplexHeatmap:::width(ht)) / 25.2 # mm to in.
  h <- as.numeric(ComplexHeatmap:::height(ht)) / 25.2 # mm to in.
  
  fname <- paste0(folder, "fits/", name, "I_MLGA_heat_", ncomp, ".pdf")
  # pdf(file = fname, width = w, height = h)
  ht
  dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
  dev.off()
  
  
  #---------Predictions and Assesment-------------#
  newx <- as.matrix(testing[, as.character(1:ncomp)])
  
  #--Number (and indices) of non-zero coefficients
  prediction.non0 <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "nonzero") # length of list elements are numbers of coefficients at lambda.1se!
  
  #--Predicted classes
  prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
  prediction <- factor(prediction, levels = levels(testing$use_celltype))
  
  #--Predicted probabilities
  prediction.prob <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "response") # array n.cells x n.factors x 1 
  prediction.prob <- as.matrix(prediction.prob[, 1:length(outcomes), 1]) %>% as.data.frame() %>% as_tibble()
  # assess if the prediction probabilities are high enough
  prediction.prob <- prediction.prob %>% 
    add_column(prediction = prediction) %>% 
    dplyr::mutate(max.prob = pmax(granulocyte,promonocyte,monocyte,macrophage,dendritic_cell,mast_cell,other)) %>% 
    dplyr::mutate(prediction_alt = if_else(max.prob > 0.55, as.character(prediction), "no_confidence"))
  # see how many cells have < 0.55 pred. prob.
  prediction.prob %>% dplyr::filter(max.prob < 0.55)
  # can try to use this alternative predicted classes with "no_confidence" class in confusion statistics
  # prediction_alt <- pull(prediction.prob, prediction_alt)
  
  # library("mixOmics") # alt. with more options: "measures" package
  confussion <- get.confusion_matrix(truth = testing$use_celltype, predicted = prediction)
  ave.ber <- get.BER(confussion)
  ave.ber
  
  pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
  confussion <- as.data.frame(confussion) %>% rownames_to_column()
  colnames(confussion) <- c("(Truth)", pred.names)
  print(confussion) # also kable it !
  
  fname <- paste0(folder, "fits/", name, "I_MLGA_confussion_", ncomp, ".csv")
  write_excel_csv(confussion, path = fname)
  
  # confussion_alt <- get.confusion_matrix(truth = testing$use_celltype, predicted = prediction_alt)
  
  # library("caret")
  confussion.c <- confusionMatrix(data=prediction, reference=testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
  confussion.c$overall[1:2] %>% round(digits = 3)
  confussionByClass <-confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3)
  # Balanced Accuracy is derived from Sensitivity & Specificity
  # F1 is derived form Precision & Recall
  
  # confussion_alt.c <- confusionMatrix(data=prediction_alt, reference=testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") # or "prec_recall", or "everything"
  
  # Tuck-in some of these stats in A_stats list, save and plot later
  # within them: nbetas, ave.sparsity, ave.ber, ave.kappa, ave.bacc, bacc
  #-nbetas:
  nbeta <- purrr::map_depth(prediction.non0, 2, .f=length) %>% unlist() %>% .[1:length(outcomes)-1]
  j <- ncomp/100
  A_stats$I$nbetas[j, ] <- nbeta
  #-ave.sparsity:
  A_stats$I$ave.sparsity[j] <- ave.spars
  #-ave.ber:
  A_stats$I$ave.ber[j] <- ave.ber
  #-ave.kappa:
  A_stats$I$ave.kappa[j] <- confussion.c$overall[2]
  #-ave.bacc:
  ave.bac <- confussion.c$overall[c(1,3:4)]
  A_stats$I$ave.bacc[j, ] <- ave.bac
  #-bacc:
  A_stats$I$bacc[j, ] <- confussionByClass[,6]
  # save this at the end of all I-VI analyses !
  
}
  
# .d
for (ncomp in c(100, 200, 300)) {
  # import a. if not loaded
  a <- get(paste0("a.", ncomp, ".d"))
  testing <- testing  <- a[-inTrain.tissue, ]
  # import cv.fit
  fname <- paste0(folder, "fits/", name, "I_MLGA_fitCV_", ncomp, "_D", ".rds")
  fit.cv <- readRDS(file = fname)
  
  # use fuctions: plot(), print(), coef(), predict()
  
  
  #---All betas and analyses are at lambda.1se (lambda.min would result in less sparsity)
  betas <- coef(fit.cv, s = "lambda.1se") # a list
  outcomes <- names(betas)
  
  betas <- purrr::map(betas, .f = as.matrix)
  betas <- do.call(cbind, betas)
  colnames(betas) <- outcomes
  
  betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble()
  betas <- betas %>% dplyr::filter_at(vars(-compID, -other), any_vars(.!=0)) # disregard "other only" != 0 coefficients
  print(betas, n=ncomp) # also kable it !
  ave.spars <- nrow(betas)-1
  
  fname <- paste0(folder, "fits/", name, "I_MLGA_betas_", ncomp, "_D", ".csv")
  write_excel_csv(betas, path = fname)
  
  # make heatplot for easier assesment of coefficients
  # library("ComplexHeatmap")
  b.mat <- as.matrix(betas[2:nrow(betas), 2:ncol(betas)])
  rownames(b.mat) <- betas$compID[2:nrow(betas)]
  colnames(b.mat) <- colnames(b.mat) %>% str_replace_all("_", " ")
  b.mat <- t(b.mat)
  
  # maybe exclude weak coefficients from visualisation?.. 0.5sd away from zero?..
  plot(density(b.mat))
  cutof <- sd(b.mat)/2
  
  # library("circlize")
  col_fun = colorRamp2(breaks = c(min(b.mat)+2, -cutof, cutof, max(b.mat)-2), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
  cols <- col_fun(seq(min(b.mat), max(b.mat), by = 2))
  ratio <- ncol(b.mat) / nrow(b.mat)
  
  #---first see how to exclude small values +/- ~ 0.5 or 0.25 sd and less !
  ht <- Heatmap(b.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
                row_names_max_width = max_text_width(rownames(b.mat)),
                column_dend_height = unit(0.7, "cm"), row_dend_width = unit(0.8, "cm"), name = "",
                clustering_distance_rows="spearman", clustering_distance_columns="spearman",
                clustering_method_rows="average", clustering_method_columns="average",
                width = unit(10, "in"), height = unit(2*10/ratio, "in"),
                show_heatmap_legend=FALSE, na_col = "white", col=cols)
  # col=cols, cluster_rows=F, cluster_columns=F, clustering_distance_rows="spearman", clustering_method_rows="ward.D", "ward.D2", "single", "complete", "average"
  ht <- draw(ht)
  w <- as.numeric(ComplexHeatmap:::width(ht)) / 25.2 # mm to in.
  h <- as.numeric(ComplexHeatmap:::height(ht)) / 25.2 # mm to in.
  
  fname <- paste0(folder, "fits/", name, "I_MLGA_heat_", ncomp, "_D", ".pdf")
  # pdf(file = fname, width = w, height = h)
  ht
  dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
  dev.off()
  
  
  #---------Predictions and Assesment-------------#
  newx <- as.matrix(testing[, as.character(1:ncomp)])
  
  #--Number (and indices) of non-zero coefficients
  prediction.non0 <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "nonzero") # length of list elements are numbers of coefficients at lambda.1se!
  
  #--Predicted classes
  prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
  prediction <- factor(prediction, levels = levels(testing$use_celltype))
  
  #--Predicted probabilities
  prediction.prob <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "response") # array n.cells x n.factors x 1 
  prediction.prob <- as.matrix(prediction.prob[, 1:length(outcomes), 1]) %>% as.data.frame() %>% as_tibble()
  # assess if the prediction probabilities are high enough
  prediction.prob <- prediction.prob %>% 
    add_column(prediction = prediction) %>% 
    dplyr::mutate(max.prob = pmax(granulocyte,promonocyte,monocyte,macrophage,dendritic_cell,mast_cell,other)) %>% 
    dplyr::mutate(prediction_alt = if_else(max.prob > 0.55, as.character(prediction), "no_confidence"))
  # see how many cells have < 0.55 pred. prob.
  prediction.prob %>% dplyr::filter(max.prob < 0.55)
  # can try to use this alternative predicted classes with "no_confidence" class in confusion statistics
  # prediction_alt <- pull(prediction.prob, prediction_alt)
  
  # library("mixOmics") # alt. with more options: "measures" package
  confussion <- get.confusion_matrix(truth = testing$use_celltype, predicted = prediction)
  ave.ber <- get.BER(confussion)
  ave.ber
  
  pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
  confussion <- as.data.frame(confussion) %>% rownames_to_column()
  colnames(confussion) <- c("(Truth)", pred.names)
  print(confussion) # also kable it !
  
  fname <- paste0(folder, "fits/", name, "I_MLGA_confussion_", ncomp, "_D", ".csv")
  write_excel_csv(confussion, path = fname)
  
  # confussion_alt <- get.confusion_matrix(truth = testing$use_celltype, predicted = prediction_alt)
  
  # library("caret")
  confussion.c <- confusionMatrix(data=prediction, reference=testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
  confussion.c$overall[1:2] %>% round(digits = 3)
  confussionByClass <-confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3)
  # Balanced Accuracy is derived from Sensitivity & Specificity
  # F1 is derived form Precision & Recall
  
  # confussion_alt.c <- confusionMatrix(data=prediction_alt, reference=testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") # or "prec_recall", or "everything"
  
  # Tuck-in some of these stats in A_stats list, save and plot later
  # within them: nbetas, ave.sparsity, ave.ber, ave.kappa, ave.bacc, bacc
  #-nbetas:
  nbeta <- purrr::map_depth(prediction.non0, 2, .f=length) %>% unlist() %>% .[1:length(outcomes)-1]
  j <- ncomp/100
  A_stats$I_d$nbetas[j, ] <- nbeta
  #-ave.sparsity:
  A_stats$I_d$ave.sparsity[j] <- ave.spars
  #-ave.ber:
  A_stats$I_d$ave.ber[j] <- ave.ber
  #-ave.kappa:
  A_stats$I_d$ave.kappa[j] <- confussion.c$overall[2]
  #-ave.bacc:
  ave.bac <- confussion.c$overall[c(1,3:4)]
  A_stats$I_d$ave.bacc[j, ] <- ave.bac
  #-bacc:
  A_stats$I_d$bacc[j, ] <- confussionByClass[,6]
  # save this at the end of all I-VI analyses !
  
}



### Ia analysis: (MLG2A)
#-----------------------------------------------------------------------------#

# for Ia IIa
outcomes <- c("macrophage_Marrow", "macrophage_Spleen", "macrophage_Lung", "macrophage_Kidney", "macrophage_Mammary_Gland", "macrophage_Limb_Muscle")

#-nbetas:
A_stats$Ia$nbetas <- matrix(0, nrow = 3, ncol = length(outcomes))
colnames(A_stats$Ia$nbetas) <- outcomes[1:length(outcomes)]
rownames(A_stats$Ia$nbetas) <- c(100,200,300)
#-ave.sparsity:
A_stats$Ia$ave.sparsity <- vector(mode = "numeric", length = 3)
names(A_stats$Ia$ave.sparsity) <- c(100,200,300)
#-ave.ber:
A_stats$Ia$ave.ber <- vector(mode = "numeric", length = 3)
names(A_stats$Ia$ave.ber) <- c(100,200,300)
#-ave.kappa:
A_stats$Ia$ave.kappa <- vector(mode = "numeric", length = 3)
names(A_stats$Ia$ave.kappa) <- c(100,200,300)
#-ave.bacc:
A_stats$Ia$ave.bacc <- matrix(0, nrow = 3, ncol = 3)
colnames(A_stats$Ia$ave.bacc) <- c("Accuracy", "loCI", "hiCI")
rownames(A_stats$Ia$ave.bacc) <- c(100,200,300)
#-bacc:
A_stats$Ia$bacc <- matrix(0, nrow = 3, ncol = length(outcomes))
colnames(A_stats$Ia$bacc) <- outcomes[1:length(outcomes)]
rownames(A_stats$Ia$bacc) <- c(100,200,300)
# duplicate for .d analysis:
A_stats$Ia_d <- A_stats$Ia


for (ncomp in c(100, 200, 300)) {
  # import a. if not loaded
  a <- get(paste0("a.", ncomp))
  testing  <- a[-inTrain.tissue, ] %>% dplyr::filter(str_detect(use_celltype_tissue, "macrophage_")) %>% dplyr::mutate(use_celltype_tissue = fct_drop(use_celltype_tissue))
  # import cv.fit
  fname <- paste0(folder, "fits/", name, "Ia_MLG2A_fitCV_", ncomp, ".rds")
  fit.cv <- readRDS(file = fname)
  
  # use fuctions: plot(), print(), coef(), predict()
  
  
  #---All betas and analyses are at lambda.1se (lambda.min would result in less sparsity)
  betas <- coef(fit.cv, s = "lambda.1se") # a list
  outcomes <- names(betas)
  
  betas <- purrr::map(betas, .f = as.matrix)
  betas <- do.call(cbind, betas)
  colnames(betas) <- outcomes
  
  betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble()
  betas <- betas %>% dplyr::filter_at(vars(-compID), any_vars(.!=0)) # disregard "other only" != 0 coefficients
  print(betas, n=ncomp) # also kable it !
  ave.spars <- nrow(betas)-1
  
  fname <- paste0(folder, "fits/", name, "Ia_MLG2A_betas_", ncomp, ".csv")
  write_excel_csv(betas, path = fname)
  
  # make heatplot for easier assesment of coefficients
  # library("ComplexHeatmap")
  b.mat <- as.matrix(betas[2:nrow(betas), 2:ncol(betas)])
  rownames(b.mat) <- betas$compID[2:nrow(betas)]
  colnames(b.mat) <- colnames(b.mat) %>% str_replace_all("_", " ")
  b.mat <- t(b.mat)
  
  # maybe exclude weak coefficients from visualisation?.. 0.5sd away from zero?..
  plot(density(b.mat))
  cutof <- sd(b.mat)/2
  
  # library("circlize")
  col_fun = colorRamp2(breaks = c(min(b.mat)+2, -cutof, cutof, max(b.mat)-2), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
  cols <- col_fun(seq(min(b.mat), max(b.mat), by = 2))
  ratio <- ncol(b.mat) / nrow(b.mat)
  
  #---first see how to exclude small values +/- ~ 0.5 or 0.25 sd and less !
  ht <- Heatmap(b.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
                row_names_max_width = max_text_width(rownames(b.mat)),
                column_dend_height = unit(0.7, "cm"), row_dend_width = unit(0.8, "cm"), name = "",
                clustering_distance_rows="spearman", clustering_distance_columns="spearman",
                clustering_method_rows="average", clustering_method_columns="average",
                width = unit(10, "in"), height = unit(2*10/ratio, "in"),
                show_heatmap_legend=FALSE, na_col = "white", col=cols)
  # col=cols, cluster_rows=F, cluster_columns=F, clustering_distance_rows="spearman", clustering_method_rows="ward.D", "ward.D2", "single", "complete", "average"
  ht <- draw(ht)
  w <- as.numeric(ComplexHeatmap:::width(ht)) / 25.2 # mm to in.
  h <- as.numeric(ComplexHeatmap:::height(ht)) / 25.2 # mm to in.
  
  fname <- paste0(folder, "fits/", name, "Ia_MLG2A_heat_", ncomp, ".pdf")
  # pdf(file = fname, width = w, height = h)
  ht
  dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
  dev.off()
  
  
  #---------Predictions and Assesment-------------#
  newx <- as.matrix(testing[, as.character(1:ncomp)])
  
  #--Number (and indices) of non-zero coefficients
  prediction.non0 <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "nonzero") # length of list elements are numbers of coefficients at lambda.1se!
  
  #--Predicted classes
  prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
  prediction <- factor(prediction, levels = levels(testing$use_celltype_tissue))
  
  #--Predicted probabilities
  prediction.prob <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "response") # array n.cells x n.factors x 1 
  prediction.prob <- as.matrix(prediction.prob[, 1:length(outcomes), 1]) %>% as.data.frame() %>% as_tibble()
  # assess if the prediction probabilities are high enough
  prediction.prob <- prediction.prob %>% 
    add_column(prediction = prediction) %>% 
    dplyr::mutate(max.prob = pmax(macrophage_Marrow,macrophage_Spleen,macrophage_Lung,macrophage_Kidney,macrophage_Mammary_Gland,macrophage_Limb_Muscle)) %>% 
    dplyr::mutate(prediction_alt = if_else(max.prob > 0.55, as.character(prediction), "no_confidence"))
  # see how many cells have < 0.55 pred. prob.
  prediction.prob %>% dplyr::filter(max.prob < 0.55)
  # can try to use this alternative predicted classes with "no_confidence" class in confusion statistics
  # prediction_alt <- pull(prediction.prob, prediction_alt)
  
  # library("mixOmics") # alt. with more options: "measures" package
  confussion <- get.confusion_matrix(truth = testing$use_celltype_tissue, predicted = prediction)
  ave.ber <- get.BER(confussion)
  ave.ber
  
  pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
  confussion <- as.data.frame(confussion) %>% rownames_to_column()
  colnames(confussion) <- c("(Truth)", pred.names)
  print(confussion) # also kable it !
  
  fname <- paste0(folder, "fits/", name, "Ia_MLG2A_confussion_", ncomp, ".csv")
  write_excel_csv(confussion, path = fname)
  
  # confussion_alt <- get.confusion_matrix(truth = testing$use_celltype, predicted = prediction_alt)
  
  # library("caret")
  confussion.c <- confusionMatrix(data=prediction, reference=testing$use_celltype_tissue, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
  confussion.c$overall[1:2] %>% round(digits = 3)
  confussionByClass <-confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3)
  # Balanced Accuracy is derived from Sensitivity & Specificity
  # F1 is derived form Precision & Recall
  
  # confussion_alt.c <- confusionMatrix(data=prediction_alt, reference=testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") # or "prec_recall", or "everything"
  
  # Tuck-in some of these stats in A_stats list, save and plot later
  # within them: nbetas, ave.sparsity, ave.ber, ave.kappa, ave.bacc, bacc
  #-nbetas:
  nbeta <- purrr::map_depth(prediction.non0, 2, .f=length) %>% unlist() %>% .[1:length(outcomes)]
  j <- ncomp/100
  A_stats$Ia$nbetas[j, ] <- nbeta
  #-ave.sparsity:
  A_stats$Ia$ave.sparsity[j] <- ave.spars
  #-ave.ber:
  A_stats$Ia$ave.ber[j] <- ave.ber
  #-ave.kappa:
  A_stats$Ia$ave.kappa[j] <- confussion.c$overall[2]
  #-ave.bacc:
  ave.bac <- confussion.c$overall[c(1,3:4)]
  A_stats$Ia$ave.bacc[j, ] <- ave.bac
  #-bacc:
  A_stats$Ia$bacc[j, ] <- confussionByClass[,6]
  # save this at the end of all I-VI analyses !
  
}


# .d
for (ncomp in c(100, 200, 300)) {
  # import a. if not loaded
  a <- get(paste0("a.", ncomp, ".d"))
  testing  <- a[-inTrain.tissue, ] %>% dplyr::filter(str_detect(use_celltype_tissue, "macrophage_")) %>% dplyr::mutate(use_celltype_tissue = fct_drop(use_celltype_tissue))
  # import cv.fit
  fname <- paste0(folder, "fits/", name, "Ia_MLG2A_fitCV_", ncomp, "_D", ".rds")
  fit.cv <- readRDS(file = fname)
  
  # use fuctions: plot(), print(), coef(), predict()
  
  
  #---All betas and analyses are at lambda.1se (lambda.min would result in less sparsity)
  betas <- coef(fit.cv, s = "lambda.1se") # a list
  outcomes <- names(betas)
  
  betas <- purrr::map(betas, .f = as.matrix)
  betas <- do.call(cbind, betas)
  colnames(betas) <- outcomes
  
  betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble()
  betas <- betas %>% dplyr::filter_at(vars(-compID), any_vars(.!=0)) # disregard "other only" != 0 coefficients
  print(betas, n=ncomp) # also kable it !
  ave.spars <- nrow(betas)-1
  
  fname <- paste0(folder, "fits/", name, "Ia_MLG2A_betas_", ncomp, "_D", ".csv")
  write_excel_csv(betas, path = fname)
  
  # make heatplot for easier assesment of coefficients
  # library("ComplexHeatmap")
  b.mat <- as.matrix(betas[2:nrow(betas), 2:ncol(betas)])
  rownames(b.mat) <- betas$compID[2:nrow(betas)]
  colnames(b.mat) <- colnames(b.mat) %>% str_replace_all("_", " ")
  b.mat <- t(b.mat)
  
  # maybe exclude weak coefficients from visualisation?.. 0.5sd away from zero?..
  plot(density(b.mat))
  cutof <- sd(b.mat)/2
  
  # library("circlize")
  col_fun = colorRamp2(breaks = c(min(b.mat)+2, -cutof, cutof, max(b.mat)-2), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
  cols <- col_fun(seq(min(b.mat), max(b.mat), by = 2))
  ratio <- ncol(b.mat) / nrow(b.mat)
  
  #---first see how to exclude small values +/- ~ 0.5 or 0.25 sd and less !
  ht <- Heatmap(b.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
                row_names_max_width = max_text_width(rownames(b.mat)),
                column_dend_height = unit(0.7, "cm"), row_dend_width = unit(0.8, "cm"), name = "",
                clustering_distance_rows="spearman", clustering_distance_columns="spearman",
                clustering_method_rows="average", clustering_method_columns="average",
                width = unit(10, "in"), height = unit(2*10/ratio, "in"),
                show_heatmap_legend=FALSE, na_col = "white", col=cols)
  # col=cols, cluster_rows=F, cluster_columns=F, clustering_distance_rows="spearman", clustering_method_rows="ward.D", "ward.D2", "single", "complete", "average"
  ht <- draw(ht)
  w <- as.numeric(ComplexHeatmap:::width(ht)) / 25.2 # mm to in.
  h <- as.numeric(ComplexHeatmap:::height(ht)) / 25.2 # mm to in.
  
  fname <- paste0(folder, "fits/", name, "Ia_MLG2A_heat_", ncomp, "_D", ".pdf")
  # pdf(file = fname, width = w, height = h)
  ht
  dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
  dev.off()
  
  
  #---------Predictions and Assesment-------------#
  newx <- as.matrix(testing[, as.character(1:ncomp)])
  
  #--Number (and indices) of non-zero coefficients
  prediction.non0 <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "nonzero") # length of list elements are numbers of coefficients at lambda.1se!
  
  #--Predicted classes
  prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
  prediction <- factor(prediction, levels = levels(testing$use_celltype_tissue))
  
  #--Predicted probabilities
  prediction.prob <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "response") # array n.cells x n.factors x 1 
  prediction.prob <- as.matrix(prediction.prob[, 1:length(outcomes), 1]) %>% as.data.frame() %>% as_tibble()
  # assess if the prediction probabilities are high enough
  prediction.prob <- prediction.prob %>% 
    add_column(prediction = prediction) %>% 
    dplyr::mutate(max.prob = pmax(macrophage_Marrow,macrophage_Spleen,macrophage_Lung,macrophage_Kidney,macrophage_Mammary_Gland,macrophage_Limb_Muscle)) %>% 
    dplyr::mutate(prediction_alt = if_else(max.prob > 0.55, as.character(prediction), "no_confidence"))
  # see how many cells have < 0.55 pred. prob.
  prediction.prob %>% dplyr::filter(max.prob < 0.55)
  # can try to use this alternative predicted classes with "no_confidence" class in confusion statistics
  # prediction_alt <- pull(prediction.prob, prediction_alt)
  
  # library("mixOmics") # alt. with more options: "measures" package
  confussion <- get.confusion_matrix(truth = testing$use_celltype_tissue, predicted = prediction)
  ave.ber <- get.BER(confussion)
  ave.ber
  
  pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
  confussion <- as.data.frame(confussion) %>% rownames_to_column()
  colnames(confussion) <- c("(Truth)", pred.names)
  print(confussion) # also kable it !
  
  fname <- paste0(folder, "fits/", name, "Ia_MLG2A_confussion_", ncomp, "_D", ".csv")
  write_excel_csv(confussion, path = fname)
  
  # confussion_alt <- get.confusion_matrix(truth = testing$use_celltype, predicted = prediction_alt)
  
  # library("caret")
  confussion.c <- confusionMatrix(data=prediction, reference=testing$use_celltype_tissue, dnn = c("Prediction", "Reference"), mode = "sens_spec")
  confussion.c$overall[1:2] %>% round(digits = 3)
  confussionByClass <-confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3)
  # Balanced Accuracy is derived from Sensitivity & Specificity
  # F1 is derived form Precision & Recall
  
  # confussion_alt.c <- confusionMatrix(data=prediction_alt, reference=testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") # or "prec_recall", or "everything"
  
  # Tuck-in some of these stats in A_stats list, save and plot later
  # within them: nbetas, ave.sparsity, ave.ber, ave.kappa, ave.bacc, bacc
  #-nbetas:
  nbeta <- purrr::map_depth(prediction.non0, 2, .f=length) %>% unlist() %>% .[1:length(outcomes)]
  j <- ncomp/100
  A_stats$Ia_d$nbetas[j, ] <- nbeta
  #-ave.sparsity:
  A_stats$Ia_d$ave.sparsity[j] <- ave.spars
  #-ave.ber:
  A_stats$Ia_d$ave.ber[j] <- ave.ber
  #-ave.kappa:
  A_stats$Ia_d$ave.kappa[j] <- confussion.c$overall[2]
  #-ave.bacc:
  ave.bac <- confussion.c$overall[c(1,3:4)]
  A_stats$Ia_d$ave.bacc[j, ] <- ave.bac
  #-bacc:
  A_stats$Ia_d$bacc[j, ] <- confussionByClass[,6]
  # save this at the end of all I-VI analyses !
  
}



### II analysis: (MEGA)
#-----------------------------------------------------------------------------#

# for I II
outcomes <- c("granulocyte", "promonocyte", "monocyte", "macrophage", "dendritic_cell", "mast_cell", "other")
#-nbetas:
A_stats$II$nbetas <- matrix(0, nrow = 3, ncol = length(outcomes)-1)
colnames(A_stats$II$nbetas) <- outcomes[1:length(outcomes)-1]
rownames(A_stats$II$nbetas) <- c(100,200,300)
#-ave.sparsity:
A_stats$II$ave.sparsity <- vector(mode = "numeric", length = 3)
names(A_stats$II$ave.sparsity) <- c(100,200,300)
#-ave.ber:
A_stats$II$ave.ber <- vector(mode = "numeric", length = 3)
names(A_stats$II$ave.ber) <- c(100,200,300)
#-ave.kappa:
A_stats$II$ave.kappa <- vector(mode = "numeric", length = 3)
names(A_stats$II$ave.kappa) <- c(100,200,300)
#-ave.bacc:
A_stats$II$ave.bacc <- matrix(0, nrow = 3, ncol = 3)
colnames(A_stats$II$ave.bacc) <- c("Accuracy", "loCI", "hiCI")
rownames(A_stats$II$ave.bacc) <- c(100,200,300)
#-bacc:
A_stats$II$bacc <- matrix(0, nrow = 3, ncol = length(outcomes))
colnames(A_stats$II$bacc) <- outcomes[1:length(outcomes)]
rownames(A_stats$II$bacc) <- c(100,200,300)
# duplicate for .d analysis:
A_stats$II_d <- A_stats$II


# Saves: coefficients table, coefficients heatmap, confussion table; and builds a list with stats 
for (ncomp in c(100, 200, 300)) {
  # import a. if not loaded
  a <- get(paste0("a.", ncomp))
  testing <- testing  <- a[-inTrain.tissue, ]
  # import cv.fit
  fname <- paste0(folder, "fits/", name, "II_MEGA_fitCV_", ncomp, ".rds")
  fit.cv <- readRDS(file = fname)
  
  # use fuctions: plot(), print(), coef(), predict()
  
  
  #---All betas and analyses are at lambda.1se (lambda.min would result in less sparsity)
  betas <- coef(fit.cv, s = "lambda.1se") # a list
  outcomes <- names(betas)
  
  betas <- purrr::map(betas, .f = as.matrix)
  betas <- do.call(cbind, betas)
  colnames(betas) <- outcomes
  
  betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble()
  betas <- betas %>% dplyr::filter_at(vars(-compID, -other), any_vars(.!=0)) # disregard "other only" != 0 coefficients
  print(betas, n=ncomp) # also kable it !
  ave.spars <- nrow(betas)-1
  
  fname <- paste0(folder, "fits/", name, "II_MEGA_betas_", ncomp, ".csv")
  write_excel_csv(betas, path = fname)
  
  # make heatplot for easier assesment of coefficients
  # library("ComplexHeatmap")
  b.mat <- as.matrix(betas[2:nrow(betas), 2:ncol(betas)])
  rownames(b.mat) <- betas$compID[2:nrow(betas)]
  colnames(b.mat) <- colnames(b.mat) %>% str_replace_all("_", " ")
  b.mat <- t(b.mat)
  
  # maybe exclude weak coefficients from visualisation?.. 0.5sd away from zero?..
  plot(density(b.mat))
  cutof <- sd(b.mat)/2
  
  # library("circlize")
  col_fun = colorRamp2(breaks = c(min(b.mat)+2, -cutof, cutof, max(b.mat)-2), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
  cols <- col_fun(seq(min(b.mat), max(b.mat), by = 2))
  ratio <- ncol(b.mat) / nrow(b.mat)
  
  #---first see how to exclude small values +/- ~ 0.5 or 0.25 sd and less !
  ht <- Heatmap(b.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
                row_names_max_width = max_text_width(rownames(b.mat)),
                column_dend_height = unit(0.7, "cm"), row_dend_width = unit(0.8, "cm"), name = "",
                clustering_distance_rows="spearman", clustering_distance_columns="spearman",
                clustering_method_rows="average", clustering_method_columns="average",
                width = unit(10, "in"), height = unit(2*10/ratio, "in"),
                show_heatmap_legend=FALSE, na_col = "white", col=cols)
  # col=cols, cluster_rows=F, cluster_columns=F, clustering_distance_rows="spearman", clustering_method_rows="ward.D", "ward.D2", "single", "complete", "average"
  ht <- draw(ht)
  w <- as.numeric(ComplexHeatmap:::width(ht)) / 25.2 # mm to in.
  h <- as.numeric(ComplexHeatmap:::height(ht)) / 25.2 # mm to in.
  
  fname <- paste0(folder, "fits/", name, "II_MEGA_heat_", ncomp, ".pdf")
  # pdf(file = fname, width = w, height = h)
  ht
  dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
  dev.off()
  
  
  #---------Predictions and Assesment-------------#
  newx <- as.matrix(testing[, as.character(1:ncomp)])
  
  #--Number (and indices) of non-zero coefficients
  prediction.non0 <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "nonzero") # length of list elements are numbers of coefficients at lambda.1se!
  
  #--Predicted classes
  prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
  prediction <- factor(prediction, levels = levels(testing$use_celltype))
  
  #--Predicted probabilities
  prediction.prob <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "response") # array n.cells x n.factors x 1 
  prediction.prob <- as.matrix(prediction.prob[, 1:length(outcomes), 1]) %>% as.data.frame() %>% as_tibble()
  # assess if the prediction probabilities are high enough
  prediction.prob <- prediction.prob %>% 
    add_column(prediction = prediction) %>% 
    dplyr::mutate(max.prob = pmax(granulocyte,promonocyte,monocyte,macrophage,dendritic_cell,mast_cell,other)) %>% 
    dplyr::mutate(prediction_alt = if_else(max.prob > 0.55, as.character(prediction), "no_confidence"))
  # see how many cells have < 0.55 pred. prob.
  prediction.prob %>% dplyr::filter(max.prob < 0.55)
  # can try to use this alternative predicted classes with "no_confidence" class in confusion statistics
  # prediction_alt <- pull(prediction.prob, prediction_alt)
  
  # library("mixOmics") # alt. with more options: "measures" package
  confussion <- get.confusion_matrix(truth = testing$use_celltype, predicted = prediction)
  ave.ber <- get.BER(confussion)
  ave.ber
  
  pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
  confussion <- as.data.frame(confussion) %>% rownames_to_column()
  colnames(confussion) <- c("(Truth)", pred.names)
  print(confussion) # also kable it !
  
  fname <- paste0(folder, "fits/", name, "II_MEGA_confussion_", ncomp, ".csv")
  write_excel_csv(confussion, path = fname)
  
  # confussion_alt <- get.confusion_matrix(truth = testing$use_celltype, predicted = prediction_alt)
  
  # library("caret")
  confussion.c <- confusionMatrix(data=prediction, reference=testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
  confussion.c$overall[1:2] %>% round(digits = 3)
  confussionByClass <-confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3)
  # Balanced Accuracy is derived from Sensitivity & Specificity
  # F1 is derived form Precision & Recall
  
  # confussion_alt.c <- confusionMatrix(data=prediction_alt, reference=testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") # or "prec_recall", or "everything"
  
  # Tuck-in some of these stats in A_stats list, save and plot later
  # within them: nbetas, ave.sparsity, ave.ber, ave.kappa, ave.bacc, bacc
  #-nbetas:
  nbeta <- purrr::map_depth(prediction.non0, 2, .f=length) %>% unlist() %>% .[1:length(outcomes)-1]
  j <- ncomp/100
  A_stats$II$nbetas[j, ] <- nbeta
  #-ave.sparsity:
  A_stats$II$ave.sparsity[j] <- ave.spars
  #-ave.ber:
  A_stats$II$ave.ber[j] <- ave.ber
  #-ave.kappa:
  A_stats$II$ave.kappa[j] <- confussion.c$overall[2]
  #-ave.bacc:
  ave.bac <- confussion.c$overall[c(1,3:4)]
  A_stats$II$ave.bacc[j, ] <- ave.bac
  #-bacc:
  A_stats$II$bacc[j, ] <- confussionByClass[,6]
  # save this at the end of all I-VI analyses !
  
}

# .d
for (ncomp in c(100, 200, 300)) {
  # import a. if not loaded
  a <- get(paste0("a.", ncomp, ".d"))
  testing <- testing  <- a[-inTrain.tissue, ]
  # import cv.fit
  fname <- paste0(folder, "fits/", name, "II_MEGA_fitCV_", ncomp, "_D", ".rds")
  fit.cv <- readRDS(file = fname)
  
  # use fuctions: plot(), print(), coef(), predict()
  
  
  #---All betas and analyses are at lambda.1se (lambda.min would result in less sparsity)
  betas <- coef(fit.cv, s = "lambda.1se") # a list
  outcomes <- names(betas)
  
  betas <- purrr::map(betas, .f = as.matrix)
  betas <- do.call(cbind, betas)
  colnames(betas) <- outcomes
  
  betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble()
  betas <- betas %>% dplyr::filter_at(vars(-compID, -other), any_vars(.!=0)) # disregard "other only" != 0 coefficients
  print(betas, n=ncomp) # also kable it !
  ave.spars <- nrow(betas)-1
  
  fname <- paste0(folder, "fits/", name, "II_MEGA_betas_", ncomp, "_D", ".csv")
  write_excel_csv(betas, path = fname)
  
  # make heatplot for easier assesment of coefficients
  # library("ComplexHeatmap")
  b.mat <- as.matrix(betas[2:nrow(betas), 2:ncol(betas)])
  rownames(b.mat) <- betas$compID[2:nrow(betas)]
  colnames(b.mat) <- colnames(b.mat) %>% str_replace_all("_", " ")
  b.mat <- t(b.mat)
  
  # maybe exclude weak coefficients from visualisation?.. 0.5sd away from zero?..
  plot(density(b.mat))
  cutof <- sd(b.mat)/2
  
  # library("circlize")
  col_fun = colorRamp2(breaks = c(min(b.mat)+2, -cutof, cutof, max(b.mat)-2), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
  cols <- col_fun(seq(min(b.mat), max(b.mat), by = 2))
  ratio <- ncol(b.mat) / nrow(b.mat)
  
  #---first see how to exclude small values +/- ~ 0.5 or 0.25 sd and less !
  ht <- Heatmap(b.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
                row_names_max_width = max_text_width(rownames(b.mat)),
                column_dend_height = unit(0.7, "cm"), row_dend_width = unit(0.8, "cm"), name = "",
                clustering_distance_rows="spearman", clustering_distance_columns="spearman",
                clustering_method_rows="average", clustering_method_columns="average",
                width = unit(10, "in"), height = unit(2*10/ratio, "in"),
                show_heatmap_legend=FALSE, na_col = "white", col=cols)
  # col=cols, cluster_rows=F, cluster_columns=F, clustering_distance_rows="spearman", clustering_method_rows="ward.D", "ward.D2", "single", "complete", "average"
  ht <- draw(ht)
  w <- as.numeric(ComplexHeatmap:::width(ht)) / 25.2 # mm to in.
  h <- as.numeric(ComplexHeatmap:::height(ht)) / 25.2 # mm to in.
  
  fname <- paste0(folder, "fits/", name, "II_MEGA_heat_", ncomp, "_D", ".pdf")
  # pdf(file = fname, width = w, height = h)
  ht
  dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
  dev.off()
  
  
  #---------Predictions and Assesment-------------#
  newx <- as.matrix(testing[, as.character(1:ncomp)])
  
  #--Number (and indices) of non-zero coefficients
  prediction.non0 <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "nonzero") # length of list elements are numbers of coefficients at lambda.1se!
  
  #--Predicted classes
  prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
  prediction <- factor(prediction, levels = levels(testing$use_celltype))
  
  #--Predicted probabilities
  prediction.prob <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "response") # array n.cells x n.factors x 1 
  prediction.prob <- as.matrix(prediction.prob[, 1:length(outcomes), 1]) %>% as.data.frame() %>% as_tibble()
  # assess if the prediction probabilities are high enough
  prediction.prob <- prediction.prob %>% 
    add_column(prediction = prediction) %>% 
    dplyr::mutate(max.prob = pmax(granulocyte,promonocyte,monocyte,macrophage,dendritic_cell,mast_cell,other)) %>% 
    dplyr::mutate(prediction_alt = if_else(max.prob > 0.55, as.character(prediction), "no_confidence"))
  # see how many cells have < 0.55 pred. prob.
  prediction.prob %>% dplyr::filter(max.prob < 0.55)
  # can try to use this alternative predicted classes with "no_confidence" class in confusion statistics
  # prediction_alt <- pull(prediction.prob, prediction_alt)
  
  # library("mixOmics") # alt. with more options: "measures" package
  confussion <- get.confusion_matrix(truth = testing$use_celltype, predicted = prediction)
  ave.ber <- get.BER(confussion)
  ave.ber
  
  pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
  confussion <- as.data.frame(confussion) %>% rownames_to_column()
  colnames(confussion) <- c("(Truth)", pred.names)
  print(confussion) # also kable it !
  
  fname <- paste0(folder, "fits/", name, "II_MEGA_confussion_", ncomp, "_D", ".csv")
  write_excel_csv(confussion, path = fname)
  
  # confussion_alt <- get.confusion_matrix(truth = testing$use_celltype, predicted = prediction_alt)
  
  # library("caret")
  confussion.c <- confusionMatrix(data=prediction, reference=testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
  confussion.c$overall[1:2] %>% round(digits = 3)
  confussionByClass <-confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3)
  # Balanced Accuracy is derived from Sensitivity & Specificity
  # F1 is derived form Precision & Recall
  
  # confussion_alt.c <- confusionMatrix(data=prediction_alt, reference=testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") # or "prec_recall", or "everything"
  
  # Tuck-in some of these stats in A_stats list, save and plot later
  # within them: nbetas, ave.sparsity, ave.ber, ave.kappa, ave.bacc, bacc
  #-nbetas:
  nbeta <- purrr::map_depth(prediction.non0, 2, .f=length) %>% unlist() %>% .[1:length(outcomes)-1]
  j <- ncomp/100
  A_stats$II_d$nbetas[j, ] <- nbeta
  #-ave.sparsity:
  A_stats$II_d$ave.sparsity[j] <- ave.spars
  #-ave.ber:
  A_stats$II_d$ave.ber[j] <- ave.ber
  #-ave.kappa:
  A_stats$II_d$ave.kappa[j] <- confussion.c$overall[2]
  #-ave.bacc:
  ave.bac <- confussion.c$overall[c(1,3:4)]
  A_stats$II_d$ave.bacc[j, ] <- ave.bac
  #-bacc:
  A_stats$II_d$bacc[j, ] <- confussionByClass[,6]
  # save this at the end of all I-VI analyses !
  
}



### IIa analysis: (MEG2A)
#-----------------------------------------------------------------------------#

# for Ia IIa
outcomes <- c("macrophage_Marrow", "macrophage_Spleen", "macrophage_Lung", "macrophage_Kidney", "macrophage_Mammary_Gland", "macrophage_Limb_Muscle")

#-nbetas:
A_stats$IIa$nbetas <- matrix(0, nrow = 3, ncol = length(outcomes))
colnames(A_stats$IIa$nbetas) <- outcomes[1:length(outcomes)]
rownames(A_stats$IIa$nbetas) <- c(100,200,300)
#-ave.sparsity:
A_stats$IIa$ave.sparsity <- vector(mode = "numeric", length = 3)
names(A_stats$IIa$ave.sparsity) <- c(100,200,300)
#-ave.ber:
A_stats$IIa$ave.ber <- vector(mode = "numeric", length = 3)
names(A_stats$IIa$ave.ber) <- c(100,200,300)
#-ave.kappa:
A_stats$IIa$ave.kappa <- vector(mode = "numeric", length = 3)
names(A_stats$IIa$ave.kappa) <- c(100,200,300)
#-ave.bacc:
A_stats$IIa$ave.bacc <- matrix(0, nrow = 3, ncol = 3)
colnames(A_stats$IIa$ave.bacc) <- c("Accuracy", "loCI", "hiCI")
rownames(A_stats$IIa$ave.bacc) <- c(100,200,300)
#-bacc:
A_stats$IIa$bacc <- matrix(0, nrow = 3, ncol = length(outcomes))
colnames(A_stats$IIa$bacc) <- outcomes[1:length(outcomes)]
rownames(A_stats$IIa$bacc) <- c(100,200,300)
# duplicate for .d analysis:
A_stats$IIa_d <- A_stats$IIa


for (ncomp in c(100, 200, 300)) {
  # import a. if not loaded
  a <- get(paste0("a.", ncomp))
  testing  <- a[-inTrain.tissue, ] %>% dplyr::filter(str_detect(use_celltype_tissue, "macrophage_")) %>% dplyr::mutate(use_celltype_tissue = fct_drop(use_celltype_tissue))
  # import cv.fit
  fname <- paste0(folder, "fits/", name, "IIa_MEG2A_fitCV_", ncomp, ".rds")
  fit.cv <- readRDS(file = fname)
  
  # use fuctions: plot(), print(), coef(), predict()
  
  
  #---All betas and analyses are at lambda.1se (lambda.min would result in less sparsity)
  betas <- coef(fit.cv, s = "lambda.1se") # a list
  outcomes <- names(betas)
  
  betas <- purrr::map(betas, .f = as.matrix)
  betas <- do.call(cbind, betas)
  colnames(betas) <- outcomes
  
  betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble()
  betas <- betas %>% dplyr::filter_at(vars(-compID), any_vars(.!=0)) # disregard "other only" != 0 coefficients
  print(betas, n=ncomp) # also kable it !
  ave.spars <- nrow(betas)-1
  
  fname <- paste0(folder, "fits/", name, "IIa_MEG2A_betas_", ncomp, ".csv")
  write_excel_csv(betas, path = fname)
  
  # make heatplot for easier assesment of coefficients
  # library("ComplexHeatmap")
  b.mat <- as.matrix(betas[2:nrow(betas), 2:ncol(betas)])
  rownames(b.mat) <- betas$compID[2:nrow(betas)]
  colnames(b.mat) <- colnames(b.mat) %>% str_replace_all("_", " ")
  b.mat <- t(b.mat)
  
  # maybe exclude weak coefficients from visualisation?.. 0.5sd away from zero?..
  plot(density(b.mat))
  cutof <- sd(b.mat)/2
  
  # library("circlize")
  col_fun = colorRamp2(breaks = c(min(b.mat)+2, -cutof, cutof, max(b.mat)-2), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
  cols <- col_fun(seq(min(b.mat), max(b.mat), by = 2))
  ratio <- ncol(b.mat) / nrow(b.mat)
  
  #---first see how to exclude small values +/- ~ 0.5 or 0.25 sd and less !
  ht <- Heatmap(b.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
                row_names_max_width = max_text_width(rownames(b.mat)),
                column_dend_height = unit(0.7, "cm"), row_dend_width = unit(0.8, "cm"), name = "",
                clustering_distance_rows="spearman", clustering_distance_columns="spearman",
                clustering_method_rows="average", clustering_method_columns="average",
                width = unit(10, "in"), height = unit(2*10/ratio, "in"),
                show_heatmap_legend=FALSE, na_col = "white", col=cols)
  # col=cols, cluster_rows=F, cluster_columns=F, clustering_distance_rows="spearman", clustering_method_rows="ward.D", "ward.D2", "single", "complete", "average"
  ht <- draw(ht)
  w <- as.numeric(ComplexHeatmap:::width(ht)) / 25.2 # mm to in.
  h <- as.numeric(ComplexHeatmap:::height(ht)) / 25.2 # mm to in.
  
  fname <- paste0(folder, "fits/", name, "IIa_MEG2A_heat_", ncomp, ".pdf")
  # pdf(file = fname, width = w, height = h)
  ht
  dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
  dev.off()
  
  
  #---------Predictions and Assesment-------------#
  newx <- as.matrix(testing[, as.character(1:ncomp)])
  
  #--Number (and indices) of non-zero coefficients
  prediction.non0 <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "nonzero") # length of list elements are numbers of coefficients at lambda.1se!
  
  #--Predicted classes
  prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
  prediction <- factor(prediction, levels = levels(testing$use_celltype_tissue))
  
  #--Predicted probabilities
  prediction.prob <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "response") # array n.cells x n.factors x 1 
  prediction.prob <- as.matrix(prediction.prob[, 1:length(outcomes), 1]) %>% as.data.frame() %>% as_tibble()
  # assess if the prediction probabilities are high enough
  prediction.prob <- prediction.prob %>% 
    add_column(prediction = prediction) %>% 
    dplyr::mutate(max.prob = pmax(macrophage_Marrow,macrophage_Spleen,macrophage_Lung,macrophage_Kidney,macrophage_Mammary_Gland,macrophage_Limb_Muscle)) %>% 
    dplyr::mutate(prediction_alt = if_else(max.prob > 0.55, as.character(prediction), "no_confidence"))
  # see how many cells have < 0.55 pred. prob.
  prediction.prob %>% dplyr::filter(max.prob < 0.55)
  # can try to use this alternative predicted classes with "no_confidence" class in confusion statistics
  # prediction_alt <- pull(prediction.prob, prediction_alt)
  
  # library("mixOmics") # alt. with more options: "measures" package
  confussion <- get.confusion_matrix(truth = testing$use_celltype_tissue, predicted = prediction)
  ave.ber <- get.BER(confussion)
  ave.ber
  
  pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
  confussion <- as.data.frame(confussion) %>% rownames_to_column()
  colnames(confussion) <- c("(Truth)", pred.names)
  print(confussion) # also kable it !
  
  fname <- paste0(folder, "fits/", name, "IIa_MEG2A_confussion_", ncomp, ".csv")
  write_excel_csv(confussion, path = fname)
  
  # confussion_alt <- get.confusion_matrix(truth = testing$use_celltype, predicted = prediction_alt)
  
  # library("caret")
  confussion.c <- confusionMatrix(data=prediction, reference=testing$use_celltype_tissue, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
  confussion.c$overall[1:2] %>% round(digits = 3)
  confussionByClass <-confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3)
  # Balanced Accuracy is derived from Sensitivity & Specificity
  # F1 is derived form Precision & Recall
  
  # confussion_alt.c <- confusionMatrix(data=prediction_alt, reference=testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") # or "prec_recall", or "everything"
  
  # Tuck-in some of these stats in A_stats list, save and plot later
  # within them: nbetas, ave.sparsity, ave.ber, ave.kappa, ave.bacc, bacc
  #-nbetas:
  nbeta <- purrr::map_depth(prediction.non0, 2, .f=length) %>% unlist() %>% .[1:length(outcomes)]
  j <- ncomp/100
  A_stats$IIa$nbetas[j, ] <- nbeta
  #-ave.sparsity:
  A_stats$IIa$ave.sparsity[j] <- ave.spars
  #-ave.ber:
  A_stats$IIa$ave.ber[j] <- ave.ber
  #-ave.kappa:
  A_stats$IIa$ave.kappa[j] <- confussion.c$overall[2]
  #-ave.bacc:
  ave.bac <- confussion.c$overall[c(1,3:4)]
  A_stats$IIa$ave.bacc[j, ] <- ave.bac
  #-bacc:
  A_stats$IIa$bacc[j, ] <- confussionByClass[,6]
  # save this at the end of all I-VI analyses !
  
}


# .d
for (ncomp in c(100, 200, 300)) {
  # import a. if not loaded
  a <- get(paste0("a.", ncomp, ".d"))
  testing  <- a[-inTrain.tissue, ] %>% dplyr::filter(str_detect(use_celltype_tissue, "macrophage_")) %>% dplyr::mutate(use_celltype_tissue = fct_drop(use_celltype_tissue))
  # import cv.fit
  fname <- paste0(folder, "fits/", name, "IIa_MEG2A_fitCV_", ncomp, "_D", ".rds")
  fit.cv <- readRDS(file = fname)
  
  # use fuctions: plot(), print(), coef(), predict()
  
  
  #---All betas and analyses are at lambda.1se (lambda.min would result in less sparsity)
  betas <- coef(fit.cv, s = "lambda.1se") # a list
  outcomes <- names(betas)
  
  betas <- purrr::map(betas, .f = as.matrix)
  betas <- do.call(cbind, betas)
  colnames(betas) <- outcomes
  
  betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble()
  betas <- betas %>% dplyr::filter_at(vars(-compID), any_vars(.!=0)) # disregard "other only" != 0 coefficients
  print(betas, n=ncomp) # also kable it !
  ave.spars <- nrow(betas)-1
  
  fname <- paste0(folder, "fits/", name, "IIa_MEG2A_betas_", ncomp, "_D", ".csv")
  write_excel_csv(betas, path = fname)
  
  # make heatplot for easier assesment of coefficients
  # library("ComplexHeatmap")
  b.mat <- as.matrix(betas[2:nrow(betas), 2:ncol(betas)])
  rownames(b.mat) <- betas$compID[2:nrow(betas)]
  colnames(b.mat) <- colnames(b.mat) %>% str_replace_all("_", " ")
  b.mat <- t(b.mat)
  
  # maybe exclude weak coefficients from visualisation?.. 0.5sd away from zero?..
  plot(density(b.mat))
  cutof <- sd(b.mat)/2
  
  # library("circlize")
  col_fun = colorRamp2(breaks = c(min(b.mat)+2, -cutof, cutof, max(b.mat)-2), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
  cols <- col_fun(seq(min(b.mat), max(b.mat), by = 2))
  ratio <- ncol(b.mat) / nrow(b.mat)
  
  #---first see how to exclude small values +/- ~ 0.5 or 0.25 sd and less !
  ht <- Heatmap(b.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
                row_names_max_width = max_text_width(rownames(b.mat)),
                column_dend_height = unit(0.7, "cm"), row_dend_width = unit(0.8, "cm"), name = "",
                clustering_distance_rows="spearman", clustering_distance_columns="spearman",
                clustering_method_rows="average", clustering_method_columns="average",
                width = unit(10, "in"), height = unit(2*10/ratio, "in"),
                show_heatmap_legend=FALSE, na_col = "white", col=cols)
  # col=cols, cluster_rows=F, cluster_columns=F, clustering_distance_rows="spearman", clustering_method_rows="ward.D", "ward.D2", "single", "complete", "average"
  ht <- draw(ht)
  w <- as.numeric(ComplexHeatmap:::width(ht)) / 25.2 # mm to in.
  h <- as.numeric(ComplexHeatmap:::height(ht)) / 25.2 # mm to in.
  
  fname <- paste0(folder, "fits/", name, "IIa_MEG2A_heat_", ncomp, "_D", ".pdf")
  # pdf(file = fname, width = w, height = h)
  ht
  dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
  dev.off()
  
  
  #---------Predictions and Assesment-------------#
  newx <- as.matrix(testing[, as.character(1:ncomp)])
  
  #--Number (and indices) of non-zero coefficients
  prediction.non0 <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "nonzero") # length of list elements are numbers of coefficients at lambda.1se!
  
  #--Predicted classes
  prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
  prediction <- factor(prediction, levels = levels(testing$use_celltype_tissue))
  
  #--Predicted probabilities
  prediction.prob <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "response") # array n.cells x n.factors x 1 
  prediction.prob <- as.matrix(prediction.prob[, 1:length(outcomes), 1]) %>% as.data.frame() %>% as_tibble()
  # assess if the prediction probabilities are high enough
  prediction.prob <- prediction.prob %>% 
    add_column(prediction = prediction) %>% 
    dplyr::mutate(max.prob = pmax(macrophage_Marrow,macrophage_Spleen,macrophage_Lung,macrophage_Kidney,macrophage_Mammary_Gland,macrophage_Limb_Muscle)) %>% 
    dplyr::mutate(prediction_alt = if_else(max.prob > 0.55, as.character(prediction), "no_confidence"))
  # see how many cells have < 0.55 pred. prob.
  prediction.prob %>% dplyr::filter(max.prob < 0.55)
  # can try to use this alternative predicted classes with "no_confidence" class in confusion statistics
  # prediction_alt <- pull(prediction.prob, prediction_alt)
  
  # library("mixOmics") # alt. with more options: "measures" package
  confussion <- get.confusion_matrix(truth = testing$use_celltype_tissue, predicted = prediction)
  ave.ber <- get.BER(confussion)
  ave.ber
  
  pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
  confussion <- as.data.frame(confussion) %>% rownames_to_column()
  colnames(confussion) <- c("(Truth)", pred.names)
  print(confussion) # also kable it !
  
  fname <- paste0(folder, "fits/", name, "IIa_MEG2A_confussion_", ncomp, "_D", ".csv")
  write_excel_csv(confussion, path = fname)
  
  # confussion_alt <- get.confusion_matrix(truth = testing$use_celltype, predicted = prediction_alt)
  
  # library("caret")
  confussion.c <- confusionMatrix(data=prediction, reference=testing$use_celltype_tissue, dnn = c("Prediction", "Reference"), mode = "sens_spec")
  confussion.c$overall[1:2] %>% round(digits = 3)
  confussionByClass <-confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3)
  # Balanced Accuracy is derived from Sensitivity & Specificity
  # F1 is derived form Precision & Recall
  
  # confussion_alt.c <- confusionMatrix(data=prediction_alt, reference=testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") # or "prec_recall", or "everything"
  
  # Tuck-in some of these stats in A_stats list, save and plot later
  # within them: nbetas, ave.sparsity, ave.ber, ave.kappa, ave.bacc, bacc
  #-nbetas:
  nbeta <- purrr::map_depth(prediction.non0, 2, .f=length) %>% unlist() %>% .[1:length(outcomes)]
  j <- ncomp/100
  A_stats$IIa_d$nbetas[j, ] <- nbeta
  #-ave.sparsity:
  A_stats$IIa_d$ave.sparsity[j] <- ave.spars
  #-ave.ber:
  A_stats$IIa_d$ave.ber[j] <- ave.ber
  #-ave.kappa:
  A_stats$IIa_d$ave.kappa[j] <- confussion.c$overall[2]
  #-ave.bacc:
  ave.bac <- confussion.c$overall[c(1,3:4)]
  A_stats$IIa_d$ave.bacc[j, ] <- ave.bac
  #-bacc:
  A_stats$IIa_d$bacc[j, ] <- confussionByClass[,6]
  # save this at the end of all I-VI analyses !
  
}



### III analysis: (MLGExp)
#-----------------------------------------------------------------------------#

# for III IV
outcomes <- c("granulocyte", "promonocyte", "monocyte", "macrophage", "dendritic_cell", "mast_cell", "other")
#-nbetas:
A_stats$III$nbetas <- matrix(0, nrow = 1, ncol = length(outcomes)-1)
colnames(A_stats$III$nbetas) <- outcomes[1:length(outcomes)-1]
#-ave.sparsity:
A_stats$III$ave.sparsity <- vector(mode = "numeric", length = 1)
#-ave.ber:
A_stats$III$ave.ber <- vector(mode = "numeric", length = 1)
#-ave.kappa:
A_stats$III$ave.kappa <- vector(mode = "numeric", length = 1)
#-ave.bacc:
A_stats$III$ave.bacc <- matrix(0, nrow = 1, ncol = 3)
colnames(A_stats$III$ave.bacc) <- c("Accuracy", "loCI", "hiCI")
#-bacc:
A_stats$III$bacc <- matrix(0, nrow = 1, ncol = length(outcomes))
colnames(A_stats$III$bacc) <- outcomes[1:length(outcomes)]



# import exp.df if not loaded
exp.testing  <- exp.df[-inTrain.tissue, ]
# import cv.fit
fname <- paste0(folder, "fits/", name, "III_MLGExp_fitCV", ".rds")
fit.cv <- readRDS(file = fname)

# use fuctions: plot(), print(), coef(), predict()


#---All betas and analyses are at lambda.1se (lambda.min would result in less sparsity)
betas <- coef(fit.cv, s = "lambda.1se") # a list
outcomes <- names(betas)

betas <- purrr::map(betas, .f = as.matrix)
betas <- do.call(cbind, betas)
colnames(betas) <- outcomes

betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble()
betas <- betas %>% dplyr::filter_at(vars(-compID, -other), any_vars(.!=0)) # disregard "other only" != 0 coefficients
print(betas, n=250) # also kable it !
ave.spars <- nrow(betas)-1

fname <- paste0(folder, "fits/", name, "III_MLGExp_betas", ".csv")
write_excel_csv(betas, path = fname)

# make heatplot for easier assesment of coefficients
# library("ComplexHeatmap")
b.mat <- as.matrix(betas[2:nrow(betas), 2:ncol(betas)])
rownames(b.mat) <- betas$compID[2:nrow(betas)]
colnames(b.mat) <- colnames(b.mat) %>% str_replace_all("_", " ")
b.mat <- t(b.mat)

# maybe exclude weak coefficients from visualisation?.. 0.5sd away from zero?..
plot(density(b.mat))
cutof <- sd(b.mat)/8

# library("circlize")
col_fun = colorRamp2(breaks = c(min(b.mat), -cutof, cutof, max(b.mat)), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
cols <- col_fun(seq(min(b.mat), max(b.mat), by = 0.002))
ratio <- ncol(b.mat) / nrow(b.mat)

#---first see how to exclude small values +/- ~ 0.5 or 0.25 sd and less !
ht <- Heatmap(b.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(b.mat)),
              column_dend_height = unit(0.7, "cm"), row_dend_width = unit(0.8, "cm"), name = "",
              clustering_distance_rows="spearman", clustering_distance_columns="spearman",
              clustering_method_rows="average", clustering_method_columns="average",
              width = unit(10, "in"), height = unit(2*10/ratio, "in"),
              show_heatmap_legend=FALSE, na_col = "white", col=cols)
# col=cols, cluster_rows=F, cluster_columns=F, clustering_distance_rows="spearman", clustering_method_rows="ward.D", "ward.D2", "single", "complete", "average"
ht <- draw(ht)
w <- as.numeric(ComplexHeatmap:::width(ht)) / 25.2 # mm to in.
h <- as.numeric(ComplexHeatmap:::height(ht)) / 25.2 # mm to in.

fname <- paste0(folder, "fits/", name, "III_MLGExp_heat", ".pdf")
# pdf(file = fname, width = w, height = h)
ht
dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
dev.off()


#---------Predictions and Assesment-------------#
newx <- as.matrix(exp.testing[, hvg])

#--Number (and indices) of non-zero coefficients
prediction.non0 <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "nonzero") # length of list elements are numbers of coefficients at lambda.1se!

#--Predicted classes
prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
prediction <- factor(prediction, levels = levels(exp.testing$use_celltype))

#--Predicted probabilities
prediction.prob <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "response") # array n.cells x n.factors x 1 
prediction.prob <- as.matrix(prediction.prob[, 1:length(outcomes), 1]) %>% as.data.frame() %>% as_tibble()
# assess if the prediction probabilities are high enough
prediction.prob <- prediction.prob %>% 
  add_column(prediction = prediction) %>% 
  dplyr::mutate(max.prob = pmax(granulocyte,promonocyte,monocyte,macrophage,dendritic_cell,mast_cell,other)) %>% 
  dplyr::mutate(prediction_alt = if_else(max.prob > 0.55, as.character(prediction), "no_confidence"))
# see how many cells have < 0.55 pred. prob.
print(prediction.prob %>% dplyr::filter(max.prob < 0.55))
# can try to use this alternative predicted classes with "no_confidence" class in confusion statistics
# prediction_alt <- pull(prediction.prob, prediction_alt)

# library("mixOmics") # alt. with more options: "measures" package
confussion <- get.confusion_matrix(truth = exp.testing$use_celltype, predicted = prediction)
ave.ber <- get.BER(confussion)
ave.ber

pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
confussion <- as.data.frame(confussion) %>% rownames_to_column()
colnames(confussion) <- c("(Truth)", pred.names)
print(confussion) # also kable it !

fname <- paste0(folder, "fits/", name, "III_MLGExp_confussion", ".csv")
write_excel_csv(confussion, path = fname)

# confussion_alt <- get.confusion_matrix(truth = testing$use_celltype, predicted = prediction_alt)

# library("caret")
confussion.c <- confusionMatrix(data=prediction, reference=exp.testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
confussion.c$overall[1:2] %>% round(digits = 3)
confussionByClass <-confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3)
# Balanced Accuracy is derived from Sensitivity & Specificity
# F1 is derived form Precision & Recall

# confussion_alt.c <- confusionMatrix(data=prediction_alt, reference=testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") # or "prec_recall", or "everything"

# Tuck-in some of these stats in A_stats list, save and plot later
# within them: nbetas, ave.sparsity, ave.ber, ave.kappa, ave.bacc, bacc
#-nbetas:
nbeta <- purrr::map_depth(prediction.non0, 2, .f=length) %>% unlist() %>% .[1:length(outcomes)-1]
A_stats$III$nbetas[1, ] <- nbeta
#-ave.sparsity:
A_stats$III$ave.sparsity[1] <- ave.spars
#-ave.ber:
A_stats$III$ave.ber[1] <- ave.ber
#-ave.kappa:
A_stats$III$ave.kappa[1] <- confussion.c$overall[2]
#-ave.bacc:
ave.bac <- confussion.c$overall[c(1,3:4)]
A_stats$III$ave.bacc[1, ] <- ave.bac
#-bacc:
A_stats$III$bacc[1, ] <- confussionByClass[,6]
# save this at the end of all I-VI analyses !



### IIIa analysis: (MLG2Exp)
#-----------------------------------------------------------------------------#

# for IIIa IVa
outcomes <- c("macrophage_Marrow", "macrophage_Spleen", "macrophage_Lung", "macrophage_Kidney", "macrophage_Mammary_Gland", "macrophage_Limb_Muscle")
#-nbetas:
A_stats$IIIa$nbetas <- matrix(0, nrow = 1, ncol = length(outcomes))
colnames(A_stats$IIIa$nbetas) <- outcomes[1:length(outcomes)]
#-ave.sparsity:
A_stats$IIIa$ave.sparsity <- vector(mode = "numeric", length = 1)
#-ave.ber:
A_stats$IIIa$ave.ber <- vector(mode = "numeric", length = 1)
#-ave.kappa:
A_stats$IIIa$ave.kappa <- vector(mode = "numeric", length = 1)
#-ave.bacc:
A_stats$IIIa$ave.bacc <- matrix(0, nrow = 1, ncol = 3)
colnames(A_stats$IIIa$ave.bacc) <- c("Accuracy", "loCI", "hiCI")
#-bacc:
A_stats$IIIa$bacc <- matrix(0, nrow = 1, ncol = length(outcomes))
colnames(A_stats$IIIa$bacc) <- outcomes[1:length(outcomes)]



# import exp.df if not loaded
exp.testing  <- exp.df[-inTrain.tissue, ] %>% dplyr::filter(str_detect(use_celltype_tissue, "macrophage_")) %>% dplyr::mutate(use_celltype_tissue = fct_drop(use_celltype_tissue))
# import cv.fit
fname <- paste0(folder, "fits/", name, "IIIa_MLG2Exp_fitCV", ".rds")
fit.cv <- readRDS(file = fname)

# use fuctions: plot(), print(), coef(), predict()


#---All betas and analyses are at lambda.1se (lambda.min would result in less sparsity)
betas <- coef(fit.cv, s = "lambda.1se") # a list
outcomes <- names(betas)

betas <- purrr::map(betas, .f = as.matrix)
betas <- do.call(cbind, betas)
colnames(betas) <- outcomes

betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble()
betas <- betas %>% dplyr::filter_at(vars(-compID), any_vars(.!=0)) # disregard "other only" != 0 coefficients
print(betas, n=250) # also kable it !
ave.spars <- nrow(betas)-1

fname <- paste0(folder, "fits/", name, "IIIa_MLG2Exp_betas", ".csv")
write_excel_csv(betas, path = fname)

# make heatplot for easier assesment of coefficients
# library("ComplexHeatmap")
b.mat <- as.matrix(betas[2:nrow(betas), 2:ncol(betas)])
rownames(b.mat) <- betas$compID[2:nrow(betas)]
colnames(b.mat) <- colnames(b.mat) %>% str_replace_all("_", " ")
b.mat <- t(b.mat)

# maybe exclude weak coefficients from visualisation?.. 0.5sd away from zero?..
plot(density(b.mat))
cutof <- sd(b.mat)/8

# library("circlize")
col_fun = colorRamp2(breaks = c(min(b.mat), -cutof, cutof, max(b.mat)), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
cols <- col_fun(seq(min(b.mat), max(b.mat), by = 0.05))
ratio <- ncol(b.mat) / nrow(b.mat)

#---first see how to exclude small values +/- ~ 0.5 or 0.25 sd and less !
ht <- Heatmap(b.mat, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 5),
              row_names_max_width = max_text_width(rownames(b.mat)),
              column_dend_height = unit(0.4, "cm"), row_dend_width = unit(0.4, "cm"), name = "",
              clustering_distance_rows="spearman", clustering_distance_columns="spearman",
              clustering_method_rows="average", clustering_method_columns="average",
              width = unit(14, "in"), height = unit(2*14/ratio, "in"),
              show_heatmap_legend=FALSE, na_col = "white", col=cols)
# col=cols, cluster_rows=F, cluster_columns=F, clustering_distance_rows="spearman", clustering_method_rows="ward.D", "ward.D2", "single", "complete", "average"
ht <- draw(ht)
w <- as.numeric(ComplexHeatmap:::width(ht)) / 25.2 # mm to in.
h <- as.numeric(ComplexHeatmap:::height(ht)) / 25.2 # mm to in.

fname <- paste0(folder, "fits/", name, "IIIa_MLG2Exp_heat", ".pdf")
# pdf(file = fname, width = w, height = h)
ht
dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
dev.off()


#---------Predictions and Assesment-------------#
newx <- as.matrix(exp.testing[, hvg])

#--Number (and indices) of non-zero coefficients
prediction.non0 <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "nonzero") # length of list elements are numbers of coefficients at lambda.1se!

#--Predicted classes
prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
prediction <- factor(prediction, levels = levels(exp.testing$use_celltype_tissue))

#--Predicted probabilities
prediction.prob <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "response") # array n.cells x n.factors x 1 
prediction.prob <- as.matrix(prediction.prob[, 1:length(outcomes), 1]) %>% as.data.frame() %>% as_tibble()
# assess if the prediction probabilities are high enough
prediction.prob <- prediction.prob %>% 
  add_column(prediction = prediction) %>% 
  dplyr::mutate(max.prob = pmax(macrophage_Marrow,macrophage_Spleen,macrophage_Lung,macrophage_Kidney,macrophage_Mammary_Gland,macrophage_Limb_Muscle)) %>% 
  dplyr::mutate(prediction_alt = if_else(max.prob > 0.55, as.character(prediction), "no_confidence"))
# see how many cells have < 0.55 pred. prob.
print(prediction.prob %>% dplyr::filter(max.prob < 0.55))
# can try to use this alternative predicted classes with "no_confidence" class in confusion statistics
# prediction_alt <- pull(prediction.prob, prediction_alt)

# library("mixOmics") # alt. with more options: "measures" package
confussion <- get.confusion_matrix(truth = exp.testing$use_celltype_tissue, predicted = prediction)
ave.ber <- get.BER(confussion)
ave.ber

pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
confussion <- as.data.frame(confussion) %>% rownames_to_column()
colnames(confussion) <- c("(Truth)", pred.names)
print(confussion) # also kable it !

fname <- paste0(folder, "fits/", name, "IIIa_MLG2Exp_confussion", ".csv")
write_excel_csv(confussion, path = fname)

# confussion_alt <- get.confusion_matrix(truth = testing$use_celltype, predicted = prediction_alt)

# library("caret")
confussion.c <- confusionMatrix(data=prediction, reference=exp.testing$use_celltype_tissue, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
confussion.c$overall[1:2] %>% round(digits = 3)
confussionByClass <-confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3)
# Balanced Accuracy is derived from Sensitivity & Specificity
# F1 is derived form Precision & Recall

# confussion_alt.c <- confusionMatrix(data=prediction_alt, reference=testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") # or "prec_recall", or "everything"

# Tuck-in some of these stats in A_stats list, save and plot later
# within them: nbetas, ave.sparsity, ave.ber, ave.kappa, ave.bacc, bacc
#-nbetas:
nbeta <- purrr::map_depth(prediction.non0, 2, .f=length) %>% unlist() %>% .[1:length(outcomes)]
A_stats$IIIa$nbetas[1, ] <- nbeta
#-ave.sparsity:
A_stats$IIIa$ave.sparsity[1] <- ave.spars
#-ave.ber:
A_stats$IIIa$ave.ber[1] <- ave.ber
#-ave.kappa:
A_stats$IIIa$ave.kappa[1] <- confussion.c$overall[2]
#-ave.bacc:
ave.bac <- confussion.c$overall[c(1,3:4)]
A_stats$IIIa$ave.bacc[1, ] <- ave.bac
#-bacc:
A_stats$IIIa$bacc[1, ] <- confussionByClass[,6]
# save this at the end of all I-VI analyses !




### IV analysis: (MEGExp)
#-----------------------------------------------------------------------------#

# having seen extreme sparsness of the above MLGExp model, it would have been better to set alpha = 0.925 or 0.9 during fit.cv
# for III IV
outcomes <- c("granulocyte", "promonocyte", "monocyte", "macrophage", "dendritic_cell", "mast_cell", "other")
#-nbetas:
A_stats$IV$nbetas <- matrix(0, nrow = 1, ncol = length(outcomes)-1)
colnames(A_stats$IV$nbetas) <- outcomes[1:length(outcomes)-1]
#-ave.sparsity:
A_stats$IV$ave.sparsity <- vector(mode = "numeric", length = 1)
#-ave.ber:
A_stats$IV$ave.ber <- vector(mode = "numeric", length = 1)
#-ave.kappa:
A_stats$IV$ave.kappa <- vector(mode = "numeric", length = 1)
#-ave.bacc:
A_stats$IV$ave.bacc <- matrix(0, nrow = 1, ncol = 3)
colnames(A_stats$IV$ave.bacc) <- c("Accuracy", "loCI", "hiCI")
#-bacc:
A_stats$IV$bacc <- matrix(0, nrow = 1, ncol = length(outcomes))
colnames(A_stats$IV$bacc) <- outcomes[1:length(outcomes)]



# import exp.df if not loaded
exp.testing  <- exp.df[-inTrain.tissue, ]
# import cv.fit
fname <- paste0(folder, "fits/", name, "IV_MEGExp_fitCV", ".rds")
fit.cv <- readRDS(file = fname)

# use fuctions: plot(), print(), coef(), predict()


#---All betas and analyses are at lambda.1se (lambda.min would result in less sparsity)
betas <- coef(fit.cv, s = "lambda.1se") # a list
outcomes <- names(betas)

betas <- purrr::map(betas, .f = as.matrix)
betas <- do.call(cbind, betas)
colnames(betas) <- outcomes

betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble()
betas <- betas %>% dplyr::filter_at(vars(-compID, -other), any_vars(.!=0)) # disregard "other only" != 0 coefficients
print(betas, n=250) # also kable it !
ave.spars <- nrow(betas)-1

fname <- paste0(folder, "fits/", name, "IV_MEGExp_betas", ".csv")
write_excel_csv(betas, path = fname)

# make heatplot for easier assesment of coefficients
# library("ComplexHeatmap")
b.mat <- as.matrix(betas[2:nrow(betas), 2:ncol(betas)])
rownames(b.mat) <- betas$compID[2:nrow(betas)]
colnames(b.mat) <- colnames(b.mat) %>% str_replace_all("_", " ")
b.mat <- t(b.mat)

# maybe exclude weak coefficients from visualisation?.. 0.5sd away from zero?..
plot(density(b.mat))
cutof <- sd(b.mat)/8

# library("circlize")
col_fun = colorRamp2(breaks = c(min(b.mat), -cutof, cutof, max(b.mat)), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
cols <- col_fun(seq(min(b.mat), max(b.mat), by = 0.002))
ratio <- ncol(b.mat) / nrow(b.mat)

#---first see how to exclude small values +/- ~ 0.5 or 0.25 sd and less !
ht <- Heatmap(b.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(b.mat)),
              column_dend_height = unit(0.7, "cm"), row_dend_width = unit(0.8, "cm"), name = "",
              clustering_distance_rows="spearman", clustering_distance_columns="spearman",
              clustering_method_rows="average", clustering_method_columns="average",
              width = unit(10, "in"), height = unit(2*10/ratio, "in"),
              show_heatmap_legend=FALSE, na_col = "white", col=cols)
# col=cols, cluster_rows=F, cluster_columns=F, clustering_distance_rows="spearman", clustering_method_rows="ward.D", "ward.D2", "single", "complete", "average"
ht <- draw(ht)
w <- as.numeric(ComplexHeatmap:::width(ht)) / 25.2 # mm to in.
h <- as.numeric(ComplexHeatmap:::height(ht)) / 25.2 # mm to in.

fname <- paste0(folder, "fits/", name, "IV_MEGExp_heat", ".pdf")
# pdf(file = fname, width = w, height = h)
ht
dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
dev.off()


#---------Predictions and Assesment-------------#
newx <- as.matrix(exp.testing[, hvg])

#--Number (and indices) of non-zero coefficients
prediction.non0 <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "nonzero") # length of list elements are numbers of coefficients at lambda.1se!

#--Predicted classes
prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
prediction <- factor(prediction, levels = levels(exp.testing$use_celltype))

#--Predicted probabilities
prediction.prob <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "response") # array n.cells x n.factors x 1 
prediction.prob <- as.matrix(prediction.prob[, 1:length(outcomes), 1]) %>% as.data.frame() %>% as_tibble()
# assess if the prediction probabilities are high enough
prediction.prob <- prediction.prob %>% 
  add_column(prediction = prediction) %>% 
  dplyr::mutate(max.prob = pmax(granulocyte,promonocyte,monocyte,macrophage,dendritic_cell,mast_cell,other)) %>% 
  dplyr::mutate(prediction_alt = if_else(max.prob > 0.55, as.character(prediction), "no_confidence"))
# see how many cells have < 0.55 pred. prob.
print(prediction.prob %>% dplyr::filter(max.prob < 0.55))
# can try to use this alternative predicted classes with "no_confidence" class in confusion statistics
# prediction_alt <- pull(prediction.prob, prediction_alt)

# library("mixOmics") # alt. with more options: "measures" package
confussion <- get.confusion_matrix(truth = exp.testing$use_celltype, predicted = prediction)
ave.ber <- get.BER(confussion)
ave.ber

pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
confussion <- as.data.frame(confussion) %>% rownames_to_column()
colnames(confussion) <- c("(Truth)", pred.names)
print(confussion) # also kable it !

fname <- paste0(folder, "fits/", name, "IV_MEGExp_confussion", ".csv")
write_excel_csv(confussion, path = fname)

# confussion_alt <- get.confusion_matrix(truth = testing$use_celltype, predicted = prediction_alt)

# library("caret")
confussion.c <- confusionMatrix(data=prediction, reference=exp.testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
confussion.c$overall[1:2] %>% round(digits = 3)
confussionByClass <-confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3)
# Balanced Accuracy is derived from Sensitivity & Specificity
# F1 is derived form Precision & Recall

# confussion_alt.c <- confusionMatrix(data=prediction_alt, reference=testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") # or "prec_recall", or "everything"

# Tuck-in some of these stats in A_stats list, save and plot later
# within them: nbetas, ave.sparsity, ave.ber, ave.kappa, ave.bacc, bacc
#-nbetas:
nbeta <- purrr::map_depth(prediction.non0, 2, .f=length) %>% unlist() %>% .[1:length(outcomes)-1]
A_stats$IV$nbetas[1, ] <- nbeta
#-ave.sparsity:
A_stats$IV$ave.sparsity[1] <- ave.spars
#-ave.ber:
A_stats$IV$ave.ber[1] <- ave.ber
#-ave.kappa:
A_stats$IV$ave.kappa[1] <- confussion.c$overall[2]
#-ave.bacc:
ave.bac <- confussion.c$overall[c(1,3:4)]
A_stats$IV$ave.bacc[1, ] <- ave.bac
#-bacc:
A_stats$IV$bacc[1, ] <- confussionByClass[,6]
# save this at the end of all I-VI analyses !



### IVa analysis: (MEG2Exp)
#-----------------------------------------------------------------------------#

# for IIIa IVa
outcomes <- c("macrophage_Marrow", "macrophage_Spleen", "macrophage_Lung", "macrophage_Kidney", "macrophage_Mammary_Gland", "macrophage_Limb_Muscle")
#-nbetas:
A_stats$IVa$nbetas <- matrix(0, nrow = 1, ncol = length(outcomes))
colnames(A_stats$IVa$nbetas) <- outcomes[1:length(outcomes)]
#-ave.sparsity:
A_stats$IVa$ave.sparsity <- vector(mode = "numeric", length = 1)
#-ave.ber:
A_stats$IVa$ave.ber <- vector(mode = "numeric", length = 1)
#-ave.kappa:
A_stats$IVa$ave.kappa <- vector(mode = "numeric", length = 1)
#-ave.bacc:
A_stats$IVa$ave.bacc <- matrix(0, nrow = 1, ncol = 3)
colnames(A_stats$IVa$ave.bacc) <- c("Accuracy", "loCI", "hiCI")
#-bacc:
A_stats$IVa$bacc <- matrix(0, nrow = 1, ncol = length(outcomes))
colnames(A_stats$IVa$bacc) <- outcomes[1:length(outcomes)]



# import exp.df if not loaded
exp.testing  <- exp.df[-inTrain.tissue, ] %>% dplyr::filter(str_detect(use_celltype_tissue, "macrophage_")) %>% dplyr::mutate(use_celltype_tissue = fct_drop(use_celltype_tissue))
# import cv.fit
fname <- paste0(folder, "fits/", name, "IVa_MEG2Exp_fitCV", ".rds")
fit.cv <- readRDS(file = fname)

# use fuctions: plot(), print(), coef(), predict()


#---All betas and analyses are at lambda.1se (lambda.min would result in less sparsity)
betas <- coef(fit.cv, s = "lambda.1se") # a list
outcomes <- names(betas)

betas <- purrr::map(betas, .f = as.matrix)
betas <- do.call(cbind, betas)
colnames(betas) <- outcomes

betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble()
betas <- betas %>% dplyr::filter_at(vars(-compID), any_vars(.!=0)) # disregard "other only" != 0 coefficients
print(betas, n=250) # also kable it !
ave.spars <- nrow(betas)-1

fname <- paste0(folder, "fits/", name, "IVa_MEG2Exp_betas", ".csv")
write_excel_csv(betas, path = fname)

# make heatplot for easier assesment of coefficients
# library("ComplexHeatmap")
b.mat <- as.matrix(betas[2:nrow(betas), 2:ncol(betas)])
rownames(b.mat) <- betas$compID[2:nrow(betas)]
colnames(b.mat) <- colnames(b.mat) %>% str_replace_all("_", " ")
b.mat <- t(b.mat)

# maybe exclude weak coefficients from visualisation?.. 0.5sd away from zero?..
plot(density(b.mat))
cutof <- sd(b.mat)/8

# library("circlize")
col_fun = colorRamp2(breaks = c(min(b.mat), -cutof, cutof, max(b.mat)), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
cols <- col_fun(seq(min(b.mat), max(b.mat), by = 0.05))
ratio <- ncol(b.mat) / nrow(b.mat)

#---first see how to exclude small values +/- ~ 0.5 or 0.25 sd and less !
ht <- Heatmap(b.mat, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 5),
              row_names_max_width = max_text_width(rownames(b.mat)),
              column_dend_height = unit(0.4, "cm"), row_dend_width = unit(0.4, "cm"), name = "",
              clustering_distance_rows="spearman", clustering_distance_columns="spearman",
              clustering_method_rows="average", clustering_method_columns="average",
              width = unit(15, "in"), height = unit(2*15/ratio, "in"),
              show_heatmap_legend=FALSE, na_col = "white", col=cols)
# col=cols, cluster_rows=F, cluster_columns=F, clustering_distance_rows="spearman", clustering_method_rows="ward.D", "ward.D2", "single", "complete", "average"
ht <- draw(ht)
w <- as.numeric(ComplexHeatmap:::width(ht)) / 25.2 # mm to in.
h <- as.numeric(ComplexHeatmap:::height(ht)) / 25.2 # mm to in.

fname <- paste0(folder, "fits/", name, "IVa_MEG2Exp_heat", ".pdf")
# pdf(file = fname, width = w, height = h)
ht
dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
dev.off()


#---------Predictions and Assesment-------------#
newx <- as.matrix(exp.testing[, hvg])

#--Number (and indices) of non-zero coefficients
prediction.non0 <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "nonzero") # length of list elements are numbers of coefficients at lambda.1se!

#--Predicted classes
prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
prediction <- factor(prediction, levels = levels(exp.testing$use_celltype_tissue))

#--Predicted probabilities
prediction.prob <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "response") # array n.cells x n.factors x 1 
prediction.prob <- as.matrix(prediction.prob[, 1:length(outcomes), 1]) %>% as.data.frame() %>% as_tibble()
# assess if the prediction probabilities are high enough
prediction.prob <- prediction.prob %>% 
  add_column(prediction = prediction) %>% 
  dplyr::mutate(max.prob = pmax(macrophage_Marrow,macrophage_Spleen,macrophage_Lung,macrophage_Kidney,macrophage_Mammary_Gland,macrophage_Limb_Muscle)) %>% 
  dplyr::mutate(prediction_alt = if_else(max.prob > 0.55, as.character(prediction), "no_confidence"))
# see how many cells have < 0.55 pred. prob.
print(prediction.prob %>% dplyr::filter(max.prob < 0.55))
# can try to use this alternative predicted classes with "no_confidence" class in confusion statistics
# prediction_alt <- pull(prediction.prob, prediction_alt)

# library("mixOmics") # alt. with more options: "measures" package
confussion <- get.confusion_matrix(truth = exp.testing$use_celltype_tissue, predicted = prediction)
ave.ber <- get.BER(confussion)
ave.ber

pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
confussion <- as.data.frame(confussion) %>% rownames_to_column()
colnames(confussion) <- c("(Truth)", pred.names)
print(confussion) # also kable it !

fname <- paste0(folder, "fits/", name, "IVa_MEG2Exp_confussion", ".csv")
write_excel_csv(confussion, path = fname)

# confussion_alt <- get.confusion_matrix(truth = testing$use_celltype, predicted = prediction_alt)

# library("caret")
confussion.c <- confusionMatrix(data=prediction, reference=exp.testing$use_celltype_tissue, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
confussion.c$overall[1:2] %>% round(digits = 3)
confussionByClass <-confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3)
# Balanced Accuracy is derived from Sensitivity & Specificity
# F1 is derived form Precision & Recall

# confussion_alt.c <- confusionMatrix(data=prediction_alt, reference=testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") # or "prec_recall", or "everything"

# Tuck-in some of these stats in A_stats list, save and plot later
# within them: nbetas, ave.sparsity, ave.ber, ave.kappa, ave.bacc, bacc
#-nbetas:
nbeta <- purrr::map_depth(prediction.non0, 2, .f=length) %>% unlist() %>% .[1:length(outcomes)]
A_stats$IVa$nbetas[1, ] <- nbeta
#-ave.sparsity:
A_stats$IVa$ave.sparsity[1] <- ave.spars
#-ave.ber:
A_stats$IVa$ave.ber[1] <- ave.ber
#-ave.kappa:
A_stats$IVa$ave.kappa[1] <- confussion.c$overall[2]
#-ave.bacc:
ave.bac <- confussion.c$overall[c(1,3:4)]
A_stats$IVa$ave.bacc[1, ] <- ave.bac
#-bacc:
A_stats$IVa$bacc[1, ] <- confussionByClass[,6]
# save this at the end of all I-VI analyses !


### V analysis: (MLTA)
#-----------------------------------------------------------------------------#


# for V
outcomes <- c("granulocyte_Marrow","promonocyte_Marrow","monocyte_Marrow","monocyte_Lung","macrophage_Marrow","macrophage_Spleen","macrophage_Lung", 
              "macrophage_Kidney","macrophage_Mammary_Gland","macrophage_Limb_Muscle","dendritic_cell_Spleen","mast_cell_Lung","other")
#-nbetas:
A_stats$V$nbetas <- matrix(0, nrow = 3, ncol = length(outcomes)-1)
colnames(A_stats$V$nbetas) <- outcomes[1:length(outcomes)-1]
rownames(A_stats$V$nbetas) <- c(100,200,300)
#-ave.sparsity:
A_stats$V$ave.sparsity <- vector(mode = "numeric", length = 3)
names(A_stats$V$ave.sparsity) <- c(100,200,300)
#-ave.ber:
A_stats$V$ave.ber <- vector(mode = "numeric", length = 3)
names(A_stats$V$ave.ber) <- c(100,200,300)
#-ave.kappa:
A_stats$V$ave.kappa <- vector(mode = "numeric", length = 3)
names(A_stats$V$ave.kappa) <- c(100,200,300)
#-ave.bacc:
A_stats$V$ave.bacc <- matrix(0, nrow = 3, ncol = 3)
colnames(A_stats$V$ave.bacc) <- c("Accuracy", "loCI", "hiCI")
rownames(A_stats$V$ave.bacc) <- c(100,200,300)
#-bacc:
A_stats$V$bacc <- matrix(0, nrow = 3, ncol = length(outcomes))
colnames(A_stats$V$bacc) <- outcomes[1:length(outcomes)]
rownames(A_stats$V$bacc) <- c(100,200,300)
# duplicate for .d analysis:
A_stats$V_d <- A_stats$V


# Saves: coefficients table, coefficients heatmap, confussion table; and builds a list with stats 
for (ncomp in c(100, 200, 300)) {
  # import a. if not loaded
  a <- get(paste0("a.", ncomp))
  testing <- testing  <- a[-inTrain.tissue, ]
  # import cv.fit
  fname <- paste0(folder, "fits/", name, "V_MLTA_fitCV_", ncomp, ".rds")
  fit.cv <- readRDS(file = fname)
  
  # use fuctions: plot(), print(), coef(), predict()
  
  
  #---All betas and analyses are at lambda.1se (lambda.min would result in less sparsity)
  betas <- coef(fit.cv, s = "lambda.1se") # a list
  outcomes <- names(betas)
  
  betas <- purrr::map(betas, .f = as.matrix)
  betas <- do.call(cbind, betas)
  colnames(betas) <- outcomes
  
  betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble()
  betas <- betas %>% dplyr::filter_at(vars(-compID, -other), any_vars(.!=0)) # disregard "other only" != 0 coefficients
  print(betas, n=ncomp) # also kable it !
  ave.spars <- nrow(betas)-1
  
  fname <- paste0(folder, "fits/", name, "V_MLTA_betas_", ncomp, ".csv")
  write_excel_csv(betas, path = fname)
  
  # make heatplot for easier assesment of coefficients
  # library("ComplexHeatmap")
  b.mat <- as.matrix(betas[2:nrow(betas), 2:ncol(betas)])
  rownames(b.mat) <- betas$compID[2:nrow(betas)]
  colnames(b.mat) <- colnames(b.mat) %>% str_replace_all("_", " ")
  b.mat <- t(b.mat)
  
  # maybe exclude weak coefficients from visualisation?.. 0.5sd away from zero?..
  plot(density(b.mat))
  cutof <- sd(b.mat)/2
  
  # library("circlize")
  col_fun = colorRamp2(breaks = c(min(b.mat)+2, -cutof, cutof, max(b.mat)-2), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
  cols <- col_fun(seq(min(b.mat), max(b.mat), by = 2))
  ratio <- ncol(b.mat) / nrow(b.mat)
  
  #---first see how to exclude small values +/- ~ 0.5 or 0.25 sd and less !
  ht <- Heatmap(b.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
                row_names_max_width = max_text_width(rownames(b.mat)),
                column_dend_height = unit(0.7, "cm"), row_dend_width = unit(0.8, "cm"), name = "",
                clustering_distance_rows="spearman", clustering_distance_columns="spearman",
                clustering_method_rows="average", clustering_method_columns="average",
                width = unit(10, "in"), height = unit(2*10/ratio, "in"),
                show_heatmap_legend=FALSE, na_col = "white", col=cols)
  # col=cols, cluster_rows=F, cluster_columns=F, clustering_distance_rows="spearman", clustering_method_rows="ward.D", "ward.D2", "single", "complete", "average"
  ht <- draw(ht)
  w <- as.numeric(ComplexHeatmap:::width(ht)) / 25.2 # mm to in.
  h <- as.numeric(ComplexHeatmap:::height(ht)) / 25.2 # mm to in.
  
  fname <- paste0(folder, "fits/", name, "V_MLTA_heat_", ncomp, ".pdf")
  # pdf(file = fname, width = w, height = h)
  ht
  dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
  dev.off()
  
  
  #---------Predictions and Assesment-------------#
  newx <- as.matrix(testing[, as.character(1:ncomp)])
  
  #--Number (and indices) of non-zero coefficients
  prediction.non0 <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "nonzero") # length of list elements are numbers of coefficients at lambda.1se!
  
  #--Predicted classes
  prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
  prediction <- factor(prediction, levels = levels(testing$use_celltype_tissue))
  
  #--Predicted probabilities
  prediction.prob <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "response") # array n.cells x n.factors x 1 
  prediction.prob <- as.matrix(prediction.prob[, 1:length(outcomes), 1]) %>% as.data.frame() %>% as_tibble()
  # assess if the prediction probabilities are high enough
  prediction.prob <- prediction.prob %>% 
    add_column(prediction = prediction) %>% 
    dplyr::mutate(max.prob = pmax(granulocyte_Marrow,promonocyte_Marrow,monocyte_Marrow,monocyte_Lung,macrophage_Marrow,macrophage_Spleen,macrophage_Lung, 
                                  macrophage_Kidney,macrophage_Mammary_Gland,macrophage_Limb_Muscle,dendritic_cell_Spleen,mast_cell_Lung,other)) %>% 
    dplyr::mutate(prediction_alt = if_else(max.prob > 0.55, as.character(prediction), "no_confidence"))
  # see how many cells have < 0.55 pred. prob.
  prediction.prob %>% dplyr::filter(max.prob < 0.55)
  # can try to use this alternative predicted classes with "no_confidence" class in confusion statistics
  # prediction_alt <- pull(prediction.prob, prediction_alt)
  
  # library("mixOmics") # alt. with more options: "measures" package
  confussion <- get.confusion_matrix(truth = testing$use_celltype_tissue, predicted = prediction)
  ave.ber <- get.BER(confussion)
  ave.ber
  
  pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
  confussion <- as.data.frame(confussion) %>% rownames_to_column()
  colnames(confussion) <- c("(Truth)", pred.names)
  print(confussion) # also kable it !
  
  fname <- paste0(folder, "fits/", name, "V_MLTA_confussion_", ncomp, ".csv")
  write_excel_csv(confussion, path = fname)
  
  # confussion_alt <- get.confusion_matrix(truth = testing$use_celltype, predicted = prediction_alt)
  
  # library("caret")
  confussion.c <- confusionMatrix(data=prediction, reference=testing$use_celltype_tissue, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
  confussion.c$overall[1:2] %>% round(digits = 3)
  confussionByClass <-confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3)
  # Balanced Accuracy is derived from Sensitivity & Specificity
  # F1 is derived form Precision & Recall
  
  # confussion_alt.c <- confusionMatrix(data=prediction_alt, reference=testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") # or "prec_recall", or "everything"
  
  # Tuck-in some of these stats in A_stats list, save and plot later
  # within them: nbetas, ave.sparsity, ave.ber, ave.kappa, ave.bacc, bacc
  #-nbetas:
  nbeta <- purrr::map_depth(prediction.non0, 2, .f=length) %>% unlist() %>% .[1:length(outcomes)-1]
  j <- ncomp/100
  A_stats$V$nbetas[j, ] <- nbeta
  #-ave.sparsity:
  A_stats$V$ave.sparsity[j] <- ave.spars
  #-ave.ber:
  A_stats$V$ave.ber[j] <- ave.ber
  #-ave.kappa:
  A_stats$V$ave.kappa[j] <- confussion.c$overall[2]
  #-ave.bacc:
  ave.bac <- confussion.c$overall[c(1,3:4)]
  A_stats$V$ave.bacc[j, ] <- ave.bac
  #-bacc:
  A_stats$V$bacc[j, ] <- confussionByClass[,6]
  # save this at the end of all I-VI analyses !
  
}

# .d
for (ncomp in c(100, 200, 300)) {
  # import a. if not loaded
  a <- get(paste0("a.", ncomp, ".d"))
  testing <- testing  <- a[-inTrain.tissue, ]
  # import cv.fit
  fname <- paste0(folder, "fits/", name, "V_MLTA_fitCV_", ncomp, "_D", ".rds")
  fit.cv <- readRDS(file = fname)
  
  # use fuctions: plot(), print(), coef(), predict()
  
  
  #---All betas and analyses are at lambda.1se (lambda.min would result in less sparsity)
  betas <- coef(fit.cv, s = "lambda.1se") # a list
  outcomes <- names(betas)
  
  betas <- purrr::map(betas, .f = as.matrix)
  betas <- do.call(cbind, betas)
  colnames(betas) <- outcomes
  
  betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble()
  betas <- betas %>% dplyr::filter_at(vars(-compID, -other), any_vars(.!=0)) # disregard "other only" != 0 coefficients
  print(betas, n=ncomp) # also kable it !
  ave.spars <- nrow(betas)-1
  
  fname <- paste0(folder, "fits/", name, "V_MLTA_betas_", ncomp, "_D", ".csv")
  write_excel_csv(betas, path = fname)
  
  # make heatplot for easier assesment of coefficients
  # library("ComplexHeatmap")
  b.mat <- as.matrix(betas[2:nrow(betas), 2:ncol(betas)])
  rownames(b.mat) <- betas$compID[2:nrow(betas)]
  colnames(b.mat) <- colnames(b.mat) %>% str_replace_all("_", " ")
  b.mat <- t(b.mat)
  
  # maybe exclude weak coefficients from visualisation?.. 0.5sd away from zero?..
  plot(density(b.mat))
  cutof <- sd(b.mat)/2
  
  # library("circlize")
  col_fun = colorRamp2(breaks = c(min(b.mat)+2, -cutof, cutof, max(b.mat)-2), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
  cols <- col_fun(seq(min(b.mat), max(b.mat), by = 2))
  ratio <- ncol(b.mat) / nrow(b.mat)
  
  #---first see how to exclude small values +/- ~ 0.5 or 0.25 sd and less !
  ht <- Heatmap(b.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
                row_names_max_width = max_text_width(rownames(b.mat)),
                column_dend_height = unit(0.7, "cm"), row_dend_width = unit(0.8, "cm"), name = "",
                clustering_distance_rows="spearman", clustering_distance_columns="spearman",
                clustering_method_rows="average", clustering_method_columns="average",
                width = unit(10, "in"), height = unit(2*10/ratio, "in"),
                show_heatmap_legend=FALSE, na_col = "white", col=cols)
  # col=cols, cluster_rows=F, cluster_columns=F, clustering_distance_rows="spearman", clustering_method_rows="ward.D", "ward.D2", "single", "complete", "average"
  ht <- draw(ht)
  w <- as.numeric(ComplexHeatmap:::width(ht)) / 25.2 # mm to in.
  h <- as.numeric(ComplexHeatmap:::height(ht)) / 25.2 # mm to in.
  
  fname <- paste0(folder, "fits/", name, "V_MLTA_heat_", ncomp, "_D", ".pdf")
  # pdf(file = fname, width = w, height = h)
  ht
  dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
  dev.off()
  
  
  #---------Predictions and Assesment-------------#
  newx <- as.matrix(testing[, as.character(1:ncomp)])
  
  #--Number (and indices) of non-zero coefficients
  prediction.non0 <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "nonzero") # length of list elements are numbers of coefficients at lambda.1se!
  
  #--Predicted classes
  prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
  prediction <- factor(prediction, levels = levels(testing$use_celltype_tissue))
  
  #--Predicted probabilities
  prediction.prob <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "response") # array n.cells x n.factors x 1 
  prediction.prob <- as.matrix(prediction.prob[, 1:length(outcomes), 1]) %>% as.data.frame() %>% as_tibble()
  # assess if the prediction probabilities are high enough
  prediction.prob <- prediction.prob %>% 
    add_column(prediction = prediction) %>% 
    dplyr::mutate(max.prob = pmax(granulocyte_Marrow,promonocyte_Marrow,monocyte_Marrow,monocyte_Lung,macrophage_Marrow,macrophage_Spleen,macrophage_Lung, 
                                  macrophage_Kidney,macrophage_Mammary_Gland,macrophage_Limb_Muscle,dendritic_cell_Spleen,mast_cell_Lung,other)) %>% 
    dplyr::mutate(prediction_alt = if_else(max.prob > 0.55, as.character(prediction), "no_confidence"))
  # see how many cells have < 0.55 pred. prob.
  prediction.prob %>% dplyr::filter(max.prob < 0.55)
  # can try to use this alternative predicted classes with "no_confidence" class in confusion statistics
  # prediction_alt <- pull(prediction.prob, prediction_alt)
  
  # library("mixOmics") # alt. with more options: "measures" package
  confussion <- get.confusion_matrix(truth = testing$use_celltype_tissue, predicted = prediction)
  ave.ber <- get.BER(confussion)
  ave.ber
  
  pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
  confussion <- as.data.frame(confussion) %>% rownames_to_column()
  colnames(confussion) <- c("(Truth)", pred.names)
  print(confussion) # also kable it !
  
  fname <- paste0(folder, "fits/", name, "V_MLTA_confussion_", ncomp, "_D", ".csv")
  write_excel_csv(confussion, path = fname)
  
  # confussion_alt <- get.confusion_matrix(truth = testing$use_celltype, predicted = prediction_alt)
  
  # library("caret")
  confussion.c <- confusionMatrix(data=prediction, reference=testing$use_celltype_tissue, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
  confussion.c$overall[1:2] %>% round(digits = 3)
  confussionByClass <-confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3)
  # Balanced Accuracy is derived from Sensitivity & Specificity
  # F1 is derived form Precision & Recall
  
  # confussion_alt.c <- confusionMatrix(data=prediction_alt, reference=testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") # or "prec_recall", or "everything"
  
  # Tuck-in some of these stats in A_stats list, save and plot later
  # within them: nbetas, ave.sparsity, ave.ber, ave.kappa, ave.bacc, bacc
  #-nbetas:
  nbeta <- purrr::map_depth(prediction.non0, 2, .f=length) %>% unlist() %>% .[1:length(outcomes)-1]
  j <- ncomp/100
  A_stats$V_d$nbetas[j, ] <- nbeta
  #-ave.sparsity:
  A_stats$V_d$ave.sparsity[j] <- ave.spars
  #-ave.ber:
  A_stats$V_d$ave.ber[j] <- ave.ber
  #-ave.kappa:
  A_stats$V_d$ave.kappa[j] <- confussion.c$overall[2]
  #-ave.bacc:
  ave.bac <- confussion.c$overall[c(1,3:4)]
  A_stats$V_d$ave.bacc[j, ] <- ave.bac
  #-bacc:
  A_stats$V_d$bacc[j, ] <- confussionByClass[,6]
  # save this at the end of all I-VI analyses !
  
}



### VI analysis: (MLTExp)
#-----------------------------------------------------------------------------#

# for VI
outcomes <- c("granulocyte_Marrow","promonocyte_Marrow","monocyte_Marrow","monocyte_Lung","macrophage_Marrow","macrophage_Spleen","macrophage_Lung", 
              "macrophage_Kidney","macrophage_Mammary_Gland","macrophage_Limb_Muscle","dendritic_cell_Spleen","mast_cell_Lung","other")
#-nbetas:
A_stats$VI$nbetas <- matrix(0, nrow = 1, ncol = length(outcomes)-1)
colnames(A_stats$VI$nbetas) <- outcomes[1:length(outcomes)-1]
#-ave.sparsity:
A_stats$VI$ave.sparsity <- vector(mode = "numeric", length = 1)
#-ave.ber:
A_stats$VI$ave.ber <- vector(mode = "numeric", length = 1)
#-ave.kappa:
A_stats$VI$ave.kappa <- vector(mode = "numeric", length = 1)
#-ave.bacc:
A_stats$VI$ave.bacc <- matrix(0, nrow = 1, ncol = 3)
colnames(A_stats$VI$ave.bacc) <- c("Accuracy", "loCI", "hiCI")
#-bacc:
A_stats$VI$bacc <- matrix(0, nrow = 1, ncol = length(outcomes))
colnames(A_stats$VI$bacc) <- outcomes[1:length(outcomes)]



# import exp.df if not loaded
exp.testing  <- exp.df[-inTrain.tissue, ]
# import cv.fit
fname <- paste0(folder, "fits/", name, "VI_MLTExp_fitCV", ".rds")
fit.cv <- readRDS(file = fname)

# use fuctions: plot(), print(), coef(), predict()


#---All betas and analyses are at lambda.1se (lambda.min would result in less sparsity)
betas <- coef(fit.cv, s = "lambda.1se") # a list
outcomes <- names(betas)

betas <- purrr::map(betas, .f = as.matrix)
betas <- do.call(cbind, betas)
colnames(betas) <- outcomes

betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble()
betas <- betas %>% dplyr::filter_at(vars(-compID, -other), any_vars(.!=0)) # disregard "other only" != 0 coefficients
print(betas, n=250) # also kable it !
ave.spars <- nrow(betas)-1

fname <- paste0(folder, "fits/", name, "VI_MLTExp_betas", ".csv")
write_excel_csv(betas, path = fname)

# make heatplot for easier assesment of coefficients
# library("ComplexHeatmap")
b.mat <- as.matrix(betas[2:nrow(betas), 2:ncol(betas)])
rownames(b.mat) <- betas$compID[2:nrow(betas)]
colnames(b.mat) <- colnames(b.mat) %>% str_replace_all("_", " ")
b.mat <- t(b.mat)

# maybe exclude weak coefficients from visualisation?.. 0.5sd away from zero?..
plot(density(b.mat))
cutof <- sd(b.mat)/8

# library("circlize")
col_fun = colorRamp2(breaks = c(min(b.mat), -cutof, cutof, max(b.mat)), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
cols <- col_fun(seq(min(b.mat), max(b.mat), by = 0.002))
ratio <- ncol(b.mat) / nrow(b.mat)

#---first see how to exclude small values +/- ~ 0.5 or 0.25 sd and less !
ht <- Heatmap(b.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(b.mat)),
              column_dend_height = unit(0.7, "cm"), row_dend_width = unit(0.8, "cm"), name = "",
              clustering_distance_rows="spearman", clustering_distance_columns="spearman",
              clustering_method_rows="average", clustering_method_columns="average",
              width = unit(10, "in"), height = unit(2*10/ratio, "in"),
              show_heatmap_legend=FALSE, na_col = "white", col=cols)
# col=cols, cluster_rows=F, cluster_columns=F, clustering_distance_rows="spearman", clustering_method_rows="ward.D", "ward.D2", "single", "complete", "average"
ht <- draw(ht)
w <- as.numeric(ComplexHeatmap:::width(ht)) / 25.2 # mm to in.
h <- as.numeric(ComplexHeatmap:::height(ht)) / 25.2 # mm to in.

fname <- paste0(folder, "fits/", name, "VI_MLTExp_heat", ".pdf")
# pdf(file = fname, width = w, height = h)
ht
dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
dev.off()


#---------Predictions and Assesment-------------#
newx <- as.matrix(exp.testing[, hvg])

#--Number (and indices) of non-zero coefficients
prediction.non0 <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "nonzero") # length of list elements are numbers of coefficients at lambda.1se!

#--Predicted classes
prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
prediction <- factor(prediction, levels = levels(exp.testing$use_celltype_tissue))

#--Predicted probabilities
prediction.prob <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "response") # array n.cells x n.factors x 1 
prediction.prob <- as.matrix(prediction.prob[, 1:length(outcomes), 1]) %>% as.data.frame() %>% as_tibble()
# assess if the prediction probabilities are high enough
prediction.prob <- prediction.prob %>% 
  add_column(prediction = prediction) %>% 
  dplyr::mutate(max.prob = pmax(granulocyte_Marrow,promonocyte_Marrow,monocyte_Marrow,monocyte_Lung,macrophage_Marrow,macrophage_Spleen,macrophage_Lung, 
                                macrophage_Kidney,macrophage_Mammary_Gland,macrophage_Limb_Muscle,dendritic_cell_Spleen,mast_cell_Lung,other)) %>% 
  dplyr::mutate(prediction_alt = if_else(max.prob > 0.55, as.character(prediction), "no_confidence"))
# see how many cells have < 0.55 pred. prob.
print(prediction.prob %>% dplyr::filter(max.prob < 0.55))
# can try to use this alternative predicted classes with "no_confidence" class in confusion statistics
# prediction_alt <- pull(prediction.prob, prediction_alt)

# library("mixOmics") # alt. with more options: "measures" package
confussion <- get.confusion_matrix(truth = exp.testing$use_celltype_tissue, predicted = prediction)
ave.ber <- get.BER(confussion)
ave.ber

pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
confussion <- as.data.frame(confussion) %>% rownames_to_column()
colnames(confussion) <- c("(Truth)", pred.names)
print(confussion) # also kable it !

fname <- paste0(folder, "fits/", name, "VI_MLTExp_confussion", ".csv")
write_excel_csv(confussion, path = fname)

# confussion_alt <- get.confusion_matrix(truth = testing$use_celltype, predicted = prediction_alt)

# library("caret")
confussion.c <- confusionMatrix(data=prediction, reference=exp.testing$use_celltype_tissue, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
confussion.c$overall[1:2] %>% round(digits = 3)
confussionByClass <-confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3)
# Balanced Accuracy is derived from Sensitivity & Specificity
# F1 is derived form Precision & Recall

# confussion_alt.c <- confusionMatrix(data=prediction_alt, reference=testing$use_celltype, dnn = c("Prediction", "Reference"), mode = "sens_spec") # or "prec_recall", or "everything"

# Tuck-in some of these stats in A_stats list, save and plot later
# within them: nbetas, ave.sparsity, ave.ber, ave.kappa, ave.bacc, bacc
#-nbetas:
nbeta <- purrr::map_depth(prediction.non0, 2, .f=length) %>% unlist() %>% .[1:length(outcomes)-1]
A_stats$VI$nbetas[1, ] <- nbeta
#-ave.sparsity:
A_stats$VI$ave.sparsity[1] <- ave.spars
#-ave.ber:
A_stats$VI$ave.ber[1] <- ave.ber
#-ave.kappa:
A_stats$VI$ave.kappa[1] <- confussion.c$overall[2]
#-ave.bacc:
ave.bac <- confussion.c$overall[c(1,3:4)]
A_stats$VI$ave.bacc[1, ] <- ave.bac
#-bacc:
A_stats$VI$bacc[1, ] <- confussionByClass[,6]
# save this at the end of all I-VI analyses !


fname <- paste0(folder, "fits/", name, "A_stats", ".rds")
saveRDS(A_stats, file = fname)

########################
########################
######-UNCOMENT-########


# Now plot these summary stats 

# Global Myelo and Macs, Consensus comps and Exp with Lasso vs. Elastic net 
#-------------------------------------------------------------------#
nComps <- c(100, 200, 300, 500) # 500 is for exp..
# gen-consensus:
## nbeta
gen.ave.nbeta.c.lasso <- c(A_stats$I$ave.sparsity, A_stats$III$ave.sparsity) # I+III
gen.ave.nbeta.c.enet <- c(A_stats$II$ave.sparsity, A_stats$IV$ave.sparsity) # II+IV
## ber
gen.ave.1.ber.c.lasso <- 1 - c(A_stats$I$ave.ber, A_stats$III$ave.ber)
gen.ave.1.ber.c.enet <- 1 - c(A_stats$II$ave.ber, A_stats$IV$ave.ber)
## acc
gen.ave.acc.c.lasso <- c(A_stats$I$ave.bacc[,1], A_stats$III$ave.bacc[,1])
gen.ave.acc.c.enet <- c(A_stats$II$ave.bacc[,1], A_stats$IV$ave.bacc[,1])
## kappa
gen.ave.kappa.c.lasso <- c(A_stats$I$ave.kappa, A_stats$III$ave.kappa)
gen.ave.kappa.c.enet <- c(A_stats$II$ave.kappa, A_stats$IV$ave.kappa)


fname <- paste0(folder, "figures/", name, "A_stats_GenConsLassoVsEnet", ".png")
png(filename = fname, width = 600, height = 500, units = "px") #  pointsize = 13
## add extra space to right margin of plot within frame
par(mar=c(6, 5, 5, 6) + 0.1)
## Plot first set of data and draw its axis
plot(nComps, gen.ave.nbeta.c.lasso, pch=16, type="b", lty=1, col="black", axes=FALSE, ylim=c(0,355), xlab="", ylab="")
points(nComps, gen.ave.nbeta.c.enet, pch=21, type="b", lty=2, col="black", ylim=c(0,355), xlab="", ylab="")
mtext(side=3, line=3, cex=1.3, "Sparsity and Predictive performance (Lasso vs. Elastic Net)")
mtext(side=3, line=1.25, cex=1.1, "Myeloid cells")
mtext(side=2, line=4, cex=1.2, "Total number of active coefficients (components or genes)", col="black")
axis(2, ylim=c(0,125), col="black", col.axis="black", las=1)  ## las=1 makes horizontal labels
box()
## Allow a second plot on the same graph
par(new=TRUE)
## Plot the second plot and put axis scale on right
plot(nComps, gen.ave.kappa.c.lasso, pch=16, type="b", lty=1, col="red", axes=FALSE, ylim=c(0.55,1), xlab="", ylab="")
points(nComps, gen.ave.kappa.c.enet, pch=21, type="b", lty=2, col="red", ylim=c(0.55,1), xlab="", ylab="")

points(nComps, gen.ave.1.ber.c.lasso, pch=16, type="b", lty=1, col="darkorange2", ylim=c(0.55,1), xlab="", ylab="")
points(nComps, gen.ave.1.ber.c.enet, pch=21, type="b", lty=2, col="darkorange2", ylim=c(0.55,1), xlab="", ylab="")

points(nComps, gen.ave.acc.c.lasso, pch=16, type="b", lty=1, col="darkred", ylim=c(0.55,1), xlab="", ylab="")
points(nComps, gen.ave.acc.c.enet, pch=21, type="b", lty=2, col="darkred", ylim=c(0.55,1), xlab="", ylab="")
## a little farther out (line=4) to make room for labels
mtext("Predictive performance statistics", side=4, col="orangered3", line=4, cex=1.25)
axis(4, ylim=c(0.55,1), col="orangered3",col.axis="orangered3", las=1)
## Draw the time axis
axis(1, at = nComps, labels = c("100 components", "200 components", "300 components", "7500 hv genes"))
## Add Legend
legend("bottom",legend=c("1-Balanced error rate","Kappa statistics","Overal accuracy","Number of coefficients","Lasso","Elastic net (alpha 0.95)"), 
       text.col=c("darkorange2","red","darkred","black","grey45","grey45"), 
       col=c("darkorange2","red","darkred","black","grey45","grey45"), 
       pch=c(16,16,16,16,16,21), lty=c(0,0,0,0,1,2), ncol = 2, inset=0.01)
dev.off()


# mac-consensus:
## nbeta
mac.ave.nbeta.c.lasso <- c(A_stats$Ia$ave.sparsity, A_stats$IIIa$ave.sparsity) # Ia+IIIa
mac.ave.nbeta.c.enet <- c(A_stats$IIa$ave.sparsity, A_stats$IVa$ave.sparsity) # IIa+IVa
## ber
mac.ave.1.ber.c.lasso <- 1 - c(A_stats$Ia$ave.ber, A_stats$IIIa$ave.ber)
mac.ave.1.ber.c.enet <- 1 - c(A_stats$IIa$ave.ber, A_stats$IVa$ave.ber)
## acc
mac.ave.acc.c.lasso <- c(A_stats$Ia$ave.bacc[,1], A_stats$IIIa$ave.bacc[,1])
mac.ave.acc.c.enet <- c(A_stats$IIa$ave.bacc[,1], A_stats$IVa$ave.bacc[,1])
## kappa
mac.ave.kappa.c.lasso <- c(A_stats$Ia$ave.kappa, A_stats$IIIa$ave.kappa)
mac.ave.kappa.c.enet <- c(A_stats$IIa$ave.kappa, A_stats$IVa$ave.kappa)


fname <- paste0(folder, "figures/", name, "A_stats_MacConsLassoVsEnet", ".png")
png(filename = fname, width = 600, height = 500, units = "px") #  pointsize = 13
## add extra space to right margin of plot within frame
par(mar=c(6, 5, 5, 6) + 0.1)
## Plot first set of data and draw its axis
plot(nComps, mac.ave.nbeta.c.lasso, pch=16, type="b", lty=1, col="black", axes=FALSE, ylim=c(0,355), xlab="", ylab="")
points(nComps, mac.ave.nbeta.c.enet, pch=21, type="b", lty=2, col="black", ylim=c(0,355), xlab="", ylab="")
mtext(side=3, line=3, cex=1.3, "Sparsity and Predictive performance (Lasso vs. Elastic Net)")
mtext(side=3, line=1.25, cex=1.1, "Macrophages")
mtext(side=2, line=4, cex=1.2, "Total number of active coefficients (components or genes)", col="black")
axis(2, ylim=c(0,125), col="black", col.axis="black", las=1)  ## las=1 makes horizontal labels
box()
## Allow a second plot on the same graph
par(new=TRUE)
## Plot the second plot and put axis scale on right
plot(nComps, mac.ave.kappa.c.lasso, pch=16, type="b", lty=1, col="red", axes=FALSE, ylim=c(0.85,1), xlab="", ylab="")
points(nComps, mac.ave.kappa.c.enet, pch=21, type="b", lty=2, col="red", ylim=c(0.85,1), xlab="", ylab="")

points(nComps, mac.ave.1.ber.c.lasso, pch=16, type="b", lty=1, col="darkorange2", ylim=c(0.85,1), xlab="", ylab="")
points(nComps, mac.ave.1.ber.c.enet, pch=21, type="b", lty=2, col="darkorange2", ylim=c(0.85,1), xlab="", ylab="")

points(nComps, mac.ave.acc.c.lasso, pch=16, type="b", lty=1, col="darkred", ylim=c(0.85,1), xlab="", ylab="")
points(nComps, mac.ave.acc.c.enet, pch=21, type="b", lty=2, col="darkred", ylim=c(0.85,1), xlab="", ylab="")
## a little farther out (line=4) to make room for labels
mtext("Predictive performance statistics", side=4, col="orangered3", line=4, cex=1.25)
axis(4, ylim=c(0.85,1), col="orangered3",col.axis="orangered3", las=1)
## Draw the time axis
axis(1, at = nComps, labels = c("100 components", "200 components", "300 components", "7500 hv genes"))
## Add Legend
legend("bottom",legend=c("1-Balanced error rate","Kappa statistics","Overal accuracy","Number of coefficients","Lasso","Elastic net (alpha 0.95)"), 
       text.col=c("darkorange2","red","darkred","black","grey45","grey45"), 
       col=c("darkorange2","red","darkred","black","grey45","grey45"), 
       pch=c(16,16,16,16,16,21), lty=c(0,0,0,0,1,2), ncol = 2, inset=0.01)
dev.off()




# Global Myelo and Macs, Diagonal comps and Exp with Lasso vs. Elastic net
#-------------------------------------------------------------------#
nComps <- c(100, 200, 300, 500) # 500 is for exp..
# gen-diagonal:
## nbeta
gen.ave.nbeta.dg.lasso <- c(A_stats$I_d$ave.sparsity, A_stats$III$ave.sparsity) # I_d+III
gen.ave.nbeta.dg.enet <- c(A_stats$II_d$ave.sparsity, A_stats$IV$ave.sparsity) # II_d+IV
## ber
gen.ave.1.ber.dg.lasso <- 1 - c(A_stats$I_d$ave.ber, A_stats$III$ave.ber)
gen.ave.1.ber.dg.enet <- 1 - c(A_stats$II_d$ave.ber, A_stats$IV$ave.ber)
## acc
gen.ave.acc.dg.lasso <- c(A_stats$I_d$ave.bacc[,1], A_stats$III$ave.bacc[,1])
gen.ave.acc.dg.enet <- c(A_stats$II_d$ave.bacc[,1], A_stats$IV$ave.bacc[,1])
## kappa
gen.ave.kappa.dg.lasso <- c(A_stats$I_d$ave.kappa, A_stats$III$ave.kappa)
gen.ave.kappa.dg.enet <- c(A_stats$II_d$ave.kappa, A_stats$IV$ave.kappa)


fname <- paste0(folder, "figures/", name, "A_stats_GenDiagLassoVsEnet", ".png")
png(filename = fname, width = 600, height = 500, units = "px") #  pointsize = 13
## add extra space to right margin of plot within frame
par(mar=c(6, 5, 5, 6) + 0.1)
## Plot first set of data and draw its axis
plot(nComps, gen.ave.nbeta.dg.lasso, pch=16, type="b", lty=1, col="black", axes=FALSE, ylim=c(0,355), xlab="", ylab="")
points(nComps, gen.ave.nbeta.dg.enet, pch=21, type="b", lty=2, col="black", ylim=c(0,355), xlab="", ylab="")
mtext(side=3, line=3, cex=1.3, "Sparsity and Predictive performance (Lasso vs. Elastic Net)")
mtext(side=3, line=1.25, cex=1.1, "Myeloid cells")
mtext(side=2, line=4, cex=1.2, "Total number of active coefficients (components or genes)", col="black")
axis(2, ylim=c(0,125), col="black", col.axis="black", las=1)  ## las=1 makes horizontal labels
box()
## Allow a second plot on the same graph
par(new=TRUE)
## Plot the second plot and put axis scale on right
plot(nComps, gen.ave.kappa.dg.lasso, pch=16, type="b", lty=1, col="red", axes=FALSE, ylim=c(0.55,1), xlab="", ylab="")
points(nComps, gen.ave.kappa.dg.enet, pch=21, type="b", lty=2, col="red", ylim=c(0.55,1), xlab="", ylab="")

points(nComps, gen.ave.1.ber.dg.lasso, pch=16, type="b", lty=1, col="darkorange2", ylim=c(0.55,1), xlab="", ylab="")
points(nComps, gen.ave.1.ber.dg.enet, pch=21, type="b", lty=2, col="darkorange2", ylim=c(0.55,1), xlab="", ylab="")

points(nComps, gen.ave.acc.dg.lasso, pch=16, type="b", lty=1, col="darkred", ylim=c(0.55,1), xlab="", ylab="")
points(nComps, gen.ave.acc.dg.enet, pch=21, type="b", lty=2, col="darkred", ylim=c(0.55,1), xlab="", ylab="")
## a little farther out (line=4) to make room for labels
mtext("Predictive performance statistics", side=4, col="orangered3", line=4, cex=1.25)
axis(4, ylim=c(0.55,1), col="orangered3",col.axis="orangered3", las=1)
## Draw the time axis
axis(1, at = nComps, labels = c("100 components", "200 components", "300 components", "7500 hv genes"))
## Add Legend
legend("bottom",legend=c("1-Balanced error rate","Kappa statistics","Overal accuracy","Number of coefficients","Lasso","Elastic net (alpha 0.95)"), 
       text.col=c("darkorange2","red","darkred","black","grey45","grey45"), 
       col=c("darkorange2","red","darkred","black","grey45","grey45"), 
       pch=c(16,16,16,16,16,21), lty=c(0,0,0,0,1,2), ncol = 2, inset=0.01)
dev.off()


# mac-diagonal:
## nbeta
mac.ave.nbeta.dg.lasso <- c(A_stats$Ia_d$ave.sparsity, A_stats$IIIa$ave.sparsity) # Ia_d+IIIa
mac.ave.nbeta.dg.enet <- c(A_stats$IIa_d$ave.sparsity, A_stats$IVa$ave.sparsity) # IIa_d+IVa
## ber
mac.ave.1.ber.dg.lasso <- 1 - c(A_stats$Ia_d$ave.ber, A_stats$IIIa$ave.ber)
mac.ave.1.ber.dg.enet <- 1 - c(A_stats$IIa_d$ave.ber, A_stats$IVa$ave.ber)
## acc
mac.ave.acc.dg.lasso <- c(A_stats$Ia_d$ave.bacc[,1], A_stats$IIIa$ave.bacc[,1])
mac.ave.acc.dg.enet <- c(A_stats$IIa_d$ave.bacc[,1], A_stats$IVa$ave.bacc[,1])
## kappa
mac.ave.kappa.dg.lasso <- c(A_stats$Ia_d$ave.kappa, A_stats$IIIa$ave.kappa)
mac.ave.kappa.dg.enet <- c(A_stats$IIa_d$ave.kappa, A_stats$IVa$ave.kappa)


fname <- paste0(folder, "figures/", name, "A_stats_MacDiagLassoVsEnet", ".png")
png(filename = fname, width = 600, height = 500, units = "px") #  pointsize = 13
## add extra space to right margin of plot within frame
par(mar=c(6, 5, 5, 6) + 0.1)
## Plot first set of data and draw its axis
plot(nComps, mac.ave.nbeta.dg.lasso, pch=16, type="b", lty=1, col="black", axes=FALSE, ylim=c(0,355), xlab="", ylab="")
points(nComps, mac.ave.nbeta.dg.enet, pch=21, type="b", lty=2, col="black", ylim=c(0,355), xlab="", ylab="")
mtext(side=3, line=3, cex=1.3, "Sparsity and Predictive performance (Lasso vs. Elastic Net)")
mtext(side=3, line=1.25, cex=1.1, "Macrophages")
mtext(side=2, line=4, cex=1.2, "Total number of active coefficients (components or genes)", col="black")
axis(2, ylim=c(0,125), col="black", col.axis="black", las=1)  ## las=1 makes horizontal labels
box()
## Allow a second plot on the same graph
par(new=TRUE)
## Plot the second plot and put axis scale on right
plot(nComps, mac.ave.kappa.dg.lasso, pch=16, type="b", lty=1, col="red", axes=FALSE, ylim=c(0.85,1), xlab="", ylab="")
points(nComps, mac.ave.kappa.dg.enet, pch=21, type="b", lty=2, col="red", ylim=c(0.85,1), xlab="", ylab="")

points(nComps, mac.ave.1.ber.dg.lasso, pch=16, type="b", lty=1, col="darkorange2", ylim=c(0.85,1), xlab="", ylab="")
points(nComps, mac.ave.1.ber.dg.enet, pch=21, type="b", lty=2, col="darkorange2", ylim=c(0.85,1), xlab="", ylab="")

points(nComps, mac.ave.acc.dg.lasso, pch=16, type="b", lty=1, col="darkred", ylim=c(0.85,1), xlab="", ylab="")
points(nComps, mac.ave.acc.dg.enet, pch=21, type="b", lty=2, col="darkred", ylim=c(0.85,1), xlab="", ylab="")
## a little farther out (line=4) to make room for labels
mtext("Predictive performance statistics", side=4, col="orangered3", line=4, cex=1.25)
axis(4, ylim=c(0.85,1), col="orangered3",col.axis="orangered3", las=1)
## Draw the time axis
axis(1, at = nComps, labels = c("100 components", "200 components", "300 components", "7500 hv genes"))
## Add Legend
legend("bottom",legend=c("1-Balanced error rate","Kappa statistics","Overal accuracy","Number of coefficients","Lasso","Elastic net (alpha 0.95)"), 
       text.col=c("darkorange2","red","darkred","black","grey45","grey45"), 
       col=c("darkorange2","red","darkred","black","grey45","grey45"), 
       pch=c(16,16,16,16,16,21), lty=c(0,0,0,0,1,2), ncol = 2, inset=0.01)
dev.off()





# Cell type breakdown Myelo and Macs, Consensus comps and Exp with Lasso
#-------------------------------------------------------------------#
# gen-consensus:
## nbeta
gen.nbeta.c.lasso <- rbind(A_stats$I$nbetas, A_stats$III$nbetas) %>% t(); colnames(gen.nbeta.c.lasso) <- c("100 c", "200 c", "300 c", "7500 hvg") # I+III (by celltype)
gen.nbeta.c.enet <- rbind(A_stats$II$nbetas, A_stats$IV$nbetas) %>% t(); colnames(gen.nbeta.c.enet) <- c("100 c", "200 c", "300 c", "7500 hvg") # II+IV (by celltype)
## acc
gen.acc.c.lasso <- rbind(A_stats$I$bacc[,-7], A_stats$III$bacc[,-7]) %>% t() %>% round(digits = 2); colnames(gen.acc.c.lasso) <- c("100 c", "200 c", "300 c", "7500 hvg")
gen.acc.c.enet <- rbind(A_stats$II$bacc[,-7], A_stats$IV$bacc[,-7]) %>% t() %>% round(digits = 2); colnames(gen.acc.c.enet) <- c("100 c", "200 c", "300 c", "7500 hvg")

library("RColorBrewer")
display.brewer.all(colorblindFriendly = TRUE)
cols <- c(brewer.pal(6, "Set2")[1], brewer.pal(6, "Purples")[c(2,4)], brewer.pal(6, "Set2")[4], brewer.pal(6, "Set2")[6], brewer.pal(7, "Reds")[6])

fname <- paste0(folder, "figures/", name, "A_stats_CellGenConsLasso", ".png")
png(filename = fname, width = 600, height = 500, units = "px") #  pointsize = 13
par(mar=c(4,6,5,3) + 0.1)
bp <- barplot(gen.nbeta.c.lasso, beside=TRUE, col = cols, ylim=c(0,80))
text(bp, gen.nbeta.c.lasso+2, labels = gen.acc.c.lasso, cex=0.7) 
mtext(side=3, line=3, cex=1.3, "Sparsity and Predictive performance (Lasso)")
mtext(side=3, line=1.25, cex=1.1, "Myeloid cells")
mtext(side=2, line=4, cex=1.1, "Number of active coefficients (components or genes)", col="black")
legend("topright", legend=c("granulocyte","promonocyte","monocyte","macrophage","dendritic cell","mast cell"), fill=cols) # bty="n",
dev.off()


# mac-consensus:
## nbeta
mac.nbeta.c.lasso <- rbind(A_stats$Ia$nbetas, A_stats$IIIa$nbetas) %>% t(); colnames(mac.nbeta.c.lasso) <- c("100 c", "200 c", "300 c", "7500 hvg") # I+III
mac.nbeta.c.enet <- rbind(A_stats$IIa$nbetas, A_stats$IVa$nbetas) %>% t(); colnames(mac.nbeta.c.enet) <- c("100 c", "200 c", "300 c", "7500 hvg") # II+IV
## acc
mac.acc.c.lasso <- rbind(A_stats$Ia$bacc, A_stats$IIIa$bacc) %>% t() %>% round(digits = 2); colnames(mac.acc.c.lasso) <- c("100 c", "200 c", "300 c", "7500 hvg")
mac.acc.c.enet <- rbind(A_stats$IIa$bacc, A_stats$IVa$bacc) %>% t() %>% round(digits = 2); colnames(mac.acc.c.enet) <- c("100 c", "200 c", "300 c", "7500 hvg")

# change colours..
cols <- c(brewer.pal(6, "Set2")[5], brewer.pal(7, "Reds")[6], brewer.pal(6, "Purples")[4], brewer.pal(6, "Set2")[2], brewer.pal(6, "Set2")[4], brewer.pal(7, "Set2")[7])

fname <- paste0(folder, "figures/", name, "A_stats_CellMacConsLasso", ".png")
png(filename = fname, width = 600, height = 500, units = "px") #  pointsize = 13
par(mar=c(4,6,5,3) + 0.1)
bp <- barplot(mac.nbeta.c.lasso, beside=TRUE, col = cols, ylim=c(0,100))
text(bp, mac.nbeta.c.lasso+2, labels = mac.acc.c.lasso, cex=0.7) 
mtext(side=3, line=3, cex=1.3, "Sparsity and Predictive performance (Lasso)")
mtext(side=3, line=1.25, cex=1.1, "Macrophages")
mtext(side=2, line=4, cex=1.1, "Number of active coefficients (components or genes)", col="black")
legend("topleft", legend=str_replace_all(rownames(mac.nbeta.dg.lasso), "_", " "), fill=cols, inset=0.01) # bty="n",
dev.off()

# no need to also plot elastic net ones, as 1 to 0.95 change did not result in any different result 
# time permitting try alpha = 0.9 or 0.8



# Cell type breakdown Myelo and Macs, Diagonal comps and Exp with Lasso vs. Elastic net
#-------------------------------------------------------------------#
# gen-diagonal:
## nbeta
gen.nbeta.dg.lasso <- rbind(A_stats$I_d$nbetas, A_stats$III$nbetas) %>% t(); colnames(gen.nbeta.dg.lasso) <- c("100 c", "200 c", "300 c", "7500 hvg") # I_d+III (by celltype)
gen.nbeta.dg.enet <- rbind(A_stats$II_d$nbetas, A_stats$IV$nbetas) %>% t(); colnames(gen.nbeta.dg.enet) <- c("100 c", "200 c", "300 c", "7500 hvg") # II_d+IV (by celltype)
## acc
gen.acc.dg.lasso <- rbind(A_stats$I_d$bacc[,-7], A_stats$III$bacc[,-7]) %>% t() %>% round(digits = 2); colnames(gen.acc.dg.lasso) <- c("100 c", "200 c", "300 c", "7500 hvg")
gen.acc.dg.enet <- rbind(A_stats$II_d$bacc[,-7], A_stats$IV$bacc[,-7]) %>% t() %>% round(digits = 2); colnames(gen.acc.dg.enet) <- c("100 c", "200 c", "300 c", "7500 hvg")


cols <- c(brewer.pal(6, "Set2")[1], brewer.pal(6, "Purples")[c(2,4)], brewer.pal(6, "Set2")[4], brewer.pal(6, "Set2")[6], brewer.pal(7, "Reds")[6])

fname <- paste0(folder, "figures/", name, "A_stats_CellGenDiagLasso", ".png")
png(filename = fname, width = 600, height = 500, units = "px") #  pointsize = 13
par(mar=c(4,6,5,3) + 0.1)
bp <- barplot(gen.nbeta.dg.lasso, beside=TRUE, col = cols, ylim=c(0,81))
text(bp, gen.nbeta.dg.lasso+2, labels = gen.acc.dg.lasso, cex=0.7) 
mtext(side=3, line=3, cex=1.3, "Sparsity and Predictive performance (Lasso)")
mtext(side=3, line=1.25, cex=1.1, "Myeloid cells")
mtext(side=2, line=4, cex=1.1, "Number of active coefficients (components or genes)", col="black")
legend("topright", legend=c("granulocyte","promonocyte","monocyte","macrophage","dendritic cell","mast cell"), fill=cols) # bty="n",
dev.off()


# mac-diagonal:
## nbeta
mac.nbeta.dg.lasso <- rbind(A_stats$Ia_d$nbetas, A_stats$IIIa$nbetas) %>% t(); colnames(mac.nbeta.dg.lasso) <- c("100 c", "200 c", "300 c", "7500 hvg") # Ia_d+III (by celltype)
mac.nbeta.dg.enet <- rbind(A_stats$IIa_d$nbetas, A_stats$IVa$nbetas) %>% t(); colnames(mac.nbeta.dg.enet) <- c("100 c", "200 c", "300 c", "7500 hvg") # IIa_d+IV (by celltype)
## acc
mac.acc.dg.lasso <- rbind(A_stats$Ia_d$bacc, A_stats$IIIa$bacc) %>% t() %>% round(digits = 2); colnames(mac.acc.dg.lasso) <- c("100 c", "200 c", "300 c", "7500 hvg")
mac.acc.dg.enet <- rbind(A_stats$IIa_d$bacc, A_stats$IVa$bacc) %>% t() %>% round(digits = 2); colnames(mac.acc.dg.enet) <- c("100 c", "200 c", "300 c", "7500 hvg")

# change colours..
cols <- c(brewer.pal(6, "Set2")[5], brewer.pal(7, "Reds")[6], brewer.pal(6, "Purples")[4], brewer.pal(6, "Set2")[2], brewer.pal(6, "Set2")[4], brewer.pal(7, "Set2")[7])

fname <- paste0(folder, "figures/", name, "A_stats_CellMacDiagLasso", ".png")
png(filename = fname, width = 600, height = 500, units = "px") #  pointsize = 13
par(mar=c(4,6,5,3) + 0.1)
bp <- barplot(mac.nbeta.dg.lasso, beside=TRUE, col = cols, ylim=c(0,100))
text(bp, mac.nbeta.dg.lasso+2, labels = mac.acc.dg.lasso, cex=0.7) 
mtext(side=3, line=3, cex=1.3, "Sparsity and Predictive performance (Lasso)")
mtext(side=3, line=1.25, cex=1.1, "Macrophages")
mtext(side=2, line=4, cex=1.1, "Number of active coefficients (components or genes)", col="black")
legend("topleft", legend=str_replace_all(rownames(mac.nbeta.dg.lasso), "_", " "), fill=cols, inset=0.01) # bty="n",
dev.off()

# no need to also plot elastic net ones, as 1 to 0.95 change did not result in any different result 
# time permitting try alpha = 0.9 or 0.8















# Tissue-Cell type breakdown Myelo, Consensus comps and Exp with Lasso only
#-------------------------------------------------------------------#
nComps <- c(100, 200, 300, 500) # 500 is for exp..
# gen-ts-consensus:
## nbeta
gen.ave.nbeta.c.lasso.ts <- c(A_stats$V$ave.sparsity, A_stats$VI$ave.sparsity) # V+VI
## ber
gen.ave.1.ber.c.lasso.ts <- 1 - c(A_stats$V$ave.ber, A_stats$VI$ave.ber)
## acc
gen.ave.acc.c.lasso.ts <- c(A_stats$V$ave.bacc[,1], A_stats$VI$ave.bacc[,1])
## kappa
gen.ave.kappa.c.lasso.ts <- c(A_stats$V$ave.kappa, A_stats$VI$ave.kappa)

fname <- paste0(folder, "figures/", name, "A_stats_GenTsConsLassoVsEnet", ".png")
png(filename = fname, width = 600, height = 500, units = "px") #  pointsize = 13
## add extra space to right margin of plot within frame
par(mar=c(6, 5, 5, 6) + 0.1)
## Plot first set of data and draw its axis
plot(nComps, gen.ave.nbeta.c.lasso.ts, pch=16, type="b", lty=1, col="black", axes=FALSE, ylim=c(0,355), xlab="", ylab="")
mtext(side=3, line=3, cex=1.3, "Sparsity and Predictive performance (Lasso)")
mtext(side=3, line=1.25, cex=1.1, "Tissue-specific Myeloid cells")
mtext(side=2, line=4, cex=1.2, "Total number of active coefficients (components or genes)", col="black")
axis(2, ylim=c(0,125), col="black", col.axis="black", las=1)  ## las=1 makes horizontal labels
box()
## Allow a second plot on the same graph
par(new=TRUE)
## Plot the second plot and put axis scale on right
plot(nComps, gen.ave.kappa.c.lasso.ts, pch=16, type="b", lty=1, col="red", axes=FALSE, ylim=c(0.55,1), xlab="", ylab="")

points(nComps, gen.ave.1.ber.c.lasso.ts, pch=16, type="b", lty=1, col="darkorange2", ylim=c(0.55,1), xlab="", ylab="")

points(nComps, gen.ave.acc.c.lasso.ts, pch=16, type="b", lty=1, col="darkred", ylim=c(0.55,1), xlab="", ylab="")
## a little farther out (line=4) to make room for labels
mtext("Predictive performance statistics", side=4, col="orangered3", line=4, cex=1.25)
axis(4, ylim=c(0.55,1), col="orangered3",col.axis="orangered3", las=1)
## Draw the time axis
axis(1, at = nComps, labels = c("100 components", "200 components", "300 components", "7500 hv genes"))
## Add Legend
legend("bottom",legend=c("1-Balanced error rate","Kappa statistics","Overal accuracy","Number of coefficients"), 
       text.col=c("darkorange2","red","darkred","black"), 
       col=c("darkorange2","red","darkred","black"), 
       pch=c(16,16,16,16), lty=c(1,1,1,1), ncol = 1, inset=0.01)
dev.off()








# Tissue-Cell type breakdown Myelo, Consensus comps and Exp with Lasso only
#-------------------------------------------------------------------#
# gen-consensus:
## nbeta
gen.nbeta.c.lasso.ts <- rbind(A_stats$V$nbetas, A_stats$VI$nbetas) %>% t(); colnames(gen.nbeta.c.lasso.ts) <- c("100 c", "200 c", "300 c", "7500 hvg") # V+VI (by celltype)
## acc
gen.acc.c.lasso.ts <- rbind(A_stats$V$bacc[,-13], A_stats$VI$bacc[,-13]) %>% t() %>% round(digits = 2); colnames(gen.acc.c.lasso.ts) <- c("100 c", "200 c", "300 c", "7500 hvg")

library("RColorBrewer")
cols <- c(brewer.pal(6, "Set2")[1], brewer.pal(6, "Purples")[c(2,4)], brewer.pal(6, "Set2")[4], brewer.pal(9, "Oranges")[3:8], brewer.pal(6, "Set2")[6], brewer.pal(8, "Dark2")[4])

fname <- paste0(folder, "figures/", name, "A_stats_CellGenTsConsLasso", ".png")
png(filename = fname, width = 600, height = 500, units = "px") #  pointsize = 13
par(mar=c(4,6,5,3) + 0.1)
bp <- barplot(gen.nbeta.c.lasso.ts, beside=TRUE, col = cols, ylim=c(0,80))
text(bp, gen.nbeta.c.lasso.ts+2, labels = gen.acc.c.lasso.ts, cex=0.5) 
mtext(side=3, line=3, cex=1.3, "Sparsity and Predictive performance (Lasso)")
mtext(side=3, line=1.25, cex=1.1, "Tissue-specific Myeloid cells")
mtext(side=2, line=4, cex=1.1, "Number of active coefficients (components or genes)", col="black")
legend("topleft", legend=str_replace_all(rownames(gen.nbeta.c.lasso.ts), "_", " "), fill=cols, inset=0.01, cex=0.85) # bty="n",
dev.off()


# mac-consensus:
## nbeta
mac.nbeta.c.lasso <- rbind(A_stats$Ia$nbetas, A_stats$IIIa$nbetas) %>% t(); colnames(mac.nbeta.c.lasso) <- c("100 c", "200 c", "300 c", "7500 hvg") # I+III
mac.nbeta.c.enet <- rbind(A_stats$IIa$nbetas, A_stats$IVa$nbetas) %>% t(); colnames(mac.nbeta.c.enet) <- c("100 c", "200 c", "300 c", "7500 hvg") # II+IV
## acc
mac.acc.c.lasso <- rbind(A_stats$Ia$bacc, A_stats$IIIa$bacc) %>% t() %>% round(digits = 2); colnames(mac.acc.c.lasso) <- c("100 c", "200 c", "300 c", "7500 hvg")
mac.acc.c.enet <- rbind(A_stats$IIa$bacc, A_stats$IVa$bacc) %>% t() %>% round(digits = 2); colnames(mac.acc.c.enet) <- c("100 c", "200 c", "300 c", "7500 hvg")

# change colours..
cols <- c(brewer.pal(6, "Set2")[5], brewer.pal(7, "Reds")[6], brewer.pal(6, "Purples")[4], brewer.pal(6, "Set2")[2], brewer.pal(6, "Set2")[4], brewer.pal(7, "Set2")[7])

fname <- paste0(folder, "figures/", name, "A_stats_CellMacConsLasso", ".png")
png(filename = fname, width = 600, height = 500, units = "px") #  pointsize = 13
par(mar=c(4,6,5,3) + 0.1)
bp <- barplot(mac.nbeta.c.lasso, beside=TRUE, col = cols, ylim=c(0,100))
text(bp, mac.nbeta.c.lasso+2, labels = mac.acc.c.lasso, cex=0.7) 
mtext(side=3, line=3, cex=1.3, "Sparsity and Predictive performance (Lasso)")
mtext(side=3, line=1.25, cex=1.1, "Macrophages")
mtext(side=2, line=4, cex=1.1, "Number of active coefficients (components or genes)", col="black")
legend("topleft", legend=str_replace_all(rownames(mac.nbeta.dg.lasso), "_", " "), fill=cols, inset=0.01) # bty="n",
dev.off()

# no need to also plot elastic net ones, as 1 to 0.95 change did not result in any different result 
# time permitting try alpha = 0.9 or 0.8














##### Find time to go from A to data frame heatmap #####















