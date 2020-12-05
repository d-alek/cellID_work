### chose_ncomp.R scipt plots two types of plots that help in deciding on the optimal number of components to select for ica runs

### INPUTS are several outputs from the icafast_extractor.r/icafast_extractor_job.sh R script and sbatch script for the cluster:
### 1. Silhouette scores (global scores per ncomp, and those per consensus component) extracted when clustering ica components 
### across 100 runs to make consensus componets, they are reflection of the stability (reproducibility) of components; 
### 2. Reconstitution errors per ncomp, they are an idication of how much information is being captured by the model, compared 
### to the original data (Error = 1 - Fit); 3. "Variance accounted for" for consensus and "diagonal" components, and null vafs 

### OUTPUTS are two type so plots: Stability vs. Error and Scree-like vs. Stability


library("here")
library("tidyverse")

spikeF_silh_stats <- readRDS(file = here("results", "spikeF_silh_stats.rds"))
spikeF_reconst_e <- readRDS(file = here("results", "spikeF_reconst_e.rds"))
spikeF_vafs <- readRDS(file = here("results", "spikeF_vafs.rds"))

sumF_silh_stats <- readRDS(file = here("results", "sumF_silh_stats.rds"))
sumF_reconst_e <- readRDS(file = here("results", "sumF_reconst_e.rds"))
sumF_vafs <- readRDS(file = here("results", "sumF_vafs.rds"))



#----------- StabilityVsError PLOT -----------------------#

pdf(file = "figures/StabilityVsError.pdf", width = 11, height = 5, pointsize = 10.5)
par(mfrow = c(1,2))

nComps <- c(50,100,150,200,250)
for(norm in c("spikeF", "sumF")){
  
  silhStat <- get(paste0(norm, "_silh_stats"))
  silhStat <- silhStat[[1]]
  reconstErr <- get(paste0(norm, "_reconst_e"))
  reconstErr_c <- reconstErr[[1]]
  reconstErr_d <- reconstErr[[2]]
  
  title <- paste0("Stability vs. Reconstitution error (", norm, " norm.)")
  # add extra space to right margin of plot
  par(mar=c(6, 4, 4, 6) + 0.1)
  # silhouette stats - left
  plot(nComps, silhStat, ylim=c(0,1), axes=FALSE, xlab="", ylab="", 
       pch=16, type="b",col="dodgerblue3", main=title)
  mtext("Component stability  (Mean silhouette width)", side=2, col="dodgerblue3", line=2.9)
  axis(2, ylim=c(0,1), col="dodgerblue3", col.axis="dodgerblue3", lwd=1.3, las=1)  ## las=1 for horizontal labels
  box(lwd=0.75)
  ## add consensus & diagonal reconst. error - right
  par(new=TRUE)
  plot(nComps, reconstErr_c, ylim=c(0.4,1), axes=FALSE, xlab="", ylab="",  
       pch=16, type="b", col="red3")
  lines(nComps, reconstErr_d, pch=21, type="b", lty = 5, col="tomato")
  mtext("Reconstitution error  (normalised)", side=4, col="red3", line=3) 
  axis(4, ylim=c(0.4,0.7), col="red3",col.axis="red3", lwd=1.2, las=1)
  # x axis
  axis(1, at = nComps)
  mtext("Number of components", side=1, col="black", line=2.5)  
  # legend
  legend("topright", legend=c("Component stability","Reconstitution error (consensus)","Reconstitution error (diagonal)"),
         text.col=c("dodgerblue3","red3","tomato"), pch=c(16,16,21), col=c("dodgerblue3","red3","tomato"), 
         cex=0.8, box.lwd=0.75, inset=0.01)
  
}
dev.off()



#----------- ScreeVsStability (consensus) PLOT -----------#

pdf(file = "figures/ScreeVsStability_c.pdf", width = 7, height = 12.5, pointsize = 8)
par(mfcol = c(5,2))

for(norm in c("spikeF", "sumF")) {
  
  vafs <- get(paste0(norm, "_vafs"))
  vaf_c <- vafs$vaf_c
  null_vaf <- vafs$null_vaf
  silhStat <- get(paste0(norm, "_silh_stats"))
  silhStat <- silhStat[[2]]
  
  for(ncomp in c(50,100,150,200,250)) {
    
    comps <- 1:ncomp
    i <- ncomp/50
    
    vaf <- vaf_c[[i]]
    vaf0 <- null_vaf[[i]]
    silh <- sort(silhStat[[i]], decreasing=TRUE)
    
    title <- paste0(ncomp, " components (", norm, " norm.)")
    par(mar = c(4, 5, 3, 5)) # default c(5.1, 4.1, 4.1, 2.1)
    plot(x=comps, y=vaf, type = "h", lwd = 2, xlab="", ylab="", axes=FALSE, ylim=c(0,0.092), main = title) 
    lines(x=comps[-1], y=vaf0[-1], col = "deeppink3", lwd = 2)
    axis(2, ylim=c(0,0.092), col="deeppink3", col.axis="deeppink3", lwd=1.2) # no space here for las=1 
    mtext("Variance accounted for (data & null)", side=2, col="deeppink3", line=2.5) 
    axis(1)
    mtext("Component number", side=1, line = 2.5) 
    box(lwd=0.5)
    par(new=TRUE)
    plot(x=comps, y=silh, type = "p", pch=16, col="dodgerblue3", xlab="", ylab="", axes=FALSE, ylim=c(-0.3,1))
    axis(4, ylim=c(-0.3,1), col="dodgerblue3", col.axis="dodgerblue3", lwd=1.3)
    mtext("Component stability (Silh. width)", side=4, col="dodgerblue3", line=2.5)
  }
}
dev.off()



#----------- ScreeVsStability (diagonal) PLOT ------------#

vaf_d <- vafs$vaf_d
# nothing new to learn with this one, so no need for plotting
# consensus and diagonal behave very similar for 50 and 100 ncomp, not after that 

#---------------------------------------------------------#

# ncomp = 100 seems like the best choice 

#---------------------------------------------------------#
