# Hierarchical fastICA 
# P. Secchi, S. Vantini, and P. Zanini (2014). Hierarchical Independent Component Analysis: a multi- resolution non-orthogonal data-driven basis

library(fastHICA)

exp.mat.sumF <- readRDS("data/processed/tm.smartseq2.sumF.mat.hvg.123.rds")
exp.mat.sumF <- as.matrix(exp.mat.sumF)
dim(exp.mat.sumF)

basis <- basis_hica(X=exp.mat.sumF, maxlev=18415, dim.subset=1024)

energy <- energy_hica(HICA.obj, maxcomp = 1, nlevel = 1, plot = FALSE)

extract18319 <- extract_hica(energy.obj=energy, comp=100, level=18319)

# Unfortunately computation never finished - to slow/demanding for this data size. 
# You can try after reducing number of samples by averaging gene exp per 10 cells !
