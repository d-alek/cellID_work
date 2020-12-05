### modules_by_genes.R 

library("tidyverse")
library("dplyr")
library("magrittr")
library("stringr")
library("ggplot2")
library("cowplot")
library("RColorBrewer")
library("SingleCellExperiment") # SCE
library("scran")
library("scater")
library("here")

library("clusterProfiler")
library("enrichplot")
library("ReactomePA")
library("msigdbr")
library("topGO")
library("igraph")


# here working with sumF normalised data only!

# import consensus ncomp=100 S matrices for spikeF and sumF versions (extracted in icafast_extractor.r/.sh)
sumF.S <- readRDS(file = here("results", "sumF_S_c.rds"))
sumF.S <- sumF.S[["100"]]  

# you will need CSE objects too for entrezids (for pathway analyses)
sumF.sc <- readRDS(file = here("data", "processed", "tm.smartseq2.sumF.hvg.123.icaDims.rds"))

# from consensus S matrix (100 components in 100 columns) make a list of 100, 
# each containing: S values, Gene Symbol, Enrezgene id, Ensemble gene id
genenames <- tibble(Symbol = rowData(sumF.sc)$Symbol, entrez = rowData(sumF.sc)$entrez, ensembl = rowData(sumF.sc)$ensembl) 
all.equal(rownames(sumF.S), genenames$Symbol)

# get S sd-s
S.sds <- apply(sumF.S, 2 , var) %>% sqrt() # as expected all are = 1 
# get S mad-s
S.mads <- apply(sumF.S, 2 , mad) # 0.26-0.66iga

S.lst <- list(all=list(), sd2=list(), sd3=list(), sd4=list(), mad4=list(), mad5=list(), mad6=list())
for(i in 1:100) { 
  i <- as.character(i)
  S.lst$all[[i]] <- genenames
  S.lst$all[[i]]$S <- sumF.S[, i]
  S.lst$all[[i]]$mad <- S.mads[i]
}

# make module genes lists filtered at 2, 3, and 4 sd
S.lst$sd2 <- purrr::map(S.lst$all, .f = dplyr::filter, abs(S) >= 2)
S.lst$sd3 <- purrr::map(S.lst$all, .f = dplyr::filter, abs(S) >= 3)
S.lst$sd4 <- purrr::map(S.lst$all, .f = dplyr::filter, abs(S) >= 4)
# make module genes lists filtered at 3, 4, and 5 mad - when assesed on plots below this approach is less optimal
# (this approach adapts to within-component variability - unlike above there's different mad for every component)
S.lst$mad4 <- purrr::map(S.lst$all, .f = dplyr::filter, abs(S) >= 4*mad)
S.lst$mad5 <- purrr::map(S.lst$all, .f = dplyr::filter, abs(S) >= 5*mad)
S.lst$mad6 <- purrr::map(S.lst$all, .f = dplyr::filter, abs(S) >= 6*mad)

saveRDS(S.lst, "results/sumF_S_modules.rds")
# S.lst <- readRDS(file = here("results", "sumF_S_modules.rds"))

###################################################################################################

## Histograms of module sizes, weights and gene reuses - based on sd cutoffs ###
pdf(file = "figures/Module.sizes.weight.reuse.SD.pdf", width = 9, height = 7) 
par(mfrow=c(3,3))
for(sd in c("sd2", "sd3", "sd4")) {
  sd.name <- str_replace(sd, "(sd)([234])", "\\2 \\1")
  
  uniq <- purrr::map(S.lst[[sd]], .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length() 
  reus <- round(sum(sapply(S.lst[[sd]], nrow)) / uniq, 1)
  
  # size
  hist(sapply(S.lst[[sd]], nrow), breaks = seq(0, 300, by = 5), col = "tomato1" , xlab = "Module size", cex.lab = 1, main =paste("Size of modules - at", sd.name), cex.main = 1.2, ylim = c(0,22))
  text(200, 18, paste0(uniq, "  (", reus, "x", ")"), cex = 1.1, col = "dodgerblue3", font = 2)
  
  # weight
  purrr::map(S.lst[[sd]], .f = dplyr::select, S) %>% purrr::map(.f = function(x) x^2) %>% purrr::map(.f = sum) %>% unlist() %>% (function(x) x / 4748) %>% 
    hist(breaks = seq(0, 1, by = 0.02), ylim = c(0,22), col = "thistle" , xlab = "Proportion of component weight captured", cex.lab = 1, main = paste("Module weight - at", sd.name), cex.main = 1.2)
  text(0.8, 18, paste0(uniq, "  (", reus, "x", ")"), cex = 1.1, col = "dodgerblue3", font = 2)
  # e.g. proportion of component's total "weight" i.e. variance captured by the module genes 
  # sum(S.lst[["all"]][["1"]]$S^2) # 4748 - comonent's weigh is always 4748 = number of genes 
  # sum(S.lst[["sd3"]][["1"]]$S^2)/4748 # 64%
  
  # gene re-use
  purrr::map(S.lst[[sd]], .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% table() %>% sort(decreasing = T) %>% table() %>% 
    barplot(xlab = "number of modules", ylab = "number of unique genes", cex.lab = 1, main = paste("Gene reuse - at", sd.name), cex.main = 1.2, ylim = c(0,1300), xlim = c(1,35), col = "seagreen")
  text(20, 1000, paste0(uniq, "  (", reus, "x", ")"), cex = 1.1, col = "dodgerblue3", font = 2)
  
}
dev.off()


## Histograms of module sizes, weights and gene reuses - based on mad cutoffs ###
pdf(file = "figures/Module.sizes.weight.reuse.MAD.pdf", width = 9, height = 7) 
par(mfrow=c(3,3))
for(mad in c("mad4", "mad5", "mad6")) {
  mad.name <- str_replace(mad, "(mad)([456])", "\\2 \\1")
  
  uniq <- purrr::map(S.lst[[mad]], .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length() 
  reus <- round(sum(sapply(S.lst[[mad]], nrow)) / uniq, 1)
  
  # size
  hist(sapply(S.lst[[mad]], nrow), breaks = seq(0, 440, by = 10), col = "tomato1" , xlab = "Module size", cex.lab = 1, main =paste("Size of modules - at", mad.name), cex.main = 1.2, ylim = c(0,22))
  text(200, 18, paste0(uniq, "  (", reus, "x", ")"), cex = 1.1, col = "dodgerblue3", font = 2)
  
  # weight
  purrr::map(S.lst[[mad]], .f = dplyr::select, S) %>% purrr::map(.f = function(x) x^2) %>% purrr::map(.f = sum) %>% unlist() %>% (function(x) x / 4748) %>% 
    hist(breaks = seq(0, 1, by = 0.02), ylim = c(0,22), col = "thistle" , xlab = "Proportion of component weight captured", cex.lab = 1, main = paste("Module weight - at", mad.name), cex.main = 1.2)
  text(0.8, 18, paste0(uniq, "  (", reus, "x", ")"), cex = 1.1, col = "dodgerblue3", font = 2)
  # e.g. proportion of component's total "weight" i.e. variance captured by the module genes 
  # sum(S.lst[["all"]][["1"]]$S^2) # 4748 - comonent's weigh is always 4748 = number of genes 
  # sum(S.lst[["sd3"]][["1"]]$S^2)/4748 # 64%
  
  # gene re-use
  purrr::map(S.lst[[mad]], .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% table() %>% sort(decreasing = T) %>% table() %>% 
    barplot(xlab = "number of modules", ylab = "number of unique genes", cex.lab = 1, main = paste("Gene reuse - at", mad.name), cex.main = 1.2, ylim = c(0,1300), xlim = c(1,35), col = "seagreen")
  text(20, 1000, paste0(uniq, "  (", reus, "x", ")"), cex = 1.1, col = "dodgerblue3", font = 2)
  
}
dev.off()

#-------------------------------------------------------------------------------------------------#

# Create by-celltype gene expression (median, mad) table for all hv. genes 
sumF.GE <- logcounts(sumF.sc) %>% t() %>% as.matrix() %>% as.data.frame() %>% rownames_to_column(var="cell")
cell.annos <- cbind(cell=sumF.sc$cell, celltype=sumF.sc$cell_ontology_class, tissue=sumF.sc$tissue) %>% as.data.frame(stringsAsFactors=F) %>% as_tibble()
sumF.GE <- left_join(sumF.GE, cell.annos, by="cell") %>% dplyr::select(cell, celltype, tissue, Atraid:Zyx)
# then can make varius summaries.. or have them ready:
GE.med <- sumF.GE %>% dplyr::select(-cell) %>% group_by(celltype) %>% summarise_all(.funs=median)
GE.mad <- sumF.GE %>% dplyr::select(-cell) %>% group_by(celltype) %>% summarise_all(.funs=mad)
# you can also use this as a basis for per celltype, per gene bubble plots (dot plots) - or just use existing functions from scater for this
# (colour = median expression; size = proportion of cells expressing, or size = inverse of MAD)


###################################################################################################

# Create and save Module Gene Expression dotplots - for all 100 modules - for each analysis level (Myelo.1-6, Endo.1-2)

# (done with 4 sd cutoff here, othervise gene lists would become too long for visualisation)
# (move to 2 or 3 sd cutoff when doing pairwise correlation anlysis to capture broader biology behind module)


# define myeloid and endothelial types 
myelo.names <- c("basophil", "brain pericyte", "classical monocyte", "granulocyte", "granulocyte monocyte progenitor cell", 
                 "granulocytopoietic cell", "Kupffer cell", "macrophage", "microglial cell", "monocyte")
endo.names <- c("endothelial cell", "endothelial cell of hepatic sinusoid", "lung endothelial cell")

# set cell groups for various levels of analysis below (this is aligned with SCE object cell order - all.equal(cell.grouping.GE$cell, sumF.sc$cell))
cell.grouping.GE <- sumF.GE %>% 
  dplyr::mutate(tissue = str_replace_all(tissue, "_", " ")) %>% 
  dplyr::mutate(tissue = if_else(str_detect(tissue, "^Brain "), "Brain", tissue)) %>% 
  tidyr::unite(celltype, tissue, col = celltype_tissue, sep = "-", remove = F) %>% 
  dplyr::mutate(myelo_celltype = if_else(celltype %in% myelo.names, celltype, "non-myeloid")) %>% 
  dplyr::mutate(myelo_celltype = if_else(myelo_celltype == "classical monocyte", "monocyte", myelo_celltype)) %>% 
  dplyr::mutate(myelo_celltype = if_else(myelo_celltype == "granulocyte monocyte progenitor cell", "GMP", myelo_celltype)) %>% 
  dplyr::mutate(macrophage = if_else(myelo_celltype != "macrophage", "other", myelo_celltype)) %>% 
  tidyr::unite(myelo_celltype, tissue, col = myelo_celltype_tissue, sep = "-", remove = F) %>% 
  dplyr::mutate(macrophage_tissue = if_else(str_detect(myelo_celltype_tissue, "^macrophage-") | 
                                              str_detect(myelo_celltype_tissue, "^Kupffer cell-") | 
                                              str_detect(myelo_celltype_tissue, "^microglial cell-"), myelo_celltype_tissue, "other")) %>% 
  dplyr::mutate(monocyte_tissue = if_else(str_detect(myelo_celltype_tissue, "^monocyte-"), myelo_celltype_tissue, "other")) %>% 
  dplyr::mutate(monocyte_subtype = if_else(str_detect(celltype_tissue, "^monocyte-Lung") | 
                                             str_detect(celltype_tissue, "^classical monocyte-Lung"), celltype_tissue, "other")) %>% 
  dplyr::mutate(endo_celltype = if_else(celltype %in% endo.names, "endothelial cell", "non-endothelial")) %>% 
  tidyr::unite(endo_celltype, tissue, col = endo_celltype_tissue, sep = "-", remove = F) %>% 
  dplyr::mutate(endo_celltype_tissue = if_else(str_detect(endo_celltype_tissue, "^non-endothelial-"), "non-endothelial", endo_celltype_tissue)) %>% 
  dplyr::select(cell, celltype, tissue, celltype_tissue, myelo_celltype, macrophage, myelo_celltype_tissue, macrophage_tissue, 
                monocyte_tissue, monocyte_subtype, endo_celltype, endo_celltype_tissue, Atraid:Zyx)


# Actually it would be usefull to attach some of this celltype groupings info to colData
all.equal(sumF.sc$cell, cell.grouping.GE$cell) # TRUE !
sumF.sc$myelo_celltype <- cell.grouping.GE$myelo_celltype
sumF.sc$myelo_celltype_tissue <- cell.grouping.GE$myelo_celltype_tissue
sumF.sc$macrophage <- cell.grouping.GE$macrophage
sumF.sc$macrophage_tissue <- cell.grouping.GE$macrophage_tissue
sumF.sc$monocyte_tissue <- cell.grouping.GE$monocyte_tissue
sumF.sc$monocyte_subtype <- cell.grouping.GE$monocyte_subtype
sumF.sc$endo_celltype <- cell.grouping.GE$endo_celltype
sumF.sc$endo_celltype_tissue <- cell.grouping.GE$endo_celltype_tissue

saveRDS(sumF.sc, "data/processed/tm.smartseq2.sumF.hvg.123.icaDims.rds")
# sumF.sc <- readRDS(file = here("data", "processed", "tm.smartseq2.sumF.hvg.123.icaDims.rds"))

#-----------------------------------------------------------------------------#
### Myeloid 1. & 2. (myeloid cell types + others)
#-----------------------------------------------------------------------------#

# set cell groupings and labels (the same info is now available from within SCE object - see above)
cell.grouping <- cell.grouping.GE$myelo_celltype 

cell.order <- unique(cell.grouping)[c(4,6,7,10,5,2,8,3,9,1)]
cell.names <- str_replace_all(cell.order, "granulocytopoietic cell", "granulocytopoietic\ncell")

# create GE dotplots per module
for (i in as.character(1:100)) { 
  
  module.genes <- S.lst$sd4[[i]] %>% arrange(S) %>% dplyr::select(Symbol) %>% pull()
  negs <- S.lst$sd4[[i]] %>% dplyr::filter(S<0) %>% nrow()
  poss <- S.lst$sd4[[i]] %>% dplyr::filter(S>0) %>% nrow()
  
  # calculare mean expressions for cell groups in focus, for genes in module in focus
  # extract max and sd from these: for the plot decide max_ave = mean.GE$p99 (99th percentile)
  mean.GE <- cell.grouping.GE %>% 
    dplyr::select(myelo_celltype, module.genes) %>% 
    group_by(myelo_celltype) %>% 
    summarise_all(mean) %>% 
    dplyr::select(-myelo_celltype) %>% 
    gather() %>% 
    summarise(max=max(value), p99=quantile(value, 0.99), sd=sd(value))
  
  
  w <- length(cell.names)*0.6+1 # [in]
  h <- length(module.genes)*0.175+1 # [in]
  fname <- paste0("figures/module_GE_dots/Myelo.1/", "Myelo.1.module.", i, ".pdf")
  pdf(file = fname, width = w, height = h) 
  p <- plotDots(sumF.sc, features=module.genes, group = I(cell.grouping), exprs_values = "logcounts",
                detection_limit = 0, low_color = "white", high_color = "royalblue4",
                max_ave = mean.GE$p99, max_detected = NULL, other_fields = list()) + 
    theme(axis.text.x = element_text(size=11, angle=60, hjust=1, vjust=1), 
          axis.text.y = element_text(size=11, colour=c(rep("red3", negs), rep("navy", poss))), 
          legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour = NA, fill = "transparent")) + 
    scale_y_discrete(paste("Module", i), limits=module.genes) + scale_x_discrete(NULL, limits=cell.order, labels=cell.names)
  print(p)
  dev.off()
  
}


#-----------------------------------------------------------------------------#
### Myeloid 3. (macrohages vs. ALL other, inc. myeloid non-macs) (binomial)
#-----------------------------------------------------------------------------#

# set cell groupings and labels
cell.grouping <- cell.grouping.GE$macrophage 

cell.order <- unique(cell.grouping)[c(2,1)]
cell.names <- cell.order

# create GE dotplots per module
for (i in as.character(1:100)) { 
  
  module.genes <- S.lst$sd4[[i]] %>% arrange(S) %>% dplyr::select(Symbol) %>% pull()
  negs <- S.lst$sd4[[i]] %>% dplyr::filter(S<0) %>% nrow()
  poss <- S.lst$sd4[[i]] %>% dplyr::filter(S>0) %>% nrow()
  
  # calculare mean expressions for cell groups in focus, for genes in module in focus
  # extract max and sd from these: for the plot decide max_ave = mean.GE$p99 (99th percentile)
  mean.GE <- cell.grouping.GE %>% 
    dplyr::select(macrophage, module.genes) %>% 
    group_by(macrophage) %>% 
    summarise_all(mean) %>% 
    dplyr::select(-macrophage) %>% 
    gather() %>% 
    summarise(max=max(value), p99=quantile(value, 0.99), sd=sd(value))
  
  
  w <- length(cell.names)*0.7+2 # [in]
  h <- length(module.genes)*0.175+1 # [in]
  fname <- paste0("figures/module_GE_dots/Myelo.3/", "Myelo.3.module.", i, ".pdf")
  pdf(file = fname, width = w, height = h) 
  p <- plotDots(sumF.sc, features=module.genes, group = I(cell.grouping), exprs_values = "logcounts",
                detection_limit = 0, low_color = "white", high_color = "royalblue4",
                max_ave = mean.GE$p99, max_detected = NULL, other_fields = list()) + 
    theme(axis.text.x = element_text(size=11, angle=60, hjust=1, vjust=1), 
          axis.text.y = element_text(size=11, colour=c(rep("red3", negs), rep("navy", poss))), 
          legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour = NA, fill = "transparent")) + 
    scale_y_discrete(paste("Module", i), limits=module.genes) + scale_x_discrete(NULL, limits=cell.order, labels=cell.names)
  print(p)
  dev.off()
  
}


#-----------------------------------------------------------------------------#
### Myeloid 4. (macrophages across tissues only, inc. Kupfer & microglia too)
#-----------------------------------------------------------------------------#

# set cell groupings and labels
cell.grouping <- cell.grouping.GE$macrophage_tissue 

cell.order <- unique(cell.grouping)[c(5,6,7,2,8,4,3,1)]
cell.names <- str_replace_all(cell.order, "-", "\n")

# create GE dotplots per module
for (i in as.character(1:100)) { 
  
  module.genes <- S.lst$sd4[[i]] %>% arrange(S) %>% dplyr::select(Symbol) %>% pull()
  negs <- S.lst$sd4[[i]] %>% dplyr::filter(S<0) %>% nrow()
  poss <- S.lst$sd4[[i]] %>% dplyr::filter(S>0) %>% nrow()
  
  # calculare mean expressions for cell groups in focus, for genes in module in focus
  # extract max and sd from these: for the plot decide max_ave = mean.GE$p99 (99th percentile)
  mean.GE <- cell.grouping.GE %>% 
    dplyr::select(macrophage_tissue, module.genes) %>% 
    group_by(macrophage_tissue) %>% 
    summarise_all(mean) %>% 
    dplyr::select(-macrophage_tissue) %>% 
    gather() %>% 
    summarise(max=max(value), p99=quantile(value, 0.99), sd=sd(value))
  
  
  w <- length(cell.names)*0.65+1 # [in]
  h <- length(module.genes)*0.175+1 # [in]
  fname <- paste0("figures/module_GE_dots/Myelo.4/", "Myelo.4.module.", i, ".pdf")
  pdf(file = fname, width = w, height = h) 
  p <- plotDots(sumF.sc, features=module.genes, group = I(cell.grouping), exprs_values = "logcounts",
                detection_limit = 0, low_color = "white", high_color = "royalblue4",
                max_ave = mean.GE$p99, max_detected = NULL, other_fields = list()) + 
    theme(axis.text.x = element_text(size=11, angle=60, hjust=1, vjust=1), 
          axis.text.y = element_text(size=11, colour=c(rep("red3", negs), rep("navy", poss))), 
          legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour = NA, fill = "transparent")) + 
    scale_y_discrete(paste("Module", i), limits=module.genes) + scale_x_discrete(NULL, limits=cell.order, labels=cell.names)
  print(p)
  dev.off()
  
}


#-----------------------------------------------------------------------------#
### Myeloid 5. (monocytes across two tissues only) (binomial)
#-----------------------------------------------------------------------------#

# set cell groupings and labels
cell.grouping <- cell.grouping.GE$monocyte_tissue 

cell.order <- unique(cell.grouping)[c(3,2,1)]
cell.names <- str_replace_all(cell.order, "-", "\n")

# create GE dotplots per module
for (i in as.character(1:100)) { 
  
  module.genes <- S.lst$sd4[[i]] %>% arrange(S) %>% dplyr::select(Symbol) %>% pull()
  negs <- S.lst$sd4[[i]] %>% dplyr::filter(S<0) %>% nrow()
  poss <- S.lst$sd4[[i]] %>% dplyr::filter(S>0) %>% nrow()
  
  # calculare mean expressions for cell groups in focus, for genes in module in focus
  # extract max and sd from these: for the plot decide max_ave = mean.GE$p99 (99th percentile)
  mean.GE <- cell.grouping.GE %>% 
    dplyr::select(monocyte_tissue, module.genes) %>% 
    group_by(monocyte_tissue) %>% 
    summarise_all(mean) %>% 
    dplyr::select(-monocyte_tissue) %>% 
    gather() %>% 
    summarise(max=max(value), p99=quantile(value, 0.99), sd=sd(value))
  
  
  w <- length(cell.names)*0.7+2 # [in]
  h <- length(module.genes)*0.175+1 # [in]
  fname <- paste0("figures/module_GE_dots/Myelo.5/", "Myelo.5.module.", i, ".pdf")
  pdf(file = fname, width = w, height = h) 
  p <- plotDots(sumF.sc, features=module.genes, group = I(cell.grouping), exprs_values = "logcounts",
                detection_limit = 0, low_color = "white", high_color = "royalblue4",
                max_ave = mean.GE$p99, max_detected = NULL, other_fields = list()) + 
    theme(axis.text.x = element_text(size=11, angle=60, hjust=1, vjust=1), 
          axis.text.y = element_text(size=11, colour=c(rep("red3", negs), rep("navy", poss))), 
          legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour = NA, fill = "transparent")) + 
    scale_y_discrete(paste("Module", i), limits=module.genes) + scale_x_discrete(NULL, limits=cell.order, labels=cell.names)
  print(p)
  dev.off()
  
}


#-----------------------------------------------------------------------------#
### Myeloid 6. (classic monocytes vs. monocytes in Lung) (binomial)
#-----------------------------------------------------------------------------#

# set cell groupings and labels
cell.grouping <- cell.grouping.GE$monocyte_subtype 

cell.order <- unique(cell.grouping)[c(3,2,1)]
cell.names <- str_replace_all(cell.order, "classical monocyte-Lung", "classical\nmonocyte Lung") %>% str_replace_all("monocyte-Lung", "monocyte\nLung")

# create GE dotplots per module
for (i in as.character(1:100)) { 
  
  module.genes <- S.lst$sd4[[i]] %>% arrange(S) %>% dplyr::select(Symbol) %>% pull()
  negs <- S.lst$sd4[[i]] %>% dplyr::filter(S<0) %>% nrow()
  poss <- S.lst$sd4[[i]] %>% dplyr::filter(S>0) %>% nrow()
  
  # calculare mean expressions for cell groups in focus, for genes in module in focus
  # extract max and sd from these: for the plot decide max_ave = mean.GE$p99 (99th percentile)
  mean.GE <- cell.grouping.GE %>% 
    dplyr::select(monocyte_subtype, module.genes) %>% 
    group_by(monocyte_subtype) %>% 
    summarise_all(mean) %>% 
    dplyr::select(-monocyte_subtype) %>% 
    gather() %>% 
    summarise(max=max(value), p99=quantile(value, 0.99), sd=sd(value))
  
  
  w <- length(cell.names)*0.7+2 # [in]
  h <- length(module.genes)*0.18+1 # [in]
  fname <- paste0("figures/module_GE_dots/Myelo.6/", "Myelo.6.module.", i, ".pdf")
  pdf(file = fname, width = w, height = h) 
  p <- plotDots(sumF.sc, features=module.genes, group = I(cell.grouping), exprs_values = "logcounts",
                detection_limit = 0, low_color = "white", high_color = "royalblue4",
                max_ave = mean.GE$p99, max_detected = NULL, other_fields = list()) + 
    theme(axis.text.x = element_text(size=11, angle=60, hjust=1, vjust=1), 
          axis.text.y = element_text(size=11, colour=c(rep("red3", negs), rep("navy", poss))), 
          legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour = NA, fill = "transparent")) + 
    scale_y_discrete(paste("Module", i), limits=module.genes) + scale_x_discrete(NULL, limits=cell.order, labels=cell.names)
  print(p)
  dev.off()
  
}


#-----------------------------------------------------------------------------#
### Endothelial 1. (endothelial cells + others) (binomial)
#-----------------------------------------------------------------------------#

# set cell groupings and labels
cell.grouping <- cell.grouping.GE$endo_celltype 

cell.order <- unique(cell.grouping)
cell.names <- cell.order

# create GE dotplots per module
for (i in as.character(1:100)) { 
  
  module.genes <- S.lst$sd4[[i]] %>% arrange(S) %>% dplyr::select(Symbol) %>% pull()
  negs <- S.lst$sd4[[i]] %>% dplyr::filter(S<0) %>% nrow()
  poss <- S.lst$sd4[[i]] %>% dplyr::filter(S>0) %>% nrow()
  
  # calculare mean expressions for cell groups in focus, for genes in module in focus
  # extract max and sd from these: for the plot decide max_ave = mean.GE$p99 (99th percentile)
  mean.GE <- cell.grouping.GE %>% 
    dplyr::select(endo_celltype, module.genes) %>% 
    group_by(endo_celltype) %>% 
    summarise_all(mean) %>% 
    dplyr::select(-endo_celltype) %>% 
    gather() %>% 
    summarise(max=max(value), p99=quantile(value, 0.99), sd=sd(value))
  
  
  w <- length(cell.names)*0.7+2 # [in]
  h <- length(module.genes)*0.18+1 # [in]
  fname <- paste0("figures/module_GE_dots/Endo.1/", "Endo.1.module.", i, ".pdf")
  pdf(file = fname, width = w, height = h) 
  p <- plotDots(sumF.sc, features=module.genes, group = I(cell.grouping), exprs_values = "logcounts",
                detection_limit = 0, low_color = "white", high_color = "royalblue4",
                max_ave = mean.GE$p99, max_detected = NULL, other_fields = list()) + 
    theme(axis.text.x = element_text(size=11, angle=60, hjust=1, vjust=1), 
          axis.text.y = element_text(size=11, colour=c(rep("red3", negs), rep("navy", poss))), 
          legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour = NA, fill = "transparent")) + 
    scale_y_discrete(paste("Module", i), limits=module.genes) + scale_x_discrete(NULL, limits=cell.order, labels=cell.names)
  print(p)
  dev.off()
  
}


#-----------------------------------------------------------------------------#
### Endothelial 2. (endothelial cells across tissues only)
#-----------------------------------------------------------------------------#

# set cell groupings and labels
cell.grouping <- cell.grouping.GE$endo_celltype_tissue 

cell.order <- unique(cell.grouping)[c(2,7,6,8,9,5,1,10,11,3)]
cell.names <- str_replace_all(cell.order, "-", "\n")

# create GE dotplots per module
for (i in as.character(1:100)) { 
  
  module.genes <- S.lst$sd4[[i]] %>% arrange(S) %>% dplyr::select(Symbol) %>% pull()
  negs <- S.lst$sd4[[i]] %>% dplyr::filter(S<0) %>% nrow()
  poss <- S.lst$sd4[[i]] %>% dplyr::filter(S>0) %>% nrow()
  
  # calculare mean expressions for cell groups in focus, for genes in module in focus
  # extract max and sd from these: for the plot decide max_ave = mean.GE$p99 (99th percentile)
  mean.GE <- cell.grouping.GE %>% 
    dplyr::select(endo_celltype_tissue, module.genes) %>% 
    group_by(endo_celltype_tissue) %>% 
    summarise_all(mean) %>% 
    dplyr::select(-endo_celltype_tissue) %>% 
    gather() %>% 
    summarise(max=max(value), p99=quantile(value, 0.99), sd=sd(value))
  
  
  w <- length(cell.names)*0.7+1 # [in]
  h <- length(module.genes)*0.18+1 # [in]
  fname <- paste0("figures/module_GE_dots/Endo.2/", "Endo.2.module.", i, ".pdf")
  pdf(file = fname, width = w, height = h) 
  p <- plotDots(sumF.sc, features=module.genes, group = I(cell.grouping), exprs_values = "logcounts",
                detection_limit = 0, low_color = "white", high_color = "royalblue4",
                max_ave = mean.GE$p99, max_detected = NULL, other_fields = list()) + 
    theme(axis.text.x = element_text(size=11, angle=60, hjust=1, vjust=1), 
          axis.text.y = element_text(size=11, colour=c(rep("red3", negs), rep("navy", poss))), 
          legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour = NA, fill = "transparent")) + 
    scale_y_discrete(paste("Module", i), limits=module.genes) + scale_x_discrete(NULL, limits=cell.order, labels=cell.names)
  print(p)
  dev.off()
  
}


#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#


###################################################################################################





###################################################################################################

# Explore Gene-Set analysis in several ways GO, Reactome pathway, Motifs enrichment.. 
# the first two are just to show groupings of genes in biologically coherent sets
# motfs work is to establish potential regulators behind module gene co-expression


### Prep data for Over-representation type analysis at 3 sd - split to + and - (alt. all in one if too small numbers in -)
#-------------------------------------------------------------------------------------------------#
# Create (+)-Major and (-)-minor modules based on the sign of S score 
S.lst <- readRDS(file = here("results", "sumF_S_modules.rds"))
S.lst <- S.lst[c("all", "sd2", "sd3", "sd4")] # mad-based modules didn't look so good 

Ma.names <- paste(1:100, rep("ma", 100), sep="_")
Mi.names <- paste(1:100, rep("mi", 100), sep="_")

# apply S +/- cutoffs to deeper levels of the list (each module)
S.lst.Ma <- S.lst %>% purrr::map_depth(.depth = 2, .f = dplyr::filter, S > 0) %>% purrr::map(.f = setNames, Ma.names)
S.lst.Mi <- S.lst %>% purrr::map_depth(.depth = 2, .f = dplyr::filter, S < 0) %>% purrr::map(.f = setNames, Mi.names)

# Unite Major (+) and minor (-) modules to create lists of 200 at each sd level
S.lst.Ma.Mi <- purrr::map2(S.lst.Ma, S.lst.Mi, .f = append)
saveRDS(S.lst.Ma.Mi, "results/sumF_S_modules_Ma_Mi.rds")
# S.lst.Ma.Mi <- readRDS(file = here("results", "sumF_S_modules_Ma_Mi.rds"))

# Extract just entezids for gene set analyses (exclude NAs)
S.lst.Ma.Mi.entezid <- S.lst.Ma.Mi %>% purrr::map_depth(.depth = 2, .f = dplyr::pull, entrez)


### Prep data for GSEA type analysis on all ~ 4750 genes scored by S and ranked + to -
#-------------------------------------------------------------------------------------------------#
# for this need ot create named vector of S scores, where names are entrezids, then sort them descendingly
S.lst.S <- S.lst %>% purrr::map_depth(.depth = 2, .f = dplyr::select, S) %>% purrr::map_depth(.depth = 2, .f = dplyr::pull, S)
S.lst.entezid <- S.lst %>% purrr::map_depth(.depth = 2, .f = dplyr::select, entrez) %>% purrr::map_depth(.depth = 2, .f = dplyr::pull, entrez)

for (level in c("all", "sd2", "sd3", "sd4")) {
  for (i in 1:100) {
    names(S.lst.S[[level]][[i]]) <- as.character(S.lst.entezid[[level]][[i]])
  }
}
# rename accordingly and sort decreasingly for GSEA
S.lst.entezid.named <- S.lst.S %>% purrr::map_depth(.depth = 2, .f = sort, decreasing = T)



########## Reactome ##########

## Over-Represenation Analyses (ORA) 
# with split +/- modules (sd3)
ora.reactome.sd3.ma.mi <- compareCluster(S.lst.Ma.Mi.entezid$sd3, fun="enrichPathway", organism="mouse", pvalueCutoff=0.1,
                                         pAdjustMethod="BH", qvalueCutoff=0.25, universe, minGSSize=5, maxGSSize=1000, readable=T)

# with bundled +/- modules (sd3)
ora.reactome.sd3.bundle <- compareCluster(S.lst.entezid$sd3, fun="enrichPathway", organism="mouse", pvalueCutoff=0.1,
                                          pAdjustMethod="BH", qvalueCutoff=0.25, universe, minGSSize=5, maxGSSize=1000, readable=T)


## GSEA (Gene set Enrichment Analyses) (all)
gsea.reactome.all <- purrr::map(.x=S.lst.entezid.named[["all"]], .f=gsePathway, 
                                organism="mouse", exponent=1, nPerm=100000,
                                minGSSize = 5, maxGSSize = 1000, pvalueCutoff = 0.1,
                                pAdjustMethod="BH", verbose=TRUE, seed=FALSE, by="fgsea") # need to make "readable"? 
gsea.reactome.all <- purrr::map(.x=gsea.reactome.all, .f=setReadable, OrgDb='org.Mm.eg.db', keyType="ENTREZID")


saveRDS(ora.reactome.sd3.ma.mi, "results/ora_gsea/reactome_ora_ma_mi.rds")
saveRDS(ora.reactome.sd3.bundle, "results/ora_gsea/reactome_ora_bundle.rds")
saveRDS(gsea.reactome.all, "results/ora_gsea/reactome_gsea.rds")
# test
colnames(ora.reactome.sd3.ma.mi@compareClusterResult) # 10: Cluster, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count
colnames(gsea.reactome.all$`52`@result) # 11: ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalues, rank, leading_edge, core_enrichment



############  GO  ############

## Over-Represenation Analyses (ORA) 
# with split +/- modules (sd3)
ora.GO.sd3.ma.mi <- compareCluster(S.lst.Ma.Mi.entezid$sd3, fun="enrichGO", OrgDb='org.Mm.eg.db', keyType="ENTREZID", ont="BP", 
                                      pvalueCutoff=0.1, pAdjustMethod="BH", universe, qvalueCutoff=0.25, minGSSize=5, maxGSSize=1000, 
                                      readable=TRUE, pool=FALSE)

# with bundled +/- modules (sd3)
ora.GO.sd3.bundle <- compareCluster(S.lst.entezid$sd3, fun="enrichGO", OrgDb='org.Mm.eg.db', keyType="ENTREZID", ont="BP", 
                                    pvalueCutoff=0.1, pAdjustMethod="BH", universe, qvalueCutoff=0.25, minGSSize=5, maxGSSize=1000, 
                                    readable=TRUE, pool=FALSE)

# to potentially simplify GO.. either run "groupGO" insted, or run simplify() or filter() or dropGO() on the uotput above


## GSEA (Gene set Enrichment Analyses) (all)
gsea.GO.all <- purrr::map(.x=S.lst.entezid.named[["all"]], .f=gseGO, 
                          OrgDb='org.Mm.eg.db', keyType="ENTREZID", ont="BP", exponent=1, nPerm=100000,
                          minGSSize = 5, maxGSSize = 1000, pvalueCutoff = 0.1,
                          pAdjustMethod="BH", verbose=TRUE, seed=FALSE, by="fgsea") # need to make "readable"? 
gsea.GO.all <- purrr::map(.x=gsea.GO.all, .f=setReadable, OrgDb='org.Mm.eg.db', keyType="ENTREZID")


saveRDS(ora.GO.sd3.ma.mi, "results/ora_gsea/GO_ora_ma_mi.rds")
saveRDS(ora.GO.sd3.bundle, "results/ora_gsea/GO_ora_bundle.rds")
saveRDS(gsea.GO.all, "results/ora_gsea/GO_gsea.rds")
# test
colnames(ora.GO.sd3.ma.mi@compareClusterResult) # 11: Cluster, ONTOLOGY, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count
colnames(gsea.GO.all$`36`@result) # 12: ONTOLOGY, ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalues, rank, leading_edge, core_enrichment


######### TF Motifs ###########

# prep gmt file (can use cran msigdbr human-mouse or BioC gskb mouse package):
# (abandoned gskb approach as it focuses on smaller experimentally validated targets, so very biased set)

#--------MSigDb-based---------#
mm.msigdb.c3.tf = msigdbr(species="Mus musculus", category="C3", subcategory="TFT") # 153,310 x 9 tibble
colnames(mm.msigdb.c3.tf)

mm.msigdb.c3.tf <- mm.msigdb.c3.tf %>% dplyr::select(gs_name, entrez_gene, gene_symbol) %>% as.data.frame()


## Over-Represenation Analyses (ORA) 
# with split +/- modules (sd3)
ora.TF.sd3.ma.mi <- map(.x=S.lst.Ma.Mi.entezid$sd3, .f = enricher, TERM2GENE=mm.msigdb.c3.tf[,1:2], 
                              pvalueCutoff=0.1, pAdjustMethod="BH", universe, minGSSize=5, maxGSSize=100000, qvalueCutoff=0.25)
# must remove NULLs before next step:
ora.TF.sd3.ma.mi[["16_mi"]] = NULL; ora.TF.sd3.ma.mi[["61_mi"]] = NULL; ora.TF.sd3.ma.mi[["64_mi"]] = NULL
ora.TF.sd3.ma.mi <- purrr::map(.x=ora.TF.sd3.ma.mi, .f=setReadable, OrgDb='org.Mm.eg.db', keyType="ENTREZID")

# with bundled +/- modules (sd3)
ora.TF.sd3.bundle <- map(.x=S.lst.entezid$sd3, .f = enricher, TERM2GENE=mm.msigdb.c3.tf[,1:2], 
                               pvalueCutoff=0.1, pAdjustMethod="BH", universe, minGSSize=5, maxGSSize=100000, qvalueCutoff=0.25)
ora.TF.sd3.bundle <- purrr::map(.x=ora.TF.sd3.bundle, .f=setReadable, OrgDb='org.Mm.eg.db', keyType="ENTREZID")


## GSEA (Gene set Enrichment Analyses) (all)
gsea.TF.all <- map(.x=S.lst.entezid.named[["all"]], .f=GSEA, TERM2GENE=mm.msigdb.c3.tf[,1:2], 
                         exponent=1, nPerm=100000, minGSSize = 5, maxGSSize = 1000, pvalueCutoff = 0.1, 
                         pAdjustMethod="BH", verbose=TRUE, seed=FALSE, by="fgsea") # need to make "readable"? 
gsea.TF.all <- purrr::map(.x=gsea.TF.all, .f=setReadable, OrgDb='org.Mm.eg.db', keyType="ENTREZID")


saveRDS(ora.TF.sd3.ma.mi, "results/ora_gsea/TF_ora_ma_mi.rds")
saveRDS(ora.TF.sd3.bundle, "results/ora_gsea/TF_ora_bundle.rds")
saveRDS(gsea.TF.all, "results/ora_gsea/TF_gsea.rds")
# test
colnames(ora.TF.sd3.ma.mi$`52_ma`@result) # 9: ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count
colnames(gsea.TF.all$`36`@result) # 11: ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalues, rank, leading_edge, core_enrichment


#-------------------------------------------------------------------------------------------------#


# NOW PLOT SOME OF THIS ORA / GSEA RESULTS !!!

ora.reactome.sd3.ma.mi <- readRDS("results/ora_gsea/reactome_ora_ma_mi.rds")
ora.reactome.sd3.bundle <- readRDS("results/ora_gsea/reactome_ora_bundle.rds")
gsea.reactome.all <- readRDS("results/ora_gsea/reactome_gsea.rds")

ora.GO.sd3.ma.mi <- readRDS("results/ora_gsea/GO_ora_ma_mi.rds")
ora.GO.sd3.bundle <- readRDS("results/ora_gsea/GO_ora_bundle.rds")
gsea.GO.all <- readRDS("results/ora_gsea/GO_gsea.rds")

S.lst <- readRDS(file = here("results", "sumF_S_modules.rds"))


################# Summaries of Reactome Pathway analyses ######################
#-----------------------------------------------------------------------------#

### 1. Some global paterns modules (64,6,51,5,35,37,41,42)
#-----------------------------------------------------------------------------#
# ORA with bundled modules - Comparative Dotplot:
# Subset ORA and GSEA CompareClusters objects as you dont want to visualise all modules at once (would be too messy):
ora.reactome.bundle.subst <- ora.reactome.sd3.bundle
subst <- ora.reactome.sd3.bundle@compareClusterResult %>% as_tibble() %>% dplyr::filter(Cluster %in% c(64,6,51,5,35,37,41,42)) %>% as.data.frame()
ora.reactome.bundle.subst@compareClusterResult <- subst
pdf(file = "figures/module_ora_gsea_tf/global.reactome.ora.bundle.pdf", width = 9, height = 11)
dotplot(ora.reactome.bundle.subst, 
        color = "p.adjust", showCategory = 10, split = NULL,
        font.size = 10, title = "Global pattern modules - Reactiome ORA", by = "geneRatio", includeAll = TRUE) + 
  scale_color_gradient(low = "red2", high = "white")
dev.off()

# ORA with +/- split modules - Comparative Dotplot:
ora.reactome.ma.mi.subst <- ora.reactome.sd3.ma.mi
mod.names <- c(str_c(c(64,6,51,5,35,37,41,42), "_ma"), str_c(c(64,6,51,5,35,37,41,42), "_mi"))
subst <- ora.reactome.sd3.ma.mi@compareClusterResult %>% as_tibble() %>% 
  dplyr::filter(Cluster %in% mod.names) %>% 
  dplyr::mutate(Cluster = str_replace_all(Cluster, "_ma", "M")) %>% dplyr::mutate(Cluster = str_replace_all(Cluster, "_mi", "m")) %>% as.data.frame()
ora.reactome.ma.mi.subst@compareClusterResult <- subst
pdf(file = "figures/module_ora_gsea_tf/global.reactome.ora.ma.mi.pdf", width = 12, height = 12)
dotplot(ora.reactome.ma.mi.subst, 
        color = "p.adjust", showCategory = 6, split = NULL,
        font.size = 10, title = "Global pattern modules - Reactiome ORA", by = "geneRatio", includeAll = TRUE) + 
  scale_color_gradient(low = "red2", high = "white")
dev.off()

# e.g. of some other ORA plotting summaries on one module (module 6 here)
S.entrezid.3sd.6 <- S.lst$sd3$`6`$entrez
S.score.3sd.6 <- S.lst$sd3$`6`$S
names(S.score.3sd.6) <- S.entrezid.3sd.6
ora.reactome.6 <- enrichPathway(S.entrezid.3sd.6, organism="mouse", pvalueCutoff=0.1,
                                pAdjustMethod="BH", qvalueCutoff=0.25, minGSSize=5, maxGSSize=1000, readable=T)
# emap
pdf(file = "figures/module_ora_gsea_tf/global.reactome.ora.emap.6.pdf", width = 10, height = 10)
emapplot(ora.reactome.6, pie_scale=1.5, layout="kk") + scale_color_gradient(low = "red2", high = "white")
dev.off()

# heatmap
pdf(file = "figures/module_ora_gsea_tf/global.reactome.ora.heatmap.6.pdf", width = 12, height = 7)
heatplot(ora.reactome.6, foldChange=S.score.3sd.6) + 
  scale_fill_gradient2("S score", low="blue", mid="white", high="red2")
dev.off()

# gene-concept network
pdf(file = "figures/module_ora_gsea_tf/global.reactome.ora.cnet.6.pdf", width = 13, height = 6)
p1 <- cnetplot(ora.reactome.6, showCategory=7, node_label="category", layout="kk") + 
  theme(legend.position="none")
p2 <- cnetplot(ora.reactome.6, showCategory=7, node_label="gene", foldChange=S.score.3sd.6, layout="kk") + 
  scale_colour_gradient2("S score", low="blue", mid="white", high="red2")
cowplot::plot_grid(p1, p2, ncol=2, labels=c("Pathways", "Genes"), scale=0.9) 
dev.off()


# GSEA - Running Score plots: no clusterprofiler for this so make a loop for all modules and select three top gene sets for each to plot 
#       (or if they are very overlapping in terms of biology, select any sensible three from the top to display)
for (i in as.character(c(64,6,51,5,35,37,41,42))) {
  
  module.gsea <- gsea.reactome.all[[i]]
  top3 <- gsea.reactome.all[[i]]@result %>% as_tibble() %>% arrange(desc(NES)) %>% as.data.frame()
  module.gsea@result <- top3
  p1 <- gsearank(module.gsea, geneSetID=1, title = module.gsea[1, "Description"]) + xlab(NULL) + annotate("text", x=3500, y=module.gsea[1, "enrichmentScore"]*0.75, label=paste0("NES = ", round(module.gsea[1,"NES"],2)), hjust=0, vjust=0, colour="red3")
  p2 <- gsearank(module.gsea, geneSetID=2, title = module.gsea[2, "Description"]) + xlab(NULL) + annotate("text", x=3500, y=module.gsea[2, "enrichmentScore"]*0.75, label=paste0("NES = ", round(module.gsea[2,"NES"],2)), hjust=0, vjust=0, colour="red3")
  p3 <- gsearank(module.gsea, geneSetID=3, title = module.gsea[3, "Description"]) + annotate("text", x=3500, y=module.gsea[3, "enrichmentScore"]*0.75, label=paste0("NES = ", round(module.gsea[3,"NES"],2)), hjust=0, vjust=0, colour="red3")
  pdf(file = paste0("figures/module_ora_gsea_tf/gsea/global.reactome.gsea.", i, ".pdf"), width = 6, height = 7)
  p <- plot_grid(p1,p2,p3, ncol=1)
  print(p)
  dev.off()
  
}

# e.g. of some other GSEA plotting summaries on one module (module 6 here)
S.entrezid.all.6 <- S.lst$all$`6`$entrez
S.score.all.6 <- S.lst$all$`6`$S
names(S.score.all.6) <- S.entrezid.all.6
# emap
pdf(file = "figures/module_ora_gsea_tf/global.reactome.gsea.emap.6.pdf", width = 9, height = 9)
emapplot(gsea.reactome.all$`6`, showCategory=25, pie_scale=1.5, layout="kk") + 
  scale_color_gradient(low = "red2", high = "white")
dev.off()

# heatmap
pdf(file = "figures/module_ora_gsea_tf/global.reactome.gsea.heatmap.6.pdf", width = 15, height = 7)
heatplot(gsea.reactome.all$`6`, foldChange=S.score.all.6) + 
  scale_fill_gradient2("S score", low="blue", mid="white", high="red2") + 
  scale_x_discrete(breaks=NULL) + xlab("Genes")
dev.off()

# gene-concept network
pdf(file = "figures/module_ora_gsea_tf/global.reactome.gsea.cnet.6.pdf", width = 14, height = 6)
p1 <- cnetplot(gsea.reactome.all$`6`, showCategory=7, node_label="category", layout="kk") + 
  theme(legend.position="none")
p2 <- cnetplot(gsea.reactome.all$`6`, showCategory=7, node_label="gene", foldChange=S.score.all.6, layout="kk") + 
  scale_colour_gradient2("S score", low="blue", mid="white", high="red2")
cowplot::plot_grid(p1, p2, ncol=2, labels=c("Pathways", "Genes"), scale=0.9) 
dev.off()



### 2. Main myeloid, B cell, pericyte type modules (15,81,72,52,40,70,71,63?,12,21,23,36,29,30,79)
#-----------------------------------------------------------------------------#
# ORA with bundled modules - Comparative Dotplot:
# Subset ORA and GSEA CompareClusters objects as you dont want to visualise all modules at once (would be too messy):
ora.reactome.bundle.subst <- ora.reactome.sd3.bundle
subst <- ora.reactome.sd3.bundle@compareClusterResult %>% as_tibble() %>% dplyr::filter(Cluster %in% c(15,81,72,52,40,70,71,63,12,21,23,36,29,30,79)) %>% as.data.frame()
ora.reactome.bundle.subst@compareClusterResult <- subst
pdf(file = "figures/module_ora_gsea_tf/blood.reactome.ora.bundle.pdf", width = 14, height = 11)
dotplot(ora.reactome.bundle.subst, 
        color = "p.adjust", showCategory = 10, split = NULL,
        font.size = 10, title = "Blood cells modules - Reactiome ORA", by = "geneRatio", includeAll = TRUE) + 
  scale_color_gradient(low = "red2", high = "white")
dev.off()

# ORA with +/- split modules - Comparative Dotplot:
ora.reactome.ma.mi.subst <- ora.reactome.sd3.ma.mi
mod.names <- c(str_c(c(15,81,72,52,40,70,71,63,12,21,23,36,29,30,79), "_ma"), str_c(c(15,81,72,52,40,70,71,63,12,21,23,36,29,30,79), "_mi"))
subst <- ora.reactome.sd3.ma.mi@compareClusterResult %>% as_tibble() %>% 
  dplyr::filter(Cluster %in% mod.names) %>% 
  dplyr::mutate(Cluster = str_replace_all(Cluster, "_ma", "M")) %>% dplyr::mutate(Cluster = str_replace_all(Cluster, "_mi", "m")) %>% as.data.frame()
ora.reactome.ma.mi.subst@compareClusterResult <- subst
pdf(file = "figures/module_ora_gsea_tf/blood.reactome.ora.ma.mi.pdf", width = 17.5, height = 12.5)
dotplot(ora.reactome.ma.mi.subst, 
        color = "p.adjust", showCategory = 6, split = NULL,
        font.size = 10, title = "Blood cells modules - Reactiome ORA", by = "geneRatio", includeAll = TRUE) + 
  scale_color_gradient(low = "red2", high = "white")
dev.off()

# e.g. of some other ORA plotting summaries on one module (module 52 here)
S.entrezid.3sd.52 <- S.lst$sd3$`52`$entrez
S.score.3sd.52 <- S.lst$sd3$`52`$S
names(S.score.3sd.52) <- S.entrezid.3sd.52
ora.reactome.52 <- enrichPathway(S.entrezid.3sd.52, organism="mouse", pvalueCutoff=0.1,
                                pAdjustMethod="BH", qvalueCutoff=0.25, minGSSize=5, maxGSSize=1000, readable=T)
# emap
pdf(file = "figures/module_ora_gsea_tf/blood.reactome.ora.emap.52.pdf", width = 8, height = 8)
emapplot(ora.reactome.52, pie_scale=1.5, layout="kk") + scale_color_gradient(low = "red2", high = "white")
dev.off()

# heatmap
pdf(file = "figures/module_ora_gsea_tf/blood.reactome.ora.heatmap.52.pdf", width = 10, height = 4)
heatplot(ora.reactome.52, foldChange=S.score.3sd.52) + 
  scale_fill_gradient2("S score", low="blue", mid="white", high="red2") + 
  theme(axis.text.x=element_text(size=8))
dev.off()

# gene-concept network
pdf(file = "figures/module_ora_gsea_tf/blood.reactome.ora.cnet.52.pdf", width = 13, height = 6)
p1 <- cnetplot(ora.reactome.52, showCategory=7, node_label="category", layout="kk") + 
  theme(legend.position="none")
p2 <- cnetplot(ora.reactome.52, showCategory=7, node_label="gene", foldChange=S.score.3sd.52, layout="kk") + 
  scale_colour_gradient2("S score", low="blue", mid="white", high="red2")
cowplot::plot_grid(p1, p2, ncol=2, labels=c("Pathways", "Genes"), scale=0.9) 
dev.off()


# GSEA - Running Score plots: no clusterprofiler for this so make a loop for all modules and select three top gene sets for each to plot 
#       (or if they are very overlapping in terms of biology, select any sensible three from the top to display)
for (i in as.character(c(15,72,52,71,12,21,23,29,79))) {
  
  module.gsea <- gsea.reactome.all[[i]]
  top3 <- gsea.reactome.all[[i]]@result %>% as_tibble() %>% arrange(desc(NES)) %>% as.data.frame()
  module.gsea@result <- top3
  p1 <- gsearank(module.gsea, geneSetID=1, title = module.gsea[1, "Description"]) + xlab(NULL) + annotate("text", x=3000, y=module.gsea[1, "enrichmentScore"]*0.75, label=paste0("NES = ", round(module.gsea[1,"NES"],2)), hjust=0, vjust=0, colour="red3")
  p2 <- gsearank(module.gsea, geneSetID=2, title = module.gsea[2, "Description"]) + xlab(NULL) + annotate("text", x=3000, y=module.gsea[2, "enrichmentScore"]*0.75, label=paste0("NES = ", round(module.gsea[2,"NES"],2)), hjust=0, vjust=0, colour="red3")
  p3 <- gsearank(module.gsea, geneSetID=3, title = module.gsea[3, "Description"]) + annotate("text", x=3000, y=module.gsea[3, "enrichmentScore"]*0.75, label=paste0("NES = ", round(module.gsea[3,"NES"],2)), hjust=0, vjust=0, colour="red3")
  pdf(file = paste0("figures/module_ora_gsea_tf/gsea/blood.reactome.gsea.", i, ".pdf"), width = 6, height = 7)
  p <- plot_grid(p1,p2,p3, ncol=1)
  print(p)
  dev.off()
  
} # 81, 40, 70, 63, 36, 30 are causing errors... ?

# e.g. of some other GSEA plotting summaries on one module (module 52 here)
S.entrezid.all.52 <- S.lst$all$`52`$entrez
S.score.all.52 <- S.lst$all$`52`$S
names(S.score.all.52) <- S.entrezid.all.52
# emap
pdf(file = "figures/module_ora_gsea_tf/blood.reactome.gsea.emap.52.pdf", width = 4, height = 4)
emapplot(gsea.reactome.all$`52`, showCategory=25, pie_scale=1.5, layout="kk") + 
  scale_color_gradient(low = "red2", high = "white")
dev.off()

# heatmap
pdf(file = "figures/module_ora_gsea_tf/blood.reactome.gsea.heatmap.52.pdf", width = 15, height = 3)
heatplot(gsea.reactome.all$`52`, foldChange=S.score.all.52) + 
  scale_fill_gradient2("S score", low="blue", mid="white", high="red2") + 
  scale_x_discrete(breaks=NULL) + xlab("Genes")
dev.off()

# gene-concept network
pdf(file = "figures/module_ora_gsea_tf/blood.reactome.gsea.cnet.52.pdf", width = 14, height = 6)
p1 <- cnetplot(gsea.reactome.all$`52`, showCategory=4, node_label="category", layout="kk") + 
  theme(legend.position="none")
p2 <- cnetplot(gsea.reactome.all$`52`, showCategory=4, node_label="gene", foldChange=S.score.all.52, layout="kk") + 
  scale_colour_gradient2("S score", low="blue", mid="white", high="red2")
cowplot::plot_grid(p1, p2, ncol=2, labels=c("Pathways", "Genes"), scale=0.9) 
dev.off()


### 3. Main Endothelial cells modules (9,66,67,82,89,69,61,94)
#-----------------------------------------------------------------------------#
# ORA with bundled modules - Comparative Dotplot:
# Subset ORA and GSEA CompareClusters objects as you dont want to visualise all modules at once (would be too messy):
ora.reactome.bundle.subst <- ora.reactome.sd3.bundle
subst <- ora.reactome.sd3.bundle@compareClusterResult %>% as_tibble() %>% dplyr::filter(Cluster %in% c(9,66,67,82,89,69,61,94)) %>% as.data.frame()
ora.reactome.bundle.subst@compareClusterResult <- subst
pdf(file = "figures/module_ora_gsea_tf/endothel.reactome.ora.bundle.pdf", width = 12, height = 9)
dotplot(ora.reactome.bundle.subst, 
        color = "p.adjust", showCategory = 10, split = NULL,
        font.size = 10, title = "Endothelial cells modules - Reactiome ORA", by = "geneRatio", includeAll = TRUE) + 
  scale_color_gradient(low = "red2", high = "white")
dev.off()

# ORA with +/- split modules - Comparative Dotplot:
ora.reactome.ma.mi.subst <- ora.reactome.sd3.ma.mi
mod.names <- c(str_c(c(9,66,67,82,89,69,61,94), "_ma"), str_c(c(9,66,67,82,89,69,61,94), "_mi"))
subst <- ora.reactome.sd3.ma.mi@compareClusterResult %>% as_tibble() %>% 
  dplyr::filter(Cluster %in% mod.names) %>% 
  dplyr::mutate(Cluster = str_replace_all(Cluster, "_ma", "M")) %>% dplyr::mutate(Cluster = str_replace_all(Cluster, "_mi", "m")) %>% as.data.frame()
ora.reactome.ma.mi.subst@compareClusterResult <- subst
pdf(file = "figures/module_ora_gsea_tf/endothel.reactome.ora.ma.mi.pdf", width = 14, height = 10)
dotplot(ora.reactome.ma.mi.subst, 
        color = "p.adjust", showCategory = 6, split = NULL,
        font.size = 10, title = "Endothelial cells modules - Reactiome ORA", by = "geneRatio", includeAll = TRUE) + 
  scale_color_gradient(low = "red2", high = "white")
dev.off()

# e.g. of some other ORA plotting summaries on one module (module 67 here)
S.entrezid.3sd.67 <- S.lst$sd3$`67`$entrez
S.score.3sd.67 <- S.lst$sd3$`67`$S
names(S.score.3sd.67) <- S.entrezid.3sd.67
ora.reactome.67 <- enrichPathway(S.entrezid.3sd.67, organism="mouse", pvalueCutoff=0.1,
                                 pAdjustMethod="BH", qvalueCutoff=0.25, minGSSize=5, maxGSSize=1000, readable=T)
# emap
pdf(file = "figures/module_ora_gsea_tf/endothel.reactome.ora.emap.67.pdf", width = 8, height = 8)
emapplot(ora.reactome.67, pie_scale=1.5, layout="kk") + scale_color_gradient(low = "red2", high = "white")
dev.off()

# heatmap
pdf(file = "figures/module_ora_gsea_tf/endothel.reactome.ora.heatmap.67.pdf", width = 10, height = 5)
heatplot(ora.reactome.67, foldChange=S.score.3sd.67) + 
  scale_fill_gradient2("S score", low="blue", mid="white", high="red2") + 
  theme(axis.text.x=element_text(size=8))
dev.off()

# gene-concept network
pdf(file = "figures/module_ora_gsea_tf/endothel.reactome.ora.cnet.67.pdf", width = 13, height = 6)
p1 <- cnetplot(ora.reactome.67, showCategory=7, node_label="category", layout="kk") + 
  theme(legend.position="none")
p2 <- cnetplot(ora.reactome.67, showCategory=7, node_label="gene", foldChange=S.score.3sd.67, layout="kk") + 
  scale_colour_gradient2("S score", low="blue", mid="white", high="red2")
cowplot::plot_grid(p1, p2, ncol=2, labels=c("Pathways", "Genes"), scale=0.9) 
dev.off()


# GSEA - Running Score plots: no clusterprofiler for this so make a loop for all modules and select three top gene sets for each to plot 
#       (or if they are very overlapping in terms of biology, select any sensible three from the top to display)
for (i in as.character(c(67,69,61))) {
  
  module.gsea <- gsea.reactome.all[[i]]
  top3 <- gsea.reactome.all[[i]]@result %>% as_tibble() %>% arrange(desc(NES)) %>% as.data.frame()
  module.gsea@result <- top3
  p1 <- gsearank(module.gsea, geneSetID=1, title = module.gsea[1, "Description"]) + xlab(NULL) + annotate("text", x=2500, y=module.gsea[1, "enrichmentScore"]*0.75, label=paste0("NES = ", round(module.gsea[1,"NES"],2)), hjust=0, vjust=0, colour="red3")
  p2 <- gsearank(module.gsea, geneSetID=2, title = module.gsea[2, "Description"]) + xlab(NULL) + annotate("text", x=2500, y=module.gsea[2, "enrichmentScore"]*0.75, label=paste0("NES = ", round(module.gsea[2,"NES"],2)), hjust=0, vjust=0, colour="red3")
  p3 <- gsearank(module.gsea, geneSetID=3, title = module.gsea[3, "Description"]) + annotate("text", x=2500, y=module.gsea[3, "enrichmentScore"]*0.75, label=paste0("NES = ", round(module.gsea[3,"NES"],2)), hjust=0, vjust=0, colour="red3")
  pdf(file = paste0("figures/module_ora_gsea_tf/gsea/endothel.reactome.gsea.", i, ".pdf"), width = 6, height = 7)
  p <- plot_grid(p1,p2,p3, ncol=1)
  print(p)
  dev.off()
  
} # 9, 66, 82, 89, 94 are causing errors... ?

# e.g. of some other GSEA plotting summaries on one module (module 67 here)
S.entrezid.all.67 <- S.lst$all$`67`$entrez
S.score.all.67 <- S.lst$all$`67`$S
names(S.score.all.67) <- S.entrezid.all.67
# emap
pdf(file = "figures/module_ora_gsea_tf/endothel.reactome.gsea.emap.67.pdf", width = 6, height = 6)
emapplot(gsea.reactome.all$`67`, showCategory=25, pie_scale=1.5, layout="kk") + 
  scale_color_gradient(low = "red2", high = "white")
dev.off()

# heatmap
pdf(file = "figures/module_ora_gsea_tf/endothel.reactome.gsea.heatmap.67.pdf", width = 8, height = 2)
heatplot(gsea.reactome.all$`67`, foldChange=S.score.all.67) + 
  scale_fill_gradient2("S score", low="blue", mid="white", high="red2") + 
  theme(axis.text.x=element_text(size=7))
dev.off()

# gene-concept network
pdf(file = "figures/module_ora_gsea_tf/endothel.reactome.gsea.cnet.67.pdf", width = 14, height = 6)
p1 <- cnetplot(gsea.reactome.all$`67`, showCategory=4, node_label="category", layout="kk") + 
  theme(legend.position="none")
p2 <- cnetplot(gsea.reactome.all$`67`, showCategory=4, node_label="gene", foldChange=S.score.all.67, layout="kk") + 
  scale_colour_gradient2("S score", low="blue", mid="white", high="red2")
cowplot::plot_grid(p1, p2, ncol=2, labels=c("Pathways", "Genes"), scale=0.9) 
dev.off()






################ Summaries of Gene Ontology (GO) analyses #####################
#-----------------------------------------------------------------------------#




# go
goplot(gsea.GO.all$`36`)






########## Summaries of Transcription Factors Motifs analyses #################
#-----------------------------------------------------------------------------#

ora.TF.sd3.ma.mi <- readRDS("results/ora_gsea/TF_ora_ma_mi.rds")
ora.TF.sd3.bundle <- readRDS("results/ora_gsea/TF_ora_bundle.rds")
gsea.TF.all <- readRDS("results/ora_gsea/TF_gsea.rds")



### 1. Some global paterns modules (64,6,51,5,35,37,41,42)
#-----------------------------------------------------------------------------#
# ORA with bundled modules - Comparative Dotplot:
# Subset ORA and GSEA CompareClusters objects as you dont want to visualise all modules at once (would be too messy):

for (i in as.character(c(c(6,51,5,35,37,41,42), 
                         c(15,81,72,52,40,70,71,63,12,21,23,36,29,30,79), 
                         c(9,66,67,82,89,69,94)))) { 
  
  ma <- paste0(i, "_ma")
  ma.name <- paste0(i, "M")
  mi <- paste0(i, "_mi")
  mi.name <- paste0(i, "m")
  
  pdf(file = paste0("figures/module_ora_gsea_tf/TF_motifs/global.tf.", i, ".pdf"), width = 22, height = 7)
  if (is.null(dim(ora.TF.sd3.ma.mi[[ma]])) | dim(ora.TF.sd3.ma.mi[[ma]])[1] == 0) { 
    p1 <- ggplot() + geom_blank() } else {
    p1 <- dotplot(ora.TF.sd3.ma.mi[[ma]], color = "p.adjust", showCategory=30, 
                  font.size=10, title=paste(ma.name, "- TF motifs ORA")) + scale_color_gradient(low = "red2", high = "white") + 
      guides(color=guide_colorbar(order=1), size=guide_legend(order=0))
  }
  if (is.null(dim(ora.TF.sd3.ma.mi[[mi]])) | dim(ora.TF.sd3.ma.mi[[mi]])[1] == 0) { 
    p2 <- ggplot() + geom_blank() } else {
    p2 <- dotplot(ora.TF.sd3.ma.mi[[mi]], color = "p.adjust", showCategory=30, 
                  font.size=10, title=paste(mi.name, "- TF motifs ORA")) + scale_color_gradient(low = "red2", high = "white") + 
      guides(color=guide_colorbar(order=1), size=guide_legend(order=0))
    }
  if (is.null(dim(ora.TF.sd3.bundle[[i]])) | dim(ora.TF.sd3.bundle[[i]])[1] == 0) { 
    p3 <- ggplot() + geom_blank() } else {
    p3 <- dotplot(ora.TF.sd3.bundle[[i]], color = "p.adjust", showCategory=30, 
                  font.size=10, title=paste(i, "- TF motifs ORA")) + scale_color_gradient(low = "red2", high = "white") + 
      guides(color=guide_colorbar(order=1), size=guide_legend(order=0))
    }
  if (is.null(dim(gsea.TF.all[[i]])) | dim(gsea.TF.all[[i]])[1] == 0) { 
    p4 <- ggplot() + geom_blank() } else {
    p4 <- dotplot(gsea.TF.all[[i]], color = "p.adjust", showCategory=30, 
                  font.size=10, title=paste(i, "- TF motifs GSEA")) + scale_color_gradient(low = "red2", high = "white") + 
      guides(color=guide_colorbar(order=1), size=guide_legend(order=0))
    }
  
  p <- plot_grid(p1,p2,p3,p4, nrow=1)
  print(p)
  dev.off()
  
  }

# 64 and 61 on its own
for (i in as.character(c(64, 61))) { 
  
  ma <- paste0(i, "_ma")
  ma.name <- paste0(i, "M")
  mi <- paste0(i, "_mi")
  mi.name <- paste0(i, "m")
  
  pdf(file = paste0("figures/module_ora_gsea_tf/TF_motifs/global.tf.", i, ".pdf"), width = 22, height = 7)
  if (is.null(dim(ora.TF.sd3.ma.mi[[ma]])) | dim(ora.TF.sd3.ma.mi[[ma]])[1] == 0) { 
    p1 <- ggplot() + geom_blank() } else {
      p1 <- dotplot(ora.TF.sd3.ma.mi[[ma]], color = "p.adjust", showCategory=30, 
                    font.size=10, title=paste(ma.name, "- TF motifs ORA")) + scale_color_gradient(low = "red2", high = "white") + 
        guides(color=guide_colorbar(order=1), size=guide_legend(order=0))
    }
  if (is.null(dim(ora.TF.sd3.ma.mi[[mi]]))) { 
    p2 <- ggplot() + geom_blank() } else {
      p2 <- dotplot(ora.TF.sd3.ma.mi[[mi]], color = "p.adjust", showCategory=30, 
                    font.size=10, title=paste(mi.name, "- TF motifs ORA")) + scale_color_gradient(low = "red2", high = "white") + 
        guides(color=guide_colorbar(order=1), size=guide_legend(order=0))
    }
  if (is.null(dim(ora.TF.sd3.bundle[[i]])) | dim(ora.TF.sd3.bundle[[i]])[1] == 0) { 
    p3 <- ggplot() + geom_blank() } else {
      p3 <- dotplot(ora.TF.sd3.bundle[[i]], color = "p.adjust", showCategory=30, 
                    font.size=10, title=paste(i, "- TF motifs ORA")) + scale_color_gradient(low = "red2", high = "white") + 
        guides(color=guide_colorbar(order=1), size=guide_legend(order=0))
    }
  if (is.null(dim(gsea.TF.all[[i]])) | dim(gsea.TF.all[[i]])[1] == 0) { 
    p4 <- ggplot() + geom_blank() } else {
      p4 <- dotplot(gsea.TF.all[[i]], color = "p.adjust", showCategory=30, 
                    font.size=10, title=paste(i, "- TF motifs GSEA")) + scale_color_gradient(low = "red2", high = "white") + 
        guides(color=guide_colorbar(order=1), size=guide_legend(order=0))
    }
  
  p <- plot_grid(p1,p2,p3,p4, nrow=1)
  print(p)
  dev.off()
  
}

# In summary this analysis should be repeated either with relaxed thresholds on adjusted p or on reduced "universe", or both !



p1 <- dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
p2 <- dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")
plot_grid(p1, p2, ncol=2)

ora.reactome.bundle.subst <- ora.reactome.sd3.bundle
subst <- ora.reactome.sd3.bundle@compareClusterResult %>% as_tibble() %>% dplyr::filter(Cluster %in% c(64,6,51,5,35,37,41,42)) %>% as.data.frame()
ora.reactome.bundle.subst@compareClusterResult <- subst
pdf(file = "figures/module_ora_gsea_tf/global.reactome.ora.bundle.pdf", width = 9, height = 11)
dotplot(ora.reactome.bundle.subst, 
        color = "p.adjust", showCategory = 10, split = NULL,
        font.size = 10, title = "Global pattern modules - Reactiome ORA", by = "geneRatio", includeAll = TRUE) + 
  scale_color_gradient(low = "red2", high = "white")
dev.off()

# ORA with +/- split modules - Comparative Dotplot:

pdf(file = "figures/module_ora_gsea_tf/global.reactome.ora.ma.mi.pdf", width = 12, height = 12)
dotplot(ora.reactome.ma.mi.subst, 
        color = "p.adjust", showCategory = 6, split = NULL,
        font.size = 10, title = "Global pattern modules - Reactiome ORA", by = "geneRatio", includeAll = TRUE) + 
  scale_color_gradient(low = "red2", high = "white")
dev.off()

# e.g. of some other ORA plotting summaries on one module (module 6 here)
S.entrezid.3sd.6 <- S.lst$sd3$`6`$entrez
S.score.3sd.6 <- S.lst$sd3$`6`$S
names(S.score.3sd.6) <- S.entrezid.3sd.6
ora.reactome.6 <- enrichPathway(S.entrezid.3sd.6, organism="mouse", pvalueCutoff=0.1,
                                pAdjustMethod="BH", qvalueCutoff=0.25, minGSSize=5, maxGSSize=1000, readable=T)
# emap
pdf(file = "figures/module_ora_gsea_tf/global.reactome.ora.emap.6.pdf", width = 10, height = 10)
emapplot(ora.reactome.6, pie_scale=1.5, layout="kk") + scale_color_gradient(low = "red2", high = "white")
dev.off()

# heatmap
pdf(file = "figures/module_ora_gsea_tf/global.reactome.ora.heatmap.6.pdf", width = 12, height = 7)
heatplot(ora.reactome.6, foldChange=S.score.3sd.6) + 
  scale_fill_gradient2("S score", low="blue", mid="white", high="red2")
dev.off()

# gene-concept network
pdf(file = "figures/module_ora_gsea_tf/global.reactome.ora.cnet.6.pdf", width = 13, height = 6)
p1 <- cnetplot(ora.reactome.6, showCategory=7, node_label="category", layout="kk") + 
  theme(legend.position="none")
p2 <- cnetplot(ora.reactome.6, showCategory=7, node_label="gene", foldChange=S.score.3sd.6, layout="kk") + 
  scale_colour_gradient2("S score", low="blue", mid="white", high="red2")
cowplot::plot_grid(p1, p2, ncol=2, labels=c("Pathways", "Genes"), scale=0.9) 
dev.off()


# GSEA - Running Score plots: no clusterprofiler for this so make a loop for all modules and select three top gene sets for each to plot 
#       (or if they are very overlapping in terms of biology, select any sensible three from the top to display)
for (i in as.character(c(64,6,51,5,35,37,41,42))) {
  
  module.gsea <- gsea.reactome.all[[i]]
  top3 <- gsea.reactome.all[[i]]@result %>% as_tibble() %>% arrange(desc(NES)) %>% as.data.frame()
  module.gsea@result <- top3
  p1 <- gsearank(module.gsea, geneSetID=1, title = module.gsea[1, "Description"]) + xlab(NULL) + annotate("text", x=3500, y=module.gsea[1, "enrichmentScore"]*0.75, label=paste0("NES = ", round(module.gsea[1,"NES"],2)), hjust=0, vjust=0, colour="red3")
  p2 <- gsearank(module.gsea, geneSetID=2, title = module.gsea[2, "Description"]) + xlab(NULL) + annotate("text", x=3500, y=module.gsea[2, "enrichmentScore"]*0.75, label=paste0("NES = ", round(module.gsea[2,"NES"],2)), hjust=0, vjust=0, colour="red3")
  p3 <- gsearank(module.gsea, geneSetID=3, title = module.gsea[3, "Description"]) + annotate("text", x=3500, y=module.gsea[3, "enrichmentScore"]*0.75, label=paste0("NES = ", round(module.gsea[3,"NES"],2)), hjust=0, vjust=0, colour="red3")
  pdf(file = paste0("figures/module_ora_gsea_tf/gsea/global.reactome.gsea.", i, ".pdf"), width = 6, height = 7)
  p <- plot_grid(p1,p2,p3, ncol=1)
  print(p)
  dev.off()
  
}

# e.g. of some other GSEA plotting summaries on one module (module 6 here)
S.entrezid.all.6 <- S.lst$all$`6`$entrez
S.score.all.6 <- S.lst$all$`6`$S
names(S.score.all.6) <- S.entrezid.all.6
# emap
pdf(file = "figures/module_ora_gsea_tf/global.reactome.gsea.emap.6.pdf", width = 9, height = 9)
emapplot(gsea.reactome.all$`6`, showCategory=25, pie_scale=1.5, layout="kk") + 
  scale_color_gradient(low = "red2", high = "white")
dev.off()

# heatmap
pdf(file = "figures/module_ora_gsea_tf/global.reactome.gsea.heatmap.6.pdf", width = 15, height = 7)
heatplot(gsea.reactome.all$`6`, foldChange=S.score.all.6) + 
  scale_fill_gradient2("S score", low="blue", mid="white", high="red2") + 
  scale_x_discrete(breaks=NULL) + xlab("Genes")
dev.off()

# gene-concept network
pdf(file = "figures/module_ora_gsea_tf/global.reactome.gsea.cnet.6.pdf", width = 14, height = 6)
p1 <- cnetplot(gsea.reactome.all$`6`, showCategory=7, node_label="category", layout="kk") + 
  theme(legend.position="none")
p2 <- cnetplot(gsea.reactome.all$`6`, showCategory=7, node_label="gene", foldChange=S.score.all.6, layout="kk") + 
  scale_colour_gradient2("S score", low="blue", mid="white", high="red2")
cowplot::plot_grid(p1, p2, ncol=2, labels=c("Pathways", "Genes"), scale=0.9) 
dev.off()



### 2. Main myeloid, B cell, pericyte type modules (15,81,72,52,40,70,71,63?,12,21,23,36,29,30,79)
#-----------------------------------------------------------------------------#
# ORA with bundled modules - Comparative Dotplot:
# Subset ORA and GSEA CompareClusters objects as you dont want to visualise all modules at once (would be too messy):

pdf(file = "figures/module_ora_gsea_tf/blood.reactome.ora.bundle.pdf", width = 14, height = 11)
dotplot(ora.reactome.bundle.subst, 
        color = "p.adjust", showCategory = 10, split = NULL,
        font.size = 10, title = "Blood cells modules - Reactiome ORA", by = "geneRatio", includeAll = TRUE) + 
  scale_color_gradient(low = "red2", high = "white")
dev.off()

# ORA with +/- split modules - Comparative Dotplot:

pdf(file = "figures/module_ora_gsea_tf/blood.reactome.ora.ma.mi.pdf", width = 17.5, height = 12.5)
dotplot(ora.reactome.ma.mi.subst, 
        color = "p.adjust", showCategory = 6, split = NULL,
        font.size = 10, title = "Blood cells modules - Reactiome ORA", by = "geneRatio", includeAll = TRUE) + 
  scale_color_gradient(low = "red2", high = "white")
dev.off()

# e.g. of some other ORA plotting summaries on one module (module 52 here)
S.entrezid.3sd.52 <- S.lst$sd3$`52`$entrez
S.score.3sd.52 <- S.lst$sd3$`52`$S
names(S.score.3sd.52) <- S.entrezid.3sd.52
ora.reactome.52 <- enrichPathway(S.entrezid.3sd.52, organism="mouse", pvalueCutoff=0.1,
                                 pAdjustMethod="BH", qvalueCutoff=0.25, minGSSize=5, maxGSSize=1000, readable=T)
# emap
pdf(file = "figures/module_ora_gsea_tf/blood.reactome.ora.emap.52.pdf", width = 8, height = 8)
emapplot(ora.reactome.52, pie_scale=1.5, layout="kk") + scale_color_gradient(low = "red2", high = "white")
dev.off()

# heatmap
pdf(file = "figures/module_ora_gsea_tf/blood.reactome.ora.heatmap.52.pdf", width = 10, height = 4)
heatplot(ora.reactome.52, foldChange=S.score.3sd.52) + 
  scale_fill_gradient2("S score", low="blue", mid="white", high="red2") + 
  theme(axis.text.x=element_text(size=8))
dev.off()

# gene-concept network
pdf(file = "figures/module_ora_gsea_tf/blood.reactome.ora.cnet.52.pdf", width = 13, height = 6)
p1 <- cnetplot(ora.reactome.52, showCategory=7, node_label="category", layout="kk") + 
  theme(legend.position="none")
p2 <- cnetplot(ora.reactome.52, showCategory=7, node_label="gene", foldChange=S.score.3sd.52, layout="kk") + 
  scale_colour_gradient2("S score", low="blue", mid="white", high="red2")
cowplot::plot_grid(p1, p2, ncol=2, labels=c("Pathways", "Genes"), scale=0.9) 
dev.off()


# GSEA - Running Score plots: no clusterprofiler for this so make a loop for all modules and select three top gene sets for each to plot 
#       (or if they are very overlapping in terms of biology, select any sensible three from the top to display)
for (i in as.character(c(15,72,52,71,12,21,23,29,79))) {
  
  module.gsea <- gsea.reactome.all[[i]]
  top3 <- gsea.reactome.all[[i]]@result %>% as_tibble() %>% arrange(desc(NES)) %>% as.data.frame()
  module.gsea@result <- top3
  p1 <- gsearank(module.gsea, geneSetID=1, title = module.gsea[1, "Description"]) + xlab(NULL) + annotate("text", x=3000, y=module.gsea[1, "enrichmentScore"]*0.75, label=paste0("NES = ", round(module.gsea[1,"NES"],2)), hjust=0, vjust=0, colour="red3")
  p2 <- gsearank(module.gsea, geneSetID=2, title = module.gsea[2, "Description"]) + xlab(NULL) + annotate("text", x=3000, y=module.gsea[2, "enrichmentScore"]*0.75, label=paste0("NES = ", round(module.gsea[2,"NES"],2)), hjust=0, vjust=0, colour="red3")
  p3 <- gsearank(module.gsea, geneSetID=3, title = module.gsea[3, "Description"]) + annotate("text", x=3000, y=module.gsea[3, "enrichmentScore"]*0.75, label=paste0("NES = ", round(module.gsea[3,"NES"],2)), hjust=0, vjust=0, colour="red3")
  pdf(file = paste0("figures/module_ora_gsea_tf/gsea/blood.reactome.gsea.", i, ".pdf"), width = 6, height = 7)
  p <- plot_grid(p1,p2,p3, ncol=1)
  print(p)
  dev.off()
  
} # 81, 40, 70, 63, 36, 30 are causing errors... ?

# e.g. of some other GSEA plotting summaries on one module (module 52 here)
S.entrezid.all.52 <- S.lst$all$`52`$entrez
S.score.all.52 <- S.lst$all$`52`$S
names(S.score.all.52) <- S.entrezid.all.52
# emap
pdf(file = "figures/module_ora_gsea_tf/blood.reactome.gsea.emap.52.pdf", width = 4, height = 4)
emapplot(gsea.reactome.all$`52`, showCategory=25, pie_scale=1.5, layout="kk") + 
  scale_color_gradient(low = "red2", high = "white")
dev.off()

# heatmap
pdf(file = "figures/module_ora_gsea_tf/blood.reactome.gsea.heatmap.52.pdf", width = 15, height = 3)
heatplot(gsea.reactome.all$`52`, foldChange=S.score.all.52) + 
  scale_fill_gradient2("S score", low="blue", mid="white", high="red2") + 
  scale_x_discrete(breaks=NULL) + xlab("Genes")
dev.off()

# gene-concept network
pdf(file = "figures/module_ora_gsea_tf/blood.reactome.gsea.cnet.52.pdf", width = 14, height = 6)
p1 <- cnetplot(gsea.reactome.all$`52`, showCategory=4, node_label="category", layout="kk") + 
  theme(legend.position="none")
p2 <- cnetplot(gsea.reactome.all$`52`, showCategory=4, node_label="gene", foldChange=S.score.all.52, layout="kk") + 
  scale_colour_gradient2("S score", low="blue", mid="white", high="red2")
cowplot::plot_grid(p1, p2, ncol=2, labels=c("Pathways", "Genes"), scale=0.9) 
dev.off()


### 3. Main Endothelial cells modules (9,66,67,82,89,69,61,94)
#-----------------------------------------------------------------------------#
# ORA with bundled modules - Comparative Dotplot:
# Subset ORA and GSEA CompareClusters objects as you dont want to visualise all modules at once (would be too messy):

pdf(file = "figures/module_ora_gsea_tf/endothel.reactome.ora.bundle.pdf", width = 12, height = 9)
dotplot(ora.reactome.bundle.subst, 
        color = "p.adjust", showCategory = 10, split = NULL,
        font.size = 10, title = "Endothelial cells modules - Reactiome ORA", by = "geneRatio", includeAll = TRUE) + 
  scale_color_gradient(low = "red2", high = "white")
dev.off()

# ORA with +/- split modules - Comparative Dotplot:

pdf(file = "figures/module_ora_gsea_tf/endothel.reactome.ora.ma.mi.pdf", width = 14, height = 10)
dotplot(ora.reactome.ma.mi.subst, 
        color = "p.adjust", showCategory = 6, split = NULL,
        font.size = 10, title = "Endothelial cells modules - Reactiome ORA", by = "geneRatio", includeAll = TRUE) + 
  scale_color_gradient(low = "red2", high = "white")
dev.off()

# e.g. of some other ORA plotting summaries on one module (module 67 here)
S.entrezid.3sd.67 <- S.lst$sd3$`67`$entrez
S.score.3sd.67 <- S.lst$sd3$`67`$S
names(S.score.3sd.67) <- S.entrezid.3sd.67
ora.reactome.67 <- enrichPathway(S.entrezid.3sd.67, organism="mouse", pvalueCutoff=0.1,
                                 pAdjustMethod="BH", qvalueCutoff=0.25, minGSSize=5, maxGSSize=1000, readable=T)
# emap
pdf(file = "figures/module_ora_gsea_tf/endothel.reactome.ora.emap.67.pdf", width = 8, height = 8)
emapplot(ora.reactome.67, pie_scale=1.5, layout="kk") + scale_color_gradient(low = "red2", high = "white")
dev.off()

# heatmap
pdf(file = "figures/module_ora_gsea_tf/endothel.reactome.ora.heatmap.67.pdf", width = 10, height = 5)
heatplot(ora.reactome.67, foldChange=S.score.3sd.67) + 
  scale_fill_gradient2("S score", low="blue", mid="white", high="red2") + 
  theme(axis.text.x=element_text(size=8))
dev.off()

# gene-concept network
pdf(file = "figures/module_ora_gsea_tf/endothel.reactome.ora.cnet.67.pdf", width = 13, height = 6)
p1 <- cnetplot(ora.reactome.67, showCategory=7, node_label="category", layout="kk") + 
  theme(legend.position="none")
p2 <- cnetplot(ora.reactome.67, showCategory=7, node_label="gene", foldChange=S.score.3sd.67, layout="kk") + 
  scale_colour_gradient2("S score", low="blue", mid="white", high="red2")
cowplot::plot_grid(p1, p2, ncol=2, labels=c("Pathways", "Genes"), scale=0.9) 
dev.off()


# GSEA - Running Score plots: no clusterprofiler for this so make a loop for all modules and select three top gene sets for each to plot 
#       (or if they are very overlapping in terms of biology, select any sensible three from the top to display)
for (i in as.character(c(67,69,61))) {
  
  module.gsea <- gsea.reactome.all[[i]]
  top3 <- gsea.reactome.all[[i]]@result %>% as_tibble() %>% arrange(desc(NES)) %>% as.data.frame()
  module.gsea@result <- top3
  p1 <- gsearank(module.gsea, geneSetID=1, title = module.gsea[1, "Description"]) + xlab(NULL) + annotate("text", x=2500, y=module.gsea[1, "enrichmentScore"]*0.75, label=paste0("NES = ", round(module.gsea[1,"NES"],2)), hjust=0, vjust=0, colour="red3")
  p2 <- gsearank(module.gsea, geneSetID=2, title = module.gsea[2, "Description"]) + xlab(NULL) + annotate("text", x=2500, y=module.gsea[2, "enrichmentScore"]*0.75, label=paste0("NES = ", round(module.gsea[2,"NES"],2)), hjust=0, vjust=0, colour="red3")
  p3 <- gsearank(module.gsea, geneSetID=3, title = module.gsea[3, "Description"]) + annotate("text", x=2500, y=module.gsea[3, "enrichmentScore"]*0.75, label=paste0("NES = ", round(module.gsea[3,"NES"],2)), hjust=0, vjust=0, colour="red3")
  pdf(file = paste0("figures/module_ora_gsea_tf/gsea/endothel.reactome.gsea.", i, ".pdf"), width = 6, height = 7)
  p <- plot_grid(p1,p2,p3, ncol=1)
  print(p)
  dev.off()
  
} # 9, 66, 82, 89, 94 are causing errors... ?

# e.g. of some other GSEA plotting summaries on one module (module 67 here)
S.entrezid.all.67 <- S.lst$all$`67`$entrez
S.score.all.67 <- S.lst$all$`67`$S
names(S.score.all.67) <- S.entrezid.all.67
# emap
pdf(file = "figures/module_ora_gsea_tf/endothel.reactome.gsea.emap.67.pdf", width = 6, height = 6)
emapplot(gsea.reactome.all$`67`, showCategory=25, pie_scale=1.5, layout="kk") + 
  scale_color_gradient(low = "red2", high = "white")
dev.off()

# heatmap
pdf(file = "figures/module_ora_gsea_tf/endothel.reactome.gsea.heatmap.67.pdf", width = 8, height = 2)
heatplot(gsea.reactome.all$`67`, foldChange=S.score.all.67) + 
  scale_fill_gradient2("S score", low="blue", mid="white", high="red2") + 
  theme(axis.text.x=element_text(size=7))
dev.off()

# gene-concept network
pdf(file = "figures/module_ora_gsea_tf/endothel.reactome.gsea.cnet.67.pdf", width = 14, height = 6)
p1 <- cnetplot(gsea.reactome.all$`67`, showCategory=4, node_label="category", layout="kk") + 
  theme(legend.position="none")
p2 <- cnetplot(gsea.reactome.all$`67`, showCategory=4, node_label="gene", foldChange=S.score.all.67, layout="kk") + 
  scale_colour_gradient2("S score", low="blue", mid="white", high="red2")
cowplot::plot_grid(p1, p2, ncol=2, labels=c("Pathways", "Genes"), scale=0.9) 
dev.off()





















# motif analysis
motif.tf.100.001 <- map(entrez.modules.ica.100.var.lfdr.001, .f = enricher, TERM2GENE = setC3tf)
motif.tf.smluni.100.001 <- map(entrez.modules.ica.100.var.lfdr.001, .f = enricher, TERM2GENE = setC3tf, universe = as.character(genenames.var$entrezgene))

motif.tf.100.001$`9.M`@result %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::select(ID)
motif.tf.100.001$`43.M`@result %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::select(ID)
motif.tf.100.001$`7.M`@result %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::select(ID)
motif.tf.100.001$`8.M`@result %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::select(ID)
motif.tf.100.001$`8.m`@result %>% dplyr::filter(p.adjust < 0.01) %>% dplyr::select(ID)

ids <- motif.tf.100.001$`9.M`@result %>% dplyr::filter(p.adjust < 0.05, ID == "RGAGGAARY_PU1_Q6") %>% dplyr::select(geneID) %>% str_split(pattern = "/") %>% unlist()
syms <- bitr(ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Mm.eg.db")
syms$SYMBOL
















###################################################################################################

# Explore correlatePairs() approach in few ways:
# A. Correlations between pairs of genes within a particular module (e.g. 52, 63, 36) 
#    - do it across all myeloid cells, across all macrophages, across all monocytes, accross endothelial cells
#    - this will create different correlation structures (visualised as networks) accross these diferent cell types 
#      so in addition to info on how much module is used in particular cell - this will add info on which parts of the module are used in certain cell types 
# B. Correlations between a particular module genes and TFs (e.g. 52, 63, 36 vs. TFs within hvg genes ) 
#    - across all myeloid cells, across all macrophages, across all monocytes, accross endothelial cells
#    - this is not to be interpreted as this TF regulates this gene (there is too much colinearity in system here for this)
#      but this will generally be able to score TF which tend to be correlated with module genes (data-driven TF exploration)
#    - this then needs to be intersected with TF motif analysis from above (knowledge-driven TF exploration) 
#      to create lists of most likely TFs behing modules in various cell types
# C. Same as 2. but with expanded list of TFs (including those not originally picked as higly-variable genes) - student project maybe...


correlation_pairs <- list(intra.module.genes = list(myeloid = list(), macrophage = list(), microglia = list(), monocyte = list(), granulocyte = list(), endothel = list()), 
                          module.genes.to.TF = list(myeloid = list(), macrophage = list(), microglia = list(), monocyte = list(), granulocyte = list(), endothel = list()) )

# 176 TFs expressed in the curent HVG set:
tf.sumF.hv <- readRDS(file = here("data", "processed", "tf_sumF_hv.rds"))


### 1. across all myeloid cells ###########################

is.myeloid <- sumF.sc$myelo_celltype %in% c("macrophage", "microglial cell", "GMP", "monocyte", 
                                            "granulocytopoietic cell", "granulocyte", "Kupffer cell", "basophil")
sce <- sumF.sc[, is.myeloid]

for (x in as.character(1:100)) { 
  
  module.x <- S.lst$sd3[[x]]$Symbol
  
  ## A. Correlations between genes within a module (intra.module.genes corr):
  # pair-wise stats
  cP <- correlatePairs(sce, subset.row=module.x, pairings=NULL)
  cP <- cP %>% as_tibble() %>% dplyr::filter(FDR <= 0.01)
  # global stats per gene 
  cG <- correlateGenes(cP)
  cG <- as_tibble(cG)
  
  # save
  correlation_pairs$intra.module.genes$myeloid[[x]][["cP"]] <- cP
  correlation_pairs$intra.module.genes$myeloid[[x]][["cG"]] <- cG
  
  
  ## B. Correlations between module genes and TFs (module.genes.to.TF corr):
  # pair-wise stats
  cP.TF <- correlatePairs(sce, subset.row=c(module.x, tf.sumF.hv), pairings=list(module.x, tf.sumF.hv))
  cP.TF <- cP.TF %>% as_tibble() %>% dplyr::filter(FDR <= 0.01)
  # global stats per gene 
  cTF <- correlateGenes(cP.TF)
  cTF <- as_tibble(cTF) %>% dplyr::filter(gene %in% tf.sumF.hv) # filters only TFs and their global scores  
  
  # save
  correlation_pairs$module.genes.to.TF$myeloid[[x]][["cP.TF"]] <- cP.TF
  correlation_pairs$module.genes.to.TF$myeloid[[x]][["cTF"]] <- cTF
  
}


### 2. across all macrophages (inc. Kupfer) ###############

is.mac <- sumF.sc$myelo_celltype %in% c("macrophage", "Kupffer cell") # 444
sce <- sumF.sc[, is.mac]

for (x in as.character(1:100)) { 
  
  module.x <- S.lst$sd3[[x]]$Symbol
  
  ## A. Correlations between genes within a module (intra.module.genes corr):
  # pair-wise stats
  cP <- correlatePairs(sce, subset.row=module.x, pairings=NULL)
  cP <- cP %>% as_tibble() %>% dplyr::filter(FDR <= 0.01)
  # global stats per gene 
  cG <- correlateGenes(cP)
  cG <- as_tibble(cG)
  
  # save
  correlation_pairs$intra.module.genes$macrophage[[x]][["cP"]] <- cP
  correlation_pairs$intra.module.genes$macrophage[[x]][["cG"]] <- cG
  
  
  ## B. Correlations between module genes and TFs (module.genes.to.TF corr):
  # pair-wise stats
  cP.TF <- correlatePairs(sce, subset.row=c(module.x, tf.sumF.hv), pairings=list(module.x, tf.sumF.hv))
  cP.TF <- cP.TF %>% as_tibble() %>% dplyr::filter(FDR <= 0.01)
  # global stats per gene 
  cTF <- correlateGenes(cP.TF)
  cTF <- as_tibble(cTF) %>% dplyr::filter(gene %in% tf.sumF.hv) # filters only TFs and their global scores  
  
  # save
  correlation_pairs$module.genes.to.TF$macrophage[[x]][["cP.TF"]] <- cP.TF
  correlation_pairs$module.genes.to.TF$macrophage[[x]][["cTF"]] <- cTF
  
}


### 3. across microglia ###################################

is.mic <- sumF.sc$myelo_celltype %in% c("microglial cell") # 600
sce <- sumF.sc[, is.mic]

for (x in as.character(1:100)) { 
  
  module.x <- S.lst$sd3[[x]]$Symbol
  
  ## A. Correlations between genes within a module (intra.module.genes corr):
  # pair-wise stats
  cP <- correlatePairs(sce, subset.row=module.x, pairings=NULL)
  cP <- cP %>% as_tibble() %>% dplyr::filter(FDR <= 0.01)
  # global stats per gene 
  cG <- correlateGenes(cP)
  cG <- as_tibble(cG)
  
  # save
  correlation_pairs$intra.module.genes$microglia[[x]][["cP"]] <- cP
  correlation_pairs$intra.module.genes$microglia[[x]][["cG"]] <- cG
  
  
  ## B. Correlations between module genes and TFs (module.genes.to.TF corr):
  # pair-wise stats
  cP.TF <- correlatePairs(sce, subset.row=c(module.x, tf.sumF.hv), pairings=list(module.x, tf.sumF.hv))
  cP.TF <- cP.TF %>% as_tibble() %>% dplyr::filter(FDR <= 0.01)
  # global stats per gene 
  cTF <- correlateGenes(cP.TF)
  cTF <- as_tibble(cTF) %>% dplyr::filter(gene %in% tf.sumF.hv) # filters only TFs and their global scores  
  
  # save
  correlation_pairs$module.genes.to.TF$microglia[[x]][["cP.TF"]] <- cP.TF
  correlation_pairs$module.genes.to.TF$microglia[[x]][["cTF"]] <- cTF
  
}


### 4. across all monocytes ###############################

is.mo <- sumF.sc$myelo_celltype %in% c("monocyte") # 409
sce <- sumF.sc[, is.mo]

for (x in as.character(1:100)) { 
  
  module.x <- S.lst$sd3[[x]]$Symbol
  
  ## A. Correlations between genes within a module (intra.module.genes corr):
  # pair-wise stats
  cP <- correlatePairs(sce, subset.row=module.x, pairings=NULL)
  cP <- cP %>% as_tibble() %>% dplyr::filter(FDR <= 0.01)
  # global stats per gene 
  cG <- correlateGenes(cP)
  cG <- as_tibble(cG)
  
  # save
  correlation_pairs$intra.module.genes$monocyte[[x]][["cP"]] <- cP
  correlation_pairs$intra.module.genes$monocyte[[x]][["cG"]] <- cG
  
  
  ## B. Correlations between module genes and TFs (module.genes.to.TF corr):
  # pair-wise stats
  cP.TF <- correlatePairs(sce, subset.row=c(module.x, tf.sumF.hv), pairings=list(module.x, tf.sumF.hv))
  cP.TF <- cP.TF %>% as_tibble() %>% dplyr::filter(FDR <= 0.01)
  # global stats per gene 
  cTF <- correlateGenes(cP.TF)
  cTF <- as_tibble(cTF) %>% dplyr::filter(gene %in% tf.sumF.hv) # filters only TFs and their global scores  
  
  # save
  correlation_pairs$module.genes.to.TF$monocyte[[x]][["cP.TF"]] <- cP.TF
  correlation_pairs$module.genes.to.TF$monocyte[[x]][["cTF"]] <- cTF
  
}


### 5. across all granulocytic ############################

is.gr <- sumF.sc$myelo_celltype %in% c("granulocytopoietic cell", "granulocyte") # 519
sce <- sumF.sc[, is.gr]

for (x in as.character(1:100)) { 
  
  module.x <- S.lst$sd3[[x]]$Symbol
  
  ## A. Correlations between genes within a module (intra.module.genes corr):
  # pair-wise stats
  cP <- correlatePairs(sce, subset.row=module.x, pairings=NULL)
  cP <- cP %>% as_tibble() %>% dplyr::filter(FDR <= 0.01)
  # global stats per gene 
  cG <- correlateGenes(cP)
  cG <- as_tibble(cG)
  
  # save
  correlation_pairs$intra.module.genes$granulocyte[[x]][["cP"]] <- cP
  correlation_pairs$intra.module.genes$granulocyte[[x]][["cG"]] <- cG
  
  
  ## B. Correlations between module genes and TFs (module.genes.to.TF corr):
  # pair-wise stats
  cP.TF <- correlatePairs(sce, subset.row=c(module.x, tf.sumF.hv), pairings=list(module.x, tf.sumF.hv))
  cP.TF <- cP.TF %>% as_tibble() %>% dplyr::filter(FDR <= 0.01)
  # global stats per gene 
  cTF <- correlateGenes(cP.TF)
  cTF <- as_tibble(cTF) %>% dplyr::filter(gene %in% tf.sumF.hv) # filters only TFs and their global scores  
  
  # save
  correlation_pairs$module.genes.to.TF$granulocyte[[x]][["cP.TF"]] <- cP.TF
  correlation_pairs$module.genes.to.TF$granulocyte[[x]][["cTF"]] <- cTF
  
}


### 6. across all endothelial cells #######################

is.endo <- sumF.sc$endo_celltype == "endothelial cell" # 1930
sce <- sumF.sc[, is.endo]

for (x in as.character(1:100)) { 
  
  module.x <- S.lst$sd3[[x]]$Symbol
  
  ## A. Correlations between genes within a module (intra.module.genes corr):
  # pair-wise stats
  cP <- correlatePairs(sce, subset.row=module.x, pairings=NULL)
  cP <- cP %>% as_tibble() %>% dplyr::filter(FDR <= 0.01)
  # global stats per gene 
  cG <- correlateGenes(cP)
  cG <- as_tibble(cG)
  
  # save
  correlation_pairs$intra.module.genes$endothel[[x]][["cP"]] <- cP
  correlation_pairs$intra.module.genes$endothel[[x]][["cG"]] <- cG
  
  
  ## B. Correlations between module genes and TFs (module.genes.to.TF corr):
  # pair-wise stats
  cP.TF <- correlatePairs(sce, subset.row=c(module.x, tf.sumF.hv), pairings=list(module.x, tf.sumF.hv))
  cP.TF <- cP.TF %>% as_tibble() %>% dplyr::filter(FDR <= 0.01)
  # global stats per gene 
  cTF <- correlateGenes(cP.TF)
  cTF <- as_tibble(cTF) %>% dplyr::filter(gene %in% tf.sumF.hv) # filters only TFs and their global scores  
  
  # save
  correlation_pairs$module.genes.to.TF$endothel[[x]][["cP.TF"]] <- cP.TF
  correlation_pairs$module.genes.to.TF$endothel[[x]][["cTF"]] <- cTF
  
}

saveRDS(correlation_pairs, "results/correlation_pairs.rds")

# correlation_pairs <- readRDS("results/correlation_pairs.rds")




# you can summarise there by listing top x TFs that appear most often ???

# for all components you can also list the TFs at "lower" level of the S component - 1 or 2 sd cutoff (module-centric)
# alt. for all expressed 176 TFs in hvg - you can summarise appearance and S score across components (TF-centric)



# Some network visualisations... start with this Tutorial:
# https://kateto.net/network-visualization

# Examples - Module 52 and 23

S.lst <- readRDS(file = here("results", "sumF_S_modules.rds"))
correlation_pairs <- readRDS(file = here("results", "correlation_pairs.rds"))

# 52 intramodular:
m.52.nodes <- S.lst$sd3[['52']]$Symbol

m.52.edges.myelo <- correlation_pairs$intra.module.genes$myeloid$`52`$cP %>% dplyr::filter(rho > 0.3)
m.52.edges.mac <- correlation_pairs$intra.module.genes$macrophage$`52`$cP %>% dplyr::filter(rho > 0.3)
m.52.edges.mic <- correlation_pairs$intra.module.genes$microglia$`52`$cP %>% dplyr::filter(rho > 0.3)
m.52.edges.mo <- correlation_pairs$intra.module.genes$monocyte$`52`$cP %>% dplyr::filter(rho > 0.3)
m.52.edges.gr <- correlation_pairs$intra.module.genes$granulocyte$`52`$cP %>% dplyr::filter(rho > 0.3)

layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] # check different layout options

# myeloid cell in general (myelo)
net.52.myelo <- graph_from_data_frame(d=m.52.edges.myelo, vertices=NULL, directed=F) # vertices=m.52.nodes (if to display unused nodes too)
modul <- cluster_optimal(net.52.myelo, weights=E(net.52.myelo)$rho)
modul.louvian <- cluster_louvain(net.52.myelo, weights=E(net.52.myelo)$rho)
modul.info <- cluster_infomap(net.52.myelo, e.weights=E(net.52.myelo)$rho, v.weights=NULL, nb.trials=100, modularity=TRUE)
V(net.52.myelo)$community <- modul$membership
colrs <- adjustcolor(c("yellowgreen", "tomato", "gold"), alpha=.7)
deg.52.myelo <- degree(net.52.myelo, normalized=F)
pdf(file = "figures/module_networks/52.myelo.pdf", width = 12, height = 10)
plot(net.52.myelo, 
     vertex.color=colrs[V(net.52.myelo)$community], vertex.size=(deg.52.myelo/2)+5, 
     vertex.label.family="Helvetica", vertex.label.cex=0.9, vertex.label.color="gray20", 
     edge.color="gray70", edge.width=(E(net.52.myelo)$rho*3)^3, main="Module 52 - myeloid cells", 
     layout=layout_components, margin=rep(-0.05, 4), asp=0) # axes=F, ylim=c(-1, 0.7), xlim=c(-1,0.6)
dev.off()
# layout_components, layout_with_lgl, layout_with_dh, layout_with_kk, layout_in_circle, layout_nicely, layout_with_graphopt, layout_with_fr, 

# macrophages (mac)
net.52.mac <- graph_from_data_frame(d=m.52.edges.mac, vertices=NULL, directed=F) # vertices=m.52.nodes (if to display unused nodes too)
# modul <- cluster_optimal(net.52.mac, weights=E(net.52.mac)$rho)
modul.louvian <- cluster_louvain(net.52.mac, weights=E(net.52.mac)$rho)
modul.info <- cluster_infomap(net.52.mac, e.weights=E(net.52.mac)$rho, v.weights=NULL, nb.trials=100, modularity=TRUE)
V(net.52.mac)$community <- modul.louvian$membership
colrs <- adjustcolor(c("yellowgreen", "tomato", "gold", "dodgerblue"), alpha=.7)
deg.52.mac <- degree(net.52.mac, normalized=F)
pdf(file = "figures/module_networks/52.mac.pdf", width = 12, height = 10)
plot(net.52.mac, 
     vertex.color=colrs[V(net.52.mac)$community], vertex.size=(deg.52.mac/4)+5, 
     vertex.label.family="Helvetica", vertex.label.cex=0.9, vertex.label.color="gray20", 
     edge.color="gray70", edge.width=(E(net.52.mac)$rho*3)^3, main="Module 52 - macrophages", 
     layout=layout_with_lgl(net.52.mac, area=vcount(net.52.mac)^3, repulserad=50*vcount(net.52.mac)*vcount(net.52.mac)^3), margin=rep(-0.1, 4), asp=0)
dev.off()
# there is need to spread nodes more.. but with that many edges the only way might be
# plotting just nodes (sized by degrees) on a sphere layaout 
pdf(file = "figures/module_networks/52.mac.nodes.pdf", width = 12, height = 12)
plot(net.52.mac, 
     vertex.color=colrs[V(net.52.mac)$community], vertex.size=(deg.52.mac/4)+5, 
     vertex.label.family="Helvetica", vertex.label.cex=0.9, vertex.label.color="gray20", 
     edge.color="gray90", edge.width=(E(net.52.mac)$rho*3)^2, main="Module 52 - macrophages", 
     layout=layout_on_sphere, margin=rep(-0.05, 4), asp=0)
dev.off()

# monocytes (mo)
net.52.mo <- graph_from_data_frame(d=m.52.edges.mo, vertices=NULL, directed=F) # vertices=m.52.nodes (if to display unused nodes too)
modul <- cluster_optimal(net.52.mo)
V(net.52.mo)$community <- modul$membership
colrs <- adjustcolor(c("yellowgreen", "tomato", "gold"), alpha=.7)
deg.52.mo <- degree(net.52.mo, normalized=F)
pdf(file = "figures/module_networks/52.mo.pdf", width = 10, height = 10)
plot(net.52.mo, 
     vertex.color=colrs[V(net.52.mo)$community], vertex.size=(deg.52.mo/2)+5, 
     vertex.label.family="Helvetica", vertex.label.cex=0.9, vertex.label.color="gray20", 
     edge.color="gray70", edge.width=(E(net.52.mo)$rho*3)^3, main="Module 52 - monocytes", 
     layout=layout_components, margin=rep(-0.05, 4), asp=0)
dev.off()

# microglia (mic)
net.52.mic <- graph_from_data_frame(d=m.52.edges.mic, vertices=NULL, directed=F) # vertices=m.52.nodes (if to display unused nodes too)
modul <- cluster_optimal(net.52.mic)
V(net.52.mic)$community <- modul$membership
colrs <- adjustcolor(c("yellowgreen"), alpha=.7)
deg.52.mic <- degree(net.52.mic, normalized=F)
pdf(file = "figures/module_networks/52.mic.pdf", width = 10, height = 10)
plot(net.52.mic, 
     vertex.color=colrs[V(net.52.mic)$community], vertex.size=(deg.52.mic/2)+5, 
     vertex.label.family="Helvetica", vertex.label.cex=0.9, vertex.label.color="gray20", 
     edge.color="gray70", edge.width=(E(net.52.mic)$rho*3)^3, main="Module 52 - microglia", 
     layout=layout_components, margin=rep(-0.05, 4), asp=0)
dev.off()

# granulocytes (gr)
net.52.gr <- graph_from_data_frame(d=m.52.edges.gr, vertices=NULL, directed=F) # vertices=m.52.nodes (if to display unused nodes too)
modul <- cluster_optimal(net.52.gr)
V(net.52.gr)$community <- modul$membership
colrs <- adjustcolor(c("yellowgreen", "tomato"), alpha=.7)
deg.52.gr <- degree(net.52.gr, normalized=F)
pdf(file = "figures/module_networks/52.gr.pdf", width = 10, height = 10)
plot(net.52.gr, 
     vertex.color=colrs[V(net.52.gr)$community], vertex.size=(deg.52.gr/2)+5, 
     vertex.label.family="Helvetica", vertex.label.cex=0.9, vertex.label.color="gray20", 
     edge.color="gray70", edge.width=(E(net.52.gr)$rho*3)^3, main="Module 52 - granulocytes", 
     layout=layout_components, margin=rep(-0.05, 4), asp=0)
dev.off()




# 52 genes & TFs:

m.52.edges.myelo <- correlation_pairs$module.genes.to.TF$myeloid$`52`$cP.TF %>% dplyr::filter(rho > 0.3) %>% dplyr::mutate(gene2 = str_c(".", gene2)) %>% dplyr::rename(weight = rho)
m.52.edges.mac <- correlation_pairs$module.genes.to.TF$macrophage$`52`$cP.TF %>% dplyr::filter(rho > 0.3) %>% dplyr::mutate(gene2 = str_c(".", gene2)) %>% dplyr::rename(weight = rho)
m.52.edges.mic <- correlation_pairs$module.genes.to.TF$microglia$`52`$cP.TF %>% dplyr::filter(rho > 0.3) %>% dplyr::mutate(gene2 = str_c(".", gene2)) %>% dplyr::rename(weight = rho)
m.52.edges.mo <- correlation_pairs$module.genes.to.TF$monocyte$`52`$cP.TF %>% dplyr::filter(rho > 0.3) %>% dplyr::mutate(gene2 = str_c(".", gene2)) %>% dplyr::rename(weight = rho)
m.52.edges.gr <- correlation_pairs$module.genes.to.TF$granulocyte$`52`$cP.TF %>% dplyr::filter(rho > 0.3) %>% dplyr::mutate(gene2 = str_c(".", gene2)) %>% dplyr::rename(weight = rho)


# myeloid cell in general (myelo)
net.52.myelo <- graph_from_data_frame(d=m.52.edges.myelo, vertices=NULL, directed=F) # vertices=m.52.nodes (if to display unused nodes too)
V(net.52.myelo)$type <- bipartite_mapping(net.52.myelo)$type
colrs <- adjustcolor(c("gold", "tomato"), alpha=.7)
deg.52.myelo <- degree(net.52.myelo, normalized=F)
modul <- cluster_optimal(net.52.myelo, weights=E(net.52.myelo)$weight)
# modul.louvian <- cluster_louvain(net.52.myelo, weights=E(net.52.myelo)$weight)

pdf(file = "figures/module_networks/52.myelo.TF.pdf", width = 12, height = 10)
plot(net.52.myelo, 
     vertex.color=colrs[V(net.52.myelo)$type+1], vertex.shape=c("circle","square")[V(net.52.myelo)$type+1], vertex.size=(deg.52.myelo/2)+5, 
     vertex.label.family="Helvetica", vertex.label.cex=0.9, vertex.label.color="gray20", 
     edge.color="gray70", edge.width=(E(net.52.myelo)$weight*3)^3, main="Module 52 - myeloid cells", 
     mark.groups=modul, mark.col=adjustcolor(c("lightsteelblue"), alpha=.2), mark.border=NA, mark.shape=1, 
     layout=layout_nicely, margin=rep(-0.05, 4), asp=0) # axes=F, ylim=c(-1, 0.7), xlim=c(-1,0.6)
dev.off()


# macrophages (mac)
net.52.mac <- graph_from_data_frame(d=m.52.edges.mac, vertices=NULL, directed=F) # vertices=m.52.nodes (if to display unused nodes too)
V(net.52.mac)$type <- bipartite_mapping(net.52.mac)$type
colrs <- adjustcolor(c("gold", "tomato"), alpha=.7)
deg.52.mac <- degree(net.52.mac, normalized=F)
# modul <- cluster_optimal(net.52.mac, weights=E(net.52.mac)$weight)
modul.louvian <- cluster_louvain(net.52.mac, weights=E(net.52.mac)$weight)

pdf(file = "figures/module_networks/52.mac.TF.pdf", width = 12, height = 10)
plot(net.52.mac, 
     vertex.color=colrs[V(net.52.mac)$type+1], vertex.shape=c("circle","square")[V(net.52.mac)$type+1], vertex.size=(deg.52.mac/4)+5, 
     vertex.label.family="Helvetica", vertex.label.cex=0.9, vertex.label.color="gray20", 
     edge.color="gray70", edge.width=(E(net.52.mac)$weight*3)^2, main="Module 52 - macrophages", 
     mark.groups=modul.louvian, mark.col=adjustcolor(c("lightsteelblue"), alpha=.2), mark.border=NA, mark.shape=1, 
     layout=layout_with_graphopt(net.52.mac, charge=0.05), margin=rep(-0.05, 4), asp=0) # axes=F, ylim=c(-1, 0.7), xlim=c(-1,0.6)
dev.off()


# monocyte (mo)
net.52.mo <- graph_from_data_frame(d=m.52.edges.mo, vertices=NULL, directed=F) # vertices=m.52.nodes (if to display unused nodes too)
V(net.52.mo)$type <- bipartite_mapping(net.52.mo)$type
colrs <- adjustcolor(c("gold", "tomato"), alpha=.7)
deg.52.mo <- degree(net.52.mo, normalized=F)
modul <- cluster_optimal(net.52.mo, weights=E(net.52.mo)$weight)
# modul.louvian <- cluster_louvain(net.52.mo, weights=E(net.52.mo)$weight)

pdf(file = "figures/module_networks/52.mo.TF.pdf", width = 10, height = 8)
plot(net.52.mo, 
     vertex.color=colrs[V(net.52.mo)$type+1], vertex.shape=c("circle","square")[V(net.52.mo)$type+1], vertex.size=(deg.52.mo/2)+5, 
     vertex.label.family="Helvetica", vertex.label.cex=0.9, vertex.label.color="gray20", 
     edge.color="gray70", edge.width=(E(net.52.mo)$weight*3)^3, main="Module 52 - monocytes", 
     mark.groups=modul, mark.col=adjustcolor(c("lightsteelblue"), alpha=.2), mark.border=NA, mark.shape=1, 
     layout=layout_as_tree, margin=rep(-0.05, 4), asp=0) # axes=F, ylim=c(-1, 0.7), xlim=c(-1,0.6)
dev.off()


# microglia (mic)
net.52.mic <- graph_from_data_frame(d=m.52.edges.mic, vertices=NULL, directed=F) # vertices=m.52.nodes (if to display unused nodes too)
V(net.52.mic)$type <- bipartite_mapping(net.52.mic)$type
colrs <- adjustcolor(c("gold", "tomato"), alpha=.7)
deg.52.mic <- degree(net.52.mic, normalized=F)
modul <- cluster_optimal(net.52.mic, weights=E(net.52.mic)$weight)
# modul.louvian <- cluster_louvain(net.52.mic, weights=E(net.52.mic)$weight)

pdf(file = "figures/module_networks/52.mic.TF.pdf", width = 8, height = 8)
plot(net.52.mic, 
     vertex.color=colrs[V(net.52.mic)$type+1], vertex.shape=c("circle","square")[V(net.52.mic)$type+1], vertex.size=(deg.52.mic/2)+5, 
     vertex.label.family="Helvetica", vertex.label.cex=0.9, vertex.label.color="gray20", 
     edge.color="gray70", edge.width=(E(net.52.mic)$weight*3)^3, main="Module 52 - microglia", 
     mark.groups=modul, mark.col=adjustcolor(c("lightsteelblue"), alpha=.2), mark.border=NA, mark.shape=1, 
     layout=layout_nicely, margin=rep(-0.05, 4), asp=0) # axes=F, ylim=c(-1, 0.7), xlim=c(-1,0.6)
dev.off()


# granulocytes (gr)
net.52.gr <- graph_from_data_frame(d=m.52.edges.gr, vertices=NULL, directed=F) # vertices=m.52.nodes (if to display unused nodes too)
V(net.52.gr)$type <- bipartite_mapping(net.52.gr)$type
colrs <- adjustcolor(c("gold", "tomato"), alpha=.7)
deg.52.gr <- degree(net.52.gr, normalized=F)
modul <- cluster_optimal(net.52.gr, weights=E(net.52.gr)$weight)
# modul.louvian <- cluster_louvain(net.52.gr, weights=E(net.52.gr)$weight)

pdf(file = "figures/module_networks/52.gr.TF.pdf", width = 12, height = 10)
plot(net.52.gr, 
     vertex.color=colrs[V(net.52.gr)$type+1], vertex.shape=c("circle","square")[V(net.52.gr)$type+1], vertex.size=(deg.52.gr/2)+5, 
     vertex.label.family="Helvetica", vertex.label.cex=0.9, vertex.label.color="gray20", 
     edge.color="gray70", edge.width=(E(net.52.gr)$weight*3)^3, main="Module 52 - granulocytes", 
     mark.groups=modul, mark.col=adjustcolor(c("lightsteelblue"), alpha=.2), mark.border=NA, mark.shape=1, 
     layout=layout_nicely, margin=rep(-0.05, 4), asp=0) # axes=F, ylim=c(-1, 0.7), xlim=c(-1,0.6)
dev.off()





# NODES	 
# vertex.color	 Node color
# vertex.frame.color	 Node border color
# vertex.shape	 One of “none”, “circle”, “square”, “csquare”, “rectangle”
# “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
# vertex.size	 Size of the node (default is 15)
# vertex.size2	 The second size of the node (e.g. for a rectangle)
# vertex.label	 Character vector used to label the nodes
# vertex.label.family	 Font family of the label (e.g.“Times”, “Helvetica”)
# vertex.label.font	 Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
# vertex.label.cex	 Font size (multiplication factor, device-dependent)
# vertex.label.dist	 Distance between the label and the vertex
# vertex.label.degree	 The position of the label in relation to the vertex, where
# 0 is right, “pi” is left, “pi/2” is below, and “-pi/2” is above
# EDGES	 
# edge.color	 Edge color
# edge.width	 Edge width, defaults to 1
# edge.arrow.size	 Arrow size, defaults to 1
# edge.arrow.width	 Arrow width, defaults to 1
# edge.lty	 Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”,
# 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
# edge.label	 Character vector used to label edges
# edge.label.family	 Font family of the label (e.g.“Times”, “Helvetica”)
# edge.label.font	 Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
# edge.label.cex	 Font size for edge labels
# edge.curved	 Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
# arrow.mode	 Vector specifying whether edges should have arrows,
# possible values: 0 no arrow, 1 back, 2 forward, 3 both
# OTHER	 
# margin	 Empty space margins around the plot, vector with length 4
# frame	 if TRUE, the plot will be framed
# main	 If set, adds a title to the plot
# sub	 If set, adds a subtitle to the plot
# asp	  Numeric, the aspect ratio of a plot (y/x).
# palette	 A color palette to use for vertex color
# rescale	 Whether to rescale coordinates to [-1,1]. Default is TRUE.





# 23
m.23.nodes <- S.lst$sd3[['23']]$Symbol










S.lst$sd2$`52`$Symbol
S.lst$sd3$`52`$Symbol
S.lst$sd4$`52`$Symbol
S.lst$all$`52` %>% dplyr::filter(Symbol == "Spi1")
S.lst$all$`52` %>% dplyr::filter(Symbol == "Spib")

sumF.S %>% as.data.frame() %>% rownames_to_column(var="Symbol") %>% dplyr::filter(Symbol == "Spi1") %>% slice(1)
# so you cam expore by tf across cponents !


S.lst$sd3$`52` %>% dplyr::filter(S >0) %>% dplyr::select(Symbol) %>% pull()
S.lst$sd3$`52` %>% dplyr::filter(S <0) %>% dplyr::select(Symbol) %>% pull()
S.lst$mad5$`52` %>% dplyr::filter(S >0) %>% dplyr::select(Symbol) %>% pull()
S.lst$mad5$`52` %>% dplyr::filter(S <0) %>% dplyr::select(Symbol) %>% pull()

S.lst$sd3$`63` %>% dplyr::filter(S >0) %>% dplyr::select(Symbol) %>% pull()
S.lst$sd3$`63` %>% dplyr::filter(S <0) %>% dplyr::select(Symbol) %>% pull()
S.lst$mad5$`63` %>% dplyr::filter(S >0) %>% dplyr::select(Symbol) %>% pull()
S.lst$mad5$`63` %>% dplyr::filter(S <0) %>% dplyr::select(Symbol) %>% pull()

S.lst$sd3$`35` %>% dplyr::filter(S >0) %>% dplyr::select(Symbol) %>% pull()
S.lst$sd3$`35` %>% dplyr::filter(S <0) %>% dplyr::select(Symbol) %>% pull()

S.lst$sd3$`15` %>% dplyr::filter(S >0) %>% dplyr::select(Symbol) %>% pull()
S.lst$sd3$`15` %>% dplyr::filter(S <0) %>% dplyr::select(Symbol) %>% pull()

S.lst$sd3$`40` %>% dplyr::filter(S >0) %>% dplyr::select(Symbol) %>% pull()
S.lst$sd3$`40` %>% dplyr::filter(S <0) %>% dplyr::select(Symbol) %>% pull()

S.lst$sd3$`70` %>% dplyr::filter(S >0) %>% dplyr::select(Symbol) %>% pull()
S.lst$sd3$`70` %>% dplyr::filter(S <0) %>% dplyr::select(Symbol) %>% pull()

S.lst$sd3$`71` %>% dplyr::filter(S >0) %>% dplyr::select(Symbol) %>% pull()
S.lst$sd3$`71` %>% dplyr::filter(S <0) %>% dplyr::select(Symbol) %>% pull()

S.lst$sd3$`72` %>% dplyr::filter(S >0) %>% dplyr::select(Symbol) %>% pull()
S.lst$sd3$`72` %>% dplyr::filter(S <0) %>% dplyr::select(Symbol) %>% pull()

S.lst$sd3$`79` %>% dplyr::filter(S >0) %>% dplyr::select(Symbol) %>% pull()
S.lst$sd3$`79` %>% dplyr::filter(S <0) %>% dplyr::select(Symbol) %>% pull()

S.lst$sd3$`6` %>% dplyr::filter(S >0) %>% dplyr::select(Symbol) %>% pull()
S.lst$sd3$`6` %>% dplyr::filter(S <0) %>% dplyr::select(Symbol) %>% pull()
















