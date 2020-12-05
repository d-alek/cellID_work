### modules_by_cells.R

library("tidyverse")
library("dplyr")
library("magrittr")
library("stringr")
library("ggplot2")
library("cowplot")
library("directlabels")
library("RColorBrewer")
library("SingleCellExperiment") # SCE
library("scran")
library("scater")
library("here")
library("caret")
library("parallel")
library("doParallel")
library("glmnet")
library("ComplexHeatmap")
library("circlize")



# import starting SCE objects for spikeF and sumF versions (need cells metadata)
spikeF.sc <- readRDS(file = here("data", "processed", "tm.smartseq2.spikeF.123.rds"))
sumF.sc <- readRDS(file = here("data", "processed", "tm.smartseq2.sumF.123.rds"))

# define myeloid and endothelial types 
myelo.names <- c("basophil", "brain pericyte", "classical monocyte", "granulocyte", "granulocyte monocyte progenitor cell", 
                 "granulocytopoietic cell", "Kupffer cell", "macrophage", "microglial cell", "monocyte")
endo.names <- c("endothelial cell", "endothelial cell of hepatic sinusoid", "lung endothelial cell") # "endothelial cell" in 10 tissues
# NOTE: I included "brain pericytes" here as "myeloid".. published, although they are not the classic myeloid cell.. interesting to explore

# colect useful cell metadata and reorganise cell annotations.. for myeloid- and endothelial-focused analysis
cellData <- colData(spikeF.sc) %>% as.data.frame() %>% as_tibble() %>% 
  dplyr::select(cell, cell_ontology_class, tissue, mouse_id, mouse_sex, plate_barcode, total_counts, total_features_by_counts) %>% 
  dplyr::mutate(tissue = str_replace_all(tissue, "_", " ")) %>% 
  dplyr::rename(celltype = cell_ontology_class) %>% 
  tidyr::unite(celltype, tissue, col = celltype_tissue, sep = " - ", remove = F) %>% 
  dplyr::mutate(myelo_celltype = if_else(celltype %in% myelo.names, celltype, "non-myeloid")) %>% 
  dplyr::mutate(myelo_celltype = if_else(myelo_celltype == "classical monocyte", "monocyte", myelo_celltype)) %>% 
  tidyr::unite(myelo_celltype, tissue, col = myelo_celltype_tissue, sep = " - ", remove = F) %>% 
  dplyr::mutate(myelo_celltype_tissue = if_else(str_detect(myelo_celltype_tissue, "^non-myeloid -"), "non-myeloid", myelo_celltype_tissue)) %>% 
  dplyr::mutate(endo_celltype = if_else(celltype %in% endo.names, "endothelial cell", "non-endothelial")) %>% 
  tidyr::unite(endo_celltype, tissue, col = endo_celltype_tissue, sep = " - ", remove = F) %>% 
  dplyr::mutate(endo_celltype_tissue = if_else(str_detect(endo_celltype_tissue, "^non-endothelial -"), "non-endothelial", endo_celltype_tissue)) %>% 
  dplyr::select(cell, celltype, tissue, celltype_tissue, myelo_celltype, myelo_celltype_tissue, endo_celltype, endo_celltype_tissue, 
                mouse_id, mouse_sex, plate_barcode, total_counts, total_features_by_counts)

# some checks
cellData %>% dplyr::count(celltype) %>% print(n=100) # 1457
cellData %>% dplyr::count(tissue) %>% print(n=100)
cellData %>% dplyr::count(celltype_tissue) %>% print(n=101)
cellData %>% dplyr::count(myelo_celltype_tissue) %>% print(n=100)
cellData %>% dplyr::count(endo_celltype_tissue) %>% print(n=100)

# import consensus ncomp=100 A matrices for spikeF and sumF versions (extracted in icafast_extractor.r/.sh)
# spikeF:
spikeF.A <- readRDS(file = here("results", "spikeF_A_c.rds"))
spikeF.A <- spikeF.A[["100"]]
rownames(spikeF.A) <- spikeF.sc$cell
# sumF:
sumF.A <- readRDS(file = here("results", "sumF_A_c.rds"))
sumF.A <- sumF.A[["100"]]
rownames(sumF.A) <- sumF.sc$cell

# make it into A.100.df and attach cell annotations and metadata
# spikeF:
spikeF.A.df <- spikeF.A %>% 
  as.data.frame() %>% rownames_to_column(var = "cell") %>% as_tibble() %>% 
  left_join(y = cellData, by = "cell")
# sumF:
sumF.A.df <- sumF.A %>% 
  as.data.frame() %>% rownames_to_column(var = "cell") %>% as_tibble() %>% 
  left_join(y = cellData, by = "cell")

# save both
saveRDS(spikeF.A.df, "results/spikeF_A_c_meta.rds") # 18420 x 113
saveRDS(sumF.A.df, "results/sumF_A_c_meta.rds") # 18420 x 113

# spikeF.A.df <- readRDS(file = here("results", "spikeF_A_c_meta.rds"))
# sumF.A.df <- readRDS(file = here("results", "sumF_A_c_meta.rds"))


#---------------------------------------------------------------------------------------#
# As a quick detour here, it is useful to attach the consensus ICA reduced dimensions to the SCE objects above, 
# for easy access later on via reducedDims slot. 
# Both A comps (sampleFactors i.e. cellFactors) and S comps (featureLoadings i.e. geneLoadings) can be attached:

# import consensus ncomp=100 S matrices for spikeF and sumF versions (extracted in icafast_extractor.r/.sh)
# spikeF:
spikeF.S <- readRDS(file = here("results", "spikeF_S_c.rds"))
spikeF.S <- spikeF.S[["100"]]
# sumF:
sumF.S <- readRDS(file = here("results", "sumF_S_c.rds"))
sumF.S <- sumF.S[["100"]]  

# create LEM objects to attach to SCE-s
# spikeF:
spikeF.ica <- LinearEmbeddingMatrix(sampleFactors = spikeF.A, featureLoadings = spikeF.S)
rownames(spikeF.ica) <- rownames(spikeF.A)
colnames(spikeF.ica) <- colnames(spikeF.A)
# sumF:
sumF.ica <- LinearEmbeddingMatrix(sampleFactors = sumF.A, featureLoadings = sumF.S)
rownames(sumF.ica) <- rownames(sumF.A)
colnames(sumF.ica) <- colnames(sumF.A)

# trim SCE objects to appropriate hvg/HVG
# spikeF:
keep.hvg.nonspike <- rowData(spikeF.sc)$is_hvg
spikeF.sc <- spikeF.sc[keep.hvg.nonspike, ]
# sumF:
keep.HVG.nonspike <- rowData(sumF.sc)$is_HVG
sumF.sc <- sumF.sc[keep.HVG.nonspike, ]

# add reducedDim to SCE
# spikeF:
reducedDim(spikeF.sc, type="ICA") <- spikeF.ica
# sumF:
reducedDim(sumF.sc, type="ICA") <- sumF.ica

# save these trimmed and reducedDim updated SCE!
saveRDS(spikeF.sc, "data/processed/tm.smartseq2.spikeF.hvg.123.icaDims.rds") # 6342 x 18420
saveRDS(sumF.sc, "data/processed/tm.smartseq2.sumF.hvg.123.icaDims.rds") # 4748 x 18420

# spikeF.sc <- readRDS(file = here("data", "processed", "tm.smartseq2.spikeF.hvg.123.icaDims.rds"))
# sumF.sc <- readRDS(file = here("data", "processed", "tm.smartseq2.sumF.hvg.123.icaDims.rds"))
#---------------------------------------------------------------------------------------#


# Do UMAP based on both attached ICA cellfactors and original logcounts.. then visualise few ICA components 

#--- UMAP using existing ICA reduced dims ------------------------------------#
set.seed(7012347)
sumF.sc <- runUMAP(sumF.sc, dimred="ICA", name="UMAP.ica", n_neighbors = 30, min_dist=0.75, spread=1.1)
# quick testing parameter adequacy
# detour: to use them for visualisation, it is easier to append ICA cellFactors to colData (for now)
ica_factor_df <- as.data.frame(sampleFactors(reducedDim(sumF.sc, "ICA")))
all.equal(rownames(colData(sumF.sc)), rownames(ica_factor_df)) # TRUE
colData(sumF.sc) <- cbind(colData(sumF.sc), ica_factor_df) 
## celltype:
nb.cols <- length(unique(sumF.sc$cell_ontology_class))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
plotReducedDim(sumF.sc, dimred="UMAP.ica", colour_by="cell_ontology_class", text_by="cell_ontology_class", 
               text_size=3, text_colour="gray20", point_alpha=0.95, point_size=2) + 
  scale_fill_manual(values = mycolors) + 
  theme(legend.position = "none")
## components:
plotReducedDim(sumF.sc, dimred="UMAP.ica", colour_by="52", text_by="cell_ontology_class", 
               text_size=3, text_colour="gray20", point_alpha=0.7, point_size=2) + 
  scale_fill_gradient2("52", low="gray90", high="red3")
## some genes:
plotReducedDim(sumF.sc, dimred="UMAP.ica", colour_by="Spi1", text_by="cell_ontology_class", 
               text_size=3, text_colour="gray20", point_alpha=0.7, point_size=2) + 
  scale_fill_gradient("Spi1", low="gray95", high="springgreen2")


#--- UMAP using logcounts and newly calc. PCA reduced dims -------------------#
sumF.sc <- runPCA(sumF.sc, dimred=NULL, exprs_values="logcounts", subset_row=1:4748, scale=TRUE, ncomponents=100, name="PCA")
set.seed(98701)
sumF.sc <- runUMAP(sumF.sc, dimred="PCA", name="UMAP.pca", n_neighbors = 40, min_dist=0.65, spread=1.1)
# quick testing parameter adequacy
## celltype:
nb.cols <- length(unique(sumF.sc$cell_ontology_class))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
plotReducedDim(sumF.sc, dimred="UMAP.pca", colour_by="cell_ontology_class", text_by="cell_ontology_class", 
               text_size=3, text_colour="gray20", point_alpha=0.95, point_size=2) + 
  scale_fill_manual(values = mycolors) + 
  theme(legend.position = "none")
## components:
plotReducedDim(sumF.sc, dimred="UMAP.pca", colour_by="52", text_by="cell_ontology_class", 
               text_size=3, text_colour="gray20", point_alpha=0.7, point_size=2) + 
  scale_fill_gradient2("52", low="gray90", high="red3")
## some genes:
plotReducedDim(sumF.sc, dimred="UMAP.pca", colour_by="Spi1", text_by="cell_ontology_class", 
               text_size=3, text_colour="gray20", point_alpha=0.7, point_size=2) + 
  scale_fill_gradient("Spi1", low="gray95", high="springgreen2")


# Prepare a big tibble for ggplot2 plotting cells flexibly in UMAP coordinates - containing: 
# cell annotations, UMAP coords, ICA coords, gene expression for all ~ 4750 genes 
cellDataUmapPlot <- colData(sumF.sc) %>% as.data.frame() %>% as_tibble() %>% 
  dplyr::select(cell, cell_ontology_class, tissue) %>% 
  dplyr::mutate(celltype=case_when(cell_ontology_class=="Brush cell of epithelium proper of large intestine"~"Brush cell of large intestine", 
                                   cell_ontology_class=="ciliated columnar cell of tracheobronchial tree"~"ciliated cell of tracheobronchial tree", 
                                   cell_ontology_class=="common lymphoid progenitor"~"CLP", 
                                   cell_ontology_class=="DN1 thymic pro-T cell"~"DN1 pro-T cell", 
                                   cell_ontology_class=="enterocyte of epithelium of large intestine"~"enterocyte of large intestine", 
                                   cell_ontology_class=="granulocyte monocyte progenitor cell"~"GMP", 
                                   cell_ontology_class=="hematopoietic precursor cell"~"hematopoietic precursor", 
                                   cell_ontology_class=="immature natural killer cell"~"immature NK cell", 
                                   cell_ontology_class=="mature natural killer cell"~"mature NK cell", 
                                   cell_ontology_class=="megakaryocyte-erythroid progenitor cell"~"MEP", 
                                   cell_ontology_class=="myofibroblast cell"~"myofibroblast", 
                                   cell_ontology_class=="natural killer cell"~"NK cell", 
                                   cell_ontology_class=="oligodendrocyte precursor cell"~"oligodendrocyte precursor", 
                                   cell_ontology_class=="pre-natural killer cell"~"pre-NK cell", 
                                   cell_ontology_class=="Slamf1-negative multipotent progenitor cell"~"Slamf1- progenitor", 
                                   cell_ontology_class=="Slamf1-positive multipotent progenitor cell"~"Slamf1+ progenitor", 
                                   cell_ontology_class=="type B pancreatic cell"~"pancreatic B cell", 
                                   TRUE ~ cell_ontology_class)) %>% 
  dplyr::mutate(tissue = str_replace_all(tissue, "_", " ")) %>% 
  tidyr::unite(celltype, tissue, col = celltype_tissue, sep = " - ", remove = F) %>% 
  dplyr::select(cell, celltype, tissue, celltype_tissue)

# only when UMAP parameters are decided, then get coordinates:
umap_ica <- reducedDim(sumF.sc, "UMAP.ica") %>% as.data.frame() %>% rownames_to_column() %>% as_tibble() %>% 
  rename(rowname="cell", V1="UMAP1.ica", V2="UMAP2.ica")
umap_pcs <- reducedDim(sumF.sc, "UMAP.pca") %>% as.data.frame() %>% rownames_to_column() %>% as_tibble() %>% 
  rename(rowname="cell", V1="UMAP1.pca", V2="UMAP2.pca")

# get ICA coordinates:
sumF.A.df <- readRDS(file = here("results", "sumF_A_c_meta.rds"))
ica_A <- sumF.A.df %>% select(cell, `1`:`100`)

# get gene expressions:
head(colnames(logcounts(sumF.sc)))
head(rownames(logcounts(sumF.sc)))
genexp <- t(as.matrix(logcounts(sumF.sc))) %>% as.data.frame() %>% rownames_to_column("cell") %>% as_tibble()

# join all in a single tibble
cellDataUmapPlot <- cellDataUmapPlot %>% 
  left_join(y=umap_ica, by="cell") %>% 
  left_join(y=umap_pcs, by="cell") %>% 
  left_join(y=ica_A, by="cell") %>% 
  left_join(y=genexp, by="cell")
write_tsv(cellDataUmapPlot, "results/cellDataUmapPlot.txt")
# cellDataUmapPlot <- read_tsv("results/cellDataUmapPlot.txt")


#--- UMAP using existing ICA reduced dims ------------------------------------#

# Then finally build your way of gglotting cells in UMAP coordinates...
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
get.medians <- gapply.fun(dl.summarize(d, x = median(x), y = median(y)))

# Cell types plot
ggplot(data=cellDataUmapPlot, aes(x=UMAP1.ica, y=UMAP2.ica, fill=celltype)) + 
  geom_point(shape=21, size=2.5, colour="grey75", stroke=0.25, alpha=0.7) + 
  directlabels::geom_dl(aes(label=celltype), method=list("get.medians", cex=0.5, alpha=0.6)) + 
  theme_bw() + 
  scale_fill_manual(values = mycolors) + 
  theme(legend.position = "none", panel.grid=element_blank())
ggsave("figures/umaps/allcells_umap_ica.pdf", width = 8, height = 8)
ggsave("figures/umaps/allcells_umap_ica.png", width = 8, height = 8)

# Organs plots
nb.tissues <- length(unique(sumF.sc$tissue))
mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.tissues)
ggplot(data=cellDataUmapPlot, aes(x=UMAP1.ica, y=UMAP2.ica, fill=tissue)) + 
  geom_point(shape=21, size=2.5, colour="grey75", stroke=0.25, alpha=0.9) + 
  theme_bw() + 
  scale_fill_manual(values = mycolors) + 
  theme(panel.grid=element_blank())
ggsave("figures/umaps/alltissues_umap_ica.pdf", width = 8.75, height = 8)
ggsave("figures/umaps/alltissues_umap_ica.png", width = 8.75, height = 8)

# Together celltype & organs
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
p1 <- ggplot(data=cellDataUmapPlot, aes(x=UMAP1.ica, y=UMAP2.ica, fill=celltype)) + 
  geom_point(shape=21, size=2.5, colour="grey75", stroke=0.2, alpha=0.7) + 
  directlabels::geom_dl(aes(label=celltype), method=list("get.medians", cex=0.5, alpha=0.6)) + 
  theme_bw() + 
  scale_fill_manual(values = mycolors) + 
  theme(legend.position = "none", panel.grid=element_blank())

nb.tissues <- length(unique(sumF.sc$tissue))
mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.tissues)
p2 <- ggplot(data=cellDataUmapPlot, aes(x=UMAP1.ica, y=UMAP2.ica, fill=tissue)) + 
  geom_point(shape=21, size=2.5, colour="grey75", stroke=0.2, alpha=0.9) + 
  theme_bw() + 
  scale_fill_manual(values = mycolors) + 
  theme(panel.grid=element_blank())

plot_grid(p1, p2, rel_widths = c(1, 1.25))
ggsave("figures/umaps/allcellstissues_umap_ica.pdf", width = 13, height = 6)
ggsave("figures/umaps/allcellstissues_umap_ica.png", width = 13, height = 6)


# Components plots
# note directlabels does not work the best here, when legend is drawn! Can use annotate() for manual labelling of few cell types
ggplot(data=cellDataUmapPlot, aes(x=UMAP1.ica, y=UMAP2.ica, fill=`63`)) + 
  geom_point(shape=21, size=2.5, colour="grey75", stroke=0.25, alpha=0.7) + 
  theme_bw() + 
  scale_fill_gradient2(low="dodgerblue1", high="red3", space="Lab") + 
  theme(panel.grid=element_blank())
ggsave("figures/umaps/comp63_umap_ica.pdf", width = 8.25, height = 7.75)
ggsave("figures/umaps/comp63_umap_ica.png", width = 8.25, height = 7.75)

# Genes plots
ggplot(data=cellDataUmapPlot, aes(x=UMAP1.ica, y=UMAP2.ica, fill=Spi1)) + 
  geom_point(shape=21, size=2.5, colour="grey75", stroke=0.25, alpha=0.8) + 
  theme_bw() + 
  scale_fill_gradient2(low="grey95", high="springgreen2", space="Lab") + 
  theme(panel.grid=element_blank())
ggsave("figures/umaps/gSpi1_umap_ica.pdf", width = 8.25, height = 7.75)
ggsave("figures/umaps/gSpi1_umap_ica.png", width = 8.25, height = 7.75)

# MAKE THEM ALL!
for (i in 1:100) {
  i <- as.character(i)
  ggplot(data=cellDataUmapPlot, aes(x=UMAP1.ica, y=UMAP2.ica, fill=get(i))) + 
    geom_point(shape=21, size=2.5, colour="grey75", stroke=0.2, alpha=0.7) + 
    theme_bw() + 
    scale_fill_gradient2(name=i, low="dodgerblue1", high="red3", space="Lab") + 
    theme(panel.grid=element_blank())
  ggsave(paste0("figures/umaps/comp" ,i, "_umap_ica.pdf"), width = 8.25, height = 7.5)
  ggsave(paste0("figures/umaps/comp" ,i, "_umap_ica.png"), width = 8.25, height = 7.5)
}



#--- UMAP using PCA reduced dims ---------------------------------------------#

# Cell types plot
ggplot(data=cellDataUmapPlot, aes(x=UMAP1.pca, y=UMAP2.pca, fill=celltype)) + 
  geom_point(shape=21, size=2.5, colour="grey75", stroke=0.2, alpha=0.7) + 
  directlabels::geom_dl(aes(label=celltype), method=list("get.medians", cex=0.5, alpha=0.6)) + 
  theme_bw() + 
  scale_fill_manual(values = mycolors) + 
  theme(legend.position = "none", panel.grid=element_blank())
ggsave("figures/umaps/allcells_umap_pca.pdf", width = 8, height = 8)
ggsave("figures/umaps/allcells_umap_pca.png", width = 8, height = 8)

# Organs plots
nb.tissues <- length(unique(sumF.sc$tissue))
mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.tissues)
ggplot(data=cellDataUmapPlot, aes(x=UMAP1.pca, y=UMAP2.pca, fill=tissue)) + 
  geom_point(shape=21, size=2.5, colour="grey75", stroke=0.2, alpha=0.9) + 
  theme_bw() + 
  scale_fill_manual(values = mycolors) + 
  theme(panel.grid=element_blank())
ggsave("figures/umaps/alltissues_umap_pca.pdf", width = 8.75, height = 8)
ggsave("figures/umaps/alltissues_umap_pca.png", width = 8.75, height = 8)

# Together celltype & organs
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
p1 <- ggplot(data=cellDataUmapPlot, aes(x=UMAP1.pca, y=UMAP2.pca, fill=celltype)) + 
  geom_point(shape=21, size=2.5, colour="grey75", stroke=0.2, alpha=0.7) + 
  directlabels::geom_dl(aes(label=celltype), method=list("get.medians", cex=0.5, alpha=0.6)) + 
  theme_bw() + 
  scale_fill_manual(values = mycolors) + 
  theme(legend.position = "none", panel.grid=element_blank())

nb.tissues <- length(unique(sumF.sc$tissue))
mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.tissues)
p2 <- ggplot(data=cellDataUmapPlot, aes(x=UMAP1.pca, y=UMAP2.pca, fill=tissue)) + 
  geom_point(shape=21, size=2.5, colour="grey75", stroke=0.2, alpha=0.9) + 
  theme_bw() + 
  scale_fill_manual(values = mycolors) + 
  theme(panel.grid=element_blank())

plot_grid(p1, p2, rel_widths = c(1, 1.25))
ggsave("figures/umaps/umaps_pca/allcellstissues_umap_pca.pdf", width = 13, height = 6)
ggsave("figures/umaps/umaps_pca/allcellstissues_umap_pca.png", width = 13, height = 6)

# Components plots
# note directlabels does not work the best here, when legend is drawn! Can use annotate() for manual labelling of few cell types
ggplot(data=cellDataUmapPlot, aes(x=UMAP1.pca, y=UMAP2.pca, fill=`63`)) + 
  geom_point(shape=21, size=2.5, colour="grey75", stroke=0.25, alpha=0.7) + 
  theme_bw() + 
  scale_fill_gradient2(low="dodgerblue1", high="red3", space="Lab") + 
  theme(panel.grid=element_blank())
ggsave("figures/umaps/umaps_pca/comp63_umap_pca.pdf", width = 8.25, height = 7.75)
ggsave("figures/umaps/umaps_pca/comp63_umap_pca.png", width = 8.25, height = 7.75)

# Genes plots
ggplot(data=cellDataUmapPlot, aes(x=UMAP1.pca, y=UMAP2.pca, fill=Spi1)) + 
  geom_point(shape=21, size=2.5, colour="grey75", stroke=0.25, alpha=0.8) + 
  theme_bw() + 
  scale_fill_gradient2(low="grey95", high="springgreen2", space="Lab") + 
  theme(panel.grid=element_blank())
ggsave("figures/umaps/umaps_pca/gSpi1_umap_pca.pdf", width = 8.25, height = 7.75)
ggsave("figures/umaps/umaps_pca/gSpi1_umap_pca.png", width = 8.25, height = 7.75)

# MAKE THEM ALL!
for (i in 1:100) {
  i <- as.character(i)
  ggplot(data=cellDataUmapPlot, aes(x=UMAP1.pca, y=UMAP2.pca, fill=get(i))) + 
    geom_point(shape=21, size=2.5, colour="grey75", stroke=0.2, alpha=0.7) + 
    theme_bw() + 
    scale_fill_gradient2(name=i, low="dodgerblue1", high="red3", space="Lab") + 
    theme(panel.grid=element_blank())
  ggsave(paste0("figures/umaps/umaps_pca/comp" ,i, "_umap_pca.pdf"), width = 8.25, height = 7.5)
  ggsave(paste0("figures/umaps/umaps_pca/comp" ,i, "_umap_pca.png"), width = 8.25, height = 7.5)
}













#---------------------------------------------------------------------------------------#
# Before going into biology, it would be good to see if any of the 100 components are significantly associated with 
# technical experimental factors such as mouse sex, mouse id, total counts per cell, and number of detected genes per cell 
# (rather than biology) - and earmark them or maybe even remove them

#-----Below is exploration of this issue on the gene level first (rather than enitre component):

# This computes % of variance of gene expression that is explained by some of these technical factors 
# (as well as celltype as biological factor) for each gene. It is based on R-squared form wrapped lm(gene.exp. ~ factor).
# It can be used to spot a problematic technical experimental factors, and genes that are most affected by them.

# saved exp_fact_var_spikeF.pdf (% variance explained by various experimental factors)
spikeF.genes.r2mat <- getVarianceExplained(spikeF.sc, exprs_values = "logcounts", 
                                           variables = c("total_counts", "total_features_by_counts", 
                                                         "plate_barcode", "mouse_id", "mouse_sex", "cell_ontology_class"))
p1 <- plotExplanatoryVariables(spikeF.genes.r2mat, variables = c("total_counts", "total_features_by_counts", 
                                                                 "plate_barcode", "mouse_id", "mouse_sex", "cell_ontology_class"), 
                               theme_size = 10)
# saved exp_fact_var_sumF.pdf (% variance explained by various experimental factors)
sumF.genes.r2mat <- getVarianceExplained(sumF.sc, exprs_values = "logcounts", 
                                         variables = c("total_counts", "total_features_by_counts", 
                                                       "plate_barcode", "mouse_id", "mouse_sex", "cell_ontology_class"))
p2 <- plotExplanatoryVariables(sumF.genes.r2mat, variables = c("total_counts", "total_features_by_counts", 
                                                               "plate_barcode", "mouse_id", "mouse_sex", "cell_ontology_class"), 
                               theme_size = 10)
pdf(file = "figures/varByFact_genes.pdf", width = 10, height = 4)
multiplot(p1 + ggtitle("Gen.exp variance - spikeF normalisation only") + scale_y_continuous(limits=c(0,2), breaks=seq(0,2, by=0.5)), 
          p2 + ggtitle("Gen.exp variance - sumF normalisation") + scale_y_continuous(limits=c(0,2), breaks=seq(0,2, by=0.5)), cols = 2)
dev.off()
# as expected, variance explained by total counts and total detected genes per cell is much smaller in sumF normalised data..
# plate bacrcode is expected to explain lots of variance as the experiment design was such (complete overlap between batch
# and biology - plates and tissues/cell types)

# just a quick manual exploration here..

# spikeF:
spikeF.genes.r2 <- spikeF.genes.r2mat %>% as.data.frame() %>% rownames_to_column(var = "GeneSymbol") %>% as_tibble()

spikeF.genes.r2 %>% dplyr::arrange(desc(total_counts)) # Ph4b max 0.1
spikeF.genes.r2 %>% dplyr::arrange(desc(total_features_by_counts)) # Rn45s max 0.5
spikeF.genes.r2 %>% dplyr::arrange(desc(mouse_id)) # max 0.15
spikeF.genes.r2 %>% dplyr::arrange(desc(mouse_sex)) # Sod1 max 0.05
spikeF.genes.r2 %>% dplyr::arrange(desc(cell_ontology_class)) # this pointless unless stratified on individual celltypes

# sumF:
sumF.genes.r2 <- sumF.genes.r2mat %>% as.data.frame() %>% rownames_to_column(var = "GeneSymbol") %>% as_tibble()

sumF.genes.r2 %>% dplyr::arrange(desc(total_counts)) # Ph4b max 0.06 !
sumF.genes.r2 %>% dplyr::arrange(desc(total_features_by_counts)) %>% print(n=15) # Ppa1 max 0.27
sumF.genes.r2 %>% dplyr::arrange(desc(mouse_id)) # max 0.18
sumF.genes.r2 %>% dplyr::arrange(desc(mouse_sex)) # Skil max 0.05
sumF.genes.r2 %>% dplyr::arrange(desc(cell_ontology_class)) # this pointless unless stratified on individual celltypes


#-----The same can be done on the component level by creating SCE with components as (reduced dimension) variables instead of genes:

# This computes % of variance of A component (cellFactor) that is explained by some of these technical factors 
# (as well as celltype as biological factor) for each component. The components explained largely by technical factors can be
# either excluded from downstream component-celltype analysis or taken with reservation if they are selected in the LASSO models below

# spikeF:
spikeF.mock.sc.A <- SingleCellExperiment(assays = list(ica.cellFactors = t(spikeF.A)), colData = cellData) # or attach all collData
spikeF.ica.r2mat <- getVarianceExplained(spikeF.mock.sc.A, exprs_values = "ica.cellFactors", 
                                         variables = c("total_counts", "total_features_by_counts", 
                                                       "plate_barcode", "mouse_id", "mouse_sex", "celltype", 
                                                       "myelo_celltype", "myelo_celltype_tissue", "endo_celltype", "endo_celltype_tissue"))
# plotExplanatoryVariables(spikeF.ica.r2mat)
p3 <- plotExplanatoryVariables(spikeF.ica.r2mat[,c(1:6)])

# sumF:
sumF.mock.sc.A <- SingleCellExperiment(assays = list(ica.cellFactors = t(sumF.A)), colData = cellData) # or attach all collData
sumF.ica.r2mat <- getVarianceExplained(sumF.mock.sc.A, exprs_values = "ica.cellFactors", 
                                         variables = c("total_counts", "total_features_by_counts", 
                                                       "plate_barcode", "mouse_id", "mouse_sex", "celltype", 
                                                       "myelo_celltype", "myelo_celltype_tissue", "endo_celltype", "endo_celltype_tissue"))
# plotExplanatoryVariables(sumF.ica.r2mat)
p4 <- plotExplanatoryVariables(sumF.ica.r2mat[,c(1:6)])
pdf(file = "figures/varByFact_components.pdf", width = 10, height = 4)
multiplot(p3 + ggtitle("ICA component variance - spikeF normalisation only") + scale_y_continuous(limits=c(0,3), breaks=seq(0,3, by=0.5)), 
          p4 + ggtitle("ICA component variance - sumF normalisation") + scale_y_continuous(limits=c(0,3), breaks=seq(0,3, by=0.5)), cols = 2)
dev.off()

# It is actually more useful to identify individual components that are associated with key technical covariates
# just a quick manual exploration here..

# spikeF:
spikeF.ica.r2 <- spikeF.ica.r2mat %>% as.data.frame() %>% rownames_to_column(var = "component") %>% as_tibble()
# technical
spikeF.ica.r2 %>% dplyr::arrange(desc(total_counts)) # 17,16 (both < 0.15)
spikeF.ica.r2 %>% dplyr::arrange(desc(total_features_by_counts)) %>% print(n=16) # 14,2,21,4,17,12,29,16,59!,41,10,18 (all > 0.4)
spikeF.ica.r2 %>% dplyr::arrange(desc(mouse_id)) # 34,44 (both < 0.4)
spikeF.ica.r2 %>% dplyr::arrange(desc(mouse_sex)) # 44 (< 0.25)
# biological (below, LASSO will better pick a combination of these per celltype)
spikeF.ica.r2 %>% dplyr::arrange(desc(myelo_celltype)) # 38,26,32,60,87,57,40,51,59! (above 0.2)
spikeF.ica.r2 %>% dplyr::arrange(desc(myelo_celltype_tissue)) # 38,26,60,32,57,87,40,59!,51
spikeF.ica.r2 %>% dplyr::arrange(desc(endo_celltype)) # 9
spikeF.ica.r2 %>% dplyr::arrange(desc(endo_celltype_tissue)) # 73,9,78,88,98

# sumF:
sumF.ica.r2 <- sumF.ica.r2mat %>% as.data.frame() %>% rownames_to_column(var = "component") %>% as_tibble()
# technical
sumF.ica.r2 %>% dplyr::arrange(desc(total_counts)) # 18 (0.07)
sumF.ica.r2 %>% dplyr::arrange(desc(total_features_by_counts)) %>% print(n=15) # 20,21,2,6,8,14,18,70,17,7,71!,47, (all > 0.2, < 0.4!)
sumF.ica.r2 %>% dplyr::arrange(desc(mouse_id)) # 33,73
sumF.ica.r2 %>% dplyr::arrange(desc(mouse_sex)) # 73 (< 0.3)
# biological (below, LASSO will better pick a combination of these per celltype)
sumF.ica.r2 %>% dplyr::arrange(desc(myelo_celltype)) # 23,15,52,40,71!,81,72,36,70,51 (above 0.2)
sumF.ica.r2 %>% dplyr::arrange(desc(myelo_celltype_tissue)) # 15,23,40,52,81,71!,36,72,51,70
sumF.ica.r2 %>% dplyr::arrange(desc(endo_celltype)) # 9, 69
sumF.ica.r2 %>% dplyr::arrange(desc(endo_celltype_tissue)) # 67,9,82,89,69,55

# as seen above, the components somewhat associated with technical factors don't seem to be be associated with celltype factors 
# with exception of comp spikeF 59 / sumF 71 (they are associated with total_features/cell.. which could be related with cell identity as well).
# I will run the first round of analysis below without removing any components (but this could be re-investigated)
# THE EASIEST WAY TO EXCLUDE COMPONENTS FROM DOWNSTREAM ANALYSIS IS VIA glmnet() / cv.glmnet() argument "exclude" !
# IT IS ALSO POSSIBLE TO INCREASE PENALTY ON SUSPICIOUS COMPONENT VIA glmnet() / cv.glmnet() argument "penalty.factor"

# IT WOULD BE GOOD TO VISUALLY SUMMARISE SOME OF TOP COMPONENT - EXP.FACTOR ASSOCIATIONS ABOVE (DOT / LOLIPOP PLOT)

#---------------------------------------------------------------------------------------#


# ... EXPLAIN ...

#------Split data to Training and Testing - taking care of celltype:tissue balance------#
library(caret)

set.seed(123)
inTraining <- createDataPartition(y = sumF.A.df$celltype_tissue, p = .85, list = FALSE) 
#---------------------------------------------------------------------------------------#


# LASSO models to run to find a combinations of components (modules) that best explain varoius cell identites 
# (each model below on two versions of data normalisation - spikeF and sumF)
# (each model below in 4 versions - lasso + deviance, lasso + misclass.error, el.net + deviance, el.net + misclass.error):


#-- Myeloid cells
### A. Components
# 1. myeloid cell types + others
# 2. myeloid cell types only ? (compare with conclusions from 1.)
# 3. macrohages vs. ALL other cells ? (compare with conclusions from 1. and 2.) (binomial)
# 4. macrophages across tissues only
# 5. monocytes across two tissues 
# 6. classic monocytes vs. monocytes (lung)

### B. Genes
# 1. myeloid cell types + others
# 2. myeloid cell types only ? (compare with conclusions from 1.)
# 3. macrohages vs. ALL other cells ? (compare with conclusions from 1. and 2.) (binomial)
# 4. macrophages across tissues
# 5. monocytes across two tissues 
# 6. classic monocytes vs. monocytes (lung)

# (Myelo.A.1.spikeF.lasso.dev).. sumF.elnet.merr..

#-- Endothelial cells
### A. Components
# 1. endothelial cells + others (binomial)
# 2. endothelial cells across tissues

### B. Genes
# 1. endothelial cells + others (binomial)
# 2. endothelial cells across tissues


# THAT IS 128 MODELS TO FIT, CROSSVALIDATE ETC...!
# PERHAPS DO THIS FOR SUMMARY STATS TO DECIDE WHICH ONE TO GO WITH 
# THEN CAN BE SIMPLIFIED TO ~ 7 MODELS BASICALLY


#-- Setup data inputs (A and Gen.exp matrices) - variables in columns for glmnet!
data <- list("A_spikeF" = spikeF.A.df, 
             "A_sumF" = sumF.A.df, 
             "GE_spikeF" = t(as.matrix(logcounts(spikeF.sc))), 
             "GE_sumF" = t(as.matrix(logcounts(sumF.sc))))


###### I run the 8 analyses below (~1000 lines of code that follows) on 8 cores in parallel #######
###################################################################################################
#-- Setup a list of CV and Prediction results (cv.glmnet objects, their betas and predicts, confussion matrix stats)

out.struct <- list("cv.fit"=list(), 
                    "betas"=list("nBetasOveral"=0, "nBetasByClass"=numeric(), "Betas"=tibble()), 
                    "predicts"=list("statsOveral"=numeric(), "statsByClass"=list("balancedAcc"=numeric(), "F1"=numeric()), 
                                    "predictClass"=character(), "predictProbs"=array(), "confussionMat"=list()))

four.methods <- list("lasso.deviance" = out.struct, 
                     "lasso.classerr" = out.struct, 
                     "elnet.deviance" = out.struct, 
                     "elnet.classerr" = out.struct)

list.struct <- list("A_spikeF" = four.methods, 
                    "A_sumF" = four.methods, 
                    "GE_spikeF" = four.methods, 
                    "GE_sumF" = four.methods) 

glmnetCV <- list("Myelo.1" = list.struct, 
                 "Myelo.2" = list.struct, 
                 "Myelo.3" = list.struct, 
                 "Myelo.4" = list.struct, 
                 "Myelo.5" = list.struct, 
                 "Myelo.6" = list.struct, 
                 "Endo.1" = list.struct, 
                 "Endo.2" = list.struct)


#-- Setup parallel 
cores <- makeCluster(detectCores(logical = FALSE), type='PSOCK')
registerDoParallel(cores)


#-----------------------------------------------------------------------------#
### Myeloid 1. (myeloid cell types + others)
#-----------------------------------------------------------------------------#

for (i in 1:4) {
  
  if (i <= 2) {
    A.df <- data[[i]] # best not to introduce factor levels here for easier tibble wrangling
    
    training <- A.df[inTraining, ]
    testing <- A.df[-inTraining, ]
    
    if (i == 1) {
      print("----- Myeloid 1 training -----")
      print(dplyr::count(training, myelo_celltype)); cat("\n")
      print("----- Myeloid 1 testing -----")
      print(dplyr::count(testing, myelo_celltype)); cat("\n")
      
      #-- Setup outcomes once - they do not change per set (e.g. "Myelo.1)
      # training:
      y <- training$myelo_celltype %>% 
        fct_relevel(c("granulocyte monocyte progenitor cell", "granulocytopoietic cell", "granulocyte", "basophil", 
                      "monocyte", "macrophage", "Kupffer cell", "microglial cell", "brain pericyte", 
                      "non-myeloid")) %>% 
        fct_recode(GMP = "granulocyte monocyte progenitor cell")
      # testing:
      truth <- testing$myelo_celltype %>% 
        fct_relevel(c("granulocyte monocyte progenitor cell", "granulocytopoietic cell", "granulocyte", "basophil", 
                      "monocyte", "macrophage", "Kupffer cell", "microglial cell", "brain pericyte", 
                      "non-myeloid")) %>% 
        fct_recode(GMP = "granulocyte monocyte progenitor cell")
    }
    
    x <- as.matrix(training[, as.character(1:100)])
    newx <- as.matrix(testing[, as.character(1:100)])
    
  } else if (i > 2) {
    GE.mat <- data[[i]]
    
    training <- GE.mat[inTraining, ]
    testing <- GE.mat[-inTraining, ]
    
    x <- training
    newx <- testing

  }

  # Fit & CV (lasso deviance)--------------------------------------------------#
  fit.cv.1 <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="deviance", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.1[[i]]$lasso.deviance$cv.fit <- fit.cv.1
  #--start plot
  pdf(file = paste0("results/fits/Myelo.1.CVplot.", names(data)[i], ".pdf"), width = 7.5, height = 6, pointsize = 10)
  par(mfrow=c(2,2))
  plot(fit.cv.1); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.1", "lasso", "deviance", "-", names(data)[i]))
  
  # Fit & CV (lasso class error)----------------------------------------------#
  fit.cv.2 <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="class", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.1[[i]]$lasso.classerr$cv.fit <- fit.cv.2
  plot(fit.cv.2); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.1", "lasso", "classerr", "-", names(data)[i]))
  
  # Fit & CV (elnet deviance)-------------------------------------------------#
  fit.cv.3 <- cv.glmnet(x=x, y=y, alpha=0.9, family="multinomial", type.measure="deviance", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.1[[i]]$elnet.deviance$cv.fit <- fit.cv.3
  plot(fit.cv.3); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.1", "elnet", "deviance", "-", names(data)[i]))
  
  # Fit & CV (elnet class error)----------------------------------------------#
  fit.cv.4 <- cv.glmnet(x=x, y=y, alpha=0.9, family="multinomial", type.measure="class", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.1[[i]]$elnet.classerr$cv.fit <- fit.cv.4
  plot(fit.cv.4); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.1", "elnet", "classerr", "-", names(data)[i]))
  
  dev.off()
  
  # Then extract coefficiants, predictions and confussion stats from the CV objects
  for (j in 1:4) {
    
    fit.cv <- get(paste0("fit.cv.", j))
    
    #-----------Betas-------------------------------------#
    betas <- coef(fit.cv, s="lambda.1se")
    outcomes <- names(betas)
    betas <- purrr::map(betas, .f = as.matrix) %>% do.call(cbind, .)
    colnames(betas) <- outcomes
    
    betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% 
      as_tibble() %>% dplyr::filter(compID != "(Intercept)") %>% dplyr::filter_at(vars(-compID), any_vars(.!=0)) 
    nbetas.overall <- nrow(betas)
    nbetas.byclass <- betas %>% summarize_at(vars(-compID), .funs = function(x) sum(x != 0)) %>% unlist()
    
    glmnetCV$Myelo.1[[i]][[j]]$betas$nBetasOveral <- nbetas.overall
    glmnetCV$Myelo.1[[i]][[j]]$betas$nBetasByClass <- nbetas.byclass
    glmnetCV$Myelo.1[[i]][[j]]$betas$Betas <- betas
    
    #-----------Predictions & Confussion stats------------#
    prediction <- predict(fit.cv, newx=newx, s="lambda.1se", type="class") # a vector
    prediction <- factor(prediction, levels=levels(truth))
    
    prediction.p <- predict(fit.cv, newx=newx, s="lambda.1se", type="response")
    
    confussion <- confusionMatrix(data=prediction, reference=truth, positive = NULL, dnn = c("Prediction", "Truth"), mode = "everything") # positive - change it for two factors
    
    confussion.overal <- confussion$overall[1:2] # just Accuracy & Kappa
    confussion.byclass <- confussion$byClass[, c("Balanced Accuracy", "F1")] # matrix!
    rownames(confussion.byclass) <- str_sub(rownames(confussion.byclass), start=8)
    
    glmnetCV$Myelo.1[[i]][[j]]$predicts$statsOveral <- confussion.overal
    glmnetCV$Myelo.1[[i]][[j]]$predicts$statsByClass$balancedAcc <- confussion.byclass[,1]
    glmnetCV$Myelo.1[[i]][[j]]$predicts$statsByClass$F1 <- confussion.byclass[,2]
    glmnetCV$Myelo.1[[i]][[j]]$predicts$predictClass <- prediction
    glmnetCV$Myelo.1[[i]][[j]]$predicts$predictProbs <- prediction.p
    glmnetCV$Myelo.1[[i]][[j]]$predicts$confussionMat <- confussion
    
  }
  
  print(paste("Myelo.1 -", names(data)[i], "... Done:", Sys.time()))
}


#-----------------------------------------------------------------------------#
### Myeloid 2. (myeloid cell types only)
#-----------------------------------------------------------------------------#

for (i in 1:4) {
  
  if (i <= 2) {
    A.df <- data[[i]] # best not to introduce factor levels here for easier tibble wrangling
    
    training <- A.df[inTraining, ] %>% dplyr::filter(myelo_celltype != "non-myeloid")
    testing <- A.df[-inTraining, ] %>% dplyr::filter(myelo_celltype != "non-myeloid")
    
    if (i == 1) {
      print("----- Myeloid 2 training -----")
      print(dplyr::count(training, myelo_celltype)); cat("\n")
      print("----- Myeloid 2 testing -----")
      print(dplyr::count(testing, myelo_celltype)); cat("\n")
      
      #-- Setup outcomes - they do not change per set (e.g. "Myelo.1)
      # training:
      y <- training$myelo_celltype %>% 
        fct_relevel(c("granulocyte monocyte progenitor cell", "granulocytopoietic cell", "granulocyte", "basophil", 
                      "monocyte", "macrophage", "Kupffer cell", "microglial cell", "brain pericyte")) %>% 
        fct_recode(GMP = "granulocyte monocyte progenitor cell")
      # testing:
      truth <- testing$myelo_celltype %>% 
        fct_relevel(c("granulocyte monocyte progenitor cell", "granulocytopoietic cell", "granulocyte", "basophil", 
                      "monocyte", "macrophage", "Kupffer cell", "microglial cell", "brain pericyte")) %>% 
        fct_recode(GMP = "granulocyte monocyte progenitor cell")
      
      #-- Setup myeloid subset cellnames for GE.mat subsetting
      inTraining.myeloid <- training$cell
      inTesting.myeloid <- testing$cell
    }
    
    x <- as.matrix(training[, as.character(1:100)])
    newx <- as.matrix(testing[, as.character(1:100)])
    
  } else if (i > 2) {
    GE.mat <- data[[i]]
    
    training <- GE.mat[inTraining.myeloid, ] # select only myeloid as above
    testing <- GE.mat[inTesting.myeloid, ]
    
    x <- training
    newx <- testing
    
  }
  
  # Fit & CV (lasso deviance)--------------------------------------------------#
  fit.cv.1 <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="deviance", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.2[[i]]$lasso.deviance$cv.fit <- fit.cv.1
  #--start plot
  pdf(file = paste0("results/fits/Myelo.2.CVplot.", names(data)[i], ".pdf"), width = 7.5, height = 6, pointsize = 10)
  par(mfrow=c(2,2))
  plot(fit.cv.1); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.2", "lasso", "deviance", "-", names(data)[i]))
  
  # Fit & CV (lasso class error)----------------------------------------------#
  fit.cv.2 <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="class", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.2[[i]]$lasso.classerr$cv.fit <- fit.cv.2
  plot(fit.cv.2); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.2", "lasso", "classerr", "-", names(data)[i]))
  
  # Fit & CV (elnet deviance)-------------------------------------------------#
  fit.cv.3 <- cv.glmnet(x=x, y=y, alpha=0.9, family="multinomial", type.measure="deviance", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.2[[i]]$elnet.deviance$cv.fit <- fit.cv.3
  plot(fit.cv.3); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.2", "elnet", "deviance", "-", names(data)[i]))
  
  # Fit & CV (elnet class error)----------------------------------------------#
  fit.cv.4 <- cv.glmnet(x=x, y=y, alpha=0.9, family="multinomial", type.measure="class", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.2[[i]]$elnet.classerr$cv.fit <- fit.cv.4
  plot(fit.cv.4); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.2", "elnet", "classerr", "-", names(data)[i]))
  
  dev.off()
  
  # Then extract coefficiants, predictions and confussion stats from the CV objects
  for (j in 1:4) {
    
    fit.cv <- get(paste0("fit.cv.", j))
    
    #-----------Betas-------------------------------------#
    betas <- coef(fit.cv, s="lambda.1se")
    outcomes <- names(betas)
    betas <- purrr::map(betas, .f = as.matrix) %>% do.call(cbind, .)
    colnames(betas) <- outcomes
    
    betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% 
      as_tibble() %>% dplyr::filter(compID != "(Intercept)") %>% dplyr::filter_at(vars(-compID), any_vars(.!=0)) 
    nbetas.overall <- nrow(betas)
    nbetas.byclass <- betas %>% summarize_at(vars(-compID), .funs = function(x) sum(x != 0)) %>% unlist()
    
    glmnetCV$Myelo.2[[i]][[j]]$betas$nBetasOveral <- nbetas.overall
    glmnetCV$Myelo.2[[i]][[j]]$betas$nBetasByClass <- nbetas.byclass
    glmnetCV$Myelo.2[[i]][[j]]$betas$Betas <- betas
    
    #-----------Predictions & Confussion stats------------#
    prediction <- predict(fit.cv, newx=newx, s="lambda.1se", type="class") # a vector
    prediction <- factor(prediction, levels=levels(truth))
    
    prediction.p <- predict(fit.cv, newx=newx, s="lambda.1se", type="response")
    
    confussion <- confusionMatrix(data=prediction, reference=truth, positive = NULL, dnn = c("Prediction", "Truth"), mode = "everything") # positive - change it for two factors
    
    confussion.overal <- confussion$overall[1:2] # just Accuracy & Kappa
    confussion.byclass <- confussion$byClass[, c("Balanced Accuracy", "F1")] # matrix!
    rownames(confussion.byclass) <- str_sub(rownames(confussion.byclass), start=8)
    
    glmnetCV$Myelo.2[[i]][[j]]$predicts$statsOveral <- confussion.overal
    glmnetCV$Myelo.2[[i]][[j]]$predicts$statsByClass$balancedAcc <- confussion.byclass[,1]
    glmnetCV$Myelo.2[[i]][[j]]$predicts$statsByClass$F1 <- confussion.byclass[,2]
    glmnetCV$Myelo.2[[i]][[j]]$predicts$predictClass <- prediction
    glmnetCV$Myelo.2[[i]][[j]]$predicts$predictProbs <- prediction.p
    glmnetCV$Myelo.2[[i]][[j]]$predicts$confussionMat <- confussion
    
  }
  
  print(paste("Myelo.2 -", names(data)[i], "... Done:", Sys.time()))
}
  

#-----------------------------------------------------------------------------#
### Myeloid 3. (macrohages vs. ALL other, inc. myeloid non-macs) (binomial)
#-----------------------------------------------------------------------------#

for (i in 1:4) {
  
  if (i <= 2) {
    A.df <- data[[i]] # best not to introduce factor levels here for easier tibble wrangling
    
    training <- A.df[inTraining, ] %>% dplyr::mutate(myelo_celltype = if_else(myelo_celltype == "macrophage", "macrophage", "other"))
    testing <- A.df[-inTraining, ] %>% dplyr::mutate(myelo_celltype = if_else(myelo_celltype == "macrophage", "macrophage", "other"))
    
    if (i == 1) {
      print("----- Myeloid 3 training -----")
      print(dplyr::count(training, myelo_celltype)); cat("\n")
      print("----- Myeloid 3 testing -----")
      print(dplyr::count(testing, myelo_celltype)); cat("\n")
      
      #-- Setup outcomes - they do not change per set (e.g. "Myelo.1)
      # training:
      y <- training$myelo_celltype %>% fct_relevel(c("other", "macrophage"))
      # testing:
      truth <- testing$myelo_celltype %>% fct_relevel(c("other", "macrophage"))
    }
    
    x <- as.matrix(training[, as.character(1:100)])
    newx <- as.matrix(testing[, as.character(1:100)])
    
  } else if (i > 2) {
    GE.mat <- data[[i]]
    
    training <- GE.mat[inTraining, ]
    testing <- GE.mat[-inTraining, ]
    
    x <- training
    newx <- testing
    
  }
  
  # Fit & CV (lasso deviance)--------------------------------------------------#
  fit.cv.1 <- cv.glmnet(x=x, y=y, alpha=1, family="binomial", type.measure="deviance", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.3[[i]]$lasso.deviance$cv.fit <- fit.cv.1
  #--start plot
  pdf(file = paste0("results/fits/Myelo.3.CVplot.", names(data)[i], ".pdf"), width = 7.5, height = 6, pointsize = 10)
  par(mfrow=c(2,2))
  plot(fit.cv.1); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.3", "lasso", "deviance", "-", names(data)[i]))
  
  # Fit & CV (lasso class error)----------------------------------------------#
  fit.cv.2 <- cv.glmnet(x=x, y=y, alpha=1, family="binomial", type.measure="class", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.3[[i]]$lasso.classerr$cv.fit <- fit.cv.2
  plot(fit.cv.2); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.3", "lasso", "classerr", "-", names(data)[i]))
  
  # Fit & CV (elnet deviance)-------------------------------------------------#
  fit.cv.3 <- cv.glmnet(x=x, y=y, alpha=0.9, family="binomial", type.measure="deviance", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.3[[i]]$elnet.deviance$cv.fit <- fit.cv.3
  plot(fit.cv.3); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.3", "elnet", "deviance", "-", names(data)[i]))
  
  # Fit & CV (elnet class error)----------------------------------------------#
  fit.cv.4 <- cv.glmnet(x=x, y=y, alpha=0.9, family="binomial", type.measure="class", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.3[[i]]$elnet.classerr$cv.fit <- fit.cv.4
  plot(fit.cv.4); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.3", "elnet", "classerr", "-", names(data)[i]))
  
  dev.off()
  
  # Then extract coefficiants, predictions and confussion stats from the CV objects
  for (j in 1:4) {
    
    fit.cv <- get(paste0("fit.cv.", j))
    
    #-----------Betas-------------------------------------#
    betas <- coef(fit.cv, s = "lambda.1se")
    betas <- as.matrix(betas)
    colnames(betas) <- levels(y)[2]
    
    betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% 
      as_tibble() %>% dplyr::filter(compID != "(Intercept)") %>% dplyr::filter_at(vars(-compID), any_vars(.!=0)) 
    nbetas.overall <- nrow(betas)
    nbetas.byclass <- betas %>% summarize_at(vars(-compID), .funs = function(x) sum(x != 0)) %>% unlist()
    
    glmnetCV$Myelo.3[[i]][[j]]$betas$nBetasOveral <- nbetas.overall
    glmnetCV$Myelo.3[[i]][[j]]$betas$nBetasByClass <- nbetas.byclass
    glmnetCV$Myelo.3[[i]][[j]]$betas$Betas <- betas
    
    #-----------Predictions & Confussion stats------------#
    prediction <- predict(fit.cv, newx=newx, s="lambda.1se", type="class") # a vector
    prediction <- factor(prediction, levels=levels(truth))
    
    prediction.p <- predict(fit.cv, newx=newx, s="lambda.1se", type="response")
    
    confussion <- confusionMatrix(data=prediction, reference=truth, positive=levels(truth)[2], dnn=c("Prediction", "Truth"), mode="everything") # positive - change it for two factors
    
    confussion.overal <- confussion$overall[1:2] # just Accuracy & Kappa
    confussion.byclass <- matrix(confussion$byClass[c("Balanced Accuracy", "F1")], nrow=1) # make it matrix for consistency!
    rownames(confussion.byclass) <- levels(truth)[2]
    
    glmnetCV$Myelo.3[[i]][[j]]$predicts$statsOveral <- confussion.overal
    glmnetCV$Myelo.3[[i]][[j]]$predicts$statsByClass$balancedAcc <- confussion.byclass[,1]
    glmnetCV$Myelo.3[[i]][[j]]$predicts$statsByClass$F1 <- confussion.byclass[,2]
    glmnetCV$Myelo.3[[i]][[j]]$predicts$predictClass <- prediction
    glmnetCV$Myelo.3[[i]][[j]]$predicts$predictProbs <- prediction.p
    glmnetCV$Myelo.3[[i]][[j]]$predicts$confussionMat <- confussion
    
  }
  
  print(paste("Myelo.3 -", names(data)[i], "... Done:", Sys.time()))
}


#-----------------------------------------------------------------------------#
### Myeloid 4. (macrophages across tissues only, inc. Kupfer & microglia too)
#-----------------------------------------------------------------------------#

for (i in 1:4) {
  
  if (i <= 2) {
    A.df <- data[[i]] # best not to introduce factor levels here for easier tibble wrangling
    
    macs.tissues <- c("macrophage - Marrow", "macrophage - Spleen", "macrophage - Limb Muscle", "macrophage - Kidney", "Kupffer cell - Liver", "macrophage - Brain Myeloid", "microglial cell - Brain Myeloid")
    
    training <- A.df[inTraining, ] %>% dplyr::filter(myelo_celltype_tissue %in% macs.tissues)
    testing <- A.df[-inTraining, ] %>% dplyr::filter(myelo_celltype_tissue %in% macs.tissues)
    
    if (i == 1) {
      print("----- Myeloid 4 training -----")
      print(dplyr::count(training, myelo_celltype_tissue)); cat("\n")
      print("----- Myeloid 4 testing -----")
      print(dplyr::count(testing, myelo_celltype_tissue)); cat("\n")
      
      #-- Setup outcomes - they do not change per set (e.g. "Myelo.1)
      # training:
      y <- training$myelo_celltype_tissue %>% 
        fct_relevel(macs.tissues) %>% 
        fct_recode(`macrophage - Brain` = "macrophage - Brain Myeloid", `microglial cell - Brain` = "microglial cell - Brain Myeloid")
      # testing:
      truth <- testing$myelo_celltype_tissue %>% 
        fct_relevel(macs.tissues) %>% 
        fct_recode(`macrophage - Brain` = "macrophage - Brain Myeloid", `microglial cell - Brain` = "microglial cell - Brain Myeloid")
      
      #-- Setup myeloid subset cellnames for GE.mat subsetting
      inTraining.macs.tissues <- training$cell
      inTesting.macs.tissues <- testing$cell
    }
    
    x <- as.matrix(training[, as.character(1:100)])
    newx <- as.matrix(testing[, as.character(1:100)])
    
  } else if (i > 2) {
    GE.mat <- data[[i]]
    
    training <- GE.mat[inTraining.macs.tissues, ] # select only macrophages as above
    testing <- GE.mat[inTesting.macs.tissues, ]
    
    x <- training
    newx <- testing
    
  }
  
  # Fit & CV (lasso deviance)--------------------------------------------------#
  fit.cv.1 <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="deviance", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.4[[i]]$lasso.deviance$cv.fit <- fit.cv.1
  #--start plot
  pdf(file = paste0("results/fits/Myelo.4.CVplot.", names(data)[i], ".pdf"), width = 7.5, height = 6, pointsize = 10)
  par(mfrow=c(2,2))
  plot(fit.cv.1); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.4", "lasso", "deviance", "-", names(data)[i]))
  
  # Fit & CV (lasso class error)----------------------------------------------#
  fit.cv.2 <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="class", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.4[[i]]$lasso.classerr$cv.fit <- fit.cv.2
  plot(fit.cv.2); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.4", "lasso", "classerr", "-", names(data)[i]))
  
  # Fit & CV (elnet deviance)-------------------------------------------------#
  fit.cv.3 <- cv.glmnet(x=x, y=y, alpha=0.9, family="multinomial", type.measure="deviance", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.4[[i]]$elnet.deviance$cv.fit <- fit.cv.3
  plot(fit.cv.3); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.4", "elnet", "deviance", "-", names(data)[i]))
  
  # Fit & CV (elnet class error)----------------------------------------------#
  fit.cv.4 <- cv.glmnet(x=x, y=y, alpha=0.9, family="multinomial", type.measure="class", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.4[[i]]$elnet.classerr$cv.fit <- fit.cv.4
  plot(fit.cv.4); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.4", "elnet", "classerr", "-", names(data)[i]))
  
  dev.off()
  
  # Then extract coefficiants, predictions and confussion stats from the CV objects
  for (j in 1:4) {
    
    fit.cv <- get(paste0("fit.cv.", j))
    
    #-----------Betas-------------------------------------#
    betas <- coef(fit.cv, s="lambda.1se")
    outcomes <- names(betas)
    betas <- purrr::map(betas, .f = as.matrix) %>% do.call(cbind, .)
    colnames(betas) <- outcomes
    
    betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% 
      as_tibble() %>% dplyr::filter(compID != "(Intercept)") %>% dplyr::filter_at(vars(-compID), any_vars(.!=0)) 
    nbetas.overall <- nrow(betas)
    nbetas.byclass <- betas %>% summarize_at(vars(-compID), .funs = function(x) sum(x != 0)) %>% unlist()
    
    glmnetCV$Myelo.4[[i]][[j]]$betas$nBetasOveral <- nbetas.overall
    glmnetCV$Myelo.4[[i]][[j]]$betas$nBetasByClass <- nbetas.byclass
    glmnetCV$Myelo.4[[i]][[j]]$betas$Betas <- betas
    
    #-----------Predictions & Confussion stats------------#
    prediction <- predict(fit.cv, newx=newx, s="lambda.1se", type="class") # a vector
    prediction <- factor(prediction, levels=levels(truth))
    
    prediction.p <- predict(fit.cv, newx=newx, s="lambda.1se", type="response")
    
    confussion <- confusionMatrix(data=prediction, reference=truth, positive = NULL, dnn = c("Prediction", "Truth"), mode = "everything") # positive - change it for two factors
    
    confussion.overal <- confussion$overall[1:2] # just Accuracy & Kappa
    confussion.byclass <- confussion$byClass[, c("Balanced Accuracy", "F1")] # matrix!
    rownames(confussion.byclass) <- str_sub(rownames(confussion.byclass), start=8)
    
    glmnetCV$Myelo.4[[i]][[j]]$predicts$statsOveral <- confussion.overal
    glmnetCV$Myelo.4[[i]][[j]]$predicts$statsByClass$balancedAcc <- confussion.byclass[,1]
    glmnetCV$Myelo.4[[i]][[j]]$predicts$statsByClass$F1 <- confussion.byclass[,2]
    glmnetCV$Myelo.4[[i]][[j]]$predicts$predictClass <- prediction
    glmnetCV$Myelo.4[[i]][[j]]$predicts$predictProbs <- prediction.p
    glmnetCV$Myelo.4[[i]][[j]]$predicts$confussionMat <- confussion
    
  }
  
  print(paste("Myelo.4 -", names(data)[i], "... Done:", Sys.time()))
}


#-----------------------------------------------------------------------------#
### Myeloid 5. (monocytes across two tissues only) (binomial)
#-----------------------------------------------------------------------------#

for (i in 1:4) {
  
  if (i <= 2) {
    A.df <- data[[i]] # best not to introduce factor levels here for easier tibble wrangling
    
    training <- A.df[inTraining, ] %>% dplyr::filter(myelo_celltype_tissue %in% c("monocyte - Marrow", "monocyte - Lung"))
    testing <- A.df[-inTraining, ] %>% dplyr::filter(myelo_celltype_tissue %in% c("monocyte - Marrow", "monocyte - Lung"))
    
    if (i == 1) {
      print("----- Myeloid 5 training -----")
      print(dplyr::count(training, myelo_celltype_tissue)); cat("\n")
      print("----- Myeloid 5 testing -----")
      print(dplyr::count(testing, myelo_celltype_tissue)); cat("\n")
      
      #-- Setup outcomes - they do not change per set (e.g. "Myelo.1)
      # training:
      y <- training$myelo_celltype_tissue %>% 
        fct_relevel(c("monocyte - Marrow", "monocyte - Lung"))
      # testing:
      truth <- testing$myelo_celltype_tissue %>% 
        fct_relevel(c("monocyte - Marrow", "monocyte - Lung"))
      
      #-- Setup myeloid subset cellnames for GE.mat subsetting
      inTraining.mo.tissues <- training$cell
      inTesting.mo.tissues <- testing$cell
    }
    
    x <- as.matrix(training[, as.character(1:100)])
    newx <- as.matrix(testing[, as.character(1:100)])
    
  } else if (i > 2) {
    GE.mat <- data[[i]]
    
    training <- GE.mat[inTraining.mo.tissues, ] # select only monocytes as above
    testing <- GE.mat[inTesting.mo.tissues, ]
    
    x <- training
    newx <- testing
    
  }
  
  # Fit & CV (lasso deviance)--------------------------------------------------#
  fit.cv.1 <- cv.glmnet(x=x, y=y, alpha=1, family="binomial", type.measure="deviance", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.5[[i]]$lasso.deviance$cv.fit <- fit.cv.1
  #--start plot
  pdf(file = paste0("results/fits/Myelo.5.CVplot.", names(data)[i], ".pdf"), width = 7.5, height = 6, pointsize = 10)
  par(mfrow=c(2,2))
  plot(fit.cv.1); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.5", "lasso", "deviance", "-", names(data)[i]))
  
  # Fit & CV (lasso class error)----------------------------------------------#
  fit.cv.2 <- cv.glmnet(x=x, y=y, alpha=1, family="binomial", type.measure="class", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.5[[i]]$lasso.classerr$cv.fit <- fit.cv.2
  plot(fit.cv.2); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.5", "lasso", "classerr", "-", names(data)[i]))
  
  # Fit & CV (elnet deviance)-------------------------------------------------#
  fit.cv.3 <- cv.glmnet(x=x, y=y, alpha=0.9, family="binomial", type.measure="deviance", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.5[[i]]$elnet.deviance$cv.fit <- fit.cv.3
  plot(fit.cv.3); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.5", "elnet", "deviance", "-", names(data)[i]))
  
  # Fit & CV (elnet class error)----------------------------------------------#
  fit.cv.4 <- cv.glmnet(x=x, y=y, alpha=0.9, family="binomial", type.measure="class", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.5[[i]]$elnet.classerr$cv.fit <- fit.cv.4
  plot(fit.cv.4); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.5", "elnet", "classerr", "-", names(data)[i]))
  
  dev.off()
  
  # Then extract coefficiants, predictions and confussion stats from the CV objects
  for (j in 1:4) {
    
    fit.cv <- get(paste0("fit.cv.", j))
    
    #-----------Betas-------------------------------------#
    betas <- coef(fit.cv, s = "lambda.1se")
    betas <- as.matrix(betas)
    colnames(betas) <- levels(y)[2]
    
    betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% 
      as_tibble() %>% dplyr::filter(compID != "(Intercept)") %>% dplyr::filter_at(vars(-compID), any_vars(.!=0)) 
    nbetas.overall <- nrow(betas)
    nbetas.byclass <- betas %>% summarize_at(vars(-compID), .funs = function(x) sum(x != 0)) %>% unlist()
    
    glmnetCV$Myelo.5[[i]][[j]]$betas$nBetasOveral <- nbetas.overall
    glmnetCV$Myelo.5[[i]][[j]]$betas$nBetasByClass <- nbetas.byclass
    glmnetCV$Myelo.5[[i]][[j]]$betas$Betas <- betas
    
    #-----------Predictions & Confussion stats------------#
    prediction <- predict(fit.cv, newx=newx, s="lambda.1se", type="class") # a vector
    prediction <- factor(prediction, levels=levels(truth))
    
    prediction.p <- predict(fit.cv, newx=newx, s="lambda.1se", type="response")
    
    confussion <- confusionMatrix(data=prediction, reference=truth, positive=levels(truth)[2], dnn=c("Prediction", "Truth"), mode="everything") # positive - change it for two factors
    
    confussion.overal <- confussion$overall[1:2] # just Accuracy & Kappa
    confussion.byclass <- matrix(confussion$byClass[c("Balanced Accuracy", "F1")], nrow=1) # make it matrix for consistency!
    rownames(confussion.byclass) <- levels(truth)[2]
    
    glmnetCV$Myelo.5[[i]][[j]]$predicts$statsOveral <- confussion.overal
    glmnetCV$Myelo.5[[i]][[j]]$predicts$statsByClass$balancedAcc <- confussion.byclass[,1]
    glmnetCV$Myelo.5[[i]][[j]]$predicts$statsByClass$F1 <- confussion.byclass[,2]
    glmnetCV$Myelo.5[[i]][[j]]$predicts$predictClass <- prediction
    glmnetCV$Myelo.5[[i]][[j]]$predicts$predictProbs <- prediction.p
    glmnetCV$Myelo.5[[i]][[j]]$predicts$confussionMat <- confussion
    
  }
  
  print(paste("Myelo.5 -", names(data)[i], "... Done:", Sys.time()))
}


#-----------------------------------------------------------------------------#
### Myeloid 6. (classic monocytes vs. monocytes in Lung) (binomial)
#-----------------------------------------------------------------------------#

for (i in 1:4) {
  
  if (i <= 2) {
    A.df <- data[[i]] # best not to introduce factor levels here for easier tibble wrangling
    
    training <- A.df[inTraining, ] %>% dplyr::filter(celltype_tissue %in% c("classical monocyte - Lung", "monocyte - Lung"))
    testing <- A.df[-inTraining, ] %>% dplyr::filter(celltype_tissue %in% c("classical monocyte - Lung", "monocyte - Lung"))
    
    if (i == 1) {
      print("----- Myeloid 6 training -----")
      print(dplyr::count(training, celltype_tissue)); cat("\n")
      print("----- Myeloid 6 testing -----")
      print(dplyr::count(testing, celltype_tissue)); cat("\n")
      
      #-- Setup outcomes - they do not change per set (e.g. "Myelo.1)
      # training:
      y <- training$celltype_tissue %>% 
        fct_relevel(c("monocyte - Lung", "classical monocyte - Lung"))
      # testing:
      truth <- testing$celltype_tissue %>% 
        fct_relevel(c("monocyte - Lung", "classical monocyte - Lung"))
      
      #-- Setup myeloid subset cellnames for GE.mat subsetting
      inTraining.class.mo <- training$cell
      inTesting.class.mo <- testing$cell
    }
    
    x <- as.matrix(training[, as.character(1:100)])
    newx <- as.matrix(testing[, as.character(1:100)])
    
  } else if (i > 2) {
    GE.mat <- data[[i]]
    
    training <- GE.mat[inTraining.class.mo, ] # select only lung monocytes as above
    testing <- GE.mat[inTesting.class.mo, ]
    
    x <- training
    newx <- testing
    
  }
  
  # Fit & CV (lasso deviance)--------------------------------------------------#
  fit.cv.1 <- cv.glmnet(x=x, y=y, alpha=1, family="binomial", type.measure="deviance", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.6[[i]]$lasso.deviance$cv.fit <- fit.cv.1
  #--start plot
  pdf(file = paste0("results/fits/Myelo.6.CVplot.", names(data)[i], ".pdf"), width = 7.5, height = 6, pointsize = 10)
  par(mfrow=c(2,2))
  plot(fit.cv.1); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.6", "lasso", "deviance", "-", names(data)[i]))
  
  # Fit & CV (lasso class error)----------------------------------------------#
  fit.cv.2 <- cv.glmnet(x=x, y=y, alpha=1, family="binomial", type.measure="class", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.6[[i]]$lasso.classerr$cv.fit <- fit.cv.2
  plot(fit.cv.2); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.6", "lasso", "classerr", "-", names(data)[i]))
  
  # Fit & CV (elnet deviance)-------------------------------------------------#
  fit.cv.3 <- cv.glmnet(x=x, y=y, alpha=0.9, family="binomial", type.measure="deviance", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.6[[i]]$elnet.deviance$cv.fit <- fit.cv.3
  plot(fit.cv.3); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.6", "elnet", "deviance", "-", names(data)[i]))
  
  # Fit & CV (elnet class error)----------------------------------------------#
  fit.cv.4 <- cv.glmnet(x=x, y=y, alpha=0.9, family="binomial", type.measure="class", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Myelo.6[[i]]$elnet.classerr$cv.fit <- fit.cv.4
  plot(fit.cv.4); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Myelo.6", "elnet", "classerr", "-", names(data)[i]))
  
  dev.off()
  
  # Then extract coefficiants, predictions and confussion stats from the CV objects
  for (j in 1:4) {
    
    fit.cv <- get(paste0("fit.cv.", j))
    
    #-----------Betas-------------------------------------#
    betas <- coef(fit.cv, s = "lambda.1se")
    betas <- as.matrix(betas)
    colnames(betas) <- levels(y)[2]
    
    betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% 
      as_tibble() %>% dplyr::filter(compID != "(Intercept)") %>% dplyr::filter_at(vars(-compID), any_vars(.!=0)) 
    nbetas.overall <- nrow(betas)
    nbetas.byclass <- betas %>% summarize_at(vars(-compID), .funs = function(x) sum(x != 0)) %>% unlist()
    
    glmnetCV$Myelo.6[[i]][[j]]$betas$nBetasOveral <- nbetas.overall
    glmnetCV$Myelo.6[[i]][[j]]$betas$nBetasByClass <- nbetas.byclass
    glmnetCV$Myelo.6[[i]][[j]]$betas$Betas <- betas
    
    #-----------Predictions & Confussion stats------------#
    prediction <- predict(fit.cv, newx=newx, s="lambda.1se", type="class") # a vector
    prediction <- factor(prediction, levels=levels(truth))
    
    prediction.p <- predict(fit.cv, newx=newx, s="lambda.1se", type="response")
    
    confussion <- confusionMatrix(data=prediction, reference=truth, positive=levels(truth)[2], dnn=c("Prediction", "Truth"), mode="everything") # positive - change it for two factors
    
    confussion.overal <- confussion$overall[1:2] # just Accuracy & Kappa
    confussion.byclass <- matrix(confussion$byClass[c("Balanced Accuracy", "F1")], nrow=1) # make it matrix for consistency!
    rownames(confussion.byclass) <- levels(truth)[2]
    
    glmnetCV$Myelo.6[[i]][[j]]$predicts$statsOveral <- confussion.overal
    glmnetCV$Myelo.6[[i]][[j]]$predicts$statsByClass$balancedAcc <- confussion.byclass[,1]
    glmnetCV$Myelo.6[[i]][[j]]$predicts$statsByClass$F1 <- confussion.byclass[,2]
    glmnetCV$Myelo.6[[i]][[j]]$predicts$predictClass <- prediction
    glmnetCV$Myelo.6[[i]][[j]]$predicts$predictProbs <- prediction.p
    glmnetCV$Myelo.6[[i]][[j]]$predicts$confussionMat <- confussion
    
  }
  
  print(paste("Myelo.6 -", names(data)[i], "... Done:", Sys.time()))
}



#-----------------------------------------------------------------------------#
### Endothelial 1. (endothelial cells + others) (binomial)
#-----------------------------------------------------------------------------#

for (i in 1:4) {
  
  if (i <= 2) {
    A.df <- data[[i]] # best not to introduce factor levels here for easier tibble wrangling
    
    training <- A.df[inTraining, ]
    testing <- A.df[-inTraining, ]
    
    if (i == 1) {
      print("----- Endothelial 1 training -----")
      print(dplyr::count(training, endo_celltype)); cat("\n")
      print("----- Endothelial 1 testing -----")
      print(dplyr::count(testing, endo_celltype)); cat("\n")
      
      #-- Setup outcomes - they do not change per set (e.g. "Myelo.1)
      # training:
      y <- training$endo_celltype %>% 
        fct_relevel(c("non-endothelial", "endothelial cell")) %>% 
        fct_recode(other = "non-endothelial")
      # testing:
      truth <- testing$endo_celltype %>% 
        fct_relevel(c("non-endothelial", "endothelial cell")) %>% 
        fct_recode(other = "non-endothelial")
    }
    
    x <- as.matrix(training[, as.character(1:100)])
    newx <- as.matrix(testing[, as.character(1:100)])
    
  } else if (i > 2) {
    GE.mat <- data[[i]]
    
    training <- GE.mat[inTraining, ]
    testing <- GE.mat[-inTraining, ]
    
    x <- training
    newx <- testing
    
  }
  
  # Fit & CV (lasso deviance)--------------------------------------------------#
  fit.cv.1 <- cv.glmnet(x=x, y=y, alpha=1, family="binomial", type.measure="deviance", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Endo.1[[i]]$lasso.deviance$cv.fit <- fit.cv.1
  #--start plot
  pdf(file = paste0("results/fits/Endo.1.CVplot.", names(data)[i], ".pdf"), width = 7.5, height = 6, pointsize = 10)
  par(mfrow=c(2,2))
  plot(fit.cv.1); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Endo.1", "lasso", "deviance", "-", names(data)[i]))
  
  # Fit & CV (lasso class error)----------------------------------------------#
  fit.cv.2 <- cv.glmnet(x=x, y=y, alpha=1, family="binomial", type.measure="class", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Endo.1[[i]]$lasso.classerr$cv.fit <- fit.cv.2
  plot(fit.cv.2); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Endo.1", "lasso", "classerr", "-", names(data)[i]))
  
  # Fit & CV (elnet deviance)-------------------------------------------------#
  fit.cv.3 <- cv.glmnet(x=x, y=y, alpha=0.9, family="binomial", type.measure="deviance", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Endo.1[[i]]$elnet.deviance$cv.fit <- fit.cv.3
  plot(fit.cv.3); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Endo.1", "elnet", "deviance", "-", names(data)[i]))
  
  # Fit & CV (elnet class error)----------------------------------------------#
  fit.cv.4 <- cv.glmnet(x=x, y=y, alpha=0.9, family="binomial", type.measure="class", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Endo.1[[i]]$elnet.classerr$cv.fit <- fit.cv.4
  plot(fit.cv.4); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Endo.1", "elnet", "classerr", "-", names(data)[i]))
  
  dev.off()
  
  # Then extract coefficiants, predictions and confussion stats from the CV objects
  for (j in 1:4) {
    
    fit.cv <- get(paste0("fit.cv.", j))
    
    #-----------Betas-------------------------------------#
    betas <- coef(fit.cv, s = "lambda.1se")
    betas <- as.matrix(betas)
    colnames(betas) <- levels(y)[2]
    
    betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% 
      as_tibble() %>% dplyr::filter(compID != "(Intercept)") %>% dplyr::filter_at(vars(-compID), any_vars(.!=0)) 
    nbetas.overall <- nrow(betas)
    nbetas.byclass <- betas %>% summarize_at(vars(-compID), .funs = function(x) sum(x != 0)) %>% unlist()
    
    glmnetCV$Endo.1[[i]][[j]]$betas$nBetasOveral <- nbetas.overall
    glmnetCV$Endo.1[[i]][[j]]$betas$nBetasByClass <- nbetas.byclass
    glmnetCV$Endo.1[[i]][[j]]$betas$Betas <- betas
    
    #-----------Predictions & Confussion stats------------#
    prediction <- predict(fit.cv, newx=newx, s="lambda.1se", type="class") # a vector
    prediction <- factor(prediction, levels=levels(truth))
    
    prediction.p <- predict(fit.cv, newx=newx, s="lambda.1se", type="response")
    
    confussion <- confusionMatrix(data=prediction, reference=truth, positive=levels(truth)[2], dnn=c("Prediction", "Truth"), mode="everything") # positive - change it for two factors
    
    confussion.overal <- confussion$overall[1:2] # just Accuracy & Kappa
    confussion.byclass <- matrix(confussion$byClass[c("Balanced Accuracy", "F1")], nrow=1) # make it matrix for consistency!
    rownames(confussion.byclass) <- levels(truth)[2]
    
    glmnetCV$Endo.1[[i]][[j]]$predicts$statsOveral <- confussion.overal
    glmnetCV$Endo.1[[i]][[j]]$predicts$statsByClass$balancedAcc <- confussion.byclass[,1]
    glmnetCV$Endo.1[[i]][[j]]$predicts$statsByClass$F1 <- confussion.byclass[,2]
    glmnetCV$Endo.1[[i]][[j]]$predicts$predictClass <- prediction
    glmnetCV$Endo.1[[i]][[j]]$predicts$predictProbs <- prediction.p
    glmnetCV$Endo.1[[i]][[j]]$predicts$confussionMat <- confussion
    
  }
  
  print(paste("Endo.1 -", names(data)[i], "... Done:", Sys.time()))
}


#-----------------------------------------------------------------------------#
### Endothelial 2. (endothelial cells across tissues only)
#-----------------------------------------------------------------------------#

for (i in 1:4) {
  
  if (i <= 2) {
    A.df <- data[[i]] # best not to introduce factor levels here for easier tibble wrangling
    
    training <- A.df[inTraining, ] %>% dplyr::filter(str_detect(endo_celltype_tissue, "^endothelial cell - ")) 
    testing <- A.df[-inTraining, ] %>% dplyr::filter(str_detect(endo_celltype_tissue, "^endothelial cell - "))
    
    if (i == 1) {
      print("----- Endothelial 2 training -----")
      print(dplyr::count(training, endo_celltype_tissue)); cat("\n")
      print("----- Endothelial 2 testing -----")
      print(dplyr::count(testing, endo_celltype_tissue)); cat("\n")
      
      #-- Setup outcomes - they do not change per set (e.g. "Myelo.1)
      # training:
      y <- training$endo_celltype_tissue %>% 
        fct_relevel(c("endothelial cell - Heart", "endothelial cell - Lung", "endothelial cell - Kidney", "endothelial cell - Liver", 
                    "endothelial cell - Pancreas", "endothelial cell - Brain Non-Myeloid", "endothelial cell - Fat", 
                    "endothelial cell - Mammary Gland", "endothelial cell - Limb Muscle", "endothelial cell - Trachea")) %>% 
        fct_recode(`endothelial cell - Brain` = "endothelial cell - Brain Non-Myeloid")
      # testing:
      truth <- testing$endo_celltype_tissue %>% 
        fct_relevel(c("endothelial cell - Heart", "endothelial cell - Lung", "endothelial cell - Kidney", "endothelial cell - Liver", 
                      "endothelial cell - Pancreas", "endothelial cell - Brain Non-Myeloid", "endothelial cell - Fat", 
                      "endothelial cell - Mammary Gland", "endothelial cell - Limb Muscle", "endothelial cell - Trachea")) %>% 
        fct_recode(`endothelial cell - Brain` = "endothelial cell - Brain Non-Myeloid")
      
      #-- Setup myeloid subset cellnames for GE.mat subsetting
      inTraining.endo.tissues <- training$cell
      inTesting.endo.tissues <- testing$cell
    }
    
    x <- as.matrix(training[, as.character(1:100)])
    newx <- as.matrix(testing[, as.character(1:100)])
    
  } else if (i > 2) {
    GE.mat <- data[[i]]
    
    training <- GE.mat[inTraining.endo.tissues, ] # select only endothelial cells as above
    testing <- GE.mat[inTesting.endo.tissues, ]
    
    x <- training
    newx <- testing
    
  }
  
  # Fit & CV (lasso deviance)--------------------------------------------------#
  fit.cv.1 <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="deviance", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Endo.2[[i]]$lasso.deviance$cv.fit <- fit.cv.1
  #--start plot
  pdf(file = paste0("results/fits/Endo.2.CVplot.", names(data)[i], ".pdf"), width = 7.5, height = 6, pointsize = 10)
  par(mfrow=c(2,2))
  plot(fit.cv.1); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Endo.2", "lasso", "deviance", "-", names(data)[i]))
  
  # Fit & CV (lasso class error)----------------------------------------------#
  fit.cv.2 <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="class", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Endo.2[[i]]$lasso.classerr$cv.fit <- fit.cv.2
  plot(fit.cv.2); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Endo.2", "lasso", "classerr", "-", names(data)[i]))
  
  # Fit & CV (elnet deviance)-------------------------------------------------#
  fit.cv.3 <- cv.glmnet(x=x, y=y, alpha=0.9, family="multinomial", type.measure="deviance", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Endo.2[[i]]$elnet.deviance$cv.fit <- fit.cv.3
  plot(fit.cv.3); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Endo.2", "elnet", "deviance", "-", names(data)[i]))
  
  # Fit & CV (elnet class error)----------------------------------------------#
  fit.cv.4 <- cv.glmnet(x=x, y=y, alpha=0.9, family="multinomial", type.measure="class", nfolds=20, nlambda=200, parallel=TRUE) 
  glmnetCV$Endo.2[[i]]$elnet.classerr$cv.fit <- fit.cv.4
  plot(fit.cv.4); mtext(side=3, line=2.25, cex=1.2, col="blue", paste("Endo.2", "elnet", "classerr", "-", names(data)[i]))
  
  dev.off()
  
  # Then extract coefficiants, predictions and confussion stats from the CV objects
  for (j in 1:4) {
    
    fit.cv <- get(paste0("fit.cv.", j))
    
    #-----------Betas-------------------------------------#
    betas <- coef(fit.cv, s="lambda.1se")
    outcomes <- names(betas)
    betas <- purrr::map(betas, .f = as.matrix) %>% do.call(cbind, .)
    colnames(betas) <- outcomes
    
    betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% 
      as_tibble() %>% dplyr::filter(compID != "(Intercept)") %>% dplyr::filter_at(vars(-compID), any_vars(.!=0)) 
    nbetas.overall <- nrow(betas)
    nbetas.byclass <- betas %>% summarize_at(vars(-compID), .funs = function(x) sum(x != 0)) %>% unlist()
    
    glmnetCV$Endo.2[[i]][[j]]$betas$nBetasOveral <- nbetas.overall
    glmnetCV$Endo.2[[i]][[j]]$betas$nBetasByClass <- nbetas.byclass
    glmnetCV$Endo.2[[i]][[j]]$betas$Betas <- betas
    
    #-----------Predictions & Confussion stats------------#
    prediction <- predict(fit.cv, newx=newx, s="lambda.1se", type="class") # a vector
    prediction <- factor(prediction, levels=levels(truth))
    
    prediction.p <- predict(fit.cv, newx=newx, s="lambda.1se", type="response")
    
    confussion <- confusionMatrix(data=prediction, reference=truth, positive = NULL, dnn = c("Prediction", "Truth"), mode = "everything") # positive - change it for two factors
    
    confussion.overal <- confussion$overall[1:2] # just Accuracy & Kappa
    confussion.byclass <- confussion$byClass[, c("Balanced Accuracy", "F1")] # matrix!
    rownames(confussion.byclass) <- str_sub(rownames(confussion.byclass), start=8)
    
    glmnetCV$Endo.2[[i]][[j]]$predicts$statsOveral <- confussion.overal
    glmnetCV$Endo.2[[i]][[j]]$predicts$statsByClass$balancedAcc <- confussion.byclass[,1]
    glmnetCV$Endo.2[[i]][[j]]$predicts$statsByClass$F1 <- confussion.byclass[,2]
    glmnetCV$Endo.2[[i]][[j]]$predicts$predictClass <- prediction
    glmnetCV$Endo.2[[i]][[j]]$predicts$predictProbs <- prediction.p
    glmnetCV$Endo.2[[i]][[j]]$predicts$confussionMat <- confussion
    
  }
  
  print(paste("Endo.2 -", names(data)[i], "... Done:", Sys.time()))
}


stopCluster(cores)

saveRDS(glmnetCV, file = "results/glmnet_cv_predict_6myelo_2endo.rds")

###################################################################################################

rm(list=setdiff(ls(), c("data", "spikeF.sc", "sumF.sc", "spikeF.A.df", "sumF.A.df")))
# glmnetCV <- readRDS(file = here("results", "glmnet_cv_predict_6myelo_2endo.rds"))


# Create matrices with the key summary stats from glmnetCV object for visualisation and assesment across all 128 models

list.struct <- list("A_spikeF"=list("nBetas"=matrix(), "F1"=matrix(), "bAcc"=matrix()), 
                    "A_sumF"=list("nBetas"=matrix(), "F1"=matrix(), "bAcc"=matrix()), 
                    "GE_spikeF"=list("nBetas"=matrix(), "F1"=matrix(), "bAcc"=matrix()),
                    "GE_sumF"=list("nBetas"=matrix(), "F1"=matrix(), "bAcc"=matrix()))

glmnetStats <- list("Myelo.1" = list.struct, 
                    "Myelo.2" = list.struct, 
                    "Myelo.3" = list.struct, 
                    "Myelo.4" = list.struct, 
                    "Myelo.5" = list.struct, 
                    "Myelo.6" = list.struct, 
                    "Endo.1" = list.struct, 
                    "Endo.2" = list.struct)

for (analysis in c("Myelo.1","Myelo.2","Myelo.3","Myelo.4","Myelo.5","Myelo.6","Endo.1","Endo.2")) {
  
  for (dataset in c("A_spikeF","A_sumF","GE_spikeF","GE_sumF")) { 
    
    # nBetas
    `Lasso-deviance` <- c(glmnetCV[[analysis]][[dataset]]$lasso.deviance$betas$nBetasOveral, glmnetCV[[analysis]][[dataset]]$lasso.deviance$betas$nBetasByClass)
    `Lasso-classerror` <- c(glmnetCV[[analysis]][[dataset]]$lasso.classerr$betas$nBetasOveral, glmnetCV[[analysis]][[dataset]]$lasso.classerr$betas$nBetasByClass)
    `El.Net-deviance` <- c(glmnetCV[[analysis]][[dataset]]$elnet.deviance$betas$nBetasOveral, glmnetCV[[analysis]][[dataset]]$elnet.deviance$betas$nBetasByClass)
    `El.Net-classerror` <- c(glmnetCV[[analysis]][[dataset]]$elnet.classerr$betas$nBetasOveral, glmnetCV[[analysis]][[dataset]]$elnet.classerr$betas$nBetasByClass)
    
    mat <- cbind(`Lasso-deviance`, `Lasso-classerror`, `El.Net-deviance`, `El.Net-classerror`)
    rownames(mat) <- c("overall", rownames(mat)[2:nrow(mat)])
    glmnetStats[[analysis]][[dataset]]$nBetas <- mat
    
    
    # predictive performance: Kappa stat (overal), F1 (per class) 
    `Lasso-deviance` <- c(glmnetCV[[analysis]][[dataset]]$lasso.deviance$predicts$statsOveral[2], glmnetCV[[analysis]][[dataset]]$lasso.deviance$predicts$statsByClass$F1)
    `Lasso-classerror` <- c(glmnetCV[[analysis]][[dataset]]$lasso.classerr$predicts$statsOveral[2], glmnetCV[[analysis]][[dataset]]$lasso.classerr$predicts$statsByClass$F1)
    `El.Net-deviance` <- c(glmnetCV[[analysis]][[dataset]]$elnet.deviance$predicts$statsOveral[2], glmnetCV[[analysis]][[dataset]]$elnet.deviance$predicts$statsByClass$F1)
    `El.Net-classerror` <- c(glmnetCV[[analysis]][[dataset]]$elnet.classerr$predicts$statsOveral[2], glmnetCV[[analysis]][[dataset]]$elnet.classerr$predicts$statsByClass$F1)
    
    mat <- cbind(`Lasso-deviance`, `Lasso-classerror`, `El.Net-deviance`, `El.Net-classerror`)
    rownames(mat) <- c("Kappa (overall)", rownames(mat)[2:nrow(mat)]) # on plot needs left side label: F1 for "per class"
    glmnetStats[[analysis]][[dataset]]$F1 <- mat
    
    
    # predictive performance: Accuracy (overal), Balanced accuracy (per class) 
    `Lasso-deviance` <- c(glmnetCV[[analysis]][[dataset]]$lasso.deviance$predicts$statsOveral[1], glmnetCV[[analysis]][[dataset]]$lasso.deviance$predicts$statsByClass$balancedAcc)
    `Lasso-classerror` <- c(glmnetCV[[analysis]][[dataset]]$lasso.classerr$predicts$statsOveral[1], glmnetCV[[analysis]][[dataset]]$lasso.classerr$predicts$statsByClass$balancedAcc)
    `El.Net-deviance` <- c(glmnetCV[[analysis]][[dataset]]$elnet.deviance$predicts$statsOveral[1], glmnetCV[[analysis]][[dataset]]$elnet.deviance$predicts$statsByClass$balancedAcc)
    `El.Net-classerror` <- c(glmnetCV[[analysis]][[dataset]]$elnet.classerr$predicts$statsOveral[1], glmnetCV[[analysis]][[dataset]]$elnet.classerr$predicts$statsByClass$balancedAcc)
    
    mat <- cbind(`Lasso-deviance`, `Lasso-classerror`, `El.Net-deviance`, `El.Net-classerror`)
    rownames(mat) <- c("Accuracy (overall)", rownames(mat)[2:nrow(mat)]) # on plot needs left side label: Balanced accuracy for "per class"
    glmnetStats[[analysis]][[dataset]]$bAcc <- mat
    
    }
}

saveRDS(glmnetStats, file = "results/glmnet_nBeta_predict_stats_6myelo_2endo.rds")

###################################################################################################

# Plot the stats.. per analysis (6 x Myeloid, 2 x Endothelial)

# this are the set of plots focussing on number of active coefficients (modules in the model) and 
# cell type prediction performance stats based on confussion matrix - currently this is plotted for
# soemwhat less forgiving overall Kappa and per cell type F1 stats.. to this can be added somewhat more optimistic 
# overall accuracy and per cell type balanced accuracy stats (they are already extracted in glmnetStat object)
for (analysis in c("Myelo.1","Myelo.2","Myelo.3","Myelo.4","Myelo.5","Myelo.6","Endo.1","Endo.2")) { 
  
  #---nBetas
  mat_A <- cbind(glmnetStats[[analysis]]$A_spikeF$nBetas, glmnetStats[[analysis]]$A_sumF$nBetas)
  getRange_A <- range(mat_A[-1, , drop=F])
  nr <- nrow(mat_A[-1, , drop=F])
  mat_GE <- cbind(glmnetStats[[analysis]]$GE_spikeF$nBetas, glmnetStats[[analysis]]$GE_sumF$nBetas)
  getRange_GE <- range(mat_GE[-1, , drop=F])
  
  haA = HeatmapAnnotation(overall = anno_barplot(mat_A[1,], bar_width=0.5, height=unit(2, "cm")), 
                          annotation_name_side="left", annotation_name_rot=0, annotation_name_offset=unit(1.25, "cm")) 
  haGE = HeatmapAnnotation(overall = anno_barplot(mat_GE[1,], bar_width=0.5, height=unit(2, "cm")), show_annotation_name=F)
  
  col_fun_A = colorRamp2(breaks = c(0, getRange_A[2]), colors = c("white", "tomato"), space = "LAB") # "LAB", "RGB"
  cols_A <- col_fun_A(seq(0, getRange_A[2], by = 1))
  col_fun_GE = colorRamp2(breaks = c(0, getRange_GE[2]), colors = c("white", "purple"), space = "LAB") # "LAB", "RGB"
  cols_GE <- col_fun_GE(seq(0, getRange_GE[2], by = 2))
  
  htA <- Heatmap(mat_A[-1, , drop=F], name="nBeta\n(modules)", cluster_rows=F, cluster_columns=F, na_col="grey20", 
                 column_split=rep(c("Modules - spikeF norm.", "Modules - sumF norm."), each=4), 
                 cluster_column_slices=F, column_title_gp = gpar(fontsize=11, fontface="bold"), 
                 row_names_side="left", row_names_max_width = max_text_width(rownames(mat_A)), column_names_rot=45, column_names_gp = gpar(fontsize=10), 
                 col=cols_A, top_annotation=haA, column_gap = unit(15, "mm"), border=TRUE, 
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   if(j == 3) { 
                     grid.text(sprintf("%g", mat_A[-1, , drop=F][i, j]), x, y, gp = gpar(fontsize = 8))
                   } else if(j == 7) {
                     grid.text(sprintf("%g", mat_A[-1, , drop=F][i, j]), x, y, gp = gpar(fontsize = 8))
                   }
                 }) 
  htGE <- Heatmap(mat_GE[-1, , drop=F], name="nBeta\n(genes)", cluster_rows=F, cluster_columns=F, na_col="grey20", 
                  column_split=rep(c("Gene - spikeF norm.", "Genes - sumF norm."), each=4), 
                  cluster_column_slices=F, column_title_gp = gpar(fontsize=11, fontface="bold"), 
                  row_names_side="left", row_names_max_width = max_text_width(rownames(mat_GE)), column_names_rot=45, column_names_gp = gpar(fontsize=10), 
                  col=cols_GE, top_annotation=haGE, column_gap = unit(15, "mm"), border=TRUE, 
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    if(j == 3) { 
                      grid.text(sprintf("%g", mat_GE[-1, , drop=F][i, j]), x, y, gp = gpar(fontsize = 8))
                    } else if(j ==7) {
                      grid.text(sprintf("%g", mat_GE[-1, , drop=F][i, j]), x, y, gp = gpar(fontsize = 8))
                    }
                  }) # I filled cells with the stats for choosen el.net-deviance approach here.. 
  # below is a version with best stat for the cell type, for the version of data.. can be changed depending on the message
  
  htList = htA + htGE
  ht1 <- draw(htList, ht_gap = unit(25, "mm"))
  fname <- paste0("figures/", analysis, ".nBetas.pdf")
  ht1
  dev.copy(pdf, fname,width = 11.75, height = 0.2756*nr+2.25) # this works better in a loop
  dev.off()

  
  #---Kappa / F1 (these must be drawn as slices of one common heatplot to keep score colour consistent)
  mat_2 <- cbind(glmnetStats[[analysis]]$A_spikeF$F1, glmnetStats[[analysis]]$A_sumF$F1, glmnetStats[[analysis]]$GE_spikeF$F1, glmnetStats[[analysis]]$GE_sumF$F1)
  getRange_2 <- range(mat_2[-1, , drop=F], na.rm = T)
  nr <- nrow(mat_2[-1, , drop=F])
  
  ha.base <- 0.95*min(mat_2[1,], na.rm=T)
  ha2 = HeatmapAnnotation(`overall Kappa` = anno_barplot(mat_2[1,], baseline=ha.base, bar_width=0.5, height=unit(2, "cm")), 
                          annotation_name_side="left", annotation_name_rot=0, annotation_name_offset=unit(1.25, "cm")) 
  
  col_fun_2 = colorRamp2(breaks = c(0, 1), colors = c("dodgerblue", "white"), space = "LAB") # "LAB", "RGB"
  cols_2 <- col_fun_2(seq(0, 1, by = 0.05))
  
  ht2 <- Heatmap(mat_2[-1, , drop=F], name="F1  score", cluster_rows=F, cluster_columns=F, na_col="grey20", 
                 column_split=factor(rep(c("Modules - spikeF norm.", "Modules - sumF norm.", "Genes - spikeF norm.", "Genes - sumF norm."), each=4), 
                                     levels=c("Modules - spikeF norm.", "Modules - sumF norm.", "Genes - spikeF norm.", "Genes - sumF norm.")), 
                 cluster_column_slices=F, column_title_gp = gpar(fontsize=11, fontface="bold"), 
                 row_names_side="left", column_names_rot=45, column_names_gp = gpar(fontsize=10), 
                 col=cols_2, top_annotation=ha2, column_gap = unit(c(15,25,15), "mm"), border=TRUE, row_names_max_width = max_text_width(rownames(mat_2)), 
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   if(j == 3 & !is.na(mat_2[-1, , drop=F][i, j])) { 
                     grid.text(sprintf("%.2f", mat_2[-1, , drop=F][i, j]), x, y, gp = gpar(fontsize = 8))
                   } else if(j == 7 & !is.na(mat_2[-1, , drop=F][i, j])) {
                     grid.text(sprintf("%.2f", mat_2[-1, , drop=F][i, j]), x, y, gp = gpar(fontsize = 8))
                   } else if(j == 11 & !is.na(mat_2[-1, , drop=F][i, j])) {
                     grid.text(sprintf("%.2f", mat_2[-1, , drop=F][i, j]), x, y, gp = gpar(fontsize = 8))
                   } else if(j == 15 & !is.na(mat_2[-1, , drop=F][i, j])) {
                     grid.text(sprintf("%.2f", mat_2[-1, , drop=F][i, j]), x, y, gp = gpar(fontsize = 8))
                   }
                 }) # I filled cells with the stats for choosen el.net-deviance approach here.. 
  # below is a version with best stat for the cell type, for the version of data.. can be changed depending on the message
  
  ht2 <- draw(ht2)
  fname <- paste0("figures/", analysis, ".F1.pdf")
  ht2
  dev.copy(pdf, fname, width = 11.75, height = 0.2756*nr+2.25) # this works better in a loop
  dev.off()
  
  
  #---Overall / per cell balanced accuracy (these must be drawn as slices of one common heatplot to keep score colour consistent)
  mat_3 <- cbind(glmnetStats[[analysis]]$A_spikeF$bAcc, glmnetStats[[analysis]]$A_sumF$bAcc, glmnetStats[[analysis]]$GE_spikeF$bAcc, glmnetStats[[analysis]]$GE_sumF$bAcc)
  getRange_3 <- range(mat_3[-1, , drop=F], na.rm = T)
  nr <- nrow(mat_3[-1, , drop=F])
  
  ha.base <- 0.95*min(mat_3[1,], na.rm=T)
  ha3 = HeatmapAnnotation(`overall Accuracy` = anno_barplot(mat_3[1,], baseline=ha.base, bar_width=0.5, height=unit(2, "cm")), 
                          annotation_name_side="left", annotation_name_rot=0, annotation_name_offset=unit(1.25, "cm")) 
  
  col_fun_3 = colorRamp2(breaks = c(0, 1), colors = c("springgreen3", "white"), space = "LAB") # "LAB", "RGB"
  cols_3 <- col_fun_3(seq(0, 1, by = 0.05))
  
  ht3 <- Heatmap(mat_3[-1, , drop=F], name="Balanced\nAccuracy", cluster_rows=F, cluster_columns=F, na_col="grey20", 
                 column_split=factor(rep(c("Modules - spikeF norm.", "Modules - sumF norm.", "Genes - spikeF norm.", "Genes - sumF norm."), each=4), 
                                     levels=c("Modules - spikeF norm.", "Modules - sumF norm.", "Genes - spikeF norm.", "Genes - sumF norm.")), 
                 cluster_column_slices=F, column_title_gp = gpar(fontsize=11, fontface="bold"), 
                 row_names_side="left", column_names_rot=45, column_names_gp = gpar(fontsize=10), 
                 col=cols_3, top_annotation=ha3, column_gap = unit(c(15,25,15), "mm"), border=TRUE, row_names_max_width = max_text_width(rownames(mat_3)), 
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   if(j == 3) { 
                     grid.text(sprintf("%.2f", mat_3[-1, , drop=F][i, j]), x, y, gp = gpar(fontsize = 8))
                   } else if(j == 7) {
                     grid.text(sprintf("%.2f", mat_3[-1, , drop=F][i, j]), x, y, gp = gpar(fontsize = 8))
                   } else if(j == 11) {
                     grid.text(sprintf("%.2f", mat_3[-1, , drop=F][i, j]), x, y, gp = gpar(fontsize = 8))
                   } else if(j == 15) {
                     grid.text(sprintf("%.2f", mat_3[-1, , drop=F][i, j]), x, y, gp = gpar(fontsize = 8))
                   }
                 }) # I filled cells with the stats for choosen el.net-deviance approach here.. 
  # below is a version with best stat for the cell type, for the version of data.. can be changed depending on the message
  
  ht3 <- draw(ht3)
  fname <- paste0("figures/", analysis, ".Accuracy.pdf")
  ht3
  dev.copy(pdf, fname, width = 11.75, height = 0.2756*nr+2.25) # this works better in a loop
  dev.off()
  
  }

# # thise are versions of cell_fun that display the best stat values per cell type, per version of data for the plots above:
# # nBeta stats
# cell_fun = function(j, i, x, y, width, height, fill) {
#   if(j <=4 & mat_A[-1, , drop=F][i, j] == min(mat_A[-1, , drop=F][i, 1:4], na.rm=TRUE)) { 
#     grid.text(sprintf("%g", mat_A[-1, , drop=F][i, j]), x, y, gp = gpar(fontsize = 8))
#   } else if(j>4 & j<=8 & mat_A[-1, , drop=F][i, j] == min(mat_A[-1, , drop=F][i, 5:8], na.rm=TRUE)) {
#     grid.text(sprintf("%g", mat_A[-1, , drop=F][i, j]), x, y, gp = gpar(fontsize = 8))
#   }
# }
# # confussion mat stats
# cell_fun = function(j, i, x, y, width, height, fill) {
#   if(j <=4 & !is.na(mat_3[-1, , drop=F][i, j]) & mat_3[-1, , drop=F][i, j] == max(mat_3[-1, , drop=F][i, 1:4], na.rm=TRUE)) { 
#     grid.text(sprintf("%.2f", mat_3[-1, , drop=F][i, j]), x, y, gp = gpar(fontsize = 8))
#   } else if(j>4 & j<=8 & !is.na(mat_3[-1, , drop=F][i, j]) & mat_3[-1, , drop=F][i, j] == max(mat_3[-1, , drop=F][i, 5:8], na.rm=TRUE)) {
#     grid.text(sprintf("%.2f", mat_3[-1, , drop=F][i, j]), x, y, gp = gpar(fontsize = 8))
#   } else if(j>8 & j<=12 & !is.na(mat_3[-1, , drop=F][i, j]) & mat_3[-1, , drop=F][i, j] == max(mat_3[-1, , drop=F][i, 9:12], na.rm=TRUE)) {
#     grid.text(sprintf("%.2f", mat_3[-1, , drop=F][i, j]), x, y, gp = gpar(fontsize = 8))
#   } else if(j>12 & j<=16 & !is.na(mat_3[-1, , drop=F][i, j]) & mat_3[-1, , drop=F][i, j] == max(mat_3[-1, , drop=F][i, 13:16], na.rm=TRUE)) {
#     grid.text(sprintf("%.2f", mat_3[-1, , drop=F][i, j]), x, y, gp = gpar(fontsize = 8))
#   }
# }



# so after this you:
# 1. review all lots and decide on which data norm - model combination to chose: as you analyse write down potential ideas, visualise it as a result ! 
# 2. plot those matching coeffiinet heatplot (and UMAPS with the key cell-type specific components!.. and maybe matching A matrix median per cell type heatplots too)
# 3. start isolting S module genes (these are same for 8 analyses), 2 data norms?
# 4. gene-set go, reactome, motifs analysis + 
# 5. module genes ~ TF correlation and A scores ~ TF lasso model
# 6. overlef and citations...

# After assesing all the stats - conclusion to:
# Plot coefficent heatmaps for all 8 analyses.. but only sumF normalised data and Elastic.net-Deviance model 
# (add Elastic.net-Class.error model only for Myelo.6 and Endo.1 analyses - much sparser modeel & better predictions in these cases)
# this is going to reveal which modules to pick for which level of analysis (level of cell type discrimination)

# When it comes to S module genes - only focus on sumF normalised version of data 


###################################################################################################

# Plot the heatmaps of acive coefficients chosen by Elastic.net-Deviance model for all 8 analyses (only sumF normalised data)
# for each save one 1-100 comp no clustering version for easier comparison across analyses, and one version with clustering 

for (analysis in c("Myelo.1","Myelo.2","Myelo.3","Myelo.4","Myelo.5","Myelo.6","Endo.1","Endo.2")) {
  
  beta.tibb <- glmnetCV[[analysis]]$A_sumF$elnet.deviance$betas$Betas
  
  beta.mat <- beta.tibb %>% as.data.frame() %>% column_to_rownames(var="compID") %>% as.matrix()
  beta.mat <- t(beta.mat)
  nr <- nrow(beta.mat)
  nc <- ncol(beta.mat)
  
  col_fun = colorRamp2(breaks = c(min(beta.mat)+1, -0.1, 0.1, max(beta.mat)-1), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
  cols <- col_fun(seq(min(beta.mat), max(beta.mat), by = 0.5))
  
  if (nr > 1) {
    
    #---version with clustering 
    hm <- Heatmap(beta.mat, row_names_gp = gpar(fontsize = 10.5), column_names_gp = gpar(fontsize = 7),
                  row_names_max_width = max_text_width(rownames(beta.mat)),
                  column_dend_height = unit(0.9, "cm"), row_dend_width = unit(0.9, "cm"), name = "",
                  clustering_distance_rows="spearman", clustering_distance_columns="spearman",
                  clustering_method_rows="average", clustering_method_columns="average",
                  width = unit(0.3*nc, "cm"), height = unit(0.6*nr, "cm"), na_col = "grey20", col=cols, border=TRUE, 
                  show_heatmap_legend=TRUE, heatmap_legend_param = list(title=expression(beta), labels_gp = gpar(fontsize = 8)))
    hm <- draw(hm)
    w <- ComplexHeatmap:::width(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
    h <- ComplexHeatmap:::height(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
    
    fname <- paste0("figures/beta_heatmaps/", analysis, ".betas.clust.pdf")
    hm
    dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
    dev.off()
    
  }
  
  #---version with 1-100 ordered coefficients (no clustering) 
  hm <- Heatmap(beta.mat, row_names_gp = gpar(fontsize = 11), row_names_side="left", column_names_gp = gpar(fontsize = 7),
                row_names_max_width = max_text_width(rownames(beta.mat)), cluster_rows=F, cluster_columns=F,
                width = unit(0.3*nc, "cm"), height = unit(0.6*nr, "cm"), na_col = "grey20", col=cols, border=TRUE, 
                show_heatmap_legend=TRUE, heatmap_legend_param = list(title=expression(beta), labels_gp = gpar(fontsize = 8)))
  hm <- draw(hm)
  w <- ComplexHeatmap:::width(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
  h <- ComplexHeatmap:::height(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
  if (nr == 1) h <- h+0.75
  
  fname <- paste0("figures/beta_heatmaps/", analysis, ".betas.ordered.pdf")
  hm
  dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
  dev.off()
  
}


###################################################################################################

# Plot UMAPS with key components:

# (the list from the notes ...)

###################################################################################################

# Create and save component's median cell usage scores (A matrix) heatmaps
# Global (celltypes & celltypes-tissues) and for each analysis level (Myelo.1-6, Endo.1-2)
# (this useful for explaining and understanding +ve & -ve values from S & A matrix, and Lasso coefficients - and connecting with GE)

# define myeloid and endothelial types 
myelo.names <- c("basophil", "brain pericyte", "classical monocyte", "granulocyte", "granulocyte monocyte progenitor cell", 
                 "granulocytopoietic cell", "Kupffer cell", "macrophage", "microglial cell", "monocyte")
endo.names <- c("endothelial cell", "endothelial cell of hepatic sinusoid", "lung endothelial cell")

# Summarise sumF.A.df by cell type, then plot overall and for myelo & endo cells
summed.A <- sumF.A.df %>% 
  dplyr::select(cell, "1":"100", celltype, tissue) %>% 
  dplyr::mutate(celltype = if_else(str_detect(celltype, "large intestine$"), str_replace_all(celltype, "large intestine$", "lg. int."), celltype)) %>% 
  dplyr::mutate(celltype = if_else(str_detect(celltype, "tracheobronchial"), str_replace_all(celltype, "tracheobronchial", "trachbronch."), celltype)) %>% 
  dplyr::mutate(celltype = if_else(str_detect(celltype, "mammary gland"), str_replace_all(celltype, "mammary gland", "mammary gl."), celltype)) %>% 
  dplyr::mutate(celltype = if_else(str_detect(celltype, "-positive"), str_replace_all(celltype, "-positive", "+"), celltype)) %>% 
  dplyr::mutate(celltype = if_else(str_detect(celltype, "-negative"), str_replace_all(celltype, "-negative", "-"), celltype)) %>% 
  dplyr::mutate(tissue = if_else(str_detect(tissue, "^Brain "), "Brain", tissue)) %>% 
  tidyr::unite(celltype, tissue, col = celltype_tissue, sep = " . ", remove = F) %>% 
  dplyr::mutate(myelo_celltype = if_else(celltype %in% myelo.names, celltype, "non-myeloid")) %>% 
  dplyr::mutate(myelo_celltype = if_else(myelo_celltype == "classical monocyte", "monocyte", myelo_celltype)) %>% 
  dplyr::mutate(myelo_celltype = if_else(myelo_celltype == "granulocyte monocyte progenitor cell", "GMP", myelo_celltype)) %>% 
  dplyr::mutate(macrophage = if_else(myelo_celltype != "macrophage", "other", myelo_celltype)) %>% 
  tidyr::unite(myelo_celltype, tissue, col = myelo_celltype_tissue, sep = " . ", remove = F) %>% 
  dplyr::mutate(macrophage_tissue = if_else(str_detect(myelo_celltype_tissue, "^macrophage . ") | 
                                              str_detect(myelo_celltype_tissue, "^Kupffer cell . ") | 
                                              str_detect(myelo_celltype_tissue, "^microglial cell . "), myelo_celltype_tissue, "other")) %>% 
  dplyr::mutate(monocyte_tissue = if_else(str_detect(myelo_celltype_tissue, "^monocyte . "), myelo_celltype_tissue, "other")) %>% 
  dplyr::mutate(monocyte_subtype = if_else(str_detect(celltype_tissue, "^monocyte . Lung") | 
                                             str_detect(celltype_tissue, "^classical monocyte . Lung"), celltype_tissue, "other")) %>% 
  dplyr::mutate(endo_celltype = if_else(celltype %in% endo.names, "endothelial cell", "non-endothelial")) %>% 
  tidyr::unite(endo_celltype, tissue, col = endo_celltype_tissue, sep = " . ", remove = F) %>% 
  dplyr::mutate(endo_celltype_tissue = if_else(str_detect(endo_celltype_tissue, "^non-endothelial . "), "non-endothelial", endo_celltype_tissue)) %>% 
  dplyr::select(cell, celltype, tissue, celltype_tissue, myelo_celltype, macrophage, myelo_celltype_tissue, macrophage_tissue, 
                monocyte_tissue, monocyte_subtype, endo_celltype, endo_celltype_tissue, "1":"100")

# Overal summaries across all cells together:
summed.A %>% summarise_at(.vars=vars("1":"100"), .funs=mean) %>% gather() %>% pull(value) %>% sort() %>% round(2) # -0.24 - 0.76
summed.A %>% summarise_at(.vars=vars("1":"100"), .funs=sd) %>% gather() %>% pull(value) %>% sort() %>% round(3) # 0.069 - 0.415

summed.A %>% summarise_at(.vars=vars("1":"100"), .funs=median) %>% gather() %>% pull(value) %>% sort() %>% round(2) # -0.26 - 0.75
summed.A %>% summarise_at(.vars=vars("1":"100"), .funs=mad) %>% gather() %>% pull(value) %>% sort() %>% round(3) # 0.045 - 0.319

# e.g. patterns of component usage
hist(summed.A$`35`, breaks=100)
hist(summed.A$`56`, breaks=100)
hist(summed.A$`63`, breaks=100)
hist(summed.A$`4`, breaks=100)
hist(summed.A$`8`, breaks=100)
hist(summed.A$`14`, breaks=100)
hist(summed.A$`23`, breaks=100)


# Summarise across cell types - in few different ways...
# Plot the heatmaps of A cell usage scores accross celltypes - for all celltypes / tissues, and for all 8 analyses (only sumF normalised data)
# for each save one 1-100 comp no clustering version for easier comparison across analyses, and one version with clustering 
# stick to terminology CELL EMBEDDINGS / GENE LOADINGS (WEIGHTS)

#-----------------------------------------------------------------------------#
### Celltype
#-----------------------------------------------------------------------------#

summed.a.median <- summed.A %>% dplyr::select(celltype, "1":"100") %>% dplyr::group_by(celltype) %>% dplyr::summarise_all(median)
summed.a.mean <- summed.A %>% dplyr::select(celltype, "1":"100") %>% dplyr::group_by(celltype) %>% dplyr::summarise_all(mean)

# e.g.
par(mfrow=c(2,2))
hist(summed.a.median$`23`, breaks=50)
hist(summed.a.median$`35`, breaks=50)
hist(summed.a.median$`56`, breaks=50)
hist(summed.a.median$`63`, breaks=50)

par(mfrow=c(2,2))
hist(summed.a.mean$`23`, breaks=50)
hist(summed.a.mean$`35`, breaks=50)
hist(summed.a.mean$`56`, breaks=50)
hist(summed.a.mean$`63`, breaks=50)


summed.a.mat <- summed.a.mean %>% as.data.frame() %>% column_to_rownames(var="celltype") %>% as.matrix()

nr <- nrow(summed.a.mat)
nc <- ncol(summed.a.mat)

col_fun = colorRamp2(breaks = c(min(summed.a.mat)+0.2, -0.05, 0.05, max(summed.a.mat)-0.2), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
cols <- col_fun(seq(min(summed.a.mat), max(summed.a.mat), by = 0.1))

#---version with clustering 
hm <- Heatmap(summed.a.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(summed.a.mat)),
              column_dend_height = unit(0.9, "cm"), row_dend_width = unit(1.1, "cm"), name = "cell type\nembeddings\n(mean)",
              clustering_distance_rows="euclidean", clustering_distance_columns="euclidean",
              clustering_method_rows="ward.D", clustering_method_columns="ward.D",
              width = unit(0.25*nc, "cm"), height = unit(0.35*nr, "cm"), na_col = "grey20", col=cols, border=TRUE, 
              show_heatmap_legend=TRUE, heatmap_legend_param = list(title="cell type\nembeddings\n(mean)", labels_gp = gpar(fontsize = 7)))
hm <- draw(hm)
w <- ComplexHeatmap:::width(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(hm) %>% convertHeight(., "inch", valueOnly = TRUE)

fname <- paste0("figures/A_heatmaps/all.celltypes.clust.pdf")
hm
dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
dev.off()


#---version with clustering of cells only
hm <- Heatmap(summed.a.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(summed.a.mat)),
              column_dend_height = unit(0.9, "cm"), row_dend_width = unit(1.1, "cm"), name = "cell type\nembeddings\n(mean)",
              clustering_distance_rows="euclidean", cluster_columns=FALSE, 
              clustering_method_rows="ward.D", 
              width = unit(0.25*nc, "cm"), height = unit(0.35*nr, "cm"), na_col = "grey20", col=cols, border=TRUE, 
              show_heatmap_legend=TRUE, heatmap_legend_param = list(title="cell type\nembeddings\n(mean)", labels_gp = gpar(fontsize = 7)))
hm <- draw(hm)
w <- ComplexHeatmap:::width(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(hm) %>% convertHeight(., "inch", valueOnly = TRUE)

fname <- paste0("figures/A_heatmaps/all.celltypes.cell.clust.pdf")
hm
dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
dev.off()



#-----------------------------------------------------------------------------#
### Celltype - Tissue
#-----------------------------------------------------------------------------#

summed.a.mean <- summed.A %>% dplyr::select(celltype_tissue, "1":"100") %>% dplyr::group_by(celltype_tissue) %>% dplyr::summarise_all(mean)
summed.a.mat <- summed.a.mean %>% as.data.frame() %>% column_to_rownames(var="celltype_tissue") %>% as.matrix()

nr <- nrow(summed.a.mat)
nc <- ncol(summed.a.mat)

col_fun = colorRamp2(breaks = c(min(summed.a.mat)+0.2, -0.05, 0.05, max(summed.a.mat)-0.2), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
cols <- col_fun(seq(min(summed.a.mat), max(summed.a.mat), by = 0.1))

#---version with clustering 
hm <- Heatmap(summed.a.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(summed.a.mat)),
              column_dend_height = unit(0.9, "cm"), row_dend_width = unit(1.1, "cm"), name = "cell type\nembeddings\n(mean)",
              clustering_distance_rows="euclidean", clustering_distance_columns="euclidean",
              clustering_method_rows="ward.D", clustering_method_columns="ward.D",
              width = unit(0.25*nc, "cm"), height = unit(0.35*nr, "cm"), na_col = "grey20", col=cols, border=TRUE, 
              show_heatmap_legend=TRUE, heatmap_legend_param = list(title="cell type\nembeddings\n(mean)", labels_gp = gpar(fontsize = 7)))
hm <- draw(hm)
w <- ComplexHeatmap:::width(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(hm) %>% convertHeight(., "inch", valueOnly = TRUE)

fname <- paste0("figures/A_heatmaps/all.celltypes-tissues.clust.pdf")
hm
dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
dev.off()


#---version with clustering of cells only
hm <- Heatmap(summed.a.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(summed.a.mat)),
              column_dend_height = unit(0.9, "cm"), row_dend_width = unit(1.1, "cm"), name = "cell type\nembeddings\n(mean)",
              clustering_distance_rows="euclidean", cluster_columns=FALSE, 
              clustering_method_rows="ward.D", 
              width = unit(0.25*nc, "cm"), height = unit(0.35*nr, "cm"), na_col = "grey20", col=cols, border=TRUE, 
              show_heatmap_legend=TRUE, heatmap_legend_param = list(title="cell type\nembeddings\n(mean)", labels_gp = gpar(fontsize = 7)))
hm <- draw(hm)
w <- ComplexHeatmap:::width(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(hm) %>% convertHeight(., "inch", valueOnly = TRUE)

fname <- paste0("figures/A_heatmaps/all.celltypes-tissues.cell.clust.pdf")
hm
dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
dev.off()



#-----------------------------------------------------------------------------#
### Myeloid 1. & 2. (myeloid cell types + others)
#-----------------------------------------------------------------------------#

summed.a.mean <- summed.A %>% dplyr::select(myelo_celltype, "1":"100") %>% dplyr::group_by(myelo_celltype) %>% dplyr::summarise_all(mean)
summed.a.mat <- summed.a.mean %>% as.data.frame() %>% column_to_rownames(var="myelo_celltype") %>% as.matrix()

nr <- nrow(summed.a.mat)
nc <- ncol(summed.a.mat)

col_fun = colorRamp2(breaks = c(min(summed.a.mat)+0.1, -0.05, 0.05, max(summed.a.mat)-0.1), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
cols <- col_fun(seq(min(summed.a.mat), max(summed.a.mat), by = 0.1))

#---version with clustering 
hm <- Heatmap(summed.a.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(summed.a.mat)),
              column_dend_height = unit(0.9, "cm"), row_dend_width = unit(1.1, "cm"), name = "cell type\nembeddings\n(mean)",
              clustering_distance_rows="euclidean", clustering_distance_columns="euclidean",
              clustering_method_rows="ward.D", clustering_method_columns="ward.D",
              width = unit(0.25*nc, "cm"), height = unit(0.4*nr, "cm"), na_col = "grey20", col=cols, border=TRUE, 
              show_heatmap_legend=TRUE, heatmap_legend_param = list(title="cell type\nembeddings\n(mean)", labels_gp = gpar(fontsize = 7)))
hm <- draw(hm)
w <- ComplexHeatmap:::width(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(hm) %>% convertHeight(., "inch", valueOnly = TRUE)

fname <- paste0("figures/A_heatmaps/Myelo.1.clust.pdf")
hm
dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
dev.off()


#---version with clustering of cells only
hm <- Heatmap(summed.a.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(summed.a.mat)),
              column_dend_height = unit(0.9, "cm"), row_dend_width = unit(1.1, "cm"), name = "cell type\nembeddings\n(mean)",
              clustering_distance_rows="euclidean", cluster_columns=FALSE, 
              clustering_method_rows="ward.D", 
              width = unit(0.25*nc, "cm"), height = unit(0.4*nr, "cm"), na_col = "grey20", col=cols, border=TRUE, 
              show_heatmap_legend=TRUE, heatmap_legend_param = list(title="cell type\nembeddings\n(mean)", labels_gp = gpar(fontsize = 7)))
hm <- draw(hm)
w <- ComplexHeatmap:::width(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(hm) %>% convertHeight(., "inch", valueOnly = TRUE)

fname <- paste0("figures/A_heatmaps/Myelo.1.cell.clust.pdf")
hm
dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
dev.off()



#-----------------------------------------------------------------------------#
### Myeloid 3. (macrohages vs. ALL other, inc. myeloid non-macs) (binomial)
#-----------------------------------------------------------------------------#

summed.a.mean <- summed.A %>% dplyr::select(macrophage, "1":"100") %>% dplyr::group_by(macrophage) %>% dplyr::summarise_all(mean)
summed.a.mat <- summed.a.mean %>% as.data.frame() %>% column_to_rownames(var="macrophage") %>% as.matrix()

nr <- nrow(summed.a.mat)
nc <- ncol(summed.a.mat)

col_fun = colorRamp2(breaks = c(min(summed.a.mat)+0.0, -0.05, 0.05, max(summed.a.mat)-0.0), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
cols <- col_fun(seq(min(summed.a.mat), max(summed.a.mat), by = 0.1))

#---version with clustering 
hm <- Heatmap(summed.a.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(summed.a.mat)),
              column_dend_height = unit(0.9, "cm"), row_dend_width = unit(1.1, "cm"), name = "cell type\nembeddings\n(mean)",
              cluster_rows=FALSE, clustering_distance_columns="euclidean",
              clustering_method_rows="ward.D", clustering_method_columns="ward.D",
              width = unit(0.25*nc, "cm"), height = unit(0.5*nr, "cm"), na_col = "grey20", col=cols, border=TRUE, 
              show_heatmap_legend=TRUE, heatmap_legend_param = list(title="cell type\nembeddings\n(mean)", labels_gp = gpar(fontsize = 7)))
hm <- draw(hm)
w <- ComplexHeatmap:::width(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(hm) %>% convertHeight(., "inch", valueOnly = TRUE)

fname <- paste0("figures/A_heatmaps/Myelo.3.clust.pdf")
hm
dev.copy(pdf, fname, width = w, height = h+0.5) # this works better in a loop
dev.off()


#---version with clustering of cells only
hm <- Heatmap(summed.a.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(summed.a.mat)),
              column_dend_height = unit(0.9, "cm"), row_dend_width = unit(1.1, "cm"), name = "cell type\nembeddings\n(mean)",
              cluster_rows=FALSE, cluster_columns=FALSE, 
              width = unit(0.25*nc, "cm"), height = unit(0.5*nr, "cm"), na_col = "grey20", col=cols, border=TRUE, 
              show_heatmap_legend=TRUE, heatmap_legend_param = list(title="cell type\nembeddings\n(mean)", labels_gp = gpar(fontsize = 7)))
hm <- draw(hm)
w <- ComplexHeatmap:::width(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(hm) %>% convertHeight(., "inch", valueOnly = TRUE)

fname <- paste0("figures/A_heatmaps/Myelo.3.cell.clust.pdf")
hm
dev.copy(pdf, fname, width = w, height = h+0.5) # this works better in a loop
dev.off()



#-----------------------------------------------------------------------------#
### Myeloid 4. (macrophages across tissues only, inc. Kupfer & microglia too)
#-----------------------------------------------------------------------------#

summed.a.mean <- summed.A %>% dplyr::select(macrophage_tissue, "1":"100") %>% dplyr::group_by(macrophage_tissue) %>% dplyr::summarise_all(mean)
summed.a.mat <- summed.a.mean %>% as.data.frame() %>% column_to_rownames(var="macrophage_tissue") %>% as.matrix()

nr <- nrow(summed.a.mat)
nc <- ncol(summed.a.mat)

col_fun = colorRamp2(breaks = c(min(summed.a.mat)+0.05, -0.05, 0.05, max(summed.a.mat)-0.05), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
cols <- col_fun(seq(min(summed.a.mat), max(summed.a.mat), by = 0.1))

#---version with clustering 
hm <- Heatmap(summed.a.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(summed.a.mat)),
              column_dend_height = unit(0.9, "cm"), row_dend_width = unit(1.1, "cm"), name = "cell type\nembeddings\n(mean)",
              clustering_distance_rows="manhattan", clustering_distance_columns="manhattan",
              clustering_method_rows="ward.D2", clustering_method_columns="ward.D2",
              width = unit(0.25*nc, "cm"), height = unit(0.4*nr, "cm"), na_col = "grey20", col=cols, border=TRUE, 
              show_heatmap_legend=TRUE, heatmap_legend_param = list(title="cell type\nembeddings\n(mean)", labels_gp = gpar(fontsize = 7)))
hm <- draw(hm)
w <- ComplexHeatmap:::width(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(hm) %>% convertHeight(., "inch", valueOnly = TRUE)

fname <- paste0("figures/A_heatmaps/Myelo.4.clust.pdf")
hm
dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
dev.off()


#---version with clustering of cells only
hm <- Heatmap(summed.a.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(summed.a.mat)),
              column_dend_height = unit(0.9, "cm"), row_dend_width = unit(1.1, "cm"), name = "cell type\nembeddings\n(mean)",
              clustering_distance_rows="manhattan", cluster_columns=FALSE, 
              clustering_method_rows="ward.D2", 
              width = unit(0.25*nc, "cm"), height = unit(0.4*nr, "cm"), na_col = "grey20", col=cols, border=TRUE, 
              show_heatmap_legend=TRUE, heatmap_legend_param = list(title="cell type\nembeddings\n(mean)", labels_gp = gpar(fontsize = 7)))
hm <- draw(hm)
w <- ComplexHeatmap:::width(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(hm) %>% convertHeight(., "inch", valueOnly = TRUE)

fname <- paste0("figures/A_heatmaps/Myelo.4.cell.clust.pdf")
hm
dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
dev.off()



#-----------------------------------------------------------------------------#
### Myeloid 5. (monocytes across two tissues only) (binomial)
#-----------------------------------------------------------------------------#

summed.a.mean <- summed.A %>% dplyr::select(monocyte_tissue, "1":"100") %>% dplyr::group_by(monocyte_tissue) %>% dplyr::summarise_all(mean)
summed.a.mat <- summed.a.mean %>% as.data.frame() %>% column_to_rownames(var="monocyte_tissue") %>% as.matrix()

nr <- nrow(summed.a.mat)
nc <- ncol(summed.a.mat)

col_fun = colorRamp2(breaks = c(min(summed.a.mat)+0.0, -0.05, 0.05, max(summed.a.mat)-0.0), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
cols <- col_fun(seq(min(summed.a.mat), max(summed.a.mat), by = 0.1))

#---version with clustering 
hm <- Heatmap(summed.a.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(summed.a.mat)),
              column_dend_height = unit(0.9, "cm"), row_dend_width = unit(1.1, "cm"), name = "cell type\nembeddings\n(mean)",
              clustering_distance_rows="euclidean", clustering_distance_columns="euclidean",
              clustering_method_rows="ward.D", clustering_method_columns="ward.D",
              width = unit(0.25*nc, "cm"), height = unit(0.5*nr, "cm"), na_col = "grey20", col=cols, border=TRUE, 
              show_heatmap_legend=TRUE, heatmap_legend_param = list(title="cell type\nembeddings\n(mean)", labels_gp = gpar(fontsize = 7)))
hm <- draw(hm)
w <- ComplexHeatmap:::width(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(hm) %>% convertHeight(., "inch", valueOnly = TRUE)

fname <- paste0("figures/A_heatmaps/Myelo.5.clust.pdf")
hm
dev.copy(pdf, fname, width = w, height = h+0.5) # this works better in a loop
dev.off()


#---version with clustering of cells only
hm <- Heatmap(summed.a.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(summed.a.mat)),
              column_dend_height = unit(0.9, "cm"), row_dend_width = unit(1.1, "cm"), name = "cell type\nembeddings\n(mean)",
              clustering_distance_rows="euclidean", cluster_columns=FALSE, 
              clustering_method_rows="ward.D", 
              width = unit(0.25*nc, "cm"), height = unit(0.5*nr, "cm"), na_col = "grey20", col=cols, border=TRUE, 
              show_heatmap_legend=TRUE, heatmap_legend_param = list(title="cell type\nembeddings\n(mean)", labels_gp = gpar(fontsize = 7)))
hm <- draw(hm)
w <- ComplexHeatmap:::width(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(hm) %>% convertHeight(., "inch", valueOnly = TRUE)

fname <- paste0("figures/A_heatmaps/Myelo.5.cell.clust.pdf")
hm
dev.copy(pdf, fname, width = w, height = h+0.5) # this works better in a loop
dev.off()



#-----------------------------------------------------------------------------#
### Myeloid 6. (classic monocytes vs. monocytes in Lung) (binomial)
#-----------------------------------------------------------------------------#

summed.a.mean <- summed.A %>% dplyr::select(monocyte_subtype, "1":"100") %>% dplyr::group_by(monocyte_subtype) %>% dplyr::summarise_all(mean)
summed.a.mat <- summed.a.mean %>% as.data.frame() %>% column_to_rownames(var="monocyte_subtype") %>% as.matrix()

nr <- nrow(summed.a.mat)
nc <- ncol(summed.a.mat)

col_fun = colorRamp2(breaks = c(min(summed.a.mat)+0.0, -0.05, 0.05, max(summed.a.mat)-0.0), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
cols <- col_fun(seq(min(summed.a.mat), max(summed.a.mat), by = 0.1))

#---version with clustering 
hm <- Heatmap(summed.a.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(summed.a.mat)),
              column_dend_height = unit(0.9, "cm"), row_dend_width = unit(1.1, "cm"), name = "cell type\nembeddings\n(mean)",
              clustering_distance_rows="euclidean", clustering_distance_columns="euclidean",
              clustering_method_rows="ward.D", clustering_method_columns="ward.D",
              width = unit(0.25*nc, "cm"), height = unit(0.5*nr, "cm"), na_col = "grey20", col=cols, border=TRUE, 
              show_heatmap_legend=TRUE, heatmap_legend_param = list(title="cell type\nembeddings\n(mean)", labels_gp = gpar(fontsize = 7)))
hm <- draw(hm)
w <- ComplexHeatmap:::width(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(hm) %>% convertHeight(., "inch", valueOnly = TRUE)

fname <- paste0("figures/A_heatmaps/Myelo.6.clust.pdf")
hm
dev.copy(pdf, fname, width = w, height = h+0.5) # this works better in a loop
dev.off()


#---version with clustering of cells only
hm <- Heatmap(summed.a.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(summed.a.mat)),
              column_dend_height = unit(0.9, "cm"), row_dend_width = unit(1.1, "cm"), name = "cell type\nembeddings\n(mean)",
              clustering_distance_rows="euclidean", cluster_columns=FALSE, 
              clustering_method_rows="ward.D", 
              width = unit(0.25*nc, "cm"), height = unit(0.5*nr, "cm"), na_col = "grey20", col=cols, border=TRUE, 
              show_heatmap_legend=TRUE, heatmap_legend_param = list(title="cell type\nembeddings\n(mean)", labels_gp = gpar(fontsize = 7)))
hm <- draw(hm)
w <- ComplexHeatmap:::width(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(hm) %>% convertHeight(., "inch", valueOnly = TRUE)

fname <- paste0("figures/A_heatmaps/Myelo.6.cell.clust.pdf")
hm
dev.copy(pdf, fname, width = w, height = h+0.5) # this works better in a loop
dev.off()



#-----------------------------------------------------------------------------#
### Endothelial 1. (endothelial cells + others) (binomial)
#-----------------------------------------------------------------------------#

summed.a.mean <- summed.A %>% dplyr::select(endo_celltype, "1":"100") %>% dplyr::group_by(endo_celltype) %>% dplyr::summarise_all(mean)
summed.a.mat <- summed.a.mean %>% as.data.frame() %>% column_to_rownames(var="endo_celltype") %>% as.matrix()

nr <- nrow(summed.a.mat)
nc <- ncol(summed.a.mat)

col_fun = colorRamp2(breaks = c(min(summed.a.mat)+0.0, -0.05, 0.05, max(summed.a.mat)-0.0), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
cols <- col_fun(seq(min(summed.a.mat), max(summed.a.mat), by = 0.1))

#---version with clustering 
hm <- Heatmap(summed.a.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(summed.a.mat)),
              column_dend_height = unit(0.9, "cm"), row_dend_width = unit(1.1, "cm"), name = "cell type\nembeddings\n(mean)",
              cluster_rows=FALSE, clustering_distance_columns="euclidean",
              clustering_method_columns="ward.D",
              width = unit(0.25*nc, "cm"), height = unit(0.5*nr, "cm"), na_col = "grey20", col=cols, border=TRUE, 
              show_heatmap_legend=TRUE, heatmap_legend_param = list(title="cell type\nembeddings\n(mean)", labels_gp = gpar(fontsize = 7)))
hm <- draw(hm)
w <- ComplexHeatmap:::width(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(hm) %>% convertHeight(., "inch", valueOnly = TRUE)

fname <- paste0("figures/A_heatmaps/Endo.1.clust.pdf")
hm
dev.copy(pdf, fname, width = w, height = h+0.5) # this works better in a loop
dev.off()


#---version with clustering of cells only
hm <- Heatmap(summed.a.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(summed.a.mat)),
              column_dend_height = unit(0.9, "cm"), row_dend_width = unit(1.1, "cm"), name = "cell type\nembeddings\n(mean)",
              cluster_rows=FALSE, cluster_columns=FALSE, 
              width = unit(0.25*nc, "cm"), height = unit(0.5*nr, "cm"), na_col = "grey20", col=cols, border=TRUE, 
              show_heatmap_legend=TRUE, heatmap_legend_param = list(title="cell type\nembeddings\n(mean)", labels_gp = gpar(fontsize = 7)))
hm <- draw(hm)
w <- ComplexHeatmap:::width(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(hm) %>% convertHeight(., "inch", valueOnly = TRUE)

fname <- paste0("figures/A_heatmaps/Endo.1.cell.clust.pdf")
hm
dev.copy(pdf, fname, width = w, height = h+0.5) # this works better in a loop
dev.off()



#-----------------------------------------------------------------------------#
### Endothelial 2. (endothelial cells across tissues only)
#-----------------------------------------------------------------------------#

summed.a.mean <- summed.A %>% dplyr::select(endo_celltype_tissue, "1":"100") %>% dplyr::group_by(endo_celltype_tissue) %>% dplyr::summarise_all(mean)
summed.a.mat <- summed.a.mean %>% as.data.frame() %>% column_to_rownames(var="endo_celltype_tissue") %>% as.matrix()

nr <- nrow(summed.a.mat)
nc <- ncol(summed.a.mat)

col_fun = colorRamp2(breaks = c(min(summed.a.mat)+0.05, -0.05, 0.05, max(summed.a.mat)-0.05), colors = c("blue", "white", "white", "red"), space = "LAB") # "LAB", "RGB"
cols <- col_fun(seq(min(summed.a.mat), max(summed.a.mat), by = 0.1))

#---version with clustering 
hm <- Heatmap(summed.a.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(summed.a.mat)),
              column_dend_height = unit(0.9, "cm"), row_dend_width = unit(1.1, "cm"), name = "cell type\nembeddings\n(mean)",
              clustering_distance_rows="canberra", clustering_distance_columns="canberra",
              clustering_method_rows="ward.D2", clustering_method_columns="ward.D2",
              width = unit(0.25*nc, "cm"), height = unit(0.4*nr, "cm"), na_col = "grey20", col=cols, border=TRUE, 
              show_heatmap_legend=TRUE, heatmap_legend_param = list(title="cell type\nembeddings\n(mean)", labels_gp = gpar(fontsize = 7)))
hm <- draw(hm)
w <- ComplexHeatmap:::width(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(hm) %>% convertHeight(., "inch", valueOnly = TRUE)

fname <- paste0("figures/A_heatmaps/Endo.2.clust.pdf")
hm
dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
dev.off()


#---version with clustering of cells only
hm <- Heatmap(summed.a.mat, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 7),
              row_names_max_width = max_text_width(rownames(summed.a.mat)),
              column_dend_height = unit(0.9, "cm"), row_dend_width = unit(1.1, "cm"), name = "cell type\nembeddings\n(mean)",
              clustering_distance_rows="euclidean", cluster_columns=FALSE, 
              clustering_method_rows="ward.D", 
              width = unit(0.25*nc, "cm"), height = unit(0.4*nr, "cm"), na_col = "grey20", col=cols, border=TRUE, 
              show_heatmap_legend=TRUE, heatmap_legend_param = list(title="cell type\nembeddings\n(mean)", labels_gp = gpar(fontsize = 7)))
hm <- draw(hm)
w <- ComplexHeatmap:::width(hm) %>% convertHeight(., "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(hm) %>% convertHeight(., "inch", valueOnly = TRUE)

fname <- paste0("figures/A_heatmaps/Endo.2.cell.clust.pdf")
hm
dev.copy(pdf, fname, width = w, height = h) # this works better in a loop
dev.off()





###################################################################################################

# Then go to mudules_by_genes.R ... for work on genes level of modules 














