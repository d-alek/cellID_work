# Testing ica-like approach awith count modelling / factorisation 

library(tidyverse)
library(dplyr)
library(magrittr)
library(stringr)
library(ggplot2)
library(RColorBrewer)

library(Seurat)
library(sctransform)

# The idea is to have both FACS and DROP data normalised in two different ways (with two options within each method):
# 1. vst_log (relative counts, scaled, log transformed) --- could substitute it with csran() type size factor, log transform
# 2. vst_rc (relative counts, scaled, NOT log transformed) --- could substitute it with csran() type size factor, NO log transform
# 3. sct_res (SCTransform, residuals output)
# 4. sct_cnt (SCTransform, normalised counts output)

#-----------------------------------------------------------------------------#
# FACS dataset
# (set.seed 707: min 500 cells per onto = 31,438 cells; exclused 0 exp, ERCC, transgene = 22,921 genes )

facs.exp <- readRDS("~/projects/icatest/data/facs.exp.707.rds") # 31,438 cells; 22,921 genes (exclused: 0 exp, ERCC, transgene)
range(facs.exp@meta.data$nCount_RNA) # 50043 16887339
range(facs.exp@meta.data$nFeature_RNA) # 502 9910

# first remove genes expressed in less than 10 cells (these stats I previously tucked in the Seurat object)
sum(facs.exp[["RNA"]]@meta.features$RawCellGExp >= 10) # 20987
keep <- facs.exp[["RNA"]]@meta.features %>% as_tibble() %>% dplyr::filter(RawCellGExp >= 10) %>% dplyr::select(GeneSymbol) %>% pull()
facs.exp <- facs.exp[keep, ]

#---1. vst_log 
facs.exp.vst.log <- NormalizeData(facs.exp, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 1e6, verbose = T)
facs.exp.vst.log <- FindVariableFeatures(facs.exp.vst.log, assay = "RNA", selection.method = "vst", nfeatures = 7500, loess.span = 0.3, clip.max = "auto", mean.cutoff = c(0.1, 10), dispersion.cutoff = c(1, Inf), verbose = T)
saveRDS(facs.exp.vst.log, file = "~/projects/icatest/data/facs.exp.7500_vst.rds")

facs.exp.mat.VST.LOG <- as.matrix(GetAssayData(facs.exp.vst.log, assay = "RNA", slot = "data"))[VariableFeatures(facs.exp.vst.log), ] # 7500 31438
saveRDS(facs.exp.mat.VST.LOG, file = "~/projects/icatest/data/facs.exp.mat.7500.VST.LOG.rds")

#---2. vst_rc
facs.exp.vst.rc <- NormalizeData(facs.exp, assay = "RNA", normalization.method = "RC", scale.factor = 1e6, verbose = T) # method = "RC" for no log
facs.exp.vst.rc <- FindVariableFeatures(facs.exp.vst.rc, assay = "RNA", selection.method = "vst", nfeatures = 7500, loess.span = 0.3, clip.max = "auto", mean.cutoff = c(0.1, 10), dispersion.cutoff = c(1, Inf), verbose = T)

facs.exp.mat.VST.RC <- as.matrix(GetAssayData(facs.exp.vst.rc, assay = "RNA", slot = "data"))[VariableFeatures(facs.exp.vst.rc), ]
saveRDS(facs.exp.mat.VST.RC, file = "~/projects/icatest/data/facs.exp.mat.7500.VST.RC.rds")

#---3. sct_res
facs.exp.sct <- SCTransform(facs.exp, variable.features.n = 7500, return.only.var.genes = TRUE, n_genes = 5000, min_cells = 10, verbose = T) # residuals are cenetered but not scaled to unit variance 
# warnings: iteration limit reached ... CONSIDER CHANGING INTERNAL PARAMETERS [ ?sctransform::vst ]
saveRDS(facs.exp.sct, file = "~/projects/icatest/data/facs.exp.7500_sct.rds") # long run (1.5-2hr), so save it

facs.exp.mat.SCT.RES <- as.matrix(GetAssayData(facs.exp.sct, assay = "SCT", slot = "scale.data")) # order of genes is different!
saveRDS(facs.exp.mat.SCT.RES, file = "~/projects/icatest/data/facs.exp.mat.7500.SCT.RES.rds") # 7500 31438

#---4. sct_cnt
facs.exp.mat.SCT.CNT <- as.matrix(GetAssayData(facs.exp.sct, assay = "SCT", slot = "counts"))[VariableFeatures(facs.exp.sct), ]
saveRDS(facs.exp.mat.SCT.CNT, file = "~/projects/icatest/data/facs.exp.mat.7500.SCT.CNT.rds") 

facs.exp.mat.SCT.LOGC <- as.matrix(GetAssayData(facs.exp.sct, assay = "SCT", slot = "data"))[VariableFeatures(facs.exp.sct), ] 
saveRDS(facs.exp.mat.SCT.LOGC, file = "~/projects/icatest/data/facs.exp.mat.7500.SCT.LOGC.rds") 

# check for overlap
length(intersect(VariableFeatures(facs.exp.vst.log), VariableFeatures(facs.exp.sct))) # 5713 (out of 7500)

# plot vst - sct comparison
p1 <- VariableFeaturePlot(facs.exp.vst.log, log = NULL) %>% LabelPoints(points = head(VariableFeatures(facs.exp.vst.log), 20), repel = TRUE, xnudge = 0, ynudge = 0) + 
  ggtitle("Log normalisation, Variance stabilising transformation (SmartSeq2)") + theme(plot.title = element_text(hjust=0.5))
p2 <- VariableFeaturePlot(facs.exp.sct, log = NULL) %>% LabelPoints(points = head(VariableFeatures(facs.exp.sct), 20), repel = TRUE, xnudge = 0, ynudge = 0) + 
  ggtitle("Regularized negative binomial regression normalisation (SmartSeq2)") + theme(plot.title = element_text(hjust=0.5))
CombinePlots(plots = list(p1, p2), ncol = 1, legend = 'none')
ggsave("~/projects/icatest/output/facs/figures/VariableFeaturePlot_facs.png", width = 6, height = 8, units = "in", scale = 1.5)



#-----------------------------------------------------------------------------#
# 10X dataset
# (set.seed 707: min 500 cells per onto = 27,705 cells; exclused 0 exp, ERCC, transgene = 19,397 genes )

drop.exp <- readRDS("~/projects/icatest/data/drop.exp.707.rds") # 27,705 cells; 19,397 genes (exclused: 0 exp, ERCC, transgene)
range(drop.exp@meta.data$nCount_RNA) # 1002 140634
range(drop.exp@meta.data$nFeature_RNA) # 501 8121

# first remove genes expressed in less than 10 cells (these stats I previously tucked in the Seurat object)
sum(drop.exp[["RNA"]]@meta.features$RawCellGExp >= 10) # 17482
keep <- drop.exp[["RNA"]]@meta.features %>% as_tibble() %>% dplyr::filter(RawCellGExp >= 10) %>% dplyr::select(GeneSymbol) %>% pull()
drop.exp <- drop.exp[keep, ]

#---1. vst_log 
drop.exp.vst.log <- NormalizeData(drop.exp, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 1e4, verbose = T)
drop.exp.vst.log <- FindVariableFeatures(drop.exp.vst.log, assay = "RNA", selection.method = "vst", nfeatures = 7500, loess.span = 0.3, clip.max = "auto", mean.cutoff = c(0.1, 10), dispersion.cutoff = c(1, Inf), verbose = T)
saveRDS(drop.exp.vst.log, file = "~/projects/icatest/data/drop.exp.7500_vst.rds")

drop.exp.mat.VST.LOG <- as.matrix(GetAssayData(drop.exp.vst.log, assay = "RNA", slot = "data"))[VariableFeatures(drop.exp.vst.log), ] # 7500 27705
saveRDS(drop.exp.mat.VST.LOG, file = "~/projects/icatest/data/drop.exp.mat.7500.VST.LOG.rds")

#---2. vst_rc
drop.exp.vst.rc <- NormalizeData(drop.exp, assay = "RNA", normalization.method = "RC", scale.factor = 1e4, verbose = T) # method = "RC" for no log
drop.exp.vst.rc <- FindVariableFeatures(drop.exp.vst.rc, assay = "RNA", selection.method = "vst", nfeatures = 7500, loess.span = 0.3, clip.max = "auto", mean.cutoff = c(0.1, 10), dispersion.cutoff = c(1, Inf), verbose = T)

drop.exp.mat.VST.RC <- as.matrix(GetAssayData(drop.exp.vst.rc, assay = "RNA", slot = "data"))[VariableFeatures(drop.exp.vst.rc), ]
saveRDS(drop.exp.mat.VST.RC, file = "~/projects/icatest/data/drop.exp.mat.7500.VST.RC.rds")

#---3. sct_res
drop.exp.sct <- SCTransform(drop.exp, variable.features.n = 7500, return.only.var.genes = TRUE, n_genes = 5000, min_cells = 10, verbose = T) # residuals are cenetered but not scaled to unit variance 
# warnings: iteration limit reached ... CONSIDER CHANGING INTERNAL PARAMETERS [ ?sctransform::vst ]
saveRDS(drop.exp.sct, file = "~/projects/icatest/data/drop.exp.7500_sct.rds") # long run (1.5hr), so save it

drop.exp.mat.SCT.RES <- as.matrix(GetAssayData(drop.exp.sct, assay = "SCT", slot = "scale.data")) # order of genes is different!
saveRDS(drop.exp.mat.SCT.RES, file = "~/projects/icatest/data/drop.exp.mat.7500.SCT.RES.rds") # 7500 27705

#---4. sct_cnt
drop.exp.mat.SCT.CNT <- as.matrix(GetAssayData(drop.exp.sct, assay = "SCT", slot = "counts"))[VariableFeatures(drop.exp.sct), ]
saveRDS(drop.exp.mat.SCT.CNT, file = "~/projects/icatest/data/drop.exp.mat.7500.SCT.CNT.rds") 

drop.exp.mat.SCT.LOGC <- as.matrix(GetAssayData(drop.exp.sct, assay = "SCT", slot = "data"))[VariableFeatures(drop.exp.sct), ] 
saveRDS(drop.exp.mat.SCT.LOGC, file = "~/projects/icatest/data/drop.exp.mat.7500.SCT.LOGC.rds") 

# check for overlap
length(intersect(VariableFeatures(drop.exp.vst.log), VariableFeatures(drop.exp.sct))) # 5636 (out of 7500)

# plot vst - sct comparison
p1 <- VariableFeaturePlot(drop.exp.vst.log, log = NULL) %>% LabelPoints(points = head(VariableFeatures(drop.exp.vst.log), 20), repel = TRUE, xnudge = 0, ynudge = 0) + 
  ggtitle("Log normalisation, Variance stabilising transformation (10X - UMI counts)") + theme(plot.title = element_text(hjust=0.5))
p2 <- VariableFeaturePlot(drop.exp.sct, log = NULL) %>% LabelPoints(points = head(VariableFeatures(drop.exp.sct), 20), repel = TRUE, xnudge = 0, ynudge = 0) + 
  ggtitle("Regularized negative binomial regression normalisation (10X - UMI counts)") + theme(plot.title = element_text(hjust=0.5))
CombinePlots(plots = list(p1, p2), ncol = 1, legend = 'none')
ggsave("~/projects/icatest/output/drop/figures/VariableFeaturePlot_drop.png", width = 6, height = 8, units = "in", scale = 1.5)



# Check for overlap between datasets:
length(Reduce(intersect, list(a = VariableFeatures(facs.exp.vst.log), b = VariableFeatures(facs.exp.sct), 
                              c = VariableFeatures(drop.exp.vst.log), d = VariableFeatures(drop.exp.sct)))) # 3296

length(intersect(VariableFeatures(facs.exp.vst.log), VariableFeatures(drop.exp.vst.log))) # 4318
length(intersect(VariableFeatures(facs.exp.vst.log), VariableFeatures(drop.exp.sct))) # 4337



# later test the overlap with those used under 3sd in 100 ncomp run ... !!!



# After ICA runs are done (on the cluster) for "icafast: 8 test dataversions"
# attach regrouped cell type annotation and other cell metadata to A matrices
# this task is different for FACS and DROP datasets, but same for all versions of FACS/DROP (same cells)!

ica_100_8test <- readRDS("~/projects/icatest/output/ica_100_8test.rds")
names(ica_100_8test)
rm( list = setdiff(ls(), "ica_100_8test") )



###############################################################################
#####                   A matrix side of things - CELLS                   #####
###############################################################################


### FACS workflow ----------------------------------------------------------###

# representative FACS dataset (for annos and metadata)
exp <- readRDS(file = "~/projects/icatest/data/facs.exp.7500_vst.rds")
exp.mat <- readRDS(file = "~/projects/icatest/data/facs.exp.mat.7500.VST.LOG.rds")
exp.df <- as.data.frame(t(exp.mat)) %>% rownames_to_column(var = "cellnames")

all.equal(colnames(exp.mat), colnames(exp)) 
all.equal(exp.df$cellnames, colnames(exp)) 
all.equal(exp.df$cellnames, rownames(exp@meta.data)) 

# representative FACS A matrix (just for annos and metadata)
A.100.df <- ica_100_8test[[1]]$M %>% as.data.frame()
colnames(A.100.df) <- as.character(1:100)
rownames(A.100.df) <- exp.df$cellnames
A.100.df <- A.100.df %>% rownames_to_column(var = "cellnames") %>% as_tibble()

# Attach annos and meta to A.100.df - then transfer them to other A data frames
# first, colect useful meta and fix fix cell annotation strings
dim(exp@meta.data); colnames(exp@meta.data)
exp@meta.data %>% as_tibble() %>% count(cell_ontology_class) %>% print(n=100)

meta <- exp@meta.data %>% 
  rownames_to_column(var = "cellnames") %>% 
  as_tibble() %>% 
  dplyr::select(cellnames, cell_ontology_class, tissue, mouse.id, mouse.sex, plate.barcode, nCount_RNA, nFeature_RNA) %>% 
  dplyr::rename(celltype = cell_ontology_class) %>% 
  dplyr::mutate(celltype = if_else(celltype == "", "unknown", celltype))

meta %>% count(celltype) %>% print(n=100)
# reorganise and delete ambigous cells annotations - ONE OPTION IS TO DO THIS BEFORE ICA RUNS ? should have sampled more endothelial cells too !
myelo.names <- c("basophil", "brain pericyte", "classical monocyte", "granulocyte", "Kupffer cell", "macrophage", "microglial cell", "monocyte")
endo.names <- c("endothelial cell", "endothelial cell of hepatic sinusoid", "lung endothelial cell") # "endothelial cell" in 10 tissues
# NOTE: I included "brain pericytes" here as "myeloid".. published, although they are not the classic myeloid cell.. interesting to explore

meta.attach <- meta %>% 
  dplyr::filter(!(celltype %in% c("unknown", "blood_cell", "leukocyte", "myeloid_cell", "professional_antigen_presenting_cell", "granulocytopoietic_cell", "granulocyte monocyte progenitor cell", "hematopoietic precursor cell"))) %>% 
  dplyr::mutate(celltype_my = if_else(celltype %in% myelo.names, "myeloid", "non myeloid")) %>% 
  dplyr::mutate(celltype_my_type = if_else(celltype %in% myelo.names, celltype, "non myeloid")) %>% 
  dplyr::mutate(celltype_my_type = case_when(celltype_my_type == "classical monocyte" ~ "monocyte", 
                                             TRUE ~ as.character(celltype_my_type))) %>% 
  tidyr::unite(celltype_my_type, tissue, col = celltype_my_type_tissue, sep = " - ", remove = F) %>% 
  dplyr::mutate(celltype_my_type_tissue = if_else(str_detect(celltype_my_type_tissue, "^non myeloid -"), "non myeloid", celltype_my_type_tissue)) %>% 
  dplyr::mutate(celltype_endo = if_else(celltype %in% endo.names, "endothelial", "non endothelial")) %>% 
  tidyr::unite(celltype_endo, tissue, col = celltype_endo_tissue, sep = " - ", remove = F) %>% 
  dplyr::mutate(celltype_endo_tissue = if_else(str_detect(celltype_endo_tissue, "^non endothelial -"), "non endothelial", celltype_endo_tissue))

# some checks
meta.attach %>% count(celltype_my); meta.attach %>% count(celltype_my_type); meta.attach %>% count(celltype_my_type_tissue) %>% print(n=100); 
meta.attach %>% count(celltype_endo); meta.attach %>% count(celltype_endo_tissue)

meta.attach <- meta.attach %>% 
  dplyr::mutate(celltype_my = factor(celltype_my, levels = c("myeloid", "non myeloid"))) %>% 
  dplyr::mutate(celltype_my_type = factor(celltype_my_type, levels = c("granulocyte", "basophil", "monocyte", "macrophage", "Kupffer cell", "microglial cell", "brain pericyte", "non myeloid"))) %>% 
  dplyr::mutate(celltype_my_type_tissue = factor(celltype_my_type_tissue, levels = c("granulocyte - Marrow", "basophil - Marrow", "monocyte - Marrow", "monocyte - Lung", 
                                                                                     "macrophage - Marrow", "macrophage - Spleen", "macrophage - Kidney", "macrophage - Limb_Muscle", "macrophage - Diaphragm", "macrophage - Brain_Myeloid",   
                                                                                     "Kupffer cell - Liver", "microglial cell - Brain_Myeloid", "brain pericyte - Brain_Non-Myeloid", "non myeloid"))) %>% 
  dplyr::mutate(celltype_endo = factor(celltype_endo, levels = c("endothelial", "non endothelial"))) %>% 
  dplyr::mutate(celltype_endo_tissue = factor(celltype_endo_tissue, levels = c("endothelial - Aorta", "endothelial - Heart", "endothelial - Lung", "endothelial - Kidney", "endothelial - Liver", "endothelial - Pancreas", "endothelial - Mammary_Gland", 
                                                                              "endothelial - Fat", "endothelial - Limb_Muscle", "endothelial - Diaphragm", "endothelial - Trachea", "endothelial - Brain_Non-Myeloid", "non endothelial")))

meta.attach <- meta.attach %>% 
  dplyr::select("cellnames", "celltype", "tissue", "celltype_my", "celltype_my_type", "celltype_my_type_tissue", "celltype_endo", "celltype_endo_tissue", "mouse.id", "mouse.sex", "plate.barcode", "nCount_RNA", "nFeature_RNA")

# make a universal DROP keep vector for all A matrices
keep <- meta.attach$cellnames

# attach annos and meta to all 8 A matrices and keep them in a list
ica_100_8test_A_facs <- list()
for(i in 1:4){
  A <- ica_100_8test[[i]]$M %>% as.data.frame()
  colnames(A) <- as.character(1:100)
  rownames(A) <- exp.df$cellnames
  A <- A %>% 
    rownames_to_column(var = "cellnames") %>% as_tibble() %>% 
    dplyr::filter(cellnames %in% keep) %>% 
    left_join(y = meta.attach, by = "cellnames")
  
  ica_100_8test_A_facs[[i]] <- A
}
names(ica_100_8test_A_facs) <- names(ica_100_8test)[1:4]
saveRDS(ica_100_8test_A_facs, "~/projects/icatest/output/test8/ica_100_8test_A_facs.rds")

rm( list = setdiff(ls(), "ica_100_8test") )



### DROP workflow ----------------------------------------------------------###

# representative DROP Seurat object (for annos and metadata)
exp <- readRDS(file = "~/projects/icatest/data/drop.exp.7500_vst.rds")
exp.mat <- readRDS(file = "~/projects/icatest/data/drop.exp.mat.7500.VST.LOG.rds")
exp.df <- as.data.frame(t(exp.mat)) %>% rownames_to_column(var = "cellnames")

all.equal(colnames(exp.mat), colnames(exp)) 
all.equal(exp.df$cellnames, colnames(exp)) 
all.equal(exp.df$cellnames, rownames(exp@meta.data)) 

# representative DROP A matrix (just for annos and metadata)
A.100.df <- ica_100_8test[[5]]$M %>% as.data.frame()
colnames(A.100.df) <- as.character(1:100)
rownames(A.100.df) <- exp.df$cellnames
A.100.df <- A.100.df %>% rownames_to_column(var = "cellnames") %>% as_tibble()

# Attach annos and meta to A.100.df - then transfer them to other A data frames
# first, colect useful meta and fix fix cell annotation strings
dim(exp@meta.data); colnames(exp@meta.data)
exp@meta.data %>% as_tibble() %>% count(cell_ontology_class) %>% print(n=100)

meta <- exp@meta.data %>% 
  rownames_to_column(var = "cellnames") %>% 
  as_tibble() %>% 
  dplyr::select(cellnames, cell_ontology_class, tissue, mouse.id, mouse.sex, nCount_RNA, nFeature_RNA) %>% 
  dplyr::rename(celltype = cell_ontology_class) %>% 
  dplyr::mutate(celltype = if_else(celltype == "", "unknown", celltype))

meta %>% count(celltype) %>% print(n=100)
# reorganise and delete ambigous cells annotations - ONE OPTION IS TO DO THIS BEFORE ICA RUNS ? should have sampled more endothelial cells too !
myelo.names <- c("alveolar macrophage", "basophil", "classical monocyte", "dendritic cell", "granulocyte", "Langerhans cell", "macrophage", "mast cell", "monocyte", "non-classical monocyte", "promonocyte")
endo.names <- c("endothelial cell", "endothelial cell of hepatic sinusoid", "kidney capillary endothelial cell", "lung endothelial cell") # "endothelial cell" in 5 tissues

meta.attach <- meta %>% 
  dplyr::filter(!(celltype %in% c("unknown", "blood cell", "leukocyte", "myeloid cell", "granulocytopoietic cell", "hematopoietic precursor cell"))) %>% 
  dplyr::mutate(celltype_my = if_else(celltype %in% myelo.names, "myeloid", "non myeloid")) %>% 
  dplyr::mutate(celltype_my_type = if_else(celltype %in% myelo.names, celltype, "non myeloid")) %>% 
  dplyr::mutate(celltype_my_type = case_when(celltype_my_type == "classical monocyte" ~ "monocyte", 
                                             celltype_my_type == "non-classical monocyte" ~ "monocyte", 
                                             celltype_my_type == "alveolar macrophage" ~ "macrophage", 
                                             TRUE ~ as.character(celltype_my_type))) %>% 
  tidyr::unite(celltype_my_type, tissue, col = celltype_my_type_tissue, sep = " - ", remove = F) %>% 
  dplyr::mutate(celltype_my_type_tissue = if_else(str_detect(celltype_my_type_tissue, "^non myeloid -"), "non myeloid", celltype_my_type_tissue)) %>% 
  dplyr::mutate(celltype_endo = if_else(celltype %in% endo.names, "endothelial", "non endothelial")) %>% 
  tidyr::unite(celltype_endo, tissue, col = celltype_endo_tissue, sep = " - ", remove = F) %>% 
  dplyr::mutate(celltype_endo_tissue = if_else(str_detect(celltype_endo_tissue, "^non endothelial -"), "non endothelial", celltype_endo_tissue))

# some checks
meta.attach %>% count(celltype_my); meta.attach %>% count(celltype_my_type); meta.attach %>% count(celltype_my_type_tissue) %>% print(n=100)
meta.attach %>% count(celltype_endo); meta.attach %>% count(celltype_endo_tissue)

meta.attach <- meta.attach %>% 
  dplyr::mutate(celltype_my = factor(celltype_my, levels = c("myeloid", "non myeloid"))) %>% 
  dplyr::mutate(celltype_my_type = factor(celltype_my_type, levels = c("granulocyte", "basophil", "promonocyte", "monocyte", "macrophage", "dendritic cell", "Langerhans cell", "mast cell", "non myeloid"))) %>% 
  dplyr::mutate(celltype_my_type_tissue = factor(celltype_my_type_tissue, levels = c("granulocyte - Marrow", "basophil - Marrow", "promonocyte - Marrow", "monocyte - Marrow", "monocyte - Lung", 
                                                                                     "macrophage - Marrow", "macrophage - Spleen", "macrophage - Lung", "macrophage - Kidney", "macrophage - Mammary_Gland", "macrophage - Limb_Muscle",  
                                                                                     "dendritic cell - Spleen", "Langerhans cell - Tongue", "mast cell - Lung", "non myeloid"))) %>% 
  dplyr::mutate(celltype_endo = factor(celltype_endo, levels = c("endothelial", "non endothelial"))) %>% 
  dplyr::mutate(celltype_endo_tissue = factor(celltype_endo_tissue, levels = c("endothelial - Heart_and_Aorta", "endothelial - Lung", "endothelial - Kidney", "endothelial - Liver", "endothelial - Mammary_Gland", 
                                                                               "endothelial - Limb_Muscle", "endothelial - Trachea", "endothelial - Bladder" , "non endothelial")))

meta.attach <- meta.attach %>% 
  dplyr::select("cellnames", "celltype", "tissue", "celltype_my", "celltype_my_type", "celltype_my_type_tissue", "celltype_endo", "celltype_endo_tissue", "mouse.id", "mouse.sex", "nCount_RNA", "nFeature_RNA")

# make a universal DROP keep vector for all A matrices
keep <- meta.attach$cellnames

# attach annos and meta to all 8 A matrices and keep them in a list
ica_100_8test_A_drop <- list()
for(i in 1:4){
  A <- ica_100_8test[[i+4]]$M %>% as.data.frame()
  colnames(A) <- as.character(1:100)
  rownames(A) <- exp.df$cellnames
  A <- A %>% 
    rownames_to_column(var = "cellnames") %>% as_tibble() %>% 
    dplyr::filter(cellnames %in% keep) %>% 
    left_join(y = meta.attach, by = "cellnames")
  
  ica_100_8test_A_drop[[i]] <- A
}
names(ica_100_8test_A_drop) <- names(ica_100_8test)[5:8]
saveRDS(ica_100_8test_A_drop, "~/projects/icatest/output/test8/ica_100_8test_A_drop.rds")



#-----------------------------------------------------------------------------#

# Now you can go on and do lassoo etc...
rm( list = setdiff(ls(), c("ica_100_8test", "ica_100_8test_A_facs")) )

#-----------------------------------------------------------------------------#


# LASSO GLM options to consider:
#-----------------------------------------------#
# 1. bi- vs. multi-nomial
# 2. pure lasso vs. "strong" elastic net
# 3. global then tissue-specific on macs
# 4. with A matrix values vs. gene expression values !
# 5. Group-Overlap-LASSO: with module genes at different sd (this one would be apart - try "grpregOverlap" (extension of "grpreg") or "mlgl" packages)!

# Below comparisons: 

### Ia: bi-nomial, level A: myeloid / non myeloid (celltype_my), .A mat values
### Ib: multi-nomial, level B: myeloid cell types (celltype_my_type), .A mat values
### Ic: multi-nomial, level C: macrophages accross tissues (celltype_my_type_tissue), .A mat values
### Inc: multi-nomial, non-chierarchical (myeloid cell types and non myeloid) (celltype_my_type), .A mat values


library("glmnet")
# try parallel
library(parallel)
detectCores() # 4
detectCores(logical = FALSE) # 2 (2 physical, but can do 4...)
library(doParallel)
cores <- makeCluster(detectCores(logical = FALSE), type='PSOCK') # fine to just put 2 or 4 in braces



############################################################################
### FACS lasso ----------------------------------------------------------###
############################################################################


#---------Split to Training and Testing datasets----------#
library(caret)

# These cell indices are same for Gen.Exp. Df and all A matrices
# Best to make them based on use_celltype_tissue variable (to keep celltype-tissue representation)!
set.seed(707)
inTrain.tissue <- createDataPartition(y = ica_100_8test_A_facs[[1]]$celltype_my_type_tissue, p = .85, list = FALSE) 
saveRDS(inTrain.tissue, file = "~/projects/icatest/output/test8/inTrain_facs.rds")
#---------------------------------------------------------#


################### Create loops for just saving fit.cv #######################

# get the right data split
inTrain.tissue <- readRDS("~/projects/icatest/output/test8/inTrain_facs.rds")

# setup a list of cv results
ica_100_8test_glmnetCV_facs <- list() # each with CVfit object for four data versions (1:"facs.vst.log" 2:"facs.vst.rc"  3:"facs.sct.res" 4:"facs.sct.cnt")

#-----------------------------------------------------------------------------#
### Ia: bi-nomial, level A: myeloid / non myeloid (celltype_my), .A mat values
#-----------------------------------------------------------------------------#

cores <- makeCluster(detectCores(logical = FALSE), type='PSOCK')
registerDoParallel(cores)

for (i in 1:4) {
  
  a <- ica_100_8test_A_facs[[i]]
  
  # global
  training <- a[ inTrain.tissue, ]
  testing  <- a[-inTrain.tissue, ]
  
  print(dplyr::count(training, celltype_my))
  print(dplyr::count(testing, celltype_my))
  
  x <- as.matrix(training[, as.character(1:100)])
  y <- training$celltype_my %>% fct_relevel(c("non myeloid", "myeloid"))
  
  # Fit & Cross-validation (lasso deviance)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="binomial", type.measure="deviance", parallel=TRUE) 
  # here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
  ica_100_8test_glmnetCV_facs$Ia.lasso.dev[[i]] <- fit.cv
  
  png(filename = paste0("~/projects/icatest/output/test8/Ia.CVplot.facs.", i, ".png"), width = 800, height = 600, units = "px") #  pointsize = 13
  par(mfrow=c(2,2))
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="blue", paste("Ia", "lasso", "deviance", "facs", i))
  
  # Fit & Cross-validation (lasso class)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="binomial", type.measure="class", parallel=TRUE) 
  ica_100_8test_glmnetCV_facs$Ia.lasso.class[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="blue", paste("Ia", "lasso", "class", "facs", i))
  
  # Fit & Cross-validation (elnet deviance)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=0.9, family="binomial", type.measure="deviance", parallel=TRUE) 
  ica_100_8test_glmnetCV_facs$Ia.elnet.dev[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="blue", paste("Ia", "elnet", "deviance", "facs", i))
  
  # Fit & Cross-validation (elnet class)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=0.9, family="binomial", type.measure="class", parallel=TRUE) 
  ica_100_8test_glmnetCV_facs$Ia.elnet.class[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="blue", paste("Ia", "elnet", "class", "facs", i))
  dev.off()
  
}

#-----------------------------------------------------------------------------#
### Ib: multi-nomial, level B: myeloid cell types (celltype_my_type), .A mat values
#-----------------------------------------------------------------------------#

for (i in 1:4) {
  
  a <- ica_100_8test_A_facs[[i]]
  
  # global
  training <- a[ inTrain.tissue, ] %>% dplyr::filter(str_detect(celltype_my_type, "^non myeloid", negate = T)) %>% dplyr::mutate(celltype_my_type = fct_drop(celltype_my_type))
  testing  <- a[-inTrain.tissue, ] %>% dplyr::filter(str_detect(celltype_my_type, "^non myeloid", negate = T)) %>% dplyr::mutate(celltype_my_type = fct_drop(celltype_my_type))
  
  print(dplyr::count(training, celltype_my_type))
  print(dplyr::count(testing, celltype_my_type))
  
  x <- as.matrix(training[, as.character(1:100)]) 
  y <- training$celltype_my_type
  
  # Fit & Cross-validation (lasso deviance)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  # here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
  ica_100_8test_glmnetCV_facs$Ib.lasso.dev[[i]] <- fit.cv
  
  png(filename = paste0("~/projects/icatest/output/test8/Ib.CVplot.facs.", i, ".png"), width = 800, height = 600, units = "px") #  pointsize = 13
  par(mfrow=c(2,2))
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="blue", paste("Ib", "lasso", "deviance", "facs", i))
  
  # Fit & Cross-validation (lasso class)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="class", type.multinomial="ungrouped", parallel=TRUE) 
  ica_100_8test_glmnetCV_facs$Ib.lasso.class[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="blue", paste("Ib", "lasso", "class", "facs", i))
  
  # Fit & Cross-validation (elnet deviance)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=0.9, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  ica_100_8test_glmnetCV_facs$Ib.elnet.dev[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="blue", paste("Ib", "elnet", "deviance", "facs", i))
  
  # Fit & Cross-validation (elnet class)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=0.9, family="multinomial", type.measure="class", type.multinomial="ungrouped", parallel=TRUE) 
  ica_100_8test_glmnetCV_facs$Ib.elnet.class[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="blue", paste("Ib", "elnet", "class", "facs", i))
  dev.off()
  
}

#-----------------------------------------------------------------------------#
### Ic: multi-nomial, level C: macrophages accross tissues (celltype_my_type_tissue), .A mat values
#-----------------------------------------------------------------------------#

for (i in 1:4) {
  
  a <- ica_100_8test_A_facs[[i]]
  
  # global
  training <- a[ inTrain.tissue, ] %>% dplyr::filter(str_detect(celltype_my_type_tissue, "^macrophage -")|str_detect(celltype_my_type_tissue, "^Kupffer cell -")|str_detect(celltype_my_type_tissue, "^microglial cell -")) %>% dplyr::mutate(celltype_my_type_tissue = fct_drop(celltype_my_type_tissue))
  testing  <- a[-inTrain.tissue, ] %>% dplyr::filter(str_detect(celltype_my_type_tissue, "^macrophage -")|str_detect(celltype_my_type_tissue, "^Kupffer cell -")|str_detect(celltype_my_type_tissue, "^microglial cell -")) %>% dplyr::mutate(celltype_my_type_tissue = fct_drop(celltype_my_type_tissue))
  
  print(dplyr::count(training, celltype_my_type_tissue))
  print(dplyr::count(testing, celltype_my_type_tissue))
  
  x <- as.matrix(training[, as.character(1:100)]) 
  y <- training$celltype_my_type_tissue
  
  # Fit & Cross-validation (lasso deviance)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  # here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
  ica_100_8test_glmnetCV_facs$Ic.lasso.dev[[i]] <- fit.cv
  
  png(filename = paste0("~/projects/icatest/output/test8/Ic.CVplot.facs.", i, ".png"), width = 800, height = 600, units = "px") #  pointsize = 13
  par(mfrow=c(2,2))
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="blue", paste("Ic", "lasso", "deviance", "facs", i))
  
  # Fit & Cross-validation (lasso class)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="class", type.multinomial="ungrouped", parallel=TRUE) 
  ica_100_8test_glmnetCV_facs$Ic.lasso.class[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="blue", paste("Ic", "lasso", "class", "facs", i))
  
  # Fit & Cross-validation (elnet deviance)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=0.9, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  ica_100_8test_glmnetCV_facs$Ic.elnet.dev[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="blue", paste("Ic", "elnet", "deviance", "facs", i))
  
  # Fit & Cross-validation (elnet class)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=0.9, family="multinomial", type.measure="class", type.multinomial="ungrouped", parallel=TRUE) 
  ica_100_8test_glmnetCV_facs$Ic.elnet.class[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="blue", paste("Ic", "elnet", "class", "facs", i))
  dev.off()
  
}

#-----------------------------------------------------------------------------#
### Inc: multi-nomial, non-chierarchical (myeloid cell types and non myeloid) (celltype_my_type), .A mat values
#-----------------------------------------------------------------------------#

for (i in 1:4) {
  
  a <- ica_100_8test_A_facs[[i]]
  
  # global
  training <- a[ inTrain.tissue, ] 
  testing  <- a[-inTrain.tissue, ] 
  
  print(dplyr::count(training, celltype_my_type))
  print(dplyr::count(testing, celltype_my_type))
  
  x <- as.matrix(training[, as.character(1:100)]) 
  y <- training$celltype_my_type
  
  # Fit & Cross-validation (lasso deviance)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  # here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
  ica_100_8test_glmnetCV_facs$Inc.lasso.dev[[i]] <- fit.cv
  
  png(filename = paste0("~/projects/icatest/output/test8/Inc.CVplot.facs.", i, ".png"), width = 800, height = 600, units = "px") #  pointsize = 13
  par(mfrow=c(2,2))
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="blue", paste("Inc", "lasso", "deviance", "facs", i))
  
  # Fit & Cross-validation (lasso class)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="class", type.multinomial="ungrouped", parallel=TRUE) 
  ica_100_8test_glmnetCV_facs$Inc.lasso.class[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="blue", paste("Inc", "lasso", "class", "facs", i))
  
  # Fit & Cross-validation (elnet deviance)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=0.9, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  ica_100_8test_glmnetCV_facs$Inc.elnet.dev[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="blue", paste("Inc", "elnet", "deviance", "facs", i))
  
  # Fit & Cross-validation (elnet class)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=0.9, family="multinomial", type.measure="class", type.multinomial="ungrouped", parallel=TRUE) 
  ica_100_8test_glmnetCV_facs$Inc.elnet.class[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="blue", paste("Inc", "elnet", "class", "facs", i))
  dev.off()
  
}
stopCluster(cores)

saveRDS(ica_100_8test_glmnetCV_facs, file = "~/projects/icatest/output/test8/ica_100_8test_glmnetCV_facs.rds")

#-----------------------------------------------------------------------------#





############### Analyses and extracting stats from fit.cv objects #############

# get the right data split
inTrain.tissue <- readRDS("~/projects/icatest/output/test8/inTrain_facs.rds")

### Ia analysis ###
#-----------------------------------------------------------------------------#
outcomes <- c("myeloid")

# setup a list of cv and prediction analyses
Ia.list.str <- list(nbetas = matrix(0, nrow = 4, ncol = length(outcomes), dimnames = list(c("VSt.log", "VSt.rc", "SCt.res", "SCt.rc"), outcomes)), 
                    ave.sparsity = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.ber = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.kappa = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.bacc = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    bacc = matrix(0, nrow = 4, ncol = length(outcomes), dimnames = list(c("VSt.log", "VSt.rc", "SCt.res", "SCt.rc"), outcomes)), 
                    Betas = list(), 
                    Confusion = list())

ica_100_8test_glmnet_stat_facs <- list(Ia = list(lasso.dev = Ia.list.str, 
                                                 lasso.class = Ia.list.str, 
                                                 elnet.dev = Ia.list.str, 
                                                 elnet.class = Ia.list.str), 
                                       Ib = list(), 
                                       Ic = list(), 
                                       Inc = list())

for (i in 1:4) {
  
  a <- ica_100_8test_A_facs[[i]]
  
  # global
  training <- a[ inTrain.tissue, ] %>% dplyr::mutate(celltype_my = fct_relevel(celltype_my, c("non myeloid", "myeloid")))
  testing  <- a[-inTrain.tissue, ] %>% dplyr::mutate(celltype_my = fct_relevel(celltype_my, c("non myeloid", "myeloid")))
  
  print(dplyr::count(training, celltype_my))
  print(dplyr::count(testing, celltype_my))
  
  newx <- as.matrix(testing[, as.character(1:100)])
  
  for (mod in c("lasso.dev", "lasso.class", "elnet.dev", "elnet.class")) {
    
    I.mod <- paste0("Ia.", mod)
    fit.cv <- ica_100_8test_glmnetCV_facs[[I.mod]][[i]]
    betas <- as.matrix(coef(fit.cv, s = "lambda.1se"))
    colnames(betas) <- "myeloid"
    betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble() %>% dplyr::filter(myeloid != 0) # change!
    
    ica_100_8test_glmnet_stat_facs$Ia[[mod]]$Betas[[i]] <- betas
    ica_100_8test_glmnet_stat_facs$Ia[[mod]]$ave.sparsity[i] <- nrow(betas)-1
    ica_100_8test_glmnet_stat_facs$Ia[[mod]]$nbetas[i, 1] <- nrow(betas)-1 # change!
    
    #---------Predictions and Assesment-------------#
    
    # Predicted classes
    prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
    prediction <- factor(prediction, levels = levels(testing$celltype_my))
    
    # library("mixOmics") # alt. with more options: "measures" package
    confussion <- get.confusion_matrix(truth = testing$celltype_my, predicted = prediction)
    BER <- get.BER(confussion)
    ica_100_8test_glmnet_stat_facs$Ia[[mod]]$ave.ber[i] <- BER
    
    pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
    confussion <- as.data.frame(confussion) %>% rownames_to_column()
    colnames(confussion) <- c("(Truth)", pred.names)
    # print(confussion) # also kable it !
    ica_100_8test_glmnet_stat_facs$Ia[[mod]]$Confusion[[i]] <- confussion
    
    # library("caret")
    confussion.c <- confusionMatrix(data=prediction, reference=testing$celltype_my, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
    confussionByClass <- confussion.c$byClass[c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3) # change!
    
    ica_100_8test_glmnet_stat_facs$Ia[[mod]]$ave.kappa[i] <- confussion.c$overall[2]
    ica_100_8test_glmnet_stat_facs$Ia[[mod]]$ave.bacc[i] <- confussion.c$overall[1]
    ica_100_8test_glmnet_stat_facs$Ia[[mod]]$bacc[i, 1] <- confussionByClass[6] # change!
    
    print(paste0("Level: Ia", "; Type: ", mod, "; Data: ", i, " - nBeta: ", nrow(betas)-1, " BER: ", round(BER, 3), " aveAcc: ", round(confussion.c$overall[1], 3), " aveKappa: ", round(confussion.c$overall[2], 3)))
    
  }
  
} 


### Ib analysis ###
#-----------------------------------------------------------------------------#
outcomes <- c("granulocyte", "basophil", "monocyte", "macrophage", "Kupffer cell", "microglial cell", "brain pericyte")

# setup a list of cv and prediction analyses
Ib.list.str <- list(nbetas = matrix(0, nrow = 4, ncol = length(outcomes), dimnames = list(c("VSt.log", "VSt.rc", "SCt.res", "SCt.rc"), outcomes)), 
                    ave.sparsity = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.ber = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.kappa = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.bacc = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    bacc = matrix(0, nrow = 4, ncol = length(outcomes), dimnames = list(c("VSt.log", "VSt.rc", "SCt.res", "SCt.rc"), outcomes)), 
                    Betas = list(), 
                    Confusion = list())

ica_100_8test_glmnet_stat_facs$Ib <- list(lasso.dev = Ib.list.str, 
                                          lasso.class = Ib.list.str, 
                                          elnet.dev = Ib.list.str, 
                                          elnet.class = Ib.list.str)  

for (i in 1:4) {
  
  a <- ica_100_8test_A_facs[[i]]
  
  # global
  training <- a[ inTrain.tissue, ] %>% dplyr::filter(str_detect(celltype_my_type, "^non myeloid", negate = T)) %>% dplyr::mutate(celltype_my_type = fct_drop(celltype_my_type))
  testing  <- a[-inTrain.tissue, ] %>% dplyr::filter(str_detect(celltype_my_type, "^non myeloid", negate = T)) %>% dplyr::mutate(celltype_my_type = fct_drop(celltype_my_type))
  
  print(dplyr::count(training, celltype_my_type))
  print(dplyr::count(testing, celltype_my_type))
  
  newx <- as.matrix(testing[, as.character(1:100)])
  
  for (mod in c("lasso.dev", "lasso.class", "elnet.dev", "elnet.class")) {
    
    I.mod <- paste0("Ib.", mod)
    fit.cv <- ica_100_8test_glmnetCV_facs[[I.mod]][[i]]
    betas <- coef(fit.cv, s = "lambda.1se")
    outcomes <- names(betas)
    betas <- purrr::map(betas, .f = as.matrix)
    betas <- do.call(cbind, betas)
    colnames(betas) <- outcomes
    betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble() %>% dplyr::filter_at(vars(-compID), any_vars(.!=0)) 
    
    ica_100_8test_glmnet_stat_facs$Ib[[mod]]$Betas[[i]] <- betas
    ica_100_8test_glmnet_stat_facs$Ib[[mod]]$ave.sparsity[i] <- nrow(betas)-1
    ica_100_8test_glmnet_stat_facs$Ib[[mod]]$nbetas[i, ] <- betas %>% summarize_at(vars(-compID), .funs = function(x) sum(x != 0)-1) %>% slice(1) %>% unlist() # change!
    
    #---------Predictions and Assesment-------------#
    
    # Predicted classes
    prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
    prediction <- factor(prediction, levels = levels(testing$celltype_my_type))
    
    # library("mixOmics") # alt. with more options: "measures" package
    confussion <- get.confusion_matrix(truth = testing$celltype_my_type, predicted = prediction)
    BER <- get.BER(confussion)
    ica_100_8test_glmnet_stat_facs$Ib[[mod]]$ave.ber[i] <- BER
    
    pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
    confussion <- as.data.frame(confussion) %>% rownames_to_column()
    colnames(confussion) <- c("(Truth)", pred.names)
    # print(confussion) # also kable it !
    ica_100_8test_glmnet_stat_facs$Ib[[mod]]$Confusion[[i]] <- confussion
    
    # library("caret")
    confussion.c <- confusionMatrix(data=prediction, reference=testing$celltype_my_type, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
    confussionByClass <- confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3) # change!
    
    ica_100_8test_glmnet_stat_facs$Ib[[mod]]$ave.kappa[i] <- confussion.c$overall[2]
    ica_100_8test_glmnet_stat_facs$Ib[[mod]]$ave.bacc[i] <- confussion.c$overall[1]
    ica_100_8test_glmnet_stat_facs$Ib[[mod]]$bacc[i, ] <- confussionByClass[, 6] # change!
    
    print(paste0("Level: Ib", "; Type: ", mod, "; Data: ", i, " - nBeta: ", nrow(betas)-1, " BER: ", round(BER, 3), " aveAcc: ", round(confussion.c$overall[1], 3), " aveKappa: ", round(confussion.c$overall[2], 3)))
    
  }
  
}


### Ic analysis ###
#-----------------------------------------------------------------------------#
outcomes <- c("macrophage - Marrow", "macrophage - Spleen", "macrophage - Kidney", "macrophage - Limb_Muscle", "macrophage - Diaphragm", "macrophage - Brain_Myeloid", "Kupffer cell - Liver", "microglial cell - Brain_Myeloid")

# setup a list of cv and prediction analyses
Ic.list.str <- list(nbetas = matrix(0, nrow = 4, ncol = length(outcomes), dimnames = list(c("VSt.log", "VSt.rc", "SCt.res", "SCt.rc"), outcomes)), 
                    ave.sparsity = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.ber = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.kappa = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.bacc = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    bacc = matrix(0, nrow = 4, ncol = length(outcomes), dimnames = list(c("VSt.log", "VSt.rc", "SCt.res", "SCt.rc"), outcomes)), 
                    Betas = list(), 
                    Confusion = list())

ica_100_8test_glmnet_stat_facs$Ic <- list(lasso.dev = Ic.list.str, 
                                          lasso.class = Ic.list.str, 
                                          elnet.dev = Ic.list.str, 
                                          elnet.class = Ic.list.str)  

for (i in 1:4) {
  
  a <- ica_100_8test_A_facs[[i]]
  
  # global
  training <- a[ inTrain.tissue, ] %>% dplyr::filter(str_detect(celltype_my_type_tissue, "^macrophage -")|str_detect(celltype_my_type_tissue, "^Kupffer cell -")|str_detect(celltype_my_type_tissue, "^microglial cell -")) %>% dplyr::mutate(celltype_my_type_tissue = fct_drop(celltype_my_type_tissue))
  testing  <- a[-inTrain.tissue, ] %>% dplyr::filter(str_detect(celltype_my_type_tissue, "^macrophage -")|str_detect(celltype_my_type_tissue, "^Kupffer cell -")|str_detect(celltype_my_type_tissue, "^microglial cell -")) %>% dplyr::mutate(celltype_my_type_tissue = fct_drop(celltype_my_type_tissue))
  
  print(dplyr::count(training, celltype_my_type_tissue))
  print(dplyr::count(testing, celltype_my_type_tissue))
  
  newx <- as.matrix(testing[, as.character(1:100)])
  
  for (mod in c("lasso.dev", "lasso.class", "elnet.dev", "elnet.class")) {
    
    I.mod <- paste0("Ic.", mod)
    fit.cv <- ica_100_8test_glmnetCV_facs[[I.mod]][[i]]
    betas <- coef(fit.cv, s = "lambda.1se")
    outcomes <- names(betas)
    betas <- purrr::map(betas, .f = as.matrix)
    betas <- do.call(cbind, betas)
    colnames(betas) <- outcomes
    betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble() %>% dplyr::filter_at(vars(-compID), any_vars(.!=0)) 
    
    ica_100_8test_glmnet_stat_facs$Ic[[mod]]$Betas[[i]] <- betas
    ica_100_8test_glmnet_stat_facs$Ic[[mod]]$ave.sparsity[i] <- nrow(betas)-1
    ica_100_8test_glmnet_stat_facs$Ic[[mod]]$nbetas[i, ] <- betas %>% summarize_at(vars(-compID), .funs = function(x) sum(x != 0)-1) %>% slice(1) %>% unlist() # change!
    
    #---------Predictions and Assesment-------------#
    
    # Predicted classes
    prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
    prediction <- factor(prediction, levels = levels(testing$celltype_my_type_tissue))
    
    # library("mixOmics") # alt. with more options: "measures" package
    confussion <- get.confusion_matrix(truth = testing$celltype_my_type_tissue, predicted = prediction)
    BER <- get.BER(confussion)
    ica_100_8test_glmnet_stat_facs$Ic[[mod]]$ave.ber[i] <- BER
    
    pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
    confussion <- as.data.frame(confussion) %>% rownames_to_column()
    colnames(confussion) <- c("(Truth)", pred.names)
    # print(confussion) # also kable it !
    ica_100_8test_glmnet_stat_facs$Ic[[mod]]$Confusion[[i]] <- confussion
    
    # library("caret")
    confussion.c <- confusionMatrix(data=prediction, reference=testing$celltype_my_type_tissue, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
    confussionByClass <- confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3) # change!
    
    ica_100_8test_glmnet_stat_facs$Ic[[mod]]$ave.kappa[i] <- confussion.c$overall[2]
    ica_100_8test_glmnet_stat_facs$Ic[[mod]]$ave.bacc[i] <- confussion.c$overall[1]
    ica_100_8test_glmnet_stat_facs$Ic[[mod]]$bacc[i, ] <- confussionByClass[, 6] # change!
    
    print(paste0("Level: Ic", "; Type: ", mod, "; Data: ", i, " - nBeta: ", nrow(betas)-1, " BER: ", round(BER, 3), " aveAcc: ", round(confussion.c$overall[1], 3), " aveKappa: ", round(confussion.c$overall[2], 3)))
    
  }
  
}


### Inc analysis ###
#-----------------------------------------------------------------------------#
outcomes <- c("granulocyte", "basophil", "monocyte", "macrophage", "Kupffer cell", "microglial cell", "brain pericyte", "non myeloid")

# setup a list of cv and prediction analyses
Inc.list.str <- list(nbetas = matrix(0, nrow = 4, ncol = length(outcomes), dimnames = list(c("VSt.log", "VSt.rc", "SCt.res", "SCt.rc"), outcomes)), 
                     ave.sparsity = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                     ave.ber = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                     ave.kappa = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                     ave.bacc = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                     bacc = matrix(0, nrow = 4, ncol = length(outcomes), dimnames = list(c("VSt.log", "VSt.rc", "SCt.res", "SCt.rc"), outcomes)), 
                     Betas = list(), 
                     Confusion = list())

ica_100_8test_glmnet_stat_facs$Inc <- list(lasso.dev = Inc.list.str, 
                                           lasso.class = Inc.list.str, 
                                           elnet.dev = Inc.list.str, 
                                           elnet.class = Inc.list.str)  

for (i in 1:4) {
  
  a <- ica_100_8test_A_facs[[i]]
  
  # global
  training <- a[ inTrain.tissue, ] 
  testing  <- a[-inTrain.tissue, ] 
  
  print(dplyr::count(training, celltype_my_type))
  print(dplyr::count(testing, celltype_my_type))
  
  newx <- as.matrix(testing[, as.character(1:100)])
  
  for (mod in c("lasso.dev", "lasso.class", "elnet.dev", "elnet.class")) {
    
    I.mod <- paste0("Inc.", mod)
    fit.cv <- ica_100_8test_glmnetCV_facs[[I.mod]][[i]]
    betas <- coef(fit.cv, s = "lambda.1se")
    outcomes <- names(betas)
    betas <- purrr::map(betas, .f = as.matrix)
    betas <- do.call(cbind, betas)
    colnames(betas) <- outcomes
    betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble() %>% dplyr::filter_at(vars(-compID), any_vars(.!=0)) 
    
    ica_100_8test_glmnet_stat_facs$Inc[[mod]]$Betas[[i]] <- betas
    ica_100_8test_glmnet_stat_facs$Inc[[mod]]$ave.sparsity[i] <- nrow(betas)-1
    ica_100_8test_glmnet_stat_facs$Inc[[mod]]$nbetas[i, ] <- betas %>% summarize_at(vars(-compID), .funs = function(x) sum(x != 0)-1) %>% slice(1) %>% unlist() # change!
    
    #---------Predictions and Assesment-------------#
    
    # Predicted classes
    prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
    prediction <- factor(prediction, levels = levels(testing$celltype_my_type))
    
    # library("mixOmics") # alt. with more options: "measures" package
    confussion <- get.confusion_matrix(truth = testing$celltype_my_type, predicted = prediction)
    BER <- get.BER(confussion)
    ica_100_8test_glmnet_stat_facs$Inc[[mod]]$ave.ber[i] <- BER
    
    pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
    confussion <- as.data.frame(confussion) %>% rownames_to_column()
    colnames(confussion) <- c("(Truth)", pred.names)
    # print(confussion) # also kable it !
    ica_100_8test_glmnet_stat_facs$Inc[[mod]]$Confusion[[i]] <- confussion
    
    # library("caret")
    confussion.c <- confusionMatrix(data=prediction, reference=testing$celltype_my_type, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
    confussionByClass <- confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3) # change!
    
    ica_100_8test_glmnet_stat_facs$Inc[[mod]]$ave.kappa[i] <- confussion.c$overall[2]
    ica_100_8test_glmnet_stat_facs$Inc[[mod]]$ave.bacc[i] <- confussion.c$overall[1]
    ica_100_8test_glmnet_stat_facs$Inc[[mod]]$bacc[i, ] <- confussionByClass[, 6] # change!
    
    print(paste0("Level: Inc", "; Type: ", mod, "; Data: ", i, " - nBeta: ", nrow(betas)-1, " BER: ", round(BER, 3), " aveAcc: ", round(confussion.c$overall[1], 3), " aveKappa: ", round(confussion.c$overall[2], 3)))
    
  }
  
}

saveRDS(ica_100_8test_glmnet_stat_facs, file = "~/projects/icatest/output/test8/ica_100_8test_glmnet_stat_facs.rds")





############################################################################
### DROP lasso ----------------------------------------------------------###
############################################################################


#---------Split to Training and Testing datasets----------#
library(caret)

# These cell indices are same for Gen.Exp. Df and all A matrices
# Best to make them based on use_celltype_tissue variable (to keep celltype-tissue representation)!
ica_100_8test_A_drop <- readRDS("~/projects/icatest/output/test8/ica_100_8test_A_drop.rds")
set.seed(707)
inTrain.tissue <- createDataPartition(y = ica_100_8test_A_drop[[1]]$celltype_my_type_tissue, p = .85, list = FALSE) 
saveRDS(inTrain.tissue, file = "~/projects/icatest/output/test8/inTrain_drop.rds")
#---------------------------------------------------------#


################### Create loops for just saving fit.cv #######################

# get the right data split
inTrain.tissue <- readRDS("~/projects/icatest/output/test8/inTrain_drop.rds")

# setup a list of cv results
ica_100_8test_glmnetCV_drop <- list() # each with CVfit object for four data versions (1:"facs.vst.log" 2:"facs.vst.rc"  3:"facs.sct.res" 4:"facs.sct.cnt")

#-----------------------------------------------------------------------------#
### Ia: bi-nomial, level A: myeloid / non myeloid (celltype_my), .A mat values
#-----------------------------------------------------------------------------#

cores <- makeCluster(detectCores(logical = FALSE), type='PSOCK')
registerDoParallel(cores)

for (i in 1:4) {
  
  a <- ica_100_8test_A_drop[[i]]
  
  # global
  training <- a[ inTrain.tissue, ]
  testing  <- a[-inTrain.tissue, ]
  
  print(dplyr::count(training, celltype_my))
  print(dplyr::count(testing, celltype_my))
  
  x <- as.matrix(training[, as.character(1:100)])
  y <- training$celltype_my %>% fct_relevel(c("non myeloid", "myeloid"))
  
  # Fit & Cross-validation (lasso deviance)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="binomial", type.measure="deviance", parallel=TRUE) 
  # here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
  ica_100_8test_glmnetCV_drop$Ia.lasso.dev[[i]] <- fit.cv
  
  png(filename = paste0("~/projects/icatest/output/test8/Ia.CVplot.drop.", i, ".png"), width = 800, height = 600, units = "px") #  pointsize = 13
  par(mfrow=c(2,2))
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="seagreen", paste("Ia", "lasso", "deviance", "drop", i))
  
  # Fit & Cross-validation (lasso class)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="binomial", type.measure="class", parallel=TRUE) 
  ica_100_8test_glmnetCV_drop$Ia.lasso.class[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="seagreen", paste("Ia", "lasso", "class", "drop", i))
  
  # Fit & Cross-validation (elnet deviance)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=0.9, family="binomial", type.measure="deviance", parallel=TRUE) 
  ica_100_8test_glmnetCV_drop$Ia.elnet.dev[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="seagreen", paste("Ia", "elnet", "deviance", "drop", i))
  
  # Fit & Cross-validation (elnet class)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=0.9, family="binomial", type.measure="class", parallel=TRUE) 
  ica_100_8test_glmnetCV_drop$Ia.elnet.class[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="seagreen", paste("Ia", "elnet", "class", "drop", i))
  dev.off()
  
}

#-----------------------------------------------------------------------------#
### Ib: multi-nomial, level B: myeloid cell types (celltype_my_type), .A mat values
#-----------------------------------------------------------------------------#

for (i in 1:4) {
  
  a <- ica_100_8test_A_drop[[i]]
  
  # global
  training <- a[ inTrain.tissue, ] %>% dplyr::filter(str_detect(celltype_my_type, "^non myeloid", negate = T)) %>% dplyr::mutate(celltype_my_type = fct_drop(celltype_my_type))
  testing  <- a[-inTrain.tissue, ] %>% dplyr::filter(str_detect(celltype_my_type, "^non myeloid", negate = T)) %>% dplyr::mutate(celltype_my_type = fct_drop(celltype_my_type))
  
  print(dplyr::count(training, celltype_my_type))
  print(dplyr::count(testing, celltype_my_type))
  
  x <- as.matrix(training[, as.character(1:100)]) 
  y <- training$celltype_my_type
  
  # Fit & Cross-validation (lasso deviance)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  # here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
  ica_100_8test_glmnetCV_drop$Ib.lasso.dev[[i]] <- fit.cv
  
  png(filename = paste0("~/projects/icatest/output/test8/Ib.CVplot.drop.", i, ".png"), width = 800, height = 600, units = "px") #  pointsize = 13
  par(mfrow=c(2,2))
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="seagreen", paste("Ib", "lasso", "deviance", "drop", i))
  
  # Fit & Cross-validation (lasso class)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="class", type.multinomial="ungrouped", parallel=TRUE) 
  ica_100_8test_glmnetCV_drop$Ib.lasso.class[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="seagreen", paste("Ib", "lasso", "class", "drop", i))
  
  # Fit & Cross-validation (elnet deviance)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=0.9, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  ica_100_8test_glmnetCV_drop$Ib.elnet.dev[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="seagreen", paste("Ib", "elnet", "deviance", "drop", i))
  
  # Fit & Cross-validation (elnet class)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=0.9, family="multinomial", type.measure="class", type.multinomial="ungrouped", parallel=TRUE) 
  ica_100_8test_glmnetCV_drop$Ib.elnet.class[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="seagreen", paste("Ib", "elnet", "class", "drop", i))
  dev.off()
  
}

#-----------------------------------------------------------------------------#
### Ic: multi-nomial, level C: macrophages accross tissues (celltype_my_type_tissue), .A mat values
#-----------------------------------------------------------------------------#

for (i in 1:4) {
  
  a <- ica_100_8test_A_drop[[i]]
  
  # global
  training <- a[ inTrain.tissue, ] %>% dplyr::filter(str_detect(celltype_my_type_tissue, "^macrophage -")|str_detect(celltype_my_type_tissue, "^Kupffer cell -")|str_detect(celltype_my_type_tissue, "^microglial cell -")) %>% dplyr::mutate(celltype_my_type_tissue = fct_drop(celltype_my_type_tissue))
  testing  <- a[-inTrain.tissue, ] %>% dplyr::filter(str_detect(celltype_my_type_tissue, "^macrophage -")|str_detect(celltype_my_type_tissue, "^Kupffer cell -")|str_detect(celltype_my_type_tissue, "^microglial cell -")) %>% dplyr::mutate(celltype_my_type_tissue = fct_drop(celltype_my_type_tissue))
  
  print(dplyr::count(training, celltype_my_type_tissue))
  print(dplyr::count(testing, celltype_my_type_tissue))
  
  x <- as.matrix(training[, as.character(1:100)]) 
  y <- training$celltype_my_type_tissue
  
  # Fit & Cross-validation (lasso deviance)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  # here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
  ica_100_8test_glmnetCV_drop$Ic.lasso.dev[[i]] <- fit.cv
  
  png(filename = paste0("~/projects/icatest/output/test8/Ic.CVplot.drop.", i, ".png"), width = 800, height = 600, units = "px") #  pointsize = 13
  par(mfrow=c(2,2))
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="seagreen", paste("Ic", "lasso", "deviance", "drop", i))
  
  # Fit & Cross-validation (lasso class)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="class", type.multinomial="ungrouped", parallel=TRUE) 
  ica_100_8test_glmnetCV_drop$Ic.lasso.class[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="seagreen", paste("Ic", "lasso", "class", "drop", i))
  
  # Fit & Cross-validation (elnet deviance)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=0.9, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  ica_100_8test_glmnetCV_drop$Ic.elnet.dev[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="seagreen", paste("Ic", "elnet", "deviance", "drop", i))
  
  # Fit & Cross-validation (elnet class)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=0.9, family="multinomial", type.measure="class", type.multinomial="ungrouped", parallel=TRUE) 
  ica_100_8test_glmnetCV_drop$Ic.elnet.class[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="seagreen", paste("Ic", "elnet", "class", "drop", i))
  dev.off()
  
}

#-----------------------------------------------------------------------------#
### Inc: multi-nomial, non-chierarchical (myeloid cell types and non myeloid) (celltype_my_type), .A mat values
#-----------------------------------------------------------------------------#

cores <- makeCluster(detectCores(logical = FALSE), type='PSOCK')
registerDoParallel(cores)

for (i in 1:4) {
  
  a <- ica_100_8test_A_drop[[i]]
  
  # global
  training <- a[ inTrain.tissue, ] 
  testing  <- a[-inTrain.tissue, ] 
  
  print(dplyr::count(training, celltype_my_type))
  print(dplyr::count(testing, celltype_my_type))
  
  x <- as.matrix(training[, as.character(1:100)]) 
  y <- training$celltype_my_type
  
  # Fit & Cross-validation (lasso deviance)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  # here: can play with type.measure = "deviance" (default), = "class", = "auc" (binomial only)
  ica_100_8test_glmnetCV_drop$Inc.lasso.dev[[i]] <- fit.cv
  
  png(filename = paste0("~/projects/icatest/output/test8/Inc.CVplot.drop.", i, ".png"), width = 800, height = 600, units = "px") #  pointsize = 13
  par(mfrow=c(2,2))
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="seagreen", paste("Inc", "lasso", "deviance", "drop", i))
  
  # Fit & Cross-validation (lasso class)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=1, family="multinomial", type.measure="class", type.multinomial="ungrouped", parallel=TRUE) 
  ica_100_8test_glmnetCV_drop$Inc.lasso.class[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="seagreen", paste("Inc", "lasso", "class", "drop", i))
  
  # Fit & Cross-validation (elnet deviance)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=0.9, family="multinomial", type.measure="deviance", type.multinomial="ungrouped", parallel=TRUE) 
  ica_100_8test_glmnetCV_drop$Inc.elnet.dev[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="seagreen", paste("Inc", "elnet", "deviance", "drop", i))
  
  # Fit & Cross-validation (elnet class)
  fit.cv <- cv.glmnet(x=x, y=y, alpha=0.9, family="multinomial", type.measure="class", type.multinomial="ungrouped", parallel=TRUE) 
  ica_100_8test_glmnetCV_drop$Inc.elnet.class[[i]] <- fit.cv
  plot(fit.cv); mtext(side=3, line=2.25, cex=1.25, col="seagreen", paste("Inc", "elnet", "class", "drop", i))
  dev.off()
  
}
stopCluster(cores)

saveRDS(ica_100_8test_glmnetCV_drop, file = "~/projects/icatest/output/test8/ica_100_8test_glmnetCV_drop.rds")

#-----------------------------------------------------------------------------#



############### Analyses and extracting stats from fit.cv objects #############

# get the right data split
inTrain.tissue <- readRDS("~/projects/icatest/output/test8/inTrain_drop.rds")

### Ia analysis ###
#-----------------------------------------------------------------------------#
outcomes <- c("myeloid")

# setup a list of cv and prediction analyses
Ia.list.str <- list(nbetas = matrix(0, nrow = 4, ncol = length(outcomes), dimnames = list(c("VSt.log", "VSt.rc", "SCt.res", "SCt.rc"), outcomes)), 
                    ave.sparsity = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.ber = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.kappa = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.bacc = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    bacc = matrix(0, nrow = 4, ncol = length(outcomes), dimnames = list(c("VSt.log", "VSt.rc", "SCt.res", "SCt.rc"), outcomes)), 
                    Betas = list(), 
                    Confusion = list())

ica_100_8test_glmnet_stat_drop <- list(Ia = list(lasso.dev = Ia.list.str, 
                                                 lasso.class = Ia.list.str, 
                                                 elnet.dev = Ia.list.str, 
                                                 elnet.class = Ia.list.str), 
                                       Ib = list(), 
                                       Ic = list(), 
                                       Inc = list())

for (i in 1:4) {
  
  a <- ica_100_8test_A_drop[[i]]
  
  # global
  training <- a[ inTrain.tissue, ] %>% dplyr::mutate(celltype_my = fct_relevel(celltype_my, c("non myeloid", "myeloid")))
  testing  <- a[-inTrain.tissue, ] %>% dplyr::mutate(celltype_my = fct_relevel(celltype_my, c("non myeloid", "myeloid")))
  
  print(dplyr::count(training, celltype_my))
  print(dplyr::count(testing, celltype_my))
  
  newx <- as.matrix(testing[, as.character(1:100)])
  
  for (mod in c("lasso.dev", "lasso.class", "elnet.dev", "elnet.class")) {
    
    I.mod <- paste0("Ia.", mod)
    fit.cv <- ica_100_8test_glmnetCV_drop[[I.mod]][[i]]
    betas <- as.matrix(coef(fit.cv, s = "lambda.1se"))
    colnames(betas) <- "myeloid"
    betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble() %>% dplyr::filter(myeloid != 0) # change!
    
    ica_100_8test_glmnet_stat_drop$Ia[[mod]]$Betas[[i]] <- betas
    ica_100_8test_glmnet_stat_drop$Ia[[mod]]$ave.sparsity[i] <- nrow(betas)-1
    ica_100_8test_glmnet_stat_drop$Ia[[mod]]$nbetas[i, 1] <- nrow(betas)-1 # change!
    
    #---------Predictions and Assesment-------------#
    
    # Predicted classes
    prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
    prediction <- factor(prediction, levels = levels(testing$celltype_my))
    
    # library("mixOmics") # alt. with more options: "measures" package
    confussion <- get.confusion_matrix(truth = testing$celltype_my, predicted = prediction)
    BER <- get.BER(confussion)
    ica_100_8test_glmnet_stat_drop$Ia[[mod]]$ave.ber[i] <- BER
    
    pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
    confussion <- as.data.frame(confussion) %>% rownames_to_column()
    colnames(confussion) <- c("(Truth)", pred.names)
    # print(confussion) # also kable it !
    ica_100_8test_glmnet_stat_drop$Ia[[mod]]$Confusion[[i]] <- confussion
    
    # library("caret")
    confussion.c <- confusionMatrix(data=prediction, reference=testing$celltype_my, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
    confussionByClass <- confussion.c$byClass[c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3) # change!
    
    ica_100_8test_glmnet_stat_drop$Ia[[mod]]$ave.kappa[i] <- confussion.c$overall[2]
    ica_100_8test_glmnet_stat_drop$Ia[[mod]]$ave.bacc[i] <- confussion.c$overall[1]
    ica_100_8test_glmnet_stat_drop$Ia[[mod]]$bacc[i, 1] <- confussionByClass[6] # change!
    
    print(paste0("Level: Ia", "; Type: ", mod, "; Data: ", i, " - nBeta: ", nrow(betas)-1, " BER: ", round(BER, 3), " aveAcc: ", round(confussion.c$overall[1], 3), " aveKappa: ", round(confussion.c$overall[2], 3)))
    
  }
  
} 


### Ib analysis ###
#-----------------------------------------------------------------------------#
outcomes <- c("granulocyte", "basophil", "promonocyte", "monocyte", "macrophage", "dendritic cell", "Langerhans cell", "mast cell")

# setup a list of cv and prediction analyses
Ib.list.str <- list(nbetas = matrix(0, nrow = 4, ncol = length(outcomes), dimnames = list(c("VSt.log", "VSt.rc", "SCt.res", "SCt.rc"), outcomes)), 
                    ave.sparsity = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.ber = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.kappa = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.bacc = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    bacc = matrix(0, nrow = 4, ncol = length(outcomes), dimnames = list(c("VSt.log", "VSt.rc", "SCt.res", "SCt.rc"), outcomes)), 
                    Betas = list(), 
                    Confusion = list())

ica_100_8test_glmnet_stat_drop$Ib <- list(lasso.dev = Ib.list.str, 
                                          lasso.class = Ib.list.str, 
                                          elnet.dev = Ib.list.str, 
                                          elnet.class = Ib.list.str)  

for (i in 1:4) {
  
  a <- ica_100_8test_A_drop[[i]]
  
  # global
  training <- a[ inTrain.tissue, ] %>% dplyr::filter(str_detect(celltype_my_type, "^non myeloid", negate = T)) %>% dplyr::mutate(celltype_my_type = fct_drop(celltype_my_type))
  testing  <- a[-inTrain.tissue, ] %>% dplyr::filter(str_detect(celltype_my_type, "^non myeloid", negate = T)) %>% dplyr::mutate(celltype_my_type = fct_drop(celltype_my_type))
  
  print(dplyr::count(training, celltype_my_type))
  print(dplyr::count(testing, celltype_my_type))
  
  newx <- as.matrix(testing[, as.character(1:100)])
  
  for (mod in c("lasso.dev", "lasso.class", "elnet.dev", "elnet.class")) {
    
    I.mod <- paste0("Ib.", mod)
    fit.cv <- ica_100_8test_glmnetCV_drop[[I.mod]][[i]]
    betas <- coef(fit.cv, s = "lambda.1se")
    outcomes <- names(betas)
    betas <- purrr::map(betas, .f = as.matrix)
    betas <- do.call(cbind, betas)
    colnames(betas) <- outcomes
    betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble() %>% dplyr::filter_at(vars(-compID), any_vars(.!=0)) 
    
    ica_100_8test_glmnet_stat_drop$Ib[[mod]]$Betas[[i]] <- betas
    ica_100_8test_glmnet_stat_drop$Ib[[mod]]$ave.sparsity[i] <- nrow(betas)-1
    ica_100_8test_glmnet_stat_drop$Ib[[mod]]$nbetas[i, ] <- betas %>% summarize_at(vars(-compID), .funs = function(x) sum(x != 0)-1) %>% slice(1) %>% unlist() # change!
    
    #---------Predictions and Assesment-------------#
    
    # Predicted classes
    prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
    prediction <- factor(prediction, levels = levels(testing$celltype_my_type))
    
    # library("mixOmics") # alt. with more options: "measures" package
    confussion <- get.confusion_matrix(truth = testing$celltype_my_type, predicted = prediction)
    BER <- get.BER(confussion)
    ica_100_8test_glmnet_stat_drop$Ib[[mod]]$ave.ber[i] <- BER
    
    pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
    confussion <- as.data.frame(confussion) %>% rownames_to_column()
    colnames(confussion) <- c("(Truth)", pred.names)
    # print(confussion) # also kable it !
    ica_100_8test_glmnet_stat_drop$Ib[[mod]]$Confusion[[i]] <- confussion
    
    # library("caret")
    confussion.c <- confusionMatrix(data=prediction, reference=testing$celltype_my_type, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
    confussionByClass <- confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3) # change!
    
    ica_100_8test_glmnet_stat_drop$Ib[[mod]]$ave.kappa[i] <- confussion.c$overall[2]
    ica_100_8test_glmnet_stat_drop$Ib[[mod]]$ave.bacc[i] <- confussion.c$overall[1]
    ica_100_8test_glmnet_stat_drop$Ib[[mod]]$bacc[i, ] <- confussionByClass[, 6] # change!
    
    print(paste0("Level: Ib", "; Type: ", mod, "; Data: ", i, " - nBeta: ", nrow(betas)-1, " BER: ", round(BER, 3), " aveAcc: ", round(confussion.c$overall[1], 3), " aveKappa: ", round(confussion.c$overall[2], 3)))
    
  }
  
}


### Ic analysis ###
#-----------------------------------------------------------------------------#
outcomes <- c("macrophage - Marrow", "macrophage - Spleen", "macrophage - Lung", "macrophage - Kidney", "macrophage - Mammary_Gland", "macrophage - Limb_Muscle")

# setup a list of cv and prediction analyses
Ic.list.str <- list(nbetas = matrix(0, nrow = 4, ncol = length(outcomes), dimnames = list(c("VSt.log", "VSt.rc", "SCt.res", "SCt.rc"), outcomes)), 
                    ave.sparsity = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.ber = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.kappa = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.bacc = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    bacc = matrix(0, nrow = 4, ncol = length(outcomes), dimnames = list(c("VSt.log", "VSt.rc", "SCt.res", "SCt.rc"), outcomes)), 
                    Betas = list(), 
                    Confusion = list())

ica_100_8test_glmnet_stat_drop$Ic <- list(lasso.dev = Ic.list.str, 
                                          lasso.class = Ic.list.str, 
                                          elnet.dev = Ic.list.str, 
                                          elnet.class = Ic.list.str)  

for (i in 1:4) {
  
  a <- ica_100_8test_A_drop[[i]]
  
  # global
  training <- a[ inTrain.tissue, ] %>% dplyr::filter(str_detect(celltype_my_type_tissue, "^macrophage -")) %>% dplyr::mutate(celltype_my_type_tissue = fct_drop(celltype_my_type_tissue))
  testing  <- a[-inTrain.tissue, ] %>% dplyr::filter(str_detect(celltype_my_type_tissue, "^macrophage -")) %>% dplyr::mutate(celltype_my_type_tissue = fct_drop(celltype_my_type_tissue))
  
  print(dplyr::count(training, celltype_my_type_tissue))
  print(dplyr::count(testing, celltype_my_type_tissue))
  
  newx <- as.matrix(testing[, as.character(1:100)])
  
  for (mod in c("lasso.dev", "lasso.class", "elnet.dev", "elnet.class")) {
    
    I.mod <- paste0("Ic.", mod)
    fit.cv <- ica_100_8test_glmnetCV_drop[[I.mod]][[i]]
    betas <- coef(fit.cv, s = "lambda.1se")
    outcomes <- names(betas)
    betas <- purrr::map(betas, .f = as.matrix)
    betas <- do.call(cbind, betas)
    colnames(betas) <- outcomes
    betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble() %>% dplyr::filter_at(vars(-compID), any_vars(.!=0)) 
    
    ica_100_8test_glmnet_stat_drop$Ic[[mod]]$Betas[[i]] <- betas
    ica_100_8test_glmnet_stat_drop$Ic[[mod]]$ave.sparsity[i] <- nrow(betas)-1
    ica_100_8test_glmnet_stat_drop$Ic[[mod]]$nbetas[i, ] <- betas %>% summarize_at(vars(-compID), .funs = function(x) sum(x != 0)-1) %>% slice(1) %>% unlist() # change!
    
    #---------Predictions and Assesment-------------#
    
    # Predicted classes
    prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
    prediction <- factor(prediction, levels = levels(testing$celltype_my_type_tissue))
    
    # library("mixOmics") # alt. with more options: "measures" package
    confussion <- get.confusion_matrix(truth = testing$celltype_my_type_tissue, predicted = prediction)
    BER <- get.BER(confussion)
    ica_100_8test_glmnet_stat_drop$Ic[[mod]]$ave.ber[i] <- BER
    
    pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
    confussion <- as.data.frame(confussion) %>% rownames_to_column()
    colnames(confussion) <- c("(Truth)", pred.names)
    # print(confussion) # also kable it !
    ica_100_8test_glmnet_stat_drop$Ic[[mod]]$Confusion[[i]] <- confussion
    
    # library("caret")
    confussion.c <- confusionMatrix(data=prediction, reference=testing$celltype_my_type_tissue, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
    confussionByClass <- confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3) # change!
    
    ica_100_8test_glmnet_stat_drop$Ic[[mod]]$ave.kappa[i] <- confussion.c$overall[2]
    ica_100_8test_glmnet_stat_drop$Ic[[mod]]$ave.bacc[i] <- confussion.c$overall[1]
    ica_100_8test_glmnet_stat_drop$Ic[[mod]]$bacc[i, ] <- confussionByClass[, 6] # change!
    
    print(paste0("Level: Ic", "; Type: ", mod, "; Data: ", i, " - nBeta: ", nrow(betas)-1, " BER: ", round(BER, 3), " aveAcc: ", round(confussion.c$overall[1], 3), " aveKappa: ", round(confussion.c$overall[2], 3)))
    
  }
  
}


### Inc analysis ###
#-----------------------------------------------------------------------------#
outcomes <- c("granulocyte", "basophil", "promonocyte", "monocyte", "macrophage", "dendritic cell", "Langerhans cell", "mast cell", "non myeloid")

# setup a list of cv and prediction analyses
Inc.list.str <- list(nbetas = matrix(0, nrow = 4, ncol = length(outcomes), dimnames = list(c("VSt.log", "VSt.rc", "SCt.res", "SCt.rc"), outcomes)), 
                    ave.sparsity = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.ber = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.kappa = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    ave.bacc = c("VSt.log"=0, "VSt.rc"=0, "SCt.res"=0, "SCt.rc"=0), 
                    bacc = matrix(0, nrow = 4, ncol = length(outcomes), dimnames = list(c("VSt.log", "VSt.rc", "SCt.res", "SCt.rc"), outcomes)), 
                    Betas = list(), 
                    Confusion = list())

ica_100_8test_glmnet_stat_drop$Inc <- list(lasso.dev = Inc.list.str, 
                                          lasso.class = Inc.list.str, 
                                          elnet.dev = Inc.list.str, 
                                          elnet.class = Inc.list.str)  

for (i in 1:4) {
  
  a <- ica_100_8test_A_drop[[i]]
  
  # global
  training <- a[ inTrain.tissue, ] 
  testing  <- a[-inTrain.tissue, ] 
  
  print(dplyr::count(training, celltype_my_type))
  print(dplyr::count(testing, celltype_my_type))
  
  newx <- as.matrix(testing[, as.character(1:100)])
  
  for (mod in c("lasso.dev", "lasso.class", "elnet.dev", "elnet.class")) {
    
    I.mod <- paste0("Inc.", mod)
    fit.cv <- ica_100_8test_glmnetCV_drop[[I.mod]][[i]]
    betas <- coef(fit.cv, s = "lambda.1se")
    outcomes <- names(betas)
    betas <- purrr::map(betas, .f = as.matrix)
    betas <- do.call(cbind, betas)
    colnames(betas) <- outcomes
    betas <- as.data.frame(betas) %>% rownames_to_column(var = "compID") %>% as_tibble() %>% dplyr::filter_at(vars(-compID), any_vars(.!=0)) 
    
    ica_100_8test_glmnet_stat_drop$Inc[[mod]]$Betas[[i]] <- betas
    ica_100_8test_glmnet_stat_drop$Inc[[mod]]$ave.sparsity[i] <- nrow(betas)-1
    ica_100_8test_glmnet_stat_drop$Inc[[mod]]$nbetas[i, ] <- betas %>% summarize_at(vars(-compID), .funs = function(x) sum(x != 0)-1) %>% slice(1) %>% unlist() # change!
    
    #---------Predictions and Assesment-------------#
    
    # Predicted classes
    prediction <- predict(fit.cv, newx = newx, s = "lambda.1se", type = "class") # a vector
    prediction <- factor(prediction, levels = levels(testing$celltype_my_type))
    
    # library("mixOmics") # alt. with more options: "measures" package
    confussion <- get.confusion_matrix(truth = testing$celltype_my_type, predicted = prediction)
    BER <- get.BER(confussion)
    ica_100_8test_glmnet_stat_drop$Inc[[mod]]$ave.ber[i] <- BER
    
    pred.names <- colnames(confussion) %>% str_replace_all("predicted.as", "p")
    confussion <- as.data.frame(confussion) %>% rownames_to_column()
    colnames(confussion) <- c("(Truth)", pred.names)
    # print(confussion) # also kable it !
    ica_100_8test_glmnet_stat_drop$Inc[[mod]]$Confusion[[i]] <- confussion
    
    # library("caret")
    confussion.c <- confusionMatrix(data=prediction, reference=testing$celltype_my_type, dnn = c("Prediction", "Reference"), mode = "sens_spec") 
    confussionByClass <- confussion.c$byClass[, c("Sensitivity", "Specificity", "Precision", "Recall", "F1", "Balanced Accuracy")] %>% round(digits = 3) # change!
    
    ica_100_8test_glmnet_stat_drop$Inc[[mod]]$ave.kappa[i] <- confussion.c$overall[2]
    ica_100_8test_glmnet_stat_drop$Inc[[mod]]$ave.bacc[i] <- confussion.c$overall[1]
    ica_100_8test_glmnet_stat_drop$Inc[[mod]]$bacc[i, ] <- confussionByClass[, 6] # change!
    
    print(paste0("Level: Inc", "; Type: ", mod, "; Data: ", i, " - nBeta: ", nrow(betas)-1, " BER: ", round(BER, 3), " aveAcc: ", round(confussion.c$overall[1], 3), " aveKappa: ", round(confussion.c$overall[2], 3)))
    
  }
  
}

saveRDS(ica_100_8test_glmnet_stat_drop, file = "~/projects/icatest/output/test8/ica_100_8test_glmnet_stat_drop.rds")



# FACS -------------------------------------------------------------------------------------------#

# 1 non myeloid  3595
# 2 myeloid       926
# [1] "Level: Ia; Type: lasso.dev; Data: 1 - nBeta: 66 BER: 0.026 aveAcc: 0.982 aveKappa: 0.946"
# [1] "Level: Ia; Type: lasso.class; Data: 1 - nBeta: 45 BER: 0.028 aveAcc: 0.981 aveKappa: 0.942"
# [1] "Level: Ia; Type: elnet.dev; Data: 1 - nBeta: 62 BER: 0.025 aveAcc: 0.983 aveKappa: 0.947"
# [1] "Level: Ia; Type: elnet.class; Data: 1 - nBeta: 52 BER: 0.026 aveAcc: 0.982 aveKappa: 0.945"
# 
# [1] "Level: Ia; Type: lasso.dev; Data: 2 - nBeta: 85 BER: 0.05 aveAcc: 0.971 aveKappa: 0.909"
# [1] "Level: Ia; Type: lasso.class; Data: 2 - nBeta: 85 BER: 0.05 aveAcc: 0.971 aveKappa: 0.909"
# [1] "Level: Ia; Type: elnet.dev; Data: 2 - nBeta: 83 BER: 0.053 aveAcc: 0.971 aveKappa: 0.908"
# [1] "Level: Ia; Type: elnet.class; Data: 2 - nBeta: 84 BER: 0.052 aveAcc: 0.971 aveKappa: 0.909"
# 
# [1] "Level: Ia; Type: lasso.dev; Data: 3 - nBeta: 77 BER: 0.031 aveAcc: 0.98 aveKappa: 0.939"
# [1] "Level: Ia; Type: lasso.class; Data: 3 - nBeta: 81 BER: 0.031 aveAcc: 0.98 aveKappa: 0.939"
# [1] "Level: Ia; Type: elnet.dev; Data: 3 - nBeta: 82 BER: 0.031 aveAcc: 0.981 aveKappa: 0.94"
# [1] "Level: Ia; Type: elnet.class; Data: 3 - nBeta: 82 BER: 0.03 aveAcc: 0.98 aveKappa: 0.94"
# 
# [1] "Level: Ia; Type: lasso.dev; Data: 4 - nBeta: 80 BER: 0.045 aveAcc: 0.972 aveKappa: 0.914"
# [1] "Level: Ia; Type: lasso.class; Data: 4 - nBeta: 86 BER: 0.043 aveAcc: 0.973 aveKappa: 0.916"
# [1] "Level: Ia; Type: elnet.dev; Data: 4 - nBeta: 76 BER: 0.046 aveAcc: 0.972 aveKappa: 0.914"
# [1] "Level: Ia; Type: elnet.class; Data: 4 - nBeta: 87 BER: 0.043 aveAcc: 0.973 aveKappa: 0.917"
# 
# 1 granulocyte        114
# 2 basophil             3
# 3 monocyte            62
# 4 macrophage          56
# 5 Kupffer cell         9
# 6 microglial cell    659
# 7 brain pericyte      23
# [1] "Level: Ib; Type: lasso.dev; Data: 1 - nBeta: 42 BER: 0.028 aveAcc: 0.994 aveKappa: 0.986"
# [1] "Level: Ib; Type: lasso.class; Data: 1 - nBeta: 75 BER: 0.03 aveAcc: 0.992 aveKappa: 0.984"
# [1] "Level: Ib; Type: elnet.dev; Data: 1 - nBeta: 74 BER: 0.03 aveAcc: 0.992 aveKappa: 0.984"
# [1] "Level: Ib; Type: elnet.class; Data: 1 - nBeta: 25 BER: 0.028 aveAcc: 0.994 aveKappa: 0.986"
# 
# [1] "Level: Ib; Type: lasso.dev; Data: 2 - nBeta: 59 BER: 0.101 aveAcc: 0.986 aveKappa: 0.97"
# [1] "Level: Ib; Type: lasso.class; Data: 2 - nBeta: 66 BER: 0.101 aveAcc: 0.986 aveKappa: 0.97"
# [1] "Level: Ib; Type: elnet.dev; Data: 2 - nBeta: 62 BER: 0.1 aveAcc: 0.987 aveKappa: 0.972"
# [1] "Level: Ib; Type: elnet.class; Data: 2 - nBeta: 61 BER: 0.101 aveAcc: 0.986 aveKappa: 0.97"
# 
# [1] "Level: Ib; Type: lasso.dev; Data: 3 - nBeta: 46 BER: 0.033 aveAcc: 0.991 aveKappa: 0.982"
# [1] "Level: Ib; Type: lasso.class; Data: 3 - nBeta: 39 BER: 0.033 aveAcc: 0.991 aveKappa: 0.982"
# [1] "Level: Ib; Type: elnet.dev; Data: 3 - nBeta: 50 BER: 0.033 aveAcc: 0.991 aveKappa: 0.982"
# [1] "Level: Ib; Type: elnet.class; Data: 3 - nBeta: 37 BER: 0.033 aveAcc: 0.991 aveKappa: 0.982"
# 
# [1] "Level: Ib; Type: lasso.dev; Data: 4 - nBeta: 65 BER: 0.095 aveAcc: 0.988 aveKappa: 0.975"
# [1] "Level: Ib; Type: lasso.class; Data: 4 - nBeta: 76 BER: 0.077 aveAcc: 0.99 aveKappa: 0.979"
# [1] "Level: Ib; Type: elnet.dev; Data: 4 - nBeta: 71 BER: 0.095 aveAcc: 0.988 aveKappa: 0.975"
# [1] "Level: Ib; Type: elnet.class; Data: 4 - nBeta: 71 BER: 0.095 aveAcc: 0.988 aveKappa: 0.975"
# 
# 1 macrophage - Marrow                25
# 2 macrophage - Spleen                 7
# 3 macrophage - Kidney                 5
# 4 macrophage - Limb_Muscle            6
# 5 macrophage - Diaphragm              4
# 6 macrophage - Brain_Myeloid          9
# 7 Kupffer cell - Liver                9
# 8 microglial cell - Brain_Myeloid   659
# [1] "Level: Ic; Type: lasso.dev; Data: 1 - nBeta: 42 BER: 0.163 aveAcc: 0.99 aveKappa: 0.943"
# [1] "Level: Ic; Type: lasso.class; Data: 1 - nBeta: 34 BER: 0.111 aveAcc: 0.993 aveKappa: 0.959"
# [1] "Level: Ic; Type: elnet.dev; Data: 1 - nBeta: 56 BER: 0.163 aveAcc: 0.99 aveKappa: 0.943"
# [1] "Level: Ic; Type: elnet.class; Data: 1 - nBeta: 31 BER: 0.111 aveAcc: 0.993 aveKappa: 0.959"
# 
# [1] "Level: Ic; Type: lasso.dev; Data: 2 - nBeta: 46 BER: 0.224 aveAcc: 0.985 aveKappa: 0.91"
# [1] "Level: Ic; Type: lasso.class; Data: 2 - nBeta: 57 BER: 0.203 aveAcc: 0.986 aveKappa: 0.918"
# [1] "Level: Ic; Type: elnet.dev; Data: 2 - nBeta: 51 BER: 0.203 aveAcc: 0.986 aveKappa: 0.918"
# [1] "Level: Ic; Type: elnet.class; Data: 2 - nBeta: 71 BER: 0.217 aveAcc: 0.985 aveKappa: 0.91"
# 
# [1] "Level: Ic; Type: lasso.dev; Data: 3 - nBeta: 35 BER: 0.136 aveAcc: 0.992 aveKappa: 0.951"
# [1] "Level: Ic; Type: lasso.class; Data: 3 - nBeta: 55 BER: 0.111 aveAcc: 0.993 aveKappa: 0.959"
# [1] "Level: Ic; Type: elnet.dev; Data: 3 - nBeta: 49 BER: 0.136 aveAcc: 0.992 aveKappa: 0.951"
# [1] "Level: Ic; Type: elnet.class; Data: 3 - nBeta: 49 BER: 0.136 aveAcc: 0.992 aveKappa: 0.951"
# 
# [1] "Level: Ic; Type: lasso.dev; Data: 4 - nBeta: 61 BER: 0.158 aveAcc: 0.989 aveKappa: 0.934"
# [1] "Level: Ic; Type: lasso.class; Data: 4 - nBeta: 80 BER: 0.08 aveAcc: 0.993 aveKappa: 0.959"
# [1] "Level: Ic; Type: elnet.dev; Data: 4 - nBeta: 55 BER: 0.203 aveAcc: 0.986 aveKappa: 0.918"
# [1] "Level: Ic; Type: elnet.class; Data: 4 - nBeta: 71 BER: 0.08 aveAcc: 0.994 aveKappa: 0.967"
# 
# 1 granulocyte        114
# 2 basophil             3
# 3 monocyte            62
# 4 macrophage          56
# 5 Kupffer cell         9
# 6 microglial cell    659
# 7 brain pericyte      23
# 8 non myeloid       3595
# [1] "Level: Inc; Type: lasso.dev; Data: 1 - nBeta: 29 BER: 0.347 aveAcc: 0.979 aveKappa: 0.938"
# [1] "Level: Inc; Type: lasso.class; Data: 1 - nBeta: 19 BER: 0.403 aveAcc: 0.975 aveKappa: 0.927"
# [1] "Level: Inc; Type: elnet.dev; Data: 1 - nBeta: 80 BER: 0.096 aveAcc: 0.988 aveKappa: 0.964"
# [1] "Level: Inc; Type: elnet.class; Data: 1 - nBeta: 77 BER: 0.112 aveAcc: 0.987 aveKappa: 0.963"
# 
# [1] "Level: Inc; Type: lasso.dev; Data: 2 - nBeta: 11 BER: 0.614 aveAcc: 0.957 aveKappa: 0.864"
# [1] "Level: Inc; Type: lasso.class; Data: 2 - nBeta: 11 BER: 0.618 aveAcc: 0.956 aveKappa: 0.863"
# [1] "Level: Inc; Type: elnet.dev; Data: 2 - nBeta: 35 BER: 0.496 aveAcc: 0.962 aveKappa: 0.884"
# [1] "Level: Inc; Type: elnet.class; Data: 2 - nBeta: 78 BER: 0.346 aveAcc: 0.97 aveKappa: 0.91"
# 
# [1] "Level: Inc; Type: lasso.dev; Data: 3 - nBeta: 45 BER: 0.193 aveAcc: 0.977 aveKappa: 0.933"
# [1] "Level: Inc; Type: lasso.class; Data: 3 - nBeta: 54 BER: 0.173 aveAcc: 0.978 aveKappa: 0.937"
# [1] "Level: Inc; Type: elnet.dev; Data: 3 - nBeta: 89 BER: 0.109 aveAcc: 0.984 aveKappa: 0.954"
# [1] "Level: Inc; Type: elnet.class; Data: 3 - nBeta: 89 BER: 0.109 aveAcc: 0.984 aveKappa: 0.954"
# 
# [1] "Level: Inc; Type: lasso.dev; Data: 4 - nBeta: 12 BER: 0.612 aveAcc: 0.957 aveKappa: 0.865"
# [1] "Level: Inc; Type: lasso.class; Data: 4 - nBeta: 11 BER: 0.619 aveAcc: 0.957 aveKappa: 0.864"
# [1] "Level: Inc; Type: elnet.dev; Data: 4 - nBeta: 53 BER: 0.456 aveAcc: 0.967 aveKappa: 0.901"
# [1] "Level: Inc; Type: elnet.class; Data: 4 - nBeta: 45 BER: 0.472 aveAcc: 0.965 aveKappa: 0.896"


# DROP -------------------------------------------------------------------------------------------#

# 1 non myeloid  3131
# 2 myeloid       550
# [1] "Level: Ia; Type: lasso.dev; Data: 1 - nBeta: 56 BER: 0.013 aveAcc: 0.992 aveKappa: 0.969"
# [1] "Level: Ia; Type: lasso.class; Data: 1 - nBeta: 69 BER: 0.014 aveAcc: 0.992 aveKappa: 0.969"
# [1] "Level: Ia; Type: elnet.dev; Data: 1 - nBeta: 60 BER: 0.013 aveAcc: 0.992 aveKappa: 0.969"
# [1] "Level: Ia; Type: elnet.class; Data: 1 - nBeta: 63 BER: 0.013 aveAcc: 0.992 aveKappa: 0.97"
# 
# [1] "Level: Ia; Type: lasso.dev; Data: 2 - nBeta: 56 BER: 0.019 aveAcc: 0.99 aveKappa: 0.961"
# [1] "Level: Ia; Type: lasso.class; Data: 2 - nBeta: 74 BER: 0.018 aveAcc: 0.99 aveKappa: 0.963"
# [1] "Level: Ia; Type: elnet.dev; Data: 2 - nBeta: 55 BER: 0.018 aveAcc: 0.99 aveKappa: 0.96"
# [1] "Level: Ia; Type: elnet.class; Data: 2 - nBeta: 59 BER: 0.019 aveAcc: 0.99 aveKappa: 0.961"
# 
# [1] "Level: Ia; Type: lasso.dev; Data: 3 - nBeta: 53 BER: 0.02 aveAcc: 0.99 aveKappa: 0.959"
# [1] "Level: Ia; Type: lasso.class; Data: 3 - nBeta: 66 BER: 0.019 aveAcc: 0.99 aveKappa: 0.962"
# [1] "Level: Ia; Type: elnet.dev; Data: 3 - nBeta: 54 BER: 0.019 aveAcc: 0.99 aveKappa: 0.96"
# [1] "Level: Ia; Type: elnet.class; Data: 3 - nBeta: 64 BER: 0.019 aveAcc: 0.99 aveKappa: 0.963"
# 
# [1] "Level: Ia; Type: lasso.dev; Data: 4 - nBeta: 50 BER: 0.02 aveAcc: 0.989 aveKappa: 0.958"
# [1] "Level: Ia; Type: lasso.class; Data: 4 - nBeta: 67 BER: 0.02 aveAcc: 0.99 aveKappa: 0.96"
# [1] "Level: Ia; Type: elnet.dev; Data: 4 - nBeta: 57 BER: 0.021 aveAcc: 0.989 aveKappa: 0.958"
# [1] "Level: Ia; Type: elnet.class; Data: 4 - nBeta: 66 BER: 0.021 aveAcc: 0.99 aveKappa: 0.959"
# 
# 1 granulocyte        108
# 2 basophil             9
# 3 promonocyte         38
# 4 monocyte           135
# 5 macrophage         246
# 6 dendritic cell       6
# 7 Langerhans cell      5
# 8 mast cell            3
# [1] "Level: Ib; Type: lasso.dev; Data: 1 - nBeta: 81 BER: 0.052 aveAcc: 0.978 aveKappa: 0.969"
# [1] "Level: Ib; Type: lasso.class; Data: 1 - nBeta: 80 BER: 0.052 aveAcc: 0.978 aveKappa: 0.969"
# [1] "Level: Ib; Type: elnet.dev; Data: 1 - nBeta: 82 BER: 0.052 aveAcc: 0.978 aveKappa: 0.969"
# [1] "Level: Ib; Type: elnet.class; Data: 1 - nBeta: 74 BER: 0.053 aveAcc: 0.976 aveKappa: 0.966"
# 
# [1] "Level: Ib; Type: lasso.dev; Data: 2 - nBeta: 86 BER: 0.097 aveAcc: 0.965 aveKappa: 0.951"
# [1] "Level: Ib; Type: lasso.class; Data: 2 - nBeta: 88 BER: 0.096 aveAcc: 0.967 aveKappa: 0.953"
# [1] "Level: Ib; Type: elnet.dev; Data: 2 - nBeta: 88 BER: 0.096 aveAcc: 0.965 aveKappa: 0.951"
# [1] "Level: Ib; Type: elnet.class; Data: 2 - nBeta: 88 BER: 0.096 aveAcc: 0.965 aveKappa: 0.951"
# 
# [1] "Level: Ib; Type: lasso.dev; Data: 3 - nBeta: 82 BER: 0.141 aveAcc: 0.965 aveKappa: 0.95"
# [1] "Level: Ib; Type: lasso.class; Data: 3 - nBeta: 94 BER: 0.097 aveAcc: 0.969 aveKappa: 0.956"
# [1] "Level: Ib; Type: elnet.dev; Data: 3 - nBeta: 82 BER: 0.141 aveAcc: 0.964 aveKappa: 0.948"
# [1] "Level: Ib; Type: elnet.class; Data: 3 - nBeta: 79 BER: 0.14 aveAcc: 0.965 aveKappa: 0.95"
# 
# [1] "Level: Ib; Type: lasso.dev; Data: 4 - nBeta: 74 BER: 0.227 aveAcc: 0.944 aveKappa: 0.919"
# [1] "Level: Ib; Type: lasso.class; Data: 4 - nBeta: 85 BER: 0.059 aveAcc: 0.967 aveKappa: 0.953"
# [1] "Level: Ib; Type: elnet.dev; Data: 4 - nBeta: 82 BER: 0.123 aveAcc: 0.958 aveKappa: 0.94"
# [1] "Level: Ib; Type: elnet.class; Data: 4 - nBeta: 90 BER: 0.059 aveAcc: 0.967 aveKappa: 0.953"
# 
# 1 macrophage - Marrow           33
# 2 macrophage - Spleen           69
# 3 macrophage - Lung             51
# 4 macrophage - Kidney           20
# 5 macrophage - Mammary_Gland    27
# 6 macrophage - Limb_Muscle      46
# [1] "Level: Ic; Type: lasso.dev; Data: 1 - nBeta: 58 BER: 0.009 aveAcc: 0.992 aveKappa: 0.99"
# [1] "Level: Ic; Type: lasso.class; Data: 1 - nBeta: 57 BER: 0.009 aveAcc: 0.992 aveKappa: 0.99"
# [1] "Level: Ic; Type: elnet.dev; Data: 1 - nBeta: 67 BER: 0.004 aveAcc: 0.996 aveKappa: 0.995"
# [1] "Level: Ic; Type: elnet.class; Data: 1 - nBeta: 54 BER: 0.018 aveAcc: 0.984 aveKappa: 0.98"
# 
# [1] "Level: Ic; Type: lasso.dev; Data: 2 - nBeta: 80 BER: 0.05 aveAcc: 0.955 aveKappa: 0.945"
# [1] "Level: Ic; Type: lasso.class; Data: 2 - nBeta: 80 BER: 0.05 aveAcc: 0.955 aveKappa: 0.945"
# [1] "Level: Ic; Type: elnet.dev; Data: 2 - nBeta: 83 BER: 0.05 aveAcc: 0.955 aveKappa: 0.945"
# [1] "Level: Ic; Type: elnet.class; Data: 2 - nBeta: 84 BER: 0.045 aveAcc: 0.959 aveKappa: 0.95"
# 
# [1] "Level: Ic; Type: lasso.dev; Data: 3 - nBeta: 69 BER: 0.04 aveAcc: 0.963 aveKappa: 0.955"
# [1] "Level: Ic; Type: lasso.class; Data: 3 - nBeta: 91 BER: 0.032 aveAcc: 0.972 aveKappa: 0.965"
# [1] "Level: Ic; Type: elnet.dev; Data: 3 - nBeta: 80 BER: 0.033 aveAcc: 0.972 aveKappa: 0.965"
# [1] "Level: Ic; Type: elnet.class; Data: 3 - nBeta: 78 BER: 0.037 aveAcc: 0.967 aveKappa: 0.96"
# 
# [1] "Level: Ic; Type: lasso.dev; Data: 4 - nBeta: 76 BER: 0.067 aveAcc: 0.943 aveKappa: 0.929"
# [1] "Level: Ic; Type: lasso.class; Data: 4 - nBeta: 81 BER: 0.071 aveAcc: 0.939 aveKappa: 0.924"
# [1] "Level: Ic; Type: elnet.dev; Data: 4 - nBeta: 84 BER: 0.071 aveAcc: 0.939 aveKappa: 0.924"
# [1] "Level: Ic; Type: elnet.class; Data: 4 - nBeta: 80 BER: 0.067 aveAcc: 0.943 aveKappa: 0.929"
# 
# 1 granulocyte        108
# 2 basophil             9
# 3 promonocyte         38
# 4 monocyte           135
# 5 macrophage         246
# 6 dendritic cell       6
# 7 Langerhans cell      5
# 8 mast cell            3
# 9 non myeloid       3131
# [1] "Level: Inc; Type: lasso.dev; Data: 1 - nBeta: 80 BER: 0.069 aveAcc: 0.991 aveKappa: 0.967"
# [1] "Level: Inc; Type: lasso.class; Data: 1 - nBeta: 89 BER: 0.048 aveAcc: 0.993 aveKappa: 0.975"
# [1] "Level: Inc; Type: elnet.dev; Data: 1 - nBeta: 81 BER: 0.069 aveAcc: 0.991 aveKappa: 0.967"
# [1] "Level: Inc; Type: elnet.class; Data: 1 - nBeta: 95 BER: 0.066 aveAcc: 0.993 aveKappa: 0.974"
# 
# [1] "Level: Inc; Type: lasso.dev; Data: 2 - nBeta: 45 BER: 0.277 aveAcc: 0.978 aveKappa: 0.917"
# [1] "Level: Inc; Type: lasso.class; Data: 2 - nBeta: 69 BER: 0.203 aveAcc: 0.981 aveKappa: 0.931"
# [1] "Level: Inc; Type: elnet.dev; Data: 2 - nBeta: 82 BER: 0.142 aveAcc: 0.984 aveKappa: 0.941"
# [1] "Level: Inc; Type: elnet.class; Data: 2 - nBeta: 87 BER: 0.142 aveAcc: 0.984 aveKappa: 0.941"
# 
# [1] "Level: Inc; Type: lasso.dev; Data: 3 - nBeta: 82 BER: 0.112 aveAcc: 0.987 aveKappa: 0.951"
# [1] "Level: Inc; Type: lasso.class; Data: 3 - nBeta: 88 BER: 0.11 aveAcc: 0.988 aveKappa: 0.956"
# [1] "Level: Inc; Type: elnet.dev; Data: 3 - nBeta: 82 BER: 0.114 aveAcc: 0.986 aveKappa: 0.948"
# [1] "Level: Inc; Type: elnet.class; Data: 3 - nBeta: 93 BER: 0.128 aveAcc: 0.989 aveKappa: 0.959"
# 
# [1] "Level: Inc; Type: lasso.dev; Data: 4 - nBeta: 32 BER: 0.381 aveAcc: 0.972 aveKappa: 0.896"
# [1] "Level: Inc; Type: lasso.class; Data: 4 - nBeta: 36 BER: 0.378 aveAcc: 0.974 aveKappa: 0.901"
# [1] "Level: Inc; Type: elnet.dev; Data: 4 - nBeta: 86 BER: 0.239 aveAcc: 0.981 aveKappa: 0.93"
# [1] "Level: Inc; Type: elnet.class; Data: 4 - nBeta: 92 BER: 0.108 aveAcc: 0.983 aveKappa: 0.938"


# Cell type breakdown --- FACS:
for(i in c("Ib", "Ic")) {
  
  for (mod in c("lasso.dev", "lasso.class", "elnet.dev", "elnet.class")) {
    
    print(paste0("Level: ", i, "  Type: ", mod, " -- nBetas:"))
    print(t(ica_100_8test_glmnet_stat_facs[[i]][[mod]]$nbetas))
    
    print(paste0("Level: ", i, "  Type: ", mod, " -- Acc:"))
    print(t(ica_100_8test_glmnet_stat_facs[[i]][[mod]]$bacc))
    print("----------")
  }
}
# [1] "Level: Ib  Type: lasso.dev -- nBetas:"
# VSt.log VSt.rc SCt.res SCt.rc
# granulocyte           9     10       6      9
# basophil              4      4       4      5
# monocyte             10     15      11     18
# macrophage           21     30      24     33
# Kupffer cell          4      7       4      6
# microglial cell      12     13      11     13
# brain pericyte        2      8       4      8
# [1] "Level: Ib  Type: lasso.dev -- Acc:"
# VSt.log VSt.rc SCt.res SCt.rc
# granulocyte       0.999  0.990   0.999  0.994
# basophil          1.000  0.833   1.000  0.833
# monocyte          0.983  0.957   0.967  0.974
# macrophage        0.973  0.970   0.971  0.970
# Kupffer cell      0.944  0.889   0.944  0.889
# microglial cell   0.998  0.998   0.998  0.997
# brain pericyte    1.000  1.000   1.000  1.000
# [1] "----------"
# [1] "Level: Ib  Type: lasso.class -- nBetas:"
# VSt.log VSt.rc SCt.res SCt.rc
# granulocyte          16     14       7     18
# basophil              5      5       2      7
# monocyte             28     22       6     26
# macrophage           34     36      19     33
# Kupffer cell         12      7       3      9
# microglial cell      15     17       9     18
# brain pericyte        2      8       3      8
# [1] "Level: Ib  Type: lasso.class -- Acc:"
# VSt.log VSt.rc SCt.res SCt.rc
# granulocyte       0.999  0.990   0.999  0.994
# basophil          1.000  0.833   1.000  0.833
# monocyte          0.975  0.958   0.967  0.982
# macrophage        0.972  0.970   0.971  0.971
# Kupffer cell      0.944  0.889   0.944  0.944
# microglial cell   0.998  0.998   0.998  0.997
# brain pericyte    1.000  1.000   1.000  1.000
# [1] "----------"
# [1] "Level: Ib  Type: elnet.dev -- nBetas:"
# VSt.log VSt.rc SCt.res SCt.rc
# granulocyte          16     13       8     13
# basophil              5      4       5      6
# monocyte             27     17      12     21
# macrophage           34     32      26     34
# Kupffer cell         10      7       6      7
# microglial cell      15     16      15     17
# brain pericyte        2      8       6      8
# [1] "Level: Ib  Type: elnet.dev -- Acc:"
# VSt.log VSt.rc SCt.res SCt.rc
# granulocyte       0.999  0.994   0.999  0.994
# basophil          1.000  0.833   1.000  0.833
# monocyte          0.975  0.958   0.967  0.974
# macrophage        0.972  0.970   0.971  0.970
# Kupffer cell      0.944  0.889   0.944  0.889
# microglial cell   0.998  0.998   0.998  0.997
# brain pericyte    1.000  1.000   1.000  1.000
# [1] "----------"
# [1] "Level: Ib  Type: elnet.class -- nBetas:"
# VSt.log VSt.rc SCt.res SCt.rc
# granulocyte           4     11       5     15
# basophil              4      4       2      6
# monocyte              3     15       6     22
# macrophage            8     29      16     33
# Kupffer cell          1      6       4      8
# microglial cell       8     15       9     18
# brain pericyte        2      8       4      8
# [1] "Level: Ib  Type: elnet.class -- Acc:"
# VSt.log VSt.rc SCt.res SCt.rc
# granulocyte       0.999  0.990   0.999  0.994
# basophil          1.000  0.833   1.000  0.833
# monocyte          0.983  0.957   0.967  0.974
# macrophage        0.973  0.970   0.971  0.970
# Kupffer cell      0.944  0.889   0.944  0.889
# microglial cell   0.998  0.998   0.998  0.997
# brain pericyte    1.000  1.000   1.000  1.000
# [1] "----------"
# [1] "Level: Ic  Type: lasso.dev -- nBetas:"
# VSt.log VSt.rc SCt.res SCt.rc
# macrophage - Marrow                   3      5       3     10
# macrophage - Spleen                   1      1       2      2
# macrophage - Kidney                   7      7       9      4
# macrophage - Limb_Muscle             16     13      10     17
# macrophage - Diaphragm               10     10      11     11
# macrophage - Brain_Myeloid           18      8       8     19
# Kupffer cell - Liver                  5      8       4     12
# microglial cell - Brain_Myeloid       9     15       9     18
# [1] "Level: Ic  Type: lasso.dev -- Acc:"
# VSt.log VSt.rc SCt.res SCt.rc
# macrophage - Marrow               1.000  0.999   1.000  1.000
# macrophage - Spleen               1.000  1.000   1.000  1.000
# macrophage - Kidney               1.000  0.797   0.900  0.799
# macrophage - Limb_Muscle          0.831  0.833   0.915  0.915
# macrophage - Diaphragm            0.624  0.749   0.749  0.874
# macrophage - Brain_Myeloid        0.944  0.776   0.944  0.832
# Kupffer cell - Liver              0.944  0.944   0.944  0.944
# microglial cell - Brain_Myeloid   0.992  0.992   0.992  0.992
# [1] "----------"
# [1] "Level: Ic  Type: lasso.class -- nBetas:"
# VSt.log VSt.rc SCt.res SCt.rc
# macrophage - Marrow                   3      7       5     12
# macrophage - Spleen                   1      1       2      4
# macrophage - Kidney                   7      8      12     17
# macrophage - Limb_Muscle             12     18      18     24
# macrophage - Diaphragm                6     15      19     23
# macrophage - Brain_Myeloid           15     15      13     43
# Kupffer cell - Liver                  4      9       9     14
# microglial cell - Brain_Myeloid       8     16      10     24
# [1] "Level: Ic  Type: lasso.class -- Acc:"
# VSt.log VSt.rc SCt.res SCt.rc
# macrophage - Marrow               1.000  1.000   1.000  1.000
# macrophage - Spleen               1.000  1.000   1.000  1.000
# macrophage - Kidney               1.000  0.797   1.000  1.000
# macrophage - Limb_Muscle          0.915  0.915   0.915  0.915
# macrophage - Diaphragm            0.749  0.749   0.749  0.874
# macrophage - Brain_Myeloid        0.944  0.776   0.944  0.944
# Kupffer cell - Liver              0.944  0.944   0.944  0.944
# microglial cell - Brain_Myeloid   0.992  0.992   0.992  0.992
# [1] "----------"
# [1] "Level: Ic  Type: elnet.dev -- nBetas:"
# VSt.log VSt.rc SCt.res SCt.rc
# macrophage - Marrow                   7      6       5     10
# macrophage - Spleen                   1      1       3      1
# macrophage - Kidney                   9      7      11      3
# macrophage - Limb_Muscle             19     16      15     16
# macrophage - Diaphragm               14     14      14     11
# macrophage - Brain_Myeloid           23     11      12     14
# Kupffer cell - Liver                  9      9       7     11
# microglial cell - Brain_Myeloid      12     16       9     16
# [1] "Level: Ic  Type: elnet.dev -- Acc:"
# VSt.log VSt.rc SCt.res SCt.rc
# macrophage - Marrow               1.000  0.999   1.000  1.000
# macrophage - Spleen               1.000  1.000   1.000  1.000
# macrophage - Kidney               1.000  0.797   0.900  0.798
# macrophage - Limb_Muscle          0.831  0.916   0.915  0.915
# macrophage - Diaphragm            0.624  0.749   0.749  0.749
# macrophage - Brain_Myeloid        0.944  0.776   0.944  0.776
# Kupffer cell - Liver              0.944  0.944   0.944  0.944
# microglial cell - Brain_Myeloid   0.992  0.992   0.992  0.992
# [1] "----------"
# [1] "Level: Ic  Type: elnet.class -- nBetas:"
# VSt.log VSt.rc SCt.res SCt.rc
# macrophage - Marrow                   4      9       5     11
# macrophage - Spleen                   1      3       3      4
# macrophage - Kidney                   6     12      11     15
# macrophage - Limb_Muscle              9     22      15     25
# macrophage - Diaphragm                4     20      14     19
# macrophage - Brain_Myeloid           14     24      12     27
# Kupffer cell - Liver                  3     13       7     15
# microglial cell - Brain_Myeloid      10     23       9     22
# [1] "Level: Ic  Type: elnet.class -- Acc:"
# VSt.log VSt.rc SCt.res SCt.rc
# macrophage - Marrow               1.000  1.000   1.000  1.000
# macrophage - Spleen               1.000  1.000   1.000  1.000
# macrophage - Kidney               1.000  0.797   0.900  1.000
# macrophage - Limb_Muscle          0.915  0.915   0.915  0.915
# macrophage - Diaphragm            0.749  0.749   0.749  0.874
# macrophage - Brain_Myeloid        0.944  0.721   0.944  0.944
# Kupffer cell - Liver              0.944  0.944   0.944  0.944
# microglial cell - Brain_Myeloid   0.992  0.992   0.992  0.992
# [1] "----------"


# Cell type breakdown --- FACS:
for(i in c("Ib", "Ic")) {
  
  for (mod in c("lasso.dev", "lasso.class", "elnet.dev", "elnet.class")) {
    
    print(paste0("Level: ", i, "  Type: ", mod, " -- nBetas:"))
    print(t(ica_100_8test_glmnet_stat_drop[[i]][[mod]]$nbetas))
    
    print(paste0("Level: ", i, "  Type: ", mod, " -- Acc:"))
    print(t(ica_100_8test_glmnet_stat_drop[[i]][[mod]]$bacc))
    print("----------")
  }
}
# [1] "Level: Ib  Type: lasso.dev -- nBetas:"
# VSt.log VSt.rc SCt.res SCt.rc
# granulocyte           9     12       9      7
# basophil              6     10       5      6
# promonocyte           9     17      10     10
# monocyte             43     42      49     39
# macrophage           48     56      51     38
# dendritic cell       12     13      11      8
# Langerhans cell       3     16       7      9
# mast cell             6      8      10      5
# [1] "Level: Ib  Type: lasso.dev -- Acc:"
# VSt.log VSt.rc SCt.res SCt.rc
# granulocyte       0.998  0.993   0.998  0.998
# basophil          1.000  1.000   1.000  1.000
# promonocyte       0.985  0.984   0.973  0.984
# monocyte          0.978  0.974   0.971  0.951
# macrophage        0.985  0.973   0.973  0.952
# dendritic cell    0.999  0.831   0.832  0.665
# Langerhans cell   1.000  1.000   1.000  1.000
# mast cell         0.832  0.832   0.666  0.499
# [1] "----------"
# [1] "Level: Ib  Type: lasso.class -- nBetas:"
# VSt.log VSt.rc SCt.res SCt.rc
# granulocyte           9     14      15     15
# basophil              6     10       6     12
# promonocyte          10     19      21     23
# monocyte             43     48      56     57
# macrophage           45     59      61     60
# dendritic cell       11     15      18     18
# Langerhans cell       3     16       8     11
# mast cell             6      8      11     11
# [1] "Level: Ib  Type: lasso.class -- Acc:"
# VSt.log VSt.rc SCt.res SCt.rc
# granulocyte       0.998  0.998   0.998  0.998
# basophil          1.000  1.000   1.000  1.000
# promonocyte       0.985  0.984   0.972  0.971
# monocyte          0.978  0.975   0.980  0.966
# macrophage        0.985  0.973   0.975  0.979
# dendritic cell    0.999  0.831   0.832  0.998
# Langerhans cell   1.000  1.000   1.000  1.000
# mast cell         0.832  0.832   0.832  0.832
# [1] "----------"
# [1] "Level: Ib  Type: elnet.dev -- nBetas:"
# VSt.log VSt.rc SCt.res SCt.rc
# granulocyte          10     18       9     13
# basophil              8     11       5     10
# promonocyte          10     20       8     17
# monocyte             46     49      52     48
# macrophage           48     62      51     54
# dendritic cell       12     17      11     18
# Langerhans cell       5     17       7     12
# mast cell             7      8      10      8
# [1] "Level: Ib  Type: elnet.dev -- Acc:"
# VSt.log VSt.rc SCt.res SCt.rc
# granulocyte       0.998  0.998   0.998  0.998
# basophil          1.000  1.000   1.000  1.000
# promonocyte       0.985  0.984   0.973  0.971
# monocyte          0.978  0.975   0.969  0.962
# macrophage        0.985  0.973   0.971  0.969
# dendritic cell    0.999  0.830   0.832  0.748
# Langerhans cell   1.000  1.000   1.000  1.000
# mast cell         0.832  0.832   0.666  0.832
# [1] "----------"
# [1] "Level: Ib  Type: elnet.class -- nBetas:"
# VSt.log VSt.rc SCt.res SCt.rc
# granulocyte          10     17       8     18
# basophil              6     11       4     19
# promonocyte          10     19       7     32
# monocyte             38     49      43     64
# macrophage           44     60      48     66
# dendritic cell       11     17      10     27
# Langerhans cell       5     17       6     14
# mast cell             6      8       8     12
# [1] "Level: Ib  Type: elnet.class -- Acc:"
# VSt.log VSt.rc SCt.res SCt.rc
# granulocyte       0.998  0.998   0.998  0.998
# basophil          1.000  1.000   1.000  1.000
# promonocyte       0.985  0.984   0.973  0.971
# monocyte          0.974  0.975   0.974  0.966
# macrophage        0.983  0.971   0.971  0.979
# dendritic cell    0.999  0.831   0.832  0.998
# Langerhans cell   1.000  1.000   1.000  1.000
# mast cell         0.832  0.832   0.666  0.832
# [1] "----------"
# [1] "Level: Ic  Type: lasso.dev -- nBetas:"
# VSt.log VSt.rc SCt.res SCt.rc
# macrophage - Marrow              5     17       6     17
# macrophage - Spleen             23     40      20     26
# macrophage - Lung                2      4       2      3
# macrophage - Kidney              8     17      10     11
# macrophage - Mammary_Gland      21     30      36     26
# macrophage - Limb_Muscle        17     31      28     30
# [1] "Level: Ic  Type: lasso.dev -- Acc:"
# VSt.log VSt.rc SCt.res SCt.rc
# macrophage - Marrow          0.985  0.980   0.998  0.980
# macrophage - Spleen          0.994  0.987   0.994  0.992
# macrophage - Lung            1.000  1.000   1.000  1.000
# macrophage - Kidney          1.000  0.996   0.996  0.993
# macrophage - Mammary_Gland   1.000  0.917   0.938  0.880
# macrophage - Limb_Muscle     0.989  0.943   0.932  0.919
# [1] "----------"
# [1] "Level: Ic  Type: lasso.class -- nBetas:"
# VSt.log VSt.rc SCt.res SCt.rc
# macrophage - Marrow              5     17      11     17
# macrophage - Spleen             25     40      34     31
# macrophage - Lung                2      4       5      4
# macrophage - Kidney              9     17      10     13
# macrophage - Mammary_Gland      19     30      46     33
# macrophage - Limb_Muscle        17     31      43     34
# [1] "Level: Ic  Type: lasso.class -- Acc:"
# VSt.log VSt.rc SCt.res SCt.rc
# macrophage - Marrow          0.985  0.980   0.985  0.980
# macrophage - Spleen          0.994  0.987   0.994  0.992
# macrophage - Lung            1.000  1.000   1.000  0.997
# macrophage - Kidney          1.000  0.996   0.998  0.993
# macrophage - Mammary_Gland   1.000  0.917   0.956  0.880
# macrophage - Limb_Muscle     0.989  0.943   0.954  0.908
# [1] "----------"
# [1] "Level: Ic  Type: elnet.dev -- nBetas:"
# VSt.log VSt.rc SCt.res SCt.rc
# macrophage - Marrow              6     20       6     21
# macrophage - Spleen             26     41      26     33
# macrophage - Lung                3      9       4      9
# macrophage - Kidney             10     18      11     13
# macrophage - Mammary_Gland      25     34      41     34
# macrophage - Limb_Muscle        24     33      35     35
# [1] "Level: Ic  Type: elnet.dev -- Acc:"
# VSt.log VSt.rc SCt.res SCt.rc
# macrophage - Marrow          1.000  0.980   0.998  0.980
# macrophage - Spleen          0.997  0.987   0.997  0.992
# macrophage - Lung            1.000  1.000   1.000  1.000
# macrophage - Kidney          1.000  0.996   0.998  0.993
# macrophage - Mammary_Gland   1.000  0.917   0.938  0.877
# macrophage - Limb_Muscle     0.989  0.943   0.954  0.908
# [1] "----------"
# [1] "Level: Ic  Type: elnet.class -- nBetas:"
# VSt.log VSt.rc SCt.res SCt.rc
# macrophage - Marrow              5     20       6     18
# macrophage - Spleen             23     41      24     28
# macrophage - Lung                3      9       4      9
# macrophage - Kidney             10     18      10     12
# macrophage - Mammary_Gland      17     35      39     29
# macrophage - Limb_Muscle        14     35      34     31
# [1] "Level: Ic  Type: elnet.class -- Acc:"
# VSt.log VSt.rc SCt.res SCt.rc
# macrophage - Marrow          0.985  0.995   0.998  0.980
# macrophage - Spleen          0.994  0.987   0.994  0.992
# macrophage - Lung            1.000  1.000   1.000  1.000
# macrophage - Kidney          0.998  0.996   0.998  0.993
# macrophage - Mammary_Gland   0.979  0.919   0.938  0.880
# macrophage - Limb_Muscle     0.978  0.943   0.943  0.919
# [1] "----------"



#-------------------- Now plot these summary stats ---------------------------#

#### FACS ####

#--- Ia level ---#
xlvls <- c(100, 200, 300, 400)
# nbeta
Ia.ave.sparsity.lasso.dev <- ica_100_8test_glmnet_stat_facs$Ia$lasso.dev$ave.sparsity
Ia.ave.sparsity.lasso.class <- ica_100_8test_glmnet_stat_facs$Ia$lasso.class$ave.sparsity
Ia.ave.sparsity.elnet.dev <- ica_100_8test_glmnet_stat_facs$Ia$elnet.dev$ave.sparsity
Ia.ave.sparsity.elnet.class <- ica_100_8test_glmnet_stat_facs$Ia$elnet.class$ave.sparsity
# ber
Ia.ave.ber.lasso.dev <- 1 - ica_100_8test_glmnet_stat_facs$Ia$lasso.dev$ave.ber
Ia.ave.ber.lasso.class <- 1 - ica_100_8test_glmnet_stat_facs$Ia$lasso.class$ave.ber
Ia.ave.ber.elnet.dev <- 1 - ica_100_8test_glmnet_stat_facs$Ia$elnet.dev$ave.ber
Ia.ave.ber.elnet.class <- 1 - ica_100_8test_glmnet_stat_facs$Ia$elnet.class$ave.ber
# acc
Ia.ave.bacc.lasso.dev <- ica_100_8test_glmnet_stat_facs$Ia$lasso.dev$ave.bacc
Ia.ave.bacc.lasso.class <- ica_100_8test_glmnet_stat_facs$Ia$lasso.class$ave.bacc
Ia.ave.bacc.elnet.dev <- ica_100_8test_glmnet_stat_facs$Ia$elnet.dev$ave.bacc
Ia.ave.bacc.elnet.class <- ica_100_8test_glmnet_stat_facs$Ia$elnet.class$ave.bacc
# kappa
Ia.ave.kappa.lasso.dev <- ica_100_8test_glmnet_stat_facs$Ia$lasso.dev$ave.kappa
Ia.ave.kappa.lasso.class <- ica_100_8test_glmnet_stat_facs$Ia$lasso.class$ave.kappa
Ia.ave.kappa.elnet.dev <- ica_100_8test_glmnet_stat_facs$Ia$elnet.dev$ave.kappa
Ia.ave.kappa.elnet.class <- ica_100_8test_glmnet_stat_facs$Ia$elnet.class$ave.kappa

###
png(filename = "~/projects/icatest/output/test8/facs_Ia.spars.vs.error.png", width=650, height =700, units="px") #  pointsize = 13
## add extra space to right margin of plot within frame
par(mfcol = c(2,1))
par(mar=c(0.5, 6.5, 5, 4) + 0.1)
## Plot first set of data and draw its axis
plot(xlvls, Ia.ave.sparsity.lasso.dev, pch=16, type="b", lty=1, col="blue", axes=FALSE, xlim=c(75,425), ylim=c(25,90), xlab="", ylab="")
points(xlvls, Ia.ave.sparsity.lasso.class, pch=15, type="p", col="blue", ylim=c(25,90), xlab="", ylab="")
points(xlvls, Ia.ave.sparsity.elnet.dev, pch=21, type="b", lty=2, col="blue", ylim=c(25,90), xlab="", ylab="")
points(xlvls, Ia.ave.sparsity.elnet.class, pch=22, type="p", col="blue", ylim=c(25,90), xlab="", ylab="")
mtext(side=3, line=3, cex=1.25, "Sparsity and Predictive performance: Lasso,  Elastic Net / Deviance, Misclassification error")
mtext(side=3, line=1.25, cex=1.1, "Myeloid vs. Non-myeloid cells")
mtext(side=2, line=3.5, cex=1.25, "Total # of active components", col="blue")
axis(side=2, at=seq(30,90,by=10), col="blue", col.axis="blue", las=1)  ## las=1 makes horizontal labels
## Add legend
legend("bottom",legend=c("1-Balanced error rate","Kappa statistics","Overal accuracy","Number of coefficients","Lasso - deviance", "Lasso - missclass.error", "El.Net - deviance", "El.Net - missclass.error"), 
       text.col=c("red","darkorange2","darkred","blue","grey45","grey45","grey45","grey45"), 
       col=c("red","darkorange2","darkred","blue","grey45","grey45","grey45","grey45"), 
       pch=c(16,16,16,16,16,15,21,22), lty=c(0,0,0,0,1,0,2,0), ncol = 2, inset=0.01)
box()

## Plot the second plot
par(mar=c(2.5, 6.5, 0.5, 4) + 0.1)
plot(xlvls, Ia.ave.kappa.lasso.dev, pch=16, type="b", lty=1, col="darkorange2", axes=FALSE, xlim=c(75,425), ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ia.ave.kappa.lasso.class, pch=15, type="p", col="darkorange2", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ia.ave.kappa.elnet.dev, pch=21, type="b", lty=2, col="darkorange2", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ia.ave.kappa.elnet.class, pch=22, type="p", col="darkorange2", ylim=c(0.8,1), xlab="", ylab="")

points(xlvls, Ia.ave.ber.lasso.dev, pch=16, type="b", lty=1, col="red", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ia.ave.ber.lasso.class, pch=15, type="p", col="red", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ia.ave.ber.elnet.dev, pch=21, type="b", lty=2, col="red", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ia.ave.ber.elnet.class, pch=22, type="p", col="red", ylim=c(0.8,1), xlab="", ylab="")

points(xlvls, Ia.ave.bacc.lasso.dev, pch=16, type="b", lty=1, col="darkred", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ia.ave.bacc.lasso.class, pch=15, type="p", col="darkred", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ia.ave.bacc.elnet.dev, pch=21, type="b", lty=2, col="darkred", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ia.ave.bacc.elnet.class, pch=22, type="p", col="darkred", ylim=c(0.8,1), xlab="", ylab="")
## a little farther out (line=4) to make room for labels
mtext(side=2, line=3.5, "Predictive performance statistics", col="orangered3", cex=1.25)
axis(side=2, at=seq(0.8,1,by=0.02), col="orangered3",col.axis="orangered3", las=1)
## Draw the time axis
axis(side=1, at = xlvls, labels = c("VSt log", "VSt RC", "SCt resid", "SCt RC"), cex=1.3, font=2)
box()
dev.off()


#--- Ib level ---#
xlvls <- c(100, 200, 300, 400)
# nbeta
Ib.ave.sparsity.lasso.dev <- ica_100_8test_glmnet_stat_facs$Ib$lasso.dev$ave.sparsity
Ib.ave.sparsity.lasso.class <- ica_100_8test_glmnet_stat_facs$Ib$lasso.class$ave.sparsity
Ib.ave.sparsity.elnet.dev <- ica_100_8test_glmnet_stat_facs$Ib$elnet.dev$ave.sparsity
Ib.ave.sparsity.elnet.class <- ica_100_8test_glmnet_stat_facs$Ib$elnet.class$ave.sparsity
# ber
Ib.ave.ber.lasso.dev <- 1 - ica_100_8test_glmnet_stat_facs$Ib$lasso.dev$ave.ber
Ib.ave.ber.lasso.class <- 1 - ica_100_8test_glmnet_stat_facs$Ib$lasso.class$ave.ber
Ib.ave.ber.elnet.dev <- 1 - ica_100_8test_glmnet_stat_facs$Ib$elnet.dev$ave.ber
Ib.ave.ber.elnet.class <- 1 - ica_100_8test_glmnet_stat_facs$Ib$elnet.class$ave.ber
# acc
Ib.ave.bacc.lasso.dev <- ica_100_8test_glmnet_stat_facs$Ib$lasso.dev$ave.bacc
Ib.ave.bacc.lasso.class <- ica_100_8test_glmnet_stat_facs$Ib$lasso.class$ave.bacc
Ib.ave.bacc.elnet.dev <- ica_100_8test_glmnet_stat_facs$Ib$elnet.dev$ave.bacc
Ib.ave.bacc.elnet.class <- ica_100_8test_glmnet_stat_facs$Ib$elnet.class$ave.bacc
# kappa
Ib.ave.kappa.lasso.dev <- ica_100_8test_glmnet_stat_facs$Ib$lasso.dev$ave.kappa
Ib.ave.kappa.lasso.class <- ica_100_8test_glmnet_stat_facs$Ib$lasso.class$ave.kappa
Ib.ave.kappa.elnet.dev <- ica_100_8test_glmnet_stat_facs$Ib$elnet.dev$ave.kappa
Ib.ave.kappa.elnet.class <- ica_100_8test_glmnet_stat_facs$Ib$elnet.class$ave.kappa

###
png(filename = "~/projects/icatest/output/test8/facs_Ib.spars.vs.error.png", width=650, height =700, units="px") #  pointsize = 13
## add extra space to right margin of plot within frame
par(mfcol = c(2,1))
par(mar=c(0.5, 6.5, 5, 4) + 0.1)
## Plot first set of data and draw its axis
plot(xlvls, Ib.ave.sparsity.lasso.dev, pch=16, type="b", lty=1, col="blue", axes=FALSE, xlim=c(75,425), ylim=c(25,90), xlab="", ylab="")
points(xlvls, Ib.ave.sparsity.lasso.class, pch=15, type="p", col="blue", ylim=c(25,90), xlab="", ylab="")
points(xlvls, Ib.ave.sparsity.elnet.dev, pch=21, type="b", lty=2, col="blue", ylim=c(25,90), xlab="", ylab="")
points(xlvls, Ib.ave.sparsity.elnet.class, pch=22, type="p", col="blue", ylim=c(25,90), xlab="", ylab="")
mtext(side=3, line=3, cex=1.25, "Sparsity and Predictive performance: Lasso,  Elastic Net / Deviance, Misclassification error")
mtext(side=3, line=1.25, cex=1.1, "Myeloid cell types")
mtext(side=2, line=3.5, cex=1.25, "Total # of active components", col="blue")
axis(side=2, at=seq(30,90,by=10), col="blue", col.axis="blue", las=1)  ## las=1 makes horizontal labels
## Add legend
legend("top",legend=c("1-Balanced error rate","Kappa statistics","Overal accuracy","Number of coefficients","Lasso - deviance", "Lasso - missclass.error", "El.Net - deviance", "El.Net - missclass.error"), 
       text.col=c("red","darkorange2","darkred","blue","grey45","grey45","grey45","grey45"), 
       col=c("red","darkorange2","darkred","blue","grey45","grey45","grey45","grey45"), 
       pch=c(16,16,16,16,16,15,21,22), lty=c(0,0,0,0,1,0,2,0), ncol = 2, inset=0.01)
box()

## Plot the second plot
par(mar=c(2.5, 6.5, 0.5, 4) + 0.1)
plot(xlvls, Ib.ave.kappa.lasso.dev, pch=16, type="b", lty=1, col="darkorange2", axes=FALSE, xlim=c(75,425), ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ib.ave.kappa.lasso.class, pch=15, type="p", col="darkorange2", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ib.ave.kappa.elnet.dev, pch=21, type="b", lty=2, col="darkorange2", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ib.ave.kappa.elnet.class, pch=22, type="p", col="darkorange2", ylim=c(0.8,1), xlab="", ylab="")

points(xlvls, Ib.ave.ber.lasso.dev, pch=16, type="b", lty=1, col="red", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ib.ave.ber.lasso.class, pch=15, type="p", col="red", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ib.ave.ber.elnet.dev, pch=21, type="b", lty=2, col="red", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ib.ave.ber.elnet.class, pch=22, type="p", col="red", ylim=c(0.8,1), xlab="", ylab="")

points(xlvls, Ib.ave.bacc.lasso.dev, pch=16, type="b", lty=1, col="darkred", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ib.ave.bacc.lasso.class, pch=15, type="p", col="darkred", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ib.ave.bacc.elnet.dev, pch=21, type="b", lty=2, col="darkred", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ib.ave.bacc.elnet.class, pch=22, type="p", col="darkred", ylim=c(0.8,1), xlab="", ylab="")
## a little farther out (line=4) to make room for labels
mtext(side=2, line=3.5, "Predictive performance statistics", col="orangered3", cex=1.25)
axis(side=2, at=seq(0.8,1,by=0.02), col="orangered3",col.axis="orangered3", las=1)
## Draw the time axis
axis(side=1, at = xlvls, labels = c("VSt log", "VSt RC", "SCt resid", "SCt RC"), cex=1.3, font=2)
box()
dev.off()


#--- Ic level ---#
xlvls <- c(100, 200, 300, 400)
# nbeta
Ic.ave.sparsity.lasso.dev <- ica_100_8test_glmnet_stat_facs$Ic$lasso.dev$ave.sparsity
Ic.ave.sparsity.lasso.class <- ica_100_8test_glmnet_stat_facs$Ic$lasso.class$ave.sparsity
Ic.ave.sparsity.elnet.dev <- ica_100_8test_glmnet_stat_facs$Ic$elnet.dev$ave.sparsity
Ic.ave.sparsity.elnet.class <- ica_100_8test_glmnet_stat_facs$Ic$elnet.class$ave.sparsity
# ber
Ic.ave.ber.lasso.dev <- 1 - ica_100_8test_glmnet_stat_facs$Ic$lasso.dev$ave.ber
Ic.ave.ber.lasso.class <- 1 - ica_100_8test_glmnet_stat_facs$Ic$lasso.class$ave.ber
Ic.ave.ber.elnet.dev <- 1 - ica_100_8test_glmnet_stat_facs$Ic$elnet.dev$ave.ber
Ic.ave.ber.elnet.class <- 1 - ica_100_8test_glmnet_stat_facs$Ic$elnet.class$ave.ber
# acc
Ic.ave.bacc.lasso.dev <- ica_100_8test_glmnet_stat_facs$Ic$lasso.dev$ave.bacc
Ic.ave.bacc.lasso.class <- ica_100_8test_glmnet_stat_facs$Ic$lasso.class$ave.bacc
Ic.ave.bacc.elnet.dev <- ica_100_8test_glmnet_stat_facs$Ic$elnet.dev$ave.bacc
Ic.ave.bacc.elnet.class <- ica_100_8test_glmnet_stat_facs$Ic$elnet.class$ave.bacc
# kappa
Ic.ave.kappa.lasso.dev <- ica_100_8test_glmnet_stat_facs$Ic$lasso.dev$ave.kappa
Ic.ave.kappa.lasso.class <- ica_100_8test_glmnet_stat_facs$Ic$lasso.class$ave.kappa
Ic.ave.kappa.elnet.dev <- ica_100_8test_glmnet_stat_facs$Ic$elnet.dev$ave.kappa
Ic.ave.kappa.elnet.class <- ica_100_8test_glmnet_stat_facs$Ic$elnet.class$ave.kappa

###
png(filename = "~/projects/icatest/output/test8/facs_Ic.spars.vs.error.png", width=650, height =700, units="px") #  pointsize = 13
## add extra space to right margin of plot within frame
par(mfcol = c(2,1))
par(mar=c(0.5, 6.5, 5, 4) + 0.1)
## Plot first set of data and draw its axis
plot(xlvls, Ic.ave.sparsity.lasso.dev, pch=16, type="b", lty=1, col="blue", axes=FALSE, xlim=c(75,425), ylim=c(25,90), xlab="", ylab="")
points(xlvls, Ic.ave.sparsity.lasso.class, pch=15, type="p", col="blue", ylim=c(25,90), xlab="", ylab="")
points(xlvls, Ic.ave.sparsity.elnet.dev, pch=21, type="b", lty=2, col="blue", ylim=c(25,90), xlab="", ylab="")
points(xlvls, Ic.ave.sparsity.elnet.class, pch=22, type="p", col="blue", ylim=c(25,90), xlab="", ylab="")
mtext(side=3, line=3, cex=1.25, "Sparsity and Predictive performance: Lasso,  Elastic Net / Deviance, Misclassification error")
mtext(side=3, line=1.25, cex=1.1, "Macrophages in different tissues")
mtext(side=2, line=3.5, cex=1.25, "Total # of active components", col="blue")
axis(side=2, at=seq(30,90,by=10), col="blue", col.axis="blue", las=1)  ## las=1 makes horizontal labels
## Add legend
legend("top",legend=c("1-Balanced error rate","Kappa statistics","Overal accuracy","Number of coefficients","Lasso - deviance", "Lasso - missclass.error", "El.Net - deviance", "El.Net - missclass.error"), 
       text.col=c("red","darkorange2","darkred","blue","grey45","grey45","grey45","grey45"), 
       col=c("red","darkorange2","darkred","blue","grey45","grey45","grey45","grey45"), 
       pch=c(16,16,16,16,16,15,21,22), lty=c(0,0,0,0,1,0,2,0), ncol = 2, inset=0.01)
box()

## Plot the second plot
par(mar=c(2.5, 6.5, 0.5, 4) + 0.1)
plot(xlvls, Ic.ave.kappa.lasso.dev, pch=16, type="b", lty=1, col="darkorange2", axes=FALSE, xlim=c(75,425), ylim=c(0.775,1), xlab="", ylab="")
points(xlvls, Ic.ave.kappa.lasso.class, pch=15, type="p", col="darkorange2", ylim=c(0.775,1), xlab="", ylab="")
points(xlvls, Ic.ave.kappa.elnet.dev, pch=21, type="b", lty=2, col="darkorange2", ylim=c(0.775,1), xlab="", ylab="")
points(xlvls, Ic.ave.kappa.elnet.class, pch=22, type="p", col="darkorange2", ylim=c(0.775,1), xlab="", ylab="")

points(xlvls, Ic.ave.ber.lasso.dev, pch=16, type="b", lty=1, col="red", ylim=c(0.775,1), xlab="", ylab="")
points(xlvls, Ic.ave.ber.lasso.class, pch=15, type="p", col="red", ylim=c(0.775,1), xlab="", ylab="")
points(xlvls, Ic.ave.ber.elnet.dev, pch=21, type="b", lty=2, col="red", ylim=c(0.775,1), xlab="", ylab="")
points(xlvls, Ic.ave.ber.elnet.class, pch=22, type="p", col="red", ylim=c(0.775,1), xlab="", ylab="")

points(xlvls, Ic.ave.bacc.lasso.dev, pch=16, type="b", lty=1, col="darkred", ylim=c(0.775,1), xlab="", ylab="")
points(xlvls, Ic.ave.bacc.lasso.class, pch=15, type="p", col="darkred", ylim=c(0.775,1), xlab="", ylab="")
points(xlvls, Ic.ave.bacc.elnet.dev, pch=21, type="b", lty=2, col="darkred", ylim=c(0.775,1), xlab="", ylab="")
points(xlvls, Ic.ave.bacc.elnet.class, pch=22, type="p", col="darkred", ylim=c(0.775,1), xlab="", ylab="")
## a little farther out (line=4) to make room for labels
mtext(side=2, line=3.5, "Predictive performance statistics", col="orangered3", cex=1.25)
axis(side=2, at=seq(0.8,1,by=0.02), col="orangered3",col.axis="orangered3", las=1)
## Draw the time axis
axis(side=1, at = xlvls, labels = c("VSt log", "VSt RC", "SCt resid", "SCt RC"), cex=1.3, font=2)
box()
dev.off()


#--- Inc level ---#
xlvls <- c(100, 200, 300, 400)
# nbeta
Inc.ave.sparsity.lasso.dev <- ica_100_8test_glmnet_stat_facs$Inc$lasso.dev$ave.sparsity
Inc.ave.sparsity.lasso.class <- ica_100_8test_glmnet_stat_facs$Inc$lasso.class$ave.sparsity
Inc.ave.sparsity.elnet.dev <- ica_100_8test_glmnet_stat_facs$Inc$elnet.dev$ave.sparsity
Inc.ave.sparsity.elnet.class <- ica_100_8test_glmnet_stat_facs$Inc$elnet.class$ave.sparsity
# ber
Inc.ave.ber.lasso.dev <- 1 - ica_100_8test_glmnet_stat_facs$Inc$lasso.dev$ave.ber
Inc.ave.ber.lasso.class <- 1 - ica_100_8test_glmnet_stat_facs$Inc$lasso.class$ave.ber
Inc.ave.ber.elnet.dev <- 1 - ica_100_8test_glmnet_stat_facs$Inc$elnet.dev$ave.ber
Inc.ave.ber.elnet.class <- 1 - ica_100_8test_glmnet_stat_facs$Inc$elnet.class$ave.ber
# acc
Inc.ave.bacc.lasso.dev <- ica_100_8test_glmnet_stat_facs$Inc$lasso.dev$ave.bacc
Inc.ave.bacc.lasso.class <- ica_100_8test_glmnet_stat_facs$Inc$lasso.class$ave.bacc
Inc.ave.bacc.elnet.dev <- ica_100_8test_glmnet_stat_facs$Inc$elnet.dev$ave.bacc
Inc.ave.bacc.elnet.class <- ica_100_8test_glmnet_stat_facs$Inc$elnet.class$ave.bacc
# kappa
Inc.ave.kappa.lasso.dev <- ica_100_8test_glmnet_stat_facs$Inc$lasso.dev$ave.kappa
Inc.ave.kappa.lasso.class <- ica_100_8test_glmnet_stat_facs$Inc$lasso.class$ave.kappa
Inc.ave.kappa.elnet.dev <- ica_100_8test_glmnet_stat_facs$Inc$elnet.dev$ave.kappa
Inc.ave.kappa.elnet.class <- ica_100_8test_glmnet_stat_facs$Inc$elnet.class$ave.kappa

###
png(filename = "~/projects/icatest/output/test8/facs_Inc.spars.vs.error.png", width=650, height =700, units="px") #  pointsize = 13
## add extra space to right margin of plot within frame
par(mfcol = c(2,1))
par(mar=c(0.5, 6.5, 5, 4) + 0.1)
## Plot first set of data and draw its axis
plot(xlvls, Inc.ave.sparsity.lasso.dev, pch=16, type="b", lty=1, col="blue", axes=FALSE, xlim=c(75,425), ylim=c(10,90), xlab="", ylab="")
points(xlvls, Inc.ave.sparsity.lasso.class, pch=15, type="p", col="blue", xlab="", ylab="")
points(xlvls, Inc.ave.sparsity.elnet.dev, pch=21, type="b", lty=2, col="blue", xlab="", ylab="")
points(xlvls, Inc.ave.sparsity.elnet.class, pch=22, type="p", col="blue", xlab="", ylab="")
mtext(side=3, line=3, cex=1.25, "Sparsity and Predictive performance: Lasso,  Elastic Net / Deviance, Misclassification error")
mtext(side=3, line=1.25, cex=1.1, "Macrophages in different tissues")
mtext(side=2, line=3.5, cex=1.25, "Total # of active components", col="blue")
axis(side=2, at=seq(30,90,by=10), col="blue", col.axis="blue", las=1)  ## las=1 makes horizontal labels
## Add legend
legend("top",legend=c("1-Balanced error rate","Kappa statistics","Overal accuracy","Number of coefficients","Lasso - deviance", "Lasso - missclass.error", "El.Net - deviance", "El.Net - missclass.error"), 
       text.col=c("red","darkorange2","darkred","blue","grey45","grey45","grey45","grey45"), 
       col=c("red","darkorange2","darkred","blue","grey45","grey45","grey45","grey45"), 
       pch=c(16,16,16,16,16,15,21,22), lty=c(0,0,0,0,1,0,2,0), ncol = 2, inset=0.01)
box()

## Plot the second plot
par(mar=c(2.5, 6.5, 0.5, 4) + 0.1)
plot(xlvls, Inc.ave.kappa.lasso.dev, pch=16, type="b", lty=1, col="darkorange2", axes=FALSE, xlim=c(75,425), ylim=c(0.38,1), xlab="", ylab="")
points(xlvls, Inc.ave.kappa.lasso.class, pch=15, type="p", col="darkorange2", xlab="", ylab="")
points(xlvls, Inc.ave.kappa.elnet.dev, pch=21, type="b", lty=2, col="darkorange2", xlab="", ylab="")
points(xlvls, Inc.ave.kappa.elnet.class, pch=22, type="p", col="darkorange2", xlab="", ylab="")

points(xlvls, Inc.ave.ber.lasso.dev, pch=16, type="b", lty=1, col="red", xlab="", ylab="")
points(xlvls, Inc.ave.ber.lasso.class, pch=15, type="p", col="red", xlab="", ylab="")
points(xlvls, Inc.ave.ber.elnet.dev, pch=21, type="b", lty=2, col="red", xlab="", ylab="")
points(xlvls, Inc.ave.ber.elnet.class, pch=22, type="p", col="red", xlab="", ylab="")

points(xlvls, Inc.ave.bacc.lasso.dev, pch=16, type="b", lty=1, col="darkred", xlab="", ylab="")
points(xlvls, Inc.ave.bacc.lasso.class, pch=15, type="p", col="darkred", xlab="", ylab="")
points(xlvls, Inc.ave.bacc.elnet.dev, pch=21, type="b", lty=2, col="darkred", xlab="", ylab="")
points(xlvls, Inc.ave.bacc.elnet.class, pch=22, type="p", col="darkred", xlab="", ylab="")
## a little farther out (line=4) to make room for labels
mtext(side=2, line=3.5, "Predictive performance statistics", col="orangered3", cex=1.25)
axis(side=2, at=seq(0.4,1,by=0.1), col="orangered3",col.axis="orangered3", las=1)
## Draw the time axis
axis(side=1, at = xlvls, labels = c("VSt log", "VSt RC", "SCt resid", "SCt RC"), cex=1.3, font=2)
box()
dev.off()




#### DROP ####

#--- Ia level ---#
xlvls <- c(100, 200, 300, 400)
# nbeta
Ia.ave.sparsity.lasso.dev <- ica_100_8test_glmnet_stat_drop$Ia$lasso.dev$ave.sparsity
Ia.ave.sparsity.lasso.class <- ica_100_8test_glmnet_stat_drop$Ia$lasso.class$ave.sparsity
Ia.ave.sparsity.elnet.dev <- ica_100_8test_glmnet_stat_drop$Ia$elnet.dev$ave.sparsity
Ia.ave.sparsity.elnet.class <- ica_100_8test_glmnet_stat_drop$Ia$elnet.class$ave.sparsity
# ber
Ia.ave.ber.lasso.dev <- 1 - ica_100_8test_glmnet_stat_drop$Ia$lasso.dev$ave.ber
Ia.ave.ber.lasso.class <- 1 - ica_100_8test_glmnet_stat_drop$Ia$lasso.class$ave.ber
Ia.ave.ber.elnet.dev <- 1 - ica_100_8test_glmnet_stat_drop$Ia$elnet.dev$ave.ber
Ia.ave.ber.elnet.class <- 1 - ica_100_8test_glmnet_stat_drop$Ia$elnet.class$ave.ber
# acc
Ia.ave.bacc.lasso.dev <- ica_100_8test_glmnet_stat_drop$Ia$lasso.dev$ave.bacc
Ia.ave.bacc.lasso.class <- ica_100_8test_glmnet_stat_drop$Ia$lasso.class$ave.bacc
Ia.ave.bacc.elnet.dev <- ica_100_8test_glmnet_stat_drop$Ia$elnet.dev$ave.bacc
Ia.ave.bacc.elnet.class <- ica_100_8test_glmnet_stat_drop$Ia$elnet.class$ave.bacc
# kappa
Ia.ave.kappa.lasso.dev <- ica_100_8test_glmnet_stat_drop$Ia$lasso.dev$ave.kappa
Ia.ave.kappa.lasso.class <- ica_100_8test_glmnet_stat_drop$Ia$lasso.class$ave.kappa
Ia.ave.kappa.elnet.dev <- ica_100_8test_glmnet_stat_drop$Ia$elnet.dev$ave.kappa
Ia.ave.kappa.elnet.class <- ica_100_8test_glmnet_stat_drop$Ia$elnet.class$ave.kappa

###
png(filename = "~/projects/icatest/output/test8/drop_Ia.spars.vs.error.png", width=650, height =700, units="px") #  pointsize = 13
## add extra space to right margin of plot within frame
par(mfcol = c(2,1))
par(mar=c(0.5, 6.5, 5, 4) + 0.1)
## Plot first set of data and draw its axis
plot(xlvls, Ia.ave.sparsity.lasso.dev, pch=16, type="b", lty=1, col="blue", axes=FALSE, xlim=c(75,425), ylim=c(25,90), xlab="", ylab="")
points(xlvls, Ia.ave.sparsity.lasso.class, pch=15, type="p", col="blue", ylim=c(25,90), xlab="", ylab="")
points(xlvls, Ia.ave.sparsity.elnet.dev, pch=21, type="b", lty=2, col="blue", ylim=c(25,90), xlab="", ylab="")
points(xlvls, Ia.ave.sparsity.elnet.class, pch=22, type="p", col="blue", ylim=c(25,90), xlab="", ylab="")
mtext(side=3, line=3, cex=1.25, "Sparsity and Predictive performance: Lasso,  Elastic Net / Deviance, Misclassification error")
mtext(side=3, line=1.25, cex=1.1, "Myeloid vs. Non-myeloid cells")
mtext(side=2, line=3.5, cex=1.25, "Total # of active components", col="blue")
axis(side=2, at=seq(30,90,by=10), col="blue", col.axis="blue", las=1)  ## las=1 makes horizontal labels
## Add legend
legend("bottom",legend=c("1-Balanced error rate","Kappa statistics","Overal accuracy","Number of coefficients","Lasso - deviance", "Lasso - missclass.error", "El.Net - deviance", "El.Net - missclass.error"), 
       text.col=c("red","darkorange2","darkred","blue","grey45","grey45","grey45","grey45"), 
       col=c("red","darkorange2","darkred","blue","grey45","grey45","grey45","grey45"), 
       pch=c(16,16,16,16,16,15,21,22), lty=c(0,0,0,0,1,0,2,0), ncol = 2, inset=0.01)
box()

## Plot the second plot
par(mar=c(2.5, 6.5, 0.5, 4) + 0.1)
plot(xlvls, Ia.ave.kappa.lasso.dev, pch=16, type="b", lty=1, col="darkorange2", axes=FALSE, xlim=c(75,425), ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ia.ave.kappa.lasso.class, pch=15, type="p", col="darkorange2", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ia.ave.kappa.elnet.dev, pch=21, type="b", lty=2, col="darkorange2", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ia.ave.kappa.elnet.class, pch=22, type="p", col="darkorange2", ylim=c(0.8,1), xlab="", ylab="")

points(xlvls, Ia.ave.ber.lasso.dev, pch=16, type="b", lty=1, col="red", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ia.ave.ber.lasso.class, pch=15, type="p", col="red", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ia.ave.ber.elnet.dev, pch=21, type="b", lty=2, col="red", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ia.ave.ber.elnet.class, pch=22, type="p", col="red", ylim=c(0.8,1), xlab="", ylab="")

points(xlvls, Ia.ave.bacc.lasso.dev, pch=16, type="b", lty=1, col="darkred", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ia.ave.bacc.lasso.class, pch=15, type="p", col="darkred", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ia.ave.bacc.elnet.dev, pch=21, type="b", lty=2, col="darkred", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ia.ave.bacc.elnet.class, pch=22, type="p", col="darkred", ylim=c(0.8,1), xlab="", ylab="")
## a little farther out (line=4) to make room for labels
mtext(side=2, line=3.5, "Predictive performance statistics", col="orangered3", cex=1.25)
axis(side=2, at=seq(0.8,1,by=0.02), col="orangered3",col.axis="orangered3", las=1)
## Draw the time axis
axis(side=1, at = xlvls, labels = c("VSt log", "VSt RC", "SCt resid", "SCt RC"), cex=1.3, font=2)
box()
dev.off()


#--- Ib level ---#
xlvls <- c(100, 200, 300, 400)
# nbeta
Ib.ave.sparsity.lasso.dev <- ica_100_8test_glmnet_stat_drop$Ib$lasso.dev$ave.sparsity
Ib.ave.sparsity.lasso.class <- ica_100_8test_glmnet_stat_drop$Ib$lasso.class$ave.sparsity
Ib.ave.sparsity.elnet.dev <- ica_100_8test_glmnet_stat_drop$Ib$elnet.dev$ave.sparsity
Ib.ave.sparsity.elnet.class <- ica_100_8test_glmnet_stat_drop$Ib$elnet.class$ave.sparsity
# ber
Ib.ave.ber.lasso.dev <- 1 - ica_100_8test_glmnet_stat_drop$Ib$lasso.dev$ave.ber
Ib.ave.ber.lasso.class <- 1 - ica_100_8test_glmnet_stat_drop$Ib$lasso.class$ave.ber
Ib.ave.ber.elnet.dev <- 1 - ica_100_8test_glmnet_stat_drop$Ib$elnet.dev$ave.ber
Ib.ave.ber.elnet.class <- 1 - ica_100_8test_glmnet_stat_drop$Ib$elnet.class$ave.ber
# acc
Ib.ave.bacc.lasso.dev <- ica_100_8test_glmnet_stat_drop$Ib$lasso.dev$ave.bacc
Ib.ave.bacc.lasso.class <- ica_100_8test_glmnet_stat_drop$Ib$lasso.class$ave.bacc
Ib.ave.bacc.elnet.dev <- ica_100_8test_glmnet_stat_drop$Ib$elnet.dev$ave.bacc
Ib.ave.bacc.elnet.class <- ica_100_8test_glmnet_stat_drop$Ib$elnet.class$ave.bacc
# kappa
Ib.ave.kappa.lasso.dev <- ica_100_8test_glmnet_stat_drop$Ib$lasso.dev$ave.kappa
Ib.ave.kappa.lasso.class <- ica_100_8test_glmnet_stat_drop$Ib$lasso.class$ave.kappa
Ib.ave.kappa.elnet.dev <- ica_100_8test_glmnet_stat_drop$Ib$elnet.dev$ave.kappa
Ib.ave.kappa.elnet.class <- ica_100_8test_glmnet_stat_drop$Ib$elnet.class$ave.kappa

###
png(filename = "~/projects/icatest/output/test8/drop_Ib.spars.vs.error.png", width=650, height =700, units="px") #  pointsize = 13
## add extra space to right margin of plot within frame
par(mfcol = c(2,1))
par(mar=c(0.5, 6.5, 5, 4) + 0.1)
## Plot first set of data and draw its axis
plot(xlvls, Ib.ave.sparsity.lasso.dev, pch=16, type="b", lty=1, col="blue", axes=FALSE, xlim=c(75,425), ylim=c(25,90), xlab="", ylab="")
points(xlvls, Ib.ave.sparsity.lasso.class, pch=15, type="p", col="blue", ylim=c(25,90), xlab="", ylab="")
points(xlvls, Ib.ave.sparsity.elnet.dev, pch=21, type="b", lty=2, col="blue", ylim=c(25,90), xlab="", ylab="")
points(xlvls, Ib.ave.sparsity.elnet.class, pch=22, type="p", col="blue", ylim=c(25,90), xlab="", ylab="")
mtext(side=3, line=3, cex=1.25, "Sparsity and Predictive performance: Lasso,  Elastic Net / Deviance, Misclassification error")
mtext(side=3, line=1.25, cex=1.1, "Myeloid cell types")
mtext(side=2, line=3.5, cex=1.25, "Total # of active components", col="blue")
axis(side=2, at=seq(30,90,by=10), col="blue", col.axis="blue", las=1)  ## las=1 makes horizontal labels
## Add legend
legend("bottom",legend=c("1-Balanced error rate","Kappa statistics","Overal accuracy","Number of coefficients","Lasso - deviance", "Lasso - missclass.error", "El.Net - deviance", "El.Net - missclass.error"), 
       text.col=c("red","darkorange2","darkred","blue","grey45","grey45","grey45","grey45"), 
       col=c("red","darkorange2","darkred","blue","grey45","grey45","grey45","grey45"), 
       pch=c(16,16,16,16,16,15,21,22), lty=c(0,0,0,0,1,0,2,0), ncol = 2, inset=0.01)
box()

## Plot the second plot
par(mar=c(2.5, 6.5, 0.5, 4) + 0.1)
plot(xlvls, Ib.ave.kappa.lasso.dev, pch=16, type="b", lty=1, col="darkorange2", axes=FALSE, xlim=c(75,425), ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ib.ave.kappa.lasso.class, pch=15, type="p", col="darkorange2", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ib.ave.kappa.elnet.dev, pch=21, type="b", lty=2, col="darkorange2", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ib.ave.kappa.elnet.class, pch=22, type="p", col="darkorange2", ylim=c(0.8,1), xlab="", ylab="")

points(xlvls, Ib.ave.ber.lasso.dev, pch=16, type="b", lty=1, col="red", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ib.ave.ber.lasso.class, pch=15, type="p", col="red", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ib.ave.ber.elnet.dev, pch=21, type="b", lty=2, col="red", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ib.ave.ber.elnet.class, pch=22, type="p", col="red", ylim=c(0.8,1), xlab="", ylab="")

points(xlvls, Ib.ave.bacc.lasso.dev, pch=16, type="b", lty=1, col="darkred", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ib.ave.bacc.lasso.class, pch=15, type="p", col="darkred", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ib.ave.bacc.elnet.dev, pch=21, type="b", lty=2, col="darkred", ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ib.ave.bacc.elnet.class, pch=22, type="p", col="darkred", ylim=c(0.8,1), xlab="", ylab="")
## a little farther out (line=4) to make room for labels
mtext(side=2, line=3.5, "Predictive performance statistics", col="orangered3", cex=1.25)
axis(side=2, at=seq(0.8,1,by=0.02), col="orangered3",col.axis="orangered3", las=1)
## Draw the time axis
axis(side=1, at = xlvls, labels = c("VSt log", "VSt RC", "SCt resid", "SCt RC"), cex=1.3, font=2)
box()
dev.off()


#--- Ic level ---#
xlvls <- c(100, 200, 300, 400)
# nbeta
Ic.ave.sparsity.lasso.dev <- ica_100_8test_glmnet_stat_drop$Ic$lasso.dev$ave.sparsity
Ic.ave.sparsity.lasso.class <- ica_100_8test_glmnet_stat_drop$Ic$lasso.class$ave.sparsity
Ic.ave.sparsity.elnet.dev <- ica_100_8test_glmnet_stat_drop$Ic$elnet.dev$ave.sparsity
Ic.ave.sparsity.elnet.class <- ica_100_8test_glmnet_stat_drop$Ic$elnet.class$ave.sparsity
# ber
Ic.ave.ber.lasso.dev <- 1 - ica_100_8test_glmnet_stat_drop$Ic$lasso.dev$ave.ber
Ic.ave.ber.lasso.class <- 1 - ica_100_8test_glmnet_stat_drop$Ic$lasso.class$ave.ber
Ic.ave.ber.elnet.dev <- 1 - ica_100_8test_glmnet_stat_drop$Ic$elnet.dev$ave.ber
Ic.ave.ber.elnet.class <- 1 - ica_100_8test_glmnet_stat_drop$Ic$elnet.class$ave.ber
# acc
Ic.ave.bacc.lasso.dev <- ica_100_8test_glmnet_stat_drop$Ic$lasso.dev$ave.bacc
Ic.ave.bacc.lasso.class <- ica_100_8test_glmnet_stat_drop$Ic$lasso.class$ave.bacc
Ic.ave.bacc.elnet.dev <- ica_100_8test_glmnet_stat_drop$Ic$elnet.dev$ave.bacc
Ic.ave.bacc.elnet.class <- ica_100_8test_glmnet_stat_drop$Ic$elnet.class$ave.bacc
# kappa
Ic.ave.kappa.lasso.dev <- ica_100_8test_glmnet_stat_drop$Ic$lasso.dev$ave.kappa
Ic.ave.kappa.lasso.class <- ica_100_8test_glmnet_stat_drop$Ic$lasso.class$ave.kappa
Ic.ave.kappa.elnet.dev <- ica_100_8test_glmnet_stat_drop$Ic$elnet.dev$ave.kappa
Ic.ave.kappa.elnet.class <- ica_100_8test_glmnet_stat_drop$Ic$elnet.class$ave.kappa

###
png(filename = "~/projects/icatest/output/test8/drop_Ic.spars.vs.error.png", width=650, height =700, units="px") #  pointsize = 13
## add extra space to right margin of plot within frame
par(mfcol = c(2,1))
par(mar=c(0.5, 6.5, 5, 4) + 0.1)
## Plot first set of data and draw its axis
plot(xlvls, Ic.ave.sparsity.lasso.dev, pch=16, type="b", lty=1, col="blue", axes=FALSE, xlim=c(75,425), ylim=c(25,90), xlab="", ylab="")
points(xlvls, Ic.ave.sparsity.lasso.class, pch=15, type="p", col="blue", ylim=c(25,90), xlab="", ylab="")
points(xlvls, Ic.ave.sparsity.elnet.dev, pch=21, type="b", lty=2, col="blue", ylim=c(25,90), xlab="", ylab="")
points(xlvls, Ic.ave.sparsity.elnet.class, pch=22, type="p", col="blue", ylim=c(25,90), xlab="", ylab="")
mtext(side=3, line=3, cex=1.25, "Sparsity and Predictive performance: Lasso,  Elastic Net / Deviance, Misclassification error")
mtext(side=3, line=1.25, cex=1.1, "Macrophages in different tissues")
mtext(side=2, line=3.5, cex=1.25, "Total # of active components", col="blue")
axis(side=2, at=seq(30,90,by=10), col="blue", col.axis="blue", las=1)  ## las=1 makes horizontal labels
## Add legend
legend("bottom",legend=c("1-Balanced error rate","Kappa statistics","Overal accuracy","Number of coefficients","Lasso - deviance", "Lasso - missclass.error", "El.Net - deviance", "El.Net - missclass.error"), 
       text.col=c("red","darkorange2","darkred","blue","grey45","grey45","grey45","grey45"), 
       col=c("red","darkorange2","darkred","blue","grey45","grey45","grey45","grey45"), 
       pch=c(16,16,16,16,16,15,21,22), lty=c(0,0,0,0,1,0,2,0), ncol = 2, inset=0.01)
box()

## Plot the second plot
par(mar=c(2.5, 6.5, 0.5, 4) + 0.1)
plot(xlvls, Ic.ave.kappa.lasso.dev, pch=16, type="b", lty=1, col="darkorange2", axes=FALSE, xlim=c(75,425), ylim=c(0.8,1), xlab="", ylab="")
points(xlvls, Ic.ave.kappa.lasso.class, pch=15, type="p", col="darkorange2", xlab="", ylab="")
points(xlvls, Ic.ave.kappa.elnet.dev, pch=21, type="b", lty=2, col="darkorange2", xlab="", ylab="")
points(xlvls, Ic.ave.kappa.elnet.class, pch=22, type="p", col="darkorange2", xlab="", ylab="")

points(xlvls, Ic.ave.ber.lasso.dev, pch=16, type="b", lty=1, col="red", xlab="", ylab="")
points(xlvls, Ic.ave.ber.lasso.class, pch=15, type="p", col="red", xlab="", ylab="")
points(xlvls, Ic.ave.ber.elnet.dev, pch=21, type="b", lty=2, col="red", xlab="", ylab="")
points(xlvls, Ic.ave.ber.elnet.class, pch=22, type="p", col="red", xlab="", ylab="")

points(xlvls, Ic.ave.bacc.lasso.dev, pch=16, type="b", lty=1, col="darkred", xlab="", ylab="")
points(xlvls, Ic.ave.bacc.lasso.class, pch=15, type="p", col="darkred", xlab="", ylab="")
points(xlvls, Ic.ave.bacc.elnet.dev, pch=21, type="b", lty=2, col="darkred", xlab="", ylab="")
points(xlvls, Ic.ave.bacc.elnet.class, pch=22, type="p", col="darkred", xlab="", ylab="")
## a little farther out (line=4) to make room for labels
mtext(side=2, line=3.5, "Predictive performance statistics", col="orangered3", cex=1.25)
axis(side=2, at=seq(0.8,1,by=0.02), col="orangered3",col.axis="orangered3", las=1)
## Draw the time axis
axis(side=1, at = xlvls, labels = c("VSt log", "VSt RC", "SCt resid", "SCt RC"), cex=1.3, font=2)
box()
dev.off()


#--- Inc level ---#
xlvls <- c(100, 200, 300, 400)
# nbeta
Inc.ave.sparsity.lasso.dev <- ica_100_8test_glmnet_stat_drop$Inc$lasso.dev$ave.sparsity
Inc.ave.sparsity.lasso.class <- ica_100_8test_glmnet_stat_drop$Inc$lasso.class$ave.sparsity
Inc.ave.sparsity.elnet.dev <- ica_100_8test_glmnet_stat_drop$Inc$elnet.dev$ave.sparsity
Inc.ave.sparsity.elnet.class <- ica_100_8test_glmnet_stat_drop$Inc$elnet.class$ave.sparsity
# ber
Inc.ave.ber.lasso.dev <- 1 - ica_100_8test_glmnet_stat_drop$Inc$lasso.dev$ave.ber
Inc.ave.ber.lasso.class <- 1 - ica_100_8test_glmnet_stat_drop$Inc$lasso.class$ave.ber
Inc.ave.ber.elnet.dev <- 1 - ica_100_8test_glmnet_stat_drop$Inc$elnet.dev$ave.ber
Inc.ave.ber.elnet.class <- 1 - ica_100_8test_glmnet_stat_drop$Inc$elnet.class$ave.ber
# acc
Inc.ave.bacc.lasso.dev <- ica_100_8test_glmnet_stat_drop$Inc$lasso.dev$ave.bacc
Inc.ave.bacc.lasso.class <- ica_100_8test_glmnet_stat_drop$Inc$lasso.class$ave.bacc
Inc.ave.bacc.elnet.dev <- ica_100_8test_glmnet_stat_drop$Inc$elnet.dev$ave.bacc
Inc.ave.bacc.elnet.class <- ica_100_8test_glmnet_stat_drop$Inc$elnet.class$ave.bacc
# kappa
Inc.ave.kappa.lasso.dev <- ica_100_8test_glmnet_stat_drop$Inc$lasso.dev$ave.kappa
Inc.ave.kappa.lasso.class <- ica_100_8test_glmnet_stat_drop$Inc$lasso.class$ave.kappa
Inc.ave.kappa.elnet.dev <- ica_100_8test_glmnet_stat_drop$Inc$elnet.dev$ave.kappa
Inc.ave.kappa.elnet.class <- ica_100_8test_glmnet_stat_drop$Inc$elnet.class$ave.kappa

###
png(filename = "~/projects/icatest/output/test8/drop_Inc.spars.vs.error.png", width=650, height =700, units="px") #  pointsize = 13
## add extra space to right margin of plot within frame
par(mfcol = c(2,1))
par(mar=c(0.5, 6.5, 5, 4) + 0.1)
## Plot first set of data and draw its axis
plot(xlvls, Inc.ave.sparsity.lasso.dev, pch=16, type="b", lty=1, col="blue", axes=FALSE, xlim=c(75,425), ylim=c(10,95), xlab="", ylab="")
points(xlvls, Inc.ave.sparsity.lasso.class, pch=15, type="p", col="blue", xlab="", ylab="")
points(xlvls, Inc.ave.sparsity.elnet.dev, pch=21, type="b", lty=2, col="blue", xlab="", ylab="")
points(xlvls, Inc.ave.sparsity.elnet.class, pch=22, type="p", col="blue", xlab="", ylab="")
mtext(side=3, line=3, cex=1.25, "Sparsity and Predictive performance: Lasso,  Elastic Net / Deviance, Misclassification error")
mtext(side=3, line=1.25, cex=1.1, "Macrophages in different tissues")
mtext(side=2, line=3.5, cex=1.25, "Total # of active components", col="blue")
axis(side=2, at=seq(30,90,by=10), col="blue", col.axis="blue", las=1)  ## las=1 makes horizontal labels
## Add legend
legend("bottom",legend=c("1-Balanced error rate","Kappa statistics","Overal accuracy","Number of coefficients","Lasso - deviance", "Lasso - missclass.error", "El.Net - deviance", "El.Net - missclass.error"), 
       text.col=c("red","darkorange2","darkred","blue","grey45","grey45","grey45","grey45"), 
       col=c("red","darkorange2","darkred","blue","grey45","grey45","grey45","grey45"), 
       pch=c(16,16,16,16,16,15,21,22), lty=c(0,0,0,0,1,0,2,0), ncol = 2, inset=0.01)
box()

## Plot the second plot
par(mar=c(2.5, 6.5, 0.5, 4) + 0.1)
plot(xlvls, Inc.ave.kappa.lasso.dev, pch=16, type="b", lty=1, col="darkorange2", axes=FALSE, xlim=c(75,425), ylim=c(0.38,1), xlab="", ylab="")
points(xlvls, Inc.ave.kappa.lasso.class, pch=15, type="p", col="darkorange2", xlab="", ylab="")
points(xlvls, Inc.ave.kappa.elnet.dev, pch=21, type="b", lty=2, col="darkorange2", xlab="", ylab="")
points(xlvls, Inc.ave.kappa.elnet.class, pch=22, type="p", col="darkorange2", xlab="", ylab="")

points(xlvls, Inc.ave.ber.lasso.dev, pch=16, type="b", lty=1, col="red", xlab="", ylab="")
points(xlvls, Inc.ave.ber.lasso.class, pch=15, type="p", col="red", xlab="", ylab="")
points(xlvls, Inc.ave.ber.elnet.dev, pch=21, type="b", lty=2, col="red", xlab="", ylab="")
points(xlvls, Inc.ave.ber.elnet.class, pch=22, type="p", col="red", xlab="", ylab="")

points(xlvls, Inc.ave.bacc.lasso.dev, pch=16, type="b", lty=1, col="darkred", xlab="", ylab="")
points(xlvls, Inc.ave.bacc.lasso.class, pch=15, type="p", col="darkred", xlab="", ylab="")
points(xlvls, Inc.ave.bacc.elnet.dev, pch=21, type="b", lty=2, col="darkred", xlab="", ylab="")
points(xlvls, Inc.ave.bacc.elnet.class, pch=22, type="p", col="darkred", xlab="", ylab="")
## a little farther out (line=4) to make room for labels
mtext(side=2, line=3.5, "Predictive performance statistics", col="orangered3", cex=1.25)
axis(side=2, at=seq(0.4,1,by=0.1), col="orangered3",col.axis="orangered3", las=1)
## Draw the time axis
axis(side=1, at = xlvls, labels = c("VSt log", "VSt RC", "SCt resid", "SCt RC"), cex=1.3, font=2)
box()
dev.off()



# The rest of the plotting e.g. cell type breakdown (nbetas vs. accuracy) can be added (below)
# generally RC data do not perform very well (higher nbetas and lower prediction accuracy)
# also VSt log data performs generally slightly better than SCt resid data (specially on UMI data?!)

### FACS ###
#-- Ib --#
library("RColorBrewer")
display.brewer.all(colorblindFriendly = TRUE)
cols <- c(brewer.pal(12, "Paired")[c(9,10)], brewer.pal(12, "Paired")[c(5,6,8)], brewer.pal(12, "Paired")[11], brewer.pal(8, "Set2")[1])

png(filename = "~/projects/icatest/output/test8/facs_Ib.spars.vs.acc.breakdown.png", width = 1200, height = 900, units = "px") #  pointsize = 13
par(mfrow=c(2,2))
for(mod in c("lasso.dev", "lasso.class", "elnet.dev", "elnet.class")) {
  par(mar=c(3,5,5,2) + 0.1)
  bp <- barplot(t(ica_100_8test_glmnet_stat_facs$Ib[[mod]]$nbetas), beside=TRUE, col = cols, ylim=c(0,52), space=c(0, 1.75))
  text(bp, t(ica_100_8test_glmnet_stat_facs$Ib[[mod]]$nbetas)+1.25, labels = round(t(ica_100_8test_glmnet_stat_facs$Ib[[mod]]$bacc), digits=2), cex=0.7) 
  mtext(side=3, line=3, cex=1.3, paste0("Sparsity and Predictive performance " , "(", mod, ")"))
  mtext(side=3, line=1.25, cex=1.1, "Myeloid cell types")
  mtext(side=2, line=3, cex=1.1, "Number of active coefficients (components)", col="black")
  legend("topright", legend=rownames(t(ica_100_8test_glmnet_stat_facs$Ib[[mod]]$nbetas)), fill=cols) # bty="n",
}
dev.off()

#-- Ic --#
cols <- c(brewer.pal(6, "Set2")[5], brewer.pal(7, "Reds")[5], brewer.pal(6, "Set2")[2], brewer.pal(8, "Set2")[c(7,8)], brewer.pal(12, "Paired")[c(1,8)], brewer.pal(12, "Paired")[11])

png(filename = "~/projects/icatest/output/test8/facs_Ic.spars.vs.acc.breakdown.png", width = 1200, height = 900, units = "px") #  pointsize = 13
par(mfrow=c(2,2))
for(mod in c("lasso.dev", "lasso.class", "elnet.dev", "elnet.class")) {
  par(mar=c(3,5,5,2) + 0.1)
  bp <- barplot(t(ica_100_8test_glmnet_stat_facs$Ic[[mod]]$nbetas), beside=TRUE, col = cols, ylim=c(0,50), space=c(0, 1.75))
  text(bp, t(ica_100_8test_glmnet_stat_facs$Ic[[mod]]$nbetas)+1.25, labels = round(t(ica_100_8test_glmnet_stat_facs$Ic[[mod]]$bacc), digits=2), cex=0.7) 
  mtext(side=3, line=3, cex=1.3, paste0("Sparsity and Predictive performance " , "(", mod, ")"))
  mtext(side=3, line=1.25, cex=1.1, "Macrophages in different tissues")
  mtext(side=2, line=3, cex=1.1, "Number of active coefficients (components)", col="black")
  legend("topleft", legend=rownames(t(ica_100_8test_glmnet_stat_facs$Ic[[mod]]$nbetas)), fill=cols, inset=0.01) # bty="n",
}
dev.off()


### DROP ###
#-- Ib --#
cols <- c(brewer.pal(12, "Paired")[c(9,10)], brewer.pal(12, "Paired")[c(7,5,6)], brewer.pal(8, "Set2")[c(6,8)], brewer.pal(8, "Dark2")[4])

png(filename = "~/projects/icatest/output/test8/drop_Ib.spars.vs.acc.breakdown.png", width = 1200, height = 900, units = "px") #  pointsize = 13
par(mfrow=c(2,2))
for(mod in c("lasso.dev", "lasso.class", "elnet.dev", "elnet.class")) {
  par(mar=c(3,5,5,2) + 0.1)
  bp <- barplot(t(ica_100_8test_glmnet_stat_drop$Ib[[mod]]$nbetas), beside=TRUE, col = cols, ylim=c(0,74), space=c(0, 1.75))
  text(bp, t(ica_100_8test_glmnet_stat_drop$Ib[[mod]]$nbetas)+1.25, labels = round(t(ica_100_8test_glmnet_stat_drop$Ib[[mod]]$bacc), digits=2), cex=0.7) 
  mtext(side=3, line=3, cex=1.3, paste0("Sparsity and Predictive performance " , "(", mod, ")"))
  mtext(side=3, line=1.25, cex=1.1, "Myeloid cell types")
  mtext(side=2, line=3, cex=1.1, "Number of active coefficients (components)", col="black")
  legend("topleft", legend=rownames(t(ica_100_8test_glmnet_stat_drop$Ib[[mod]]$nbetas)), fill=cols, inset=0.01) # bty="n",
}
dev.off()

#-- Ic --#
cols <- c(brewer.pal(6, "Set2")[5], brewer.pal(7, "Reds")[5], brewer.pal(6, "Purples")[4], brewer.pal(6, "Set2")[2], brewer.pal(6, "Set2")[4], brewer.pal(7, "Set2")[7])

png(filename = "~/projects/icatest/output/test8/drop_Ic.spars.vs.acc.breakdown.png", width = 1100, height = 900, units = "px") #  pointsize = 13
par(mfrow=c(2,2))
for(mod in c("lasso.dev", "lasso.class", "elnet.dev", "elnet.class")) {
  par(mar=c(3,5,5,2) + 0.1)
  bp <- barplot(t(ica_100_8test_glmnet_stat_drop$Ic[[mod]]$nbetas), beside=TRUE, col = cols, ylim=c(0,58), space=c(0, 1.75))
  text(bp, t(ica_100_8test_glmnet_stat_drop$Ic[[mod]]$nbetas)+1.25, labels = round(t(ica_100_8test_glmnet_stat_drop$Ic[[mod]]$bacc), digits=2), cex=0.7) 
  mtext(side=3, line=3, cex=1.3, paste0("Sparsity and Predictive performance " , "(", mod, ")"))
  mtext(side=3, line=1.25, cex=1.1, "Macrophages in different tissues")
  mtext(side=2, line=3, cex=1.1, "Number of active coefficients (components)", col="black")
  legend("topleft", legend=rownames(t(ica_100_8test_glmnet_stat_drop$Ic[[mod]]$nbetas)), fill=cols, inset=0.01) # bty="n",
}
dev.off()





###############################################################################
#####                   S matrix side of things - GENES                   #####
###############################################################################

rm( list = setdiff(ls(), "ica_100_8test") )

##---------------------------------------------------------------------------##
## FACS ##
##---------------------------------------------------------------------------##
ica_100_8test_S_facs <- list()
for(i in 1:4){
  S <- ica_100_8test[[i]]$S %>% as.data.frame()
  colnames(S) <- as.character(1:100)
  S <- S %>% rownames_to_column(var = "Symbol") %>% as_tibble()
  
  S.lst <- list(all=list(), sd2=list(), sd3=list(), sd4=list())
  for(j in 1:100) { 
    j <- as.character(j)
    S.lst$all[[j]] <- S[, c("Symbol", j)]
    names(S.lst$all[[j]]) <- c("Symbol", "S")
  }
  
  S.lst$sd2 <- purrr::map(S.lst$all, .f = dplyr::filter, abs(S) >= 2)
  S.lst$sd3 <- purrr::map(S.lst$all, .f = dplyr::filter, abs(S) >= 3)
  S.lst$sd4 <- purrr::map(S.lst$all, .f = dplyr::filter, abs(S) >= 4)
  
  ica_100_8test_S_facs[[i]] <- S.lst
}
names(ica_100_8test_S_facs) <- names(ica_100_8test)[1:4]
saveRDS(ica_100_8test_S_facs, "~/projects/icatest/output/test8/ica_100_8test_S_facs.rds")


### Histograms of module sizes ###
png(filename = "~/projects/icatest/output/test8/figures/facs_module.sizes.png", width = 900, height = 800, units = "px") #  pointsize = 13
par(mfrow=c(4,3))
for(mod in c("facs.vst.log", "facs.vst.rc", "facs.sct.res", "facs.sct.cnt")) {
  i <- case_when(mod == "facs.vst.log" ~ "VSt  log", 
                 mod == "facs.vst.rc" ~ "VSt  RC", 
                 mod == "facs.sct.res" ~ "SCt  res", 
                 mod == "facs.sct.cnt" ~ "SCt  RC")
  
  for(sd in c("sd2", "sd3", "sd4")) {
    j <- case_when(sd == "sd2" ~ "2 sd", 
                   sd == "sd3" ~ "3 sd", 
                   sd == "sd4" ~ "4 sd")
    
    uniq <- purrr::map(ica_100_8test_S_facs[[mod]][[sd]], .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length() 
    reus <- round(sum(sapply(ica_100_8test_S_facs[[mod]][[sd]], nrow)) / purrr::map(ica_100_8test_S_facs[[mod]][[sd]], .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length(), 1)
    
    hist(sapply(ica_100_8test_S_facs[[mod]][[sd]], nrow), breaks = seq(0, 490, by = 10), col = "tomato1" , xlab = "Module size", cex.lab = 1.2, main =paste("100 modules @", j), cex.main = 1.3, ylim = c(0,40))
    mtext(side=3, line=0.05, cex=1, paste("-", i, "-"))
    text(400, 20, paste0(uniq, "  (", reus, "x", ")"), cex = 1.25, col = "blue", font = 2)
  }
}
dev.off()


### Barplot of genes re-use by modules ###
png(filename = "~/projects/icatest/output/test8/figures/facs_module.gene.reuse.png", width = 900, height = 800, units = "px") #  pointsize = 13
par(mfrow=c(4,3))
for(mod in c("facs.vst.log", "facs.vst.rc", "facs.sct.res", "facs.sct.cnt")) {
  i <- case_when(mod == "facs.vst.log" ~ "VSt  log", 
                 mod == "facs.vst.rc" ~ "VSt  RC", 
                 mod == "facs.sct.res" ~ "SCt  res", 
                 mod == "facs.sct.cnt" ~ "SCt  RC")
  
  for(sd in c("sd2", "sd3", "sd4")) {
    j <- case_when(sd == "sd2" ~ "2 sd", 
                   sd == "sd3" ~ "3 sd", 
                   sd == "sd4" ~ "4 sd")
    
    uniq <- purrr::map(ica_100_8test_S_facs[[mod]][[sd]], .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length() 
    reus <- round(sum(sapply(ica_100_8test_S_facs[[mod]][[sd]], nrow)) / purrr::map(ica_100_8test_S_facs[[mod]][[sd]], .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length(), 1)
    
    purrr::map(ica_100_8test_S_facs[[mod]][[sd]], .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% table() %>% sort(decreasing = T) %>% table() %>% 
      barplot(xlab = "number of modules", ylab = "number of unique genes", cex.lab = 1.2, main = paste("100 modules @", j), cex.main = 1.3, ylim = c(0,1500), xlim = c(1,55), col = "seagreen")
    mtext(side=3, line=0.05, cex=1, paste("-", i, "-"))
    text(25, 800, paste0(uniq, "  (", reus, "x", ")"), cex = 1.25, col = "blue", font = 2)
  }
}
dev.off()


### Histogram of module weights ###
png(filename = "~/projects/icatest/output/test8/figures/facs_module.weight.png", width = 900, height = 800, units = "px") #  pointsize = 13
par(mfrow=c(4,3))
for(mod in c("facs.vst.log", "facs.vst.rc", "facs.sct.res", "facs.sct.cnt")) {
  i <- case_when(mod == "facs.vst.log" ~ "VSt  log", 
                 mod == "facs.vst.rc" ~ "VSt  RC", 
                 mod == "facs.sct.res" ~ "SCt  res", 
                 mod == "facs.sct.cnt" ~ "SCt  RC")
  
  for(sd in c("sd2", "sd3", "sd4")) {
    j <- case_when(sd == "sd2" ~ "2 sd", 
                   sd == "sd3" ~ "3 sd", 
                   sd == "sd4" ~ "4 sd")
    
    uniq <- purrr::map(ica_100_8test_S_facs[[mod]][[sd]], .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length() 
    reus <- round(sum(sapply(ica_100_8test_S_facs[[mod]][[sd]], nrow)) / purrr::map(ica_100_8test_S_facs[[mod]][[sd]], .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length(), 1)
    
    purrr::map(ica_100_8test_S_facs[[mod]][[sd]], .f = dplyr::select, S) %>% purrr::map(.f = function(x) x^2) %>% purrr::map(.f = sum) %>% unlist() %>% (function(x) x / 7500) %>% 
      hist(breaks = seq(0, 1, by = 0.025), ylim = c(0,50), col = "thistle" , xlab = "Module weight (proportion of total comp. weight)", cex.lab = 1.2, main = paste("100 modules @", j), cex.main = 1.3)
    mtext(side=3, line=0.05, cex=1, paste("-", i, "-"))
    text(0.25, 25, paste0(uniq, "  (", reus, "x", ")"), cex = 1.25, col = "blue", font = 2)
  }
}
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














##---------------------------------------------------------------------------##
## DROP ##
##---------------------------------------------------------------------------##
ica_100_8test_S_drop <- list()
for(i in 1:4){
  S <- ica_100_8test[[i+4]]$S %>% as.data.frame()
  colnames(S) <- as.character(1:100)
  S <- S %>% rownames_to_column(var = "Symbol") %>% as_tibble()
  
  S.lst <- list(all=list(), sd2=list(), sd3=list(), sd4=list())
  for(j in 1:100) { 
    j <- as.character(j)
    S.lst$all[[j]] <- S[, c("Symbol", j)]
    names(S.lst$all[[j]]) <- c("Symbol", "S")
  }
  
  S.lst$sd2 <- purrr::map(S.lst$all, .f = dplyr::filter, abs(S) >= 2)
  S.lst$sd3 <- purrr::map(S.lst$all, .f = dplyr::filter, abs(S) >= 3)
  S.lst$sd4 <- purrr::map(S.lst$all, .f = dplyr::filter, abs(S) >= 4)
  
  ica_100_8test_S_drop[[i]] <- S.lst
}
names(ica_100_8test_S_drop) <- names(ica_100_8test)[5:8]
saveRDS(ica_100_8test_S_drop, "~/projects/icatest/output/test8/ica_100_8test_S_drop.rds")


### Histograms of module sizes ###
png(filename = "~/projects/icatest/output/test8/figures/drop_module.sizes.png", width = 900, height = 800, units = "px") #  pointsize = 13
par(mfrow=c(4,3))
for(mod in c("drop.vst.log", "drop.vst.rc", "drop.sct.res", "drop.sct.cnt")) {
  i <- case_when(mod == "drop.vst.log" ~ "VSt  log", 
                 mod == "drop.vst.rc" ~ "VSt  RC", 
                 mod == "drop.sct.res" ~ "SCt  res", 
                 mod == "drop.sct.cnt" ~ "SCt  RC")
  
  for(sd in c("sd2", "sd3", "sd4")) {
    j <- case_when(sd == "sd2" ~ "2 sd", 
                   sd == "sd3" ~ "3 sd", 
                   sd == "sd4" ~ "4 sd")
    
    uniq <- purrr::map(ica_100_8test_S_drop[[mod]][[sd]], .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length() 
    reus <- round(sum(sapply(ica_100_8test_S_drop[[mod]][[sd]], nrow)) / purrr::map(ica_100_8test_S_drop[[mod]][[sd]], .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length(), 1)
    
    hist(sapply(ica_100_8test_S_drop[[mod]][[sd]], nrow), breaks = seq(0, 490, by = 10), col = "tomato1" , xlab = "Module size", cex.lab = 1.2, main =paste("100 modules @", j), cex.main = 1.3, ylim = c(0,40))
    mtext(side=3, line=0.05, cex=1, paste("-", i, "-"))
    text(400, 20, paste0(uniq, "  (", reus, "x", ")"), cex = 1.25, col = "blue", font = 2)
  }
}
dev.off()


### Barplot of genes re-use by modules ###
png(filename = "~/projects/icatest/output/test8/figures/drop_module.gene.reuse.png", width = 900, height = 800, units = "px") #  pointsize = 13
par(mfrow=c(4,3))
for(mod in c("drop.vst.log", "drop.vst.rc", "drop.sct.res", "drop.sct.cnt")) {
  i <- case_when(mod == "drop.vst.log" ~ "VSt  log", 
                 mod == "drop.vst.rc" ~ "VSt  RC", 
                 mod == "drop.sct.res" ~ "SCt  res", 
                 mod == "drop.sct.cnt" ~ "SCt  RC")
  
  for(sd in c("sd2", "sd3", "sd4")) {
    j <- case_when(sd == "sd2" ~ "2 sd", 
                   sd == "sd3" ~ "3 sd", 
                   sd == "sd4" ~ "4 sd")
    
    uniq <- purrr::map(ica_100_8test_S_drop[[mod]][[sd]], .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length() 
    reus <- round(sum(sapply(ica_100_8test_S_drop[[mod]][[sd]], nrow)) / purrr::map(ica_100_8test_S_drop[[mod]][[sd]], .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length(), 1)
    
    purrr::map(ica_100_8test_S_drop[[mod]][[sd]], .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% table() %>% sort(decreasing = T) %>% table() %>% 
      barplot(xlab = "number of modules", ylab = "number of unique genes", cex.lab = 1.2, main = paste("100 modules @", j), cex.main = 1.3, ylim = c(0,1500), xlim = c(1,55), col = "seagreen")
    mtext(side=3, line=0.05, cex=1, paste("-", i, "-"))
    text(25, 800, paste0(uniq, "  (", reus, "x", ")"), cex = 1.25, col = "blue", font = 2)
  }
}
dev.off()


### Histogram of module weights ###
png(filename = "~/projects/icatest/output/test8/figures/drop_module.weight.png", width = 900, height = 800, units = "px") #  pointsize = 13
par(mfrow=c(4,3))
for(mod in c("drop.vst.log", "drop.vst.rc", "drop.sct.res", "drop.sct.cnt")) {
  i <- case_when(mod == "drop.vst.log" ~ "VSt  log", 
                 mod == "drop.vst.rc" ~ "VSt  RC", 
                 mod == "drop.sct.res" ~ "SCt  res", 
                 mod == "drop.sct.cnt" ~ "SCt  RC")
  
  for(sd in c("sd2", "sd3", "sd4")) {
    j <- case_when(sd == "sd2" ~ "2 sd", 
                   sd == "sd3" ~ "3 sd", 
                   sd == "sd4" ~ "4 sd")
    
    uniq <- purrr::map(ica_100_8test_S_drop[[mod]][[sd]], .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length() 
    reus <- round(sum(sapply(ica_100_8test_S_drop[[mod]][[sd]], nrow)) / purrr::map(ica_100_8test_S_drop[[mod]][[sd]], .f = dplyr::select, Symbol) %>% flatten() %>% unlist() %>% unique() %>% length(), 1)
    
    purrr::map(ica_100_8test_S_drop[[mod]][[sd]], .f = dplyr::select, S) %>% purrr::map(.f = function(x) x^2) %>% purrr::map(.f = sum) %>% unlist() %>% (function(x) x / 7500) %>% 
      hist(breaks = seq(0, 1, by = 0.025), ylim = c(0,50), col = "thistle" , xlab = "Module weight (proportion of total comp. weight)", cex.lab = 1.2, main = paste("100 modules @", j), cex.main = 1.3)
    mtext(side=3, line=0.05, cex=1, paste("-", i, "-"))
    text(0.25, 25, paste0(uniq, "  (", reus, "x", ")"), cex = 1.25, col = "blue", font = 2)
  }
}
dev.off()




#-------------------------------------------------------------------------------------------------#

# Test expression of genes which are part of key segregating modules accross cell types (Seurat)
# try one of key modules composition at each level Ia, b, c, nc; all for lasso-deviance & elnemt-classerror = 8 plots for facs

ica_100_8test_A_facs <- readRDS("~/projects/icatest/output/test8/ica_100_8test_A_facs.rds")
ica_100_8test_glmnet_stat_facs <- readRDS("~/projects/icatest/output/test8/ica_100_8test_glmnet_stat_facs.rds")
ica_100_8test_S_facs <- readRDS("~/projects/icatest/output/test8/ica_100_8test_S_facs.rds")


ica_100_8test_glmnet_stat_facs$Ia$lasso.dev$Betas[[1]] %>% print(n=100) # 2- (microglia strong), 50- (pericite intermed), 42-,94+ (both mac strong, 42 bm mac strong)
ica_100_8test_glmnet_stat_facs$Ia$elnet.class$Betas[[1]] %>% print(n=100)

ica_100_8test_glmnet_stat_facs$Ib$lasso.dev$Betas[[1]] %>% print(n=100) # 94- mac & mo opposite dir!
ica_100_8test_glmnet_stat_facs$Ib$elnet.class$Betas[[1]] %>% print(n=100)

ica_100_8test_glmnet_stat_facs$Ic$lasso.dev$Betas[[1]] %>% print(n=100) # 64 spl mac
ica_100_8test_glmnet_stat_facs$Ic$elnet.class$Betas[[1]] %>% print(n=100)

ica_100_8test_glmnet_stat_facs$Inc$lasso.dev$Betas[[1]] %>% print(n=100) # 21 kupfer, global mac opposite slightly
ica_100_8test_glmnet_stat_facs$Inc$elnet.class$Betas[[1]] %>% print(n=100)

facs.exp.7500_vst <- readRDS("~/projects/icatest/data/facs.exp.7500_vst.rds")
keeps <- ica_100_8test_A_facs$facs.vst.log %>% pull(var = cellnames)
facs.vst.log <- subset(facs.exp.7500_vst, cells = keeps)

myelo.names <- c("basophil", "brain pericyte", "classical monocyte", "granulocyte", "Kupffer cell", "macrophage", "microglial cell", "monocyte")

celltype <- facs.vst.log@meta.data$cell_ontology_class
cell.meta <- ica_100_8test_A_facs$facs.vst.log %>% dplyr::select(cellnames, celltype, celltype_my, celltype_my_type, celltype_my_type_tissue)
all.equal(celltype, cell.meta$celltype)

celltype.Ia 
celltype.Ib 
celltype.Ic
celltype.Inc


### 2 ###--------------------------------------------------
# module 2 - Ia level
modulegenes <- ica_100_8test_S_facs$facs.vst.log$sd4[[2]] %>% arrange(S) %>% pull(Symbol)

Idents(facs.vst.log) <- factor(cell.meta$celltype_my, levels = c("non myeloid", "myeloid"))
png(filename = "~/projects/icatest/output/test8/figures/dotplot.mod2.Ia.png", width = length(modulegenes)*12.5+200, height = length(levels(Idents(facs.vst.log)))*60+100, units = "px") #  pointsize = 13
DotPlot(facs.vst.log, features = rev(modulegenes), cols = c("palegreen", "blue")) + 
  RotatedAxis() + 
  FontSize(x.text = 9, y.text = 11, x.title = FALSE, y.title = FALSE, main = 12) + 
  ggtitle("Module 2") + 
  theme(legend.text=element_text(size=10), legend.title=element_text(size=9))
dev.off()

# module 2 - Ib level
Idents(facs.vst.log) <- cell.meta$celltype_my_type
png(filename = "~/projects/icatest/output/test8/figures/dotplot.mod2.Ib.png", width = length(modulegenes)*12.5+200, height = length(levels(Idents(facs.vst.log)))*40+100, units = "px")
DotPlot(facs.vst.log, features = rev(modulegenes), cols = c("palegreen", "blue")) + 
  RotatedAxis() + 
  FontSize(x.text = 9, y.text = 11, x.title = FALSE, y.title = FALSE, main = 12) + 
  ggtitle("Module 2") + 
  theme(legend.text=element_text(size=10), legend.title=element_text(size=9))
dev.off()

# module 2 - Ic level
Idents(facs.vst.log) <- cell.meta$celltype_my_type_tissue
png(filename = "~/projects/icatest/output/test8/figures/dotplot.mod2.Ic.png", width = length(modulegenes)*12.5+200, height = length(levels(Idents(facs.vst.log)))*40+100, units = "px")
DotPlot(facs.vst.log, features = rev(modulegenes), cols = c("palegreen", "blue")) + 
  RotatedAxis() + 
  FontSize(x.text = 9, y.text = 11, x.title = FALSE, y.title = FALSE, main = 12) + 
  ggtitle("Module 2") + 
  theme(legend.text=element_text(size=10), legend.title=element_text(size=9))
dev.off()

ica_100_8test_S_facs$facs.vst.log$sd4[[2]] %>% arrange(S) %>% print(n=100)


### 94 ###--------------------------------------------------
# module 94 - Ia level
modulegenes <- ica_100_8test_S_facs$facs.vst.log$sd4[[94]] %>% arrange(S) %>% pull(Symbol)

Idents(facs.vst.log) <- factor(cell.meta$celltype_my, levels = c("non myeloid", "myeloid"))
png(filename = "~/projects/icatest/output/test8/figures/dotplot.mod94.Ia.png", width = length(modulegenes)*12.5+200, height = length(levels(Idents(facs.vst.log)))*60+100, units = "px")
DotPlot(facs.vst.log, features = rev(modulegenes), cols = c("palegreen", "blue")) + 
  RotatedAxis() + 
  FontSize(x.text = 9, y.text = 11, x.title = FALSE, y.title = FALSE, main = 12) + 
  ggtitle("Module 94") + 
  theme(legend.text=element_text(size=10), legend.title=element_text(size=9))
dev.off()

# module 94 - Ib level
Idents(facs.vst.log) <- cell.meta$celltype_my_type
png(filename = "~/projects/icatest/output/test8/figures/dotplot.mod94.Ib.png", width = length(modulegenes)*12.5+200, height = length(levels(Idents(facs.vst.log)))*40+100, units = "px")
DotPlot(facs.vst.log, features = rev(modulegenes), cols = c("palegreen", "blue")) + 
  RotatedAxis() + 
  FontSize(x.text = 9, y.text = 11, x.title = FALSE, y.title = FALSE, main = 12) + 
  ggtitle("Module 94") + 
  theme(legend.text=element_text(size=10), legend.title=element_text(size=9))
dev.off()

# module 94 - Ic level
Idents(facs.vst.log) <- cell.meta$celltype_my_type_tissue
png(filename = "~/projects/icatest/output/test8/figures/dotplot.mod94.Ic.png", width = length(modulegenes)*12.5+200, height = length(levels(Idents(facs.vst.log)))*40+100, units = "px")
DotPlot(facs.vst.log, features = rev(modulegenes), cols = c("palegreen", "blue")) + 
  RotatedAxis() + 
  FontSize(x.text = 9, y.text = 11, x.title = FALSE, y.title = FALSE, main = 12) + 
  ggtitle("Module 94") + 
  theme(legend.text=element_text(size=10), legend.title=element_text(size=9))
dev.off()

ica_100_8test_S_facs$facs.vst.log$sd4[[94]] %>% arrange(S) %>% print(n=100)


### 42 ###--------------------------------------------------
# module 42 - Ia level
modulegenes <- ica_100_8test_S_facs$facs.vst.log$sd4[[42]] %>% arrange(S) %>% pull(Symbol)

Idents(facs.vst.log) <- factor(cell.meta$celltype_my, levels = c("non myeloid", "myeloid"))
png(filename = "~/projects/icatest/output/test8/figures/dotplot.mod42.Ia.png", width = length(modulegenes)*12.5+200, height = length(levels(Idents(facs.vst.log)))*60+100, units = "px")
DotPlot(facs.vst.log, features = rev(modulegenes), cols = c("palegreen", "blue")) + 
  RotatedAxis() + 
  FontSize(x.text = 9, y.text = 11, x.title = FALSE, y.title = FALSE, main = 12) + 
  ggtitle("Module 42") + 
  theme(legend.text=element_text(size=10), legend.title=element_text(size=9))
dev.off()

# module 42 - Ib level
Idents(facs.vst.log) <- cell.meta$celltype_my_type
png(filename = "~/projects/icatest/output/test8/figures/dotplot.mod42.Ib.png", width = length(modulegenes)*12.5+200, height = length(levels(Idents(facs.vst.log)))*40+100, units = "px")
DotPlot(facs.vst.log, features = rev(modulegenes), cols = c("palegreen", "blue")) + 
  RotatedAxis() + 
  FontSize(x.text = 9, y.text = 11, x.title = FALSE, y.title = FALSE, main = 12) + 
  ggtitle("Module 42") + 
  theme(legend.text=element_text(size=10), legend.title=element_text(size=9))
dev.off()

# module 42 - Ic level
Idents(facs.vst.log) <- cell.meta$celltype_my_type_tissue
png(filename = "~/projects/icatest/output/test8/figures/dotplot.mod42.Ic.png", width = length(modulegenes)*12.5+200, height = length(levels(Idents(facs.vst.log)))*40+100, units = "px")
DotPlot(facs.vst.log, features = rev(modulegenes), cols = c("palegreen", "blue")) + 
  RotatedAxis() + 
  FontSize(x.text = 9, y.text = 11, x.title = FALSE, y.title = FALSE, main = 12) + 
  ggtitle("Module 42") + 
  theme(legend.text=element_text(size=10), legend.title=element_text(size=9))
dev.off()

ica_100_8test_S_facs$facs.vst.log$sd4[[42]] %>% arrange(S) %>% print(n=100)


### 64 ###--------------------------------------------------
# module 64 - Ia level
modulegenes <- ica_100_8test_S_facs$facs.vst.log$sd4[[64]] %>% arrange(S) %>% pull(Symbol)

Idents(facs.vst.log) <- factor(cell.meta$celltype_my, levels = c("non myeloid", "myeloid"))
png(filename = "~/projects/icatest/output/test8/figures/dotplot.mod64.Ia.png", width = length(modulegenes)*12.5+200, height = length(levels(Idents(facs.vst.log)))*55+100, units = "px")
DotPlot(facs.vst.log, features = rev(modulegenes), cols = c("palegreen", "blue")) + 
  RotatedAxis() + 
  FontSize(x.text = 9, y.text = 11, x.title = FALSE, y.title = FALSE, main = 12) + 
  ggtitle("Module 64") + 
  theme(legend.text=element_text(size=10), legend.title=element_text(size=9))
dev.off()

# module 64 - Ib level
Idents(facs.vst.log) <- cell.meta$celltype_my_type
png(filename = "~/projects/icatest/output/test8/figures/dotplot.mod64.Ib.png", width = length(modulegenes)*12.5+200, height = length(levels(Idents(facs.vst.log)))*40+100, units = "px")
DotPlot(facs.vst.log, features = rev(modulegenes), cols = c("palegreen", "blue")) + 
  RotatedAxis() + 
  FontSize(x.text = 9, y.text = 11, x.title = FALSE, y.title = FALSE, main = 12) + 
  ggtitle("Module 64") + 
  theme(legend.text=element_text(size=10), legend.title=element_text(size=9))
dev.off()

# module 64 - Ic level
Idents(facs.vst.log) <- cell.meta$celltype_my_type_tissue
png(filename = "~/projects/icatest/output/test8/figures/dotplot.mod64.Ic.png", width = length(modulegenes)*12.5+200, height = length(levels(Idents(facs.vst.log)))*40+100, units = "px")
DotPlot(facs.vst.log, features = rev(modulegenes), cols = c("palegreen", "blue")) + 
  RotatedAxis() + 
  FontSize(x.text = 9, y.text = 11, x.title = FALSE, y.title = FALSE, main = 12) + 
  ggtitle("Module 64") + 
  theme(legend.text=element_text(size=10), legend.title=element_text(size=9))
dev.off()

ica_100_8test_S_facs$facs.vst.log$sd4[[64]] %>% arrange(S) %>% print(n=110)


### 21 ###--------------------------------------------------
# module 21 - Ia level
modulegenes <- ica_100_8test_S_facs$facs.vst.log$sd4[[21]] %>% arrange(S) %>% pull(Symbol)

Idents(facs.vst.log) <- factor(cell.meta$celltype_my, levels = c("non myeloid", "myeloid"))
png(filename = "~/projects/icatest/output/test8/figures/dotplot.mod21.Ia.png", width = length(modulegenes)*12.5+200, height = length(levels(Idents(facs.vst.log)))*55+100, units = "px")
DotPlot(facs.vst.log, features = rev(modulegenes), cols = c("palegreen", "blue")) + 
  RotatedAxis() + 
  FontSize(x.text = 9, y.text = 11, x.title = FALSE, y.title = FALSE, main = 12) + 
  ggtitle("Module 21") + 
  theme(legend.text=element_text(size=10), legend.title=element_text(size=9))
dev.off()

# module 21 - Ib level
Idents(facs.vst.log) <- cell.meta$celltype_my_type
png(filename = "~/projects/icatest/output/test8/figures/dotplot.mod21.Ib.png", width = length(modulegenes)*12.5+200, height = length(levels(Idents(facs.vst.log)))*40+100, units = "px")
DotPlot(facs.vst.log, features = rev(modulegenes), cols = c("palegreen", "blue")) + 
  RotatedAxis() + 
  FontSize(x.text = 9, y.text = 11, x.title = FALSE, y.title = FALSE, main = 12) + 
  ggtitle("Module 21") + 
  theme(legend.text=element_text(size=10), legend.title=element_text(size=9))
dev.off()

# module 21 - Ic level
Idents(facs.vst.log) <- cell.meta$celltype_my_type_tissue
png(filename = "~/projects/icatest/output/test8/figures/dotplot.mod21.Ic.png", width = length(modulegenes)*12.5+200, height = length(levels(Idents(facs.vst.log)))*40+100, units = "px")
DotPlot(facs.vst.log, features = rev(modulegenes), cols = c("palegreen", "blue")) + 
  RotatedAxis() + 
  FontSize(x.text = 9, y.text = 11, x.title = FALSE, y.title = FALSE, main = 12) + 
  ggtitle("Module 21") + 
  theme(legend.text=element_text(size=10), legend.title=element_text(size=9))
dev.off()

ica_100_8test_S_facs$facs.vst.log$sd4[[21]] %>% arrange(S) %>% print(n=110)














library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genes = modulegenes
genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = genes, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
cat(genes$HGNC.symbol, sep = ",")
cat(genes$HGNC.symbol, sep = "\n")


library("msigdbr")




















# pCMF approach
#-----------------------------------------------------------------------------#
devtools::install_github("gdurif/pCMF", subdir="pkg", ref="prod_no_omp") # could not figure out "prod" option 
library(pCMF)









# mitochondrial genes 
# MGI symbol:
# mt-Tf
# mt-Rnr1
# mt-Tv
# mt-Rnr2
# mt-Tl1
# mt-Nd1
# mt-Ti
# mt-Tq
#	mt-Tm
# mt-Nd2
# mt-Tw
#	mt-Ta
# mt-Tn
# mt-Tc
# mt-Ty
# mt-Co1
# mt-Ts1
























