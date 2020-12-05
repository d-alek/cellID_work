# modules_to_TFs.R

library("tidyverse")
library("dplyr")
library("magrittr")
library("stringr")
library("ggplot2")
library("RColorBrewer")
library("SingleCellExperiment") # SCE
library("here")

# Extract the list of human transcription factors from manually currated collection of 1620
# (The Human Transcription Factors: Cell 172, February 8, 2018)
human.tf <- read_tsv(file = here("data", "external", "human_TF.txt"))

# Download mouse-human homology matches from MGI resources
# (Human and Mouse Homology Classes with Sequence information: http://www.informatics.jax.org/faq/ORTH_dload.shtml)
homologenes <- read_tsv(file = here("data", "external", "HOM_MouseHumanSequence.rpt.txt"))

# Extract matching mouse TF gene symbols + entrezids (1303 of them)
human.tf.homologene.id <- homologenes %>% 
  dplyr::filter(`Common Organism Name` == "human" & Symbol %in% human.tf$Name) %>% 
  dplyr::select(`HomoloGene ID`, Symbol)

mouse.tf <- homologenes %>% 
  dplyr::filter(`Common Organism Name` == "mouse, laboratory" & `HomoloGene ID` %in% human.tf.homologene.id$`HomoloGene ID`)
saveRDS(mouse.tf, "data/external/mouse_TF.rds")

# Check for overlap in entire TM dataset, in spikeF and sumF normalised highly-variable gene subsets
#-- full set 20698 genes
sumF.sc <- readRDS(file = here("data", "processed", "tm.smartseq2.sumF.123.rds")) # 20698

length(intersect(mouse.tf$Symbol, rowData(sumF.sc)$Symbol)) # 1257
length(intersect(mouse.tf$`EntrezGene ID`, rowData(sumF.sc)$entrez)) # 1257
length(intersect(mouse.tf$Symbol, rowData(sumF.sc)$ID)) # 1202 (old gene symbols from TM)
# save them:
tf.expressed <- intersect(mouse.tf$Symbol, rowData(sumF.sc)$Symbol)
saveRDS(tf.expressed, "data/processed/tf_expressed.rds")
# tf.expressed <- readRDS(file = here("data", "processed", "tf_expressed.rds"))

#-- 4.5-6.5k HVG subsets
spikeF.sc <- readRDS(file = here("data", "processed", "tm.smartseq2.spikeF.hvg.123.icaDims.rds")) # 6342
sumF.sc <- readRDS(file = here("data", "processed", "tm.smartseq2.sumF.hvg.123.icaDims.rds")) # 4748

# spikeF:
tf.spikeF.hv <- intersect(mouse.tf$Symbol, rowData(spikeF.sc)$Symbol) # 259 (selecting 30% hv genes, retains only 20% of TFs) (same # if entrezid used for matching!)
saveRDS(tf.spikeF.hv, "data/processed/tf_spikeF_hv.rds")

# sumF:
tf.sumF.hv <- intersect(mouse.tf$Symbol, rowData(sumF.sc)$Symbol) # 176 (selecting 23% hv genes, retains only 14% of TFs) (same # if entrezid used for matching!)
saveRDS(tf.sumF.hv, "data/processed/tf_sumF_hv.rds")



# Try to explain A scores by TFs (A score ~ Tfs Lasso / El.Net models, similar to the Celltype ~ A score)



