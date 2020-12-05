# Spec score calculation per gene, per cell type 
# Measuring cell identity in noisy biological systems. 2011. Birnbaum KD, Kussel E. Nucleic Acids Research. 39. 9093-917. 

# Spec score can be defined as an average amount of info. the gene's expression level provides about the cell's identity 
# it is related to Mutual Information (where MI = average Spec score over all cell types under analysis)
# but Spec scores provide more specific information for each cell type - gene combination

# tables and plots are made for exploration of genes with top Spec score for each cell type or 
# exploration of cell types in which a particular gene has highest Spec scores
# ..these can be used to filter ICA module genes for increased cell-specificity

# some potential improvements:
# if testing of various numbers of bins for discretisation of gene expression is needed, maybe better make this into a function
# would be good to develop permutation-based null distribution to asses significance of Spec scores.. but computationally expensive

library("tidyverse")
library("dplyr")
library("magrittr")
library("stringr")
library("SingleCellExperiment") # SCE
library("scran")
library("scater")
library("ggplot2")
library("RColorBrewer")
library("viridis")
library("pheatmap")

sumF.sc <- readRDS("data/processed/tm.smartseq2.sumF.123.rds")
sumF.A.df <- readRDS("results/sumF_A_c_meta.rds")

sum(rowData(sumF.sc)$bio_sum >= 0.00969) # 7501
sum(rowData(sumF.sc)$is_HVG == T) # 4748

# pick 7505 most variable genes (for which to calculate Spec score) and fit in tibble along with cell type annotations (cells x genes)
keeps <- rowData(sumF.sc)$bio_sum >= 0.00969 | rowData(sumF.sc)$is_HVG == T # 7505
sumF.sc.7500 <- sumF.sc[keeps, ]
tm.7500 <- t(logcounts(sumF.sc.7500)) %>% as.matrix() %>% as.data.frame() %>% as_tibble()
tm.7500 <- bind_cols(sumF.A.df[, c("cell", "celltype", "tissue", "celltype_tissue")], tm.7500) 

# the idea is to average gene expression per each 10 cells of a particular cell type, so that Spec estimates are not too noisy
dplyr::count(tm.7500, celltype) %>% print(n = Inf)
dplyr::count(tm.7500, celltype_tissue) %>% print(n = Inf)

tm.7500.summ <- dplyr::count(tm.7500, celltype) %>% 
  mutate(n.10 = as.integer(round(n/10)))

# create cell grouping vector per cell type and add it to tibble for group_by()
set.seed(123)
cell_group <- map(.x = pull(tm.7500.summ, n), .f = function(x) sample(rep(1:round(x/10), len=x)))
names(cell_group) <- pull(tm.7500.summ, celltype)

cell_group_vec <- flatten_int(cell_group)
tm.7500 <- arrange(tm.7500, celltype) %>% add_column(cell_group=cell_group_vec, .after="celltype")

# get average gene expression for each pool of 10 cells of the same type 
tm.7500.binned <- tm.7500 %>% 
  group_by(celltype, cell_group) %>% 
  summarise(across('0610005C13Rik':'Zyx', mean))

tm.7500.binned %>% ungroup() %>% pull(Wfdc21) %>% density() %>% plot()
tm.7500.binned %>% ungroup() %>% pull(Spi1) %>% hist(breaks=50)

# bin gene expression into equally spaced bins
tm.7500.bins5 <- tm.7500.binned %>% 
  ungroup() %>% 
  mutate(across('0610005C13Rik':'Zyx', cut, breaks=5, labels=c("1","2","3","4","5")), .keep="unused") 
  #mutate(across('0610005C13Rik':'Zyx', fct_collapse, "5"=c("5","6")), .keep="unused")

tm.7500.bins5 %>% pull(Wfdc21) %>% table()
tm.7500.bins5 %>% pull(Spib) %>% table()

# maybe try various numbers of bins.. easier to make it into a function later 
tm.7500.bins10
tm.7500.bins20
tm.7500.bins40

# calculate p(gen.exp.level | cell.type) ---------------------------------------------------------#
tm.7500.bins5.levelp <- tm.7500.bins5 %>% 
  group_by(celltype) %>% 
  summarise(across('0610005C13Rik':'Zyx', fct_count, sort=FALSE, prop=TRUE)) %>% 
  ungroup()

# fix the weird data structure this created
fix_colnames <- colnames(tm.7500.bins5)[3:7507]

bins5.levelp <- tm.7500.bins5.levelp["celltype"]
for (i in fix_colnames) {
  adds <- tm.7500.bins5.levelp[i] %>% flatten_dfc() %>% select(p)
  names(adds) <- i
  bins5.levelp <- bind_cols(bins5.levelp, adds)
}
bins5.levelp <- bins5.levelp %>% add_column(levels = as_factor(rep(c(1,2,3,4,5), 76)), .after="celltype")  

# add p(cell.type)
celltype_p <- tm.7500.bins5 %>% 
  group_by(celltype) %>% 
  summarise(celltype.p=n()/nrow(tm.7500.bins5)) %>% 
  pull("celltype.p")
celltype_p <- rep(celltype_p, each=5)
bins5.levelp <- bins5.levelp %>% add_column(celltype.p = celltype_p, .after="celltype")

# calculate p(gen.exp.level | cell.type) * p(cell.type) ------------------------------------------#
bins5.levelp.p <- modify2(.x=select(bins5.levelp, '0610005C13Rik':'Zyx'), .y=select(bins5.levelp, 'celltype.p'), .f= ~ .x * .y)
bins5.levelp.p <- bins5.levelp.p %>% 
  add_column(celltype = pull(bins5.levelp, 'celltype'), levels = pull(bins5.levelp, 'levels'), .before="0610005C13Rik")

# calculate sum[p(gen.exp.level | cell.type) * p(cell.type)] over cells - denominator 
bins5.levelp.p.denom <- bins5.levelp.p %>% 
  group_by(levels) %>% 
  summarise(across('0610005C13Rik':'Zyx', sum))

# calculate p(cell.type | gen.exp.level) ---------------------------------------------------------#
bins5.typep <- modify2(.x=select(bins5.levelp.p, '0610005C13Rik':'Zyx'), .y=select(bins5.levelp.p.denom, '0610005C13Rik':'Zyx'), .f= ~ .x / .y)
bins5.typep <- bins5.typep %>% 
  add_column(celltype = pull(bins5.levelp.p, 'celltype'), levels = pull(bins5.levelp.p, 'levels'), .before="0610005C13Rik")
bins5.typep[is.na(bins5.typep)] <- 0

# calculate Log p(cell.type | gen.exp.level) 
ntypes <- length(unique(bins5.typep$celltype))
bins5.typep.log <- bins5.typep %>% modify_at(.at=c(-1,-2), .f=log, base=ntypes)
bins5.typep.log[bins5.typep.log == -Inf] <- 0

# calculate gen.exp.level info -------------------------------------------------------------------#
# I(level) = 1 + Sum[p(cell.type | gen.exp.level) * Log p(cell.type | gen.exp.level)] over cell types
bins5.typep.prod <- modify2(.x=select(bins5.typep, '0610005C13Rik':'Zyx'), .y=select(bins5.typep.log, '0610005C13Rik':'Zyx'), .f= ~ .x * .y)
bins5.typep.prod <- bins5.typep.prod %>% 
  add_column(celltype = pull(bins5.typep, 'celltype'), levels = pull(bins5.typep, 'levels'), .before="0610005C13Rik")
i.level <- bins5.typep.prod %>% 
  group_by(levels) %>% 
  summarise(across('0610005C13Rik':'Zyx', sum)) %>% 
  modify_at(.at=-1, .f= ~ . + 1)

# calculate Informativeness of gene for each cell type - Spec ------------------------------------#
# Spec = Sum[p(gen.exp.level | cell.type) * I(level)] over gene exp levels
bins5.levelp.I.prod <- modify2(.x=select(bins5.levelp, '0610005C13Rik':'Zyx'), .y=select(i.level, '0610005C13Rik':'Zyx'), .f= ~ .x * .y)
bins5.levelp.I.prod <- bins5.levelp.I.prod %>% 
  add_column(celltype = pull(bins5.levelp, 'celltype'), levels = pull(bins5.levelp, 'levels'), .before="0610005C13Rik")
Spec <- bins5.levelp.I.prod %>% 
  group_by(celltype) %>% 
  summarise(across('0610005C13Rik':'Zyx', sum))
Spec %>% print(n=100)  
write_tsv(Spec, "results/spec_bin5_123.txt")
# Spec <- read_tsv("results/spec_bin5_123.txt")

Spec %>% select(celltype, Spi1) %>% arrange(desc(Spi1)) %>% print(n=77)
Spec %>% select(celltype, Il7r) %>% arrange(desc(Il7r)) %>% print(n=77)
Spec %>% select(celltype, Gata1) %>% arrange(desc(Gata1)) %>% print(n=77)
Spec %>% select(celltype, Gata2) %>% arrange(desc(Gata2)) %>% print(n=77)
Spec %>% select(celltype, Mrc1) %>% arrange(desc(Mrc1)) %>% print(n=77) 
Spec %>% select(celltype, Aadac) %>% arrange(desc(Aadac)) %>% print(n=77)
Spec %>% select(celltype, Elf3) %>% arrange(desc(Elf3)) %>% print(n=77)

# transposed version for easier exploration 
Spec.t <- Spec %>% select(-1) %>% as.matrix() %>% t() # matrix
colnames(Spec.t) <- Spec %>% pull(celltype)

Spec.tt <- Spec.t %>% as.data.frame(row.names=names(Spec)[-1]) %>% as_tibble(rownames="genes") 

# check some cell types
names(Spec.tt)
Spec.tt %>% select(genes, macrophage) %>% arrange(desc(macrophage)) %>% print(n=100)
Spec.tt %>% select(genes, monocyte, 'classical monocyte') %>% arrange(desc(monocyte)) %>% print(n=100)
Spec.tt %>% select(genes, monocyte, 'classical monocyte') %>% arrange(desc(`classical monocyte`)) %>% print(n=100)
Spec.tt %>% select(genes, granulocyte) %>% arrange(desc(granulocyte)) %>% print(n=100)
Spec.tt %>% select(genes, basophil) %>% arrange(desc(basophil)) %>% print(n=100)
Spec.tt %>% select(genes, hepatocyte) %>% arrange(desc(hepatocyte)) %>% print(n=100)
Spec.tt %>% select(genes, `epithelial cell`) %>% arrange(desc(`epithelial cell`)) %>% print(n=100)

# plots: Spec scores per cell type 
spec_cell_long <- Spec.tt %>% 
  gather(-genes, key='celltype', value='spec') %>% 
  group_by(celltype) %>% 
  mutate(spec.rank=rank(-spec, ties.method="random")) %>% 
  ungroup()

# ggplot(spec_cell_long, aes(x=spec.rank, y=spec)) +
#   geom_line(aes(group=factor(celltype), colour=factor(celltype))) +
#   theme_bw() +
#   theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), plot.title=element_text(hjust=0.5), legend.position="none") + 
#   labs(title="Spec scores per cell type, ranked across genes", x="Spec score rank", y="Spec score") + 
#   geom_label(data=spec_cell_long %>% filter((spec>=0.5 & spec.rank>=800) | (spec<=0.3 & spec.rank<=175)), aes(label=celltype), check_overlap = TRUE)

ggplot(spec_cell_long, aes(x=spec.rank, y=spec)) +
  geom_line(aes(group=factor(celltype), colour=factor(celltype))) +
  theme_bw() +
  theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), plot.title=element_text(hjust=0.5), legend.position="none") + 
  labs(title="Spec scores per cell type, ranked", x="Spec score rank", y="Spec score") + 
  annotate(geom="label", label="hepatocyte", x=1500, y=0.45, label.size=0.1)
ggsave("figures/spec_percell_bin5.pdf", width = 6, height = 4.5)

# plots: Spec scores per gene
spec_gene_long <- Spec %>% 
  gather(-celltype, key='gene', value='spec') %>% 
  group_by(gene) %>% 
  mutate(spec.rank=rank(-spec, ties.method="random")) %>% 
  ungroup()

ggplot(spec_gene_long, aes(x=spec.rank, y=spec)) +
  geom_line(aes(group=factor(gene)), size=0.05, colour="royalblue", alpha=1/10) +
  theme_bw() + 
  theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), plot.title=element_text(hjust=0.5)) + 
  labs(title="Spec scores per gene, ranked across cell types", x="Spec score rank", y="Spec score")
ggsave("figures/spec_pergene_bin5.pdf", width = 6, height = 4.5)

# heatplot of all
# first trim off genes with max Spec <= 0.5 ~ 5050 genes
rowmax <- Spec %>% summarise(across('0610005C13Rik':'Zyx', max)) %>% unlist(rowmax)
rowmax <- tibble(genes=names(rowmax), row.max=unname(rowmax))
Spec.tt.trim <- left_join(Spec.tt, rowmax, by="genes") %>% filter(row.max >= 0.5) %>% select(-row.max)
Spec.trim.mat <- Spec.tt.trim %>% select(-1) %>% as.matrix() %>% t()
colnames(Spec.trim.mat) <- Spec.tt.trim %>% pull(genes)

pdf(file = "figures/spec_heatmap_bin5.pdf", width = 11, height = 8)
pheatmap(Spec.trim.mat, color=inferno(10), 
         clustering_distance_rows="correlation", 
         clustering_distance_cols="correlation", 
         clustering_method="ward.D", 
         show_colnames=F, treeheight_col=0, fontsize_row=8)
dev.off()

pdf(file = "figures/spec_heatmap_bin5_2.pdf", width = 11, height = 8)
pheatmap(Spec.trim.mat, color=inferno(10), 
         clustering_distance_rows="correlation", 
         clustering_distance_cols="correlation", 
         clustering_method="ward.D2", 
         show_colnames=F, treeheight_col=0, fontsize_row=8)
dev.off()
