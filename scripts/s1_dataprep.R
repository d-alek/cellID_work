# Tabula Muris - data taken from BioC - data prep, gene annotation, normalisation

# Content #####################################################################
# 1. Adding gene annotations 
# 2. Quality control on the cells
# 3. Gene-level expression metrics
# 4. Normalisation of cell-specific biases
# 5. Modeling technical trend & Detecting HVG
# 6. Select useful subsets of cells for downsteram ICA etc.
# 7. Denoising the expression values using the PCA

# various dataprep & normalisation checks figures were saved manually in figures/dataprep folder

# SmartSeq2 data ##############################################################

library(tidyverse)
library(dplyr)
library(magrittr)
library(stringr)
library(ggplot2)
library(RColorBrewer)

library(ExperimentHub)
library(SingleCellExperiment)
library(TabulaMurisData)

library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)

library(scater)
library(scran)

eh <- ExperimentHub()
query(eh, "TabulaMurisData") # dated @ April 2019


tm.smartseq2 <- TabulaMurisSmartSeq2() # alt. tm.smartseq2 <- eh[["EH1618"]]

### 1. Adding gene annotations 

# Below: Reorganising TM genes metadata 
# Issue: Important to get entrez gene ids etc. for gene set analyses
# but difficult to get unique gene mappings as TM gave only NON-UNIQUE Symbols
# in addition some of these Symbols are now old mgi symbols, some not even that
# this now works and is saved .. but the code could be somewhat simplified 
###############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# explore what's added 
names(colData(tm.smartseq2))
head(colData(tm.smartseq2))
head(colnames(tm.smartseq2))

length(unique(colData(tm.smartseq2)$plate_barcode)) # 247 diff plates!
unique(colData(tm.smartseq2)$mouse_id) # 8 mice, 10 mixes = 10 levels
unique(colData(tm.smartseq2)$tissue) # 18 tissues
unique(colData(tm.smartseq2)$cell_ontology_class) # remove NAs 

head(rowData(tm.smartseq2)$ID)
head(rowData(tm.smartseq2)$Symbol)
head(rownames(tm.smartseq2)) # all same


keytypes(org.Mm.eg.db) # get: "SYMBOL", "ENSEMBL", "ENTREZID", "GENENAME"!

# After A LOT of gene nomeclature x 3 exploration - this is the best, although fairly manual, way to get updated Gene Symbols and matched to NCBI & ENS unique ids
# (other options using Biomart or org.Mm.eg.db do not result in as good matching of troubling non-unique Symbols)

# extract TM "gene symbols"
write.csv(rownames(tm.smartseq2), row.names = F, "smartseq2.ids.csv")
# then use MGI's batch query with imput type "All Symbols/Synonims/Homologs" and  asking for "Nomenclature", "Entrez Gene ID", "Ensembl ID"... to get smartseq2.MGI.batch file
tm.smartseq2.MGI.ids <- read.table("smartseq2.MGI.batch.txt", header = T, sep = "\t", na.strings = "", as.is = T, quote = "", strip.white = T, comment.char = "")
tm.smartseq2.MGI.ids <- as_tibble(tm.smartseq2.MGI.ids) %>% distinct(.keep_all = TRUE)
# unique(tm.smartseq2.MGI.ids$Input.Type)
# unique(tm.smartseq2.MGI.ids$Chr) # .. there is "MT"


#-1: first get those that match MGI "current symbol" type - then extract their matching entrezgene id (unique)
current.symbol <- tm.smartseq2.MGI.ids %>% 
  dplyr::select(Input, Input.Type, Symbol, Entrez.Gene.ID) %>% 
  dplyr::filter(Input.Type == "current symbol") %>% 
  dplyr::arrange(Symbol, Entrez.Gene.ID) %>% 
  distinct(Input, .keep_all = TRUE) # 20943 to 20939 (this way you keep entrezgene id with smaller value)

# current.symbol %>% arrange(Symbol, Entrez.Gene.ID) %>% dplyr::filter(Symbol == lead(Symbol))
length(unique(current.symbol$Entrez.Gene.ID))-1; sum(is.na(current.symbol$Entrez.Gene.ID)) # 20914 + 25 = 20939
sum(duplicated(current.symbol$Entrez.Gene.ID)); sum(duplicated(current.symbol$Symbol)) # 24; 0
current.symbol %>% dplyr::filter(is.na(Entrez.Gene.ID)) %>% print(n=25) # 25 unmapped, can try to find them in other databases... ...
# this leaves 23433-20939=2494 TM "gene symbols" unchecked

#-2: so get those that match MGI "old symbol" type - then extract their matching entrezgene id (unique)
old.symbol <- tm.smartseq2.MGI.ids %>% 
  dplyr::select(Input, Input.Type, Symbol, Entrez.Gene.ID) %>% 
  dplyr::filter((!Input %in% current.symbol$Input) & (!Symbol %in% current.symbol$Symbol) & Input.Type == "old symbol") %>% 
  dplyr::arrange(Input, Entrez.Gene.ID) %>% 
  distinct(Input, .keep_all = TRUE) %>% 
  distinct(Symbol, .keep_all = TRUE)# 2266 to 2244 (this way you keep entrezgene id with smaller value)

# old.symbol %>% arrange(Symbol, Entrez.Gene.ID) %>% dplyr::filter(Symbol == lead(Symbol)) %>% print(n=30)
length(unique(old.symbol$Entrez.Gene.ID))-1; sum(is.na(old.symbol$Entrez.Gene.ID)) # 2243 + 1 = 2244
sum(duplicated(old.symbol$Entrez.Gene.ID)); sum(duplicated(old.symbol$Symbol)) # 0; 0
old.symbol %>% dplyr::filter(is.na(Entrez.Gene.ID)) %>% print(n=25) # 1 unmapped, can try to find it in other databases... ...
# this leaves 2494-2244=250 TM "gene symbols" unchecked


set1 <- current.symbol$Input # 20939 unique (25 no entrez id)
set1.leftover <- setdiff(rownames(tm.smartseq2), set1)# 2494 unique
set2 <- intersect(set1.leftover, old.symbol$Input) # 2244 unique (1 no entrez id)
set2.leftover <- setdiff(set1.leftover, set2) # 250 unique

# 20939 + 2244 + 250 = 23433 !!!

#-3: It is now expected these 250 are not in Input category for "current symbol" or "old symbol" as these are done already?
# if they are - then to be removed as they are then product of non-uniqueness of symbols:
mystery.symbol <- tm.smartseq2.MGI.ids %>% 
  dplyr::select(Input, Input.Type, Symbol, Entrez.Gene.ID) %>% 
  dplyr::filter((Input %in% set2.leftover) & (!Input.Type %in% c("current symbol", "old symbol"))) # 225

leftover <- mystery.symbol %>% 
  dplyr::select(Input) %>% distinct(Input) %>% 
  dplyr::filter(!str_detect(Input, "^ERCC-")) %>% dplyr::filter(!str_detect(Input, "_transgene")) %>% pull() # 124

left.leftover <- leftover[!leftover %in% c(current.symbol$Symbol, old.symbol$Symbol)] #124

test2 <- mapIds(org.Mm.eg.db, keys=left.leftover, keytype="ALIAS", column="SYMBOL", multiVals = "asNA")
test2 <- test2[!is.na(test2)] # 55 more 1:1 matches, two of them present in symbols
test2 <- test2[!test2 %in% c(current.symbol$Symbol, old.symbol$Symbol)] # this removes Clca1" "Klra4" duplicates..
test2.rev <- mapIds(org.Mm.eg.db, keys=test2, keytype="SYMBOL", column="ENTREZID", multiVals = "asNA") # 53 1:1 reverse matches 
intersect(test2.rev, c(current.symbol$Entrez.Gene.ID, old.symbol$Entrez.Gene.ID)) # 0

# now fill in these in "mystery.symbol" table
# out of 250: 225 available: 92 ERCC + 4 transgene (delete) + 53 = 145 (so 250-145=105 will be lost form the genelist)
# some of these acc. to mgi & ncbi symbols are recently discontinued to lack of evidence.. but add them anyway, some are current symbols somehow missed in mgi barch search
resolved.mystery.symbol <- tibble(Input = names(test2), Input.Type = rep("mystery symbol", 53), Symbol = test2, Entrez.Gene.ID = as.integer(test2.rev))
ercc.symbol <- mystery.symbol %>% dplyr::filter(str_detect(Input, "^ERCC-"))

#-4: you can also go back and check non-mgi entrezid matches for the missed 26 above in steps -1 and -2 and add them if unique so far:
revise1 <- current.symbol %>% dplyr::filter(is.na(Entrez.Gene.ID))
revise2 <- old.symbol %>% dplyr::filter(is.na(Entrez.Gene.ID))
revise <- bind_rows(revise1, revise2)

test3 <- mapIds(org.Mm.eg.db, keys=revise$Symbol, keytype="SYMBOL", column="ENTREZID", multiVals = "asNA")
test3 <- test3[!is.na(test3)] # 7 1:1 matches, the rest remains NA - add them below

Reduce(intersect, list(current.symbol$Input, old.symbol$Input, resolved.mystery.symbol$Input))
Reduce(intersect, list(current.symbol$Symbol, old.symbol$Symbol, resolved.mystery.symbol$Symbol))
Reduce(intersect, list(current.symbol$Entrez.Gene.ID, old.symbol$Entrez.Gene.ID, resolved.mystery.symbol$Entrez.Gene.ID)) # all clear - good to go

all.symbols <- bind_rows(current.symbol, old.symbol, resolved.mystery.symbol, ercc.symbol)
new <- tibble(Symbol = names(test3), Entrez.Gene.ID = as.integer(test3))

intersect(all.symbols$Entrez.Gene.ID, new$Entrez.Gene.ID) # all clear - good to go

all.symbols <- all.symbols %>% left_join(new, by = "Symbol") %>% 
  dplyr::mutate(Entrez.Gene.ID.x = if_else(!is.na(Entrez.Gene.ID.y), Entrez.Gene.ID.y, Entrez.Gene.ID.x)) %>% 
  dplyr::select(-Entrez.Gene.ID.y) %>% 
  rename(Entrez.Gene.ID = Entrez.Gene.ID.x)

# For the record these 105 (as mentioned above) are going to be removed from the original TM rownames:
setdiff(rownames(tm.smartseq2), all.symbols$Input) # 105 (deleted) + 23328 (in all.symbols) = 23433 (in rownames TM)

# # back-check entrez id to gene symbols consistency
# symbol.back.names.asNA <- mapIds(org.Mm.eg.db, keys=as.character(all.symbols$Entrez.Gene.ID), keytype="ENTREZID", column="SYMBOL", multiVals = "asNA")
# # ALL ARE 1:1, ALL 23217 for non-missing entrez ids
# all.equal(unname(symbol.back.names.asNA), all.symbols %>% dplyr::filter(!is.na(Entrez.Gene.ID)) %>% dplyr::select(Symbol) %>% pull()) # 182 mismatches - pretty good, these are mainly old-new mismatches


# now you can get ENS names and full gene names 
sum(is.na(all.symbols$Entrez.Gene.ID)) # 111
ensembl.names.asNA <- mapIds(org.Mm.eg.db, keys=as.character(all.symbols$Entrez.Gene.ID), keytype="ENTREZID", column="ENSEMBL", multiVals = "asNA")
ensembl.names.filt <- mapIds(org.Mm.eg.db, keys=as.character(all.symbols$Entrez.Gene.ID), keytype="ENTREZID", column="ENSEMBL", multiVals = "filter")
ensembl.names.1st <- mapIds(org.Mm.eg.db, keys=sort(unique(as.character(all.symbols$Entrez.Gene.ID))), keytype="ENTREZID", column="ENSEMBL") # as vector if NAs removed (with sort() here)

length(sort(unique(as.character(all.symbols$Entrez.Gene.ID)))) # 23217 nonNA input (23217 + 111 = 23328)

length(ensembl.names.asNA) # 23217 
sum(is.na(ensembl.names.asNA)) # 826 (497 no matches + 329 dups matches)
sum(duplicated(sort(ensembl.names.asNA))) # 49 NON-UNIQUE resulting ENS ids ! (not counting NAs) 
length(unique(sort(ensembl.names.asNA))) # 22342 (not counting NAs) (22342 + 826 + 49 = 23217)

length(ensembl.names.filt) # 22888 (23217 - 329 = 22888) : non-unique excluded
sum(is.na(ensembl.names.filt)) # 497 : no-matches NAed
sum(duplicated(sort(ensembl.names.filt))) # 49 NON-UNIQUE resulting ENS ids ! (not counting NAs) 
length(unique(sort(ensembl.names.filt))) # 22342 (not counting NAs) (22342 + 497 + 49 = 22888)

# here you see how many just do not have ENS match (not worrying about uniqueness of match)
length(ensembl.names.1st) # 23217 
sum(is.na(ensembl.names.1st)) # 497 DO NOT HAVE MATCH AT ALL
sum(duplicated(sort(ensembl.names.1st))) # 52 NON-UNIQUE resulting ENS ids ! (not counting NAs) - so cannot us ethis approach
length(unique(sort(ensembl.names.1st))) # 22668 (not counting NAs) (22668 + 497 + 52 = 23217) + 111 NAs = 23328

# 497 does not have ncbi-ens match at all / 22720 does have some match:
# 329 has multiple ncbi-ens matches / 22391 has unique ncbi-ens matches - that is OK, accept the first match (the default option on mapIds())
# 49 duplicate ens results even when only 22391 unique ncbi-ens matches considered (goes to 52 if additional 329 multiple ncbi-ens matches included)

# IN SHORT GO WITH ensembl.names.1st APPROACH - DEFAULT (just keep in mind there are 52 ENS dupicates)

# locate and check ENS duplication:
# many of these are corresponding to one normal Gene Symbol, one Gm.. / ..Rik type of preliminary nomenclature symbol
# many of them are genes in indentical regions - just happened to have different Symbol and Entrez id, but same ENS id
dups <- sort(ensembl.names.1st)[duplicated(sort(ensembl.names.1st))]
length(unique(dups)) # 49... makes 52

for (i in 1:49) {
  dups.i <- ensembl.names.1st[ensembl.names.1st %in% dups[i]]
  print(unname(dups[i]))
  print(all.symbols[all.symbols$Entrez.Gene.ID %in% names(dups.i), ])
  cat("\n")
}
# so leave them for now and here is the list if needs to be checked later 


# FINALLY
# Add ENS ids to set of gene annotations (22720 all matches (52 of them ENS duplicates) + 111 entrez NAs + 497 ensemble NAs)
new.ens <- tibble(Entrez.Gene.ID = as.integer(names(ensembl.names.1st)), ENS.Gene.ID = ensembl.names.1st)
all.symbols <- all.symbols %>% left_join(new.ens, by = "Entrez.Gene.ID")

sum(is.na(all.symbols$ENS.Gene.ID)) # 608 = 111 + 497
length(unique(all.symbols$ENS.Gene.ID)) # 22669 (less than 22720 as there are some duplicates)

# ADD GENE NAMES & FEATURE TYPE TOO
# extract TM "new gene symbols"
write.csv(all.symbols$Symbol, row.names = F, "smartseq2.new.symbols.csv")
# then use MGI's batch query with imput type "All Symbols/Synonims/Homologs" and  asking for "Nomenclature", "Entrez Gene ID", "Ensembl ID"... to get smartseq2.MGI.batch file
tm.smartseq2.MGI.new.symbol <- read.table("smartseq2.MGI.new.symbol.batch.txt", header = T, sep = "\t", na.strings = "", as.is = T, quote = "", strip.white = T, comment.char = "")
tm.smartseq2.MGI.new.symbol <- as_tibble(tm.smartseq2.MGI.new.symbol) %>% distinct(.keep_all = TRUE) %>% dplyr::filter(Input != "x")
unique(tm.smartseq2.MGI.new.symbol$Input.Type) # "current symbol" NA
length(unique(tm.smartseq2.MGI.new.symbol$Symbol)) # 23222-1=23221
sum(is.na(tm.smartseq2.MGI.new.symbol$Symbol)) # 17
# 23221+17=23238 all good
sort(tm.smartseq2.MGI.new.symbol$Name)[duplicated(sort(tm.smartseq2.MGI.new.symbol$Name))] # 11
tm.smartseq2.MGI.new.symbol[is.na(tm.smartseq2.MGI.new.symbol$Name), ] # 17

# first transfer ERCC names from all.symbols Input to Symbols:
all.symbols <- all.symbols %>% dplyr::mutate(Symbol = if_else(is.na(Symbol), Input, Symbol))
# then attach "name" and "feature type"
new.names <- tm.smartseq2.MGI.new.symbol %>% dplyr::select(Symbol, Name, Feature.Type)
all.symbols <- all.symbols %>% left_join(new.names, by = "Symbol")

# checks...
sum(is.na(all.symbols$Input)); sum(duplicated(all.symbols$Input)); length(unique(all.symbols$Input))
sum(is.na(all.symbols$Symbol)); sum(duplicated(all.symbols$Symbol)); length(unique(all.symbols$Symbol))
sum(is.na(all.symbols$Entrez.Gene.ID)); sum(duplicated(all.symbols$Entrez.Gene.ID)); length(unique(all.symbols$Entrez.Gene.ID))-1
sum(is.na(all.symbols$ENS.Gene.ID)); sum(duplicated(all.symbols$ENS.Gene.ID)); length(unique(all.symbols$ENS.Gene.ID))-1
sum(is.na(all.symbols$Name)); sum(duplicated(all.symbols$Name)); length(unique(all.symbols$Name))-1

write_tsv(all.symbols, "all.symbols.tsv")

# FINALLY Add as metadata 

# explore what's added already
names(colData(tm.smartseq2))
head(colData(tm.smartseq2))
head(colnames(tm.smartseq2))

head(rowData(tm.smartseq2)$ID)
head(rowData(tm.smartseq2)$Symbol)
head(rownames(tm.smartseq2)) # all same

# essentially keep original TM rownames = ID = Symbol as "ID" in rowData
# set new rownames to new extracted Symbol, same Symbol column in duplicated in rowData

keep.genes <- rownames(tm.smartseq2) %in% all.symbols$Input
tm.smartseq2 <- tm.smartseq2[keep.genes, ]

m <- match(rownames(tm.smartseq2), all.symbols$Input)
all.equal(rownames(tm.smartseq2), all.symbols$Input[m])
all.symbols <- all.symbols[m,] # so they are now ordered same as in sce object
all.equal(rownames(tm.smartseq2), all.symbols$Input)
all.equal(rowData(tm.smartseq2)$ID, all.symbols$Input)

rownames(tm.smartseq2) <- all.symbols$Symbol
rowData(tm.smartseq2)$Symbol <- all.symbols$Symbol
rowData(tm.smartseq2)$entrez <- all.symbols$Entrez.Gene.ID
rowData(tm.smartseq2)$ensembl <- all.symbols$ENS.Gene.ID
rowData(tm.smartseq2)$Name <- all.symbols$Name
rowData(tm.smartseq2)$FeatureType <- all.symbols$Feature.Type
all.equal(rownames(tm.smartseq2), rowData(tm.smartseq2)$Symbol)
names(rowData(tm.smartseq2))

saveRDS(file="BioC_TM_smartseq2_GenesMeta.rds", tm.smartseq2)
tm.smartseq2 <- readRDS("BioC_TM_smartseq2_GenesMeta.rds")


#-------------------------------------------------------------------------------------------------#
# below is a useful exercise in BiomaRt.. but still retrieves lots of "old symbols" by MGI standards
library(biomaRt)
?useMart
listMarts() # ENSEMBL_MART_ENSEMBL
listDatasets(mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")) # mmusculus_gene_ensembl

# 1) select a mart and data set     
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host="www.ensembl.org")
listFilters(mart=mart, what=c("name", "description"))
searchFilters(mart=mart, pattern="symbol") # mgi_symbol
searchFilters(mart=mart, pattern="entrezgene") # entrezgene_accession
searchFilters(mart=mart, pattern="external") # external_gene_name

listAttributes(mart=mart, what=c("name","description","page"), page="feature_page")
searchAttributes(mart=mart, pattern="symbol") # mgi_symbol
searchAttributes(mart=mart, pattern="entrezgene") # entrezgene_id, entrezgene_accession, entrezgene_description
searchAttributes(mart=mart, pattern="ensembl") # ensembl_gene_id

# 2) run a biomart query using the getBM() function and specify the attributes and filter arguments    
# MGI SYMBOLS to ENTREZID and BACK------------------------#
bm.mgi_symbol_entrezgene <- biomaRt::getBM(attributes = c("mgi_symbol", "entrezgene_id", "entrezgene_accession"),   
                                           filters    = "mgi_symbol",     
                                           values     = rownames(tm.smartseq2),       
                                           mart       = mart)
sum(is.na(bm.mgi_symbol_entrezgene$mgi_symbol)); sum(duplicated(bm.mgi_symbol_entrezgene$mgi_symbol)); length(unique(bm.mgi_symbol_entrezgene$mgi_symbol)) # 0; 129; 20868 / =20997
sum(is.na(bm.mgi_symbol_entrezgene$entrezgene_id)); length(unique(bm.mgi_symbol_entrezgene$entrezgene_id)) # 2402; 18565 / =20997
sum(bm.mgi_symbol_entrezgene$entrezgene_accession==""); length(unique(bm.mgi_symbol_entrezgene$entrezgene_accession)) # 2402; 18565 / =20997
# so you have 18565-129=18436 unique matchings there.. but when tested on some examples, lots of them stil lead to old symbols mapping

entreznames <- unique(bm.mgi_symbol_entrezgene$entrezgene_id)
# now from those picked unique entrezgenes to mgi.. is it going to give you current gene symbols (e.g. H1f8 instead of H1foo)?
bm.entrezgene_mgi_symbol <- biomaRt::getBM(attributes = c("entrezgene_id", "entrezgene_accession", "mgi_symbol"),   
                                           filters    = "entrezgene_id",     
                                           values     = entreznames,       
                                           mart       = mart)
sum(is.na(bm.entrezgene_mgi_symbol$entrezgene_id)); sum(duplicated(bm.entrezgene_mgi_symbol$entrezgene_id)); length(unique(bm.entrezgene_mgi_symbol$entrezgene_id)) # 0; 110; 18564 / =18674
sum(is.na(bm.entrezgene_mgi_symbol$mgi_symbol)); sum(duplicated(bm.entrezgene_mgi_symbol$mgi_symbol)); length(unique(bm.entrezgene_mgi_symbol$mgi_symbol)) # 0; 137; 18537 / =18674
# BUT NO THEY STILL HAVE TENDENCY TO OLD SYMBOLS !

# MGI SYMBOLS to ENSEMBL and BACK ------------------------#
bm.mgi_symbol_ensembl <- biomaRt::getBM(attributes = c("mgi_symbol", "ensembl_gene_id"),   
                                           filters    = "mgi_symbol",     
                                           values     = rownames(tm.smartseq2),       
                                           mart       = mart)
sum(is.na(bm.mgi_symbol_ensembl$mgi_symbol)); sum(duplicated(bm.mgi_symbol_ensembl$mgi_symbol)); length(unique(bm.mgi_symbol_entrezgene$mgi_symbol)) # 0; 325; 20868 / =21193
sum(is.na(bm.mgi_symbol_ensembl$ensembl_gene_id)); sum(duplicated(bm.mgi_symbol_ensembl$ensembl_gene_id)); length(unique(bm.mgi_symbol_ensembl$ensembl_gene_id)) # 0; 4; 21189 / =21193

ensemblnames <- unique(bm.mgi_symbol_ensembl$ensembl_gene_id)
bm.ensembl_mgi_symbol <- biomaRt::getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),   
                                           filters    = "ensembl_gene_id",     
                                           values     = ensemblnames,       
                                           mart       = mart)
sum(is.na(bm.ensembl_mgi_symbol$ensembl_gene_id)); sum(duplicated(bm.ensembl_mgi_symbol$ensembl_gene_id)); length(unique(bm.ensembl_mgi_symbol$ensembl_gene_id)) # 0; 5; 20869 / =21194
sum(is.na(bm.ensembl_mgi_symbol$mgi_symbol)); sum(duplicated(bm.ensembl_mgi_symbol$mgi_symbol)); length(unique(bm.ensembl_mgi_symbol$mgi_symbol)) # 0; 325; 20869 / =21194
# SAME, THEY STILL HAVE TENDENCY TO OLD SYMBOLS !

#-------------------------------------------------------------------------------------------------#
# an alternative example of how to add feature metadata quickly via scater-biomaRt.. of course does not work well with non-unique filters :)
marrow.BioC <- getBMFeatureAnnos(marrow.BioC, 
                                 filters = c("mgi_symbol"), 
                                 attributes = c("mgi_symbol", "ensembl_gene_id", "entrezgene", "mgi_description", 
                                                "gene_biotype", "chromosome_name", "start_position", "end_position", "transcript_length"),  
                                 biomart = "ENSEMBL_MART_ENSEMBL", 
                                 dataset = "mmusculus_gene_ensembl", 
                                 host = "www.ensembl.org")
#-------------------------------------------------------------------------------------------------#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###############################################################################


### 2. Quality control on the cells ----------------------------------------###

## 2.a Remove non-annotated cells now, and vaguely annotated cells LATER ON (to not affect the number of cells per plate too much)
drop.cells <- is.na(colData(tm.smartseq2)$cell_ontology_class) 
# few NA cells are present on most plates, but 14 plates have more than 100, 7 plates have more than 200 of them (could affect multiblock normalisation if removed now)
# drop.cells <- colData(tm.smartseq2)$cell_ontology_class %in% c("blood cell", "myeloid cell", "professional antigen presenting cell", "leukocyte", "lymphocyte")
tm.smartseq2 <- tm.smartseq2[, !drop.cells] # 8981 NA cells removed (these are mostly cells removed on GC ground by TM, and some non-annotated)

# explore celltype - tissue categories
tibble(celltype=tm.smartseq2$cell_ontology_class) %>% 
  dplyr::count(celltype) %>% dplyr::arrange(celltype) %>% print(n=115)

tibble(celltype=tm.smartseq2$cell_ontology_class, tissue=tm.smartseq2$tissue) %>% 
  dplyr::count(celltype,tissue) %>% dplyr::arrange(celltype,tissue) %>% print(n=115)

# check if any plates have too few cells - also check this after QC
tibble(plate=tm.smartseq2$plate_barcode) %>% 
  dplyr::count(plate) %>% dplyr::arrange(n) %>% print(n=10) # 5 plates with less than 20 cells..


## 2.b Calculate QC metrics 
tm.smartseq2 <- calculateQCMetrics(tm.smartseq2, use_spikes = TRUE)

# create tissue-celltype category to subset on prior to quality cutoffs (ideally some small but similar cell types above could be fused here to avoid small groups)
tm.smartseq2$tissue_celltype <- paste(tm.smartseq2$tissue, tm.smartseq2$cell_ontology_class, sep = "-")

# before 
p1 <- plotColData(tm.smartseq2, y="pct_counts_ERCC", x="tissue_celltype") + scale_x_discrete(labels = NULL)
p2 <- plotColData(tm.smartseq2, y="total_counts", x="tissue_celltype") + scale_x_discrete(labels = NULL)
p3 <- plotColData(tm.smartseq2, y="total_features_by_counts", x="tissue_celltype") + scale_x_discrete(labels = NULL)
multiplot(p1,p2,p3, cols=1)
range(tm.smartseq2$total_counts)
range(tm.smartseq2$total_features_by_counts)
range(tm.smartseq2$pct_counts_ERCC)


# 2.c Make a cut (at 4 nmads because data has been already trimmed by TM in a global way)
libsize.drop <- isOutlier(tm.smartseq2$total_counts, nmads=4, type="lower", log=TRUE, batch=tm.smartseq2$tissue_celltype)
feature.drop <- isOutlier(tm.smartseq2$total_features_by_counts, nmads=4, type="lower", log=TRUE, batch=tm.smartseq2$tissue_celltype)
spike.drop <- isOutlier(tm.smartseq2$pct_counts_ERCC, nmads=4, type="higher", batch=tm.smartseq2$tissue_celltype) | tm.smartseq2$pct_counts_ERCC > 50

keep <- !(libsize.drop | feature.drop | spike.drop) # loose 2166/44779 non-NA cells = 4.8%
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), BySpike=sum(spike.drop), Remaining=sum(keep))

tm.smartseq2 <- tm.smartseq2[, keep]

# after
p1 <- plotColData(tm.smartseq2, y="pct_counts_ERCC", x="tissue_celltype") + scale_x_discrete(labels = NULL)
p2 <- plotColData(tm.smartseq2, y="total_counts", x="tissue_celltype") + scale_x_discrete(labels = NULL)
p3 <- plotColData(tm.smartseq2, y="total_features_by_counts", x="tissue_celltype") + scale_x_discrete(labels = NULL)
multiplot(p1,p2,p3, cols=1)
range(tm.smartseq2$total_counts)
range(tm.smartseq2$total_features_by_counts)
range(tm.smartseq2$pct_counts_ERCC)

# check again if any plates have too few cells
tibble(plate=tm.smartseq2$plate_barcode) %>% 
  dplyr::count(plate) %>% dplyr::arrange(n) %>% print(n=20) # now 11 plates with less than 20 HQ cells, 3 with less than 10 HQ ceels..(MAA000924 :7, MAA001856 :7, MAA000538 :9)

tibble(plate=tm.smartseq2$plate_barcode, celltype=tm.smartseq2$cell_ontology_class) %>% 
  dplyr::filter(plate %in% c("MAA000924", "MAA001856", "MAA000538", "B003284")) %>% dplyr::arrange(plate) %>% print(n=100)

# it's safe to exclude the few plates with 10 HQ cells or less to get better texchical fit below (unfortunately loosing 3 out of 25 rare Bergmann glial cell):
plate.drop <- tm.smartseq2$plate_barcode %in% c("MAA000924", "MAA001856", "MAA000538", "B003284")
tm.smartseq2 <- tm.smartseq2[, !plate.drop]



### 3. Gene-level expression metrics ---------------------------------------###

## 3.a Highly-expressed features
plotHighestExprs(tm.smartseq2, n=25) + theme(axis.text=element_text(size=11), axis.title=element_text(size=14))


## 3.b Filter-out low abundance genes

# by average counts across all cells
ave.counts <- calcAverage(tm.smartseq2, use_size_factors=FALSE)
hist(log10(ave.counts), breaks=100, main="", col="grey80", xlab=expression(Log[10]~"average count")); abline(v = log10(0.005), col = "red")

demo.keep <- ave.counts >= 0.005 | grepl("^ERCC-", names(ave.counts))
summary(demo.keep)

# by number of expressing cells
num.cells <- nexprs(tm.smartseq2, byrow=TRUE)
smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells", xlab=expression(Log[10]~"average count")); abline(v = log10(0.005), h = 10, col = "red")

to.keep <- num.cells > 10 | grepl("^ERCC-", names(num.cells))
summary(to.keep)

comb.keep <- demo.keep & to.keep
summary(comb.keep)

tm.smartseq2 <- tm.smartseq2[comb.keep,] # 20698 genes kept for now

summary(isSpike(tm.smartseq2, "ERCC"))



### 4. Normalisation of cell-specific biases -------------------------------###

## 4.a Normalising based on spike-in size factors only 
# i.e. preserving the cell type difference in total RNA content (possibly more meaningful for TM)!!!

tm.smartseq2.spikeF <- computeSpikeFactors(tm.smartseq2, type="ERCC", general.use=TRUE)
head(sizeFactors(tm.smartseq2.spikeF))
head(sizeFactors(tm.smartseq2.spikeF, "ERCC")) # same as general size factors
summary(sizeFactors(tm.smartseq2.spikeF))

plot(tm.smartseq2.spikeF$total_counts/1e6, sizeFactors(tm.smartseq2.spikeF), log="xy",
     xlab="Library size (millions)", ylab="Size factor", pch=21, col="dodgerblue3", main="Size Factors - based on spike-ins only")

# Normalise logcounts & counts:
# tm.smartseq2.spikeF.N <- normalize(tm.smartseq2.spikeF) # in logcounts slot -NO NEED TO DO THIS AS YOU NEED MULTIBATCH NORMALISE HERE
# tm.smartseq2.spikeF.N <- normalize(tm.smartseq2.spikeF, return_log = FALSE) # in normcounts slot


## 4.b Normalising based on combination of deconvolution and spike-in size factors 
# i.e. erasing the cell type difference in total RNA content 

# deconvolution factors - 1st take to get scaling input
sumF.1st <- computeSumFactors(tm.smartseq2, sizes=seq(21, 101, 5), cluster=tm.smartseq2$cell_ontology_class, max.cluster.size=4500, positive=TRUE, min.mean=1, sf.out=TRUE)
head(sumF.1st)
summary(sumF.1st) # if run with "positive=FALSE" results in 20 slightly negative factors - not bad, but could have gone more stringent on cell QC..
# deconvolution factors - 2nd take, with scaling input from first run, to get improved 
tm.smartseq2.sumF <- computeSumFactors(tm.smartseq2, sizes=seq(21, 101, 5), cluster=tm.smartseq2$cell_ontology_class, max.cluster.size=4500, positive=TRUE, min.mean=1, scaling=sumF.1st) 
head(sizeFactors(tm.smartseq2.sumF))
summary(sizeFactors(tm.smartseq2.sumF))
# separate spike-in factors
tm.smartseq2.sumF <- computeSpikeFactors(tm.smartseq2.sumF, type="ERCC", general.use=FALSE)
head(sizeFactors(tm.smartseq2.sumF, "ERCC"))
summary(sizeFactors(tm.smartseq2.sumF, "ERCC"))

plot(tm.smartseq2.sumF$total_counts/1e6, sizeFactors(tm.smartseq2.sumF), log="xy",
     xlab="Library size (millions)", ylab="Size factor", pch=21, col="tomato3", main="Size Factors - based on deconvolution method")

# Normalise logcounts & counts:
# tm.smartseq2.sumF.N <- normalize(tm.smartseq2.sumF) # in logcounts slot -NO NEED TO DO THIS AS YOU NEED MULTIBATCH NORMALISE HERE
# tm.smartseq2.sumF <- normalize(tm.smartseq2.sumF, return_log = FALSE) # in normcounts slot


# comparing spike-in & deconvolution size factors
plot(sizeFactors(tm.smartseq2.spikeF), sizeFactors(tm.smartseq2.sumF), log="xy",
     xlab="Size factor (spike-in)", ylab="Size factor (deconvolution)", pch=21, col="slategray", main="Size Factors - comparison")

# Now the question is if the GENERAL normalize() is good enough here.. or to use multiBlockNorm() to account for plates too
# on one side - need to take care of plates in trendVar() below.. 
# on the other side - normalise() size factors are formed based on clusters of cell types - and they overlap with plates !!!


### 5. Modeling technical trend & Detecting HVG

# Two options here: 
# A) use normalise() then trendVar(block=..) & decomposeVar() - fits a single trend to the plate-specific means and variances of spike-in transcripts.. but trend might not be same between plates!
# B) use multiBlockNorm() then multiBlockVar() - fits plate-specific trends, then combines - having in mind large number of plates in TM this approach is more reasonable here!?
# so go with option B

names(rowData(tm.smartseq2.spikeF.BN))
anyNA(rowData(tm.smartseq2.spikeF.BN)$log10_mean_counts)
summary(rowData(tm.smartseq2.spikeF.BN)$log10_mean_counts)
sum(rowData(tm.smartseq2.spikeF.BN)$is_feature_control_ERCC)


## 5.a Trends based on spike-in size factors normalisation 
tm.smartseq2.spikeF.BN <- multiBlockNorm(tm.smartseq2.spikeF, block=as.factor(tm.smartseq2.spikeF$plate_barcode)) # logcounts slot
# no change in logcounts from normalise() !! - this is EXPECTED as spike-in factors are not rescaled here as there is no gene-based factors in this version  
# tm.smartseq2.spikeF.BN <- multiBlockNorm(tm.smartseq2.spikeF, block=as.factor(tm.smartseq2.spikeF$plate_barcode), return_log = FALSE) # normcounts slot

#-- loess (spike)
comb.spikeF.var.out.loess <- multiBlockVar(tm.smartseq2.spikeF.BN, block=as.factor(tm.smartseq2.spikeF.BN$plate_barcode), 
                                           trend.args=list(parametric=TRUE, method="loess", loess.args=list(span=0.3)), assay.type="logcounts")
set.seed(101)
par(mfrow=c(5,5))
is.spike <- isSpike(tm.smartseq2.spikeF.BN)
for (plate in sample(levels(as.factor(tm.smartseq2.spikeF.BN$plate_barcode)), 25)) {
  cur.out <- comb.spikeF.var.out.loess$per.block[[plate]]
  plot(cur.out$mean, cur.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
       ylab="Variance of log-expression", ylim=c(0, 40), main=plate)
  curve(metadata(cur.out)$trend(x), col="dodgerblue", lwd=2, add=TRUE)
  points(cur.out$mean[is.spike], cur.out$total[is.spike], col="red", pch=16)
}
sum(comb.spikeF.var.out.loess$FDR <= 0.05, na.rm = T) # 8876
sum(comb.spikeF.var.out.loess$FDR <= 0.05 & comb.spikeF.var.out.loess$bio > 0, na.rm = T) # 8109

hvg.out.1 <- comb.spikeF.var.out.loess[,1:6][which(comb.spikeF.var.out.loess$FDR <= 0.05),]
hist(hvg.out.1$bio, breaks = 200, col = "grey80")

hvg.out.1 <- comb.spikeF.var.out.loess[,1:6][which(comb.spikeF.var.out.loess$FDR <= 0.05 & comb.spikeF.var.out.loess$bio > 0),]
hvg.out.1 <- hvg.out.1[order(hvg.out.1$bio, decreasing=TRUE),] 
summary(hvg.out.1$bio) # 1st Qu 0.231825
hist(hvg.out.1$bio, breaks = 200, col = "grey80"); abline(v=summary(hvg.out.1$bio)[c(2,3)], col="red")
hvg.out.1 <- hvg.out.1[which(hvg.out.1$bio >= summary(hvg.out.1$bio)[2]), ]; nrow(hvg.out.1) # 6082
# write.table(file="hsc_hvg.tsv", hvg.out, sep="\t", quote=FALSE, col.names=NA)
plotExpression(tm.smartseq2.spikeF.BN, features=rownames(hvg.out.1)[1:20]) + theme(axis.text=element_text(size=11), axis.title=element_text(size=14))

#-- spline (spike)
comb.spikeF.var.out.spline <- multiBlockVar(tm.smartseq2.spikeF.BN, block=as.factor(tm.smartseq2.spikeF.BN$plate_barcode), 
                                           trend.args=list(parametric=TRUE, method = "spline"), assay.type="logcounts")
set.seed(101)
par(mfrow=c(5,5))
is.spike <- isSpike(tm.smartseq2.spikeF.BN)
for (plate in sample(levels(as.factor(tm.smartseq2.spikeF.BN$plate_barcode)), 25)) {
  cur.out <- comb.spikeF.var.out.spline$per.block[[plate]]
  plot(cur.out$mean, cur.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
       ylab="Variance of log-expression", ylim=c(0, 40), main=plate)
  curve(metadata(cur.out)$trend(x), col="dodgerblue", lwd=2, add=TRUE)
  points(cur.out$mean[is.spike], cur.out$total[is.spike], col="red", pch=16)
}
sum(comb.spikeF.var.out.spline$FDR <= 0.05, na.rm = T) # 8934
sum(comb.spikeF.var.out.spline$FDR <= 0.05 & comb.spikeF.var.out.spline$bio > 0, na.rm = T) # 7842

hvg.out.2 <- comb.spikeF.var.out.spline[,1:6][which(comb.spikeF.var.out.spline$FDR <= 0.05),]
hist(hvg.out.2$bio, breaks = 200, col = "grey80")

hvg.out.2 <- comb.spikeF.var.out.spline[,1:6][which(comb.spikeF.var.out.spline$FDR <= 0.05 & comb.spikeF.var.out.spline$bio > 0),]
hvg.out.2 <- hvg.out.2[order(hvg.out.2$bio, decreasing=TRUE),] 
summary(hvg.out.2$bio) # 1st Qu 0.212451
hist(hvg.out.2$bio, breaks = 200, col = "grey80"); abline(v=summary(hvg.out.2$bio)[c(2,3)], col="red")
hvg.out.2 <- hvg.out.2[which(hvg.out.2$bio >= summary(hvg.out.2$bio)[2]), ]; nrow(hvg.out.2) # 5881
# write.table(file="hsc_hvg.tsv", hvg.out, sep="\t", quote=FALSE, col.names=NA)
plotExpression(tm.smartseq2.spikeF.BN, features=rownames(hvg.out.2)[1:50]) + theme(axis.text=element_text(size=11), axis.title=element_text(size=14))

# unite spikein ones 
hvg.spikeNorm <- union(rownames(hvg.out.1), rownames(hvg.out.2))

# compare
length(union(rownames(hvg.out.1), rownames(hvg.out.2))) # 6246
length(intersect(rownames(hvg.out.1), rownames(hvg.out.2))) # 5717
length(setdiff(rownames(hvg.out.1), rownames(hvg.out.2))) # 365
length(setdiff(rownames(hvg.out.2), rownames(hvg.out.1))) # 164

summary(hvg.out.1[setdiff(rownames(hvg.out.1), rownames(hvg.out.2)), ]$bio)
summary(hvg.out.2[setdiff(rownames(hvg.out.2), rownames(hvg.out.1)), ]$bio)



## 5.b Trends based on deconvolution size factors normalisation 
tm.smartseq2.sumF.BN <- multiBlockNorm(tm.smartseq2.sumF, block=as.factor(tm.smartseq2.sumF$plate_barcode)) # logcounts slot
# here logcounts are different from normalise() !!
# tm.smartseq2.sumF.BN <- multiBlockNorm(tm.smartseq2.sumF, block=as.factor(tm.smartseq2.sumF$plate_barcode), return_log = FALSE) # normcounts slot

#-- loess (deconvoulution)
comb.sumF.var.out.loess <- multiBlockVar(tm.smartseq2.sumF.BN, block=as.factor(tm.smartseq2.sumF.BN$plate_barcode), 
                                           trend.args=list(parametric=TRUE, loess.args=list(span=0.4)), assay.type="logcounts")
set.seed(101)
par(mfrow=c(5,5))
is.spike <- isSpike(tm.smartseq2.sumF.BN)
for (plate in sample(levels(as.factor(tm.smartseq2.sumF.BN$plate_barcode)), 25)) {
  cur.out <- comb.sumF.var.out.loess$per.block[[plate]]
  plot(cur.out$mean, cur.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
       ylab="Variance of log-expression", ylim=c(0, 40), main=plate)
  curve(metadata(cur.out)$trend(x), col="dodgerblue", lwd=2, add=TRUE)
  points(cur.out$mean[is.spike], cur.out$total[is.spike], col="red", pch=16)
}
sum(comb.sumF.var.out.loess$FDR <= 0.05, na.rm = T) # 6946
sum(comb.sumF.var.out.loess$FDR <= 0.05 & comb.sumF.var.out.loess$bio > 0, na.rm = T) # 6022

hvg.out.3 <- comb.sumF.var.out.loess[,1:6][which(comb.sumF.var.out.loess$FDR <= 0.05),]
hist(hvg.out.3$bio, breaks = 200, col = "grey80")

hvg.out.3 <- comb.sumF.var.out.loess[,1:6][which(comb.sumF.var.out.loess$FDR <= 0.05 & comb.sumF.var.out.loess$bio > 0),]
hvg.out.3 <- hvg.out.3[order(hvg.out.3$bio, decreasing=TRUE),] 
summary(hvg.out.3$bio) # 1st Qu 0.165212
hist(hvg.out.3$bio, breaks = 200, col = "grey80"); abline(v=summary(hvg.out.3$bio)[c(2,3)], col="red")
hvg.out.3 <- hvg.out.3[which(hvg.out.3$bio >= summary(hvg.out.3$bio)[2]), ]; nrow(hvg.out.3) # 4516
# write.table(file="hsc_hvg.tsv", hvg.out, sep="\t", quote=FALSE, col.names=NA)
plotExpression(tm.smartseq2.sumF.BN, features=rownames(hvg.out.3)[1:50]) + theme(axis.text=element_text(size=11), axis.title=element_text(size=14))


#-- spline (deconvoulution)
comb.sumF.var.out.spline <- multiBlockVar(tm.smartseq2.sumF.BN, block=as.factor(tm.smartseq2.sumF.BN$plate_barcode), 
                                            trend.args=list(parametric=TRUE, method = "spline"), assay.type="logcounts")
set.seed(101)
par(mfrow=c(5,5))
is.spike <- isSpike(tm.smartseq2.sumF.BN)
for (plate in sample(levels(as.factor(tm.smartseq2.sumF.BN$plate_barcode)), 25)) {
  cur.out <- comb.sumF.var.out.spline$per.block[[plate]]
  plot(cur.out$mean, cur.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
       ylab="Variance of log-expression", ylim=c(0, 40), main=plate)
  curve(metadata(cur.out)$trend(x), col="dodgerblue", lwd=2, add=TRUE)
  points(cur.out$mean[is.spike], cur.out$total[is.spike], col="red", pch=16)
}
sum(comb.sumF.var.out.spline$FDR <= 0.05, na.rm = T) # 7065
sum(comb.sumF.var.out.spline$FDR <= 0.05 & comb.sumF.var.out.spline$bio > 0, na.rm = T) # 5960

hvg.out.4 <- comb.sumF.var.out.spline[,1:6][which(comb.sumF.var.out.spline$FDR <= 0.05),]
hist(hvg.out.4$bio, breaks = 200, col = "grey80")

hvg.out.4 <- comb.sumF.var.out.spline[,1:6][which(comb.sumF.var.out.spline$FDR <= 0.05 & comb.sumF.var.out.spline$bio > 0),]
hvg.out.4 <- hvg.out.4[order(hvg.out.4$bio, decreasing=TRUE),] 
summary(hvg.out.4$bio) # 1st Qu 0.157864
hist(hvg.out.4$bio, breaks = 200, col = "grey80"); abline(v=summary(hvg.out.4$bio)[c(2,3)], col="red")
hvg.out.4 <- hvg.out.4[which(hvg.out.4$bio >= summary(hvg.out.4$bio)[2]), ]; nrow(hvg.out.4) # 4470
# write.table(file="hsc_hvg.tsv", hvg.out, sep="\t", quote=FALSE, col.names=NA)
plotExpression(tm.smartseq2.sumF.BN, features=rownames(hvg.out.4)[1:50]) + theme(axis.text=element_text(size=11), axis.title=element_text(size=14))

# unite deconvo ones 
hvg.sumNorm <- union(rownames(hvg.out.3), rownames(hvg.out.4))

# compare
length(union(rownames(hvg.out.3), rownames(hvg.out.4))) # 4748
length(intersect(rownames(hvg.out.3), rownames(hvg.out.4))) # 4238
length(setdiff(rownames(hvg.out.3), rownames(hvg.out.4))) # 278
length(setdiff(rownames(hvg.out.4), rownames(hvg.out.3))) # 232

summary(hvg.out.3[setdiff(rownames(hvg.out.3), rownames(hvg.out.4)), ]$bio)
summary(hvg.out.4[setdiff(rownames(hvg.out.4), rownames(hvg.out.3)), ]$bio)

#---#---#---#

# unite/intersect all
hvg <- union(hvg.spikeNorm, hvg.sumNorm) # 6342 - variability in total RNA across cell types was RETAINED here
HVG <- hvg.sumNorm # 4748 - variability in total RNA across cell types was REMOVED here

length(setdiff(hvg.spikeNorm, hvg.sumNorm)) # 1594
length(setdiff(hvg.sumNorm, hvg.spikeNorm)) # 96

# Atach hvg and trends to RowData
hvg.keep.spikelike <- rownames(tm.smartseq2.spikeF.BN) %in% hvg
rowData(tm.smartseq2.spikeF.BN)$is_hvg <- hvg.keep.spikelike
rowData(tm.smartseq2.spikeF.BN)$bio_spike <- comb.spikeF.var.out.loess[1:6]$bio
rowData(tm.smartseq2.spikeF.BN)$tech_spike <- comb.spikeF.var.out.loess[1:6]$tech

HVG.keep.sumlike <- rownames(tm.smartseq2.sumF.BN) %in% HVG
rowData(tm.smartseq2.sumF.BN)$is_HVG <- HVG.keep.sumlike
rowData(tm.smartseq2.sumF.BN)$bio_sum <- comb.sumF.var.out.loess[1:6]$bio
rowData(tm.smartseq2.sumF.BN)$tech_sum<- comb.sumF.var.out.loess[1:6]$tech

# Save the two versions (spikeF and sumF)
saveRDS(tm.smartseq2.spikeF.BN, "data/tm.smartseq2.spikeF.BN.rds")
saveRDS(tm.smartseq2.sumF.BN, "data/tm.smartseq2.sumF.BN.rds")

# tm.smartseq2.spikeF.BN <- readRDS("tm.smartseq2.spikeF.BN.rds")
# tm.smartseq2.sumF.BN <- readRDS("tm.smartseq2.sumF.BN.rds")
# these are normalised saved two versions out of which different subsets of cells can be extracted


### 6. Select useful subsets of cells for downsteram ICA etc.

## 6.a Remove vaguely annotated cells NOW

# spikeF
drop.cells <- colData(tm.smartseq2.spikeF.BN)$cell_ontology_class %in% c("blood cell", "myeloid cell", "professional antigen presenting cell", "leukocyte", "lymphocyte")
tm.smartseq2.spikeF.BN <- tm.smartseq2.spikeF.BN[, !drop.cells] # 2023 cells removed 
# sumF
drop.cells <- colData(tm.smartseq2.sumF.BN)$cell_ontology_class %in% c("blood cell", "myeloid cell", "professional antigen presenting cell", "leukocyte", "lymphocyte")
tm.smartseq2.sumF.BN <- tm.smartseq2.sumF.BN[, !drop.cells] # 2023 cells removed 

## 6.b Create a smaller data subset (123) with up to 300 cells per celltype_tissue (600 for microglia)

# get celltype-tissue count column
celltype.tissue <- tibble(cellID = tm.smartseq2.spikeF.BN$cell, 
                          celltype_tissue = tm.smartseq2.spikeF.BN$tissue_celltype) %>% 
  dplyr::group_by(celltype_tissue) %>% 
  dplyr::mutate(n = n()) %>% 
  dplyr::ungroup()

# sample
takeall <- celltype.tissue %>% dplyr::filter(n <= 300) %>% 
  dplyr::select(cellID) %>% pull() # 6420
set.seed(123)
samp300 <- celltype.tissue %>% dplyr::filter(n > 300 & celltype_tissue != "Brain_Myeloid-microglial cell") %>% 
  group_by(celltype_tissue) %>% sample_n(300) %>% ungroup() %>% 
  dplyr::select(cellID) %>% pull() # 11400
set.seed(123)
microgli <- celltype.tissue %>% dplyr::filter(celltype_tissue == "Brain_Myeloid-microglial cell") %>% 
  sample_n(600) %>% ungroup() %>% 
  dplyr::select(cellID) %>% pull() # 600

Reduce(f = intersect, x = c(takeall, samp300, microgli)) # 0

keep.123 <- sort(c(takeall, samp300, microgli))
# spikeF
tm.smartseq2.spikeF.123 <- tm.smartseq2.spikeF.BN[, keep.123]
tibble(tissue_celltype=tm.smartseq2.spikeF.123$tissue_celltype) %>% dplyr::count(tissue_celltype) %>% dplyr::arrange(tissue_celltype) %>% print(n=115)
saveRDS(tm.smartseq2.spikeF.123, "data/tm.smartseq2.spikeF.123.rds")

# sumF
tm.smartseq2.sumF.123 <- tm.smartseq2.sumF.BN[, keep.123]
saveRDS(tm.smartseq2.sumF.123, "data/tm.smartseq2.sumF.123.rds")


## 6.c Select HVG, and save as matrices for ICA (excluding spike-ins)

# spikeF - variability in total RNA across cell types was RETAINED here
keep.hvg.nonspike <- rowData(tm.smartseq2.spikeF.123)$is_hvg
tm.spikeF.123.hvg <- tm.smartseq2.spikeF.123[keep.hvg.nonspike, ]
tm.smartseq2.spikeF.mat.hvg.123 <- logcounts(tm.spikeF.123.hvg) # expand it with as.matrix() when needed
saveRDS(tm.smartseq2.spikeF.mat.hvg.123, "data/tm.smartseq2.spikeF.mat.hvg.123.rds")
# exp.mat.spikeF <- readRDS("data/tm.smartseq2.spikeF.mat.hvg.123.rds")
# exp.mat.spikeF <- as.matrix(exp.mat.spikeF)

# sumF - variability in total RNA across cell types was REMOVED here
keep.HVG.nonspike <- rowData(tm.smartseq2.sumF.123)$is_HVG
tm.sumF.123.hvg <- tm.smartseq2.sumF.123[keep.HVG.nonspike, ]
tm.smartseq2.sumF.mat.hvg.123 <- logcounts(tm.sumF.123.hvg) # expand it with as.matrix() when needed
saveRDS(tm.smartseq2.sumF.mat.hvg.123, "data/tm.smartseq2.sumF.mat.hvg.123.rds")
# exp.mat.sumF <- readRDS("data/tm.smartseq2.sumF.mat.hvg.123.rds")
# exp.mat.sumF <- as.matrix(exp.mat.sumF)


### 7. Denoising the expression values using the PCA - SOMETHING TO  MAYBE TRY ON THE CLUSTER, OTHERWISE NEVER FINISHES ON MY COMP !

# here you need to chose which "technical" trend from above to use: loess or spline (so youre denoising according to the modelled tech trend, to remove technical component) 
# you can make two sets of denoised data - for spike-in and for deconvolution methods - but do only with loess trend 
# primary interest is to  get "lowrank" appoximations of denoised data.. to atempt ica with this 

## 7.a Get denoising for spike-in hvg appoach 
tm.smartseq2.spikeF.123.pca <- denoisePCA(tm.smartseq2.spikeF.123, technical=rowData(tm.smartseq2.spikeF.123)$tech_spike, assay.type="logcounts", 
                                     value="pca", min.rank=5, max.rank=100, get.spikes=FALSE, sce.out=TRUE) # also value="pca"; technical=rowData(tm.smartseq2.spikeF.BN)$tech_spike; technical=comb.spikeF.var.out.loess

tm.smartseq2.spikeF.123.pcaa <- denoisePCA(tm.smartseq2.spikeF.123.pca, technical=rowData(tm.smartseq2.spikeF.123.pca)$tech_spike, assay.type="logcounts", 
                                     value="lowrank", min.rank=5, max.rank=100, get.spikes=FALSE, sce.out=TRUE) # also value="pca"; technical=rowData(tm.smartseq2.spikeF.BN)$tech_spike; technical=comb.spikeF.var.out.loess


## 7.b Get denoising for deconvolution hvg appoach 
tm.smartseq2.sumF.123.pca <- denoisePCA(tm.smartseq2.sumF.123, technical=rowData(tm.smartseq2.sumF.123)$tech_sum, assay.type="logcounts", 
                                   value="pca", min.rank=5, max.rank=100, get.spikes=FALSE, sce.out=TRUE) # also value="pca"; technical=rowData(tm.smartseq2.sumF.BN)$tech_sum; technical=comb.sumF.var.out.loess

tm.smartseq2.sumF.123.pcaa <- denoisePCA(tm.smartseq2.sumF.123.pca, technical=rowData(tm.smartseq2.sumF.BN.pca)$tech_sum, assay.type="logcounts", 
                                   value="lowrank", min.rank=5, max.rank=150, get.spikes=FALSE, sce.out=TRUE) # also value="pca"; technical=rowData(tm.smartseq2.sumF.BN)$tech_sum; technical=comb.sumF.var.out.loess



# next: ica steps - continue to the scripts in the "cluster_scripts" folder 










