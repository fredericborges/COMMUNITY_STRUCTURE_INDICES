# Author: Frédéric Borges
# frederic.borges@univ-lorraine.fr
# date: 2023/05/15

# VARIABLES----
path2abundance <- file.path("data",abundance_table)      # path to the adundance table file
path2metad <- file.path("data",metadata_file)  # path to the metadata file
col_name_sample <- "Sample_name"

# LIRBARIES----
library(tidyverse)
library(dada2)
library(Biostrings)
library(DECIPHER)
library(phangorn)
library(phyloseq)

df <- read.delim(file = file.path(path2abundance),
                 check.names = F, 
                 dec = ",")

# SEQTAB.NOCHIM----
seqtab.nochim <- df %>% 
  select(seed_sequence,12:ncol(df)) %>% 
  column_to_rownames("seed_sequence") %>% 
  as.matrix() %>% 
  t()
# TAXA----
taxa <- df %>% 
  select(blast_taxonomy, seed_sequence)
ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family","Genus","Species")
taxa[ranks] <-str_split_fixed(taxa$blast_taxonomy, ';', 7) 
taxa <- taxa %>% 
  select(-1) %>% 
  column_to_rownames("seed_sequence") %>% 
  as.matrix()
# METADATA----
metad <- read.delim(file = path2metad, check.names = F)
metad <- metad %>% 
  column_to_rownames(col_name_sample)
# TREE----
sequences <- getSequences(rownames(taxa))
names(sequences) <- sequences 
# this propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)
# Then we build a neighbour-joining tree then fit a maximum likelihood 
# tree using the neighbour-joining tree as a starting point
phang_align <- phyDat(as(alignment, 'matrix'), type='DNA')
dm <- dist.ml(phang_align)
treeNJ <- NJ(dm)  
# note, tip order != sequence order
fit = pml(treeNJ, data=phang_align)
## negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
# l’étape qui suit requiert un temps de calcul important
fitGTR <- optim.pml(fitGTR, model='GTR', optInv=TRUE, optGamma=TRUE,
                    rearrangement = 'stochastic',
                    control = pml.control(trace = 0))
tree <- fitGTR$tree
# PHYLOSEQ OBJECT----
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = F),
  tax_table(taxa), 
  sample_data(metad),
  phy_tree(fitGTR$tree)
)

abundance_table <- cbind(t(seqtab.nochim), taxa)
abundance_table <- as.data.frame(abundance_table)
abundance_table <- tibble::rownames_to_column(abundance_table, "Sequence")

save(ps, file = file.path("output", name_phyloseq_file))
rm(phang_align,tree,treeNJ,df,alignment,fit,col_name_sample,dm,name_phyloseq_file,path2abundance, path2metad, ranks,sequences)

