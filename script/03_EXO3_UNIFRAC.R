#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# EXAMPLE 3 : BETA DIVERSITY - UNIFRAC----
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Frédéric Borges
# frederic.borges@univ-lorraine.fr
# date: 2023/05/15


# Clear the environment
rm(list=ls())
# LIBRARIES----
library(tidyverse)  # manipulation de données
library(phyloseq)   # analyse de microbiome
library(ape)        # manipulation d'arbre
library(vegan)      # analyse écologique de communauté 
library(ggtree)     # tree
source(file.path("script","graphical_methods.R"))
# DATA PREPARATION----
# The objective is to format the data for the following analysis. 
# Declare the variables that will be used in the script "data_preparation.R" 
abundance_table <- "exo3_16S.txt"
metadata_file <- "exo3_metadata_16S.txt"
name_phyloseq_file <- "phyloseq.Rdata"
# script to prepare the data, including the creation of the phyloseq object ps
source(file.path("script","data_preparation.R"), echo = T)
# Phyloseq object
# Look at the phyloseq object
ps
# DICTIONNARY, COLOR----
## We create a "dictionary" for translation and order the categories 
## as we want 
dictionary = c("s1","s2","s3")
# colors
my_palette <- c('#67001f','#f4a582','#4393c3')
names(my_palette) <- dictionary
# ABUNDANCE TABLE----
abundance_table %>% select(Phylum,Genus,sam1,sam2,sam3)
# PLOT COMPOSITION----
p1 <- plot_composition(ps, "Kingdom", "Bacteria", "Genus", fill = "Genus") + 
  facet_grid(~Step, scales = "free_x", space = "free_x")+
  theme_bw()
# ALPHA DIVRESITY----
p2 <- plot_richness(ps, x = "Step", color = "Step",
              measures = c("Observed", "Shannon", "InvSimpson")) + 
  geom_point(aes(group = Step))+## add one boxplot per type of food
  scale_color_manual(values = my_palette)+
  theme_bw()

gridExtra::grid.arrange(p1,p2, 
                        ncol = 1, nrow = 2)

# BETA DIVERSITY----
# Distance 
# Jaccard
dist <- phyloseq::distance(ps, method = "cc")
dist # see the distance matrix in the console
# Bray-Curtis
dist <- phyloseq::distance(ps, method = "bray")
dist # see the distance matrix in the console
# Show tree
ggtree(ps) + 
  geom_nodelab(aes(label=label), hjust=-0.05, size=3.5)+ # tree scaffold
  geom_point(aes(x=x+hjust, color=Step, shape=Phylum, 
                 size=Abundance), na.rm=TRUE)+ # add the symbols wich size depends on abundance
  geom_tiplab(aes(label=Genus), hjust=-0.35)+ # add the name of the OTUs
  scale_size_continuous(trans=log_trans(5))+
  scale_color_manual(values = my_palette)+
  theme(legend.position="right") + hexpand(.4)


# Unifrac
dist <- phyloseq::distance(ps, method = "unifrac")
dist # see the distance matrix in the console
abundance_table %>% select(Phylum,Genus,sam1,sam2,sam3)

# REFERENCES
sessionInfo()
