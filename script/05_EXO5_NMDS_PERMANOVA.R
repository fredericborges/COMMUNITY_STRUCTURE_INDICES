#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# EXAMPLE 5 : BETA DIVERSITY - NMDS, PERMANOVA -----
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Frédéric Borges
# frederic.borges@univ-lorraine.fr
# date: 2023/05/15


# NMDS stands for Non Metric Multidimensional Scaling. 
# To understand how it works, I suggest this great tutorial from the Rob K statistics chanel:
# https://www.youtube.com/watch?v=YS89_2n9Tlw&t=515s

# Clear the environment
rm(list=ls())
# LIBRARIES----
library(tidyverse)  # manipulation de données
library(phyloseq)   # analyse de microbiome
library(ape)        # manipulation d'arbre
library(vegan)      # analyse écologique de communauté 
library(ggtree)     # tree
library(ggrepel)    # to label the points in the scatter plots
source(file.path("script","graphical_methods.R"))
# DATA PREPARATION----
# The objective is to format the data for the following analysis. 
# Declare the variables that will be used in the script "data_preparation.R" 
abundance_table <- "exo5_16S_v2.txt"
metadata_file <- "exo5_metadata_16S.txt"
name_phyloseq_file <- "phyloseq.Rdata"
# script to prepare the data, including the creation of the phyloseq object ps
source(file.path("script","data_preparation.R"), echo = T)
# Phyloseq object
# Look at the phyloseq object
ps
# DICTIONNARY, COLOR----
## We create a "dictionary" for translation and order the categories 
## as we want 
dictionary = c("low","high")
# colors
my_palette <- c('#67001f','#4393c3')
names(my_palette) <- dictionary
# PLOT COMPOSITION----
abundance_table %>% select(Phylum,Genus,sam1,sam2,sam3,sam4,sam5,sam6, sam7, sam8)
p1 <- plot_composition(ps, "Kingdom", "Bacteria", "Genus", fill = "Genus") + 
  facet_grid(~pH, scales = "free_x", space = "free_x")+
  theme_bw()
p1 # The barplots show slight differences between samples according to the variable pH. 

# ALPHA DIVRESITY----
sample_sums(ps) %>% range() # range of read numbers in the overall samples
ps_rare <- rarefy_even_depth(ps, rngseed = 20190124) # rarefaction
sample_sums(ps_rare) %>% range() # range after rarefaction
# plot de alpha diversity indices
plot_richness(ps_rare, x = "pH", color = "pH",
              measures = c("Observed", "Shannon", "InvSimpson")) + 
  geom_boxplot(aes(group = pH))+## add one boxplot per type of food
  scale_color_manual(values = my_palette)+
  theme_bw()
#Anova
div_data <- cbind(estimate_richness(ps_rare),  ## diversity indices
                  sample_data(ps_rare)         ## covariates
)

# Can you interprete the anova analysis ?
model <- aov(Shannon ~ pH, data = div_data)
anova(model)

model <- aov(InvSimpson ~ pH, data = div_data)
anova(model)

# BETA DIVERSITY----
# Compute the distances 
dist <- phyloseq::distance(ps_rare, method = "bray") # important : don't use distance() but phyloseq::distance() or you may get an error
# NMDS
set.seed(13)
nmds <-metaMDS(dist,k=2, trymax = 100) # k = 2 dimensions
nmds$stress # check that the stress value is below 0.2
stressplot(nmds) # the data should follow a line 
# GRAPHICAL REPRESENTATION
metad2 <- metad %>% 
  rownames_to_column("samples")
scores(nmds) %>% 
  as_tibble(rownames= "samples") %>% 
  inner_join(., metad2, by = "samples") %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = pH))+
  geom_point(size= 4)+
  scale_color_manual(values = my_palette)+
  geom_label_repel(aes(label = rownames(scores(nmds))), 
                  max.overlaps = 100,
                   size= 2.8,
                   min.segment.length = 0,
                   box.padding = 0.5,
                   label.padding = 0.15,
                   label.r = 0.13,
                   label.size =0.4)+
  theme_bw()
# The samples are well grouped according to pH level. 
# HIERARCHICAL CLUSTERING
plot_clust(ps_rare, dist = "bray", method = "ward.D2", color = "pH",
           palette = my_palette, 
           title = "Clustering of samples (Bray Curtis + Ward)\nsamples colored by pH")
# The samples are well grouped according to pH level. 
# PERMANOVA
# Hypothesis : the pH and beta-diversity are linked.
metadata <- sample_data(ps) %>% as("data.frame")
model <- vegan::adonis(dist ~ pH, data = metadata, permutations = 999)
model$aov.tab

# Exercice : Change the code in order to see the link between pH and beta diversity by using the weighted-unifrac index; 
# Hint : be careful, you need to change the code in several places.


# REFERENCES
sessionInfo()


