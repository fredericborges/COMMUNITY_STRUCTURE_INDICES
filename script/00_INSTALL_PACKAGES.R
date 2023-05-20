# Installation if required

pack <- c("BiocManager",
          "tidyverse",
          "ape",
          "vegan",
          "phyloseq",
          "ggtree")

for(i in 1:length(pack)){
  if(!require(pack[i], character.only = TRUE)){
    install.packages(pack[i], character.only = TRUE)
    library(pack[i], character.only = TRUE)
  }else{
    message(paste(pack[i],"is already installed."))
  }
}

