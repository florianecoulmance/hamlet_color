#!/usr/bin/env Rscript
# by: Floriane Coulmance: 24/02/2021
# usage:
# Rscript --vanilla gxp_plots.R <data_path> <figure_path>
# -------------------------------------------------------------------------------------------------------------------
# data_path in : $BASE_DIR/outputs/7_gxp/$DATASET
# figure_path in : $BASE_DIR/outputs/7_gxp/$DATASET/figures/
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed library
library(ggplot2)
library(tidyverse)
library(hypogen)
library(dplyr)
library(plyr)


# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[7:8]
print(args)

data_path <- as.character(args[1]) # Path to mean coverage data table
# data_path <- "/Users/fco/Desktop/BREMEN_OP/ibd/coverage_table"
figure_path <- as.character(args[2]) # Path to the figure folder


setwd('/user/doau0129/chapter1_2/outputs/gxp')

# df_gxp_snout <- read.table(file = 'snout.lm.50k.5k.txt.gz', sep = '\t', header = TRUE)
# df_gxp_bars <- read.table(file = 'bars_body.lm.50k.5k.txt.gz', sep = '\t', header = TRUE)
# 

files <- list.files(pattern = "lmm.50k.5k.txt.gz")
traits <- list("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

for(i in 1:length(traits)) {
  
  pdf(paste0(figure_path,"plot_",traits[i],".pdf",spe=""))
  
  par(mfrow=c(2,2))
  
  for(j in 1:length(files)){
    old_trait <- traitfile
    
    traitfile <- regmatches(x=files[j],gregexpr('[.]',files[j]),invert=TRUE)[[1]][1]
    print(traitfile)
      
    }
    
    if(traitfile==traits[i]) {
      data <- read_tsv(files[j]) %>% left_join(hypo_chrom_start) %>% mutate(GPOS = MID_POS + GSTART)
      name <- tools::file_path_sans_ext(files[j])
      name <- tools::file_path_sans_ext(name)
      print(name)
      
      
      p <- ggplot(data = data, aes(x = GPOS, y = AVG_p_wald))+
        geom_hypo_LG()+
        geom_point(size = .2)+
        scale_fill_hypo_LG_bg()+
        scale_x_hypo_LG()+
        theme_hypo() + ggtitle(name)
      
      print(p)
       
    }
  
  }
  dev.off()
  
}


