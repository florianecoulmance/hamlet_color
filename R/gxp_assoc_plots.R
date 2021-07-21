#!/usr/bin/env Rscript
# by: Floriane Coulmance: 24/02/2021
# usage:
# Rscript gxp_plots.R <data_path> <figure_path>
# -------------------------------------------------------------------------------------------------------------------
# data_path in : $BASE_DIR/outputs/7_gxp/$DATASET
# figure_path in : $BASE_DIR/outputs/7_gxp/$DATASET/figures/
# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed library
library(ggplot2)
library(tidyverse)
library(hypogen)
library(hypoimg)
library(dplyr)
library(plyr)


# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6:7]
print(args)

data_path <- as.character(args[1]) # Path to mean coverage data table
#data_path <- "/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/chapter1/"
figure_path <- as.character(args[2]) # Path to the figure folder
#figure_path <- "/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/chapter1/"

# -------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


plotgwas <- function(dataset,path) {
  p <- ggplot() + facet_wrap(RUN~., ncol = 1, dir = 'v', strip.position="right") +
    geom_hypo_LG() +
    geom_point(data = dataset, aes(x = GPOS, y = AVG_P), size = .1) +
    scale_fill_hypo_LG_bg() +
    scale_x_hypo_LG(name = "Linkage Groups") +
    scale_y_continuous(name = expression(italic('-log(p_value)'))) +
    theme_hypo() +
    theme(legend.position = 'none',
          axis.title.x = element_text(),
          axis.text.x.top= element_text(colour = 'darkgray'))
  
  hypo_save(filename = paste0(path,"assoc50k_plink_plots.png"), type="cairo",
            plot = p,
            width = 8,
            height = 8)
}


concat_files <- function(f,p) {
  count = 0
  l <- list()
  for (file in f) {
    count = count + 1 
    run_files <- file %>%
      str_sub(.,end=-19) %>%
      str_replace(.,pattern = '([a-z]{3})-([a-z]{3})-([a-z]{3})', '\\2\\1-\\3\\1')
    print(run_files)
    d <- read.table(paste0(p,file), header=TRUE)
    d$RUN <- run_files
    #assign(run_files,d)
    l[[count]]=assign(run_files,d)
    print(head(d))
    #print(l)
  }
  return(l)
}



# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


files <- list.files(data_path, pattern = "assoc.50k.5k.txt.gz")
print(files)
traits <- list("PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9")
print(traits)
files_l <- concat_files(files,data_path)
names(files_l) <- traits
files <- bind_rows(files_l, .id = 'id') %>% left_join(hypo_chrom_start) %>% mutate(GPOS = MID_POS + GSTART)
files$range <- do.call(paste, c(files[c("CHROM", "BIN_START", "BIN_END")], sep="_"))
plotgwas(files,figure_path)




