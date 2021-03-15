#!/usr/bin/env Rscript
# by: Floriane Coulmance: 24/02/2021
# usage:
# Rscript --vanilla gxp_plots.R <data_path> <figure_path> <analysis>
# -------------------------------------------------------------------------------------------------------------------
# data_path in : $BASE_DIR/outputs/7_gxp/$DATASET
# figure_path in : $BASE_DIR/outputs/7_gxp/$DATASET/figures/
# analysis in : "univariate_gemma", "multivariate_plink"
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

# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6:8]
print(args)

data_path <- as.character(args[1]) # Path to mean coverage data table
#data_path <- "/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/chapter1/"
figure_path <- as.character(args[2]) # Path to the figure folder
#figure_path <- "/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/chapter1/"
#analysis <- "univariate_gemma"
analysis <- as.character(args[3])


# -------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------
plotadgwas <- function(dataset,path,analysis,y) {
  p <- ggplot() + facet_wrap(RUN~., ncol = 1,dir = 'v', strip.position="right") +
    geom_hypo_LG() +
    geom_point(data = dataset, aes(x = GPOS, y = y), size = .1) +
    scale_fill_hypo_LG_bg() +
    scale_x_hypo_LG(name = "Linkage Groups") +
    scale_y_continuous(name = expression(italic('-log10(p_value)'))) +
    theme_hypo() +
    theme(legend.position = 'none',
          axis.title.x = element_text(),
          axis.text.x.top= element_text(colour = 'darkgray'))
  
  hypo_save(filename = paste0(path,analysis,".png"),
            plot = p,
            width = 8,
            height = 8)
}


concat_files_gem <- function(f,p) {
  count = 0
  for (file in f) {
    count = count + 1 
    run_files <- file %>%
      str_sub(.,end=-19) %>%
      str_replace(.,pattern = '([a-z]{3})-([a-z]{3})-([a-z]{3})', '\\2\\1-\\3\\1')
    print(run_files)
    d <- read.table(paste0(p,file), header=TRUE)
    d$RUN <- run_files
    assign(run_files,d)
    print(head(d))
  }
}

concat_files_plk <- function(f) {
  count = 0
  for (file in f) {
    count = count + 1 
    run_files <- file %>%
      str_sub(.,end=-23) %>%
      str_replace(.,pattern = '([a-z]{3})-([a-z]{3})-([a-z]{3})', '\\2\\1-\\3\\1')
    print(run_files)
    d <- read.table(file, header=TRUE)
    d$RUN <- run_files
    assign(run_files,d)
  }
}
# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

#setwd(data_path)

if (analysis == "univariate_gemma"){
  files <- list.files(data_path, pattern = "lmm.50k.5k.txt.gz")
  print(files)
  traits <- list("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
  print(traits)
  data <- concat_files_gem(files,data_path)
  print("here1")
  print(data)
  print(head(PC1))
  print("here2")
  files_l <- list(PC1,PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)
  print(files_l)
  print("here3")
  
} else {
  
  f <- list.files(data_path, pattern = "mvplink.50k.5k.txt.gz")
  
  if (analysis == "multivariate_plink_PC1") {
    files <- list.files(f, pattern = "PC1-")
    traits <- list("PC1", "PC1-2", "PC1-3", "PC1-4", "PC1-5", "PC1-6", "PC1-7", "PC1-8", "PC1-9", "PC1-10")
    concat_files_plk(files)
    print("here1")
    print(head(PC1))
    print("here2")
    files_l <- list(PC1, PC1_2, PC1_3, PC1_4, PC1_5, PC1_6, PC1_7, PC1_8, PC1_9, PC1_10)
    print("here3")

  } else if (analysis == "multivariate_plink_PC2"){
    files <- list.files(f, pattern = "PC2-")
    traits <- list("PC2", "PC2-3", "PC2-4", "PC2-5", "PC2-6", "PC2-7", "PC2-8", "PC2-9", "PC2-10")
    concat_files_plk(files)
    files_l <- list(PC2, PC2_3, PC2_4, PC2_5, PC2_6, PC2_7, PC2_8, PC2_9, PC2_10)
    
  } else if (analysis == "multivariate_plink_PC3"){
    files <- list.files(f, pattern = "PC3-")
    traits <- list("PC3", "PC3-4", "PC3-5", "PC3-6", "PC3-7", "PC3-8", "PC3-9", "PC3-10")
    concat_files_plk(files)
    files_l <- list(PC3, PC3_4, PC3_5, PC3_6, PC3_7, PC3_8, PC3_9, PC3_10)
    
  } else if (analysis == "multivariate_plink_PC4"){
    files <- list.files(f, pattern = "PC4-")
    traits <- list("PC4", "PC4-5", "PC4-6", "PC4-7", "PC4-8", "PC4-9", "PC4-10")
    concat_files_plk(files)
    files_l <- list(PC4, PC4_5, PC4_6, PC4_7, PC4_8, PC4_9, PC4_10)
    
  } else if (analysis == "multivariate_plink_PC5"){
    files <- list.files(f, pattern = "PC5-")
    traits <- list("PC5", "PC5-6", "PC5-7", "PC5-8", "PC5-9", "PC5-10")
    concat_files_plk(files)
    files_l <- list(PC5, PC5_6, PC5_7, PC5_8, PC5_9, PC5_10)
    
  } else if (analysis == "multivariate_plink_PC6"){
    files <- list.files(f, pattern = "PC6-")
    traits <- list("PC6", "PC6-7", "PC6-8", "PC6-9", "PC6-10")
    concat_files_plk(files)
    files_l <- list(PC6, PC6_7, PC6_8, PC6_9, PC6_10)
    
  } else if (analysis == "multivariate_plink_PC7"){
    files <- list.files(f, pattern = "PC7-")
    traits <- list("PC7", "PC7-8", "PC7-9", "PC7-10")
    concat_files_plk(files)
    files_l <- list(PC7, PC7_8, PC7_9, PC7_10)
    
  } else if (analysis == "multivariate_plink_PC8"){
    files <- list.files(f, pattern = "PC8-")
    traits <- list("PC8", "PC8-9", "PC8-10")
    concat_files_plk(files)
    files_l <- list(PC8, PC8_9, PC8_10)
    
    
  } else if (analysis == "multivariate_plink_PC9"){
    files <- list.files(f, pattern = "PC9-")
    traits <- list("PC9","PC9-10")
    concat_files_plk(files)
    files_l <- list(PC9, PC9_10)
    
  } else if (analysis == "multivariate_plink_byPCs"){
    files <- list.files("PC1.mvplink.50k.5k.txt.gz", "PC2.mvplink.50k.5k.txt.gz", "PC3.mvplink.50k.5k.txt.gz", "PC4.mvplink.50k.5k.txt.gz", "PC5.mvplink.50k.5k.txt.gz", "PC6.mvplink.50k.5k.txt.gz", "PC7.mvplink.50k.5k.txt.gz", "PC8.mvplink.50k.5k.txt.gz", "PC9.mvplink.50k.5k.txt.gz", "PC10.mvplink.50k.5k.txt.gz")
    traits <- list("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
    concat_files_plk(files)
    files_l <- list(PC1,PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)
    
  }
  
}

names(files_l) <- traits
print(files_l)
files <- bind_rows(files_l, .id = 'id') %>% left_join(hypo_chrom_start) %>% mutate(GPOS = MID_POS + GSTART)

if (analysis == "univariate_gemma"){
  plotadgwas(files,figure_path,analysis,"AVG_p_wald")
  
} else {
  plotadgwas(files,figure_path,analysis,"AVG_P")  
}
  




