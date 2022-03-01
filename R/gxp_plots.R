#!/usr/bin/env Rscript
# by: Floriane Coulmance: 24/02/2021
# usage:
# Rscript gxp_plots.R <data_path> <figure_path> <analysis>
# -------------------------------------------------------------------------------------------------------------------
# data_path in : $BASE_DIR/outputs/7_gxp/$DATASET
# figure_path in : $BASE_DIR/outputs/7_gxp/figures/$DATASET/
# analysis in : "univariate_gemma", "multivariate_plink"
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
args = args[6:8]
print(args)

data_path <- as.character(args[1]) # Path to mean coverage data table
#data_path <- "/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/chapter1/"
figure_path <- as.character(args[2]) # Path to the figure folder
#figure_path <- "/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/chapter1/"
#analysis <- "multivariate_plink_PC1"
analysis <- as.character(args[3])


# -------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


plotgwas_gem <- function(dataset,path,analysis) {
  p <- ggplot() + facet_wrap(RUN~., ncol = 1, dir = 'v', strip.position="right") +
    geom_hypo_LG() +
    geom_point(data = dataset, aes(x = GPOS, y = AVG_p_wald), size = .1) +
    scale_fill_hypo_LG_bg() +
    scale_x_hypo_LG(name = "Linkage Groups") +
    scale_y_continuous(name = expression(italic('-log(p-Wald)'))) +
    theme_hypo() +
    theme(legend.position = 'none',
          axis.title.x = element_text(),
          axis.text.x.top= element_text(colour = 'darkgray'))
  
  hypo_save(filename = paste0(path,analysis,"_lmm_50k5k.png"), type="cairo",
            plot = p,
            width = 8,
            height = 8)
}


plotgwas_plk <- function(dataset,path,analysis) {
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
  
  hypo_save(filename = paste0(path,analysis,"_50k5k.png"), type="cairo",
            plot = p,
            width = 8,
            height = 8)
}


plotgwas_mvp <- function(dataset,path,analysis) {
  p <- ggplot() + facet_wrap(RUN~., ncol = 1, dir = 'v', strip.position="right") +
    geom_hypo_LG() +
    geom_point(data = dataset, aes(x = GPOS, y = AVG_P), size = .1) +
    scale_fill_hypo_LG_bg() +
    scale_x_hypo_LG(name = "Linkage Groups") +
    scale_y_continuous(name = expression(italic('-log(p-value)'))) +
    theme_hypo() +
    theme(legend.position = 'none',
          axis.title.x = element_text(),
          axis.text.x.top= element_text(colour = 'darkgray'))
  
  
  hypo_save(filename = paste0(path,analysis,".png"), type="cairo",
            plot = p,
            width = 8,
            height = 8)
}


concat_files_gem <- function(f,p) {
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


concat_files_plk <- function(f,p) {
  count = 0
  l <- list()
  for (file in f) {
    count = count + 1 
    run_files <- file %>%
      str_sub(.,end=-21) %>%
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


concat_files_mul <- function(f,p) {
  count = 0
  l <- list()
  for (file in f) {
    count = count + 1 
    run_files <- file %>%
      str_sub(.,end=-23) %>%
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


if (analysis == "univariate_gemma") {

  files <- list.files(data_path, pattern = "lmm.50k.5k.txt.gz")
  print(files)
  traits <- list("PC1", "PC10", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9")
  print(traits)
  files_l <- concat_files_gem(files,data_path)
  # print("here1")
  # print(head(PC1))
  # print("here2")
  # files_l <- list(PC1,PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)
  # print(files_l)
  # print("here3")
  
} else if (analysis == "univariate_plink") {

  files <- list.files(data_path, pattern = "assoc.50k.5k.txt.gz")
  print(files)
  traits <- list("PC1", "PC10", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9")
  print(traits)
  files_l <- concat_files_plk(files,data_path)

} else {  
  
  f <- list.files(data_path, pattern = "mvplink.50k.5k.txt.gz")
  print(f)  
  if (analysis == "multivariate_plink_PC1") {
    files <- grep("PC1_", f, value = TRUE)
    print(files)
    traits <- list("PC1_10", "PC1_2", "PC1_3", "PC1_4", "PC1_5", "PC1_6", "PC1_7", "PC1_8", "PC1_9") #"PC1_11", "PC1_12", "PC1_13", "PC1_14", "PC1_15", 
    print(traits)
    files_l <- concat_files_mul(files,data_path)

  } else if (analysis == "multivariate_plink_PC2"){
    files <- grep("PC2_", f, value = TRUE)
    print(files)
    traits <- list("PC2_10", "PC2_3", "PC2_4", "PC2_5", "PC2_6", "PC2_7", "PC2_8", "PC2_9") #"PC2_11", "PC2_12", "PC2_13", "PC2_14", "PC2_15", 
    print(traits)
    files_l <- concat_files_mul(files,data_path)
    
  } else if (analysis == "multivariate_plink_PC3"){
    files <- grep("PC3_", f, value = TRUE)
    print(files)
    traits <- list("PC3_10", "PC3_4", "PC3_5", "PC3_6", "PC3_7", "PC3_8", "PC3_9") #"PC3_11", "PC3_12", "PC3_13", "PC3_14", "PC3_15", 
    print(traits)
    files_l <- concat_files_mul(files,data_path)
    
  } else if (analysis == "multivariate_plink_PC4"){
    files <- grep("PC4_", f, value = TRUE)
    print(files)
    traits <- list("PC4_10", "PC4_5", "PC4_6", "PC4_7", "PC4_8", "PC4_9") #"PC4_11", "PC4_12", "PC4_13", "PC4_14", "PC4_15", 
    print(traits)
    files_l <- concat_files_mul(files,data_path)
    
  } else if (analysis == "multivariate_plink_PC5"){
    files <- grep("PC5_", f, value = TRUE)
    print(files)
    traits <- list("PC5_10", "PC5_6", "PC5_7", "PC5_8", "PC5_9") #"PC5_11", "PC5_12", "PC5_13", "PC5_14", "PC5_15", 
    print(traits)
    files_l <- concat_files_mul(files, data_path)
    
  } else if (analysis == "multivariate_plink_PC6"){
    files <- grep("PC6_", f, value = TRUE)
    print(files)
    traits <- list("PC6_10", "PC6_7", "PC6_8", "PC6_9") #"PC6_11", "PC6_12", "PC6_13", "PC6_14", "PC6_15", 
    print(traits)
    files_l <- concat_files_mul(files, data_path)
    
  } else if (analysis == "multivariate_plink_PC7"){
    files <- grep("PC7_", f, value = TRUE)
    print(files)
    traits <- list("PC7_10", "PC7_8", "PC7_9") #"PC7_11", "PC7_12", "PC7_13", "PC7_14", "PC7_15", 
    print(traits)
    files_l <- concat_files_mul(files, data_path)
    
  } else if (analysis == "multivariate_plink_PC8"){
    files <- grep("PC8_", f, value = TRUE)
    print(files)
    traits <- list("PC8_10", "PC8_9") #"PC8_11", "PC8_12", "PC8_13", "PC8_14", "PC8_15",
    print(traits)
    files_l <- concat_files_mul(files, data_path)
     
  } else if (analysis == "multivariate_plink_PC9"){
    files <- grep("PC9_", f, value = TRUE)
    print(files)
    traits <- list("PC9_10") #, "PC9_11", "PC9_12", "PC9_13", "PC9_14", "PC9_15")
    print(traits)
    files_l <- concat_files_mul(files, data_path)

  # } else if (analysis == "multivariate_plink_PC10"){
  #   files <- grep("PC10_", f, value = TRUE)
  #   print(files)
  #   traits <- list("PC10_11", "PC10_12", "PC10_13", "PC10_14", "PC10_15")
  #   print(traits)
  #   files_l <- concat_files_plk(files, data_path)

  # } else if (analysis == "multivariate_plink_PC11"){
  #   files <- grep("PC11_", f, value = TRUE)
  #   print(files)
  #   traits <- list("PC11_12", "PC11_13", "PC11_14", "PC11_15")
  #   print(traits)
  #   files_l <- concat_files_plk(files, data_path)

  # } else if (analysis == "multivariate_plink_PC12"){
  #   files <- grep("PC12_", f, value = TRUE)
  #   print(files)
  #   traits <- list("PC12_13", "PC12_14", "PC12_15")
  #   print(traits)
  #   files_l <- concat_files_plk(files, data_path)

  # } else if (analysis == "multivariate_plink_PC13"){
  #   files <- grep("PC13_", f, value = TRUE)
  #   print(files)
  #   traits <- list("PC13_14", "PC13_15")
  #   print(traits)
  #   files_l <- concat_files_plk(files, data_path)

  # } else if (analysis == "multivariate_plink_PC14"){
  #   files <- grep("PC14_", f, value = TRUE)
  #   print(files)
  #   traits <- list("PC14_15")
  #   print(traits)
  #   files_l <- concat_files_plk(files, data_path)
    
  } else if (analysis == "multivariate_plink_byPCs"){
    files <- list("PC1.mvplink.50k.5k.txt.gz", "PC2.mvplink.50k.5k.txt.gz", "PC3.mvplink.50k.5k.txt.gz", "PC4.mvplink.50k.5k.txt.gz", "PC5.mvplink.50k.5k.txt.gz", "PC6.mvplink.50k.5k.txt.gz", "PC7.mvplink.50k.5k.txt.gz", "PC8.mvplink.50k.5k.txt.gz", "PC9.mvplink.50k.5k.txt.gz", "PC10.mvplink.50k.5k.txt.gz") #, "PC11.mvplink.50k.5k.txt.gz", "PC12.mvplink.50k.5k.txt.gz", "PC13.mvplink.50k.5k.txt.gz", "PC14.mvplink.50k.5k.txt.gz", "PC15.mvplink.50k.5k.txt.gz"
    print(files)
    traits <- list("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10") #, "PC11", "PC12", "PC13", "PC14", "PC15"
    print(traits)
    files_l <- concat_files_mul(files, data_path)
    
  }
  
}


names(files_l) <- traits
files <- bind_rows(files_l, .id = 'id') %>% left_join(hypo_chrom_start) %>% mutate(GPOS = MID_POS + GSTART)
files$range <- do.call(paste, c(files[c("CHROM", "BIN_START", "BIN_END")], sep="_"))


if (analysis == "univariate_gemma"){
  plotgwas_gem(files,figure_path,analysis)

} else if (analysis == "univariate_plink") {
  plotgwas_plk(files,figure_path,analysis)
  
} else {
  plotgwas_mvp(files,figure_path,analysis)  
}





