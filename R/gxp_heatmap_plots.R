#!/usr/bin/env Rscript
# by: Floriane Coulmance: 01/03/2022
# usage:
# Rscript gxp_heatmap_plots.R <data_path> <figure_path> <dataset> <metadata> <im_path> <1stPC> <2ndPC>
# -------------------------------------------------------------------------------------------------------------------
# data_path in : $BASE_DIR/outputs/7_gxp/$DATASET
# figure_path in : $BASE_DIR/figures/7_gxp/figures/$TYPE/$COLOR_SPACE/$DATASET
# dataset in : LAB_fullm_54off_59on
# metadata in : $BASE_DIR/metadata/
# im_path in : $BASE_DIR/images/$TYPE/$COLOR_SPACE/$DATASET/
# 1stPC : in 5 possible PC
# 2ndPC : in 5 possible PC
# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed library
library(ggplot2)
library(patchwork)
library(ggtext)
library(tidyverse)
library(hypogen)
library(hypoimg)
library(dplyr)
library(plyr)
library(stringr)
library(stringi)
library(png)
library(ggpubr)
library(scales)


# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6:10]
print(args)

data_path <- as.character(args[1]) # Path to mean coverage data table
print(data_path)
figure_path <- as.character(args[2]) # Path to the figure folder
print(figure_path)
dataset <- as.character(args[3])
print(dataset)
metadata <- as.character(args[4])
print(metadata)
im_path <- as.character(args[5])
print(im_path)

# data_path <- "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/"
# figure_path <- "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/"
# dataset <- "LAB_fullm_left_54off_59on"
# metadata <- "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/metadata/"
# im_path <- "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/images/continuous/LAB/LAB_fullm_left_54off_59on/"


# -------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


concat_files_gem <- function(f,p) {
  count = 0
  l <- list()
  for (file in f) {
    count = count + 1 
    run_files <- file %>% 
                 str_match(., "[.]\\s*(.*?)\\s*.50k") %>%
                 .[,2]
    print(run_files)
    d <- read.table(paste0(p,file), header=TRUE)
    d$RUN <- run_files
    print(colnames(d))
    if ("AVG_p_wald" %in% colnames(d)) {
      #d %>% rename(AVG_p_wald = AVG_P)
      names(d)[names(d) == "AVG_p_wald"] <- "AVG_P"
    }
    #assign(run_files,d)
    l[[count]]=assign(run_files,d)
    print(head(d))
    #print(l)
  }
  
  return(l)
  
}


plot_g_h <- function(table,path,prefix,tr) {
  p <- ggplot() + facet_wrap(RUN~., ncol = 1, dir = 'v', strip.position="right") +
    geom_hypo_LG() +
    geom_point(data = table, aes(x = GPOS, y = AVG_P), size = .001) +
    scale_fill_hypo_LG_bg() +
    scale_x_hypo_LG(name = "Linkage Groups") +
    scale_y_continuous(name = expression(italic('-log(p-value)'))) +
    theme_hypo() +
    theme(legend.position = 'none',
          axis.title.x = element_text(),
          axis.text.x.top= element_text(colour = 'darkgray', size = 8),
          plot.margin = unit(c(0,0,0,0), "cm")) 
    

    # im <- paste0(figure_path+dataset+"_"+trait+".png")
    img <- readPNG(paste0(path,prefix,"_",tr,".png"))
    g <- rasterGrob(img, interpolate=TRUE)

    plot <- ggarrange(p, g, ncol = 2, nrow = 1, widths = c(1.5, 1), heights = c(0.5, 1.5)) + ggtitle(tr) + theme(plot.title = element_text(hjust = 0.8))
    
    return(plot) 

}


univariate_plots <- function(path_univariate, trait_list, fig_path, dat) {
  
  # Function to create a dataframe with the 3 univariate analyses 
  # per PCs of the phenotype PCA
  
  count = 0
  print(count)
  
  l <- list()
  print(l)
  
  for (trait in trait_list) {
    print(trait)
    count = count + 1
    print(count)
    string <- paste0(trait,".lmm.50k.5k.txt.gz|",trait,".assoc.50k.5k.txt.gz|",trait,".mvplink.50k.5k.txt.gz")
    print(string)
    f <- list.files(path_univariate, pattern = string)
    print(head(f))
    # f <- f[grepl(paste0(trait,"."), names(f))]  
    # print(f)
    model <- list("PLINK", "GEMMA", "MV PLINK")
    print(model)
    files_l <- concat_files_gem(f,path_univariate)
    names(files_l) <- model
    #print(head(files_l))
    
    f <- bind_rows(files_l, .id = 'id') %>% left_join(hypo_chrom_start) %>% mutate(GPOS = MID_POS + GSTART)
    f$range <- do.call(paste, c(f[c("CHROM", "BIN_START", "BIN_END")], sep="_"))
    #print(head(f))
    
    l[[count]]=assign(trait,f)
  }
  
  # Create set of plots for each phenotype PCs and arrange them in 1 plot
  
  p <- vector('list', 10)
  print(p)
  
  count = 0
  print(count)
  
  for (d in l) {
    #print(d)
    print(head(d))
    count = count + 1
    print(count)
    trait <- trait_list[count]
    print(trait)
    p[[count]] <- plot_g_h(d,fig_path,dat,trait)
  }
  
  plot_final <- ggarrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]],
                          p[[6]], p[[7]], p[[8]], p[[9]], p[[10]], 
                          ncol = 2, nrow = 5,  align = "v",
                          common.legend = TRUE)
  
  
  hypo_save(filename = paste0(figure_path,"univariate_comparison.png"), type="cairo",
            plot = plot_final,
            width = 20,
            height = 23)
  
}


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Determine the list of traits to analyse
traits <- list("PC1", "PC2", "PC3", "PC4", "PC5")
print(traits)

# Analyses of 10 univariate PCs and their heatmaps
univariate_plots(data_path, traits, figure_path, dataset)
