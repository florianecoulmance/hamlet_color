#!/usr/bin/env Rscript
# by: Floriane Coulmance: 01/03/2022
# usage:
# Rscript gxp_heatmap_plots.R <data_path> <figure_path> <dataset>
# -------------------------------------------------------------------------------------------------------------------
# data_path in : $BASE_DIR/outputs/7_gxp/$DATASET
# figure_path in : $BASE_DIR/figures/7_gxp/figures/$TYPE/$COLOR_SPACE/$DATASET
# dataset in : LAB_fullm_54off_59on
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
library(stringr)
library(png)
library(grid)


# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6:7]
print(args)

data_path <- as.character(args[1]) # Path to mean coverage data table
print(data_path)
figure_path <- as.character(args[2]) # Path to the figure folder
print(figure_path)


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
    #assign(run_files,d)
    l[[count]]=assign(run_files,d)
    print(head(d))
    #print(l)
  }
  return(l)
}


plot_g_h <- function(table,path,prefix,trait) {
  p <- ggplot() + facet_wrap(RUN~., ncol = 1, dir = 'v', strip.position="right") +
    geom_hypo_LG() +
    geom_point(data = table, aes(x = GPOS, y = AVG_p_wald), size = .1) +
    scale_fill_hypo_LG_bg() +
    scale_x_hypo_LG(name = "Linkage Groups") +
    scale_y_continuous(name = expression(italic('-log(p-Wald)'))) +
    theme_hypo() +
    theme(legend.position = 'none',
          axis.title.x = element_text(),
          axis.text.x.top= element_text(colour = 'darkgray'))

    # im <- paste0(figure_path+dataset+"_"+trait+".png")
    img <- readPNG(paste0(path+prefix+"_"+trait+".png"))
    g <- rasterGrob(img, interpolate=TRUE)

    plot <- p+g 


}


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


traits <- list("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

count = 0
l <- list()
for (trait in traits) {
    print(trait)
    count = count + 1 

    f <- list.files(data_path, pattern = "$trait.|lmm.50k.5k.txt.gz|assoc.50k.5k.txt.gz|mvplink.50k.5k.txt.gz")
    print(f)
    model <- list("GEMMA", "PLINK", "MV PLINK")
    print(model)
    files_l <- concat_files_gem(f,data_path)
    names(files_l) <- model
    f <- bind_rows(files_l, .id = 'id') %>% left_join(hypo_chrom_start) %>% mutate(GPOS = MID_POS + GSTART)
    f$range <- do.call(paste, c(files[c("CHROM", "BIN_START", "BIN_END")], sep="_"))

    l[[count]]=assign(trait,f)

}

print(head(l))

p <- vector('list', len(l))
print(p)


for (d in l) {
  print(d)
  print(head(d))

  p[[d]] <- plot_g_h(d,figure_path,dataset,d)

}


plot_final <- grid.arrange(p[[1]], p[[2]], ..., p[[N]], nrow = 1, ncol =2)

hypo_save(filename = paste0(figure_path,"univariate_comparison.png"), type="cairo",
          plot = plot_final,
          width = 8,
          height = 8)


