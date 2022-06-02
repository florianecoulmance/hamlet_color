#!/usr/bin/env Rscript
# by: Floriane Coulmance: 30/05/2022
# usage:
# Rscript .R <data_path> <figure_path>
# -------------------------------------------------------------------------------------------------------------------
# data_path in : $BASE_DIR/outputs/7_gxp/$DATASET
# figure_path in : $BASE_DIR/outputs/7_gxp/figures/$DATASET/
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
library(png)
library(ggpubr)


# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6:7]
print(args)

data_path <- as.character(args[1]) # Path to mean coverage data table
# data_path <- "/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/chapter1/"
figure_path <- as.character(args[2]) # Path to the figure folder
# figure_path <- "/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/chapter1/figures/7_gxp/continuous/LAB/LAB_fullm_54off_59on/"


# -------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------



# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


file_gwas <- read.table(paste0(data_path, "PC1_5.mvplink.50k.5k.txt.gz"), header = TRUE) %>%
             left_join(hypo_chrom_start) %>%
             mutate(GPOS = MID_POS + GSTART)
file_gwas$range <- do.call(paste, c(file_gwas[c("CHROM", "BIN_START", "BIN_END")], sep="_"))

zoom1 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG04_2.png"))
zoom2 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG12_2.png"))
zoom3 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG12_3.png"))

plotgwas_mvp <- function(dataset,path,analysis) {
  p <- ggplot() +
    geom_hypo_LG() +
    geom_point(data = file_gwas, aes(x = GPOS, y = AVG_P), size = .1) +
    scale_fill_hypo_LG_bg() +
    scale_x_hypo_LG(name = "Linkage Groups") +
    scale_y_continuous(name = expression(bolditalic(-log[10](p)))) +
    theme_hypo() +
    theme(legend.position = 'none',
          axis.title.x = element_text(),
          axis.text.x.top= element_text(colour = 'darkgray'))
  
  g1 <- rasterGrob(zoom1, interpolate = T)
  g2 <- rasterGrob(zoom2, interpolate = T)
  g3 <- rasterGrob(zoom3, interpolate = T)
  
  g <- ggarrange(p, ggarrange(g1, g2, g3, ncol = 3, widths = c(1,1,1), heights = c(7,7,7)), ncol = 1, nrow = 2, align = "v",
                 widths = c(2, 2), heights = c(1, 7)) 
  
  
  hypo_save(filename = paste0(figure_path,"PC1_5_gwas_zoom.pdf"),
            plot = g,
            width = 8,
            height = 5)
}
