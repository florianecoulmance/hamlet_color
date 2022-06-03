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

threshold <- 2.5
thresh <- file_gwas[file_gwas[, "AVG_P"] >= threshold,]
print(thresh)

p <- ggplot() +
     geom_hypo_LG() +
     geom_vline(data = thresh, aes(xintercept = GPOS), color = "red", size = 0.3) +
     geom_point(data = file_gwas, aes(x = GPOS, y = AVG_P), size = .1) +
     scale_fill_hypo_LG_bg() +
     scale_x_hypo_LG(name = "Linkage Groups") +
     scale_y_continuous(name = expression(bolditalic(-log[10](p)))) +
     theme_hypo() +
     ggtitle("a") +
     theme(plot.title = element_text(hjust = 0, vjust = -4, size = 13, face = "bold"),
           legend.position = 'none',
           axis.title.x = element_text(),
           axis.text.x.top= element_text(colour = 'darkgray'),
           axis.title.y = ggplot2::element_text(angle = 90, size = 6))
  
g1 <- rasterGrob(zoom1, interpolate = T)
g2 <- rasterGrob(zoom2, interpolate = T)
g3 <- rasterGrob(zoom3, interpolate = T)

g_plot1 <- ggplot() +
           annotation_custom(g1) +
           ggtitle("b") +
           theme(plot.title = element_text(hjust = 0.15, vjust = -2, size = 13, face = "bold"),
                 plot.background = element_blank(),
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(-2, 0, 0, -0.65), "cm"))


g_plot2 <- ggplot() +
          annotation_custom(g2) +
          ggtitle("c") +
          theme(plot.title = element_text(hjust = 0.15, vjust = -2, size = 13, face = "bold"),
                plot.background = element_blank(),
                panel.background = element_blank(),
                panel.grid = element_blank(),
                panel.border = element_blank(),
                plot.margin = unit(c(-2, 0, 0, -0.3), "cm"))

g_plot3 <- ggplot() +
           annotation_custom(g3) +
           ggtitle("d") +
           theme(plot.title = element_text(hjust = 0.15, vjust = -2, size = 13, face = "bold"),
                 plot.background = element_blank(),
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(-2, 0, 0, 0), "cm"))

g <- ggarrange(p, ggarrange(g_plot1, g_plot2, g_plot3, ncol = 3, widths = c(1,1,1)), ncol = 1, nrow = 2, align = "hv",
                 widths = c(2, 2), heights = c(1, 3)) 
  

# try <- ggarrange(g_plot1, g_plot2, g_plot3, ncol = 3)


hypo_save(filename = paste0(figure_path,"PC1_5_gwas_zoom.pdf"),
            plot = g,
            width = 6,
            height = 5)

