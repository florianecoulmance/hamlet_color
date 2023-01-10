#!/usr/bin/env Rscript
# by: Floriane Coulmance: 30/05/2022
# usage:
# Rscript gxp_pc1_5_3peaks.R <data_path> <figure_path>
# -------------------------------------------------------------------------------------------------------------------
# data_path in : $BASE_DIR/outputs/7_gxp/$DATASET
# figure_path in : $BASE_DIR/outputs/7_gxp/figures/$DATASET/
# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------

# Clear the work space
rm(list = ls())

# Load needed library
library(ggplot2)
library(patchwork)
# library(tidyverse)
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
zoom2 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG12_3.png"))
zoom3 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG12_4.png"))

chart1 <- readPNG(paste0(figure_path,"PC1_5/LG04_11599216_pie.png"))
chart2 <- readPNG(paste0(figure_path,"PC1_5/LG12_20232954_pie.png"))
chart3 <- readPNG(paste0(figure_path,"PC1_5/LG12_22227968_pie.png"))

threshold <- 2.5
thresh <- file_gwas[file_gwas[, "AVG_P"] >= threshold,]
print(thresh)

p <- ggplot() +
     geom_hypo_LG() +
     geom_vline(data = thresh, aes(xintercept = GPOS), color = "red", alpha = 0.5, size = 0.4) +
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

pie1 <- rasterGrob(chart1, interpolate = T)
pie2 <- rasterGrob(chart2, interpolate = T)
pie3 <- rasterGrob(chart3, interpolate = T)

g_plot1 <- ggplot() +
           annotation_custom(g1) +
           ggtitle("b") +
           theme(plot.title = element_text(hjust = 0.15, vjust = -2, size = 13, face = "bold"),
                 plot.background = element_blank(),
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(-2, 0, 0, -0.95), "cm"))


g_plot2 <- ggplot() +
          annotation_custom(g2) +
          ggtitle("c") +
          theme(plot.title = element_text(hjust = 0.12, vjust = -2, size = 13, face = "bold"),
                plot.background = element_blank(),
                panel.background = element_blank(),
                panel.grid = element_blank(),
                panel.border = element_blank(),
                plot.margin = unit(c(-2, 0, 0, -0.4), "cm"))

g_plot3 <- ggplot() +
           annotation_custom(g3) +
           ggtitle("d") +
           theme(plot.title = element_text(hjust = 0.11, vjust = -2, size = 13, face = "bold"),
                 plot.background = element_blank(),
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(-2, 0, 0, 0), "cm"))

pie_1 <- ggplot() + 
         annotation_custom(pie1) +
         # ggtitle("e") +
         theme(plot.title = element_text(hjust = 0.11, vjust = -2, size = 13, face = "bold"),
               plot.background = element_blank(),
               panel.background = element_blank(),
               panel.grid = element_blank(),
               panel.border = element_blank(),
               plot.margin = unit(c(-2, 0, 0, 0), "cm"))

pie_2 <- ggplot() + 
         annotation_custom(pie2) +
         # ggtitle("f") +
         theme(plot.title = element_text(hjust = 0.11, vjust = -2, size = 13, face = "bold"),
               plot.background = element_blank(),
               panel.background = element_blank(),
               panel.grid = element_blank(),
               panel.border = element_blank(),
               plot.margin = unit(c(-2, 0, 0, 0), "cm"))

pie_3 <- ggplot() + 
         annotation_custom(pie3) +
         # ggtitle("g") +
         theme(plot.title = element_text(hjust = 0.11, vjust = -2, size = 13, face = "bold"),
               plot.background = element_blank(),
               panel.background = element_blank(),
               panel.grid = element_blank(),
               panel.border = element_blank(),
               plot.margin = unit(c(-2, 0, 0, 0), "cm"))

g <- ggarrange(p, ggarrange(g_plot1, g_plot2, g_plot3, ncol = 3, widths = c(1,1,1)), ggarrange(pie_1, pie_2, pie_3, ncol = 3, widths = c(1,1,1)), ncol = 1, nrow = 3, align = "hv",
                 widths = c(2, 2), heights = c(1, 3)) 
  

# try <- ggarrange(g_plot1, g_plot2, g_plot3, ncol = 3)


hypo_save(filename = paste0(figure_path,"PC1_5_gwas_zoom.pdf"),
            plot = g,
            width = 8.3,
            height = 10)
          
            # width = 6,
            # height = 5)



# zoom4 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG03_1.png"))
# zoom5 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG04_1.png"))
# zoom6 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG05_1.png"))
# zoom7 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG06_1.png"))
# zoom8 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG09_1.png"))
# zoom9 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG09_2.png"))
# zoom10 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG12_1.png"))
# zoom11 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG19_1.png"))
# zoom12 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG20_1.png"))
# zoom13 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG23_1.png"))


# g4 <- rasterGrob(zoom4, interpolate = T)
# g5 <- rasterGrob(zoom5, interpolate = T)
# g6 <- rasterGrob(zoom6, interpolate = T)
# g7 <- rasterGrob(zoom7, interpolate = T)
# g8 <- rasterGrob(zoom8, interpolate = T)
# g9 <- rasterGrob(zoom9, interpolate = T)
# g10 <- rasterGrob(zoom10, interpolate = T)
# g11 <- rasterGrob(zoom11, interpolate = T)
# g12 <- rasterGrob(zoom12, interpolate = T)
# g13 <- rasterGrob(zoom13, interpolate = T)

# g_plot4 <- ggplot() +
#   annotation_custom(g4) +
#   theme(plot.title = element_text(hjust = 0.15, vjust = -2, size = 13, face = "bold"),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         plot.margin = unit(c(-2, 0, 0, -0.95), "cm"))


# g_plot5 <- ggplot() +
#   annotation_custom(g5) +
#   theme(plot.title = element_text(hjust = 0.12, vjust = -2, size = 13, face = "bold"),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         plot.margin = unit(c(-2, 0, 0, -0.4), "cm"))

# g_plot6 <- ggplot() +
#   annotation_custom(g6) +
#   theme(plot.title = element_text(hjust = 0.11, vjust = -2, size = 13, face = "bold"),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         plot.margin = unit(c(-2, 0, 0, 0), "cm"))

# g_plot7 <- ggplot() +
#   annotation_custom(g7) +
#   theme(plot.title = element_text(hjust = 0.11, vjust = -2, size = 13, face = "bold"),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         plot.margin = unit(c(-2, 0, 0, 0), "cm"))

# g_plot8 <- ggplot() +
#   annotation_custom(g8) +
#   theme(plot.title = element_text(hjust = 0.11, vjust = -2, size = 13, face = "bold"),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         plot.margin = unit(c(-2, 0, 0, 0), "cm"))

# g_plot9 <- ggplot() +
#   annotation_custom(g9) +
#   theme(plot.title = element_text(hjust = 0.11, vjust = -2, size = 13, face = "bold"),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         plot.margin = unit(c(-2, 0, 0, 0), "cm"))

# g_plot10 <- ggplot() +
#   annotation_custom(g10) +
#   theme(plot.title = element_text(hjust = 0.11, vjust = -2, size = 13, face = "bold"),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         plot.margin = unit(c(-2, 0, 0, 0), "cm"))

# g_plot11 <- ggplot() +
#   annotation_custom(g11) +
#   theme(plot.title = element_text(hjust = 0.11, vjust = -2, size = 13, face = "bold"),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         plot.margin = unit(c(-2, 0, 0, 0), "cm"))

# g_plot12 <- ggplot() +
#   annotation_custom(g12) +
#   theme(plot.title = element_text(hjust = 0.11, vjust = -2, size = 13, face = "bold"),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         plot.margin = unit(c(-2, 0, 0, 0), "cm"))

# g_plot13 <- ggplot() +
#   annotation_custom(g13) +
#   theme(plot.title = element_text(hjust = 0.11, vjust = -2, size = 13, face = "bold"),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         plot.margin = unit(c(-2, 0, 0, 0), "cm"))

# g <- ggarrange(g_plot4, g_plot5, g_plot6,
#                g_plot7, g_plot8, g_plot9,
#                g_plot10, g_plot11, g_plot12,
#                g_plot13, ncol = 3, widths = c(1,1), nrow = 4, align = "hv",
#                labels = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j")) 

# hypo_save(filename = paste0(figure_path,"PC1_5_gwas_zoom_others.pdf"),
#           plot = g,
#           width = 10,
#           height = 21)



# zoom14 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG08_1.png"))
# zoom15 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG08_2.png"))
# zoom16 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG08_3.png"))
# zoom17 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG08_4.png"))
# zoom18 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG08_5.png"))
# zoom19 <- readPNG(paste0(figure_path,"PC1_5/PC1_5_LG08_6.png"))


# g14 <- rasterGrob(zoom14, interpolate = T)
# g15 <- rasterGrob(zoom15, interpolate = T)
# g16 <- rasterGrob(zoom16, interpolate = T)
# g17 <- rasterGrob(zoom17, interpolate = T)
# g18 <- rasterGrob(zoom18, interpolate = T)
# g19 <- rasterGrob(zoom19, interpolate = T)


# g_plot14 <- ggplot() +
#   annotation_custom(g14) +
#   theme(plot.title = element_text(hjust = 0.11, vjust = -2, size = 13, face = "bold"),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         plot.margin = unit(c(-2, 0, 0, 0), "cm"))

# g_plot15 <- ggplot() +
#   annotation_custom(g15) +
#   theme(plot.title = element_text(hjust = 0.11, vjust = -2, size = 13, face = "bold"),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         plot.margin = unit(c(-2, 0, 0, 0), "cm"))

# g_plot16 <- ggplot() +
#   annotation_custom(g16) +
#   theme(plot.title = element_text(hjust = 0.11, vjust = -2, size = 13, face = "bold"),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         plot.margin = unit(c(-2, 0, 0, 0), "cm"))

# g_plot17 <- ggplot() +
#   annotation_custom(g17) +
#   theme(plot.title = element_text(hjust = 0.11, vjust = -2, size = 13, face = "bold"),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         plot.margin = unit(c(-2, 0, 0, 0), "cm"))

# g_plot18 <- ggplot() +
#   annotation_custom(g18) +
#   theme(plot.title = element_text(hjust = 0.11, vjust = -2, size = 13, face = "bold"),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         plot.margin = unit(c(-2, 0, 0, 0), "cm"))

# g_plot19 <- ggplot() +
#   annotation_custom(g19) +
#   theme(plot.title = element_text(hjust = 0.11, vjust = -2, size = 13, face = "bold"),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         plot.margin = unit(c(-2, 0, 0, 0), "cm"))

# g <- ggarrange(g_plot14, g_plot15, g_plot16,
#                g_plot17, g_plot18, g_plot19, 
#                ncol = 3, widths = c(1,1), nrow = 2, align = "hv",
#                labels = c("a", "b", "c", "d", "e", "f")) 

# hypo_save(filename = paste0(figure_path,"PC1_5_gwas_zoom_LG08.pdf"),
#           plot = g,
#           width = 10,
#           height = 11)


