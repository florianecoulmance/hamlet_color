#!/usr/bin/env Rscript
# by: Floriane Coulmance: 05/04/2022
# usage:
# Rscript gxp_heatmap_plots.R <data_path> <figure_path> <dataset> <metadata> <im_path>
# -------------------------------------------------------------------------------------------------------------------
# data_path in : $BASE_DIR/outputs/7_gxp/$DATASET
# figure_path in : $BASE_DIR/figures/7_gxp/figures/$TYPE/$COLOR_SPACE/$DATASET
# dataset in : LAB_fullm_54off_59on
# metadata in : $BASE_DIR/metadata/
# im_path in : $BASE_DIR/images/$TYPE/$COLOR_SPACE/$DATASET/
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
library(ggimage)


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


plot_pca <- function(pcf, pcs, data, center_points, variance, file_pc1, file_pc2, im1, im2, fig_path, dat) {
  
  # Function to plot PCA from dataframe and centroids of groups 
  print(as.numeric(stri_sub(pcf, -1)))
  print(typeof(as.numeric(stri_sub(pcf, -1))))
  
  # print(im1)
  # print(im2)
  # print(center_points)
  
  p <- ggplot(data,aes(x=.data[[paste0(pcf,".x")]],y=.data[[paste0(pcs,".x")]])) + #, label = sample
    geom_point(aes(colour = spec), size = 3) +
    # geom_text(hjust=0, vjust=-1, size=1) +
    scale_color_manual(values=c("nig" = '#FF0033', "chl" = '#9900CC', "abe" = '#996600', "gut" = '#0000FF',
                                "gum" = '#FF00FF', "ran" = '#666699', "gem" = '#CC0000', "may" = '#FF9933',
                                "ind" = '#66CCFF', "pue" = '#FFCC00', "flo" = '#33FFCC', "tan" = '#333333',
                                "uni" = '#66CC00'),
                       labels = c("nig" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_nigricans.l.cairo.png' width='90' /><br>*H. nigricans*",
                                  "chl" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_chlorurus.l.cairo.png' width='90' /><br>*H. chlorurus*",
                                  "abe" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_aberrans.l.cairo.png' width='90' /><br>*H. aberrans*",
                                  "gut" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_guttavarius.l.cairo.png' width='90' /><br>*H. guttavarius*",
                                  "gum" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_gumigutta.l.cairo.png' width='90' /><br>*H. gummigutta*",
                                  "ran" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_randallorum.l.cairo.png' width='90' /><br>*H. randallorum*",
                                  "gem" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_gemma.l.cairo.png' width='90' /><br>*H. gemma*",
                                  "may" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_maya.l.cairo.png' width='90' /><br>*H. maya*",
                                  "ind" = "<img src='//user/doau0129/work/chapter2/ressources/logos/H_indigo.l.cairo.png' width='90' /><br>*H. indigo*",
                                  "pue" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_puella.l.cairo.png' width='90' /><br>*H. puella*",
                                  "flo" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_floridae.l.cairo.png' width='90' /><br>*H. floridae*",
                                  "tan" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_sp.l.cairo.png' width='90' /><br>*Hypoplectrus* sp.",
                                  "uni" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_unicolor.l.cairo.png' width='90' /><br>*H. unicolor*"),
                       breaks = c("nig", "chl", "abe", "gut", "gum", "ran", "gem", "may", "ind", "pue",
                                  "flo", "tan", "uni")) +
    scale_shape_manual(values = c(16,3), labels = c(Off = "flash OFF", On = "flash ON")) +
    geom_segment(aes(x=.data[[paste0(pcf,".y")]], y=.data[[paste0(pcs,".y")]], xend=.data[[paste0(pcf,".x")]], yend=.data[[paste0(pcs,".x")]], colour=spec), size = 0.1) +
    geom_image(data=center_points, aes(image = link), size = 0.07, asp = 1) +
    geom_point(data=center_points,aes(colour = spec), size=20, alpha = 0.1) +
    theme(legend.position="none",
          legend.title=element_blank(),
          legend.box = "vertical",
          legend.text =  element_markdown(size = 20),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          text = element_text(size=20),
          legend.key=element_blank(),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          plot.margin = unit(c(0, 0, 0, 0.1), "cm")) +
    # guides(color = guide_legend(nrow = 13)) +
    scale_x_continuous(position = "top",labels = unit_format(unit = "k", scale = 1e-3)) +
    scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3)) +
    labs(x = paste0(pcf,", var =  ", format(round(variance$X0[as.numeric(stri_sub(pcf, -1))] * 100, 1), nsmall = 1), " %") ,
         y = paste0(pcs,", var = ", format(round(variance$X0[as.numeric(stri_sub(pcs, -1))] * 100, 1), nsmall = 1), " %")) +
    guides(color = guide_legend(nrow = 2))
  #ggtitle(paste0("PCA ", dat))
  
  # Additional plots of GWAS corresponding to each PCA axis
  pc1 <- ggplot() + geom_hypo_LG() +
    geom_point(data = file_pc1, aes(x = GPOS, y = AVG_P), size = .1) +
    scale_fill_hypo_LG_bg() +
    # scale_x_reverse(name = "Linkage Group", expand = c(0, 0), breaks = (hypo_karyotype$GSTART + hypo_karyotype$GEND)/2,
    #                 labels = 1:24, position = "top") +
    scale_x_hypo_LG(name = "Linkage Group") +
    scale_y_continuous(name = expression(italic('-log(p-value)')),
                       limits = c(0,3),
                       position="left") +
    theme_hypo() +
    theme(legend.position = 'none',
          legend.title=element_blank(),
          legend.box = "vertical",
          legend.text =  element_markdown(size = 30),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(size = 12),
          axis.text.x= element_text(size = 12, colour = 'darkgray'),
          plot.margin = unit(c(0.1, 0, 0.1, 0.4), "cm")) #+
    # rotate()
  
  pc2 <- ggplot() + geom_hypo_LG() +
    geom_point(data = file_pc2, aes(x = GPOS, y = AVG_P), size = .1) +
    scale_fill_hypo_LG_bg() +
    scale_x_reverse(name = "Linkage Group", expand = c(0, 0), breaks = (hypo_karyotype$GSTART + hypo_karyotype$GEND)/2,
                    labels = 1:24, position = "top") +
    # scale_x_hypo_LG(name = "Linkage Group") +
    scale_y_continuous(name = expression(italic('-log(p-value)')),
                       limits = c(0,3),
                       position = "right") +
    theme_hypo() +
    theme(legend.position = 'none',
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12, colour = 'darkgray'),
          plot.margin = unit(c(0, 0.1, 0, 0.1), "cm")) +
    # rotate()
    coord_flip()
  
  # Make the images heatmaps as plots
  g1 <- rasterGrob(im1[1:558,,], interpolate = T)
  g2 <- rasterGrob(im2[1:558,,], interpolate = T)
  g3 <- rasterGrob(im1[560:684,,], interpolate = T)
  
  # if ((pcf=="PC1")&(pcs=="PC5")) {
    
    # Arranging the plot
    g <- ggarrange(p, pc2, pc1, ggarrange(g1, g2, g3, nrow = 3, labels = c("", "e", ""), font.label = list(size = 30), widths = c(4,4,1), heights = c(4,4,1)), ncol = 2, nrow = 2, align = "h",
                   widths = c(3, 1), heights = c(3, 1), labels = c("a", "b", "c", "d"), font.label = list(size = 30), common.legend = TRUE, legend="bottom")
  # } else {
  #   
  #   # Arranging the plot
  #   g <- ggarrange(p, pc2, pc1, ggarrange(g2, g1, g3, nrow = 3, widths = c(4,4,1), heights = c(4,4,1)), ncol = 2, nrow = 2, align = "h",
  #                  widths = c(2, 1), heights = c(2, 1))
  #   
  # }
 
  return(g)
  
}


pca_analysis <- function(pc_first, pc_second, pca_pheno, var, im) {
  
  # Function that opens all necessary files and create PCA + univariate plots
  
  
  # Combine PCs info with image info and calculate centroids for each species group 
  meta_table <- merge(pca_pheno, im, by = 'im') # Merge image metadata file to PCA results table 
  print(meta_table)
  centroids <- aggregate(cbind(pca_pheno[[pc_first]],pca_pheno[[pc_second]])~spec,pca_pheno,mean) # Create centroid table for PC1 PC2 for each of the species group
  # print(centroids)
  centroids <- centroids %>% setNames(., nm = c("spec", pc_first, pc_second))
  # print(centroids)
  meta_table_centroid <- merge(meta_table, centroids, by = 'spec') # Merge centroid table with the image and PCA data table
  # print(meta_table_centroid)
  # centroids["PC1.x"] <- centroids["PC1"] # Create matching columns to meta_table_centroid in centroids table to be used in plots
  # centroids["PC2.x"] <- centroids["PC2"]
  centroids <- centroids %>% mutate(x = centroids[[pc_first]], y = centroids[[pc_second]]) %>% setNames(., nm = c("spec", pc_first, pc_second, paste0(pc_first,".x"), paste0(pc_second,".x")))
  # print(centroids)
  
  df <- data.frame(spec = c("nig", "chl", "abe", "gut", "gum", "ran", "gem", "may", "ind", "pue", "flo", "tan", "uni"),
                   link = c("/user/doau0129/work/chapter2/ressources/logos/H_nigricans.l.cairo.png",
                            "/user/doau0129/work/chapter2/ressources/logos/H_chlorurus.l.cairo.png",
                            "/user/doau0129/work/chapter2/ressources/logos/H_aberrans.l.cairo.png",
                            "/user/doau0129/work/chapter2/ressources/logos/H_guttavarius.l.cairo.png",
                            "/user/doau0129/work/chapter2/ressources/logos/H_gumigutta.l.cairo.png",
                            "/user/doau0129/work/chapter2/ressources/logos/H_randallorum.l.cairo.png",
                            "/user/doau0129/work/chapter2/ressources/logos/H_gemma.l.cairo.png",
                            "/user/doau0129/work/chapter2/ressources/logos/H_maya.l.cairo.png",
                            "/user/doau0129/work/chapter2/ressources/logos/H_indigo.l.cairo.png",
                            "/user/doau0129/work/chapter2/ressources/logos/H_puella.l.cairo.png",
                            "/user/doau0129/work/chapter2/ressources/logos/H_floridae.l.cairo.png",
                            "/user/doau0129/work/chapter2/ressources/logos/H_sp.l.cairo.png",
                            "/user/doau0129/work/chapter2/ressources/logos/H_unicolor.l.cairo.png"))
  
  centroids <- merge(centroids, df, by = 'spec') # Merge centroid table with the image and PCA data table
  # print(centroids)
  
  meta_table_centroid <- merge(meta_table_centroid, df, by = 'spec') # Merge centroid table with the image and PCA data table
  print(meta_table_centroid)
  
  # Get the PC1 GWAS plot for the univariate GWAS done with MVPLINK
  PC1 <- list(paste0(pc_first,".mvplink.50k.5k.txt.gz"))
  print(PC1)
  f1 <- read.table(paste0(data_path,PC1), header=TRUE)
  f1 <- f1 %>% left_join(hypo_chrom_start) %>% mutate(GPOS = MID_POS + GSTART)
  f1$range <- do.call(paste, c(f1[c("CHROM", "BIN_START", "BIN_END")], sep="_"))
  
  # Get the PC2 GWAS plot for the univariate GWAS done with MVPLINK
  PC2 <- list(paste0(pc_second,".mvplink.50k.5k.txt.gz"))
  print(PC2)
  f2 <- read.table(paste0(data_path,PC2), header=TRUE)
  f2 <- f2 %>% left_join(hypo_chrom_start) %>% mutate(GPOS = MID_POS + GSTART)
  f2$range <- do.call(paste, c(f2[c("CHROM", "BIN_START", "BIN_END")], sep="_"))
  
  img1 <- readPNG(paste0("/user/doau0129/work/chapter1/figures/7_gxp/continuous/LAB/",dataset,"/",dataset,"_",pc_first,".png"))
  img2 <- readPNG(paste0("/user/doau0129/work/chapter1/figures/7_gxp/continuous/LAB/",dataset,"/",dataset,"_",pc_second,".png"))
  # img1 <- readPNG(paste0("/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/figures/7_gxp/continuous/LAB/LAB_fullm_54off_59on/LAB_fullm_54off_59on_",pc_first,".png"))
  # img2 <- readPNG(paste0("/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/figures/7_gxp/continuous/LAB/LAB_fullm_54off_59on/LAB_fullm_54off_59on_",pc_second,".png"))
  
  # Plot the phenotype PCA and save it as figure
  p <- plot_pca(pc_first, pc_second, meta_table_centroid, centroids, var, f1, f2, img1, img2, figure_path, dataset)
  
  # Save the plot
  hypo_save(filename = paste0(figure_path,pc_first,"_",pc_second,"_univariate_gwas.pdf"),
            # type = "cairo",
            plot = p,
            width = 19,
            height = 22)
  
  return(p)
  
}


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Analyses of image PCA and gwas of PC1 and PC2
# Open all needed files
pheno_PC <- read.csv(paste0(metadata, dataset, "_PCs.csv"), sep = ";") #%>% 
            # mutate(spec = ifelse(sample == "PL17_23nigpue", "tan", spec), 
            #        sample = gsub("PL17_23nigpue", "PL17_23tanpue", sample))
print(pheno_PC)

# pca_pheno <- read.csv("/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/metadata/LAB_fullm_left_54off_59on_PCs.csv", sep = ";")
pheno_var <- read.csv(paste0(im_path, dataset, "_var.csv"), sep = ",")
# var <- read.csv("/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/images/continuous/LAB/LAB_fullm_left_54off_59on/LAB_fullm_left_54off_59on_var.csv", sep = ",")
immeta <- read.table(file = paste0(metadata, "image_metadata.tsv"), sep = '\t', header = TRUE)
# im <- read.table(file = "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/metadata/image_metadata.tsv", sep = '\t', header = TRUE)

# Create a column with image name without suffix
immeta["im"] <- gsub('.{4}$', '', immeta$image)
print(immeta$im)

# p1 <- pca_analysis("PC1", "PC2", pheno_PC, pheno_var, immeta)
p2 <- pca_analysis("PC1", "PC3", pheno_PC, pheno_var, immeta)
# p3 <- pca_analysis("PC1", "PC4", pheno_PC, pheno_var, immeta)
p4 <- pca_analysis("PC1", "PC5", pheno_PC, pheno_var, immeta)
# p5 <- pca_analysis("PC2", "PC3", pheno_PC, pheno_var, immeta)
p6 <- pca_analysis("PC2", "PC4", pheno_PC, pheno_var, immeta)
# p7 <- pca_analysis("PC2", "PC5", pheno_PC, pheno_var, immeta)
# p8 <- pca_analysis("PC3", "PC4", pheno_PC, pheno_var, immeta)
# p9 <- pca_analysis("PC3", "PC5", pheno_PC, pheno_var, immeta)
# p10 <- pca_analysis("PC4", "PC5", pheno_PC, pheno_var, immeta)
# p11 <- pca_analysis("PC1", "PC6", pheno_PC, pheno_var, immeta)

# plot <- ((p1 | p2 | p3) / (p5 | p6 | p7) / (p8 | p9 | p10))

# plot <- ggarrange(p1, p3, p5, p7, p8, p9, p10, ncol = 3, nrow = 3, common.legend = TRUE, legend = "bottom")

# hypo_save(filename = paste0(figure_path,"pca_univariate_gwas.png"),
#           type = "cairo",
#           plot = plot,
#           width = 30,
#           height = 31)

### histogram 
names <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15")
pheno_var_5 <- pheno_var[0:15,]
pheno_var_5["PCs"] <- names
pheno_var_5["percent"] <- pheno_var_5$X0*100

# var <- merge(pheno_var_5, pheno_cum_5, by="X")

v <- ggplot(data=pheno_var_5, aes(x=X, y=percent)) +
  geom_bar(stat="identity", fill = "#CC79A7") +
  scale_x_continuous(breaks=0:14,
                     labels=c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15")) +
  labs(x = "", y = "Variance (%)") +
  geom_text(aes(label=paste0(round(percent, 1),"%")), vjust=-1, color="black", size=8) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=35))
  
hypo_save(filename = paste0(figure_path,"var_histogram.pdf"),
          plot = v,
          width = 20,
          height = 12)

# p <- ggplot(data,aes(x=.data[[paste0(pcf,".x")]],y=.data[[paste0(pcs,".x")]],color=spec)) +
#   geom_point(size = 3) +
#   scale_color_manual(values=c("nig" = '#FF0033', "chl" = '#9900CC', "abe" = '#996600', "gut" = '#0000FF',
#                               "gum" = '#FF00FF', "ran" = '#666699', "gem" = '#CC0000', "may" = '#FF9933',
#                               "ind" = '#66CCFF', "pue" = '#FFCC00', "flo" = '#33FFCC', "tan" = '#333333',
#                               "uni" = '#66CC00'),
#                      labels = c("nig" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_nigricans.l.cairo.png' width='60' /><br>*H. nigricans*",
#                                 "chl" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_chlorurus.l.cairo.png' width='60' /><br>*H. chlorurus*",
#                                 "abe" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_aberrans.l.cairo.png' width='60' /><br>*H. aberrans*",
#                                 "gut" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_guttavarius.l.cairo.png' width='60' /><br>*H. guttavarius*",
#                                 "gum" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_gumigutta.l.cairo.png' width='60' /><br>*H. gummigutta*",
#                                 "ran" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_randallorum.l.cairo.png' width='60' /><br>*H. randallorum*",
#                                 "gem" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_gemma.l.cairo.png' width='60' /><br>*H. gemma*",
#                                 "may" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_maya.l.cairo.png' width='60' /><br>*H. maya*",
#                                 "ind" = "<img src='//user/doau0129/work/chapter2/ressources/logos/H_indigo.l.cairo.png' width='60' /><br>*H. indigo*",
#                                 "pue" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_puella.l.cairo.png' width='60' /><br>*H. puella*",
#                                 "flo" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_floridae.l.cairo.png' width='60' /><br>*H. floridae*",
#                                 "tan" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_tan.l.cairo.png' width='60' /><br>*Tan hamlet*",
#                                 "uni" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_unicolor.l.cairo.png' width='60' /><br>*H. unicolor*"),
#                      breaks = c("nig", "chl", "abe", "gut", "gum", "ran", "gem", "may", "ind", "pue",
#                                 "flo", "tan", "uni")) +
#   # scale_shape_manual(values = c(16,3), labels = c(Off = "flash OFF", On = "flash ON")) +
#   geom_point(data=center_points,size=7) +
#   # geom_image(aes(image = link), size = 0.15) +
#   geom_segment(aes(x=.data[[paste0(pcf,".y")]], y=.data[[paste0(pcs,".y")]], xend=.data[[paste0(pcf,".x")]], yend=.data[[paste0(pcs,".x")]], colour=spec), size = 0.1) +
#   theme(legend.position="none",legend.title=element_blank(),
#         legend.box = "vertical", legend.text =  element_markdown(size = 10),
#         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),
#         text = element_text(size=10), legend.key=element_blank(), axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) +
#   guides(color = guide_legend(nrow = 13)) +
#   scale_x_continuous(position = "top",labels = unit_format(unit = "k", scale = 1e-3)) +
#   scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3)) +
#   labs(x = paste0(pcf,", var =  ", format(round(variance$X0[as.numeric(stri_sub(pcf, -1))] * 100, 1), nsmall = 1), " %") ,
#        y = paste0(pcs,", var = ", format(round(variance$X0[as.numeric(stri_sub(pcs, -1))] * 100, 1), nsmall = 1), " %")) #+
# #ggtitle(paste0("PCA ", dat))
