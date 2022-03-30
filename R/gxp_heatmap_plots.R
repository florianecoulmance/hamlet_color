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


# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6:12]
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
pc_first <- as.character(args[6])
print(pc_first)
pc_second <- as.character(args[7])
print(pc_second)

# data_path <- "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/"
# figure_path <- "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/"
# dataset <- "LAB_fullm_PC1-2"
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


plot_pca <- function(data, center_points, variance, file_pc1, file_pc2, fig_path, dat) {
  
  # Function to plot PCA from dataframe and centroids of groups 
  print(as.numeric(stri_sub(pc_first, -1)))
  print(typeof(as.numeric(stri_sub(pc_first, -1))))

  p <- ggplot(data,aes(x=.data[[paste0(pc_first,".x")]],y=.data[[paste0(pc_second,".x")]],color=spec)) +
    geom_point(size = 3) +
    scale_color_manual(values=c("nig" = '#FF0033', "chl" = '#9900CC', "abe" = '#996600', "gut" = '#0000FF',
                                "gum" = '#FF00FF', "ran" = '#666699', "gem" = '#CC0000', "may" = '#FF9933',
                                "ind" = '#66CCFF', "pue" = '#FFCC00', "flo" = '#33FFCC', "tan" = '#333333',
                                "uni" = '#66CC00'),
                       labels = c("nig" = "<img src='/user/doau0129/work/chapter1/ressources/logos/H_nigricans.l.cairo.png' width='70' /><br>*H. nigricans*",
                                  "chl" = "<img src='/user/doau0129/work/chapter1/ressources/logos/H_chlorurus.l.cairo.png' width='70' /><br>*H. chlorurus*",
                                  "abe" = "<img src='/user/doau0129/work/chapter1/ressources/logos/H_aberrans.l.cairo.png' width='70' /><br>*H. aberrans*",
                                  "gut" = "<img src='/user/doau0129/work/chapter1/ressources/logos/H_guttavarius.l.cairo.png' width='70' /><br>*H. guttavarius*",
                                  "gum" = "<img src='/user/doau0129/work/chapter1/ressources/logos/H_gumigutta.l.cairo.png' width='70' /><br>*H. gummigutta*",
                                  "ran" = "<img src='/user/doau0129/work/chapter1/ressources/logos/H_randallorum.l.cairo.png' width='70' /><br>*H. randallorum*",
                                  "gem" = "<img src='/user/doau0129/work/chapter1/ressources/logos/H_gemma.l.cairo.png' width='70' /><br>*H. gemma*",
                                  "may" = "<img src='/user/doau0129/work/chapter1/ressources/logos/H_maya.l.cairo.png' width='70' /><br>*H. maya*",
                                  "ind" = "<img src='/user/doau0129/work/chapter1/ressources/logos/H_indigo.l.cairo.png' width='70' /><br>*H. indigo*",
                                  "pue" = "<img src='/user/doau0129/work/chapter1/ressources/logos/H_puella.l.cairo.png' width='70' /><br>*H. puella*",
                                  "flo" = "<img src='/user/doau0129/work/chapter1/ressources/logos/H_floridae.l.cairo.png' width='70' /><br>*H. floridae*",
                                  "tan" = "<img src='/user/doau0129/work/chapter1/ressources/logos/H_tan.l.cairo.png' width='70' /><br>*Tan hamlet*",
                                  "uni" = "<img src='/user/doau0129/work/chapter1/ressources/logos/H_unicolor.l.cairo.png' width='70' /><br>*H. unicolor*"),
                       breaks = c("nig", "chl", "abe", "gut", "gum", "ran", "gem", "may", "ind", "pue",
                                  "flo", "tan", "uni")) +
    # scale_shape_manual(values = c(16,3), labels = c(Off = "flash OFF", On = "flash ON")) +
    geom_point(data=center_points,size=7) +
    geom_segment(aes(x=.data[[paste0(pc_first,".y")]], y=.data[[paste0(pc_second,".y")]], xend=.data[[paste0(pc_first,".x")]], yend=.data[[paste0(pc_second,".x")]], colour=spec), size = 0.1) +
    theme(legend.position="left",legend.title=element_blank(),
          legend.box = "vertical", legend.text =  element_markdown(size = 15),
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),
          text = element_text(size=20), legend.key=element_blank(),
          legend.key.size = unit(0.7, 'cm'), plot.margin = unit(c(0,0,0.5,0.1), "cm")) +
    guides(color = guide_legend(nrow = 13)) +
    scale_x_continuous(position = "top") +
    labs(x = paste0(pc_first,", var =  ", format(round(variance$X0[as.numeric(stri_sub(pc_first, -1))] * 100, 1), nsmall = 1), " %") ,
         y = paste0(pc_second,", var = ", format(round(variance$X0[as.numeric(stri_sub(pc_second, -1))] * 100, 1), nsmall = 1), " %")) #+
  #ggtitle(paste0("PCA ", dat))
  
  
  # Additional plots of GWAS corresponding to each PCA axis
  pc1 <- ggplot() + geom_hypo_LG() +
    geom_point(data = file_pc1, aes(x = GPOS, y = AVG_P), size = .1) +
    scale_fill_hypo_LG_bg() +
    scale_x_hypo_LG(name = "Linkage Groups") +
    scale_y_continuous(name = expression(italic('-log(p-Wald)'))) +
    theme_hypo() +
    theme(legend.position = 'none',
          axis.title.x = element_text(),
          axis.text.x.top= element_text(colour = 'darkgray'),
          plot.margin = unit(c(0.5,0,0.5,0), "cm"))
  
  pc2 <- ggplot() + geom_hypo_LG() +
    geom_point(data = file_pc2, aes(x = GPOS, y = AVG_P), size = .1) +
    scale_fill_hypo_LG_bg() +
    # scale_x_hypo_LG(name = "Linkage Groups") +
    scale_y_continuous(name = expression(italic('-log(p-Wald)')), position = "top") +
    theme_hypo() +
    theme(legend.position = 'none',
          axis.title.x = element_text(),
          axis.text.x.top= element_text(colour = 'darkgray'),
          plot.margin = unit(c(0,0.5,0,0.1), "cm")) +
    rotate() + 
    scale_x_reverse(name = "Linkage Groups", expand = c(0, 0), breaks = (hypo_karyotype$GSTART + hypo_karyotype$GEND)/2,
                    labels = 1:24, position = "top")
  
  
  # Arranging the plot
  ggarrange(p, pc2, pc1, NULL, 
            ncol = 2, nrow = 2,  align = "hv",
            widths = c(3, 1), heights = c(3, 1),
            common.legend = TRUE, legend = "left")
  
}


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# # Determine the list of traits to analyse
# traits <- list("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
# print(traits)

# # Analyses of 10 univariate PCs and their heatmaps
# univariate_plots(data_path, traits, figure_path, dataset)


# Analyses of image PCA and gwas of PC1 and PC2
# Open all needed files
pca_pheno <- read.csv(paste0(metadata, dataset, "_PCs.csv"), sep = ";")
# pca_pheno <- read.csv("/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/metadata/LAB_fullm_left_54off_59on_PCs.csv", sep = ";")
var <- read.csv(paste0(im_path, dataset, "_var.csv"), sep = ",")
# var <- read.csv("/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/images/continuous/LAB/LAB_fullm_left_54off_59on/LAB_fullm_left_54off_59on_var.csv", sep = ",")
im <- read.table(file = paste0(metadata, "image_metadata.tsv"), sep = '\t', header = TRUE)
# im <- read.table(file = "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/metadata/image_metadata.tsv", sep = '\t', header = TRUE)

# Create a column with image name without suffix
im["im"] <- gsub('.{4}$', '', im$image)

# Combine PCs info with image info and calculate centroids for each species group 
meta_table <- merge(pca_pheno, im, by = 'im') # Merge image metadata file to PCA results table 
print(meta_table)
centroids <- aggregate(cbind(pca_pheno[[pc_first]],pca_pheno[[pc_second]])~spec,pca_pheno,mean) # Create centroid table for PC1 PC2 for each of the species group
print(centroids)
centroids <- centroids %>% setNames(., nm = c("spec", pc_first, pc_second))
print(centroids)
meta_table_centroid <- merge(meta_table, centroids, by = 'spec') # Merge centroid table with the image and PCA data table
print(meta_table_centroid)
# centroids["PC1.x"] <- centroids["PC1"] # Create matching columns to meta_table_centroid in centroids table to be used in plots
# centroids["PC2.x"] <- centroids["PC2"]
centroids <- centroids %>% mutate(x = centroids[[pc_first]], y = centroids[[pc_second]]) %>% setNames(., nm = c("spec", pc_first, pc_second, paste0(pc_first,".x"), paste0(pc_second,".x")))
print(centroids)

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

# Plot the phenotype PCA and save it as figure
p <- plot_pca(meta_table_centroid, centroids, var, f1, f2, figure_path, dataset)

# Save the plot
hypo_save(filename = paste0(figure_path,pc_first,"_",pc_second,"_univariate_gwas.pdf"),
          plot = p,
          width = 12.5,
          height = 12.5)
