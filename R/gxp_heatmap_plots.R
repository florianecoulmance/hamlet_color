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
print(traits)

count = 0
l <- list()
print(count)
print(l)

for (trait in traits) {
    print(trait)
    count = count + 1
    print(count)

    f <- list.files(data_path, pattern = "lmm.50k.5k.txt.gz|assoc.50k.5k.txt.gz|mvplink.50k.5k.txt.gz")
    print(f)
    f <- f[grep(trait, names(f))]  
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

print(p)

plot_final <- grid.arrange(p[[1]], p[[2]], ..., p[[N]], nrow = 1, ncol =2)

hypo_save(filename = paste0(figure_path,"univariate_comparison.png"), type="cairo",
          plot = plot_final,
          width = 8,
          height = 8)






img <- readPNG(paste0(figure_path+dataset+"_pca.png"))
pca_plot <- rasterGrob(img, interpolate=TRUE)


PC1 <- list("PC1.mvplink.50k.5k.txt.gz")
f1 <- read.table(paste0(data_path,PC1), header=TRUE)
f1 <- bind_rows(PC1, .id = 'id') %>% left_join(hypo_chrom_start) %>% mutate(GPOS = MID_POS + GSTART)
f1$range <- do.call(paste, c(files[c("CHROM", "BIN_START", "BIN_END")], sep="_"))

pc1 <- ggplot() + facet_wrap(RUN~., ncol = 1, dir = 'v', strip.position="right") +
    geom_hypo_LG() +
    geom_point(data = f1, aes(x = GPOS, y = AVG_p_wald), size = .1) +
    scale_fill_hypo_LG_bg() +
    scale_x_hypo_LG(name = "Linkage Groups") +
    scale_y_continuous(name = expression(italic('-log(p-Wald)'))) +
    theme_hypo() +
    theme(legend.position = 'none',
          axis.title.x = element_text(),
          axis.text.x.top= element_text(colour = 'darkgray'))


PC2 <- list("PC2.mvplink.50k.5k.txt.gz")
f2 <- read.table(paste0(data_path,PC2), header=TRUE)
f2 <- bind_rows(PC2, .id = 'id') %>% left_join(hypo_chrom_start) %>% mutate(GPOS = MID_POS + GSTART)
f2$range <- do.call(paste, c(files[c("CHROM", "BIN_START", "BIN_END")], sep="_"))

pc2 <- ggplot() + facet_wrap(RUN~., ncol = 1, dir = 'v', strip.position="right") +
    geom_hypo_LG() +
    geom_point(data = f2, aes(x = GPOS, y = AVG_p_wald), size = .1) +
    scale_fill_hypo_LG_bg() +
    scale_x_hypo_LG(name = "Linkage Groups") +
    scale_y_continuous(name = expression(italic('-log(p-Wald)'))) +
    theme_hypo() +
    theme(legend.position = 'none',
          axis.title.x = element_text(),
          axis.text.x.top= element_text(colour = 'darkgray'))

# Arranging the plot
ggarrange(pca_plot, NULL, pc1, pc2, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)