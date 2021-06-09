# by: Floriane Coulmance: 06/05/2021
# usage:
# Rscript pcs_plots.R <path_phenotypes> <color_space> <metadata_path> <pca_file> <variance_file> <im_file_path>
# -------------------------------------------------------------------------------------------------------------------
# path_phenotypes in : $BASE_DIR/outputs/images/$DATASET
# metadata_path : image metadata file path $BASE_DIR/metadata/
# pca_file in : csv files with PCA results in $BASE_DIR/outputs/images/$DATASET/<color_space>/
# variance_file in : variance corresponding to PCA file in $BASE_DIR/outputs/images/$DATASET/<color_space>/
# im_file_path : $BASE_DIR/metadata/image_metadata.tsv
# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed library
library(stringr)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggtext)
library(hypoimg)


# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6:11]
print(args)

path_phenotypes <- as.character(args[1])
metadata_path <- as.character(args[2])
pca_file <- as.character(args[3])
variance_file <- as.character(args[4])
im_file_path <- as.character(args[5])
k <- as.character(args[6])
# path_phenotypes <- "/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/left_54off_59on/" # "$BASE_DIR/images/AB_fullm_left_54off_59on/"
# color_space <- "AB"
# metadata_path <- "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/metadata/"
# pca_file <- "AB_PCs_all_fullm.csv"
# variance_file <- "AB_var_all_fullm.csv"
# im_file_path <- "/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/metadata/image_metadata.tsv"


# -------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


plot_pca <- function(data, center_points, variance, pheno_path, K) {
  
  # Function to plot PCA from dataframe and centroids of groups 
  
  p <- ggplot(data,aes(x=PC1.x,y=PC2.x,color=spec))
  p <- p + geom_point(size = 1) 
  p <- p + scale_color_manual(values=c("nig" = '#FF0033', "chl" = '#9900CC', "abe" = '#996600', "gut" = '#0000FF', "gum" = '#FF00FF', "ran" = '#666699', "gem" = '#CC0000', "may" = '#FF9933', "ind" = '#66CCFF', "pue" = '#FFCC00', "flo" = '#33FFCC', "tan" = '#333333', "uni" = '#66CC00'),
                              labels = c("nig" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_nigricans.l.cairo.png' width='25' /><br>*H. nigricans*",
                                         "chl" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_chlorurus.l.cairo.png' width='25' /><br>*H. chlorurus*",
                                         "abe" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_aberrans.l.cairo.png' width='25' /><br>*H. aberrans*",
                                         "gut" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_guttavarius.l.cairo.png' width='25' /><br>*H. guttavarius*",
                                         "gum" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_gumigutta.l.cairo.png' width='25' /><br>*H. gummigutta*",
                                         "ran" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_randallorum.l.cairo.png' width='25' /><br>*H. randallorum*",
                                         "gem" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_gemma.l.cairo.png' width='25' /><br>*H. gemma*",
                                         "may" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_maya.l.cairo.png' width='25' /><br>*H. maya*",
                                         "ind" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_indigo.l.cairo.png' width='25' /><br>*H. indigo*",
                                         "pue" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_puella.l.cairo.png' width='25' /><br>*H. puella*",
                                         "flo" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_floridae.l.cairo.png' width='25' /><br>*H. floridae*",
                                         "tan" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_tan.l.cairo.png' width='25' /><br>*Tan hamlet*",
                                         "uni" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_unicolor.l.cairo.png' width='25' /><br>*H. unicolor*"),
                              breaks = c("nig", "chl", "abe", "gut", "gum", "ran", "gem", "may", "ind", "pue", "flo", "tan", "uni"))
  p <- p + geom_point(data=center_points,size=0.1)
  p <- p + geom_segment(aes(x=PC1.y, y=PC2.y, xend=PC1.x, yend=PC2.x, colour=spec), size = 0.1)
  p <- p + theme(legend.position="bottom",legend.title=element_blank(),legend.box = "vertical",legend.text =  element_markdown(size = 6))
  p <- p + guides(color = guide_legend(nrow = 1))
  p <- p + labs(x = paste0("PC1, var =  ", format(round(variance$X0[1] * 100, 1), nsmall = 1), " %") , y = paste0("PC2, var = ", format(round(variance$X0[2] * 100, 1), nsmall = 1), " %"))
  p <- p + ggtitle(paste0("PCA k", K, " color6, binary images"))
  
  
  hypo_save(filename = paste0(pheno_path,"k", K, "_color6_pca.png"),
            plot = p,
            width = 12,
            height = 8)
}


write_metadata_gxp <- function(PCs, path_meta, K) {
  
  #  Make usable metadata for gxp pipeline
  
  PCs["sample"] <- gsub("\\-.*","",PCs$images)
  
  PCs <- PCs %>% mutate(sample = gsub("por", "pue", sample),
                        sample = gsub("AG9RX", "AG9RX_", sample),
                        sample = gsub("PL17_0", "PL17_", sample),
                        sample = gsub("PL17_1uniboc", "PL17_01uniboc", sample),
                        sample = gsub("PL17_2pueboc", "PL17_02pueboc", sample),
                        sample = gsub("PL17_4pueboc", "PL17_04pueboc", sample),
                        sample = gsub("PL17_5pueboc", "PL17_05pueboc", sample),
                        sample = gsub("PL17_160pueflo", "PL17_160floflo", sample),
                        sample = gsub("PL17_125tanbel", "PL17_125ranbel", sample),
                        sample = gsub("PL17_62puepue", "PL17_62chlpue", sample))
  
  PCs["geo"] <- str_sub(PCs$sample,-3,-1)
  PCs["spec"] <- str_sub(PCs$sample,-6,-4)
  PCs["im"] <- gsub('.{4}$', '', PCs$images)
  PCs$X <- NULL
  
  write.table(PCs,paste0(path_meta,"k", K, "_color6_PCs.csv"),sep=";",row.names=FALSE, quote = FALSE)
  
  return(PCs)
  
}  


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------

# Open all needed files
pca_results <- read.csv(paste0(path_phenotypes,pca_file))
var <- read.csv(paste0(path_phenotypes,variance_file), sep = ",")
im <- read.table(file = im_file_path, sep = '\t', header = TRUE)

im["im"] <- gsub('.{4}$', '', im$image) # <- create a column with image name without suffix

PC_table <- write_metadata_gxp(pca_results, metadata_path, k) # <- create a table with sample name for further gxp analysis

meta_table <- merge(PC_table, im, by = 'im') # <- merge image metadata file to PCA results table 
centroids <- aggregate(cbind(PC1,PC2)~spec,PC_table,mean) # <- create centroid table for PC1 PC2 for each of the species group
meta_table_centroid <- merge(meta_table, centroids, by = 'spec') # <- merge centroid table with the image and PCA data table
centroids["PC1.x"] <- centroids["PC1"] # <- create matching columns to meta_table_centroid in centroids table
centroids["PC2.x"] <- centroids["PC2"]

plot_pca(meta_table_centroid, centroids, var, path_phenotypes, k) # <- plot PCA and save it