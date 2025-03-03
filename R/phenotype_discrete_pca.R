# by: Floriane Coulmance: 06/05/2021
# usage:
# Rscript phenotype_discrete_pca.R <path_phenotypes> <metadata_path> <pca_file> <variance_file> <im_file_path> <label> <figure_path>
# -------------------------------------------------------------------------------------------------------------------
# path_phenotypes in : $BASE_DIR/images/discrete/$COLOR_SPACE/$DATASET/$folder/
# metadata_path : image metadata file path $BASE_DIR/metadata/
# pca_file in : table of PCs in $BASE_DIR/images/discrete/$COLOR_SPACE/$DATASET/$folder/
# variance_file in : PCA variance table in $BASE_DIR/images/discrete/$COLOR_SPACE/$DATASET/$folder/
# im_file_path : $BASE_DIR/metadata/image_metadata.tsv
# label : label for file creation that correspond to the $DATASET and $folder considered
# figure_path : $BASE_DIR/figures/7_gxp/discrete/$COLOR_SPACE/$DATASET/$folder/
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
args = args[6:12]
print(args)

path_phenotypes <- as.character(args[1]) # Path to phenotype PCA files folder
print(path_phenotypes)
metadata_path <- as.character(args[2]) # Path to metadata folder to output phenotype metadata
print(metadata_path)
pca_file <- as.character(args[3]) # Name of the PCs table
print(pca_file)
variance_file <- as.character(args[4]) # Name of the variance table
print(variance_file)
im_file_path <- as.character(args[5]) # Path to all image metadata table
print(im_file_path)
label <- as.character(args[6]) # files label to give
print(label)
figure_path <- as.character(args[7]) # Path to figure folder
print(figure_path)


# -------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


plot_pca <- function(data, center_points, variance, fig_path, lbel) {
  
  # Function to plot PCA from dataframe and centroids of groups 
  
  p <- ggplot(data,aes(x=PC1.x,y=PC2.x,color=spec)) +
       geom_point(size = 1) +
       scale_color_manual(values=c("nig" = '#FF0033', "chl" = '#9900CC', "abe" = '#996600', "gut" = '#0000FF', "gum" = '#FF00FF', "ran" = '#666699', "gem" = '#CC0000', "may" = '#FF9933', "ind" = '#66CCFF', "pue" = '#FFCC00', "flo" = '#33FFCC', "tan" = '#333333', "uni" = '#66CC00'),    
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
                          breaks = c("nig", "chl", "abe", "gut", "gum", "ran", "gem", "may", "ind", "pue", "flo", "tan", "uni")) +
       geom_point(data=center_points,size=0.1) +
       geom_segment(aes(x=PC1.y, y=PC2.y, xend=PC1.x, yend=PC2.x, colour=spec), size = 0.1) +
       theme(legend.position="bottom",legend.title=element_blank(),legend.box = "vertical",legend.text =  element_markdown(size = 6)) +
       guides(color = guide_legend(nrow = 1)) +
       labs(x = paste0("PC1, var =  ", format(round(variance$X0[1] * 100, 1), nsmall = 1), " %") , y = paste0("PC2, var = ", format(round(variance$X0[2] * 100, 1), nsmall = 1), " %")) +
       ggtitle(paste0("PCA ", lbel))
  
  # Save figure
  hypo_save(filename = paste0(fig_path, lbel, "_pca.pdf"),
            plot = p,
            width = 12,
            height = 8)

}


write_metadata_gxp <- function(PCs, path_meta, K) {
  
  #  Make usable metadata for gxp pipeline
  
  PCs["sample"] <- gsub("\\-.*","",PCs$images)
  
  # Make sure sample names and labels are matching with genotyping data, correct naming errors
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
  
  # Reorganise and arrange table
  PCs["geo"] <- str_sub(PCs$sample,-3,-1)
  PCs["spec"] <- str_sub(PCs$sample,-6,-4)
  PCs["im"] <- gsub('.{4}$', '', PCs$images)
  PCs$X <- NULL
  
  # Save phenotype metadata table to metadata folder
  write.table(PCs,paste0(path_meta, lbel, "_PCs.csv"),sep=";",row.names=FALSE, quote = FALSE)
  
  # Return the table in a R variable
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

# Create a column with image name without suffix
im["im"] <- gsub('.{4}$', '', im$image) #Create a column with image name without suffix

# Format table as a metadata table reusable in the GWAS pipeline
PC_table <- write_metadata_gxp(pca_results, metadata_path, label) # Create a table with sample name for further gxp analysis

# Combine PCs info with image info and calculate centroids for each species group 
meta_table <- merge(PC_table, im, by = 'im') # Merge image metadata file to PCA results table 
centroids <- aggregate(cbind(PC1,PC2)~spec,PC_table,mean) # Create centroid table for PC1 PC2 for each of the species group
meta_table_centroid <- merge(meta_table, centroids, by = 'spec') # Merge centroid table with the image and PCA data table
centroids["PC1.x"] <- centroids["PC1"] # Create matching columns to meta_table_centroid in centroids table
centroids["PC2.x"] <- centroids["PC2"]

# Plot the phenotype PCA and save it as figure
plot_pca(meta_table_centroid, centroids, var, figure_path, label) # <- plot PCA and save it