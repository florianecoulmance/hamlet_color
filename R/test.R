# Clear the work space
rm(list = ls())

# Load needed library
library(stringr)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggtext)
library(hypoimg)


pca_results <- read.csv('/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/metadata/LAB_bodym_left_54off_59on_PCs.csv', sep = ";")
barred <- pca_results[ which( pca_results$spec == "pue" | pca_results$spec == "ind" | pca_results$spec == "flo") , ]
no_bar <- pca_results %>% filter(!pca_results$spec == "pue", !pca_results$spec == "ind", !pca_results$spec == "flo")


write.table(barred,"/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/metadata/LAB_bodym_left_54off_59on_38barred_PCs.csv",sep=";",row.names=FALSE, quote = FALSE)
write.table(no_bar,"/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/metadata/LAB_bodym_left_54off_59on_75unbarred_PCs.csv",sep=";",row.names=FALSE, quote = FALSE)
