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



vroom::vroom("/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/chlpue_puepue.fst.tsv.gz", delim = '\t') %>%
  # filter(!is.na(WEIR_AND_COCKERHAM_FST)) %>% # remove na values from the weir column
  mutate(FST_RANK = rank(-WEIR_AND_COCKERHAM_FST, ties.method = "random")) %>% # create column with rank of the Fst values 
  select(CHROM, POS, WEIR_AND_COCKERHAM_FST, FST_RANK) %>% 
  filter(FST_RANK < number)  %>% # filter 800 positions
  select(CHROM, POS) %>% # select column to keep
  write_tsv(path = "Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/800SNPs.snps") # write output to file





