#!/usr/bin/env Rscript
# by: Floriane Coulmance: 01/05/2020
# usage:
# Rscript --vanilla pca.R <vcf_file> <output_path> <out_prefix>
# -------------------------------------------------------------------------------------------------------------------
# vcf_file in : 1) $BASE_DIR/ 2) 3) 4) 
# output_path in : $BASE_DIR/pca/
# out_prefix in : 1) 2) 3) 4)
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed library
library(gdsfmt)
library(SNPRelate)
library(stringi)
library(tidyverse)


# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[7:9]
print(args)

vcf_file <- as.character(args[1]) # Path to vcf file
output_path <- as.character(args[2]) # Path to the figure folder
out_prefix < as.character(args[3])


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS


# Read the files and convert to pca format
ld_threshold <- 1
gds_file <- str_c(output_path,out_prefix,".gds")
snpgdsVCF2GDS(vcf.fn=vcf_file, out.fn=gds_file, method="biallelic.only")
genofile <- snpgdsOpen(gds_file)
snpset <- snpgdsLDpruning(genofile, ld.threshold = ld_threshold, autosome.only = FALSE)
snpset.id <- unlist(snpset)


# perfor pca
pca <- snpgdsPCA(genofile, snp.id = snpset.id, num.thread = 2, autosome.only = FALSE)

pca_tib <- pca$eigenvect %>%
  as_tibble() %>%
  set_names(nm = str_pad(1:length(names(.)),2,pad = '0') %>%
              str_c("PC",.)) %>%
  mutate(sample = pca$sample.id,
         spec = str_sub(sample,-6,-4),
         loc = str_sub(sample,-3,-1),
         pop = str_c(spec,loc))

explained_var <- tibble(variation = pca$varprop*100) %>%
  mutate(PC = str_c("PC", str_pad(row_number(),
                                  width = 2,
                                  pad = "0"))) %>%
  select(PC, variation) %>%
  filter(!is.na(variation))


# Save important files for plots
save(pca,file = str_c(output_path,out_prefix,"_pca.RData"))
pca_tib %>% write_tsv(path = str_c(output_path,out_prefix,"_eigenvectors.tsv"))
explained_var %>% write_tsv(path = str_c(output_path,out_prefix,"_eixplained_var.tsv"))


# Close files
snpgdsClose(genofile)
system(str_c("rm ", gds_file))
showfile.gds(closeall=TRUE)
