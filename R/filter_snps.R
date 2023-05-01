#!/usr/bin/env Rscript
# by: Floriane Coulmance: 02/10/2020
# usage:
# Rscript filter_snps.R in.weir.fst.gz 800 out_prefix
# -------------------------------------------------------------------------------------------------------------------
# in.weir.fst.gz in : 
# 800 : number of SNPs to select
# out_prefix : the population pairs prefix
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed library
# library(GenomicOriginsScripts)
library(tidyverse)
library(vroom)


# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args <- commandArgs(trailingOnly = FALSE)
args = args[6:8]
print(args)
in_file <- as.character(args[1]) # Path to input file
nr_snps <- as.numeric(args[2]) # Number of SNPs
out_prefix <- as.character(args[3]) # Path+prefix for new filename
# in_file <- "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1_2/chlpue_puepue.fst.tsv.gz"
# nr_snps <- 800
# out_prefix <- "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1_2/chlpue_puepue"


# -------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


snp_random <- function(file, number, prefix) {

  # Takes input vcf Fst file and select a number of random SNP positions

  vroom::vroom(file, delim = '\t') %>%
  filter(!is.na(WEIR_AND_COCKERHAM_FST)) %>% # remove na values from the weir column
  mutate(FST_RANK = rank(-WEIR_AND_COCKERHAM_FST, ties.method = "random")) %>% # create column with rank of the Fst values 
  select(CHROM, POS, WEIR_AND_COCKERHAM_FST, FST_RANK) %>% 
  filter(FST_RANK < number)  %>% # filter 800 positions
  select(CHROM, POS) %>% # select column to keep
  write_tsv(path = str_c(prefix, "_", number, "SNPs",'.snps' )) # write output to file

}


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


snp_random(in_file, nr_snps, out_prefix) # execute the function to create new file with 800 randomly selected SNPs

