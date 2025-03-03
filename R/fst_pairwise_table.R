#!/usr/bin/env Rscript
# by: Floriane Coulmance: 21/06/2022
# usage:
# Rscript fst_pairwise_table.R <fst_glob> <figure_path>
# -------------------------------------------------------------------------------------------------------------------
# fst_glob in : $BASE_DIR/outputs/8_fst/fst_globals_pop_nowindow.txt
# figure_path : $BASE_DIR/figures/8_fst/
# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------

# Clear the work space
rm(list = ls())

# Load needed library
library(GenomicOriginsScripts)
library(hypoimg)
library(hypogen)
library(vroom)
library(tidyverse)
library(ggplot2)
library(ggtext)
library(ggpubr)
library(gridExtra)


# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6:7]
print(args)

fst_glob <- as.character(args[1]) # Path to the global fst file
# fst_glob <- "/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/chapter1/fst_globals_pop_nowindow.txt"
figure_path <- as.character(args[2])
# figure_path <- "/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/chapter1/figures/8_fst/

# -------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Open global Fst
globals <- vroom(fst_glob, delim = '\t',
                 col_names = c('path','mean','weighted')) %>%
           mutate(pop1 = str_sub(path,-6,-1),
                  pop2 = str_sub(path,-13,-8),
                  weighted = round(weighted,3)) %>%
           select(pop1, pop2, weighted)

order <- c("puebel", "nigbel", "indbel", "maybel", "pueboc", "nigboc", "puepue", "unipue", "chlpue", "uniflo")

cmp_glob <- pivot_wider(globals, names_from = pop2, values_from = weighted) %>%
            slice(match(order, pop1)) %>%
            # column_to_rownames(var="pop1") %>%
            select(pop1, puebel, nigbel, indbel, maybel, pueboc, nigboc, puepue, unipue, chlpue, uniflo)


cmp_glob$indbel[which(cmp_glob$pop1 == "maybel")] <- cmp_glob$maybel[which(cmp_glob$pop1 == "indbel")]
cmp_glob$puepue[which(cmp_glob$pop1 == "unipue")] <- cmp_glob$unipue[which(cmp_glob$pop1 == "puepue")]
cmp_glob$puepue[which(cmp_glob$pop1 == "uniflo")] <- cmp_glob$uniflo[which(cmp_glob$pop1 == "puepue")]
cmp_glob$unipue[which(cmp_glob$pop1 == "uniflo")] <- cmp_glob$uniflo[which(cmp_glob$pop1 == "unipue")]
cmp_glob$chlpue[which(cmp_glob$pop1 == "uniflo")] <- cmp_glob$uniflo[which(cmp_glob$pop1 == "chlpue")]

cmp_glob$maybel[which(cmp_glob$pop1 == "indbel")] <- NA
cmp_glob$unipue[which(cmp_glob$pop1 == "puepue")] <- NA
cmp_glob$uniflo[which(cmp_glob$pop1 == "puepue")] <- NA
cmp_glob$uniflo[which(cmp_glob$pop1 == "unipue")] <- NA
cmp_glob$uniflo[which(cmp_glob$pop1 == "chlpue")] <- NA

cmp_glob <- sapply(cmp_glob, as.character)
cmp_glob[is.na(cmp_glob)] <- " "
cmp_glob[cmp_glob<0] <- " _ "
colnames(cmp_glob)[1] <- "POPULATIONS"


# fills <- rep("grey", each=nrow(cmp_glob[1:10,1]))

tt <- ttheme_default() #core=list(bg_params=list(fill=fills, col = "gray56"))
tt$core$fg_params <- list(fontface=matrix(c(2,1,1,1,1,1,1,1,1,1,1), ncol=ncol(cmp_glob),nrow=nrow(cmp_glob),byrow=TRUE))

pdf(paste0(figure_path,"fst_pairwise_table_nowindow.pdf"), height = 3.2, width = 8.5)
# p<-tableGrob(cmp_glob)
grid.table(cmp_glob, theme=tt)
# grid.table(cmp_glob[1:10,1], theme=tt)

dev.off()

