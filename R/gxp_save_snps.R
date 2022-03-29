# by: Floriane Coulmance: 15/04/2021
# usage:
# Rscript gxp_zooms.R <data_path> <data_file> <figure_path> <pc>
# -------------------------------------------------------------------------------------------------------------------
# data_path in : $BASE_DIR/outputs/7_gxp/$DATASET
# data_file in : *mvplink.50k.5k.txt.gz
# figure_path in : $BASE_DIR/outputs/figures/7_gxp/$DATASET/
# pc in : 55 names of univariate and multivariate PCs association analysis
# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed library
library(broom, lib.loc=.libPaths()[-1])
library(dbplyr, lib.loc=.libPaths()[-1])
# library(broom)
# library(dbplyr)
library(GenomicOriginsScripts)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(hypogen)
library(hypoimg)
library(dplyr)
library(furrr)
library(ggraph)
library(tidygraph)
library(ggtext)
library(ape)
library(igraph)


# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6:9]
print(args)

data_path <- as.character(args[1]) # Path to data folder with GWAS files
data_file <- as.character(args[2]) # Name of the GWAS 50kb windowed-averaged file
figure_path <- as.character(args[3]) # Path to the figure folder
pc <- as.character(args[4]) # Name of the multivariate GWAS
# data_path <- "/Volumes/FLO/PhD/1_CHAPTER1/1_GENETICS/chapter1/"
# data_file <- "PC1_10.mvplink.50k.5k.txt.gz"
# figure_path <- "/Volumes/FLO/PhD/1_CHAPTER1/1_GENETICS/chapter1/"
# pc <- "PC1_10"
# threshold <- 1.7


# -------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


prep_file <- function(f,path) {
  
  # Function to open and prepare file's columns 
  # relative to the type of GWAS file (50k window-averaged or snp-by-snp)

  # Identify window averaged or snp-by-snp file
  if (grepl(".50k.5k.txt.gz", f, fixed = TRUE)) {
    pattern <- "50k"
  } else if (grepl(".logarithm.txt.gz", f, fixed = TRUE)) {
    pattern <- "logarithm"
  }
  
  # Prepare file and needed columns
  if (pattern == "50k") {
    run_files <- f %>% str_sub(.,end=-23) %>%
                       str_replace(.,pattern = '([a-z]{3})-([a-z]{3})-([a-z]{3})', '\\2\\1-\\3\\1')
    
    d <- read.table(paste0(path,f), header=TRUE) %>% left_join(hypo_chrom_start) %>% 
                                                     mutate(RUN = run_files,
                                                            LOG_P = AVG_P,
                                                            GPOS = MID_POS + GSTART)

  } else if (pattern == "logarithm") {
    run_files <- f %>% str_sub(.,end=-26) %>%
                       str_replace(.,pattern = '([a-z]{3})-([a-z]{3})-([a-z]{3})', '\\2\\1-\\3\\1')
    
    d <- read.table(paste0(path,f), header=TRUE) %>% left_join(hypo_chrom_start) %>%
                                                     mutate(RUN = run_files,
                                                            LOG_P = P,
                                                            MID_POS = POS,
                                                            BIN_START = POS,
                                                            BIN_END = POS,
                                                            RANGE = paste(CHROM, POS, sep="_"))

  } else {
    print("wrong file")
  }
  
  # print(run_files)
  # print(head(d))
  return(d)
  
}


threshold_table <- function(f) {
  
  # Create a table with regions of high association signal to plot

  assoc <- f  %>% select(CHROM, BIN_START, BIN_END, LOG_P, RUN) %>%
                  setNames(., nm = c('chrom', 'start', 'end', 'log_p', 'run') ) %>%
                  left_join(.,hypogen::hypo_chrom_start, by = c('chrom' = 'CHROM')) %>%
                  mutate(gpos = (start+end)/2 + GSTART) %>%
                  group_by(chrom) %>%
                  mutate(check = gpos > lag(gpos,default = 0) + 50000,
                         id = cumsum(check),
                         gid = str_c(chrom,'_',id)) %>%
                  group_by(gid) %>%
                  dplyr::summarise(chrom = chrom[1],
                                   start = min(start),
                                   end = max(end),
                                   gstart = min(start)+GSTART[1],
                                   gend = max(end)+GSTART[1]) %>%
                  mutate(gpos = (gstart+gend)/2)
  
  return(assoc)  
  
}


save_snp <- function (lg, start, end, outlier_id, ...) {
  
  # Create table with most associated SNPs
  
  data = gxp_snp %>% filter(CHROM == lg, MID_POS > start + 1500 & MID_POS < end - 1500) %>%
                     filter(row_number(desc(LOG_P)) <= 1) %>%
                     dplyr::summarise(CHROM = CHROM, POS = POS, LOG_P = LOG_P, WEIGHTS = WEIGHTS,
                                      RUN = RUN, RANGE = RANGE, GSTART = GSTART, MID_POS = MID_POS,
                                      BIN_START = BIN_START, BIN_END = BIN_END) %>%
                     mutate(ID = outlier_id)
  
}


save_snp_table <- function(table, path, name, t, chromosome) {
  
  # Save the most associated SNPs of regions of high association
  
  a <- table %>% filter(outlier_id %in% outlier_pick) %>%
                 pmap(save_snp) %>%
                 do.call(rbind, .)

  write.table(a, file = paste0(path, name, "/", name, "_", chromosome, ".snp.txt"), sep = " ",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  
}


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Read file and prepare useful columns
gxp_data <- prep_file(data_file,data_path)
data_snp <- paste0(pc,".mvplink.logarithm.txt.gz") 
gxp_snp <- prep_file(data_snp,data_path)

# Make a list of LG (Linkage Groups aka chromosome) names
chr <- (unique(as.vector(gxp_data[c("CHROM")])))
# print(chr["CHROM"])
# print(chr[[1]][1])

# Loop over each LG in the hamlet genome
for (i in 1:length(chr[[1]])) {
  # print(i)
  chrom <- chr[[1]][i] # Get the chromosome name
  # print(chrom)
  
  # Change threshold for LG08 (because messy chromosome)
  if (chrom == "LG08") {
    threshold <- 1.9
  } else {
    threshold <- 1.5
  }
  
  # print(threshold)
  
  # Select LG regions where the association signal is the highest
  thresh <- gxp_data[gxp_data[,"CHROM"] == chrom,]
  thresh <- thresh[thresh[, "LOG_P"] >= threshold,]
  # print(thresh)
  
  # Create table with regions of interest
  region_table <- threshold_table(thresh) %>% 
                  setNames(., nm = c("outlier_id", "lg", "start", "end", "gstart", "gend", "gpos")) %>%
                  mutate(heatmap = paste0(figure_path,outlier_id,"_heatmaps.png"))
  
  # List outlier regions' names
  outlier_pick <- region_table$outlier_id
  # print(outlier_pick)
  
  nb <- length(outlier_pick) # count nb of regions 
  # print(nb)
  
  region_label_tibbles <- tibble(outlier_id = outlier_pick, label = letters[1:nb]) # prepare panels letters for plots

  # Plot and save as table GWAS peaks for the LG
  if (nb != 0) {
    save_snp_table(region_table, data_path, pc, threshold, chrom)
  } else {
    print(paste0(chrom, " does not have peaks"))
  }
  
}
