#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla table_fst_outliers.R multi_fst.50k.tsv.gz
# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed library
library(GenomicOriginsScripts)
library(vroom)


# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6:7]
print(args)

fst_file <- as.character(args[1])
path <- as.character(args[2])


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# outliers <-  vroom::vroom(fst_file,delim = '\t') %>%
#   select(CHROM, BIN_START, BIN_END, N_VARIANTS, WEIGHTED_FST) %>%
#   setNames(., nm = c('chrom', 'start', 'end', 'n_snps', 'fst') ) %>%
#   left_join(.,hypogen::hypo_chrom_start, by = c('chrom' = 'CHROM')) %>%
#   mutate(gpos = (start+end)/2 + GSTART) %>%
#   group_by(chrom) %>%
#   mutate(norm_fst = fst-mean(fst)) %>%
#   ungroup() %>%
#   filter(norm_fst >= quantile(norm_fst, .998)) %>%
#   group_by(chrom) %>%
#   mutate(check = gpos > lag(gpos,default = 0) + 50000,
#          id = cumsum(check),
#          gid = str_c(chrom,'_',id)) %>%
#   group_by(gid) %>%
#   summarise(chrom = chrom[1],
#             start = min(start),
#             end = max(end),
#             gstart = min(start)+GSTART[1],
#             gend = max(end)+GSTART[1]) %>%
#   mutate(gpos = (gstart+gend)/2)

# outliers <-  vroom::vroom(fst_file,delim = ' ') %>%
#   select(CHROM, BIN_START, BIN_END, AVG_P, RUN) %>%
#   setNames(., nm = c('chrom', 'start', 'end', 'avg_p', 'run') ) %>%
#   left_join(.,hypogen::hypo_chrom_start, by = c('chrom' = 'CHROM')) %>%
#   mutate(gpos = (start+end)/2 + GSTART) %>%
#   group_by(chrom) %>%
#   mutate(check = gpos > lag(gpos,default = 0) + 50000,
#          id = cumsum(check),
#          gid = str_c(chrom,'_',id)) %>%
#   group_by(gid) %>%
#   summarise(chrom = chrom[1],
#             start = min(start),
#             end = max(end),
#             gstart = min(start)+GSTART[1],
#             gend = max(end)+GSTART[1]) %>%
#   mutate(gpos = (gstart+gend)/2)


outliers <-  vroom::vroom(fst_file,delim = ' ') %>%
  select(CHROM, BIN_START, BIN_END, P, RUN) %>%
  setNames(., nm = c('chrom', 'start', 'end', 'p', 'run') ) %>%
  left_join(.,hypogen::hypo_chrom_start, by = c('chrom' = 'CHROM')) %>%
  mutate(gpos = (start+end)/2 + GSTART) %>%
  group_by(chrom) %>%
  mutate(check = gpos > lag(gpos,default = 0) + 50000,
         id = cumsum(check),
         gid = str_c(chrom,'_',id)) %>%
  group_by(gid) %>%
  summarise(chrom = chrom[1],
            start = min(start),
            end = max(end),
            gstart = min(start)+GSTART[1],
            gend = max(end)+GSTART[1]) %>%
  mutate(gpos = (gstart+gend)/2)

# write_tsv(outliers, paste0(path,'10kb_AVGP>3.tsv'))
# write_tsv(outliers, paste0(path,'fst_outliers_998.tsv'))
write_tsv(outliers, paste0(path,'snp_P13.tsv'))
