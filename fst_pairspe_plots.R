rm(list = ls())

#needed libraries
library(ggplot2)
library(tidyverse)
library(hypogen)
library(hypoimg)
library(usethis)
library(devtools)
library(dplyr)
library(plyr)
library(ggpubr)
library(vroom)
library(GenomicOriginsScripts)



#modified function from GenomicOriginsScripts to get 2 flags around 1 species
anno_pair_spe <- function (spe,left, right, circle_color = NA, circle_fill_left = "white",
                           circle_fill_right = "lightgray", circle_lwd = 0.5, plot_names = FALSE,
                           plot_name_size = 3, font_family = "sans", ...) {
  nr_left <- which(hypo_flag$geo %>% str_sub(.,1,3) == left)
  nr_right <- which(hypo_flag$geo %>% str_sub(.,1,3) == right)
  nr_spe <- which((hypo_img$spec %>% str_sub(.,1,3)) == spe)
  p <- ggplot() +
    theme_void() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.4, 0.38)) +
    annotation_custom(hypo_flag$flag[[nr_right]], xmin = 0.5, xmax = 1.02, ymin = -Inf, ymax = Inf) +
    annotation_custom(hypo_flag$flag[[nr_left]], xmin = -1.02, xmax = -0.5, ymin = -Inf, ymax = Inf)+
    annotation_custom(hypo_img$l[[nr_spe]], xmin = -0.3, xmax = .3, ymin = -Inf, ymax = Inf) +
    coord_cartesian(xlim = c(-1.1,1.1))
  return(p)
}


#set all needed files in right format for plots
setwd('/Users/fco/Desktop/BREMEN_OP/chapter1_2/FST/pairwise_species/')
files <- dir(pattern = "50k.windowed.weir.fst.gz")

#list of names of files corresponding to fst combinations
run_files <- files %>%
  str_sub(.,end=-26) %>%
  str_replace(.,pattern = '([a-z]{3})-([a-z]{3})-([a-z]{3})', '\\2\\1-\\3\\1')


#table assembling all the data tables
data <- purrr::pmap(tibble(file = str_c(files),run = run_files),hypo_import_windows) %>% bind_rows() %>%
  set_names(., nm = c('CHROM', 'BIN_START', 'BIN_END', 'N_VARIANTS', 'WEIGHTED_FST', 'MEAN_FST', 'GSTART', 'POS', 'GPOS', 'run')) %>% separate(run, into =c("loc","pop1","pop2")) %>% mutate(run = str_c(pop1,loc,'-',pop2,loc),run = fct_reorder(run,WEIGHTED_FST))

#subset for test : run faster                                                                                                                                                                                          
dat1 <- data %>% filter(CHROM == "LG01")

#transform the run names correspponding to fst associations
globals_file <- read_tsv("../fst_globals.txt",col_names=c("names","weighted_fst","mean_fst"))
globals_file$names <- gsub("^.{0,44}", "", globals_file$names)
globals_file <- globals_file %>% separate(names, c("loc","pop1", "pop2"), sep = "_", extra = "merge")
globals_file$run <- paste(globals_file$pop1,globals_file$pop2,sep="_")
globals_file <- globals_file[, c(1,6,4,5)]
globals_file <- globals_file[-c(10:21), ] 
write.table(globals_file,"../fst_globals2.txt",sep="\t",row.names=FALSE,col.names=FALSE)

globals <- vroom::vroom("../fst_globals2.txt", delim = '\t',
                        col_names = c('loc','run','mean','weighted')) %>%
  separate(run, into = c('pop1','pop2')) %>%
  mutate(run = str_c(pop1,loc,'-',pop2,loc),
         run = fct_reorder(run,weighted))


# set the illustrations corresponding to each associations of pairwise fsts
grobs <- list(GenomicOriginsScripts::anno_pair_flag(loc = "bel", left = "ind", right = "may") %>%
                ggplotGrob(),
              GenomicOriginsScripts::anno_pair_flag(loc = "bel", left = "ind", right = "nig") %>%
                ggplotGrob(),
              GenomicOriginsScripts::anno_pair_flag(loc = "bel", left = "ind", right = "pue") %>%
                ggplotGrob(),
              GenomicOriginsScripts::anno_pair_flag(loc = "bel", left = "may", right = "nig") %>%
                ggplotGrob(),
              GenomicOriginsScripts::anno_pair_flag(loc = "bel", left = "may", right = "pue") %>%
                ggplotGrob(),
              GenomicOriginsScripts::anno_pair_flag(loc = "bel", left = "nig", right = "pue") %>%
                ggplotGrob(),
              GenomicOriginsScripts::anno_pair_flag(loc = "pan", left = "nig", right = "pue") %>%
                ggplotGrob(),
              GenomicOriginsScripts::anno_pair_flag(loc = "pan", left = "nig", right = "uni") %>%
                ggplotGrob(),
              GenomicOriginsScripts::anno_pair_flag(loc = "pan", left = "pue", right = "uni") %>%
                ggplotGrob(),
              GenomicOriginsScripts::anno_pair_flag(loc = "pue", left = "chl", right = "pue") %>%
                ggplotGrob(),
              GenomicOriginsScripts::anno_pair_flag(loc = "pue", left = "chl", right = "uni") %>%
                ggplotGrob(),
              GenomicOriginsScripts::anno_pair_flag(loc = "pue", left = "pue", right = "uni") %>%
                ggplotGrob())

grob_list <- tibble(run = globals$run, grob = grobs)



#plot
sc_ax <- scales::cbreaks(c(0,max(globals$weighted)),
                         scales::pretty_breaks(4))


p <- ggplot() + facet_wrap(run~., ncol = 1,dir = 'v') +
  geom_vline(data = hypogen::hypo_karyotype,
             aes(xintercept = GEND),
             color = hypo_clr_lg) +
  
  geom_hypo_grob(data = grob_list,
                 aes(grob = grob, x = .9,y = .7),
                 angle = 0, height = .5, width = .16) +
  
  geom_point(data = data, aes(x = GPOS, y = WEIGHTED_FST),
             size=.2) +
  
  scale_x_hypo_LG(sec.axis =  sec_axis(~ ./hypo_karyotype$GEND[23],
                                       breaks = (sc_ax$breaks/max(globals$weighted)),
                                       labels = sprintf("%.2f", sc_ax$breaks))) +
  
  scale_y_continuous(name = expression(italic('F'[ST])),
                     limits = c(-.1,1),
                     breaks = c(0,.5,1)) +
  theme_hypo() +
  
  theme(strip.text = element_blank(),
        legend.position = 'none',
        axis.title.x = element_text(),
        axis.text.x.bottom = element_text(colour = 'darkgray'))

hypo_save(filename = 'Pairwise_spec_Fst.png',
          plot = p,
          width = 8,
          height = 12)
