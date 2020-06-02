rm(list = ls())


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



setwd('/Users/fco/Desktop/BREMEN_OP/chapter1_2/FST/')

files <- dir(pattern = "50k.windowed.weir.fst.gz")

run_files <- files %>%
  str_sub(.,1,11) %>%
  str_replace(.,pattern = '([a-z]{3})-([a-z]{3})-([a-z]{3})', '\\2\\1-\\3\\1')

data <- purrr::pmap(tibble(file = str_c(files),run = run_files),hypo_import_windows) %>% bind_rows() %>%
  set_names(., nm = c('CHROM', 'BIN_START', 'BIN_END', 'N_VARIANTS', 'WEIGHTED_FST', 'MEAN_FST', 'GSTART', 'POS', 'GPOS', 'run')) %>% separate(run, into =c("loc","pop1","pop2")) %>% mutate(run = str_c(pop1,loc,'-',pop2,loc),run = fct_reorder(run,WEIGHTED_FST))
                                                                                                                                                                                             

globals_file <- read_tsv("fst_globals.txt",col_names=c("names","weighted_fst","mean_fst"))
globals_file$names <- gsub("^.{0,44}", "", globals_file$names)
globals_file <- globals_file %>% separate(names, c("loc","pop1", "pop2"), sep = "_", extra = "merge")
globals_file$run <- paste(globals_file$pop1,globals_file$pop2,sep="_")
globals_file <- globals_file[, c(1,6,4,5)]
globals_file <- globals_file[-c(15), ] 
write.table(globals_file,"fst_globals2.txt",sep="\t",row.names=FALSE,col.names=FALSE)

globals <- vroom::vroom("fst_globals2.txt", delim = '\t',
                        col_names = c('loc','run','mean','weighted')) %>%
  separate(run, into = c('pop1','pop2')) %>%
  mutate(run = str_c(pop1,loc,'-',pop2,loc),
         run = fct_reorder(run,weighted))

runs <- globals %>%
  group_by(run) %>%
  count() %>%
  ungroup() %>%
  #select(-n) %>%
  mutate(loc = str_sub(run,4,6),
         right_short = str_sub(run,1,3),
         left_short = str_sub(run,8,10)) %>%
  mutate(left = left_short,
         right = right_short)

grobs <- list(hypo_anno_pair('indigo','maya',
                             circle_color = 'black',plot_names = TRUE) %>%
                ggplotGrob(),
              hypo_anno_pair('indigo','nigricans',
                             circle_color = 'black',plot_names = TRUE) %>%
                ggplotGrob(),
              hypo_anno_pair('indigo','puella',
                             circle_color = 'black',plot_names = TRUE) %>%
                ggplotGrob(),
              hypo_anno_pair('maya','nigricans',
                             circle_color = 'black',plot_names = TRUE) %>%
                ggplotGrob(),
              hypo_anno_pair('maya','puella',
                             circle_color = 'black',plot_names = TRUE) %>%
                ggplotGrob(),
              hypo_anno_pair('nigricans','puella',
                             circle_color = 'black',plot_names = TRUE) %>%
                ggplotGrob(),
              hypo_anno_pair('nigricans','puella',
                             circle_color = 'black',plot_names = TRUE) %>%
                ggplotGrob(),
              hypo_anno_pair('nigricans','unicolor',
                             circle_color = 'black',plot_names = TRUE) %>%
                ggplotGrob(),
              hypo_anno_pair('puella','unicolor',
                             circle_color = 'black',plot_names = TRUE) %>%
                ggplotGrob(),
              hypo_anno_flag_pair(left= 'belize', right = 'puerto_rico',
                                 flag_lwd = 1, flag_line_color = 'black',
                                 circle_color = 'black', plot_names = TRUE,
                                 plot_name_size = 5) %>%
                ggplotGrob(),
              hypo_anno_flag_pair(left= 'panama', right = 'belize',
                                  flag_lwd = 1, flag_line_color = 'black',
                                  circle_color = 'black', plot_names = TRUE,
                                  plot_name_size = 5) %>%
                ggplotGrob(),
              hypo_anno_flag_pair(left= 'panama', right = 'puerto_rico',
                                  flag_lwd = 1, flag_line_color = 'black',
                                  circle_color = 'black', plot_names = TRUE,
                                  plot_name_size = 5) %>%
                ggplotGrob(),
              hypo_anno_flag_pair(left= 'usa', right = 'belize',
                                  flag_lwd = 1, flag_line_color = 'black',
                                  circle_color = 'black', plot_names = TRUE,
                                  plot_name_size = 5) %>%
                ggplotGrob(),
              hypo_anno_flag_pair(left= 'usa', right = 'panama',
                                  flag_lwd = 1, flag_line_color = 'black',
                                  circle_color = 'black', plot_names = TRUE,
                                  plot_name_size = 5) %>%
                ggplotGrob(),
              hypo_anno_flag_pair(left= 'belize', right = 'puerto_rico',
                                  flag_lwd = 1, flag_line_color = 'black',
                                  circle_color = 'black', plot_names = TRUE,
                                  plot_name_size = 5) %>%
                ggplotGrob(),
              hypo_anno_flag_pair(left= 'panama', right = 'belize',
                                  flag_lwd = 1, flag_line_color = 'black',
                                  circle_color = 'black', plot_names = TRUE,
                                  plot_name_size = 5) %>%
                ggplotGrob(),
              hypo_anno_flag_pair(left= 'panama', right = 'puerto_rico',
                                  flag_lwd = 1, flag_line_color = 'black',
                                  circle_color = 'black', plot_names = TRUE,
                                  plot_name_size = 5) %>%
                ggplotGrob(),
              hypo_anno_flag_pair(left= 'usa', right = 'belize',
                                  flag_lwd = 1, flag_line_color = 'black',
                                  circle_color = 'black', plot_names = TRUE,
                                  plot_name_size = 5) %>%
                ggplotGrob(),
              hypo_anno_flag_pair(left= 'usa', right = 'panama',
                                  flag_lwd = 1, flag_line_color = 'black',
                                  circle_color = 'black', plot_names = TRUE,
                                  plot_name_size = 5) %>%
                ggplotGrob(),
              hypo_anno_flag_pair(left= 'usa', right = 'puerto_rico',
                                  flag_lwd = 1, flag_line_color = 'black',
                                  circle_color = 'black', plot_names = TRUE,
                                  plot_name_size = 5) %>%
                ggplotGrob(),
              hypo_anno_pair('chlorurus','puella',
                             circle_color = 'black',plot_names = TRUE) %>%
                ggplotGrob(),
              hypo_anno_pair('chlorurus','unicolor',
                             circle_color = 'black',plot_names = TRUE) %>%
                ggplotGrob(),
              hypo_anno_pair('puella','unicolor',
                             circle_color = 'black',plot_names = TRUE) %>%
                ggplotGrob())

grob_list <- tibble(RUN = runs$run, grob = grobs)



sc_ax <- scales::cbreaks(c(0,max(globals$weighted)),
                         scales::pretty_breaks(4))

p <- ggplot() + facet_wrap(. ~ run, ncol = 1,dir = 'v') +
  geom_vline(data = hypogen::hypo_karyotype,
             aes(xintercept = GEND),
             color = hypo_clr_lg) +
  # 
  # geom_hypo_grob(data = grob_list,
  #                aes(grob = grob, x = .9,y = .7),
  #                angle = 0, height = .5, width = .16) +
  
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


print(p)

