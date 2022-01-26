# by: Floriane Coulmance: 01/05/2020
# usage:
# Rscript fst_plots.R <fst_directory_path> <global_fst_file_path>
# -------------------------------------------------------------------------------------------------------------------
# fst_directory_path : $BASE_DIR/outputs/8_fst/
# global_fst_file_path : $BASE_DIR/outputs/8_fst/fst_globals.txt
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed library
library(ggplot2)
library(tidyverse)
library(hypogen)
library(hypoimg)
library(dplyr)
library(plyr)
library(ggpubr)
library(vroom)


# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6:8]
print(args)

fst_folder <- as.character(args[1]) # Path to fst folder
print(fst_folder)
fst_global <- as.character(args[2]) # Path to the global fst file
print(fst_global)
path_figures <- as.character(args[3]) # Path to the figure folder
print(path_figures)


# -------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


pairwise_table <- function(file_list) {

  # Takes a list of files and combine their content into 1 table

  run_files <- file_list %>%
             str_sub(.,1,11) %>%
             str_replace(.,pattern = '([a-z]{3})-([a-z]{3})-([a-z]{3})', '\\2\\1-\\3\\1')

  data <- purrr::pmap(tibble(file = str_c(file_list),run = run_files),hypo_import_windows) %>%
          bind_rows() %>%
          set_names(., nm = c('CHROM', 'BIN_START', 'BIN_END', 'N_VARIANTS', 'WEIGHTED_FST', 'MEAN_FST', 'GSTART', 'POS', 'GPOS', 'run')) %>%
          separate(run, into =c("loc","pop1","pop2")) %>% 
          mutate(run = str_c(pop1,loc,'-',pop2,loc),run = fct_reorder(run,WEIGHTED_FST))

  return(data)

}


fst_global_table <- function(file) {
  
  # Format global fst table with run columns, pair of species and location 
  
  global1 <- file %>% 
    mutate(names = gsub("^.{0,44}", "", names)) %>% 
    separate(names, c("dataset","pair1", "pair2"), sep = "_", extra = "merge") %>%
    mutate(run = str_c(pair1,dataset,'-',pair2,dataset),
           run = fct_reorder(run,weighted_fst)) %>%
    relocate(dataset, run, pair1, pair2, weighted_fst, mean_fst) %>%
    filter(weighted_fst != "NaN" & mean_fst != "NaN")
  
  return(global1)
  
}


legend_creation <- function(table_global, param) {
  
  # Creates the legend for each panel with logos
  
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
  
  grob_list <- tibble(dataset = global1$dataset, RUN = global1$run, grob = grobs)
  
  list_logos <- if(param=="loc")
                grob_list %>% filter(dataset != "nig" & dataset != "pue")
                else
                grob_list %>% filter(dataset != "bel" & dataset != "boc" & dataset != "puer")
  
  return(list_logos)
  
}


fst_plots <- function(table_fst, table_global, list_grob, path, prefix) {
  
  # Take fst data table and global fst table to make a faceted plot
  
  sc_ax <- scales::cbreaks(c(0,max(table_global$weighted)),
                           scales::pretty_breaks(4))
  
  p <- ggplot() +
    facet_wrap(. ~ run, ncol = 1,dir = 'v') +
    geom_vline(data = hypogen::hypo_karyotype,
               aes(xintercept = GEND),
               color = hypo_clr_lg) +
    geom_hypo_grob(data = list_grob,
                   aes(grob = grob, x = .9,y = .7),
                   angle = 0, height = .5, width = .16) +
    geom_point(data = table_fst,
               aes(x = GPOS, y = WEIGHTED_FST),
               size=.2) +
    scale_x_hypo_LG(sec.axis =  sec_axis(~ ./hypo_karyotype$GEND[23],
                                         breaks = (sc_ax$breaks/max(table_global$weighted)),
                                         labels = sprintf("%.2f", sc_ax$breaks))) +
    scale_y_continuous(name = expression(italic('F'[ST])),
                       limits = c(-.1,1),
                       breaks = c(0,.5,1)) +
    theme_hypo() +
    theme(strip.text = element_blank(),
          legend.position = 'none',
          axis.title.x = element_text(),
          axis.text.x.bottom = element_text(colour = 'darkgray'))
  
  hypo_save(filename = paste0(path,prefix,"_fst.pdf"),
            plot = f,
            width = 26,
            height = 40)
  
}


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------

# Locate all the files that will be used in plots
files <- dir(fst_folder, pattern = '.50k.windowed.weir.fst.gz')
print(files)

# Locate the files either corresponding to a location fst or a species fst
files_loc <- files[grepl("^bel_|^boc_|^puer_", files)]
print(files_loc)
files_spec <- files[grepl("^nig_|^pue_", files)]
print(files_spec)

# Create tables for pairwise fst comparisons
pairwise_loc <- pairwise_table(files_loc)
print(pairwise_loc)
pairwise_spec <- pairwise_table(files_spec)
print(pairwise_spec)

# Open the global fst file
global <- read_tsv(fst_global,col_names=c("names","weighted_fst","mean_fst"))
print(global)

# Format global fst table
global_table <- fst_global_table(global)
print(global_table)

# Create a list of logos to use in faceted plot in next step
logos_loc <- legend_creation(global_table, "loc")
print(logos_loc)
logos_spec <- legend_creation(global_table, "spec")
print(logos_spec)

# Create fst plots
fst_plots(pairwise_loc, global_table, logos_loc, "loc")
fst_plots(pairwise_spec, global_table, logos_spec, "spec")


