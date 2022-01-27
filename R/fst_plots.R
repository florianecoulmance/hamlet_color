# by: Floriane Coulmance: 01/05/2020
# usage:
# Rscript fst_plots.R <fst_directory_path> <global_fst_file_path> <figure_path>
# -------------------------------------------------------------------------------------------------------------------
# fst_directory_path : $BASE_DIR/outputs/8_fst/
# global_fst_file_path : $BASE_DIR/outputs/8_fst/fst_globals.txt
# figure_path : $BASE_DIR/figures/fst/
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed library
library(ggplot2)
library(tidyverse)
library(hypogen)
library(hypoimg)
library(stringr)
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


pairwise_table <- function(path,file_list) {

  # Takes a list of files and combine their content into 1 table

  # Create list of dataset names
  run_files <- file_list %>%
               str_extract(., "[^.]+") %>% # extract string before first occurence of "." in file name
               str_replace(.,pattern = '([a-z]{3})-([a-z]{3})-([a-z]{3})', '\\2\\1-\\3\\1') # replace _ by -

  # Combine data from all dataset and keep track of dataset for each row
  data <- purrr::pmap(tibble(file = str_c(path,file_list),run = run_files),hypo_import_windows) %>%
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


subset_global <- function(table, param) {

  # Takes a global fst table and subset according to datasets wanted
  
  global1 <- if(param=="loc") {
                table %>% filter(dataset != "nig" & dataset != "pue")
             } else {
                table %>% filter(dataset != "bel" & dataset != "boc" & dataset != "puer")
             }

  return(global1)

}


anno_pair_spe <- function (spe, left, right,
                           circle_color = NA,
                           circle_fill_left = "white",
                           circle_fill_right = "lightgray",
                           circle_lwd = 0.5, plot_names = FALSE,
                           plot_name_size = 3,
                           font_family = "sans", ...) {
  
  # Modified function from GenomicOriginsScripts to get 2 flags around 1 species

  nr_left <- which(hypo_flag$geo %>% str_sub(.,1,3) == left) # select logo to be at the left of the plot in database
  
  nr_right <- which(hypo_flag$geo %>% str_sub(.,1,3) == right) # select logo to be at the right of the plot in database
  
  nr_spe <- which((hypo_img$spec %>% str_sub(.,1,3)) == spe) # select the logo of the species to be in the middle, from database
  
  # Combine each components into one big logo plot
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


legend_creation <- function(table_global, param) {
  
  # Creates the legend for each panel with logos
  
  # Create plot for each dataset into and store plots into a list
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
                anno_pair_spe(spe ="nig",left= "bel",right="pue") %>%
                  ggplotGrob(),
                anno_pair_spe(spe ="nig",left= "pan",right="bel") %>%
                  ggplotGrob(),
                anno_pair_spe(spe ="nig",left= "pan",right="pue") %>%
                  ggplotGrob(),
                anno_pair_spe(spe ="nig",left= "usa",right="bel") %>%
                  ggplotGrob(),
                anno_pair_spe(spe ="nig",left= "usa",right="pan") %>%
                  ggplotGrob(),
                anno_pair_spe(spe ="pue",left= "bel",right="pue") %>%
                  ggplotGrob(),
                anno_pair_spe(spe ="pue",left= "pan",right="bel") %>%
                  ggplotGrob(),
                anno_pair_spe(spe ="pue",left= "pan",right="pue") %>%
                  ggplotGrob(),
                anno_pair_spe(spe ="pue",left= "usa",right="bel") %>%
                  ggplotGrob(),
                anno_pair_spe(spe ="pue",left= "usa",right="pan") %>%
                  ggplotGrob(),
                anno_pair_spe(spe ="pue",left= "usa",right="pue") %>%
                  ggplotGrob(),
                GenomicOriginsScripts::anno_pair_flag(loc = "pue", left = "chl", right = "pue") %>%
                  ggplotGrob(),
                GenomicOriginsScripts::anno_pair_flag(loc = "pue", left = "chl", right = "uni") %>%
                  ggplotGrob(),
                GenomicOriginsScripts::anno_pair_flag(loc = "pue", left = "pue", right = "uni") %>%
                  ggplotGrob())
  
  # Put into a table the dataset names and the associated plot
  grob_list <- tibble(dataset = table_global$dataset, run = table_global$run, grob = grobs)
  
  # Subset according to parameter
  list_logos <- if(param=="loc") {
                  grob_list %>% filter(dataset != "nig" & dataset != "pue")
                } else {
                  grob_list %>% filter(dataset != "bel" & dataset != "boc" & dataset != "puer")
                }

  return(list_logos)
  
}


rescale_fst <- function (fst) {

  # Modified function from GenomicOriginsScripts used internaly by the GenomicOriginsScripts::fst_bar_row_run function

  start <- 0
  end <- 1
  fst_max <- max(global_table$weighted_fst)
  scales::rescale(fst, from = c(0, fst_max), to = c(start, end))

}


fst_plots <- function(table_fst, table_global, list_grob, path, prefix) {
  
  # Take fst data table and global fst table to make a faceted plot
  
  sc_ax <- scales::cbreaks(c(0,max(table_global$weighted_fst)), scales::pretty_breaks(4)) # format x axis preliminary to plot creation
  
  # Create the plot
  p <- ggplot() +
       facet_wrap(.~run, as.table = TRUE, ncol = 1,dir = 'v') +
       geom_rect(data = table_global %>% 
                        select(weighted_fst, run) %>%
                        setNames(., nm = c('fst', 'run')) %>%
                        pmap(., GenomicOriginsScripts::fst_bar_row_run) %>%
                        bind_rows() %>%
                        mutate(run = fct_reorder(run,xmax_org)) %>%
                        mutate(xmax = xmax * hypo_karyotype$GEND[23]),
                 aes(xmin = 0, xmax = xmax,
                 ymin = -Inf, ymax = Inf),
                 color = rgb(1,1,1,0),
                 fill = GenomicOriginsScripts::clr_below) +
       geom_vline(data = hypogen::hypo_karyotype,
                  aes(xintercept = GEND),
                  color = hypo_clr_lg) +
       geom_point(data = table_fst,
                  aes(x = GPOS, y = WEIGHTED_FST),
                  size=.2) +
       geom_hypo_grob(data = list_grob,
                      aes(grob = grob, x = .9,y = .7),
                      angle = 0, height = .5, width = .16) +
       scale_x_hypo_LG(sec.axis =  sec_axis(~ ./hypo_karyotype$GEND[23],
                                            breaks = (sc_ax$breaks/max(table_global$weighted_fst)),
                                            labels = sprintf("%.2f", sc_ax$breaks),
                                            name = expression(Genomic~position/~Genome~wide~weighted~italic(F[ST])))) +
       scale_y_continuous(name = expression(italic('F'[ST])),
                          limits = c(-.1,1),
                          breaks = c(0,.5,1)) +
       theme_hypo() +
       theme(strip.text = element_blank(),
             legend.position = 'none',
             axis.title.x = element_text(),
             axis.text.x.bottom = element_text(colour = 'darkgray'))
  
  # Save plot into png image
  hypo_save(filename = paste0(path,prefix,"_fst.png"),
            plot = p,
            width = 8,
            height = 12,
            dpi = 600,
            type = "cairo")
  
}


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Locate all the files that will be used in plots
files <- dir(fst_folder, pattern = '.50k.windowed.weir.fst.gz')
print(files)

# Remove empty files from list of files
for(f in files){
  print(f)
  x <- system(paste0("wc -l ",fst_folder,f),intern=TRUE)
  print(x)
  if(grepl("^1",x)){
  files <- files[files != f]
  }
}

# print(files)

# Locate the files either corresponding to a location fst or a species fst
files_loc <- files[grepl("^bel_|^boc_|^puer_", files)]
# print(files_loc)
files_spec <- files[grepl("^nig_|^pue_", files)]
# print(files_spec)

# Create tables for pairwise fst comparisons
pairwise_loc <- pairwise_table(fst_folder,files_loc)
# print(pairwise_loc)
pairwise_spec <- pairwise_table(fst_folder,files_spec)
# print(pairwise_spec)

# Open the global fst file
global <- read_tsv(fst_global,col_names=c("names","weighted_fst","mean_fst"))
# print(global)

# Format global fst table
global_table <- fst_global_table(global)
# print(global_table)

# Subset global_table according to the different fst plots we want
global_loc <- subset_global(global_table, "loc")
print(global_loc)
global_spec <- subset_global(global_table, "spec")
print(global_spec)

# Create a list of logos to use in faceted plot in next step
logos_loc <- legend_creation(global_table, "loc")
# print(logos_loc)
logos_spec <- legend_creation(global_table, "spec")
# print(logos_spec)

# Create fst plots
fst_plots(pairwise_loc, global_loc, logos_loc, path_figures, "loc")
fst_plots(pairwise_spec, global_spec, logos_spec, path_figures, "spec")
