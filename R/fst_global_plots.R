#!/usr/bin/env Rscript
# by: Floriane Coulmance: 08/06/2022
# usage:
# Rscript .R <multi_fst.50k> <fst_outlier> <global_fst>
# -------------------------------------------------------------------------------------------------------------------
# multi_fst.50k in : $BASE_DIR/outputs/8_fst/$DATASET/multi_fst.50k.tsv.gz
# fst_outlier in : $BASE_DIR/outputs/8_fst/$DATASET/fst_outliers_998.tsv
# global_fst in : $BASE_DIR/outputs/8_fst/$DATASET/fst_globals.txt
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

# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6:8]
print(args)

fst_50 <- as.character(args[1]) # Path to mean over 50k joint fst
# fst_50 <- "/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/chapter1/multi_fst.50k.tsv.gz"
fst_out <- as.character(args[2]) # Path to the outlier fst file
# fst_out <- "/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/chapter1/fst_outliers_998.tsv"
fst_glob <- as.character(args[3]) # Path to the global fst file
# fst_glob <- "/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/chapter1/fst_globals.txt"


# -------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


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
       annotation_custom(hypo_flag$flag[[nr_right]], xmin = 0.3, xmax = 0.9, ymin = -Inf, ymax = Inf) +
       annotation_custom(hypo_flag$flag[[nr_left]], xmin = -0.9, xmax = -0.3, ymin = -Inf, ymax = Inf) +
       annotation_custom(hypo_img$l[[nr_spe]], xmin = -0.6, xmax = 0.6, ymin = -Inf, ymax = Inf) +
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
  grob_list <- tibble(label = table_global$label, run = table_global$run, grob = grobs)
  
  # Subset according to parameter
  list_logos <- if(param=="loc") {
    grob_list %>% filter(run != "nig" & run != "pue")
  } else {
    grob_list %>% filter(run != "bel" & run != "boc" & run != "puer")
  }
  
  return(list_logos)
  
}


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Open 50k windowed Fst
fst_data <- vroom(fst_50, delim = '\t') %>%
  select(CHROM, BIN_START, BIN_END, N_VARIANTS, WEIGHTED_FST) %>%
  setNames(., nm = c('CHROM', 'BIN_START', 'BIN_END', 'n_snps', 'fst') ) %>%
  add_gpos() %>%
  select(GPOS, fst) %>%
  setNames(., nm = c('GPOS','value')) %>%
  mutate(window = str_c('bold(',project_case('a'),'):joint~italic(F[ST])'))

# Open global Fst
globals <- vroom(fst_glob, delim = '\t',
                        col_names = c('path','mean','weighted')) %>%
  mutate(label = str_c(str_sub(path,-11,-9),'_',str_sub(path,-7,-5),'-',str_sub(path,-3,-1)),
         label = fct_reorder(label,weighted),
         run = str_c((str_sub(label,1,3))),
         pair1 = str_c(str_sub(label,5,7)),
         pair2 = str_c(str_sub(label,9,11)),
         pair2 = ifelse(run == "nig" & pair2 == "pue", "puer", pair2),
         pair2 = ifelse(run == "pue" & pair2 == "pue","puer", pair2),
         label = gsub("uer", "puer", label),
         run = gsub("uer", "puer", run),
         pair_1 = str_c(run,"_",pair1),
         pair_2 = str_c(run,"_",pair2)) %>%
  filter(weighted != "NaN")

# Import fst outliers
outliers <-  vroom(fst_out, delim = '\t')

# Create the list of logos to use in plot 
logos_loc <- legend_creation(globals, "loc")
logos_spec <- legend_creation(globals, "spec")
# print(logos_spec)
# print(logos_loc)

# Merge the logos for 1 location species comparison to Fst table
d_loc <- merge(globals, logos_loc, by = c("label")) #%>% arrange(weighted)

# Merge the logos for 1 species locations comparison to Fst table
d_spec <- merge(globals, logos_spec, by = c("label")) #%>% arrange(weighted)

# Merge the complete list of comparisons logos to the Fst table
all <- rbind(logos_loc,logos_spec)
d_all <- merge(globals, all, by=c("label")) 


# Create the whole genome Fst plot
p1 <- ggplot() +
      ggtitle("a") +
      # add gray/white LGs background
      geom_hypo_LG()+
      # the red highlights for the outlier regions are added
      geom_vline(data = outliers,
                 aes(xintercept = gpos),
                 color = outlr_clr) +
      # the fst, delta dxy and gxp data is plotted
      geom_point(data = fst_data,
                 aes(x = GPOS, y = value),
                 size = plot_size,
                 color = plot_clr) +
      # setting the scales
      scale_y_continuous(name = expression(bolditalic(joint(F[ST]))),
                                expand = c(0, 0),
                                limits = c(0, 0.5)) +
      scale_x_continuous(name = "Linkage Group",
                         expand = c(0, 0),
                         breaks = (hypo_karyotype$GSTART + hypo_karyotype$GEND)/2,
                         labels = 1:24,
                         position = "top") +
      scale_fill_hypo_LG_bg() +
      # theme_hypo() +
      theme(plot.background = element_blank(), 
            panel.background = element_blank(),
            panel.grid = element_blank(), 
            panel.border = element_blank(),
            axis.line = element_line(),
            strip.background = element_rect(fill = NA, color = hypo_clr_lg),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.key = element_rect(fill = "transparent",  color = NA),
            plot.title = element_text(size = 30, face = "bold"),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 20),
            plot.margin = unit(c(0.5, 0.5, 0.1, 0.5), "cm"))



p_loc <- ggplot(d_loc,
                aes(x=0, y=-2)) +
         ggtitle("b") +
         facet_wrap(.~label,
                    as.table = TRUE,
                    ncol = 3,
                    dir = 'v') +
         geom_text(aes(label = weighted),
                   color = "darkgrey",
                   nudge_x = 0,
                   nudge_y = 2,
                   size = 8) +
         geom_text(aes(label = mean),
                   color = "#5B9E2E",
                   nudge_x = 0,
                   size = 8) +
         geom_hypo_grob(data = d_loc,
                        aes(grob = grob, x = 0.2, y = 0.6),
                        angle = 0,
                        height = 9,
                        width = 0.9) +
         scale_x_continuous(limits = c(-1,1),
                            name = "Between species",
                            position = "top") +
         scale_y_continuous(limits = c(-2.8,10)) +
         theme_hypo() +
         theme(strip.text = element_blank(),
               legend.position = 'none',
               legend.title=element_blank(),
               plot.title = element_text(hjust = -0.09, size = 30, face = "bold"),
               axis.title.x.top = element_text(size = 20),
               # axis.title.x.top = element_blank(),
               axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               panel.grid.major = element_blank(),
               axis.title = element_blank(),
               axis.ticks = element_blank(),
               axis.line = element_blank(),
               panel.border = element_blank(),
               plot.margin = unit(c(-1.5, 0, 0.5, 5), "cm"))


p_spec <- ggplot(d_spec,
                 aes(x=0, y=-2)) +
          ggtitle("c") +
          facet_wrap(.~label,
                     as.table = TRUE,
                     ncol = 3,
                     dir = 'v') +
          geom_text(aes(label = weighted),
                    color = "darkgrey",
                    nudge_x = 0,
                    nudge_y = 2,
                    size = 8) +
          geom_text(aes(label = mean),
                    color = "#5B9E2E",
                    nudge_x = 0,
                    size = 8) +
          geom_hypo_grob(data = d_spec,
                         aes(grob = grob, x = 0.2, y = 0.6),
                         angle = 0,
                         height = 8,
                         width = 0.8) +
          scale_x_continuous(limits = c(-1,1),
                             name = "Between locations",
                             position = "top") +
          scale_y_continuous(limits = c(-2.8,10)) +
                             theme_hypo() +
          theme(strip.text = element_blank(),
                legend.position = 'none',
                legend.title=element_blank(),
                plot.title = element_text(hjust = -0.05, size = 30, face = "bold"),
                axis.title.x.top = element_text(size = 20),
                # axis.title.x.top = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                panel.grid.major = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank(),
                panel.border = element_blank(),
                plot.margin = unit(c(-1.5, 0, 0.5, 6), "cm"))

g <- ggarrange(p1, ggarrange(p_loc, p_spec, ncol = 2, widths = c(1,1)), ncol = 1, nrow = 2, align = "h",
               widths = c(2, 2), heights = c(1, 4)) 



hypo_save(filename = "/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/chapter1/fst_global.pdf",
          plot = g,
          width = 28,
          height = 12)





