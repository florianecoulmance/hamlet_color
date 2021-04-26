# by: Floriane Coulmance: 15/04/2021
# usage:
# Rscript gxp_plots.R <data_path> <data_file> <figure_path> <pc>
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

data_path <- as.character(args[1]) # Path to mean coverage data table
data_file <- as.character(args[2])
# data_path <- "/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/chapter1/"
# data_file <- "PC1_8.mvplink.50k.5k.txt.gz"
figure_path <- as.character(args[3]) # Path to the figure folder
# figure_path <- "/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/chapter1/"
# pc <- "PC1_8"
pc <- as.character(args[4])


# -------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


prep_file <- function(f,path) {
  
  # Function to open and prepare columns 
  # relative to the type of association file (window averaged or snp-by-snp)


  # Identify window averaged or snp-by-snp file
  if (grepl(".50k.5k.txt.gz", f, fixed = TRUE)) {
  pattern <- "50k"
  } else if (grepl(".logarithm.txt.gz", f, fixed = TRUE)) {
  pattern <- "logarithm"
  }
  
  # Prepare file and needed columns
  if (pattern == "50k"){
    run_files <- f %>%
        str_sub(.,end=-23) %>%
        str_replace(.,pattern = '([a-z]{3})-([a-z]{3})-([a-z]{3})', '\\2\\1-\\3\\1')
    d <- read.table(paste0(path,f), header=TRUE) %>% left_join(hypo_chrom_start) %>% mutate(RUN = run_files, LOG_P = AVG_P, GPOS = MID_POS + GSTART)

  } else if (pattern == "logarithm") {
    run_files <- f %>%
        str_sub(.,end=-26) %>%
        str_replace(.,pattern = '([a-z]{3})-([a-z]{3})-([a-z]{3})', '\\2\\1-\\3\\1')
    d <- read.table(paste0(path,f), header=TRUE) %>% left_join(hypo_chrom_start) %>% mutate(RUN = run_files, LOG_P = P, MID_POS = POS, BIN_START = POS, BIN_END = POS, RANGE = paste(CHROM, POS, sep="_"))

  } else {
    print("wrong file")

  }
  
  print(run_files)
  print(head(d))
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


plot_panel_gxp <- function (lg, start, end, trait, ...) {
  # Create gxp plot for windowed results
  
  ggplot2::ggplot() +
  geom_rect(inherit.aes = FALSE, data = tibble(start = start,end = end), aes(xmin = start, xmax = end),
            ymin = -Inf, ymax = Inf, fill = rgb(0.9, 0.9, 0.9, 0.3), color = rgb(0.9, 0.9, 0.9, 0.9)) + 
  geom_line(data = gxp_data %>% filter(CHROM == lg, MID_POS > start - window_buffer * 1.25, MID_POS < end + window_buffer * 1.25), aes(x = MID_POS, y = LOG_P, color = RUN), size = 0.2) +
  scale_color_manual(name = "GxP Trait", values = gxp_clr) +
  scale_x_continuous(name = lg, expand = c(0, 0), position = "top") +
  scale_y_continuous(name = expression(bolditalic(-log[10](p))),expand = c(0, 0)) +
  guides(color = guide_legend(keyheight = unit(3,"pt"), keywidth = unit(20, "pt"), override.aes = list(size = 2))) +
  theme_panels()
  
}


plot_panel_gxp_snp <- function (lg, start, end, trait, ...) {
  # Plot gxp results snp-by-snp
  
  ggplot2::ggplot() +
    geom_rect(inherit.aes = FALSE, data = tibble(start = start,end = end), aes(xmin = start, xmax = end),
              ymin = -Inf, ymax = Inf, fill = rgb(0.9, 0.9, 0.9, 0.3), color = rgb(0.9, 0.9, 0.9, 0.9)) + 
    geom_line(data = gxp_snp %>% filter(CHROM == lg, MID_POS > start - window_buffer * 1.25, MID_POS < end + window_buffer * 1.25), aes(x = MID_POS, y = LOG_P, color = RUN, label = RANGE), size = 0.2) +
    geom_text_repel(data = gxp_snp %>% filter(CHROM == lg, MID_POS > start - window_buffer * 1.25, MID_POS < end + window_buffer * 1.25) %>%  filter(row_number(desc(LOG_P)) <= 5), aes(x = MID_POS, y = LOG_P, color = RUN, label = RANGE), segment.color = 'transparent', color = "Red", size=1, hjust = -1, nudge_x = -1.5) +
    scale_color_manual(name = "GxP Trait", values = gxp_clr) +
    scale_x_continuous(name = lg, expand = c(0, 0), position = "top") +
    scale_y_continuous(name = expression(bolditalic(-log[10](p))),expand = c(0, 0)) +
    guides(color = guide_legend(keyheight = unit(3,"pt"), keywidth = unit(20, "pt"), override.aes = list(size = 2))) +
    theme_panels()
  
}


save_snp <- function (lg, start, end, trait, ...) {
  # Create table with most associated SNPs of the most associated 5okb regions
  
  data = gxp_snp %>% filter(CHROM == lg, MID_POS > start - window_buffer * 1.25, MID_POS < end + window_buffer * 1.25) %>%
         filter(row_number(desc(LOG_P)) <= 5) %>%
         dplyr::summarise(CHROM = CHROM, POS = POS, LOG_P = LOG_P, WEIGHTS = WEIGHTS, RUN = RUN, RANGE = RANGE, GSTART = GSTART, MID_POS = MID_POS, BIN_START = BIN_START, BIN_END = BIN_END)
  
}


save_snp_table <- function(table, path, name) {
  # Save the most associated snps of regions of high association
  
  a <- table %>% filter(outlier_id %in% outlier_pick) %>%
    pmap(save_snp, cool_genes = cool_genes) %>% do.call(rbind, .)
  write.table(a, file = paste0(path,name,".snp.txt"), sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)
}


plot_curt <- function (loc = "bel", outlier_id, outlier_nr, lg, start, end, 
                       cool_genes, text = TRUE, label, trait, ...) 
{
  p_g <- plot_panel_anno(lg = lg, outlier_id = outlier_id, 
                         label = label, start = start, end = end, genes = cool_genes)

  p_gxp <- plot_panel_gxp(lg = lg, start = start, end = end,
                          trait = pcs)
  
  p_snp <- plot_panel_gxp_snp(lg = lg, start = start, end = end,
                     trait = pcs)

  if (text) {
    p_curtain <- cowplot::plot_grid(p_g, p_gxp, p_snp, ncol = 1, align = "v", 
                                    rel_heights = c(1, rep(0.8, 7)))
  }
  else {
    p_curtain <- cowplot::plot_grid(p_g + no_title(), p_gxp + 
                                      no_title(), p_snp + no_title(), ncol = 1, align = "v", 
                                    rel_heights = c(1, rep(0.8, 7)))
  }
  p_curtain
}


plot_regions <- function(region, outlier_list, label, path, name) {
  # Create plot for regions with highest associations

  # Create the plot
  p_single <- region %>%
    filter(outlier_id %in% outlier_list) %>%
    left_join(label) %>%
    dplyr::mutate(outlier_nr = row_number(), text = ifelse(outlier_nr == 1,TRUE,FALSE)) %>%
    pmap(plot_curt, cool_genes = cool_genes) %>%
    cowplot::plot_grid(plotlist = ., nrow = row,
                       rel_heights = heights,
                       labels = letters[1:length(outlier_list)] %>% project_case())
  
  # Save the file in .png file
  hypo_save(filename = paste0(path, name, ".png"),
            plot = p_single,
            width = 12,
            height = 8)  
}


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Read file and prepare useful columns
gxp_data <- prep_file(data_file,data_path)
data_snp <- paste0(pc,".mvplink.logarithm.txt.gz") 
gxp_snp <- prep_file(data_snp,data_path)

# Select regions where the association signal is the highest
thresh <- gxp_data[gxp_data[, "LOG_P"] >= 1.2 ,]

# Create table with regions of interest
region_table <- threshold_table(thresh) %>%
  setNames(., nm = c("outlier_id", "lg", "start", "end", "gstart", "gend", "gpos"))

outlier_pick <- region_table$outlier_id # find outlier regions of interest

outlier_pick <- outlier_pick[lapply(outlier_pick,function(x) length(grep("LG08",x,value=FALSE))) == 0]

print(outlier_pick)
nb <- length(outlier_pick) # count nb of regions 
region_label_tibbles <- tibble(outlier_id = outlier_pick, label = letters[1:nb]) # prepare panels letters for plots

# Set parameters for plots
cool_genes <-  c('arl3','kif16b','cdx1','hmcn2',
                 'sox10','smarca4',
                 'rorb',
                 'alox12b','egr1',
                 'ube4b','casz1',
                 'hoxc8a','hoxc9','hoxc10a',
                 'hoxc13a','rarga','rarg',
                 'snai1','fam83d','mafb','sws2abeta','sws2aalpha','sws2b','lws','grm8')

# set parameter of window buffer for plot
window_buffer <- 2.5*10^5

# set color parameter for plot line
gxp_clr <- c("#5B9E2E") %>%
  darken(factor = .95)

# Determine the dimension and grid layout regarding the nb of region to plot 
if (nb <= 4) {
  row = 1
  heights = c(1, 1)
} else if (nb <= 8 & nb > 4) {
  row = 2
  heights = c(2, 2)
} else if (nb <= 12 & nb > 8) {
  row = 3
  heights = c(3, 3)
} else if (nb <= 16 & nb > 12) {
  row = 4
  heights = c(4, 4)
} else if (nb <= 20 & nb > 16) {
  row = 5
  heights = c(5, 5)
} else if (nb <= 24 & nb > 20) {
  row = 6
  heights = c(6, 6)
} else if (nb <= 28 & nb > 24) {
  row = 7
  heights = c(7, 7)
} else {
  print("nb of snp too elevated, change -log(p_val) threshold")
}

# Call the plot function with parameters and files of interest
plot_regions(region_table, outlier_pick, region_label_tibbles, figure_path, pc)

# Create a summary table with most associated snps with their respective weights
save_snp_table(region_table, data_path, pc)



