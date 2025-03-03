# by: Floriane Coulmance: 15/04/2021
# usage:
# Rscript gxp_zooms.R <data_path> <data_file> <figure_path> <pc> <chromosome>
# -------------------------------------------------------------------------------------------------------------------
# data_path in : $BASE_DIR/outputs/7_gxp/$DATASET
# data_file in : *mvplink.50k.5k.txt.gz
# figure_path in : $BASE_DIR//figures/7_gxp/$DATASET/\${NAME}
# pc in : 55 names of univariate and multivariate PCs association analysis
# chromosome in one of the 24 hamlets Linkage Group
# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed library
# library(broom, lib.loc=.libPaths()[-1])
# library(dbplyr, lib.loc=.libPaths()[-1])
library(broom)
library(dbplyr)
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
library(png)


# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6:10]
print(args)

data_path <- as.character(args[1]) # Path to data folder with GWAS files
data_file <- as.character(args[2]) # Name of the GWAS 50kb windowed-averaged file
figure_path <- as.character(args[3]) # Path to the figure folder
pc <- as.character(args[4]) # Name of the multivariate GWAS
chrom <- as.character(args[5]) # Name of the considered chromosome
# data_path <- "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/"
# data_file <- "PC1_5.mvplink.50k.5k.txt.gz"
# figure_path <- "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/figures/7_gxp/continuous/LAB/LAB_fullm_54off_59on/PC1_5/"
# pc <- "PC1_5"
# chrom <- "LG12"


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


custom_annoplot_flo <- function (..., searchLG, xrange, genes_of_interest = c(), genes_of_sec_interest = c(),
                                 anno_rown = 3, width = 0.1, gene_color = 'darkgray', start, end) {
  
  # Extract annotation information of specific LG range and plot the specific region
  
  # Import annotation data for a specific LG (chromosome) and range
  df_list <- hypogen::hypo_annotation_get(searchLG = searchLG, xrange = xrange,
                                          genes_of_interest = genes_of_interest,
                                          genes_of_sec_interest = genes_of_sec_interest,
                                          anno_rown = anno_rown)
  # Plot the annotation
  ggplot2::ggplot() +
    # add exon 
    ggplot2::geom_rect(data = df_list[[2]],
                       aes(xmin = ps, xmax = pe, ymin = yl - (width/2),
                           ymax = yl + (width/2), group = Parent),
                       fill = alpha(gene_color,.6), col = gene_color, lwd = 0.9) +
    # add outlier area
    ggplot2::geom_rect(inherit.aes = FALSE,
                       data = tibble::tibble(start = start, end = end),
                       aes(xmin = start, xmax = end),
                       ymin = -Inf, ymax = Inf,
                       fill=rgb(1,1,1,.3),color = rgb(1,1,1,.9)) +
    # add gene direction if known
    ggplot2::geom_segment(data = (df_list[[1]] %>% dplyr::filter(strand %in% c("+", "-"))),
                          aes(x = ps, xend = pe, y = yl, yend = yl, group = Parent),
                          lwd = 0.9, arrow = arrow(length = unit(1.5,"pt"), type = "closed"),
                          color = clr_genes) +
    # add gene extent if direction is unknowns
    ggplot2::geom_segment(data = (df_list[[1]] %>% dplyr::filter(!strand %in% c("+", "-"))),
                          aes(x = ps, xend = pe, y = yl, yend = yl, group = Parent),
                          lwd = 0.9, color = clr_genes) +
    # add gene label
    ggplot2::geom_text(data = df_list[[1]],
                       aes(x = labelx, label = gsub("hpv1g000000", ".", label), y = yl - 0.5), size = 7)
  
}


plot_panel_anno_flo <- function(outlier_id, label, lg, start, end, genes = c(),...)  {
  
  # Plot an annotation for outlier region of interest
  
  # Create plot title
  ttle <- stringr::str_sub(outlier_id,1,4) #%>% stringr::str_c(.,' (',project_inv_case(label),')')
  
  # Create plot of region annotation
  p <- custom_annoplot_flo(searchLG = lg,
                           xrange = c(start-window_buffer*1.25,end+window_buffer*1.25),
                           genes_of_interest = genes,
                           anno_rown = 6, start = start, end = end) +
    # layout x ayis
    ggplot2::scale_x_continuous(name = ttle,
                                position = 'top',
                                expand = c(0,0),
                                limits = c(start-window_buffer*1.25, end+window_buffer*1.25),
                                labels = ax_scl) +
    # layout y ayis
    ggplot2::scale_y_continuous(name = expression(bolditalic(Genes)), expand = c(0,.4)) +
    # use same boundaries for all panels
    ggplot2::coord_cartesian(xlim = c(start-window_buffer, end+window_buffer)) +
    # special panel layout for annotation panel
    hypogen::theme_hypo() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = rgb(.9,.9,.9),
                                                            color = rgb(1,1,1,0)),
                   legend.position = 'none',
                   axis.title.y = ggplot2::element_text(size = 24,vjust = 6),
                   axis.title.x = ggplot2::element_text(size = 28),
                   axis.line = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.x = element_text(size = 24))
  
  # Use correct greek symbols in labels if needed
  if(outlier_id == 'LG17_1') {
    p$layers[[5]]$data$label <- p$layers[[5]]$data$label %>%
      stringr::str_replace(., 'alpha', "\u03B1") %>%
      stringr::str_replace(.,  'beta', "\u03B2")
  }
  
  p
  
}


# plot_panel_gxp <- function (lg, start, end, trait, ...) {
#   
#   # Create gxp plot for 50k windowed GWAS results
#   
#   ggplot2::ggplot() +
#   geom_rect(inherit.aes = FALSE, data = tibble(start = start,end = end),
#             aes(xmin = start, xmax = end),
#             ymin = -Inf, ymax = Inf, fill = rgb(0.9, 0.9, 0.9, 0.3),
#             color = rgb(0.9, 0.9, 0.9, 0.9)) + 
#   geom_line(data = gxp_data %>% filter(CHROM == lg, MID_POS > start - window_buffer * 1.25, MID_POS < end + window_buffer * 1.25),
#             aes(x = MID_POS, y = LOG_P, color = RUN),
#             size = 0.7) +
#   scale_color_manual(name = "GxP Trait", values = gxp_clr) +
#   scale_x_continuous(name = lg, expand = c(0, 0), position = "top") +
#   scale_y_continuous(name = expression(bolditalic(-log[10](p))), expand = c(0, 0)) +
#   guides(color = guide_legend(keyheight = unit(3,"pt"), keywidth = unit(20, "pt"),
#                               override.aes = list(size = 2))) +
#     theme_hypo() + hypogen::theme_hypo_anno_extra() + 
#     ggplot2::theme(legend.position = "none",
#                    axis.title.y = ggplot2::element_text(angle = 90, size = 26), 
#                    axis.line.y = ggplot2::element_line(size = plot_lwd), 
#                    axis.ticks.y = ggplot2::element_line(size = plot_lwd), 
#                    axis.title.x = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank(), 
#                    axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
#                    axis.text.y = ggplot2::element_text(size = 24),
#                    plot.margin = ggplot2::margin(t = 1, r = 1, b = 3, 
#                                                  l = 1), ...)
# }


plot_panel_gxp <- function (lg, start, end, trait, ...) {
  
  # Create gxp plot for 50k windowed GWAS results
  
  ggplot2::ggplot() +
    geom_rect(inherit.aes = FALSE, data = tibble(start = start,end = end),
              aes(xmin = start, xmax = end),
              ymin = -Inf, ymax = Inf, fill = rgb(0.9, 0.9, 0.9, 0.3),
              color = rgb(0.9, 0.9, 0.9, 0.9)) + 
    geom_point(data = gxp_snp %>% filter(CHROM == lg, MID_POS > start - window_buffer * 1.25, MID_POS < end + window_buffer * 1.25),
               aes(x = MID_POS, y = LOG_P, color = "black"),
               size = 1, stroke = 0.2) +
    geom_point(data = gxp_snp %>% filter(CHROM == lg, MID_POS > start + 1500 & MID_POS < end - 1500)
               %>% filter(row_number(desc(LOG_P)) <= 1),
               color ="red",
               aes(x=MID_POS, y=LOG_P, label = RANGE), size = 4) +
    geom_text_repel(data = gxp_snp %>% filter(CHROM == lg, MID_POS > start + 1500 & MID_POS < end - 1500)
                    %>% filter(row_number(desc(LOG_P)) <= 1),
                    aes(x = MID_POS, y = LOG_P, color = RUN, label = RANGE),
                    segment.color = 'red', min.segment.length = 0.1, segment.size = 0.05,
                    color = "Red", size=10, hjust = -1, nudge_x = -1.5) +
    geom_line(data = gxp_data %>% filter(CHROM == lg, MID_POS > start - window_buffer * 1.25, MID_POS < end + window_buffer * 1.25),
              aes(x = MID_POS, y = LOG_P/0.20, color = RUN),
              size = 0.7, color = gxp_clr) +
    # scale_color_manual(name = "GxP Trait", values = gxp_clr) +
    scale_x_continuous(name = lg, expand = c(0, 0), position = "top") +
    scale_y_continuous(name = expression(bolditalic(-log[10](p))), expand = c(0, 0), limits = c(0, 25), sec.axis = sec_axis(~ . * 0.20)) +
    guides(color = guide_legend(keyheight = unit(3,"pt"), keywidth = unit(20, "pt"),
                                override.aes = list(size = 2))) +
    scale_color_manual(name = "GxP Trait", values = c("black", "red", gxp_clr)) +
    # scale_x_continuous(name = lg, expand = c(0, 0), position = "bottom") +
    # scale_y_continuous(name = expression(bolditalic(-log[10](p))),expand = c(0, 0)) +
    # guides(color = guide_legend(keyheight = unit(3,"pt"), keywidth = unit(20, "pt"),
    #                             override.aes = list(size = 2))) +
    theme_hypo() + hypogen::theme_hypo_anno_extra() + 
    ggplot2::theme(legend.position = "none",
                   axis.title.y = ggplot2::element_text(angle = 90, size = 26), 
                   axis.line.y = ggplot2::element_line(size = plot_lwd), 
                   axis.ticks.y = ggplot2::element_line(size = plot_lwd), 
                   axis.title.x = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank(), 
                   axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 24),
                   axis.line.y.right = element_line(color = gxp_clr),
                   axis.text.y.right = element_text(color = gxp_clr),
                   plot.margin = ggplot2::margin(t = 1, r = 1, b = 3, 
                                                 l = 1), ...)
  
}


plot_panel_gxp_snp <- function (lg, start, end, trait, ...) {
  
  # Create gxp plot for snp-by-snp GWAS results
  
  ggplot2::ggplot() +
  geom_rect(inherit.aes = FALSE, data = tibble(start = start,end = end),
            aes(xmin = start, xmax = end),
            ymin = -Inf, ymax = Inf, fill = rgb(0.9, 0.9, 0.9, 0.3),
            color = rgb(0.9, 0.9, 0.9, 0.9)) + 
  geom_point(data = gxp_snp %>% filter(CHROM == lg, MID_POS > start - window_buffer * 1.25, MID_POS < end + window_buffer * 1.25),
             aes(x = MID_POS, y = LOG_P, color = gxp_clr),
             size = 1, stroke = 0.2) +
  geom_point(data = gxp_snp %>% filter(CHROM == lg, MID_POS > start + 1500 & MID_POS < end - 1500)
                            %>% filter(row_number(desc(LOG_P)) <= 1),
             color ="red",
             aes(x=MID_POS, y=LOG_P, label = RANGE), size = 4) +
  geom_text_repel(data = gxp_snp %>% filter(CHROM == lg, MID_POS > start + 1500 & MID_POS < end - 1500)
                                 %>% filter(row_number(desc(LOG_P)) <= 1),
                  aes(x = MID_POS, y = LOG_P, color = RUN, label = RANGE),
                  segment.color = 'red', min.segment.length = 0.1, segment.size = 0.05,
                  color = "Red", size=10, hjust = -1, nudge_x = -1.5) +
  scale_color_manual(name = "GxP Trait", values = c(gxp_clr, "red")) +
  scale_x_continuous(name = lg, expand = c(0, 0), position = "bottom") +
  scale_y_continuous(name = expression(bolditalic(-log[10](p))),expand = c(0, 0)) +
  guides(color = guide_legend(keyheight = unit(3,"pt"), keywidth = unit(20, "pt"),
                              override.aes = list(size = 2))) +
  theme_hypo() + hypogen::theme_hypo_anno_extra() + 
  ggplot2::theme(legend.position = "none",
                 axis.title.y = ggplot2::element_text(angle = 90, size = 26), 
                 axis.line.y = ggplot2::element_line(size = plot_lwd), 
                 axis.ticks.y = ggplot2::element_line(size = plot_lwd), 
                 axis.title.x = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank(), 
                 axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_text(size = 24),
                 plot.margin = ggplot2::margin(t = 1, r = 1, b = 3, 
                                                 l = 1), ...)
  
}


plot_curt <- function (outlier_id, outlier_nr, lg, start, end, text = TRUE, label, trait, heatmap, ...) {
  
  # Create zoom plots into GWAS peak regions for one particular LG, with all necessary 

  print(heatmap)

  # Annotation pannel
  p_g <- plot_panel_anno_flo(lg = lg, outlier_id = outlier_id, label = label,
                             start = start, end = end, genes = cool_genes)

  # Pannel for GWAS peak zoomed and 50kb averaged windows
  p_gxp <- plot_panel_gxp(lg = lg, start = start, end = end, trait = pcs)
  
  # Pannel for GWAS SNPs zoomed
  #p_snp <- plot_panel_gxp_snp(lg = lg, start = start, end = end, trait = pcs)
  
  #Pannel for heatmap
  img <- readPNG(heatmap)
  g <- rasterGrob(img[1:540,,], interpolate=TRUE, width = 1)
  g2 <- rasterGrob(img[665:820,,], interpolate=TRUE)
  g_plot <- ggplot() + annotation_custom(g) + theme(plot.background = element_blank(), 
                                                    panel.background = element_blank(),
                                                    panel.grid = element_blank(), 
                                                    panel.border = element_blank(),
                                                    plot.margin = unit(c(0, 0, 0, 0), "cm"))
  g2_plot <- ggplot() + annotation_custom(g2) + theme(plot.background = element_blank(), 
                                                      panel.background = element_blank(),
                                                      panel.grid = element_blank(), 
                                                      panel.border = element_blank(),
                                                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  # Assemble all panels
  if (text) {
    p_curtain <- cowplot::plot_grid(p_g, p_gxp, g_plot, g2_plot, ncol = 1, align = "v",
                                    rel_heights = c(1, rep(0.8,2), 0.3), axis="tblr")
  }
  else {
    p_curtain <- cowplot::plot_grid(p_g + no_title(), p_gxp + no_title(), g, g2,
                                    ncol = 1, align = "v", rel_heights = c(1, rep(0.8, 7)), axis="tblr")
  }
  

  hypo_save(filename = paste0(figure_path, pc, "_", outlier_id, ".png"),
            plot = p_curtain,
            width = 11.2,
            height = 16,
            type = 'cairo')

  print(p_curtain)
  # p_curtain

}


plot_regions <- function(region, outlier_list, label, path, name, chromosome) {
  
  # Create plot for regions with highest associations within each LG

  # Create the plot
  p_single <- region %>% filter(outlier_id %in% outlier_list) %>%
                         left_join(label) %>%
                         dplyr::mutate(outlier_nr = row_number()) %>% #, 
                                       #text = ifelse(outlier_nr == 1,TRUE,FALSE)) %>%
                         pmap(plot_curt, cool_genes = cool_genes) %>%
                         cowplot::plot_grid(plotlist = ., nrow = row, rel_heights = heights,
                                            labels = letters[1:length(outlier_list)] %>%
                                                     project_case())
  
  # Save the file in .png file
  hypo_save(filename = paste0(path, name, "_", chromosome, ".png"),
            plot = p_single,
            width = 32,
            # height = 18,
            height = 14.5,
            type = 'cairo')
  
  p_single
  
}


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Set general variables ---------------------------------------------------------------------------------------------

# Set parameters for plots
cool_genes <-  c('arl3','kif16b','cdx1','hmcn2', 'sox10','smarca4', 'rorb', 'alox12b','egr1', 
                 'ube4b','casz1', 'hoxc8a','hoxc9','hoxc10a', 'hoxc13a','rarga','rarg', 'snai1',
                 'fam83d','mafb','sws2abeta','sws2aalpha','sws2b','lws','grm8')

# set parameter of window buffer for plot
window_buffer <- 2.5*10^5

# set color parameter for plot line
gxp_clr <- c("#5B9E2E") %>% darken(factor = .95)


# Analysis ----------------------------------------------------------------------------------------------------------

# Read file and prepare useful columns
gxp_data <- prep_file(data_file,data_path)
data_snp <- paste0(pc,".mvplink.logarithm.txt.gz") 
gxp_snp <- prep_file(data_snp,data_path)


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
print(thresh)

# Create table with regions of interest
region_table <- threshold_table(thresh) %>% 
                mutate(heatmap = paste0(figure_path,gid,"_heatmaps.png")) %>%
                setNames(., nm = c("outlier_id", "lg", "start", "end", "gstart", "gend", "gpos", "heatmap"))

print(region_table["heatmap"][[1]])
print(region_table)

# List outlier regions' names
outlier_pick <- region_table$outlier_id
# print(outlier_pick)

nb <- length(outlier_pick) # count nb of regions 
# print(nb)

region_label_tibbles <- tibble(outlier_id = outlier_pick, label = letters[1:nb]) # prepare panels letters for plots

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

# Plot and save as table GWAS peaks for the LG
if (nb != 0) {
  plot_regions(region_table, outlier_pick, region_label_tibbles, figure_path, pc, chrom)
} else {
  print("this chromosome does not have peaks")
}



