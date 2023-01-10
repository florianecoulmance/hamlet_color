# by: Floriane Coulmance: 20/10/2022
# usage:
# Rscript snp_alleles_plots.R <allele_file> <figure_path>
# --------------------------------------------------------------------------------------------------------------------
# allele_file in : $BASE_DIR/outputs/7_gxp/LAB_fullm_54off_59on/PC1_5/PC1_5_alleles.txt
# figure_path in : $BASE_DIR/figures/7_gxp/continuous/LAB/LAB_fullm_54off_59on/PC1_5/
# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed library
library(broom, lib.loc=.libPaths()[-1])
library(dbplyr, lib.loc=.libPaths()[-1])
library(GenomicOriginsScripts)
library(ggplot2)
library(tidyverse)
library(hypogen)
library(hypoimg)
library(reshape2)
library(scatterpie)
library(ggtext)

# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6:7]
print(args)

allele_file <- as.character(args[1]) # Path to list of alleles file
figure_path <- as.character(args[2]) # Path to figure folder
# allele_file <- "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/PC1_5_alleles.txt"
# figure_path <- "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/figures/7_gxp/continuous/LAB/LAB_fullm_54off_59on/PC1_5/"


# -------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


pie_plots <- function(dat) {
  
  # Take a dataframe and makes a pie chart 
  
  p <- ggplot(data = dat) +
    geom_scatterpie(aes(x=x_pos, y=y_pos, group=allele, r=radius),
                    data=dat, cols=colnames(dat)[2:14], color=NA, alpha=.8) +
    geom_text(aes(y = y_pos + 3, x = x_pos, 
                  label = paste0(label,"\n",allele,"\nnb individuals : ",sum))) +
    scale_fill_manual(values=c("abe"="#E5E5A1",
                               "aff"="#FFB1D8",
                               "chl"="#8B4513",
                               "flo"="#A39CCF",
                               "gem"="#1E90FF",
                               "gum"="#FAA342",
                               "gut"="#FFEA00",
                               "ind"="#0A0082",
                               "may"="#7EA7C2",
                               "nig"="#1E1E1E",
                               "pue"="#E17366",
                               "ran"="#7EBDB3",
                               "uni"="#CCCCCC",
                               "tan"="#D2B48C"),
                      labels = c("abe" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_aberrans.l.cairo.png' width='40' /><br>*H. aberrans*",
                                 "aff" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_tan.l.cairo.png' width='40' /><br>*H. affinis*",
                                 "chl" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_chlorurus.l.cairo.png' width='40' /><br>*H. chlorurus*",
                                 "flo" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_floridae.l.cairo.png' width='40' /><br>*H. floridae*",
                                 "gem" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_gemma.l.cairo.png' width='40' /><br>*H. gemma*",
                                 "gum" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_gumigutta.l.cairo.png' width='40' /><br>*H. gummigutta*",
                                 "gut" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_guttavarius.l.cairo.png' width='40' /><br>*H. guttavarius*",
                                 "ind" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_indigo.l.cairo.png' width='40' /><br>*H. indigo*",
                                 "may" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_maya.l.cairo.png' width='40' /><br>*H. maya*",
                                 "nig" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_nigricans.l.cairo.png' width='40' /><br>*H. nigricans*",
                                 "pue" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_puella.l.cairo.png' width='40' /><br>*H. puella*",
                                 "ran" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_randallorum.l.cairo.png' width='40' /><br>*H. randallorum*",
                                 "uni" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_unicolor.l.cairo.png' width='40' /><br>*H. unicolor*",
                                 "tan" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_tan.l.cairo.png' width='40' /><br>*H. affinis*")) + 
    
    theme(legend.title=element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.grid = element_blank(),
          legend.position="none",
          legend.text =  element_markdown(size = 10),
          legend.key=element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.border = element_blank()) +
    guides(fill = guide_legend(nrow = 1))

  
}


save_plots <- function(plot, label, path_fig) {
  
  # Takes a plot and save it with label name and the figure path
  
  hypo_save(filename = paste0(path_fig,label,"_pie.png"),
            type = "cairo",
            plot = plot,
            width = 12,
            height = 8.5)
  
}


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------

allele_df <- read_tsv(allele_file) %>% 
             mutate(ID = paste0(CHROM,"_",POS))

allele_t <- dcast(melt(allele_df, id.vars = "ID"), variable ~ ID) %>% 
            slice_tail(n = 113)

names <- colnames(allele_t[,2:20])
test <- allele_t %>%
  mutate_at(vars(names), ~ gsub(":.*","",.)) %>%
  mutate_at(vars(names), ~ gsub("\\|","/",.)) %>%
  mutate(geo = str_sub(variable, -3, -1),
         spec = str_sub(variable, -6, -4))


for (i in names(test[,2:20])) {
  print(i)
  df <- test[c("variable", "spec", "geo",i)] %>% 
        setNames(., nm = c("SampleID", "spec", "geo", "allele"))
  print(df)
  pie_dat <- dcast(df,allele~spec, fill = 0, value.var="spec", fun.aggregate = length) %>%
             mutate(label = case_when(endsWith(allele, "0/0") ~ "homozygous reference",
                                      endsWith(allele, "0/1") ~ "heterozygous",
                                      endsWith(allele, "1/1") ~ "homozygous alternative"),
                    x_pos = case_when(endsWith(allele, "0/0") ~ 0,
                                      endsWith(allele, "0/1") ~ 0,
                                      endsWith(allele, "1/1") ~ 0),
                    y_pos = case_when(endsWith(allele, "0/0") ~ 28,
                                      endsWith(allele, "0/1") ~ 12,
                                      endsWith(allele, "1/1") ~ 4),
                    sum = as.numeric(rowSums(.[,2:14])),
                    radius = log(sum)/1.8)
  print(pie_dat)
  assign(i, pie_dat)
  
  p <- pie_plots(pie_dat)
  save_plots(p,i,figure_path)
  
}


