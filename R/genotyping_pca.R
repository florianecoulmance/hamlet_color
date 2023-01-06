# by: Floriane Coulmance: 01/05/2020
# usage:
# Rscript pca.R <vcf_file> <output_path> <out_prefix>
# -------------------------------------------------------------------------------------------------------------------
# vcf_file in : 1) $BASE_DIR/ 2) 3) 4) 
# output_path in : $BASE_DIR/pca/
# out_prefix in : 1) 2) 3) 4)
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed library
library(gdsfmt)
library(SNPRelate)
library(tidyverse)
library(stringr)
library(stringi)
library(ggplot2)
library(ggpubr)
library(ggtext)
library(hypoimg)


# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6:9]
print(args)

vcf_file <- as.character(args[1]) # Path to vcf file
print(vcf_file)
output_path <- as.character(args[2]) # Path to the output folder
print(output_path)
out_prefix <- as.character(args[3]) # String for naming outputs
print(out_prefix)
figure_path <- as.character(args[4]) # Path to the figure folder
print(figure_path)

# output_path <- "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/"
# out_prefix <- "snp_filterd"
# figure_path <- "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/figures/genotyping_pca/"

# -------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


genotyping_pca_files <- function(vcf,path,prefix) {

  # Takes a vcf file with genotype SNPs and create pca and associated files

  ld_threshold <- 0.2 # set threshold of linkage desiquilibrium
  maf_threshold <- NaN

  # Read the files and convert to pca format
  gds_file <- str_c(path,prefix,ld_threshold,maf_threshold,".gds")
  snpgdsVCF2GDS(vcf.fn=vcf, out.fn=gds_file, method="biallelic.only")
  genofile <- snpgdsOpen(gds_file)
  snpset <- snpgdsLDpruning(genofile, ld.threshold = ld_threshold, maf=maf_threshold, method="corr", autosome.only = FALSE)
  snpset.id <- unlist(snpset)

  # perform pca
  pca <- snpgdsPCA(genofile, snp.id = snpset.id, num.thread = 2, autosome.only = FALSE)

  # get eigenvectors
  pca_tib <- pca$eigenvect %>%
             as_tibble() %>% 
             set_names(nm = str_pad(1:length(names(.)),2,pad = '0') %>%
             str_c("PC",.)) %>%
             mutate(sample = pca$sample.id,
                    spec = str_sub(sample,-6,-4),
                    loc = str_sub(sample,-3,-1),
                    pop = str_c(spec,loc))

  # get variance table
  explained_var <- tibble(variation = pca$varprop*100) %>%
  mutate(PC = str_c("PC", str_pad(row_number(),
                                  width = 2,
                                  pad = "0"))) %>%
  select(PC, variation) %>%
  filter(!is.na(variation))

  # Save important files for plots
  save(pca,file = str_c(path,prefix,"_pca.RData"))
  pca_tib %>% write_tsv(path = str_c(path,prefix,"_eigenvectors.tsv"))
  explained_var %>% write_tsv(path = str_c(path,prefix,"_eixplained_var.tsv"))

  # Close files
  snpgdsClose(genofile)
  system(str_c("rm ", gds_file))
  showfile.gds(closeall=TRUE)

}


genotyping_pca_plots <- function(path, prefix, pathfigure) {

  # create pca plots from pca files previously created
  
  # load and read the necessary files
  load(str_c(path,prefix,"_pca.RData"))
  variance <- read_tsv(str_c(path,prefix,"_eixplained_var.tsv"))
  names(variance) <- c("index", "variation") 
  print(variance)

  # modifications in pca files
  data <- data.frame(pca[["eigenvect"]])
  names(data) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15")

  id <- pca[["sample.id"]]
  data["id"] <- id
  print(data)
#   data$id <- gsub("PL17_108nigbel", "PL17_108indbel", data$id)
#   data$id <- gsub("PL17_111indbel", "PL17_111nigbel", data$id)
  print(data[26:36,])
  data["geo"] <- stri_sub(data$id,-3,-1)
  data["spec"] <- stri_sub(data$id,-6,-4)

  # create a list of species logo to integrate to plots
  logos_spec <- c("abe" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_aberrans.l.cairo.png' width='120' /><br>*H. aberrans*",
    "chl" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_chlorurus.l.cairo.png' width='120' /><br>*H. chlorurus*",
    "flo" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_floridae.l.cairo.png' width='120' /><br>*H. floridae*",
    "gem" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_gemma.l.cairo.png' width='120' /><br>*H. gemma*",
    "gum" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_gumigutta.l.cairo.png' width='120' /><br>*H. gummiguta*",
    "gut" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_guttavarius.l.cairo.png' width='120' /><br>*H. guttavarius*",
    "ind" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_indigo.l.cairo.png' width='120' /><br>*H. indigo*",
    "may" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_maya.l.cairo.png' width='120' /><br>*H. maya*",
    "nig" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_nigricans.l.cairo.png' width='120' /><br>*H. nigricans*",
    "pue" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_puella.l.cairo.png' width='120' /><br>*H. puella*",
    "ran" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_randallorum.l.cairo.png' width='120' /><br>*H. randallorum*",
    "tan" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_sp.l.cairo.png' width='120' /><br>*Hypoplectrus* sp.",
    "uni" = "<img src='/user/doau0129/work/chapter2/ressources/logos/H_unicolor.l.cairo.png' width='120' /><br>*H. unicolor*")

  # create a list of location logo to integrate to plots
  logos_loc <- c(bel = "*Belize*",  #<img src='/user/doau0129/work/chapter2/ressources/logos/belize.png' width='100' /><br>
                boc = "*Panama*", #<img src='/user/doau0129/work/chapter2/ressources/logos/pan.png' width='100' /><br>
                flo = "*Florida*", #<img src='/user/doau0129/work/chapter2/ressources/logos/us.png' width='100' /><br>
                pue = "*Puerto Rico*") #<img src='/user/doau0129/work/chapter2/ressources/logos/puer.png' width='100' /><br>

  # create list of colors to use in plots 
  spec_colors <- c("nig" = '#FF0033', "chl" = '#9900CC', "abe" = '#996600', "gut" = '#0000FF',
                   "gum" = '#FF00FF', "ran" = '#666699', "gem" = '#CC0000', "may" = '#FF9933',
                   "ind" = '#66CCFF', "pue" = '#FFCC00', "flo" = '#33FFCC', "tan" = '#333333',
                   "uni" = '#66CC00')

  # PC1 vs. PC2
  p1 <- ggplot(data,aes(x=PC1,y=PC2,color=spec)) + geom_point(size = 7, aes(shape=geo)) + 
        labs(x = paste0("PC1, var =  ", format(round(variance$variation[1], 1), nsmall = 1), " %") , y = paste0("PC2, var = ", format(round(variance$variation[2], 1), nsmall = 1), " %"))
  p1 <- p1 + scale_color_manual(values=spec_colors, labels = logos_spec) +
        scale_shape_manual(values=c(16,17,15,8), labels = logos_loc) +
        scale_x_continuous(position = "top") +
        theme(legend.position="bottom",legend.title=element_blank(),legend.box = "vertical",legend.text =  element_markdown(size = 30),
                  panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.9),
                  text = element_text(size=30), legend.key=element_blank()) +
        guides(color = guide_legend(nrow = 1))

  # PC1 vs PC3
  p2 <- ggplot(data,aes(x=PC1,y=PC3,color=spec)) + geom_point(size = 7, aes(shape=geo)) +
        labs(x = paste0("PC1, var =  ", format(round(variance$variation[1], 1), nsmall = 1), " %") , y = paste0("PC3, var = ", format(round(variance$variation[3], 1), nsmall = 1), " %"))
  p2 <- p2 + scale_color_manual(values=spec_colors, labels = logos_spec) +
        scale_shape_manual(values=c(16,17,15,8), labels = logos_loc) +
        scale_x_continuous(position = "top") +
        theme(legend.position="bottom",legend.title=element_blank(),legend.box = "vertical",legend.text =  element_markdown(size = 30),
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),
          text = element_text(size=30), legend.key=element_blank()) +
        guides(color = guide_legend(nrow = 1))

  # PC2 vs PC3
  p3 <- ggplot(data,aes(x=PC2,y=PC3,color=spec)) + geom_point(size = 7, aes(shape=geo)) +
        labs(x = paste0("PC2, var =  ", format(round(variance$variation[2], 1), nsmall = 1), " %") , y = paste0("PC3, var = ", format(round(variance$variation[3], 1), nsmall = 1), " %"))
  p3 <- p3 + scale_color_manual(values=spec_colors, labels = logos_spec) +
        scale_shape_manual(values=c(16,17,15,8), labels = logos_loc) +
        scale_x_continuous(position = "top") +
        theme(legend.position="bottom",legend.title=element_blank(),legend.box = "vertical",legend.text =  element_markdown(size = 30),
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),
          text = element_text(size=30), legend.key=element_blank()) +
        guides(color = guide_legend(nrow = 1))
  
  # PC1 vs PC2 zoom 1
  p4 <- ggplot(data,aes(x=PC1,y=PC2,color=spec)) + geom_point(size = 7, aes(shape=geo)) + 
        labs(x = paste0("PC1, var =  ", format(round(variance$variation[1], 1), nsmall = 1), " %") , y = paste0("PC2, var = ", format(round(variance$variation[2], 1), nsmall = 1), " %"))
  p4 <- p4 + xlim(0, 0.22) + ylim(-0.03,-0.01) 
  p4 <- p4 + scale_color_manual(values=spec_colors, labels = logos_spec) +
        scale_shape_manual(values=c(16,17,15,8), labels = logos_loc) +
        scale_x_continuous(position = "top") +
        theme(legend.position="bottom",legend.title=element_blank(),legend.box = "vertical",legend.text =  element_markdown(size = 30),
                  panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.9),
                  text = element_text(size=30), legend.key=element_blank()) +
        guides(color = guide_legend(nrow = 1))

  # PC1 vs PC2 zoom 2
  p5 <- ggplot(data,aes(x=PC1,y=PC2,color=spec)) + geom_point(size = 7, aes(shape=geo)) + 
        labs(x = paste0("PC1, var =  ", format(round(variance$variation[1], 1), nsmall = 1), " %") , y = paste0("PC2, var = ", format(round(variance$variation[2], 1), nsmall = 1), " %"))
  p5 <- p5 + xlim(-0.08, -0.03) + ylim(-0.01,0.015) 
  p5 <- p5 + scale_color_manual(values=spec_colors, labels = logos_spec) +
        scale_shape_manual(values=c(16,17,15,8), labels = logos_loc) +
        scale_x_continuous(position = "top") +
        theme(legend.position="bottom",legend.title=element_blank(),legend.box = "vertical",legend.text =  element_markdown(size = 30),
                  panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.9),
                  text = element_text(size=30), legend.key=element_blank()) +
        guides(color = guide_legend(nrow = 1))

  # PC1 vs PC2 zoom 3
  p6 <- ggplot(data,aes(x=PC1,y=PC2,color=spec, label = id)) + geom_point(size = 7, aes(shape=geo)) + 
        geom_text(hjust=0, vjust=-1, size=1) +
        labs(x = paste0("PC1, var =  ", format(round(variance$variation[1], 1), nsmall = 1), " %") , y = paste0("PC2, var = ", format(round(variance$variation[2], 1), nsmall = 1), " %"))
  p6 <- p6 + scale_color_manual(values=spec_colors, labels = logos_spec) +
        scale_shape_manual(values=c(16,17,15,8), labels = logos_loc) +
        theme(legend.position="bottom",legend.title=element_blank(),legend.box = "vertical",legend.text =  element_markdown(size = 30),
                  panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.9),
                  text = element_text(size=30), legend.key=element_blank()) +
        guides(color = guide_legend(nrow = 1))
  p6 <- p6 + xlim(0.02, 0.04) + ylim(0.27,0.31) #+ scale_x_continuous(position = "top")

  # PC1 vs PC2 zoom 4
  p7 <- ggplot(data,aes(x=PC1,y=PC2,color=spec, label = id)) + geom_point(size = 7, aes(shape=geo)) +
        geom_text(hjust=0, vjust=-1, size=1) +
        labs(x = paste0("PC1, var =  ", format(round(variance$variation[1], 1), nsmall = 1), " %") , y = paste0("PC2, var = ", format(round(variance$variation[2], 1), nsmall = 1), " %"))
  p7 <- p7 + scale_color_manual(values=spec_colors, labels = logos_spec) +
        scale_shape_manual(values=c(16,17,15,8), labels = logos_loc) +
        theme(legend.position="bottom",legend.title=element_blank(),legend.box = "vertical",legend.text =  element_markdown(size = 30),
                  panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.9),
                  text = element_text(size=30), legend.key=element_blank()) +
        guides(color = guide_legend(nrow = 1))
  p7 <- p7 + xlim(-0.04, 0.015) + ylim(-0.08,0.025) #+ scale_x_continuous(position = "top")
  

  # PC1 vs PC2 zoom 5
  p8 <- ggplot(data,aes(x=PC1,y=PC2,color=spec, label = id)) + geom_point(size = 7, aes(shape=geo)) + 
        geom_text(hjust=0, vjust=-1, size=1) +
        labs(x = paste0("PC1, var =  ", format(round(variance$variation[1], 1), nsmall = 1), " %") , y = paste0("PC2, var = ", format(round(variance$variation[2], 1), nsmall = 1), " %"))
  p8 <- p8 + scale_color_manual(values=spec_colors, labels = logos_spec) +
        scale_shape_manual(values=c(16,17,15,8), labels = logos_loc) +
        theme(legend.position="bottom",legend.title=element_blank(),legend.box = "vertical",legend.text =  element_markdown(size = 30),
                  panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.9),
                  text = element_text(size=30), legend.key=element_blank()) +
        guides(color = guide_legend(nrow = 1))
  p8 <- p8 + xlim(-0.005, 0.015) + ylim(-0.08,0.025) #+ scale_x_continuous(position = "top")
  

  # Arrange all plots in pae
  f <- ggarrange(p1, p2, p3, common.legend = TRUE, legend="bottom", ncol = 3, nrow = 1, labels = c("a", "b", "c"), font.label = list(size = 30))
  g <- if(grepl("casz1",prefix)) ggarrange(p1, p4, p5, common.legend = TRUE, legend="bottom", ncol = 3, nrow = 1) else ggarrange(p1, p6, p7, p8, common.legend = TRUE, legend="bottom", ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"), font.label = list(size = 30))

  # Save figures
  hypo_save(filename = paste0(pathfigure,"ld0.2_nomaf_corr/",prefix,"_pca_swap_genofile.pdf"),
          plot = f,
          width = 40,
          height = 18)

  hypo_save(filename = paste0(pathfigure,"ld0.2_nomaf_corr/",prefix,"_Lpca_zoom_swap_genofile.pdf"),
            plot = g,
            width = 40,
            height = 30)

}


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------

# Perform PCA
genotyping_pca_files(vcf_file,output_path,out_prefix)

# Create pca plots
genotyping_pca_plots(output_path, out_prefix, figure_path)

