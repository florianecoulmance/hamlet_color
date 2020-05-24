rm(list = ls())


library(gdsfmt)
library(SNPRelate)
library(stringi)
library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)
library(GWASTools)

setwd('/Users/fco/Desktop/BREMEN_OP/chapter1_2/GxP/')
df_gxp_snout <- read.table(file = 'snout.lm.50k.5k.txt.gz', sep = '\t', header = TRUE)
df_gxp_bars <- read.table(file = 'bars_body.lm.50k.5k.txt.gz', sep = '\t', header = TRUE)

manhattanPlot(df_gxp_snout$AVG_p_score, df_gxp_snout$CHROM, 
              ylim = NULL, trunc.lines = TRUE,
              signif = 5e-8, thinThreshold=NULL, pointsPerBin=10000, col=NULL)

manhattanPlot(df_gxp_bars$AVG_p_lrt, df_gxp_bars$CHROM, 
              ylim = NULL, trunc.lines = TRUE,
              signif = 5e-8, thinThreshold=NULL, pointsPerBin=10000, col=NULL)


