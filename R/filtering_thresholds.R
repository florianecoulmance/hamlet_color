#!/usr/bin/env Rscript
# by: Floriane Coulmance: 01/05/2020
# usage:
# Rscript --vanilla pca.R <metrics_table> <figures_path>
# -------------------------------------------------------------------------------------------------------------------
# metrics_table in : $BASE_DIR/outputs/6_genotyping/6_1_snp/raw_var_sites.table.txt
# figures_path in : $BASE_DIR/figures/
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed library
library(ggplot2)
library(gridExtra)


# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[7:8]
print(args)

metrics_table <- as.character(args[1]) # Path to metrics table
figures_path <- as.character(args[2]) # Path to the figure folder


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS
metrics <- read.table(metrics_table, sep='', header = TRUE)


mq <- ggplot(metrics, aes(x=MQ)) + 
  geom_density(color="darkblue", fill="lightblue") + xlim(50, 70) + geom_vline(xintercept=c(57.5,62.2), linetype="dotted")
#mq

#m <- metrics[which.max(metrics$MQ),]
#print(m)

qd <- ggplot(metrics, aes(x=QD)) + 
  geom_density(color="darkblue", fill="lightblue") + xlim(0, 50) + geom_vline(xintercept=c(4), linetype="dotted")

#qd

logFS <- log(metrics$FS)
fs <- ggplot(metrics, aes(x=logFS)) + 
  geom_density(color="darkblue", fill="lightblue") + xlim(-5, 15) + geom_vline(xintercept=c(1.77), linetype="dotted")
#fs

mqrs <- ggplot(metrics, aes(x=MQRankSum)) + 
  geom_density(color="darkblue", fill="lightblue") + xlim(-1, +1) + geom_vline(xintercept=c(-0.2,0.2), linetype="dotted")
#mqrs

rprs <- ggplot(metrics, aes(x=ReadPosRankSum)) + 
  geom_density(color="darkblue", fill="lightblue") + xlim(-5, +5) + geom_vline(xintercept=c(-2,2), linetype="dotted")
#rprs

figure <- grid.arrange(mq, qd, fs, mqrs, rprs, ncol = 1, nrow = 5)

ggsave(file=paste0(figures_path,"filtering_thresholds.pdf"), figure)
