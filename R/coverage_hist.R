#!/usr/bin/env Rscript
# by: Floriane Coulmance: 02/10/2020
# usage:
# Rscript --vanilla coverage_hist.R <data_path> <figure_path> <list_path> 
# -------------------------------------------------------------------------------------------------------------------
# data_path in : $BASE_DIR/outputs/coverage/coverage_table
# figure_path in : $BASE_DIR/figures/
# list_path in : $BASE_DIR/outputs/lof/
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed library
library(ggplot2)
library(tidyverse)


# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[7:9]
print(args)

data_path <- as.character(args[1]) # Path to mean coverage data table
# data_path <- "/Users/fco/Desktop/BREMEN_OP/ibd/coverage_table"
figure_path <- as.character(args[2]) # Path to the figure folder
# figure_path <- "/Users/fco/Desktop/BREMEN_OP/ibd/figures/"
list_path <- as.character(args[3]) # Path to the folder that will contain the list of genetic files to remove
# list_path <- "/Users/fco/Desktop/BREMEN_OP/ibd/"


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


coverage <- read.table(data_path,sep = "")
mean <- mean(coverage$V2)
lab <- round(mean,digits=2)
lab <- paste0("mean coverage : ",lab,"x")


p <- ggplot(coverage, aes(x=V2, y=reorder(V1, V2))) + geom_bar(stat = "identity", fill = "#DDB9F5") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        axis.line.y = element_blank()) + 
  scale_x_continuous(expand = c(0,0)) +
  geom_vline(xintercept= mean(coverage$V2), size=0.3, label=mean) +
  geom_vline(xintercept=15, linetype="dashed", color = "red", size=0.3) +
  geom_text(aes(x=mean(coverage$V2), label=lab, y=20), size=1.5, hjust=-0, vjust=28) + 
  geom_text(aes(x=15, label="15x", y=20), color = "red", size=1.5, hjust=-0.3, vjust=28) +
  labs(x = "mean coverage per bases")

pdf(paste0(figure_path,"coverage_histogram.pdf"))
print(p)
dev.off()

to_remove_table <- coverage %>%
  filter(V2 < 15.00)

fileConn <- file(paste0(list_path,"b_remove.fofn"))
writeLines(to_remove_table$V1, fileConn)
close(fileConn)
        

  
