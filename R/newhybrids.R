#!/usr/bin/env Rscript
# by: Floriane Coulmance: 02/10/2020
# usage:
# Rscript newhybrids.R <folder_path> 
# -------------------------------------------------------------------------------------------------------------------
# folder_path : $BASE_DIR/outputs/8_hybrids/\${POP_UN}_\${POP_DEUX}_nhinput/
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed library
library(parallelnewhybrid)


# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6]
print(args)

folder_path <- as.character(args[1]) # Path to mean coverage data table


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Run NEWHYBRID software, specify where to find the software and the necessary arguments
parallelnh_LINUX(folder.data = folder_path,
                 where.NH = "/gss/work/doau0129/software/newhybrids/",
                 burnin = 1000000,
                 sweeps = 10000000)
