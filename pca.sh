#!/bin/bash

#SBATCH --job-name=pca
#SBATCH --partition=carl.p
#SBATCH --output=/user/doau0129/work/chapter1_2/logs/pca_%A_%a.out
#SBATCH --error=/user/doau0129/work/chapter1_2/logs/pca_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-02:00:00

module load hpc-env/6.4
module load R/3.6.1-intel-2018a

Rscript R/test.R
