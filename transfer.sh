#!/bin/bash
#SBATCH --partition=carl.p
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90G
#SBATCH --time=4-00:00:00


rsync -r * /user/doau0129/data/chapter1/
#rsync -r figures/ /user/doau0129/data/chapter1/
#rsync -r sh/ /user/doau0129/data/chapter1/
#rsync -r R/ /user/doau0129/data/chapter1/
#rsync -r metadata/ /user/doau0129/data/chapter1/
#rsync -r images/ /user/doau0129/data/chapter1/
#rsync -r python/ /user/doau0129/data/chapter1/

