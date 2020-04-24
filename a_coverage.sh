#!/bin/bash

#SBATCH --job-name=a_coverage
#SBATCH --partition=carl.p
#SBATCH --array=1-117
#SBATCH --output=/user/doau0129/work/chapter1_2/logs/a_coverage_%A_%a.out
#SBATCH --error=/user/doau0129/work/chapter1_2/logs/a_coverage_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=02:00:00

BASE_DIR=/user/doau0129/work/chapter1_2/


INPUT_COV=$BASE_DIR/outputs/listoffiles/duplicates_coverage.fofn

COV=$(cat ${INPUT_COV} | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo $COV


sample1=${COV%.*}
sample=${sample1%.*}
echo $sample

gatk --java-options "-Xmx30G" \
  CollectWgsMetrics \
  -I=$BASE_DIR/outputs/03_mark_duplicates/mark_duplicates/duplicates/${COV} \
  -O=$BASE_DIR/outputs/a_coverage/${sample}.wgsmetrics.txt \
  -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa.gz
