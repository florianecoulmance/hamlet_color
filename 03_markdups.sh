#!/bin/bash

#SBATCH --job-name=03_3_markdups
#SBATCH --partition=carl.p
#SBATCH --array=1-4
#SBATCH --output=/user/doau0129/work/chapter1_2/logs/03_3_markdups_%A_%a.out
#SBATCH --error=/user/doau0129/work/chapter1_2/logs/03_3_markdups_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=05:00:00

BASE_DIR=/user/doau0129/work/chapter1_2/


INPUT_TAGSINTER=$BASE_DIR/outputs/listoffiles/tags_intermediate3.fofn

TAGSINTER=$(cat ${INPUT_TAGSINTER} | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo $TAGSINTER

sample1=${TAGSINTER%.*}
sample=${sample1%.*}
echo $sample

gatk --java-options "-Xmx30G" \
  MarkDuplicates \
  -I=$BASE_DIR/outputs/03_mark_duplicates/tags_intermediate/${TAGSINTER} \
  -O=$BASE_DIR/outputs/03_mark_duplicates/mark_duplicates/duplicates/${sample}.dedup.bam \
  -M=$BASE_DIR/outputs/03_mark_duplicates/mark_duplicates/metrics/${sample}.dedup.metrics.txt \
  -MAX_FILE_HANDLES=1000  \
  -TMP_DIR=$BASE_DIR/temp_files
