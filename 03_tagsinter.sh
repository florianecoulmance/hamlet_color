#!/bin/bash

#SBATCH --job-name=03_2_tagsinter
#SBATCH --partition=carl.p
#SBATCH --array=1-117
#SBATCH --output=/user/doau0129/work/chapter1_2/logs/03_2_tagsinter_%A_%a.out
#SBATCH --error=/user/doau0129/work/chapter1_2/logs/03_2_tagsinter_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=2-00:00:00

BASE_DIR=/user/doau0129/work/chapter1_2/


INPUT_SORTSAM=$BASE_DIR/outputs/listoffiles/sort_sam.fofn

SORTSAM=$(cat ${INPUT_SORTSAM} | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo $SORTSAM

sample1=${SORTSAM%.*}
sample=${sample1%.*}
echo $sample

gatk --java-options "-Xmx30G" \
  SetNmAndUqTags \
  --INPUT=$BASE_DIR/outputs/03_mark_duplicates/sort_sam/${SORTSAM} \
  --OUTPUT=$BASE_DIR/outputs/03_mark_duplicates/tags_intermediate/${sample}.intermediate.bam \
  --CREATE_INDEX=true \
  --CREATE_MD5_FILE=true \
  -TMP_DIR=$BASE_DIR/temp_files \
  --REFERENCE_SEQUENCE=$BASE_DIR/ressources/HP_genome_unmasked_01.fa.gz
