#!/bin/bash

#SBATCH --job-name=01_mark_adapters
#SBATCH --partition=carl.p
#SBATCH --array=1-117
#SBATCH --output=/user/doau0129/work/chapter1_2/logs/01_mark_adapters_%A_%a.out
#SBATCH --error=/user/doau0129/work/chapter1_2/logs/01_mark_adapters_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=04:00:00

BASE_DIR=/user/doau0129/work/chapter1_2/


mkdir $BASE_DIR/outputs/01_adapters/
mkdir $BASE_DIR/outputs/01_adapters/adapters/
mkdir $BASE_DIR/outputs/01_adapters/metrics/
mkdir $BASE_DIR/outputs/listoffiles/

INPUT_UBAMS=$BASE_DIR/outputs/listoffiles/ubams.fofn

UBAMS=$(cat ${INPUT_UBAMS} | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo $UBAMS

sample1=${UBAMS%.*}
sample=${sample1%.*}
echo $sample

gatk --java-options "-Xmx18G" \
   MarkIlluminaAdapters \
   -I=$BASE_DIR/outputs/00_ubams/${UBAMS} \
   -O=$BASE_DIR/outputs/01_adapters/adapters/${sample}.adapter.bam \
   -M=$BASE_DIR/outputs/01_adapters/metrics/${sample}.adapter.metrics.txt \
   -TMP_DIR=$BASE_DIR/temp_files
