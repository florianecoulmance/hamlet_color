#!/bin/bash

#SBATCH --job-name=01_mark_adapters_ubams
#SBATCH --partition=carl.p
#SBATCH --array=1-117
#SBATCH --output=/user/doau0129/work/chapter1_2/logs/01_mark_adapters_ubams_%A_%a.out
#SBATCH --error=/user/doau0129/work/chapter1_2/logs/01_mark_adapters_ubams_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=04:00:00

BASE_DIR=/Users/fco/Desktop/BREMEN_OP/chapter1_2/


mkdir $BASE_DIR/outputs/01_adapters/
mkdir $BASE_DIR/outputs/01_adapters/adapters/
mkdir $BASE_DIR/outputs/01_adapters/metrics/
mkdir $BASE_DIR/outputs/listoffiles/

ls -1 $BASE_DIR/outputs/00_ubams/ > $BASE_DIR/outputs/listoffiles/ubams.fofn

INPUT_UBAMS=$BASE_DIR/outputs/listoffiles/ubams.fofn

UBAMS=$(cat $INPUT_BAMS | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)

sample=${UBAMS%.*}
echo $sample

# gatk --java-options "-Xmx18G" \
#   MarkIlluminaAdapters \
#   -I=$BASE_DIR/outputs/00_ubams/${UBAMS} \
#   -O=$BASE_DIR/outputs/01_adapters/adapters/${sample}.adapter.bam \
#   -M=$BASE_DIR/outputs/01_adapters/metrics/${sample}.adapter.metrics.txt \
#   -TMP_DIR=\$BASE_DIR/temp_files
