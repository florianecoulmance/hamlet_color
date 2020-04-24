#!/bin/bash

#SBATCH --job-name=03_1_sortsam
#SBATCH --partition=carl.p
#SBATCH --array=1-117
#SBATCH --output=/user/doau0129/work/chapter1_2/logs/03_1_sortsam_%A_%a.out
#SBATCH --error=/user/doau0129/work/chapter1_2/logs/03_1_sortsam_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=2-00:00:00

BASE_DIR=/user/doau0129/work/chapter1_2/


INPUT_MAPPED_BAM=$BASE_DIR/outputs/listoffiles/mapped_bams.fofn

MAPPED_BAM=$(cat ${INPUT_MAPPED_BAM} | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo $MAPPED_BAM

sample1=${MAPPED_BAM%.*}
sample=${sample1%.*}
echo $sample

gatk --java-options "-Xmx30G" \
		SortSam \
		-I=$BASE_DIR/outputs/02_merged_map/bam_merge_alignement/${MAPPED_BAM} \
		-O=$BASE_DIR/outputs/03_mark_duplicates/sort_sam/${sample}.sorted.sam \
		--SORT_ORDER="coordinate" \
		--CREATE_INDEX=false \
		--CREATE_MD5_FILE=false \
		-TMP_DIR=$BASE_DIR/temp_files
