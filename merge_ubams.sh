#!/bin/bash

#SBATCH --job-name=00_readgroups_ubams
#SBATCH --partition=carl.p
#SBATCH --array=1-117
#SBATCH --output=/user/doau0129/work/chapter1_2/logs/%x_%A_%a.o.log
#SBATCH --error=/user/doau0129/work/chapter1_2/logs/%x_%A_%a.e.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=01:00:00

BASE_DIR=/user/doau0129/work/chapter1_2/


ls -1 $BASE_DIR/outputs/01_adapters/adapters/ > $BASE_DIR/outputs/listoffiles/adapters.bam.fofn

INPUT_ADAPT_BAMS=$BASE_DIR/outputs/listoffiles/adapters.bam.fofn

ADAPT_BAM=$(cat ${INPUT_ADAPT_BAMS} | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo $ADAPT_BAM

sample1=${ADAPT_BAM%.*}
sample=${sample1%.*}
echo $sample

gatk --java-options "-Xmx68G" \
    SamToFastq \
    -I=${ADAPT_BAM} \
    -FASTQ=$BASE_DIR/outputs/02_merged_mapped/fatsq_adapters/${sample}.adapter.fq \
    -INTERLEAVE=true \
    -NON_PF=true \
    -TMP_DIR=$BASE_DIR/temp_files | \
bwa mem -M -t 8 -p $BASE_DIR/ressources/HP_genome_unmasked_01.fa $BASE_DIR/outputs/02_merged_mapped/fatsq_adapters/${sample}.adapter.fq > $BASE_DIR/outputs/02_merged_mapped/sam_align/${sample}.aligned.sam | \
gatk --java-options "-Xmx68G" \
    MergeBamAlignment \
    --VALIDATION_STRINGENCY SILENT \
    --EXPECTED_ORIENTATIONS FR \
    --ATTRIBUTES_TO_RETAIN X0 \
    -ALIGNED_BAM=$BASE_DIR/outputs/02_merged_mapped/sam_align/${sample}.aligned.sam \
    -UNMAPPED_BAM=$BASE_DIR/outputs/00_ubams/${sample}.ubam.bam \
    -OUTPUT=$BASE_DIR/outputs/02_merged_mapped/bam_merge_map/${sample}.mapped.bam \
    --REFERENCE_SEQUENCE=$BASE_DIR/ressources/HP_genome_unmasked_01.fa.gz \
    -PAIRED_RUN true \
    --SORT_ORDER "unsorted" \
    --IS_BISULFITE_SEQUENCE false \
    --ALIGNED_READS_ONLY false \
    --CLIP_ADAPTERS false \
    --MAX_RECORDS_IN_RAM 2000000 \
    --ADD_MATE_CIGAR true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --UNMAP_CONTAMINANT_READS true \
    -TMP_DIR=$BASE_DIR/temp_files
