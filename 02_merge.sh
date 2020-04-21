#!/bin/bash

#SBATCH --job-name=02_3_merge	
#SBATCH --partition=carl.p
#SBATCH --array=1-117
#SBATCH --output=/user/doau0129/work/chapter1_2/logs/02_3_merge_%A_%a.out
#SBATCH --error=/user/doau0129/work/chapter1_2/logs/02_3_merge_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=75G
#SBATCH --time=2-00:00:00

BASE_DIR=/user/doau0129/work/chapter1_2/


INPUT_ADAPT_BAMS=$BASE_DIR/outputs/listoffiles/samalign.fofn

ADAPT_BAM=$(cat ${INPUT_ADAPT_BAMS} | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo $ADAPT_BAM

sample1=${ADAPT_BAM%.*}
sample=${sample1%.*}
echo $sample

#gatk --java-options "-Xmx68G" \
#    SamToFastq \
#    -I=$BASE_DIR/outputs/01_adapters/adapters/${ADAPT_BAM} \
#    -FASTQ=$BASE_DIR/outputs/02_merged_map/fatsq_adapters/${sample}.adapterfq \
#    -INTERLEAVE=true \
#    -NON_PF=true \
#    -TMP_DIR=$BASE_DIR/temp_files
#bwa mem -M -t 8 -p $BASE_DIR/ressources/HP_genome_unmasked_01.fa $BASE_DIR/outputs/02_merged_map/fatsq_adapters/${sample}.adapterfq > $BASE_DIR/outputs/02_merged_map/sam_align/${sample}.aligned.sam
gatk --java-options "-Xmx68G" \
    MergeBamAlignment \
    --VALIDATION_STRINGENCY SILENT \
    --EXPECTED_ORIENTATIONS FR \
    --ATTRIBUTES_TO_RETAIN X0 \
    -ALIGNED_BAM=$BASE_DIR/outputs/02_merged_map/sam_align/${sample}.aligned.sam \
    -UNMAPPED_BAM=$BASE_DIR/outputs/00_ubams/${sample}.ubam.bam \
    -OUTPUT=$BASE_DIR/outputs/02_merged_map/bam_merge_alignement/${sample}.mapped.bam \
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
