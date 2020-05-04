#!/bin/bash

#SBATCH --job-name=05_genlikely
#SBATCH --partition=carl.p
#SBATCH --array=1-113
#SBATCH --output=/user/doau0129/work/chapter1_2/logs/05_genlikely_%A_%a.out
#SBATCH --error=/user/doau0129/work/chapter1_2/logs/05_genlikely_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=45G
#SBATCH --time=3-00:00:00

BASE_DIR=/user/doau0129/work/chapter1_2/


INPUT_DUPLI=$BASE_DIR/outputs/listoffiles/duplicates_hapcaller.fofn

DUPLI=$(cat ${INPUT_DUPLI} | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo $DUPLI

sample1=${DUPLI%.*}
sample=${sample1%.*}
echo $sample

gatk --java-options "-Xmx35g" HaplotypeCaller  \
     -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
     -I=$BASE_DIR/outputs/03_mark_duplicates/mark_duplicates/duplicates/${DUPLI} \
     -O=$BASE_DIR/outputs/05_genlikely/${sample}.g.vcf.gz \
     -ERC GVCF
