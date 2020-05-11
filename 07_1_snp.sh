#!/bin/bash

#SBATCH --job-name=07_1_snp
#SBATCH --partition=carl.p
#SBATCH --output=/user/doau0129/work/chapter1_2/logs/07_1_snp_%A_%a.out
#SBATCH --error=/user/doau0129/work/chapter1_2/logs/07_1_snp_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90G
#SBATCH --time=4-00:00:00

BASE_DIR=/user/doau0129/work/chapter1_2/

gatk --java-options "-Xmx85g" \
     GenotypeGVCFs \
     -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
     -V=$BASE_DIR/outputs/06_cohort_genotyping/cohort.g.vcf.gz \
     -O=$BASE_DIR/07_genotyping/07_1_raw_snp/intermediate.vcf.gz

 gatk --java-options "-Xmx85G" \
     SelectVariants \
     -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
     -V=$BASE_DIR/07_genotyping/07_1_raw_snp/intermediate.vcf.gz \
     --select-type-to-include=SNP \
     -O=$BASE_DIR/outputs/07_genotyping/07_1_raw_snp/raw_var_sites.vcf.gz

 rm $BASE_DIR/07_genotyping/07_1_raw_snp/intermediate.*
