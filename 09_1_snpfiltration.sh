#!/bin/bash

#SBATCH --job-name=09_1_filtration
#SBATCH --partition=carl.p
#SBATCH --output=/user/doau0129/work/chapter1_2/logs/09_1_filtration_%A_%a.out
#SBATCH --error=/user/doau0129/work/chapter1_2/logs/09_1_filtration_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --time=12:00:00

BASE_DIR=/user/doau0129/work/chapter1_2/


gatk --java-options "-Xmx75G" \
       VariantFiltration \
       -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
       -V=$BASE_DIR/outputs/07_genotyping/07_1_raw_snp/raw_var_sites.vcf.gz \
       -O=$BASE_DIR/outputs/09_1_snpfiltration/intermediate.vcf.gz \
       --filter-expression "QD < 4.0" \
       --filter-name "filter_QD" \
       --filter-expression "FS > 60.0" \
       --filter-name "filter_FS" \
       --filter-expression "MQ < 57.2 || MQ > 62.2" \
       --filter-name "filter_MQ" \
       --filter-expression "MQRankSum < -0.2 || MQRankSum > 0.2" \
       --filter-name "filter_MQRankSum" \
       --filter-expression "ReadPosRankSum < -2.0 || ReadPosRankSum > 2.0 " \
       --filter-name "filter_ReadPosRankSum"

gatk --java-options "-Xmx75G" \
       SelectVariants \
       -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
       -V=$BASE_DIR/outputs/09_1_snpfiltration/intermediate.vcf.gz \
       -O=$BASE_DIR/outputs/09_1_snpfiltration/intermediate.filterd.vcf.gz \
       --exclude-filtered

vcftools \
       --gzvcf $BASE_DIR/outputs/09_1_snpfiltration/intermediate.filterd.vcf.gz \
       --max-missing-count 17 \
       --max-alleles 2 \
       --stdout  \
       --recode | \
       bgzip > $BASE_DIR/outputs/09_1_snpfiltration/filterd_bi-allelic.vcf.gz

tabix -p vcf $BASE_DIR/outputs/09_1_snpfiltration/filterd_bi-allelic.vcf.gz

rm $BASE_DIR/outputs/09_1_snpfiltration/intermediate.*
