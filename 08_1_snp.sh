#!/bin/bash

#SBATCH --job-name=08_1_snp
#SBATCH --partition=carl.p
#SBATCH --output=/user/doau0129/work/chapter1_2/logs/08_1_snp_%A_%a.out
#SBATCH --error=/user/doau0129/work/chapter1_2/logs/08_1_snp_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=06:00:00

BASE_DIR=/user/doau0129/work/chapter1_2/

gatk --java-options "-Xmx25G" \
       VariantsToTable \
       --variant=$BASE_DIR/outputs/07_genotyping/07_1_raw_snp/raw_var_sites.vcf.gz \
       --output=$BASE_DIR/outputs/08_1_variants_metrics/raw_var_sites.table.txt \
       -F=CHROM -F=POS -F=MQ \
       -F=QD -F=FS -F=MQRankSum -F=ReadPosRankSum \
       --show-filtered
