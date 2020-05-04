#!/bin/bash

#SBATCH --job-name=06_combine
#SBATCH --partition=carl.p
#SBATCH --output=/user/doau0129/work/chapter1_2/logs/06_combine_%A_%a.out
#SBATCH --error=/user/doau0129/work/chapter1_2/logs/06_combine_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90G
#SBATCH --time=4-00:00:00

BASE_DIR=/user/doau0129/work/chapter1_2/

ls -1 $BASE_DIR/outputs/05_genlikely/*.g.vcf.gz > genlikely_individual0.fofn

awk '{print "-V", $1}' $BASE_DIR/outputs/listoffiles/genlikely_individual0.fofn > $BASE_DIR/outputs/listoffiles/genlikely_individual.fofn

rm $BASE_DIR/outputs/listoffiles/genlikely_individual0.fofn

INPUT_GEN=$(cat test3)
echo $INPUT_GEN

gatk --java-options "-Xmx85g" \
      CombineGVCFs \
      -R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
      $INPUT_GEN \
      -O=cohort.g.vcf.gz
