#!/bin/bash

#SBATCH --job-name=08_2_all
#SBATCH --partition=carl.p
#SBATCH --output=/user/doau0129/work/chapter1_2/logs/08_2_all_%A_%a.out
#SBATCH --error=/user/doau0129/work/chapter1_2/logs/08_2_all_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=45G
#SBATCH --time=3-00:00:00

BASE_DIR=/user/doau0129/work/chapter1_2/

ls -1 $BASE_DIR/outputs/07_2_all/all_sites.* > $BASE_DIR/outputs/listoffiles/all_sites_LG0.fofn

awk '{print "-I", $1}' $BASE_DIR/outputs/listoffiles/all_sites_LG0.fofn > $BASE_DIR/outputs/listoffiles/all_sites_LG.fofn

rm $BASE_DIR/outputs/listoffiles/all_sites_LG0.fofn

INPUT_GEN=$(cat $BASE_DIR/outputs/listoffiles/all_sites_LG.fofn)
echo $INPUT_GEN

gatk --java-options "-Xmx85g" \
       GatherVcfs \
       $INPUT \
       -O=$BASE_DIR/outputs/08_2_merge/all_sites.vcf.gz
