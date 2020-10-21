#!/bin/bash
# by: Floriane Coulmance: 01/04/2020
# usage:
# gxp.sh -i <PATH> -j <JOB_ID>
# ------------------------------------------------------------------------------
# PATH corresponds to the path to the base directory, all outputs and necessary
# folder will be created by the script
# JOB_ID corresponds string ids from where you want  the script to be ran
# ------------------------------------------------------------------------------



# ********** Allow to enter bash options **********
# -------------------------------------------------

while getopts i:j: option
do
case "${option}"
in
i) BASE_DIR=${OPTARG};; # get the base directory path
j) JID_RES=${OPTARG};; # get the jobid from which you want to resume
esac
done



# ********* Create necessary repositories *********
# -------------------------------------------------

# Repo for job logs
mkdir $BASE_DIR/logs/

# Repo for figures
mkdir $BASE_DIR/figures/

# Outputs repo
mkdir $BASE_DIR/outputs/
mkdir $BASE_DIR/outputs/0_ubam/
mkdir $BASE_DIR/outputs/1_adapters/
mkdir $BASE_DIR/outputs/1_adapters/adapters/
mkdir $BASE_DIR/outputs/1_adapters/metrics/
mkdir $BASE_DIR/outputs/2_align/
mkdir $BASE_DIR/outputs/2_align/1_fatsq/
mkdir $BASE_DIR/outputs/2_align/2_align/
mkdir $BASE_DIR/outputs/2_align/3_merge/
mkdir $BASE_DIR/outputs/3_duplicates/
mkdir $BASE_DIR/outputs/3_duplicates/1_sort/
mkdir $BASE_DIR/outputs/3_duplicates/2_tag/
mkdir $BASE_DIR/outputs/3_duplicates/3_mark/
mkdir $BASE_DIR/outputs/3_duplicates/3_mark/duplicates/
mkdir $BASE_DIR/outputs/3_duplicates/3_mark/metrics/
mkdir $BASE_DIR/outputs/3_duplicates/3_mark/duplicates/removed/
mkdir $BASE_DIR/outputs/4_likelihood/
mkdir $BASE_DIR/outputs/5_cohort/
mkdir $BASE_DIR/outputs/6_genotyping/
mkdir $BASE_DIR/outputs/6_genotyping/6_1_snp/
mkdir $BASE_DIR/outputs/6_genotyping/6_2_all/

# Repo for coverage analysis
mkdir $BASE_DIR/outputs/coverage/

# Repo for pca analysis files
mkdir $BASE_DIR/outputs/pca/

# Annex folder for list of files
mkdir $BASE_DIR/outputs/lof/



# ********* Jobs creation *************************
# -------------------------------------------------

# --------------------------- PREPARATION -------------------------------------#

# ------------------------------------------------------------------------------
# Job 0 creates single unaligned bam files from forward and reverse sequencing
# files for each of the 49 samples considered in this study

jobfile0=0_ubam.tmp # temp file
cat > $jobfile0 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=0_ubam
#SBATCH --partition=carl.p
#SBATCH --array=2-118
#SBATCH --output=$BASE_DIR/logs/0_ubam_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/0_ubam_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=04:00:00


INPUT_META=$BASE_DIR/metadata/metadata_gxp_ben_floridae_complete
LINES=\$(cat \${INPUT_META} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)

IFS=";" read -r -a array <<< "\$LINES"
label=\${array[1]}
company=\${array[7]}
frwdfile=\${array[9]}
flowcellidfrwd=\${array[15]}
lanefrwd=\${array[16]}
revfile=\${array[24]}

echo -e "----------------------------"
echo -e "Label:\t\${label}\nFwd:\t\${frwdfile}\nRev:\t\${revfile}"
echo -e "Flowcell:\t\${flowcellidfrwd}\nLane:\t\${lanefrwd}"
echo -e "Read group:\t\${flowcellidfrwd}.\${lanefrwd}\nCompany:\t\${company}"

gatk --java-options "-Xmx20G" \
    FastqToSam \
    -SM=\${label} \
    -F1=$BASE_DIR/data/\${frwdfile} \
    -F2=$BASE_DIR/data/\${revfile} \
    -O=$BASE_DIR/outputs/0_ubam/\${label}.\${lanefrwd}.ubam.bam \
    -RG=\${label}.\${lanefrwd} \
    -LB=\${label}".lib1" \
    -PU=\${flowcellidfrwd}.\${lanefrwd} \
    -PL=Illumina \
    -CN=\${company} \
    --TMP_DIR=$BASE_DIR/temp_files

#Make a file of files
ls -1 $BASE_DIR/outputs/0_ubam/ > $BASE_DIR/outputs/lof/0_ubam.fofn


EOA



# ------------------------------------------------------------------------------
# Job 1 marks adapters and creates the corresponding files
# in /1_adapters/adapters and the metric files in /1_adapters/metrics

jobfile1=1_adapters.tmp # temp file
cat > $jobfile1 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=1_adapters
#SBATCH --partition=carl.p
#SBATCH --array=1-117
#SBATCH --output=$BASE_DIR/logs/1_adapters_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/1_adapters_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=04:00:00


INPUT_UBAMS=$BASE_DIR/outputs/lof/0_ubam.fofn
UBAMS=\$(cat \${INPUT_UBAMS} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \$UBAMS

sample1=\${UBAMS%.*}
sample=\${sample1%.*}
echo \$sample

gatk --java-options "-Xmx18G" \
   MarkIlluminaAdapters \
   -I=$BASE_DIR/outputs/0_ubam/\${UBAMS} \
   -O=$BASE_DIR/outputs/1_adapters/adapters/\${sample}.adapter.bam \
   -M=$BASE_DIR/outputs/1_adapters/metrics/\${sample}.adapter.metrics.txt \
   -TMP_DIR=$BASE_DIR/temp_files

ls -1 $BASE_DIR/outputs/1_adapters/adapters/ > $BASE_DIR/outputs/lof/1_adapters.fofn


EOA



# ------------------------------------------------------------------------------
# Job 2 uses reference genome to align the reads

jobfile2=2_align.tmp # temp file
cat > $jobfile2 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=2_align
#SBATCH --partition=mpcb.p
#SBATCH --array=1-117
#SBATCH --output=$BASE_DIR/logs/2_align_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/2_align_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=30G
#SBATCH --time=2-00:00:00


INPUT_ADAPT_BAMS=$BASE_DIR/outputs/lof/1_adapters.fofn
ADAPT_BAM=\$(cat \${INPUT_ADAPT_BAMS} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \$ADAPT_BAM

sample1=\${ADAPT_BAM%.*}
sample=\${sample1%.*}
echo \$sample

gatk --java-options "-Xmx68G" \
    SamToFastq \
    -I=$BASE_DIR/outputs/1_adapters/adapters/\${ADAPT_BAM} \
    -FASTQ=/dev/stdout \
    -INTERLEAVE=true \
    -NON_PF=true \
    -TMP_DIR=$BASE_DIR/temp_files | \
bwa mem -M -t 8 -p $BASE_DIR/ressources/HP_genome_unmasked_01.fa /dev/stdin |
gatk --java-options "-Xmx68G" \
    MergeBamAlignment \
    --VALIDATION_STRINGENCY SILENT \
    --EXPECTED_ORIENTATIONS FR \
    --ATTRIBUTES_TO_RETAIN X0 \
    -ALIGNED_BAM=/dev/stdin \
    -UNMAPPED_BAM=$BASE_DIR/outputs/0_ubam/\${sample}.ubam.bam \
    -OUTPUT=$BASE_DIR/outputs/2_align/\${sample}.mapped.bam \
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

ls -1 $BASE_DIR/outputs/2_align/ > $BASE_DIR/outputs/lof/2_align.fofn


EOA



# ------------------------------------------------------------------------------
# Job 3 tag

jobfile3=3_tag.tmp # temp file
cat > $jobfile3 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=3_tag
#SBATCH --partition=carl.p
#SBATCH --array=1-117
#SBATCH --output=$BASE_DIR/logs/3_tag_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/3_tag_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=3-00:00:00


INPUT_MAPPED_BAM=$BASE_DIR/outputs/lof/2_align.fofn
MAPPED_BAM=\$(cat \${INPUT_MAPPED_BAM} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \$MAPPED_BAM

sample1=\${MAPPED_BAM%.*}
sample=\${sample1%.*}
echo \$sample

gatk --java-options "-Xmx30G" \
		SortSam \
		-I=$BASE_DIR/outputs/2_align/3_merge/\${MAPPED_BAM} \
		-O=/dev/stdout \
		--SORT_ORDER="coordinate" \
		--CREATE_INDEX=false \
		--CREATE_MD5_FILE=false \
		-TMP_DIR=$BASE_DIR/temp_files | \
gatk --java-options "-Xmx30G" \
    SetNmMdAndUqTags \
    --INPUT=/dev/stdin \
    --OUTPUT=$BASE_DIR/outputs/3_duplicates/2_tag/\${sample}.intermediate.bam \
    --CREATE_INDEX=false \
    --CREATE_MD5_FILE=false \
    -TMP_DIR=$BASE_DIR/temp_files \
    --REFERENCE_SEQUENCE=$BASE_DIR/ressources/HP_genome_unmasked_01.fa.gz

ls -1 $BASE_DIR/outputs/3_duplicates/2_tag/*.bam |xargs -n1 basename > $BASE_DIR/outputs/lof/3_tag.fofn

EOA



# ------------------------------------------------------------------------------
# Job 4 mark duplicates

jobfile4=4_mark.tmp # temp file
cat > $jobfile4 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=4_mark
#SBATCH --partition=carl.p
#SBATCH --array=1-117
#SBATCH --output=$BASE_DIR/logs/4_mark_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/4_mark_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=15:00:00


INPUT_TAGSINTER=$BASE_DIR/outputs/lof/3_tag.fofn
TAGSINTER=\$(cat \${INPUT_TAGSINTER} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \$TAGSINTER

sample1=\${TAGSINTER%.*}
sample=\${sample1%.*}
echo \$sample

gatk --java-options "-Xmx30G" \
  MarkDuplicates \
  -I=$BASE_DIR/outputs/3_duplicates/2_tag/\${TAGSINTER} \
  -O=$BASE_DIR/outputs/3_duplicates/3_mark/duplicates/\${sample}.dedup.bam \
  -M=$BASE_DIR/outputs/3_duplicates/3_mark/metrics/\${sample}.dedup.metrics.txt \
  -MAX_FILE_HANDLES=1000  \
  -TMP_DIR=$BASE_DIR/temp_files

ls -1 $BASE_DIR/outputs/3_duplicates/3_mark/duplicates/*.bam > $BASE_DIR/outputs/lof/4_duplicates.fofn


EOA


# ------------------------------------------------------------------------------
# Job 5 creates an index file for the aligned samples with adapters and marked
# duplicates

jobfile5=5_index.tmp # temp file
cat > $jobfile5 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=5_index
#SBATCH --partition=carl.p
#SBATCH --array=1-117
#SBATCH --output=$BASE_DIR/logs/5_index_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/5_index_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=01:00:00


INPUT_DUPLI=$BASE_DIR/outputs/lof/4_duplicates.fofn
DUPLI=\$(cat \${INPUT_DUPLI} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \$DUPLI

sample1=\${DUPLI%.*}
sample=\${sample1%.*}
echo \$sample

gatk --java-options "-Xmx30G" \
  BuildBamIndex \
  -INPUT=\${DUPLI}


EOA



# --------------------------- COVERAGE ----------------------------------------#

# ------------------------------------------------------------------------------
# Job a calculates coverage for each of the samples

jobfilea=a_coverage.tmp # temp file
cat > $jobfilea <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=a_coverage
#SBATCH --partition=carl.p
#SBATCH --array=1-117
#SBATCH --output=$BASE_DIR/logs/a_coverage_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/a_coverage_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=02:00:00


INPUT_COV=$BASE_DIR/outputs/lof/3_3_duplicates.fofn
COV=\$(cat \${INPUT_COV} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \$COV

sample1=\${COV%.*}
sample=\${sample1%.*}
echo \$sample

gatk --java-options "-Xmx30G" \
  CollectWgsMetrics \
  -I=$BASE_DIR/outputs/3_duplicates/3_mark/duplicates/\${COV} \
  -O=$BASE_DIR/outputs/coverage/\${sample}.wgsmetrics.txt \
  -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa.gz


EOA



# ------------------------------------------------------------------------------
# Job b creates the mean coverage table and the histogram and runs the R script
# for the coverage table, histogram and returns the path+names of the sample
# files that have a too low coverage to be analysed

jobfileb=b_histo.tmp # temp file
cat > $jobfileb <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=b_histo
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/b_histo_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/b_histo_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=00:15:00

module load  hpc-env/8.3
module load R/4.0.2-foss-2019b

FILES=$BASE_DIR/outputs/coverage/*
echo \$FILES

for f in \$FILES
do
  echo \$f
  name_file=\${f##*/}
  echo \$name_file
  sample1=\${name_file%.*}
  sample2=\${sample1%.*}
  sample=\${sample2%.*}
  echo \$sample
  cov=\$(cat \$f | awk 'FNR == 8 {print \$2}')
  echo \$cov

  echo "\$sample \$cov" >> $BASE_DIR/outputs/coverage/coverage_table
done

Rscript --vanilla $BASE_DIR/R/coverage_hist.R $BASE_DIR/outputs/coverage/coverage_table $BASE_DIR/figures/ $BASE_DIR/outputs/lof/

while read line
do
  echo \$line
  ls $BASE_DIR/outputs/3_duplicates/3_mark/duplicates/ | grep \$line >> $BASE_DIR/outputs/lof/c_remove.fofn
done < $BASE_DIR/outputs/lof/b_remove.fofn


EOA



# ------------------------------------------------------------------------------
# Job c removes the files that have too low coverage from the analysis

jobfilec=c_remove.tmp # temp file
cat > $jobfilec <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=c_remove
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/c_remove_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/c_remove_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=02:00:00


for file in \$(cat $BASE_DIR/outputs/lof/c_remove.fofn)
do
  echo \$file
  mv  $BASE_DIR/outputs/3_duplicates/3_mark/duplicates/\$file $BASE_DIR/outputs/3_duplicates/3_mark/duplicates/removed/
done

ls -1 $BASE_DIR/outputs/3_duplicates/3_mark/duplicates/*.bam |xargs -n1 basename > $BASE_DIR/outputs/lof/3_3_new_duplicates.fofn


EOA



# --------------------------- GENOTYPING --------------------------------------#

# Create array size based on files not removed at previous step
# SIZE=$(wc $BASE_DIR/outputs/lof/3_3_new_duplicates.fofn | awk '{print $1}')
# echo $SIZE

# ------------------------------------------------------------------------------
# Job 6 calculates genotype likelihoods

jobfile6=6_likelihood.tmp # temp file
cat > $jobfile6 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=6_likelihood
#SBATCH --partition=carl.p
#SBATCH --array=1-$SIZE
#SBATCH --output=$BASE_DIR/logs/6_likelihood_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/6_likelihood_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=3-00:00:00


INPUT_DUPLI=$BASE_DIR/outputs/lof/3_3_new_duplicates.fofn
DUPLI=\$(cat \${INPUT_DUPLI} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \$DUPLI

sample1=\${DUPLI%.*}
sample=\${sample1%.*}
echo \$sample

gatk --java-options "-Xmx35g" HaplotypeCaller  \
     -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
     -I=$BASE_DIR/outputs/3_duplicates/3_mark/duplicates/\${DUPLI} \
     -O=$BASE_DIR/outputs/4_likelihood/\${sample}.g.vcf.gz \
     -ERC GVCF

ls -1 $BASE_DIR/outputs/4_likelihood/*.g.vcf.gz > $BASE_DIR/outputs/lof/4_likelihood.fofn


EOA



# ------------------------------------------------------------------------------
# Job 7 combines all the samples'gvcf files into one cohort one

jobfile7=7_cohort.tmp # temp file
cat > $jobfile7 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=7_cohort
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/7_cohort_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/7_cohort_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --time=4-00:00:00


awk '{print "-V", $1}' $BASE_DIR/outputs/lof/4_likelihood.fofn > $BASE_DIR/outputs/lof/4_likelihood2.fofn
INPUT_GEN=$(cat $BASE_DIR/outputs/lof/4_likelihood2.fofn)
echo \${INPUT_GEN}

gatk --java-options "-Xmx85g" \
      CombineGVCFs \
      -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
      \${INPUT_GEN} \
      -O=$BASE_DIR/outputs/5_cohort/cohort.g.vcf.gz


EOA



# --------------------------- SNP ONLY ----------------------------------------#

# ------------------------------------------------------------------------------
# Job 8 genotype and select just SNPs variants

jobfile8=8_snp.tmp # temp file
cat > $jobfile8 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=8_snp
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/8_snp_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/8_snp_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90G
#SBATCH --time=4-00:00:00


gatk --java-options "-Xmx85g" \
     GenotypeGVCFs \
     -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
     -V=$BASE_DIR/outputs/5_cohort/cohort.g.vcf.gz \
     -O=$BASE_DIR/6_genotyping/6_1_snp/intermediate.vcf.gz

gatk --java-options "-Xmx85G" \
     SelectVariants \
     -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
     -V=$BASE_DIR/6_genotyping/6_1_snp/intermediate.vcf.gz \
     --select-type-to-include=SNP \
     -O=$BASE_DIR/outputs/6_genotyping/6_1_snp/raw_var_sites.vcf.gz

rm $BASE_DIR/6_genotyping/6_1_snp/intermediate.*


EOA



# ------------------------------------------------------------------------------
# Job 9 creates metric table to look before filtering

jobfile9=9_snp_metrics.tmp # temp file
cat > $jobfile9 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=9_snp_metrics
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/9_snp_metrics_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/9_snp_metrics_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=06:00:00


gatk --java-options "-Xmx25G" \
       VariantsToTable \
       --variant=$BASE_DIR/outputs/6_genotyping/6_1_snp/raw_var_sites.vcf.gz \
       --output=$BASE_DIR/outputs/6_genotyping/6_1_snp/raw_var_sites.table.txt \
       -F=CHROM -F=POS -F=MQ \
       -F=QD -F=FS -F=MQRankSum -F=ReadPosRankSum \
       --show-filtered

Rscript --vanilla $BASE_DIR/R/filtering_thresholds.R $BASE_DIR/outputs/6_genotyping/6_1_snp/raw_var_sites.table.txt $BASE_DIR/figures/

EOA



# ------------------------------------------------------------------------------
# Job 10 filter the SNP variants based on metrics from previous step, and
# creates the vcffiles for casz1 gene SNPs only

jobfile10=10_snp_filter.tmp # temp file
cat > $jobfile10 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=10_snp_filter
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/10_snp_filter_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/10_snp_filter_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --time=12:00:00


gatk --java-options "-Xmx75G" \
       VariantFiltration \
       -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
       -V=$BASE_DIR/outputs/6_genotyping/6_1_snp/raw_var_sites.vcf.gz \
       -O=$BASE_DIR/outputs/6_genotyping/6_1_snp/intermediate.vcf.gz \
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
       -V=$BASE_DIR/outputs/6_genotyping/6_1_snp/intermediate.vcf.gz \
       -O=$BASE_DIR/outputs/6_genotyping/6_1_snp/intermediate.filterd.vcf.gz \
       --exclude-filtered

vcftools \
       --gzvcf $BASE_DIR/outputs/6_genotyping/6_1_snp/intermediate.filterd.vcf.gz \
       --max-missing-count 17 \
       --max-alleles 2 \
       --stdout  \
       --recode | \
       bgzip > $BASE_DIR/outputs/6_genotyping/6_1_snp/filterd_bi-allelic.vcf.gz

tabix -p vcf $BASE_DIR/outputs/6_genotyping/6_1_snp/filterd_bi-allelic.vcf.gz

rm $BASE_DIR/outputs/6_genotyping/6_1_snp/intermediate.*

echo -e "$BASE_DIR/outputs/6_genotyping/6_1_snp/filterd_bi-allelic.vcf.gz" > $BASE_DIR/outputs/lof/snp_all.fofn

EOA



# --------------------------- ALL SITES ---------------------------------------#

# ------------------------------------------------------------------------------
# Job 11 genotype and select all variants types

jobfile11=11_all.tmp # temp file
cat > $jobfile11 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=11_all
#SBATCH --partition=carl.p
#SBATCH --array=1-24
#SBATCH --output=$BASE_DIR/logs/11_all_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/11_all_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90G
#SBATCH --time=4-00:00:00


gatk --java-options "-Xmx85g" \
    GenotypeGVCFs \
    -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
    -L=LG${SLURM_ARRAY_TASK_ID} \
    -V=$BASE_DIR/outputs/5_cohort/cohort.g.vcf.gz  \
    -O=$BASE_DIR/6_genotyping/6_2_all/intermediate.vcf.gz \
    --include-non-variant-sites=true

gatk --java-options "-Xmx85G" \
    SelectVariants \
    -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
    -V=$BASE_DIR/6_genotyping/6_2_all/intermediate.vcf.gz \
    --select-type-to-exclude=INDEL \
    -O=$BASE_DIR/6_genotyping/6_2_all/all_sites.LG${SLURM_ARRAY_TASK_ID}.vcf.gz

rm $BASE_DIR/6_genotyping/6_2_all/intermediate.*


EOA



# ------------------------------------------------------------------------------
# Job 12 generates one file from all LG files genotypes for all variant sites

jobfile12=12_all_gather.tmp # temp file
cat > $jobfile12 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=12_all_gather
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/12_all_gather_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/12_all_gather_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=45G
#SBATCH --time=3-00:00:00


ls -1 $BASE_DIR/outputs/6_genotyping/6_2_all/all_sites.* > $BASE_DIR/outputs/lof/14_all.fofn
awk '{print "-I", $1}' $BASE_DIR/outputs/lof/14_all.fofn > $BASE_DIR/outputs/lof/14_all2.fofn

INPUT_GEN=$(cat $BASE_DIR/outputs/lof/14_all2.fofn)
echo \${INPUT_GEN}

gatk --java-options "-Xmx85g" \
       GatherVcfs \
       \${INPUT_GEN} \
       -O=$BASE_DIR/outputs/6_genotyping/6_2_all/all_sites.vcf.gz


EOA



# ------------------------------------------------------------------------------
# Job 13 generates one file from all LG files genotypes for all variant sites

jobfile13=13_all_filter.tmp # temp file
cat > $jobfile13 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=13_all_filter
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/13_all_filter_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/13_all_filter_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=150G
#SBATCH --time=30:00:00


tabix -p vcf $BASE_DIR/outputs/6_genotyping/6_2_all/all_sites.vcf.gz

gatk --java-options "-Xmx75G" \
       VariantFiltration \
       -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
       -V=$BASE_DIR/outputs/6_genotyping/6_2_all/all_sites.vcf.gz \
       -O=$BASE_DIR/outputs/6_genotyping/6_2_all/intermediate.vcf.gz \
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
       --QUIET true &> var_filt.log

gatk --java-options "-Xmx75G" \
       SelectVariants \
       -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
       -V=$BASE_DIR/outputs/6_genotyping/6_2_all/intermediate.vcf.gz \
       -O=$BASE_DIR/outputs/6_genotyping/6_2_all/intermediate.filterd.vcf.gz \
       --exclude-filtered \
       --QUIET true \
       --verbosity ERROR  &> var_select.log

vcftools \
       --gzvcf $BASE_DIR/outputs/6_genotyping/6_2_all/intermediate.filterd.vcf.gz \
       --max-missing-count 17 \
       --stdout  \
       --recode | \
       bgzip > $BASE_DIR/outputs/6_genotyping/6_2_all/filterd.allBP.vcf.gz

rm $BASE_DIR/outputs/6_genotyping/6_2_all/intermediate.*

echo -e "\n$BASE_DIR/outputs/6_genotyping/6_1_snp/filterd.allBP.vcf.gz" >> $BASE_DIR/outputs/lof/snp_all.fofn



EOA



# ------------------------------------------------------------------------------
# Job 14 filter the SNP variants based on metrics from previous step, and
# creates the vcffiles for casz1 gene SNPs only

jobfile17=14_changes.tmp # temp file
cat > $jobfile14 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=14_changes
#SBATCH --partition=carl.p
#SBATCH --array=1-2
#SBATCH --output=$BASE_DIR/logs/14_changes_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/14_changes_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90G
#SBATCH --time=4-00:00:00


INPUT_GENO=$BASE_DIR/outputs/lof/snp_all.fofn
GENO=\$(cat \${INPUT_GENO} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
PATH=\${GENO%/*}
echo \$GENO
echo \$PATH

if [ "\${SLURM_ARRAY_TASK_ID}" == "1" ]
then
  PREFIX="snp"
elif [ "\${SLURM_ARRAY_TASK_ID}" == "2" ]
  PREFIX="all"
fi

echo \${PREFIX}

echo "PL17_35puepue PL17_35indpue" > $BASE_DIR/outputs/lof/change_sample.txt

bcftools \
       reheader \
       --samples $BASE_DIR/outputs/lof/change_sample.txt \
       -o \${PATH}/\${PREFIX}_filterd.vcf.gz \
       \${GENO}

LG=\$(zless $BASE_DIR/annotations/HP.annotation.named.LG12.gff.gz | grep -w gene | grep -i casz1 | awk '{print \$1}')
START=\$(zless $BASE_DIR/annotations/HP.annotation.named.LG12.gff.gz | grep -w gene | grep -i casz1 | awk '{print \$4}')
END=\$(zless $BASE_DIR/annotations/HP.annotation.named.LG12.gff.gz | grep -w gene | grep -i casz1 | awk '{print \$5}')
echo \$LG
echo \$START
echo \$END

vcftools \
      --vcf \${PATH}/\${PREFIX}_filterd.vcf.gz \
      --chr \${LG} \
      --from-bp \${START} \
      --to-bp \${END} \
      --recode \
      --out \${PATH}/\${PREFIX}_filterd_casz1.vcf.gz

echo -e "\${PATH}/\${PREFIX}_filterd.vcf.gz\n\${PATH}/\${PREFIX}_filterd_casz1.vcf.gz\n" >> $BASE_DIR/outputs/lof/17_pca.txt


EOA



# ------------------------------------------------------------------------------
# Job d creates the PCA files input for the plots

jobfiled=d_pca.tmp # temp file
cat > $jobfiled <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=d_pca
#SBATCH --partition=carl.p
#SBATCH --array=1-4
#SBATCH --output=$BASE_DIR/logs/d_pca_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/d_pca_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-02:00:00


INPUT_PCA=$BASE_DIR/outputs/lof/17_pca.fofn
PCA=\$(cat \${INPUT_PCA} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
FILE=\${*/%PCA}
VCF=\${FILE%.*}
PREFIX=\${VCF%.*}
echo \$PCA
echo \$FILE
echo \$VCF
echo \$PREFIX

module load  hpc-env/8.3
module load R/4.0.2-foss-2019b

Rscript --vanilla $BASE_DIR/R/pca.R \${PCA} $BASE_DIR/pca/ \${PREFIX}


EOA


# ********** Schedule the job launching ***********
# -------------------------------------------------

# ------------------------------------------------------------------------------
# DATA PREPARATION

if [ "$JID_RES" = "jid1" ] || [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jida" ] || [ "$JID_RES" = "jidb" ] || [ "$JID_RES" = "jidc" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ] || [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jidd" ];
then
  echo "*****   0_ubam          : DONE         **"
else
  jid0=$(sbatch ${jobfile0})
fi


if [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jida" ] || [ "$JID_RES" = "jidb" ] || [ "$JID_RES" = "jidc" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ] || [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jidd" ]
then
  echo "*****   1_adapters      : DONE         **"
elif [ "$JID_RES" = jid1 ]
then
  jid1=$(sbatch ${jobfile1})
else
  jid1=$(sbatch --dependency=afterok:${jid0##* } ${jobfile1})
fi


if [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jida" ] || [ "$JID_RES" = "jidb" ] || [ "$JID_RES" = "jidc" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ] || [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jidd" ]
then
  echo "*****   2_align         : DONE         **"
elif [ "$JID_RES" = jid2 ]
then
  jid2=$(sbatch ${jobfile2})
else
  jid2=$(sbatch --dependency=afterok:${jid1##* } ${jobfile2})
fi


if [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jida" ] || [ "$JID_RES" = "jidb" ] || [ "$JID_RES" = "jidc" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ] || [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jidd" ]
then
  echo "*****   3_tag           : DONE         **"
elif [ "$JID_RES" = "jid3" ]
then
  jid3=$(sbatch ${jobfile3})
else
  jid3=$(sbatch --dependency=afterok:${jid2##* } ${jobfile3})
fi


if [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jida" ] || [ "$JID_RES" = "jidb" ] || [ "$JID_RES" = "jidc" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ] || [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jidd" ]
then
  echo "*****   4_mark          : DONE         **"
elif [ "$JID_RES" = "jid4" ]
then
  jid4=$(sbatch ${jobfile4})
else
  jid4=$(sbatch --dependency=afterok:${jid3##* } ${jobfile4})
fi


if [ "$JID_RES" = "jida" ] || [ "$JID_RES" = "jidb" ] || [ "$JID_RES" = "jidc" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ] || [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jidd" ]
then
  echo "*****   5_index         : DONE         **"
elif [ "$JID_RES" = "jid5" ]
then
  jid5=$(sbatch ${jobfile5})
else
  jid5=$(sbatch --dependency=afterok:${jid4##* } ${jobfile5})
fi


if [ "$JID_RES" = "jidb" ] || [ "$JID_RES" = "jidc" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ] || [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jidd" ]
then
  echo "*****   a_coverage      : DONE         **"
elif [ "$JID_RES" = "jida" ]
then
  jida=$(sbatch ${jobfilea})
else
  jida=$(sbatch --dependency=afterok:${jid5##* } ${jobfilea})
fi


if [ "$JID_RES" = "jidc" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ] || [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jidd" ]
then
  echo "*****   b_remove      : DONE         **"
elif [ "$JID_RES" = "jidb" ]
then
  jidb=$(sbatch ${jobfileb})
else
  jidb=$(sbatch --dependency=afterok:${jida##* } ${jobfileb})
fi


if [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ] || [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jidd" ]
then
  echo "*****   c_histo         : DONE         **"
elif [ "$JID_RES" = "jidc" ]
then
  jidc=$(sbatch ${jobfilec})
else
  jidc=$(sbatch --dependency=afterok:${jidb##* } ${jobfilec})
fi


# ------------------------------------------------------------------------------
# GENOTYPING

if [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ] || [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jidd" ]
then
  echo "*****   6_likelihood    : DONE         **"
elif [ "$JID_RES" = "jid6" ]
then
  # Create array size based on files not removed at previous step
  SIZE=$(wc $BASE_DIR/outputs/lof/3_3_new_duplicates.fofn | awk '{print $1}')
  echo $SIZE
  jid6=$(sbatch ${jobfile6})
else
  # Create array size based on files not removed at previous step
  SIZE=$(wc $BASE_DIR/outputs/lof/3_3_new_duplicates.fofn | awk '{print $1}')
  echo $SIZE
  jid6=$(sbatch --dependency=afterok:${jidc##* } ${jobfile6})
fi


if [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ] || [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jidd" ]
then
  echo "*****   7_cohort        : DONE         **"
elif [ "$JID_RES" = "jid7" ]
then
  jid7=$(sbatch ${jobfile7})
else
  jid7=$(sbatch --dependency=afterok:${jid6##* } ${jobfile7})
fi


# ------------------------------------------------------------------------------
# FILTER FOR SNPs ONLY

if [ "$JID_RES" = "jid9" ] || [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jidd" ]
then
  echo "*****   8_snp           : DONE         **"
elif [ "$JID_RES" = "jid8" ]
then
  jid8=$(sbatch ${jobfile8})
else
  jid8=$(sbatch --dependency=afterok:${jid7##* } ${jobfile8})
fi


if [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jidd" ]
then
  echo "*****   9_snp_metrics   : DONE         **"
elif [ "$JID_RES" = "jid9" ]
then
  jid9=$(sbatch ${jobfile9})
else
  jid9=$(sbatch --dependency=afterok:${jid8##* } ${jobfile9})
fi


# ------------------------------------------------------------------------------
# FILTER FOR ALL VARIANTS

if [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jidd" ]
then
  echo "*****   11_all           : DONE         **"
elif [ "$JID_RES" = "jid11" ]
then
  jid11=$(sbatch ${jobfile11})
else
  jid11=$(sbatch --dependency=afterok:${jid7##* } ${jobfile11})
fi


if [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jidd" ]
then
  echo "*****   12_all_gather    : DONE         **"
elif [ "$JID_RES" = "jid12" ]
then
  jid12=$(sbatch ${jobfile12})
else
  jid12=$(sbatch --dependency=afterok:${jid11##* } ${jobfile12})
fi


# ------------------------------------------------------------------------------
# LAST CHANGES AND PCA

if [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jidd" ]
then
  echo "*****   10_snp_filter    : DONE         **"
  echo "*****   13_all_filter    : DONE         **"
elif [ "$JID_RES" = "jid10" ]
then
  jid10=$(sbatch ${jobfile10})
elif [ "$JID_RES" = "jid13" ]
then
  jid13=$(sbatch ${jobfile13})
else
  jid10=$(sbatch --dependency=afterok:${jid9##* } ${jobfile10})
  jid13=$(sbatch --dependency=afterok:${jid12##* } ${jobfile13})
fi


if [ "$JID_RES" = "jidd" ];
then
  echo "*****   14_changes       : DONE         **"
elif [ "$JID_RES" = "jid14" ]
then
  jid14=$(sbatch ${jobfile14})
else
  jid14=$(sbatch --dependency=afterok:${jid13##* } ${jobfile14})
fi


if [ "$JID_RES" = "jidd" ];
then
  jidd=$(sbatch ${jobfiled})
else
  jidd=$(sbatch --dependency=afterok:${jid14##* } ${jobfiled})
fi



# ******** Remove useless files and folders *******
# -------------------------------------------------

rm *tmp
# rm -r $BASE_DIR/outputs/lof/
