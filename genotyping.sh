#!/bin/bash

# Allow to enter arguments
while getopts i:j: option
do
case "${option}"
in
i) BASE_DIR=${OPTARG};; # get the base directory path
j) JID_RES=${OPTARG};; # get the jobid from which you want to resume
esac
done



#author: Floriane (floriane.coulmance@cri-paris.org)
#version: 0.0.1



# Create the output repositories where will be stored files from each steps
mkdir $BASE_DIR/outputs/
mkdir $BASE_DIR/logs/



# ---------------- DATA PREPARATION ---------------- #



  # *** SPLIT SAMPLES *** #

# Create the repo for split_samples
mkdir $BASE_DIR/outputs/00_ubams/


jobfile1=00_split_samples.tmp # temp file
cat > $jobfile1 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=00_readgroups_ubams
#SBATCH --partition=carl.p
#SBATCH --array=2-118
#SBATCH --output=$BASE_DIR/logs/00_readgroups_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/00_readgroups_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=01:00:00

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
    -O=$BASE_DIR/outputs/00_ubams/\${label}.\${lanefrwd}.ubam.bam \
    -RG=\${label}.\${lanefrwd} \
    -LB=\${label}".lib1" \
    -PU=\${flowcellidfrwd}.\${lanefrwd} \
    -PL=Illumina \
    -CN=\${company} \
    --TMP_DIR=$BASE_DIR/temp_files

EOA



  # *** MARK ADAPTERS *** #

#Create necessary folders for the mark adapters step
mkdir $BASE_DIR/outputs/01_adapters/
mkdir $BASE_DIR/outputs/01_adapters/adapters/
mkdir $BASE_DIR/outputs/01_adapters/metrics/
mkdir $BASE_DIR/outputs/listoffiles/

#Create the jobs that mark adapters
jobfile2=01_mark_adapters.tmp # temp file
cat > $jobfile2 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=01_mark_adapters
#SBATCH --partition=carl.p
#SBATCH --array=1-117
#SBATCH --output=$BASE_DIR/logs/01_mark_adapters_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/01_mark_adapters_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=04:00:00

#Make a file of files
ls -1 $BASE_DIR/outputs/00_ubams/ > $BASE_DIR/outputs/listoffiles/01_ubams.fofn

INPUT_UBAMS=$BASE_DIR/outputs/listoffiles/01_ubams.fofn

UBAMS=\$(cat \${INPUT_UBAMS} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \$UBAMS

sample1=\${UBAMS%.*}
sample=\${sample1%.*}
echo \$sample

gatk --java-options "-Xmx18G" \
   MarkIlluminaAdapters \
   -I=$BASE_DIR/outputs/00_ubams/\${UBAMS} \
   -O=$BASE_DIR/outputs/01_adapters/adapters/\${sample}.adapter.bam \
   -M=$BASE_DIR/outputs/01_adapters/metrics/\${sample}.adapter.metrics.txt \
   -TMP_DIR=$BASE_DIR/temp_files

EOA



  # *** ALIGN TO REFERENCE GENOME *** #

mkdir $BASE_DIR/outputs/02_align_reference/

#---- Convert previous created sam to fastaq files ----#

mkdir $BASE_DIR/outputs/02_align_reference/02_1_fatsq_adapters/

jobfile3=02_1_samtofq.tmp # temp file
cat > $jobfile3 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=02_1_samtofq
#SBATCH --partition=carl.p
#SBATCH --array=1-117
#SBATCH --output=$BASE_DIR/logs/02_1_samtofq_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/02_1_samtofq_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=75G
#SBATCH --time=03:00:00

ls -1 $BASE_DIR/outputs/01_adapters/adapters/ > $BASE_DIR/outputs/listoffiles/02_1_adapterbams.fofn

INPUT_ADAPT_BAMS=$BASE_DIR/outputs/listoffiles/02_1_adapterbams.fofn

ADAPT_BAM=\$(cat \${INPUT_ADAPT_BAMS} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \$ADAPT_BAM

sample1=\${ADAPT_BAM%.*}
sample=\${sample1%.*}
echo \$sample

gatk --java-options "-Xmx68G" \
    SamToFastq \
    -I=$BASE_DIR/outputs/01_adapters/adapters/\${ADAPT_BAM} \
    -FASTQ=$BASE_DIR/outputs/02_align_reference/02_1_fatsq_adapters/\${sample}.adapterfq \
    -INTERLEAVE=true \
    -NON_PF=true \
    -TMP_DIR=$BASE_DIR/temp_files

EOA


#---- Align to reference genome ----#

mkdir $BASE_DIR/outputs/02_align_reference/02_2_samalign/

jobfile4=02_2_samalign.tmp # temp file
cat > $jobfile4 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=02_2_samalign
#SBATCH --partition=mpcb.p
#SBATCH --array=1-117
#SBATCH --output=$BASE_DIR/logs/02_2_samalign_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/02_2_samalign_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=15G
#SBATCH --time=2-00:00:00

ls -1 $BASE_DIR/outputs/02_align_reference/02_1_fatsq_adapters/ > $BASE_DIR/outputs/listoffiles/02_2_fastq.fofn

INPUT_ADAPT_BAMS=$BASE_DIR/outputs/listoffiles/02_2_fastq.fofn

ADAPT_BAM=\$(cat \${INPUT_ADAPT_BAMS} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \$ADAPT_BAM

sample=\${ADAPT_BAM%.*}
echo \$sample

bwa mem -M -t 8 -p $BASE_DIR/ressources/HP_genome_unmasked_01.fa $BASE_DIR/outputs/02_align_reference/02_1_fatsq_adapters/\${ADAPT_BAM} > $BASE_DIR/outputs/02_align_reference/02_2_samalign/\${sample}.aligned.sam

EOA


#---- Merge aligned bam with marked adqapters to ubams ----#

mkdir $BASE_DIR/outputs/02_align_reference/02_3_merge/

jobfile5=02_3_merge.tmp # temp file
cat > $jobfile5 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=02_3_merge
#SBATCH --partition=carl.p
#SBATCH --array=1-117
#SBATCH --output=$BASE_DIR/logs/02_3_merge_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/02_3_merge_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=75G
#SBATCH --time=05:00:00

ls -1 $BASE_DIR/outputs/02_align_reference/02_2_samalign/ > $BASE_DIR/outputs/listoffiles/02_3_aligned.fofn

INPUT_ALIGNEMENT=$BASE_DIR/outputs/listoffiles/02_3_aligned.fofn

INPUT_AL=\$(cat \${INPUT_ALIGNEMENT} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \$INPUT_AL

sample1=\${INPUT_AL%.*}
sample=\${sample1%.*}
echo \$sample


gatk --java-options "-Xmx68G" \
    MergeBamAlignment \
    --VALIDATION_STRINGENCY SILENT \
    --EXPECTED_ORIENTATIONS FR \
    --ATTRIBUTES_TO_RETAIN X0 \
    -ALIGNED_BAM=$BASE_DIR/outputs/02_align_reference/02_2_samalign/\${INPUT_AL} \
    -UNMAPPED_BAM=$BASE_DIR/outputs/00_ubams/\${sample}.ubam.bam \
    -OUTPUT=$BASE_DIR/outputs/02_align_reference/02_3_merge/\${sample}.mapped.bam \
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

EOA



  # *** MARK DUPLICATES *** #

mkdir $BASE_DIR/outputs/03_mark_duplicates/

#---- Sort the previously created sam files ----#

mkdir $BASE_DIR/outputs/03_mark_duplicates/03_1_sort_sam/

jobfile6=03_1_sortsam.tmp # temp file
cat > $jobfile6 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=03_1_sortsam
#SBATCH --partition=carl.p
#SBATCH --array=1-117
#SBATCH --output=$BASE_DIR/logs/03_1_sortsam_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/03_1_sortsam_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=2-00:00:00

ls -1 $BASE_DIR/outputs/02_align_reference/02_3_merge/ > $BASE_DIR/outputs/listoffiles/03_1_map_merge_bams.fofn

INPUT_MAPPED_BAM=$BASE_DIR/outputs/listoffiles/03_1_map_merge_bams.fofn

MAPPED_BAM=\$(cat \${INPUT_MAPPED_BAM} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \$MAPPED_BAM

sample1=\${MAPPED_BAM%.*}
sample=\${sample1%.*}
echo \$sample

gatk --java-options "-Xmx30G" \
		SortSam \
		-I=$BASE_DIR/outputs/02_align_reference/02_3_merge/\${MAPPED_BAM} \
		-O=$BASE_DIR/outputs/03_mark_duplicates/03_1_sort_sam/\${sample}.sorted.sam \
		--SORT_ORDER="coordinate" \
		--CREATE_INDEX=false \
		--CREATE_MD5_FILE=false \
		-TMP_DIR=$BASE_DIR/temp_files

EOA


#---- Process that put tags? ----#

mkdir $BASE_DIR/outputs/03_mark_duplicates/03_2_tags_intermediate/

jobfile7=03_2_tagsinter.tmp # temp file
cat > $jobfile7 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=03_2_tagsinter
#SBATCH --partition=carl.p
#SBATCH --array=1-117
#SBATCH --output=$BASE_DIR/logs/03_2_tagsinter_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/03_2_tagsinter_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=10:00:00

ls -1 $BASE_DIR/outputs/03_mark_duplicates/03_1_sort_sam/ > $BASE_DIR/outputs/listoffiles/03_2_sort_sam.fofn

INPUT_SORTSAM=$BASE_DIR/outputs/listoffiles/03_2_sort_sam.fofn

SORTSAM=\$(cat \${INPUT_SORTSAM} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \$SORTSAM

sample1=\${SORTSAM%.*}
sample=\${sample1%.*}
echo \$sample

gatk --java-options "-Xmx30G" \
  SetNmAndUqTags \
  --INPUT=$BASE_DIR/outputs/03_mark_duplicates/03_1_sort_sam/\${SORTSAM} \
  --OUTPUT=$BASE_DIR/outputs/03_mark_duplicates/03_2_tags_intermediate/\${sample}.intermediate.bam \
  --CREATE_INDEX=false \
  --CREATE_MD5_FILE=false \
  -TMP_DIR=$BASE_DIR/temp_files \
  --REFERENCE_SEQUENCE=$BASE_DIR/ressources/HP_genome_unmasked_01.fa.gz

EOA


#---- Mark duplicates ----#

mkdir $BASE_DIR/outputs/03_mark_duplicates/03_3_mark_duplicates/
mkdir $BASE_DIR/outputs/03_mark_duplicates/03_3_mark_duplicates/duplicates/
mkdir $BASE_DIR/outputs/03_mark_duplicates/03_3_mark_duplicates/metrics/

jobfile8=03_3_markdups.tmp # temp file
cat > $jobfile8 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=03_3_markdups
#SBATCH --partition=carl.p
#SBATCH --array=1-117
#SBATCH --output=$BASE_DIR/logs/03_3_markdups_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/03_3_markdups_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=05:00:00

ls -1 $BASE_DIR/outputs/03_mark_duplicates/03_2_tags_intermediate/*.bam |xargs -n1 basename > $BASE_DIR/outputs/listoffiles/03_3_tags_intermediate.fofn

INPUT_TAGSINTER=$BASE_DIR/outputs/listoffiles/03_3_tags_intermediate.fofn

TAGSINTER=\$(cat \${INPUT_TAGSINTER} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \$TAGSINTER

sample1=\${TAGSINTER%.*}
sample=\${sample1%.*}
echo \$sample

gatk --java-options "-Xmx30G" \
  MarkDuplicates \
  -I=$BASE_DIR/outputs/03_mark_duplicates/03_2_tags_intermediate/\${TAGSINTER} \
  -O=$BASE_DIR/outputs/03_mark_duplicates/03_3_mark_duplicates/duplicates/\${sample}.dedup.bam \
  -M=$BASE_DIR/outputs/03_mark_duplicates/03_3_mark_duplicates/metrics/\${sample}.dedup.metrics.txt \
  -MAX_FILE_HANDLES=1000  \
  -TMP_DIR=$BASE_DIR/temp_files

EOA



  # *** CREATE INDEX FILE FOR THE ONES WITH DUPLICATES MARKED *** #

mkdir $BASE_DIR/outputs/04_index/

jobfile9=04_index.tmp # temp file
cat > $jobfile9 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=04_index
#SBATCH --partition=carl.p
#SBATCH --array=1-117
#SBATCH --output=$BASE_DIR/logs/04_index_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/04_index_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=01:00:00

ls -1 $BASE_DIR/outputs/03_mark_duplicates/03_3_mark_duplicates/duplicates/ > $BASE_DIR/outputs/listoffiles/04_duplicates.fofn

INPUT_DUPLI=$BASE_DIR/outputs/listoffiles/04_duplicates.fofn

DUPLI=\$(cat \${INPUT_DUPLI} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \$DUPLI

sample1=\${DUPLI%.*}
sample=\${sample1%.*}
echo \$sample

gatk --java-options "-Xmx30G" \
  BuildBamIndex \
  -INPUT=$BASE_DIR/outputs/03_mark_duplicates/03_3_mark_duplicates/duplicates/\${DUPLI}

EOA



  # *** CREATE COVERAGE STATISTIC FILE FOR EACH SAMPLE *** #

mkdir $BASE_DIR/outputs/a_coverage/

jobfile10=a_coverage.tmp # temp file
cat > $jobfile10 <<EOA # generate the job file
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

ls -1 $BASE_DIR/outputs/03_mark_duplicates/03_3_mark_duplicates/duplicates/*.bam|xargs -n1 basename > $BASE_DIR/outputs/listoffiles/a_duplicates_coverage.fofn

INPUT_COV=$BASE_DIR/outputs/listoffiles/a_duplicates_coverage.fofn

COV=\$(cat \${INPUT_COV} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \$COV


sample1=\${COV%.*}
sample=\${sample1%.*}
echo \$sample

gatk --java-options "-Xmx30G" \
  CollectWgsMetrics \
  -I=$BASE_DIR/outputs/03_mark_duplicates/mark_duplicates/duplicates/\${COV} \
  -O=$BASE_DIR/outputs/a_coverage/\${sample}.wgsmetrics.txt \
  -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa.gz

EOA



jobfile101=a_coverage_table.tmp # temp file
cat > $jobfile101 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=a_coverage_table
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/a_coverage_table_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/a_coverage_table_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=00:15:00

FILES=/user/doau0129/work/chapter1_2/outputs/a_coverage/*

for f in $FILES
do
  echo $f
  name_file=${f##*/}
  echo $name_file
  sample1=${name_file%.*}
  sample2=${sample1%.*}
  sample=${sample2%.*}
  echo $sample
  cov=$(cat $f | awk 'FNR == 8 {print $2}')
  echo $cov

  echo "$sample $cov" >> $BASE_DIR/outputs/a_coverage/coverage_table
done

EOA




#Remove all the files with a coverage lower than x15 from the analysis
mkdir $BASE_DIR/outputs/b_removed_from_analysis/


# ---------------- VARIANT CALLING ---------------- #



  # *** CALCULATE GENOTYPE LIKELYHOODS *** #

mkdir $BASE_DIR/outputs/05_genlikely/

jobfile11=05_genlikely.tmp # temp file
cat > $jobfile11 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=05_genlikely
#SBATCH --partition=carl.p
#SBATCH --array=1-113
#SBATCH --output=$BASE_DIR/logs/05_genlikely_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/05_genlikely_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=3-00:00:00

mv $BASE_DIR/outputs/03_mark_duplicates/mark_duplicates/duplicates/AG9RX_47pueboc.1.dedup.bam $BASE_DIR/outputs/b_removed_from_analysis/
mv $BASE_DIR/outputs/03_mark_duplicates/mark_duplicates/duplicates/AG9RX_47pueboc.1.dedup.bai $BASE_DIR/outputs/b_removed_from_analysis/
mv $BASE_DIR/outputs/03_mark_duplicates/mark_duplicates/duplicates/PL17_79abepue.3.dedup.bam $BASE_DIR/outputs/b_removed_from_analysis/
mv $BASE_DIR/outputs/03_mark_duplicates/mark_duplicates/duplicates/PL17_79abepue.3.dedup.bai $BASE_DIR/outputs/b_removed_from_analysis/
mv $BASE_DIR/outputs/03_mark_duplicates/mark_duplicates/duplicates/PL17_101maybel.1.dedup.bam $BASE_DIR/outputs/b_removed_from_analysis/
mv $BASE_DIR/outputs/03_mark_duplicates/mark_duplicates/duplicates/PL17_101maybel.1.dedup.bai $BASE_DIR/outputs/b_removed_from_analysis/
mv $BASE_DIR/outputs/03_mark_duplicates/mark_duplicates/duplicates/PL17_97indbel.3.dedup.bam $BASE_DIR/outputs/b_removed_from_analysis/
mv $BASE_DIR/outputs/03_mark_duplicates/mark_duplicates/duplicates/PL17_97indbel.3.dedup.bai $BASE_DIR/outputs/b_removed_from_analysis/


ls -1 $BASE_DIR/outputs/03_mark_duplicates/mark_duplicates/duplicates/*.bam |xargs -n1 basename > $BASE_DIR/outputs/listoffiles/05_duplicates_hapcaller.fofn

INPUT_DUPLI=$BASE_DIR/outputs/listoffiles/05_duplicates_hapcaller.fofn

DUPLI=\$(cat \${INPUT_DUPLI} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \$DUPLI

sample1=\${DUPLI%.*}
sample=\${sample1%.*}
echo \$sample

gatk --java-options "-Xmx35g" HaplotypeCaller  \
     -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
     -I=$BASE_DIR/outputs/03_mark_duplicates/mark_duplicates/duplicates/\${DUPLI} \
     -O=$BASE_DIR/outputs/05_genlikely/\${sample}.g.vcf.gz \
     -ERC GVCF

EOA


#
#   # *** COMBINE gvcf FILES - ALL SAMPLES TOGETHER *** #
#

mkdir $BASE_DIR/outputs/06_cohort_genotyping/

jobfile12=06_combine_gvcf.tmp # temp file
cat > $jobfile12 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=06_combine
#SBATCH --partition=carl.p
#SBATCH --output=/user/doau0129/work/chapter1_2/logs/06_combine_%A_%a.out
#SBATCH --error=/user/doau0129/work/chapter1_2/logs/06_combine_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --time=4-00:00:00

ls -1 $BASE_DIR/outputs/05_genlikely/*.g.vcf.gz > $BASE_DIR/outputs/listoffiles/genlikely_individual0.fofn

awk '{print "-V", $1}' $BASE_DIR/outputs/listoffiles/genlikely_individual0.fofn > $BASE_DIR/outputs/listoffiles/genlikely_individual.fofn

rm $BASE_DIR/outputs/listoffiles/genlikely_individual0.fofn

INPUT_GEN=$(cat $BASE_DIR/outputs/listoffiles/genlikely_individual.fofn)
echo \${INPUT_GEN}

gatk --java-options "-Xmx85g" \
      CombineGVCFs \
      -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
      \${INPUT_GEN} \
      -O=$BASE_DIR/outputs/06_cohort_genotyping/cohort.g.vcf.gz


EOA
#
#
#
#   # *** ALL SAMPLES ARE JOINTLY GENOTYPED *** #
#
# #GenotypeGVCFs#

mkdir $BASE_DIR/outputs/07_genotyping/
mkdir $BASE_DIR/outputs/07_1_raw_snp/

jobfile13=07_1_genotype.tmp # temp file
cat > $jobfile13 <<EOA # generate the job file
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

EOA
#
#
# #SelectVariants#
#

mkdir $BASE_DIR/outputs/07_genotyping/07_2_all_sites/

jobfile14=07_2_all.tmp # temp file
cat > $jobfile14 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=07_2_all
#SBATCH --partition=carl.p
#SBATCH --array=1-24
#SBATCH --output=/user/doau0129/work/chapter1_2/logs/07_2_all_%A_%a.out
#SBATCH --error=/user/doau0129/work/chapter1_2/logs/07_2_all_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90G
#SBATCH --time=4-00:00:00

gatk --java-options "-Xmx85g" \
    GenotypeGVCFs \
    -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
    -L=LG${SLURM_ARRAY_TASK_ID} \
    -V=$BASE_DIR/outputs/06_cohort_genotyping/cohort.g.vcf.gz  \
    -O=$BASE_DIR/07_genotyping/07_2_all_sites/intermediate.vcf.gz \
    --include-non-variant-sites=true

gatk --java-options "-Xmx85G" \
    SelectVariants \
    -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
    -V=$BASE_DIR/07_genotyping/07_2_all_sites/intermediate.vcf.gz \
    --select-type-to-exclude=INDEL \
    -O=$BASE_DIR/07_genotyping/07_2_all_sites/all_sites.LG${SLURM_ARRAY_TASK_ID}.vcf.gz

rm $BASE_DIR/07_genotyping/07_2_all_sites/intermediate.*

EOA
#
#
#
#   # *** COLLECT THE METRICS OF RAW GENOTYPES IN TABLE *** #
#

mkdir $BASE_DIR/outputs/08_1_variants_metrics/

jobfile15=08_1_snp.tmp # temp file
cat > $jobfile15 <<EOA # generate the job file
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

gatk --java-options "-Xmx25G" \
       VariantsToTable \
       --variant=$BASE_DIR/outputs/07_genotyping/07_1_raw_snp/raw_var_sites.vcf.gz \
       --output=$BASE_DIR/outputs/08_1_variants_metrics/raw_var_sites.table.txt \
       -F=CHROM -F=POS -F=MQ \
       -F=QD -F=FS -F=MQRankSum -F=ReadPosRankSum \
       --show-filtered

EOA
#


mkdir $BASE_DIR/outputs/08_2_merge/

jobfile16=08_2_all.tmp # temp file
cat > $jobfile16 <<EOA # generate the job file
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

ls -1 $BASE_DIR/outputs/07_2_all/all_sites.* > $BASE_DIR/outputs/listoffiles/all_sites_LG0.fofn

awk '{print "-I", $1}' $BASE_DIR/outputs/listoffiles/all_sites_LG0.fofn > $BASE_DIR/outputs/listoffiles/all_sites_LG.fofn

rm $BASE_DIR/outputs/listoffiles/all_sites_LG0.fofn

INPUT_GEN=$(cat $BASE_DIR/outputs/listoffiles/all_sites_LG.fofn)
echo \${INPUT_GEN}

gatk --java-options "-Xmx85g" \
       GatherVcfs \
       \${INPUT_GEN} \
       -O=$BASE_DIR/outputs/08_2_merge/all_sites.vcf.gz

EOA





#
#   # *** TAGGING & FILTERING OF THE GENOTYPES BASED ON METRICS TABLE *** #
#
# #VariantFiltration#
#

mkdir $BASE_DIR/outputs/09_1_snpfiltration/

jobfile17=09_1_filtration.tmp # temp file
cat > $jobfile17 <<EOA # generate the job file
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

EOA
#
#
# #SelectVariants again...#
#

mkdir $BASE_DIR/outputs/09_2_allfiltration/

jobfile18=09_2_filtration.tmp # temp file
cat > $jobfile18 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=09_2_filtration
#SBATCH --partition=carl.p
#SBATCH --output=/user/doau0129/work/chapter1_2/logs/09_2_filtration_%A_%a.out
#SBATCH --error=/user/doau0129/work/chapter1_2/logs/09_2_filtration_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=150G
#SBATCH --time=30:00:00

tabix -p vcf $BASE_DIR/outputs/08_2_merge/all_sites.vcf.gz

gatk --java-options "-Xmx75G" \
       VariantFiltration \
       -R=$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
       -V=$BASE_DIR/outputs/08_2_merge/all_sites.vcf.gz \
       -O=$BASE_DIR/outputs/09_2_allfiltration/intermediate.vcf.gz \
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
       -V=$BASE_DIR/outputs/09_2_allfiltration/intermediate.vcf.gz \
       -O=$BASE_DIR/outputs/09_2_allfiltration/intermediate.filterd.vcf.gz \
       --exclude-filtered \
       --QUIET true \
       --verbosity ERROR  &> var_select.log

vcftools \
       --gzvcf $BASE_DIR/outputs/09_2_allfiltration/intermediate.filterd.vcf.gz \
       --max-missing-count 17 \
       --stdout  \
       --recode | \
       bgzip > $BASE_DIR/outputs/09_2_allfiltration/filterd.allBP.vcf.gz

rm $BASE_DIR/outputs/09_2_allfiltration/intermediate.*


EOA


if [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ] || [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jid15" ] || [ "$JID_RES" = "jid16" ] || [ "$JID_RES" = "jid17" ] || [ "$JID_RES" = "jid18" ];
then
  echo "00_split_samples DONE                   **"
else
  jid1=$(sbatch ${jobfile1})
fi



if [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ] || [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jid15" ] || [ "$JID_RES" = "jid16" ] || [ "$JID_RES" = "jid17" ] || [ "$JID_RES" = "jid18" ];
then
  echo "*****     01_mark_adapters DONE          **"
elif [ "$JID_RES" = jid2 ]
then
  jid2=$(sbatch ${jobfile2})
else
  jid2=$(sbatch --dependency=afterok:${jid1##* } ${jobfile2})
fi



if [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ] || [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jid15" ] || [ "$JID_RES" = "jid16" ] || [ "$JID_RES" = "jid17" ] || [ "$JID_RES" = "jid18" ];
then
  echo "*****     02_1_samtofq DONE              **"
elif [ "$JID_RES" = "jid3" ]
then
  jid3=$(sbatch ${jobfile3})
else
  jid3=$(sbatch --dependency=afterok:${jid2##* } ${jobfile3})
fi



if [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ] || [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jid15" ] || [ "$JID_RES" = "jid16" ] || [ "$JID_RES" = "jid17" ] || [ "$JID_RES" = "jid18" ];
then
  echo "*****     02_2_samalign DONE             **"
elif [ "$JID_RES" = "jid4" ]
then
  jid4=$(sbatch ${jobfile4})
else
  jid4=$(sbatch --dependency=afterok:${jid3##* } ${jobfile4})
fi



if [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ] || [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jid15" ] || [ "$JID_RES" = "jid16" ] || [ "$JID_RES" = "jid17" ] || [ "$JID_RES" = "jid18" ];
then
  echo "*****     02_3_merge DONE                **"
elif [ "$JID_RES" = "jid5" ]
then
  jid5=$(sbatch ${jobfile5})
else
  jid5=$(sbatch --dependency=afterok:${jid4##* } ${jobfile5})
fi



if [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ] || [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jid15" ] || [ "$JID_RES" = "jid16" ] || [ "$JID_RES" = "jid17" ] || [ "$JID_RES" = "jid18" ];
then
  echo "*****     03_1_sortsam DONE              **"
elif [ "$JID_RES" = "jid6" ]
then
  jid6=$(sbatch ${jobfile6})
else
  jid6=$(sbatch --dependency=afterok:${jid5##* } ${jobfile6})
fi



if [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ] || [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jid15" ] || [ "$JID_RES" = "jid16" ] || [ "$JID_RES" = "jid17" ] || [ "$JID_RES" = "jid18" ];
then
  echo "*****     03_2_tagsinter DONE            **"
elif [ "$JID_RES" = "jid7" ]
then
  jid7=$(sbatch ${jobfile7})
else
  jid7=$(sbatch --dependency=afterok:${jid6##* } ${jobfile7})
fi



if [ "$JID_RES" = "jid9" ] || [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jid15" ] || [ "$JID_RES" = "jid16" ] || [ "$JID_RES" = "jid17" ] || [ "$JID_RES" = "jid18" ];
then
  echo "*****     03_3_markdups DONE             **"
elif [ "$JID_RES" = "jid8" ]
then
  jid8=$(sbatch ${jobfile8})
else
  jid8=$(sbatch --dependency=afterok:${jid7##* } ${jobfile8})
fi



if [ "$JID_RES" = "jid10" ] || [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jid15" ] || [ "$JID_RES" = "jid16" ] || [ "$JID_RES" = "jid17" ] || [ "$JID_RES" = "jid18" ];
then
  echo "*****     04_index DONE                  **"
elif [ "$JID_RES" = "jid9" ]
then
  jid9=$(sbatch ${jobfile9})
else
  jid9=$(sbatch --dependency=afterok:${jid8##* } ${jobfile9})
fi

if [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jid15" ] || [ "$JID_RES" = "jid16" ] || [ "$JID_RES" = "jid17" ] || [ "$JID_RES" = "jid18" ];
then
  echo "*****     a_coverage DONE                **"
elif [ "$JID_RES" = "jid10" ]
then
  jid10=$(sbatch ${jobfile10})
else
  jid10=$(sbatch --dependency=afterok:${jid9##* } ${jobfile10})
fi

jid101=$(sbatch --dependency=afterok:${jid10##* } ${jobfile101})

if [ "$JID_RES" = "jid12" ] || [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jid15" ] || [ "$JID_RES" = "jid16" ] || [ "$JID_RES" = "jid17" ] || [ "$JID_RES" = "jid18" ];
then
  echo "*****     05_genlikely DONE              **"
elif [ "$JID_RES" = "jid11" ]
then
  jid11=$(sbatch ${jobfile11})
else
  jid11=$(sbatch --dependency=afterok:${jid10##* } ${jobfile11})
fi



# if [ "$JID_RES" = "jid13" ] || [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jid15" ] || [ "$JID_RES" = "jid16" ] || [ "$JID_RES" = "jid17" ] || [ "$JID_RES" = "jid18" ];
# then
#   echo "*****     06_combine_gvcf DONE           **"
# elif [ "$JID_RES" = "jid12" ]
# then
#   jid12=$(sbatch ${jobfile12})
# else
#   jid12=$(sbatch --dependency=afterok:${jid11##* } ${jobfile12})
# fi
#
#
#
# if [ "$JID_RES" = "jid14" ] || [ "$JID_RES" = "jid15" ] || [ "$JID_RES" = "jid16" ] || [ "$JID_RES" = "jid17" ] || [ "$JID_RES" = "jid18" ];
# then
#   echo "*****     07_1_genotype DONE             **"
# elif [ "$JID_RES" = "jid13" ]
# then
#   jid13=$(sbatch ${jobfile13})
# else
#   jid13=$(sbatch --dependency=afterok:${jid12##* } ${jobfile13})
# fi
#
#
#
# if [ "$JID_RES" = "jid15" ] || [ "$JID_RES" = "jid16" ] || [ "$JID_RES" = "jid17" ] || [ "$JID_RES" = "jid18" ];
# then
#   echo "*****     07_2_select_variants DONE      **"
# elif [ "$JID_RES" = "jid14" ]
# then
#   jid14=$(sbatch ${jobfile14})
# else
#   jid14=$(sbatch --dependency=afterok:${jid13##* } ${jobfile14})
# fi
#
#
#
# if [ "$JID_RES" = "jid16" ] || [ "$JID_RES" = "jid17" ] || [ "$JID_RES" = "jid18" ];
# then
#   echo "*****     08_variant_table DONE          **"
# elif [ "$JID_RES" = "jid15" ]
# then
#   jid15=$(sbatch ${jobfile15})
# else
#   jid15=$(sbatch --dependency=afterok:${jid14##* } ${jobfile15})
# fi
#
#
#
# if [ "$JID_RES" = "jid17" ] || [ "$JID_RES" = "jid18" ];
# then
#   echo "*****     09_1_filter_variant DONE       **"
# elif [ "$JID_RES" = "jid16" ]
# then
#   jid16=$(sbatch ${jobfile16})
# else
#   jid16=$(sbatch --dependency=afterok:${jid15##* } ${jobfile16})
# fi
#
#
#
# if [ "$JID_RES" = "jid18" ];
# then
#   echo "*****     09_2_select_variant DONE       **"
# elif [ "$JID_RES" = "jid17" ]
# then
  jid17=$(sbatch ${jobfile17})
else
   jid17=$(sbatch --dependency=afterok:${jid16##* } ${jobfile17})
fi



jid18=$(sbatch --dependency=afterok:${jid17##* } ${jobfile18})
