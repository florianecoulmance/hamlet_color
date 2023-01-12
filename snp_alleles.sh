#!/bin/bash
# by: Floriane Coulmance: 17/10/2022
# usage:
# sbatch fst.sh -i <PATH> -j <JOB_ID> 
# ------------------------------------------------------------------------------
# PATH corresponds to the path to the base directory, all outputs and necessary
# folder will be created by the script
# JOB_ID corresponds to ids of jobs from where you want the analysis to be ran
# within the pipeline
# ------------------------------------------------------------------------------



# ********** Allow to enter bash options **********
# -------------------------------------------------

while getopts i:j:k: option
do
case "${option}"
in
i) BASE_DIR=${OPTARG};; # get the base directory path
j) JID_RES=${OPTARG};; # get the jobid from which you want to resume
esac
done



# ********* Create necessary repositories *********
# -------------------------------------------------

# Repo for gxp outputs
# mkdir $BASE_DIR/outputs/8_fst/



# ********* Jobs creation *************************
# -------------------------------------------------

# ------------------------------------------------------------------------------
# Job 0 prepare snp position table 

jobfile0=0_table.tmp # generate a temp file that will be launched
cat > $jobfile0 <<EOA # indicate that EOA is the end of the file
#!/bin/bash
#SBATCH --job-name=0_keep_species                                                               # set the jobname
#SBATCH --partition=carl.p                                                                      # set the cluster partition to use
#SBATCH --output=$BASE_DIR/logs/0_table_%A_%a.out                                        # send the job output file to the log folder
#SBATCH --error=$BASE_DIR/logs/0_table_%A_%a.err                                         # send the job error file to the log folder
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G                                                                       # set the estimated memory needed for the job to run
#SBATCH --time=02:30:00                                                                         # set the estimated amount of time for the job to run


INPUT=$BASE_DIR/outputs/7_gxp/continuous/LAB/LAB_fullm_54off_59on/PC1_5/*.txt                                # input the genotyping file with SNPs only created with the genotyping.sh pipeline
echo \${INPUT}

awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' \${INPUT} > $BASE_DIR/outputs/7_gxp/continuous/LAB/LAB_fullm_54off_59on/PC1_5/all.txt 

awk '!seen[\$0]++' $BASE_DIR/outputs/7_gxp/continuous/LAB/LAB_fullm_54off_59on/PC1_5/all.txt > $BASE_DIR/outputs/7_gxp/continuous/LAB/LAB_fullm_54off_59on/PC1_5/PC1_5_snp_all.txt


EOA



# ------------------------------------------------------------------------------
# Job 1 prepare the snp file with allelle information

jobfile1=1_allele.tmp # generate a temp file that will be launched
cat > $jobfile1 <<EOA # indicate that EOA is the end of the file
#!/bin/bash
#SBATCH --job-name=1_allele                                                               # set the jobname
#SBATCH --partition=carl.p                                                                      # set the cluster partition to use
#SBATCH --output=$BASE_DIR/logs/1_allele_%A_%a.out                                        # send the job output file to the log folder
#SBATCH --error=$BASE_DIR/logs/1_allele_%A_%a.err                                         # send the job error file to the log folder
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G                                                                       # set the estimated memory needed for the job to run
#SBATCH --time=02:30:00                                                                         # set the estimated amount of time for the job to run


INPUT=$BASE_DIR/outputs/7_gxp/continuous/LAB/LAB_fullm_54off_59on/PC1_5/PC1_5_snp_all.txt

sed 1d \${INPUT} | while read -r line;
    do
        echo \${line}
        chrom=\$(echo \${line} | awk '{print \$1}')
        echo \${chrom}
        pos=\$(echo \${line} | awk '{print \$2}')
        echo \${pos}
        bcftools filter -r \${chrom}:\${pos} $BASE_DIR/outputs/6_genotyping/6_1_snp/snp_filterd.vcf.gz | grep -v '##' >> $BASE_DIR/outputs/7_gxp/continuous/LAB/LAB_fullm_54off_59on/PC1_5/intermediate.txt
done

awk '!seen[\$0]++' $BASE_DIR/outputs/7_gxp/continuous/LAB/LAB_fullm_54off_59on/PC1_5/intermediate.txt | sed 's|#||g' > $BASE_DIR/outputs/7_gxp/continuous/LAB/LAB_fullm_54off_59on/PC1_5/PC1_5_alleles.txt


EOA



# ------------------------------------------------------------------------------
# Job 2 plot pie chart for each SNP

jobfile2=2_pie.tmp # generate a temp file that will be launched
cat > $jobfile2 <<EOA # indicate that EOA is the end of the file
#!/bin/bash
#SBATCH --job-name=2_pie                                                               # set the jobname
#SBATCH --partition=carl.p                                                                      # set the cluster partition to use
#SBATCH --output=$BASE_DIR/logs/2_pie_%A_%a.out                                        # send the job output file to the log folder
#SBATCH --error=$BASE_DIR/logs/2_pie_%A_%a.err                                         # send the job error file to the log folder
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G                                                                       # set the estimated memory needed for the job to run
#SBATCH --time=02:30:00                                                                         # set the estimated amount of time for the job to run


ml hpc-env/8.3
ml R-bundle-Bioconductor/3.12-foss-2019b-R-4.0.2
ml FriBidi
ml HarfBuzz

INPUT=$BASE_DIR/outputs/7_gxp/continuous/LAB/LAB_fullm_54off_59on/PC1_5/PC1_5_alleles.txt

Rscript $BASE_DIR/R/snp_alleles_plots.R \${INPUT} $BASE_DIR/figures/7_gxp/continuous/LAB/LAB_fullm_54off_59on/PC1_5/

EOA


# ------------------------------------------------------------------------------
# Job 3 

jobfile3=3_snp_anno.tmp # generate a temp file that will be launched
cat > $jobfile3 <<EOA # indicate that EOA is the end of the file
#!/bin/bash
#SBATCH --job-name=3_snp_anno                                                               # set the jobname
#SBATCH --partition=carl.p                                                                      # set the cluster partition to use
#SBATCH --array=0-10
#SBATCH --output=$BASE_DIR/logs/3_snp_anno_%A_%a.out                                        # send the job output file to the log folder
#SBATCH --error=$BASE_DIR/logs/3_snp_anno_%A_%a.err                                         # send the job error file to the log folder
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G                                                                       # set the estimated memory needed for the job to run
#SBATCH --time=02:30:00                                                                         # set the estimated amount of time for the job to run


chr=(03 04 05 06 07 08 09 12 19 20 23)
LG=\${chr[\${SLURM_ARRAY_TASK_ID}]}
echo \${LG}
CHROM="LG\${LG}"
echo \${CHROM}

INPUT=$BASE_DIR/outputs/7_gxp/continuous/LAB/LAB_fullm_54off_59on/PC1_5/PC1_5_\${CHROM}.snp.txt


awk 'BEGIN { OFS="\t" } {print \$6, \$1, \$2, \$2, \$3="+"}' \${INPUT} | \
sed 's/RANGE\tCHROM\tPOS\tPOS\t+/Unique Peak ID\tchromosome\tstarting position\tending position\tStrand/g' > $BASE_DIR/outputs/7_gxp/continuous/LAB/LAB_fullm_54off_59on/PC1_5/PC1_5_\${CHROM}_homer.txt

/user/doau0129/miniconda3/bin/annotatePeaks.pl $BASE_DIR/outputs/7_gxp/continuous/LAB/LAB_fullm_54off_59on/PC1_5/PC1_5_\${CHROM}_homer.txt /user/doau0129/data/ref_genome/HP_genome_unmasked_01.fa.gz -gtf /user/doau0129/data/annotations/HP.annotation.named.\${CHROM}.gff.gz -annStats $BASE_DIR/outputs/7_gxp/continuous/LAB/LAB_fullm_54off_59on/PC1_5/PC1_5_\${CHROM}_homer_output_annStats.txt > $BASE_DIR/outputs/7_gxp/continuous/LAB/LAB_fullm_54off_59on/PC1_5/PC1_5_\${CHROM}_homer_output.txt



EOA





# ********** Schedule the job launching ***********
# -------------------------------------------------

if [ "$JID_RES" = "jid1" ] || [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ];
then
  echo "*****   0_table    : DONE         **"
else
  jid0=$(sbatch ${jobfile0})
fi


if [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ];
then
  echo "*****   1_alleles  : DONE         **"
elif [ "$JID_RES" = jid1 ]
then
  jid1=$(sbatch ${jobfile1})
else
  jid1=$(sbatch --dependency=afterok:${jid0##* } ${jobfile1})
fi

if [ "$JID_RES" = "jid3" ];
then
  echo "*****   2_pie      : DONE         **"
elif [ "$JID_RES" = jid2 ]
then
  jid2=$(sbatch ${jobfile2})
else
  jid2=$(sbatch --dependency=afterok:${jid1##* } ${jobfile2})
fi

if [ "$JID_RES" = "jid3" ];
then
  jid3=$(sbatch ${jobfile3})
else
  jid3=$(sbatch --dependency=afterok:${jid2##* } ${jobfile3})
fi






# ******** Remove useless files and folders *******
# -------------------------------------------------

rm *tmp
