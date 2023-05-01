#!/bin/bash
# by: Floriane Coulmance: 06/04/2021
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
mkdir $BASE_DIR/outputs/8_fst/

# Repo for gxp figures
mkdir $BASE_DIR/figures/8_fst/



# ********* Jobs creation *************************
# -------------------------------------------------

# ------------------------------------------------------------------------------
# Job 0 generates pairwise comparison labels for population with 7 or more individuals

jobfile0=0_pairs.tmp # temp file
cat > $jobfile0 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=0_pairs
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/0_pairs_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/0_pairs_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=04:30:00


INPUT=$BASE_DIR/outputs/6_genotyping/6_1_snp/snp_filterd.vcf.gz

for k in puebel nigbel maybel indbel pueboc nigboc uniflo unipue puepue chlpue;
do
  echo \$k 
  vcfsamplenames \${INPUT} | grep \${k} > $BASE_DIR/outputs/8_fst/\$k.pop

  for j in puebel nigbel maybel indbel pueboc nigboc uniflo unipue puepue chlpue;
  do
    echo \$k \$j \$(tr ' ' '\n' <<<"\$k \$j" | sort | tr -d '\n') >> $BASE_DIR/outputs/8_fst/pairwise_comparisons.txt;
  done;
done

awk '!seen[\$3]++' $BASE_DIR/outputs/8_fst/pairwise_comparisons.txt | awk '{print \$1, \$2}' > $BASE_DIR/outputs/8_fst/pairwise_comparison1.txt


EOA



# ------------------------------------------------------------------------------
# Job 1 creates 1 file per population with individuals' labels and calculates 
# pairwise Fst between 2 populations

jobfile1=1_fst.tmp # temp file
cat > $jobfile1 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=1_fst
#SBATCH --partition=carl.p
#SBATCH --array=1-55
#SBATCH --output=$BASE_DIR/logs/1_fst_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/1_fst_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=04:30:00


INPUT=$BASE_DIR/outputs/8_fst/pairwise_comparison1.txt
PAIR=\$(cat \${INPUT} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)                               # create 1 job per combination (job array)
echo \${PAIR}                                                                                      # print the combination of this job

POP1=\$(cut -d' ' -f1 <<< \$PAIR)                                                                  # find the first population
POP2=\$(cut -d' ' -f2 <<< \$PAIR)                                                                  # find the 2nd population
echo \${POP1}
echo \${POP2}

VCF=$BASE_DIR/outputs/6_genotyping/6_1_snp/snp_filterd.vcf.gz                                        # input the SNP genotyping file
echo \${VCF}

FILE1=$BASE_DIR/outputs/8_fst/\${POP1}.pop                                                     # input the 2 corresponding population files
FILE2=$BASE_DIR/outputs/8_fst/\${POP2}.pop                                                     # this will be used as an input for FST calculations

# use VCFTOOLS to calculate pairwise FST
vcftools --gzvcf \${VCF} --weir-fst-pop \${FILE1} --weir-fst-pop \${FILE2} --out $BASE_DIR/outputs/8_fst/\${POP1}_\${POP2}.nowindow 2> $BASE_DIR/outputs/8_fst/\${POP1}_\${POP2}_nowindow.log


EOA



# ------------------------------------------------------------------------------
# Job 2 retrieve weighted Fst from output file for each pairwise Fst comparison 

jobfile2=2_weighted.tmp # temp file
cat > $jobfile2 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=2_weighted
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/2_weighted_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/2_weighted_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=04:30:00


# find lines where mean and weighted FST estimates are reported
zless $BASE_DIR/outputs/8_fst/*_nowindow.log | \
    grep -E 'Weir and Cockerham|--out' | \
    grep -A 8 nowindow | \
    sed '/^--/d; s/^.*--out //g; s/.nowindow//g; /^Output/d; s/Weir and Cockerham //g; s/ Fst estimate: /\t/g' | \
    paste - - - | \
    cut -f 1,3,5 | \
    sed 's/^\\(...\\)-/\\1\\t/g' > $BASE_DIR/outputs/8_fst/fst_globals_pop_nowindow.txt            # create global FST file as a table with entry for each pairwise comparison (24 rows in total) and columns for name of the comparison, mean FST, weighted FST (3 columns) 


EOA



# ------------------------------------------------------------------------------
# Job 3 creates the joint weighted Fst for all population with 7 or more individuals

jobfile3=3_global.tmp # temp file
cat > $jobfile3 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=3_global
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/3_global_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/3_global_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=04:30:00


INPUT=$BASE_DIR/outputs/6_genotyping/6_1_snp/snp_filterd.vcf.gz

for k in abebel abeboc abepue chlpue floflo gemflo gumboc gutpue indbel indpue maybel nigbel nigboc nigflo puebel pueboc pueflo puepue ranbel tanpue uniboc uniflo unipue;
do
  vcfsamplenames \${INPUT} | grep \${k} > $BASE_DIR/outputs/8_fst/\$k.pop
  vcfsamplenames \${INPUT} | grep \${k} >> $BASE_DIR/outputs/8_fst/all_individuals.pop
done 

POP="--weir-fst-pop $BASE_DIR/outputs/8_fst/chlpue.pop \
   --weir-fst-pop $BASE_DIR/outputs/8_fst/indbel.pop \
   --weir-fst-pop $BASE_DIR/outputs/8_fst/maybel.pop \
   --weir-fst-pop $BASE_DIR/outputs/8_fst/nigbel.pop \
   --weir-fst-pop $BASE_DIR/outputs/8_fst/nigboc.pop \
   --weir-fst-pop $BASE_DIR/outputs/8_fst/puebel.pop \
   --weir-fst-pop $BASE_DIR/outputs/8_fst/pueboc.pop \
   --weir-fst-pop $BASE_DIR/outputs/8_fst/puepue.pop \
   --weir-fst-pop $BASE_DIR/outputs/8_fst/uniflo.pop \
   --weir-fst-pop $BASE_DIR/outputs/8_fst/unipue.pop"

echo "\${POP}"

vcftools --gzvcf \${INPUT} \${POP} --out $BASE_DIR/outputs/8_fst/all_pop 2> $BASE_DIR/outputs/8_fst/all_pop.log


# FILE=$BASE_DIR/outputs/8_fst/all_individuals.pop

# # use VCFTOOLS to calculate pairwise FST
# vcftools --gzvcf \${VCF} --weir-fst-pop \${FILE} --weir-fst-pop \${FILE} --out $BASE_DIR/outputs/8_fst/all_individuals 2> $BASE_DIR/outputs/8_fst/all_individuals.log


EOA


# ------------------------------------------------------------------------------
# Job 4 create table of pairwise comparison for manuscript

jobfile4=4_table.tmp # temp file
cat > $jobfile4 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=4_table
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/4_table_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/4_table_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=04:30:00


module load hpc-env/8.3                                                                           # load necessary environment for necessary version of R 
module load R/4.0.2-foss-2019b
ml FriBidi
ml HarfBuzz

INPUT=$BASE_DIR/outputs/8_fst/fst_globals_pop_nowindow.txt  

Rscript $BASE_DIR/R/fst_pairwise_table.R \${INPUT} $BASE_DIR/figures/8_fst/


EOA



# ********** Schedule the job launching ***********
# -------------------------------------------------

if [ "$JID_RES" = "jid1" ] || [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ];
then
  echo "*****   0_pairs       : DONE         **"
else
  jid0=$(sbatch ${jobfile0})
fi


if [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ];
then
  echo "*****   1_fst         : DONE         **"
elif [ "$JID_RES" = jid1 ]
then
  jid1=$(sbatch ${jobfile1})
else
  jid1=$(sbatch --dependency=afterok:${jid0##* } ${jobfile1})
fi


if [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ];
then
  echo "*****   2_weighted    : DONE         **"
elif [ "$JID_RES" = jid2 ]
then
  jid2=$(sbatch ${jobfile2})
else
  jid2=$(sbatch --dependency=afterok:${jid1##* } ${jobfile2})
fi


if [ "$JID_RES" = "jid4" ];
then
  echo "*****   3_global      : DONE         **"
elif [ "$JID_RES" = "jid3" ]
then
  jid3=$(sbatch ${jobfile3})
else
  jid3=$(sbatch --dependency=afterok:${jid2##* } ${jobfile3})
fi


if [ "$JID_RES" = "jid4" ];
then
  jid4=$(sbatch ${jobfile4})
else
  jid4=$(sbatch --dependency=afterok:${jid3##* } ${jobfile4})
fi



# ******** Remove useless files and folders *******
# -------------------------------------------------

rm *tmp
