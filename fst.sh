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



# ********* Jobs creation *************************
# -------------------------------------------------

# --------------------------- PREPARATION -------------------------------------#

# ------------------------------------------------------------------------------
# Job 0 create files with list of individuals for 5 datasets to use for FST

jobfile0=0_keep_species.tmp # generate a temp file that will be launched
cat > $jobfile0 <<EOA # indicate that EOA is the end of the file
#!/bin/bash
#SBATCH --job-name=0_keep_species                                                               # set the jobname
#SBATCH --partition=carl.p                                                                      # set the cluster partition to use
#SBATCH --array=1-5                                                                             # set the array numbers
#SBATCH --output=$BASE_DIR/logs/0_keep_species_%A_%a.out                                        # send the job output file to the log folder
#SBATCH --error=$BASE_DIR/logs/0_keep_species_%A_%a.err                                         # send the job error file to the log folder
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G                                                                       # set the estimated memory needed for the job to run
#SBATCH --time=02:30:00                                                                         # set the estimated amount of time for the job to run


INPUT=$BASE_DIR/outputs/6_genotyping/6_1_snp/snp_filterd.vcf.gz                                 # input the genotyping file with SNPs only created with the genotyping.sh pipeline

tr=(nig pue bel boc flo puer)                                                                   # list of 2 species and 3 locations that have enough data points to be analysed for FST 
printf "%s\n" "\${tr[@]}" > $BASE_DIR/outputs/lof/fst_species.fofn                              # put this list with each entry in a new line in a file

INPUT_SP=$BASE_DIR/outputs/lof/fst_species.fofn                                                 # use the previous file as input and create a job for each line of this file
SP=\$(cat \${INPUT_SP} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \${SP}                                                                                     # print the entry to the output file

if [ \${SP} = "nig" ];                                                                          # find samples to keep in genotyping file
then
  echo \${SP}
  vcfsamplenames \${INPUT} | grep nig > $BASE_DIR/outputs/8_fst/\${SP}.pop             # find samples corresponding to H. nigricans for each location
fi

if [ \${SP} = "pue" ];                                                                          # find samples corresponding to H. puella for each location
then
  echo \${SP}
  vcfsamplenames \${INPUT} | grep pue | grep -v abe | grep -v chl | grep -v gem | grep -v gum | grep -v gut | grep -v ind | grep -v may | grep -v nig | grep -v ran | grep -v tan | grep -v uni > $BASE_DIR/outputs/8_fst/\${SP}.pop
fi

if [ \${SP} = "bel" ];                                                                          # find all species within Belize that have enough entries
then
  echo \${SP}
  vcfsamplenames \${INPUT} | grep bel | grep -v abe | grep -v chl | grep -v flo | grep -v gem | grep -v gum | grep -v gut | grep -v ran | grep -v tan | grep -v uni > $BASE_DIR/outputs/8_fst/\${SP}.pop
fi

if [ \${SP} = "boc" ];                                                                          # find all species within Bocas (Panama) that have enough entries
then
  echo \${SP}
  vcfsamplenames \${INPUT} | grep boc | grep -v abe | grep -v chl | grep -v flo | grep -v gem | grep -v gum | grep -v gut | grep -v ran | grep -v tan | grep -v ind | grep -v may > $BASE_DIR/outputs/8_fst/\${SP}.pop
fi

if [ \${SP} = "puer" ];                                                                         # find all species within Puerto Rico that have enough entries
then
  echo \${SP}
  vcfsamplenames \${INPUT} | grep pue | grep -v abe | grep -v nig | grep -v flo | grep -v gem | grep -v gum | grep -v gut | grep -v ran | grep -v ind | grep -v may | grep -v bel | grep -v boc  > $BASE_DIR/outputs/8_fst/\${SP}.pop
fi

vcftools --gzvcf \${INPUT} \                                                                    # create SNP dataset for each of the 5 datasets above
      --keep $BASE_DIR/outputs/8_fst/\${SP}.pop \                                      # keep samples corresponding to each dataset
      --mac 3 \   
      --recode \
      --stdout | gzip > $BASE_DIR/outputs/8_fst/\${SP}.vcf.gz                          # output as vcf file


EOA



# ------------------------------------------------------------------------------
# Job 1 get samples labels, species abbreviation and location abbreviation for
# each sample into a table, this is a prliminary step

jobfile1=1_all_hamlets.tmp # temp file
cat > $jobfile1 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=1_all_hamlets
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/1_all_hamlets_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/1_all_hamlets_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=02:30:00


INPUT_VCF=$BASE_DIR/outputs/6_genotyping/6_1_snp/snp_filterd.vcf.gz                             # input the SNP dataset 
echo \${INPUT_VCF}

vcfsamplenames \${INPUT_VCF} | \                                                                # get all samples from input file and modify to get a table that is tab separated
       grep "abe\\|gum\\|ind\\|may\\|nig\\|pue\\|ran\\|uni\\|gem\\|flo\\|chl\\|tan\\|gut" | \
       awk '{print \$1"\\t"\$1}' | \
       sed 's/\\t.*\\(...\\)\\(...\\)\$/\\t\\1\\t\\2/g' > $BASE_DIR/outputs/8_fst/hamlets.pop.txt


EOA



# ------------------------------------------------------------------------------
# Job 2 calculates global differentiation (FST) over all the hamlet populations

jobfile2=2_multi.tmp # temp file
cat > $jobfile2 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=2_multi
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/2_multi_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/2_multi_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=15:00:00


INPUT_VCF=$BASE_DIR/outputs/6_genotyping/6_1_snp/snp_filterd.vcf.gz                             # input the SNP genotyping file
echo \${INPUT_VCF}


INPUT_POP=$BASE_DIR/outputs/8_fst/hamlets.pop.txt                                      # input the list of samples table created at the previous step
echo \${INPUT_POP}

awk '{print \$1"\\t"\$2\$3}' \${INPUT_POP} > $BASE_DIR/outputs/8_fst/pop.txt           # from samples file modify and create an intermediate file with hamlet population names

for k in abebel abeboc abepue chlpue floflo gemflo gumboc gutpue indbel indpue maybel nigbel nigboc nigflo nigpue puebel pueboc pueflo puepue ranbel tanpue uniboc uniflo unipue; do
  grep \$k pop.txt | cut -f 1 > $BASE_DIR/outputs/fst/pop.\$k.txt
  done                                                                                          # separate the different populations into different files


POP="--weir-fst-pop $BASE_DIR/outputs/fst/pop.abebel.txt \                                      # create string for vcftools arguments including all the population files
   --weir-fst-pop $BASE_DIR/outputs/fst/pop.abeboc.txt \
   --weir-fst-pop $BASE_DIR/outputs/fst/pop.abepue.txt \
   --weir-fst-pop $BASE_DIR/outputs/fst/pop.chlpue.txt \
   --weir-fst-pop $BASE_DIR/outputs/fst/pop.floflo.txt \
   --weir-fst-pop $BASE_DIR/outputs/fst/pop.gemflo.txt \
   --weir-fst-pop $BASE_DIR/outputs/fst/pop.gumboc.txt \
   --weir-fst-pop $BASE_DIR/outputs/fst/pop.gutpue.txt \
   --weir-fst-pop $BASE_DIR/outputs/fst/pop.indbel.txt \
   --weir-fst-pop $BASE_DIR/outputs/fst/pop.indpue.txt \
   --weir-fst-pop $BASE_DIR/outputs/fst/pop.maybel.txt \
   --weir-fst-pop $BASE_DIR/outputs/fst/pop.nigbel.txt \
   --weir-fst-pop $BASE_DIR/outputs/fst/pop.nigboc.txt \
   --weir-fst-pop $BASE_DIR/outputs/fst/pop.nigflo.txt \
  --weir-fst-pop $BASE_DIR/outputs/fst/pop.nigpue.txt \
  --weir-fst-pop $BASE_DIR/outputs/fst/pop.puebel.txt \
  --weir-fst-pop $BASE_DIR/outputs/fst/pop.pueboc.txt \
  --weir-fst-pop $BASE_DIR/outputs/fst/pop.pueflo.txt \
  --weir-fst-pop $BASE_DIR/outputs/fst/pop.puepue.txt \
  --weir-fst-pop $BASE_DIR/outputs/fst/pop.ranbel.txt \
  --weir-fst-pop $BASE_DIR/outputs/fst/pop.tanpue.txt \
  --weir-fst-pop $BASE_DIR/outputs/fst/pop.uniboc.txt \
  --weir-fst-pop $BASE_DIR/outputs/fst/pop.uniflo.txt \
  --weir-fst-pop $BASE_DIR/outputs/fst/pop.unipue.txt"


  # fst by SNP
     # ----------
vcftools --gzvcf \${INPUT_VCF} \                                                                # use VCFTOOLS to calculate FST by SNP
      \$POP \
     --stdout  2> multi_fst_snp.log | \
     gzip > $BASE_DIR/outputs/8_fst/multi_fst_snp.tsv.gz

     # fst 50kb window
     # ---------------
vcftools --gzvcf \${INPUT_VCF} \                                                               # use VCFTOOLS to calculate FST by 50kb SNP windows of the genome
      \$POP \
     --fst-window-step 5000 \
     --fst-window-size 50000 \
     --stdout  2> multi_fst.50k.log | \
     gzip > $BASE_DIR/outputs/8_fst/multi_fst.50k.tsv.gz

     # fst 10kb window
     # ---------------
vcftools --gzvcf \${INPUT_VCF} \                                                               # use VCFTOOLS to calcultae FST by 10kb SNP windows of the genome
      \$POP \   
     --fst-window-step 1000 \
     --fst-window-size 10000 \
     --stdout  2> multi_fst.10k.log | \
     gzip > $BASE_DIR/outputs/8_fst/multi_fst.10k.tsv.gz


Rscript \$BASE_DIR/R/table_fst_outliers.R multi_fst.50k.tsv.gz $BASE_DIR/outputs/8_fst/ # creates table of outlier positions/regions from FST values _ find original R code here : https://github.com/k-hench/hamlet_radiation/blob/master/R/table_fst_outliers.R, output is fst_outliers_998.tsv


EOA



# ------------------------------------------------------------------------------
# Job 3 create all combinations of pairwise comparisons from the 5 datasets

jobfile3=3_prep_pairwise.tmp # temp file
cat > $jobfile3 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=3_prep_pairwise
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/3_prep_pairwise_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/3_prep_pairwise_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100M
#SBATCH --time=00:00:55


printf 'nig\n%.0s' {1..6} > col1                                                               # all pairwise location combinations for H. nigricans 
col2=(flo flo flo boc boc bel)
printf "%s\n" "\${col2[@]}" > col2
paste -d " " col1 col2 > col12
col3=(boc bel pue bel pue pue)
printf "%s\n" "\${col3[@]}" > col3
paste -d " " col12 col3 > nig

printf 'pue\n%.0s' {1..6} > col11                                                              # all pairwise location combinations for H. puella 
col22=(flo flo flo boc boc bel)
printf "%s\n" "\${col22[@]}" > col22
paste -d " " col11 col22 > col1122
col33=(boc bel pue bel pue pue)
printf "%s\n" "\${col33[@]}" > col33
paste -d " " col1122 col33 > pue


printf 'bel\n%.0s' {1..6} > col111                                                             # all pairwise species combinations for Belize 
col222=(ind ind ind may may nig)
printf "%s\n" "\${col222[@]}" > col222
paste -d " " col111 col222 > col111222
col333=(may nig pue nig pue pue)
printf "%s\n" "\${col333[@]}" > col333
paste -d " " col111222 col333 > bel


printf 'boc\n%.0s' {1..3} > col1111                                                            # all pairwise species combinations for Bocas (Panama)
col2222=(nig nig pue)
printf "%s\n" "\${col2222[@]}" > col2222
paste -d " " col1111 col2222 > col11112222
col3333=(pue uni uni)
printf "%s\n" "\${col3333[@]}" > col3333
paste -d " " col11112222 col3333 > boc


printf 'puer\n%.0s' {1..3} > col11111                                                          # all pairwise species combinations for Puerto Rico
col22222=(chl chl pue)
printf "%s\n" "\${col22222[@]}" > col22222
paste -d " " col11111 col22222 > col1111122222
col33333=(pue uni uni)
printf "%s\n" "\${col33333[@]}" > col33333
paste -d " " col1111122222 col33333 > puer

awk 'NF' nig pue bel boc puer > $BASE_DIR/outputs/lof/fst_pairwise.fofn                        # combine the 5 datasets combinations into 1 file
rm nig                                                                                         # remove all unecessary files and variables
rm pue
rm bel
rm boc
rm puer
rm col*


EOA



# ------------------------------------------------------------------------------
# Job 4 caculate pairwise FST for the 5 datasets (2 species, 3 locations)

jobfile4=4_pairwise_fst.tmp # temp file
cat > $jobfile4 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=4_pairwise_fst
#SBATCH --partition=carl.p
#SBATCH --array=1-24
#SBATCH --output=$BASE_DIR/logs/4_pairwise_fst_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/4_pairwise_fst_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=04:30:00


INPUT_PW=$BASE_DIR/outputs/lof/fst_pairwise.fofn                                               # input the file with lines corresponding to different pairwise combinations for the 5 datasets
PW=\$(cat \${INPUT_PW} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)                          # create 1 job per combination (job array)
echo \${PW}                                                                                    # print the combination of this job

SP=\$(cut -d' ' -f1 <<< \$PW)                                                                  # find the dataset (either a species or a location)
LOC1=\$(cut -d' ' -f2 <<< \$PW)                                                                # find the 1st pairwise component (either a location or a species)
LOC2=\$(cut -d' ' -f3 <<< \$PW)                                                                # find the 2nd pairwise component (either a location or a species)
echo \${SP}                                                                                    # print each 
echo \${LOC1}
echo \${LOC2}

POP=$BASE_DIR/outputs/8_fst/\${SP}.pop                                                # input the list of samples correponding to the dataset (either a species or a location)
echo \${POP}

VCF=$BASE_DIR/outputs/8_fst/\${SP}.vcf.gz                                             # input the SNP genotyping file correponding to the dataset (with only the sample listed in the dataset)
echo \${VCF}


grep \${LOC1} \${POP} > $BASE_DIR/outputs/8_fst/\${SP}_\${LOC1}_\${LOC2}_pop1.txt     # create subfiles for each dataset with list of samples correponding to 1 component of pairwise combination
grep \${LOC2} \${POP} > $BASE_DIR/outputs/8_fst/\${SP}_\${LOC1}_\${LOC2}_pop2.txt     # this will be used as an input for FST calculations


vcftools --gzvcf \${VCF} \                                                                     # use VCFTOOLS to calculate pairwise FST
      --weir-fst-pop $BASE_DIR/outputs/8_fst/\${SP}_\${LOC1}_\${LOC2}_pop1.txt \
      --weir-fst-pop $BASE_DIR/outputs/8_fst/\${SP}_\${LOC1}_\${LOC2}_pop2.txt \
      --fst-window-step 5000 \                                                                 # use SNP windows of 50kb
      --fst-window-size 50000 \
      --out $BASE_DIR/outputs/8_fst/\${SP}_\${LOC1}_\${LOC2}.50k 2> $BASE_DIR/outputs/8_fst/\${SP}_\${LOC1}_\${LOC2}.50k.log

vcftools --gzvcf \${VCF} \                  
      --weir-fst-pop $BASE_DIR/outputs/8_fst/\${SP}_\${LOC1}_\${LOC2}_pop1.txt \
      --weir-fst-pop $BASE_DIR/outputs/8_fst/\${SP}_\${LOC1}_\${LOC2}_pop2.txt \
      --fst-window-size 10000 \                                                                # use SNP windows of 10kb
      --fst-window-step 1000 \
      --out $BASE_DIR/outputs/8_fst/\${SP}_\${LOC1}_\${LOC2}.10k 2> $BASE_DIR/outputs/8_fst/\${SP}_\${LOC1}_\${LOC2}.10k.log


gzip $BASE_DIR/outputs/8_fst/\${SP}_\${LOC1}_\${LOC2}.50k.windowed.weir.fst           # zip both windowed vcf files created 
gzip $BASE_DIR/outputs/8_fst/\${SP}_\${LOC1}_\${LOC2}.10k.windowed.weir.fst


EOA



# ------------------------------------------------------------------------------
# Job 5 create file of global FST per pairwise comparisons

jobfile5=5_global.tmp # temp file
cat > $jobfile5 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=5_global
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/5_global_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/5_global_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=04:30:00


cd $BASE_DIR/outputs/8_fst/                                                           # move to folder where log files from FST calculation can be found

cat $BASE_DIR/outputs/8_fst/\*.50k.log | \                                            # find lines where mean and weighted FST estimates are reported
    grep -E 'Weir and Cockerham|--out' | \
    grep -A 3 50k | \
    sed '/^--/d; s/^.*--out //g; s/.50k//g; /^Output/d; s/Weir and Cockerham //g; s/ Fst estimate: /\t/g' | \
    paste - - - | \
    cut -f 1,3,5 | \
    sed 's/^\\(...\\)-/\\1\\t/g' > $BASE_DIR/outputs/8_fst/fst_globals.txt            # create global FST file as a table with entry for each pairwise comparison (24 rows in total) and columns for name of the comparison, mean FST, weighted FST (3 columns) 


EOA



# ------------------------------------------------------------------------------
# Job 6 create plots for FST

jobfile6=6_plots.tmp # temp file
cat > $jobfile6 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=6_plots
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/6_plots_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/6_plots_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=04:30:00


ml hpc-env/8.3                                                                        # need to load environments that have the righ ressources for plots 
ml R/4.0.2-foss-2019b
ml FriBidi
ml HarfBuzz

mkdir $BASE_DIR/figures/fst/                                               # create an additional folder for the specific plots that will be created

Rscript $BASE_DIR/R/fst_plots.R $BASE_DIR/outputs/8_fst/ $BASE_DIR/outputs/8_fst/fst_globals.txt $BASE_DIR/figures/fst/


EOA



# ------------------------------------------------------------------------------
# Job 7 test

jobfile7=7_test.tmp # temp file
cat > $jobfile7 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=7_test
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/7_test_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/7_test_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=04:30:00


for k in puebel nigbel maybel indbel pueboc nigboc uniflo unipue puepue chlpue;
do
for j in puebel nigbel maybel indbel pueboc nigboc uniflo unipue puepue chlpue;
do
echo \$k \$j >> $BASE_DIR/outputs/8_fst/pairwise_comparisons.txt;
done;
done

PAIR=$BASE_DIR/outputs/8_fst/pairwise_comparisons.txt

IFS=$'\n'
for i in \$(cat \${PAIR});
do
  echo "\$i"
  a="\${i:0:6}"
  echo "\$a
  b="\${i:7:12}"
  echo "\$b

  for m in \$(cat \${PAIR});
  do
    echo "\$m"
    c="\${m:0:6}
    echo "\$c
    d="\${m:7:12}
    echo "\$d

    if [ "\$a" = "\$d" ] && [ "\$b" = "\$c" ];
    then
      echo "same pair"
      :
    else
      echo "\$a \$b" >> pairwise_comparison1.txt
    fi
  done
done


EOA



# ********** Schedule the job launching ***********
# -------------------------------------------------

if [ "$JID_RES" = "jid1" ] || [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ];
then
  echo "*****   0_keep_species  : DONE         **"
else
  jid0=$(sbatch ${jobfile0})
fi


if [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ];
then
  echo "*****   1_all_hamlets   : DONE         **"
elif [ "$JID_RES" = jid1 ]
then
  jid1=$(sbatch ${jobfile1})
else
  jid1=$(sbatch --dependency=afterok:${jid0##* } ${jobfile1})
fi


if [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ];
then
  echo "*****   2_multi         : DONE         **"
elif [ "$JID_RES" = jid2 ]
then
  jid2=$(sbatch ${jobfile2})
else
  jid2=$(sbatch --dependency=afterok:${jid1##* } ${jobfile2})
fi


if [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ]  || [ "$JID_RES" = "jid6" ];
then
  echo "*****   3_prep_pairwise : DONE         **"
elif [ "$JID_RES" = jid3 ]
then
  jid3=$(sbatch ${jobfile3})
else
  jid3=$(sbatch --dependency=afterok:${jid2##* } ${jobfile3})
fi


if [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ];
then
  echo "*****   4_pairwise_fst  : DONE         **"
elif [ "$JID_RES" = jid4 ]
then
  jid4=$(sbatch ${jobfile4})
else
  jid4=$(sbatch --dependency=afterok:${jid3##* } ${jobfile4})
fi

if [ "$JID_RES" = "jid6" ];
then
  echo "*****   5_global        : DONE         **"
elif [ "$JID_RES" = "jid5" ]
then
  jid5=$(sbatch ${jobfile5})
else
  jid5=$(sbatch --dependency=afterok:${jid4##* } ${jobfile5})
fi


if [ "$JID_RES" = "jid6" ];
then
  jid6=$(sbatch ${jobfile6})
else
  jid6=$(sbatch --dependency=afterok:${jid5##* } ${jobfile6})
fi


if [ "$JID_RES" = "jid7" ];
then
  jid7=$(sbatch ${jobfile7})
else
  jid7=$(sbatch --dependency=afterok:${jid6##* } ${jobfile7})
fi


# ******** Remove useless files and folders *******
# -------------------------------------------------

rm *tmp
