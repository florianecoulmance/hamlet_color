#!/bin/bash
# by: Floriane Coulmance: 01/04/2020
# usage:
# sbatch hybrids.sh -i <PATH> -j <JOB_ID>
# ------------------------------------------------------------------------------
# PATH corresponds to the path to the base directory, all outputs and necessary
# folder will be created by the script
# JOB_ID corresponds to ids of jobs from where you want the analysis to be ran
# within the pipeline
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
mkdir $BASE_DIR/outputs/9_newhyb/

# Annex folder for list of files
mkdir $BASE_DIR/outputs/lof/



# ********* Jobs creation *************************
# -------------------------------------------------

# --------------------------- PREPARATION -------------------------------------#

# ------------------------------------------------------------------------------
# Job 0 creates the species-location pairs

jobfile0=0_pairs.tmp # temp file
cat > $jobfile0 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=0_pairs
#SBATCH --partition=carl.p
#SBATCH --array=0-11
#SBATCH --output=$BASE_DIR/logs/0_pairs_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/0_pairs_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=04:00:00
 

INPUT=$BASE_DIR/outputs/6_genotyping/6_1_snp/snp_filterd.vcf.gz                                     # input the genotyping file with SNPs created by step 14 of genotyping.sh pipeline

tr=(indbel/maybel indbel/nigbel indbel/puebel maybel/nigbel maybel/puebel nigbel/puebel           # create list of species pairwise comparison keeping populations that have >5 individuals represented
nigboc/pueboc nigboc/uniboc pueboc/uniboc chlpue/puepue chlpue/unipue puepue/unipue)

PAIR=\${tr[\${SLURM_ARRAY_TASK_ID}]}                                                              # create 1 job (array job) per species pairwise comparison 
echo \${PAIR}                                               
POP_UN=\${PAIR%/*}                                                                                # extract both population of the pairwise comparison 
echo \${POP_UN}                                                                                   # print the pairwise comparison label and both population label
POP_DEUX=\${PAIR##*/}
echo \${POP_DEUX}

vcfsamplenames \${INPUT} | grep \${POP_UN} > $BASE_DIR/outputs/9_newhyb/\${POP_UN}.pop           # create list of individuals per population of the pairwise comparison 
vcfsamplenames \${INPUT} | grep \${POP_DEUX} > $BASE_DIR/outputs/9_newhyb/\${POP_DEUX}.pop       # put the list of individuals in 1 file for both of the populations

vcftools --gzvcf \${INPUT} --weir-fst-pop $BASE_DIR/outputs/9_newhyb/\${POP_UN}.pop --weir-fst-pop $BASE_DIR/outputs/9_newhyb/\${POP_DEUX}.pop --stdout > $BASE_DIR/outputs/9_newhyb/\${POP_UN}_\${POP_DEUX}.fst.tsv

ls -1 $BASE_DIR/outputs/9_newhyb/*.fst.tsv > $BASE_DIR/outputs/lof/8_fst.fofn                 # add all pairwise Fst file names and path to a file


EOA



# ------------------------------------------------------------------------------
# Job 1 creates the files with information on the position of the 800 most 
# differentiated SNPs

jobfile1=1_snp.tmp # temp file
cat > $jobfile1 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=1_snp
#SBATCH --partition=carl.p
#SBATCH --array=1-12
#SBATCH --output=$BASE_DIR/logs/1_snp_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/1_snp_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=00:15:00


module load hpc-env/8.3                                                                           # load the environment for the latest version of R on the Oldenburg Cluster
module load R/4.0.2-foss-2019b

INPUT_PAIRS=$BASE_DIR/outputs/lof/8_fst.fofn                                                      # input file with fst file names and path
PAIR=\$(cat \${INPUT_PAIRS} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)                        # create a job for each filename in the input file and extract the path to the file corresponding to the job array number
echo \${PAIR}

# PREFIX=\${PAIR##*/} 
# echo \${PREFIX}

PREFIX2=\${PAIR%%.*}                                                                              # create prefix (with path and future file name without the extension)
echo \${PREFIX2}

Rscript $BASE_DIR/R/filter_snps.R \${PAIR} 800 \${PREFIX2}                                        # execute R code that will select 800 random positions in Fst file


EOA



# ------------------------------------------------------------------------------
# Job 2 creates the genotyping files with only the 800 most differentiated SNPs

jobfile2=2_vcf.tmp # temp file
cat > $jobfile2 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=2_vcf
#SBATCH --partition=carl.p
#SBATCH --array=1-12
#SBATCH --output=$BASE_DIR/logs/2_vcf_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/2_vcf_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=03:00:00


INPUT=$BASE_DIR/outputs/6_genotyping/6_1_snp/snp_filterd.vcf.gz                                 # input the genotyping file with SNPs created by step 14 of genotyping.sh pipeline
echo \${INPUT}

INPUT_PAIRS=$BASE_DIR/outputs/lof/8_fst.fofn                                                      # input the list of pairwise Fst files created at job0
PAIR=\$(cat \${INPUT_PAIRS} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)                        # create a job per Fst file
echo \${PAIR}

PREFIX=\${PAIR##*/}                                                                               # retrieve file prefix including both population names of the pairwise Fst
echo \${PREFIX}
PREFIX=\${PREFIX%%.*}
echo \${PREFIX}

POP_UN=\${PREFIX%_*}                                                                              # retrieve both population names separately 
echo \${POP_UN}
POP_DEUX=\${PREFIX##*_}
echo \${POP_DEUX}

vcftools --gzvcf \${INPUT} --keep $BASE_DIR/outputs/9_newhyb/\${POP_UN}.pop --keep $BASE_DIR/outputs/9_newhyb/\${POP_DEUX}.pop --thin 5000 --out $BASE_DIR/outputs/9_newhyb/newHyb.\${POP_UN}_\${POP_DEUX} --positions $BASE_DIR/outputs/9_newhyb/\${PREFIX}_800SNPs.snps --recode                                                                                    # from the genotyping file, retrieve all individuals from both populations
                                                                                # and keep the 800 most differentiated positions based on the position file created at job1
                                                                                                              # output with recode.vcf extension

EOA



# ------------------------------------------------------------------------------
# Job 3 creates the files with 800 most differentiated SNPs

jobfile3=3_format.tmp # temp file
cat > $jobfile3 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=3_format
#SBATCH --partition=carl.p
#SBATCH --array=1-12
#SBATCH --output=$BASE_DIR/logs/3_format_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/3_format_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=03:00:00


INPUT_PAIRS=$BASE_DIR/outputs/lof/8_fst.fofn                                                      # input the file of list of pairwise Fst files                                    
PAIR=\$(cat \${INPUT_PAIRS} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)                        # create a job array per Fst file
echo \${PAIR}

PREFIX=\${PAIR##*/}                                                                               # retrieve pairwise comparison prefix
echo \${PREFIX}
PREFIX=\${PREFIX%%.*}
echo \${PREFIX}

POP_UN=\${PREFIX%_*}                                                                              # retrieve both population names of the pairwise comparison
echo \${POP_UN}
POP_DEUX=\${PREFIX##*_}
echo \${POP_DEUX}

grep '#' $BASE_DIR/outputs/9_newhyb/newHyb.\${POP_UN}_\${POP_DEUX}.recode.vcf > $BASE_DIR/outputs/9_newhyb/newHyb.\${POP_UN}_\${POP_DEUX}.80SNPs.vcf               # retrieve all the content of the vcf file that starts with # (info at the beginning and column names)
           

grep -v '#' $BASE_DIR/outputs/9_newhyb/newHyb.\${POP_UN}_\${POP_DEUX}.recode.vcf | shuf -n 80 | sort -k 1 -k2 >> $BASE_DIR/outputs/9_newhyb/newHyb.\${POP_UN}_\${POP_DEUX}.80SNPs.vcf              # select randomly 80 SNPs out of the 800 most differentiated positions
            # add those 80 SNPs it to the previously created vcf files

grep '#CHROM' $BASE_DIR/outputs/9_newhyb/newHyb.\${POP_UN}_\${POP_DEUX}.80SNPs.vcf | cut -f 10- | sed 's/\\t/\\n/g' > $BASE_DIR/outputs/9_newhyb/newHyb.\${POP_UN}_\${POP_DEUX}.80SNPs_individuals.txt           # from the newly 80SNPs dataset created, make a list of individuals' name
                                                                             # keep the indiviuals' name in a new txt file
        

java -Xmx1024m -Xms512M -jar /user/doau0129/miniconda3/share/pgdspider-2.1.1.5-0/PGDSpider2-cli.jar -inputfile $BASE_DIR/outputs/9_newhyb/newHyb.\${POP_UN}_\${POP_DEUX}.80SNPs.vcf -inputformat VCF -outputfile $BASE_DIR/outputs/9_newhyb/newHyb.\${POP_UN}_\${POP_DEUX}.80SNPs.txt -outputformat NEWHYBRIDS -spid $BASE_DIR/ressources/vcf2nh.spid

# use PGDSpider2 to create a newhybrid formatted file to be able to use with NEWHYBRID
# input the vcf file with 80 SNPs
# output the formatted file as a txt file ready to be used in NEWHYBRID


EOA



# ------------------------------------------------------------------------------
# Job 4 creates the files with 800 most differentiated SNPs

jobfile4=4_nh.tmp # temp file
cat > $jobfile4 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=4_nh
#SBATCH --partition=carl.p
#SBATCH --array=1-12
#SBATCH --output=$BASE_DIR/logs/4_nh_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/4_nh_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=3-03:00:00


module load hpc-env/8.3                                                                           # load necessary environment for necessary version of R 
module load R/4.0.2-foss-2019b

INPUT_PAIRS=$BASE_DIR/outputs/lof/8_fst.fofn                                                      # input the file with list of pairwise FST files
PAIR=\$(cat \${INPUT_PAIRS} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)                        # create 1 job per pairwise FST file
echo \${PAIR}

PREFIX=\${PAIR##*/}                                                                               # find the prefix corresponding to pairwise comparison in file name
echo \${PREFIX}
PREFIX=\${PREFIX%%.*}
echo \${PREFIX}

POP_UN=\${PREFIX%_*}                                                                              # find both the populations name of the pairwise Fst prefix
echo \${POP_UN}
POP_DEUX=\${PREFIX##*_}
echo \${POP_DEUX}

mkdir -p $BASE_DIR/outputs/9_newhyb/\${POP_UN}_\${POP_DEUX}_nhinput                              # create an output folder for each pairwise comparison
scp -r $BASE_DIR/outputs/9_newhyb/newHyb.\${POP_UN}_\${POP_DEUX}.80SNPs.txt $BASE_DIR/outputs/8_hybrids/\${POP_UN}_\${POP_DEUX}_nhinput/ # copy the NEWHYBRID formatted file with 80 SNP positions to the new output folder
scp -r $BASE_DIR/outputs/9_newhyb/newHyb.\${POP_UN}_\${POP_DEUX}.80SNPs_individuals.txt $BASE_DIR/outputs/8_hybrids/\${POP_UN}_\${POP_DEUX}_nhinput/ # copy the list of individual from both population of the pairwise comparison to the new output folder

Rscript $BASE_DIR/R/newhybrids.R $BASE_DIR/outputs/9_newhyb/\${POP_UN}_\${POP_DEUX}_nhinput/     # execute the Rscript that performs NEWHYBRIDS


EOA


# ------------------------------------------------------------------------------
# Job 5 make hybrid plot

jobfile5=5_plots.tmp # temp file
cat > $jobfile5 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=5_plots
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/5_plots_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/5_plots_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=3-03:00:00


module load hpc-env/8.3                                                                           # load necessary environment for necessary version of R 
module load R/4.0.2-foss-2019b
ml FriBidi
ml HarfBuzz

mkdir $BASE_DIR/figures/9_newhyb/                                                                 # create figure directory for NEWHYBRID plot

Rscript $BASE_DIR/R/plot_hybrids.R $BASE_DIR/outputs/9_newhyb/ $BASE_DIR/figures/9_newhyb/        # run Rscript for newhybrid plot


EOA



# ********** Schedule the job launching ***********
# -------------------------------------------------

# ------------------------------------------------------------------------------
# DATA PREPARATION

if [ "$JID_RES" = "jid1" ] || [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ];
then
  echo "*****   0_pairs         : DONE         **"
else
  jid0=$(sbatch ${jobfile0})
fi


if [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ];
then
  echo "*****   1_snp           : DONE         **"
elif [ "$JID_RES" = jid1 ]
then
  jid1=$(sbatch ${jobfile1})
else
  jid1=$(sbatch --dependency=afterok:${jid0##* } ${jobfile1})
fi


if [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ]  || [ "$JID_RES" = "jid5" ];
then
  echo "*****   2_align_fastq   : DONE         **"
elif [ "$JID_RES" = jid2 ]
then
  jid2=$(sbatch ${jobfile2})
else
  jid2=$(sbatch --dependency=afterok:${jid1##* } ${jobfile2})
fi


if [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ];
then
  echo "*****   3_format        : DONE         **"
elif [ "$JID_RES" = "jid3" ]
then
  jid3=$(sbatch ${jobfile3})
else
  jid3=$(sbatch --dependency=afterok:${jid2##* } ${jobfile3})
fi


if [ "$JID_RES" = "jid5" ];
then
  echo "*****   4_nh            : DONE         **"
elif [ "$JID_RES" = "jid4" ]
then
  jid4=$(sbatch ${jobfile4})
else
  jid4=$(sbatch --dependency=afterok:${jid3##* } ${jobfile4})
fi


if [ "$JID_RES" = "jid5" ];
then
  jid5=$(sbatch ${jobfile5})
else
  jid5=$(sbatch --dependency=afterok:${jid4##* } ${jobfile5})
fi



# ******** Remove useless files and folders *******
# -------------------------------------------------

rm *tmp
# rm -r $BASE_DIR/outputs/lof/
