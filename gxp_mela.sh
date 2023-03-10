#!/bin/bash
# by: Floriane Coulmance: 08/03/2023
# usage:
# sbatch gxp_mela.sh -i <PATH> -j <JOB_ID> -k <DATASET>
# ------------------------------------------------------------------------------
# PATH corresponds to the path to the base directory, all outputs and necessary
# folder will be created by the script
# JOB_ID corresponds to ids of jobs from where you want the analysis to be ran
# within the pipeline
# DATASET : can be looking like LAB_fullm_54off_59on (for continuous
# phenotype analysis) or like RGB_fullm_k7_color0 (for discrete phenotype analysis),
# it should be prefix of metadata file
# ------------------------------------------------------------------------------



# ********** Allow to enter bash options **********
# -------------------------------------------------

while getopts i:j:k: option
do
case "${option}"
in
i) BASE_DIR=${OPTARG};; # get the base directory path
j) JID_RES=${OPTARG};; # get the jobid from which you want to resume
k) DATASET=${OPTARG};; # get the dataset type
esac
done



# ********* Create necessary variables *********
# -------------------------------------------------

# Get color space variable from dataset name
COLOR_SPACE="${DATASET%%_*}"
echo $COLOR_SPACE

if [[ $DATASET == *"k"* ]];
then
  TYPE="discrete"
else
  TYPE="continuous"
fi

echo $TYPE



# ********* Create necessary repositories *********
# -------------------------------------------------

# Repo for gxp outputs
mkdir $BASE_DIR/outputs/7_gxp/

# Repo for the corresponding dataset
mkdir $BASE_DIR/outputs/7_gxp/$TYPE/
mkdir $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/
mkdir $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/

# Repo for the figures
mkdir $BASE_DIR/figures/
mkdir $BASE_DIR/figures/7_gxp/
mkdir $BASE_DIR/figures/7_gxp/$TYPE/
mkdir $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/
mkdir $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/


# ********* Jobs creation *************************
# -------------------------------------------------

# --------------------------- PREPARATION -------------------------------------#

# ------------------------------------------------------------------------------
# Job 0 convert genotyping file to plink format 

jobfile0=0_convert_plink.tmp # temp file
cat > $jobfile0 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=0_convert_plink
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/0_convert_plink_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/0_convert_plink_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=06:00:00

ml hpc-env/8.3                                                                        # need to load environments that have the righ ressources for plots 
ml R/4.0.2-foss-2019b
#ml R-bundle-Bioconductor/3.12-foss-2019b-R-4.0.2
ml FriBidi
ml HarfBuzz 


INPUT_VCF=$BASE_DIR/outputs/6_genotyping/6_1_snp/snp_filterd.vcf.gz                    # Input the bi allelic genotyping file

# bcftools view -S $BASE_DIR/outputs/6_genotyping/6_1_snp/pue_nig_uni_list -o $BASE_DIR/outputs/6_genotyping/6_1_snp/snp_filterd_pnu.vcf.gz \${INPUT_VCF}
# tabix -p vcf $BASE_DIR/outputs/6_genotyping/6_1_snp/snp_filterd_pnu.vcf.gz                                        # create index for the file just created

INPUT_PNU=$BASE_DIR/outputs/6_genotyping/6_1_snp/snp_filterd_pnu.vcf.gz

# Perform genotyping PCA
FILE=\${INPUT_VCF##*/}                                                                      # extract prefix for further analysis and naming of files
PREFIX=\${FILE%%.*}
echo \$INPUT_PNU
echo \$FILE
echo \$PREFIX

TEST=\${PREFIX}
echo \$TEST

# Rscript $BASE_DIR/R/genotyping_pca.R \${INPUT_PNU} $BASE_DIR/outputs/pca/ \${PREFIX} $BASE_DIR/figures/genotyping_pca/             # run the R script for plots
# Rscript $BASE_DIR/R/genotyping_pca.R \${INPUT_VCF} $BASE_DIR/outputs/pca/ \${PREFIX} $BASE_DIR/figures/genotyping_pca/             # run the R script for plots

# Convert the genotyping file to plink format
vcftools \
      --gzvcf \${INPUT_PNU} \
      --plink \
      --out $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/GxP_plink

# Convert to hapmap / not mandatory to run
plink \
      --file $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/GxP_plink \
      --recode12 \
      --out $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/hapmap

# Convert plink genotyping file to binary files to be used in GWAS
plink \
    --noweb \
    --file $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/GxP_plink \
    --make-bed \
    --out $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/GxP_plink_binary

# Save a copy of the binary .fam file
cp $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/GxP_plink_binary.fam \
   $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/GxP_plink_binary_sauvegarde.fam


EOA



# ------------------------------------------------------------------------------
# Job 1 prepare the phenotyping files for GEMMA and MVPLINK analysis

jobfile1=1_prep.tmp # temp file
cat > $jobfile1 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=1_prep.tmp
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/1_prep_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/1_prep_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500M
#SBATCH --time=04:00:00


# PREPARATION OF THE MERGED GENOTYPES AND PHENOTYPES ----------------------------

fam=$BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/GxP_plink_binary.fam                             # Path to the genotyping plink binary file
pheno=$BASE_DIR/metadata/${DATASET}_PCs.csv                                           # Path to the phenotype file
echo \${pheno}

tr=(PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10)                # Create file with list of PCs for further job array steps
printf "%s\n" "\${tr[@]}" > $BASE_DIR/outputs/lof/pcs.fofn

sort -k1 \${fam} > $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/GxP_plink_binary_sorted.fam                                                                    # Sort the .fam file on the individuals' labels column

srted=$BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/GxP_plink_binary_sorted.fam

# Merge genotyping binary .fam file and phenotype file
awk -F ";" '{print \$17,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' \${pheno} | \
sort -k1 > $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/pheno_intermediate1

join \${srted} $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/pheno_intermediate1 | \
awk -F " " '{print \$1,\$2,\$3,\$4,\$5,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16}' \
> $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/pheno_intermediate

# Format the merged file for use in GWAS
echo -e 'label Within_family_ID ID_father ID_mother Sex PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10' \
> $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/pheno_table.fam && cat $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/pheno_intermediate \
>> $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/pheno_table.fam


# PREPARATION OF THE FILES FOR MVPLINK -------------------------------------------

MAP=$BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/GxP_plink.map                                    # Input the plink .map file   
PED=$BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/GxP_plink.ped                                    # Input the plink .ped file
cp \${MAP} $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/gwas_multi.map                            # Copy both files to new files with different names, necessarry for MVPLINK analysis
cp \${PED} $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/gwas_multi.ped

P=$BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/pheno_table.fam                                    # Input the just created genotype + phenotype file

# Select needed column for MVPLINK analysis and output to .phen file with same naming as previous step
awk -F " " '{print \$1,\$2,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15}' \${P} \
> $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/gwas_multi.phen

var="FID IID PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10"                                # Add column name row to file 
sed -i "1s/.*/\${var}/" $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/gwas_multi.phen


EOA



# ------------------------------------------------------------------------------
# Job 2 run MVPLINK

jobfile2=2_mvplink.tmp # temp file
cat > $jobfile2 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=2_mvplink.tmp
#SBATCH --partition=carl.p
#SBATCH --array=1-14
#SBATCH --output=$BASE_DIR/logs/2_mvplink_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/2_mvplink_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90G
#SBATCH --time=1-00:00:00


INPUT_MV=$BASE_DIR/outputs/lof/mvplink.fofn                                           # This files of list of combination of multivariate analysis has to be created outside this pipeline

MV=\$(cat \${INPUT_MV} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)                 # Create a job for all the possible phenotypes and the associated .fam file with just one phenotype at a time
echo \${MV}

NAME="\$(cut -d ' ' -f 1 <<<"\${MV}")"
echo \${NAME}

COL="\$(cut -d ' ' -f 2- <<<"\${MV}")"
echo \${COL}

# P=\$(echo  /user/doau0129/work/software/plink.multivariate --file gwas_multi --mqfam --mult-pheno gwas_multi.phen --pheno-name \${COL})
# echo \${P}

# Run the MVPLINK command on the combinations of multivariate traits
/user/doau0129/work/software/plink.multivariate --file $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/gwas_multi --mqfam --mult-pheno $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/gwas_multi.phen --pheno-name \${COL} --out $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/\${NAME}


EOA



# ------------------------------------------------------------------------------
# Job 3 format MVPLINK results

jobfile3=3_slider.tmp # temp file
cat > $jobfile3 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=3_slider.tmp
#SBATCH --partition=carl.p
#SBATCH --array=1-14
#SBATCH --output=$BASE_DIR/logs/3_slider_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/3_slider_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90G
#SBATCH --time=1-00:00:00


body() {                                                                              # Function to read and print to each line
	IFS= read -r header
	printf '%s\n' "\$header"
	"\$@"
}

INPUT_MV=$BASE_DIR/outputs/lof/mvplink.fofn                                           # Input list of multivariate combinations (created outside the pipeline)

MV=\$(cat \${INPUT_MV} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)                 # Create a job for all the possible phenotypes and the associated .fam file with just one phenotype at a time
echo \${MV}

NAME="\$(cut -d ' ' -f 1 <<<"\${MV}")"                                                # Find the name of the combination
echo \${NAME}

FILE=$BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/\${NAME}.mqfam.total                            # Input output of mvplink analysis
echo \${FILE}

# Pre clean the output file
awk '{print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8}' \${FILE} | \
cut -d ' ' -f 2,3,6-9 | \
sed 's/SNP BP /CHROM POS /g' | \
awk '{sub(/\:.*$/,"",\$1); print \$0}' | \
awk '{if (\$3!="NA"){ print}}' | \
body sort -k1,1 -k2,2n | \
gzip > $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/\${NAME}.mvplink.txt.gz

win5=50000                                                                            # Create different sets of parameters
step5=5000
echo \${win5}
echo \${step5}

win1=10000
step1=1000
echo \${win1}
echo \${step1}

# Run the genome average of MVPLINK output with different set of parameters
$BASE_DIR/sh/mvplink_slider.sh $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/\${NAME}.mvplink.txt.gz \${win5} \${step5}
$BASE_DIR/sh/mvplink_slider.sh $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/\${NAME}.mvplink.txt.gz \${win1} \${step1}

# No average just logarithmic of GWAS values
$BASE_DIR/sh/mvplink_log.sh $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/\${NAME}.mvplink.txt.gz


EOA



# ------------------------------------------------------------------------------
# Job 4 GxP plots

jobfile4=4_plots.tmp # temp file
cat > $jobfile4 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=4_plots.tmp
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/4_plots_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/4_plots_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=22G
#SBATCH --time=04:00:00


ml hpc-env/8.3
ml R/4.0.2-foss-2019b
ml FriBidi
ml HarfBuzz

ls -1 $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/*.mvplink.50k.5k.txt.gz > $BASE_DIR/outputs/lof/mvplink_50k.fofn
ls -1 $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/*.assoc.50k.5k.txt.gz > $BASE_DIR/outputs/lof/assoc_50k.fofn

rm $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/\*.tmp

# Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ univariate_gemma
# Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ multivariate_plink_PC1
# Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$DATASET/ $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ multivariate_plink_PC2
# Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$DATASET/ $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ multivariate_plink_PC3
# Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$DATASET/ $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ multivariate_plink_PC4
# Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$DATASET/ $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ multivariate_plink_PC5
# Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$DATASET/ $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ multivariate_plink_PC6
# Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$DATASET/ $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ multivariate_plink_PC7
# Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$DATASET/ $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ multivariate_plink_PC8
# Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$DATASET/ $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ multivariate_plink_PC9
# Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ multivariate_plink_byPCs

# Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ univariate_plink

# Rscript $BASE_DIR/R/gxp_heatmap_plots.R $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ $DATASET $BASE_DIR/metadata/ $BASE_DIR/images/$TYPE/$COLOR_SPACE/$DATASET/

Rscript $BASE_DIR/R/gxp_pca_univariate.R $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ $DATASET $BASE_DIR/metadata/ $BASE_DIR/images/$TYPE/$COLOR_SPACE/$DATASET/


EOA



# ------------------------------------------------------------------------------
# Job 5 create zoom plot on region of interest

jobfile5=5_snps.tmp # temp file
cat > $jobfile5 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=5_snps.tmp
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/5_snps_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/5_snps_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=04:00:00


ml hpc-env/8.3
ml R-bundle-Bioconductor/3.12-foss-2019b-R-4.0.2
ml FriBidi
ml HarfBuzz


# Input the multivariate analysis of choice
AVG=$BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/PC1_5.mvplink.50k.5k.txt.gz

# Extract name of the multivariate analysis, color space and mask names
B=\$(basename "\${AVG}")
echo \${B}
NAME=\${B%%.*}
echo \${NAME}

EFF=${DATASET%%_*}
echo \${EFF}

C=${DATASET#*_}
echo \${C}
M=\${C%%_*}
echo \${M}

# Based on mask name, input mask file
if [ "\${M}" = "pnudctm" ];
then
  MASK=$BASE_DIR/ressources/images/pue_nig_uni_dct_mask.tif
elif [ "\${M}" = "pnufullm" ]
then
  MASK=$BASE_DIR/ressources/images/pue_nig_uni_full_mask.tif
else
  echo "please verify your dataset folder spelling in the specified -k parameter"
fi

echo \${MASK}

# Create appropriate subfolders 
mkdir $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/\${NAME}/
mkdir $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/\${NAME}/

# Run the R and python plot analyses
Rscript $BASE_DIR/R/gxp_save_snps.R $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ \${B} $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/\${NAME}/ \${NAME}


EOA



# ------------------------------------------------------------------------------
# Job 6 create heatmaps for most associated SNPs

jobfile6=6_heat.tmp # temp file
cat > $jobfile6 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=6_heat.tmp
#SBATCH --partition=carl.p
#SBATCH --array=0-23
#SBATCH --output=$BASE_DIR/logs/6_heat_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/6_heat_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=04:00:00


ml hpc-env/8.3
ml R/4.0.2-foss-2019b
ml FriBidi
ml HarfBuzz

# Input the multivariate analysis of choice
AVG=$BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/PC1_5.mvplink.50k.5k.txt.gz

# Extract name of the multivariate analysis, color space and mask names
B=\$(basename "\${AVG}")
echo \${B}
NAME=\${B%%.*}
echo \${NAME}

EFF=${DATASET%%_*}
echo \${EFF}

C=${DATASET#*_}
echo \${C}
M=\${C%%_*}
echo \${M}

# Based on mask name, input mask file
if [ "\${M}" = "pnudctm" ];
then
  MASK=$BASE_DIR/ressources/images/pue_nig_uni_dct_mask.tif
elif [ "\${M}" = "pnufullm" ]
then
  MASK=$BASE_DIR/ressources/images/pue_nig_uni_full_mask.tif
else
  echo "please verify your dataset folder spelling in the specified -k parameter"
fi

echo \${MASK}

chr=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24)

# Create a job for all the possible phenotypes and the associated .fam file with just one phenotype at a time
LG=\${chr[\${SLURM_ARRAY_TASK_ID}]}
echo \${LG}

CHROM="LG\${LG}"
FILE_N="\${NAME}_LG\${LG}"
echo \${FILE_N}

# python3 $BASE_DIR/python/plot_snp_heatmap.py $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/\${NAME}/ \${FILE_N}.snp.txt $BASE_DIR/images/$TYPE/$COLOR_SPACE/$DATASET/${DATASET}_modifiedImage.csv \${MASK} $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/\${NAME}/ \${EFF}
Rscript $BASE_DIR/R/gxp_zooms.R $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ \${B} $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/\${NAME}/ \${NAME} \${CHROM}
Rscript $BASE_DIR/R/gxp_pc1_5_3peaks.R $BASE_DIR/outputs/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/ $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/

EOA



# ********** Schedule the job launching ***********
# -------------------------------------------------

if [ "$JID_RES" = "jid1" ] || [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ];
then
  echo "*****   0_convert_plink : DONE         **"
else
  jid0=$(sbatch ${jobfile0})
fi


if [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ];
then
  echo "*****   1_prep          : DONE         **"
elif [ "$JID_RES" = jid1 ]
then
  jid1=$(sbatch ${jobfile1})
else
  jid1=$(sbatch --dependency=afterok:${jid0##* } ${jobfile1})
fi


if [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ];
then
  echo "*****   2_gemma         : DONE         **"
elif [ "$JID_RES" = "jid2" ]
then
  jid2=$(sbatch ${jobfile2})
else
  jid2=$(sbatch --dependency=afterok:${jid1##* } ${jobfile2})
fi


if [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ];
then
  echo "*****   3_windows       : DONE         **"
elif [ "$JID_RES" = "jid3" ]
then
  jid3=$(sbatch ${jobfile3})
else
  jid3=$(sbatch --dependency=afterok:${jid2##* } ${jobfile3})
fi


if [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ];
then
  echo "*****   4_mvplink       : DONE         **"
elif [ "$JID_RES" = "jid4" ]
then
  jid4=$(sbatch ${jobfile4})
else
  jid4=$(sbatch --dependency=afterok:${jid3##* } ${jobfile4})
fi


if [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ];
then
  echo "*****   5_slider        : DONE         **"
elif [ "$JID_RES" = "jid5" ]
then
  jid5=$(sbatch ${jobfile5})
else
  jid5=$(sbatch --dependency=afterok:${jid4##* } ${jobfile5})
fi
 

if [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ];
then
  echo "*****   6_plots         : DONE         **"
elif [ "$JID_RES" = "jid6" ]
then
  jid6=$(sbatch ${jobfile6})
else
  jid6=$(sbatch --dependency=afterok:${jid5##* } ${jobfile6})
fi
 

if [ "$JID_RES" = "jid8" ];
then
  echo "*****   7_zooms         : DONE         **"
elif [ "$JID_RES" = "jid7" ]
then
  jid7=$(sbatch ${jobfile7})
else
  jid7=$(sbatch --dependency=afterok:${jid6##* } ${jobfile7})
fi
 

if [ "$JID_RES" = "jid8" ];
then
  jid8=$(sbatch ${jobfile8})
else
  jid8=$(sbatch --dependency=afterok:${jid7##* } ${jobfile8})
fi




# ******** Remove useless files and folders *******
# -------------------------------------------------

rm *tmp
# rm -r $BASE_DIR/outputs/lof/









# # ------------------------------------------------------------------------------
# # Job 7 create zoom plot on region of interest

# jobfile7=7_zooms.tmp # temp file
# cat > $jobfile7 <<EOA # generate the job file
# #!/bin/bash
# #SBATCH --job-name=7_zooms.tmp
# #SBATCH --partition=carl.p
# #SBATCH --array=1-55
# #SBATCH --output=$BASE_DIR/logs/7_zooms_%A_%a.out
# #SBATCH --error=$BASE_DIR/logs/7_zooms_%A_%a.err
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=1
# #SBATCH --mem-per-cpu=40G
# #SBATCH --time=04:00:00


# ml hpc-env/8.3
# ml R-bundle-Bioconductor/3.12-foss-2019b-R-4.0.2
# ml FriBidi
# ml HarfBuzz


# ls -1 $BASE_DIR/outputs/7_gxp/$DATASET/*.mvplink.50k.5k.txt.gz > $BASE_DIR/outputs/lof/$DATASET_mvplink_50k.fofn

# INPUT_AVG=$BASE_DIR/outputs/lof/$DATASET_mvplink_50k.fofn

# #Create a job for all the possible phenotypes and the associated .fam file with just one phenotype at a time
# AVG=\$(cat \${INPUT_AVG} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
# echo \${AVG}

# B=\$(basename "\${AVG}")
# echo \${B}
# NAME=\${B%%.*}
# echo \${NAME}
# EFF=${DATASET%%_*}
# echo \${EFF}

# C=${DATASET#*_}
# echo \${C}
# M=\${C%%_*}
# echo \${M}


# if [ "\${M}" = "fullm" ];
# then
#   MASK=$BASE_DIR/ressources/full_mask.tif
# elif [ "\${M}" = "bodym" ]
# then
#   MASK=$BASE_DIR/ressources/body_mask.tif
# else
#   echo "please verify your dataset folder spelling in the specified -k parameter"
# fi


# echo \${MASK}


# mkdir $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/\${NAME}/
# mkdir $BASE_DIR/outputs/7_gxp/$DATASET/\${NAME}/

# Rscript $BASE_DIR/R/gxp_zooms.R $BASE_DIR/outputs/7_gxp/$DATASET/ \${B} $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/\${NAME}/ \${NAME}
# python3 $BASE_DIR/python/plot_snp_heatmap.py $BASE_DIR/outputs/7_gxp/$DATASET/\${T}/ \${NAME}.snp.txt $BASE_DIR/images/$DATASET/${DATASET}_modifiedImage.csv \${MASK} $BASE_DIR/figures/7_gxp/$TYPE/$COLOR_SPACE/$DATASET/\${T}/\${NAME}/ \${EFF}


# EOA
