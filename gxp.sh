#!/bin/bash
# by: Floriane Coulmance: 17/02/2021
# usage:
# gxp.sh -i <PATH> -j <JOB_ID>
# ------------------------------------------------------------------------------
# PATH corresponds to the path to the base directory, all outputs and necessary
# folder will be created by the script
# JOB_ID corresponds string ids from where you want  the script to be ran
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



# ********* Create necessary repositories *********
# -------------------------------------------------

# Repo for gxp outputs
mkdir $BASE_DIR/outputs/7_gxp/

# Repo for the corresponding dataset
mkdir $BASE_DIR/outputs/7_gxp/$DATASET/



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
#SBATCH --time=02:30:00


INPUT_BI=$BASE_DIR/outputs/6_genotyping/6_1_snp/snp_filterd.vcf.gz

vcftools \
      --gzvcf \${INPUT_BI} \
      --plink \
      --out $BASE_DIR/outputs/7_gxp/$DATASET/GxP_plink

plink \
      --file $BASE_DIR/outputs/7_gxp/$DATASET/GxP_plink \
      --recode12 \
      --out $BASE_DIR/outputs/7_gxp/$DATASET/hapmap

EOA



# ------------------------------------------------------------------------------
# Job 1 convert genotyping file to plink binary format

jobfile1=1_plink_binary.tmp # temp file
cat > $jobfile1 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=1_plink_binary
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/1_plink_binary_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/1_plink_binary_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=02:30:00


plink \
    --noweb \
    --file $BASE_DIR/outputs/7_gxp/$DATASET/GxP_plink \
    --make-bed \
    --out $BASE_DIR/outputs/7_gxp/$DATASET/GxP_plink_binary

cp $BASE_DIR/outputs/7_gxp/$DATASET/GxP_plink_binary.fam $BASE_DIR/outputs/7_gxp/$DATASET/GxP_plink_binary_sauvegarde.fam

EOA



# ------------------------------------------------------------------------------
# Job 2 prepare files for GEMMA analysis

jobfile2=2_prep.tmp # temp file
cat > $jobfile2 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=2_prep.tmp
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/2_prep_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/2_prep_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500M
#SBATCH --time=00:15:00


fam=$BASE_DIR/outputs/7_gxp/$DATASET/GxP_plink_binary.fam
pheno=$BASE_DIR/metadata/${DATASET}_PCs.csv
echo \${pheno}

tr=(PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10)
printf "%s\n" "\${tr[@]}" > $BASE_DIR/outputs/lof/pcs.fofn

#Create joint phenotype and .fam file with all phenotypes
awk -F ";" '{print \$17,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' \${pheno} > $BASE_DIR/outputs/7_gxp/$DATASET/pheno_intermediate1
sort -k1 $BASE_DIR/outputs/7_gxp/$DATASET/pheno_intermediate1 > $BASE_DIR/outputs/7_gxp/$DATASET/pheno_intermediate2
sort -k1 \${fam}
join \${fam} $BASE_DIR/outputs/7_gxp/$DATASET/pheno_intermediate2 > $BASE_DIR/outputs/7_gxp/$DATASET/pheno_intermediate3
awk -F " " '{print \$1,\$2,\$3,\$4,\$5,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16}' $BASE_DIR/outputs/7_gxp/$DATASET/pheno_intermediate3 > $BASE_DIR/outputs/7_gxp/$DATASET/pheno_intermediate4
echo -e 'label Within_family_ID ID_father ID_mother Sex PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10' > $BASE_DIR/outputs/7_gxp/$DATASET/pheno_table.fam && cat $BASE_DIR/outputs/7_gxp/$DATASET/pheno_intermediate4 >> $BASE_DIR/outputs/7_gxp/$DATASET/pheno_table.fam


EOA



# ------------------------------------------------------------------------------
# Job 3 run univariate GEMMA

jobfile3=3_gemma.tmp # temp file
cat > $jobfile3 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=3_gemma.tmp
#SBATCH --partition=carl.p
#SBATCH --array=1-10
#SBATCH --output=$BASE_DIR/logs/3_gemma_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/3_gemma_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=04:00:00


body() {
	IFS= read -r header
	printf '%s\n' "\$header"
	"\$@"
}

fam=$BASE_DIR/outputs/7_gxp/$DATASET/GxP_plink_binary.fam
INPUT_TR=$BASE_DIR/outputs/lof/pcs.fofn

#Create a job for all the possible phenotypes and the associated .fam file with just one phenotype at a time
TRAITS=\$(cat \${INPUT_TR} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \${TRAITS}

BASE_NAME=\$(echo  \${fam} | sed 's/.fam//g')
echo \${BASE_NAME}

echo \${BASE_NAME}_\${TRAITS}

cp \${BASE_NAME}.bed \${BASE_NAME}_\${TRAITS}.bed
cp \${BASE_NAME}.bim \${BASE_NAME}_\${TRAITS}.bim
cp \${BASE_NAME}.log \${BASE_NAME}_\${TRAITS}.log
cp \${BASE_NAME}.nosex \${BASE_NAME}_\${TRAITS}.nosex

awk -v t="\${TRAITS}" 'NR==1 {for (i=1; i<=NF; i++) {f[\$i] = i}} {print \$(f["label"]), \$(f["Within_family_ID"]), \$(f["ID_father"]), \$(f["ID_mother"]), \$(f["Sex"]), \$(f[t])}' $BASE_DIR/outputs/7_gxp/$DATASET/pheno_table.fam >  $BASE_DIR/outputs/7_gxp/$DATASET/pheno_header_\${TRAITS}.fam
sed '1d' $BASE_DIR/outputs/7_gxp/$DATASET/pheno_header_\${TRAITS}.fam > \${BASE_NAME}_\${TRAITS}.fam

mkdir $BASE_DIR/output
mkdir $BASE_DIR/output/$DATASET/

  # 2) create relatedness matrix of samples using gemma
gemma -bfile \${BASE_NAME}_\${TRAITS} -gk 1 -o /$DATASET/\${TRAITS}

  # 3) fit linear model using gemma (-lm)
gemma -bfile \${BASE_NAME}_\${TRAITS} -lm 4 -miss 0.1 -notsnp -o /$DATASET/\${TRAITS}.lm

  # 4) fit linear mixed model using gemma (-lmm)
gemma -bfile \${BASE_NAME}_\${TRAITS} -k output/$DATASET/\${TRAITS}.cXX.txt -lmm 4 -o /$DATASET/\${TRAITS}.lmm

  # 5) reformat output
sed 's/\\trs\\t/\\tCHROM\\tPOS\\t/g; s/\\([0-2][0-9]\\):/\\1\\t/g' output/$DATASET/\${TRAITS}.lm.assoc.txt | \
      cut -f 2,3,9-14 | body sort -k1,1 -k2,2n | gzip > $BASE_DIR/outputs/7_gxp/$DATASET/\${TRAITS}.lm.GxP.txt.gz
sed 's/\\trs\\t/\\tCHROM\\tPOS\\t/g; s/\\([0-2][0-9]\\):/\\1\\t/g' output/$DATASET/\${TRAITS}.lmm.assoc.txt | \
      cut -f 2,3,8-10,13-15 | body sort -k1,1 -k2,2n | gzip > $BASE_DIR/outputs/7_gxp/$DATASET/\${TRAITS}.lmm.GxP.txt.gz

# rm \${BASE_NAME}_\${TRAITS}.bed
# rm \${BASE_NAME}_\${TRAITS}.bim
# rm \${BASE_NAME}_\${TRAITS}.log
# rm \${BASE_NAME}_\${TRAITS}.nosex



EOA



# ------------------------------------------------------------------------------
# Job 4 sliding windows for graphs

jobfile4=4_windows.tmp # temp file
cat > $jobfile4 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=4_windows.tmp
#SBATCH --partition=carl.p
#SBATCH --array=1-10
#SBATCH --output=$BASE_DIR/logs/4_windows_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/4_windows_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=22G
#SBATCH --time=02:30:00

INPUT_TR=$BASE_DIR/outputs/lof/pcs.fofn

#Create a job for all the possible phenotypes and the associated .fam file with just one phenotype at a time
TRAITS=\$(cat \${INPUT_TR} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \${TRAITS}

lm=$BASE_DIR/outputs/7_gxp/$DATASET/\${TRAITS}.lm.GxP.txt.gz
echo \${lm}
lmm=$BASE_DIR/outputs/7_gxp/$DATASET/\${TRAITS}.lmm.GxP.txt.gz
echo \${lmm}

win5=50000
step5=5000
echo \${win5}
echo \${step5}

win1=10000
step1=1000
echo \${win1}
echo \${step1}


$BASE_DIR/sh/gxp_slider.sh \${lm} \${win5} \${step5}
$BASE_DIR/sh/gxp_slider.sh \${lm} \${win1} \${step1}
$BASE_DIR/sh/gxp_slider.sh \${lmm} \${win5} \${step5}
$BASE_DIR/sh/gxp_slider.sh \${lmm} \${win1} \${step1}

EOA



# # ------------------------------------------------------------------------------
# # Job 5 multivariate gemma

# jobfile5=5_multi.tmp # temp file
# cat > $jobfile5 <<EOA # generate the job file
# #!/bin/bash
# #SBATCH --job-name=5_multi.tmp
# #SBATCH --partition=carl.p
# #SBATCH --array=1-45
# #SBATCH --output=$BASE_DIR/logs/5_multi_%A_%a.out
# #SBATCH --error=$BASE_DIR/logs/5_multi_%A_%a.err
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=1
# #SBATCH --mem-per-cpu=22G
# #SBATCH --time=04:00:00


# body() {
#         IFS= read -r header
#         printf '%s\n' "\$header"
#         "\$@"
# }


# INPUT_M=$BASE_DIR/outputs/lof/multi.fofn
# fam=$BASE_DIR/outputs/7_gxp/$DATASET/GxP_plink_binary.fam

# #Create a job for all the possible phenotypes and the associated .fam file with just one phenotype at a time
# MULTI=\$(cat \${INPUT_M} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
# echo \${MULTI}

# NAME="\$(cut -d ' ' -f 1 <<<"\${MULTI}")"
# echo \${NAME}

# COL="\$(cut -d ' ' -f 2- <<<"\${MULTI}")"
# echo \${COL}
# IFS=' ' read -r -a C <<<\${COL}
# echo \${C[@]}
# for (( i = 0 ; i < \${#C[@]} ; i++ )) do (( C[\$i]=\${C[\$i]} - 5 )) ; done
# echo \${C[@]}
# echo \${C}

# COL2="\$(echo \${COL} | sed -e 's/ /, $/g')"
# echo \${COL2}

# string="$"
# COLUMNS="\${string}\${COL2}"
# echo \${COLUMNS}

# string2="\\\$1, \\\$2, \\\$3, \\\$4, \\\$5, "
# COLUMNS2="\${string2}\${COLUMNS}"
# echo \${COLUMNS2}

# PR="{print \${COLUMNS2}}"
# echo "\${PR}"

# BASE_NAME=\$(echo  \${fam} | sed 's/.fam//g')
# echo \${BASE_NAME}

# echo \${BASE_NAME}_\${NAME}

# mv \${fam} \$BASE_NAME-old.fam
# cp \${BASE_NAME}-old.fam \${fam}
# cp \${BASE_NAME}.bed \${BASE_NAME}_\${NAME}.bed
# cp \${BASE_NAME}.bim \${BASE_NAME}_\${NAME}.bim
# cp \${BASE_NAME}.log \${BASE_NAME}_\${NAME}.log
# cp \${BASE_NAME}.nosex \${BASE_NAME}_\${NAME}.nosex

# awk -f <(echo "\${PR}") $BASE_DIR/outputs/7_gxp/$DATASET/pheno_table.fam >  $BASE_DIR/outputs/7_gxp/$DATASET/pheno_header_\${NAME}.fam
# sed '1d' $BASE_DIR/outputs/7_gxp/$DATASET/pheno_header_\${NAME}.fam > \${BASE_NAME}_\${NAME}.fam

#   # 1) create relatedness matrix of samples using gemma
# gemma -bfile \${BASE_NAME}_\${NAME} -gk 1 -o /$DATASET/\${NAME}

#   # 2) fit multivariate linear mixed model using gemma (-lmm)
# gemma -bfile \${BASE_NAME}_\${NAME} -k output/$DATASET/\${NAME}.cXX.txt -lmm 4 -n \${C[@]} -o /$DATASET/\${NAME}.lmm

#   # 3) reformat output
# sed 's/\\trs\\t/\\tCHROM\\tPOS\\t/g; s/\\([0-2][0-9]\\):/\\1\\t/g' output/$DATASET/\${NAME}.lmm.assoc.txt | \
#       cut -f 2,3,8-10,13-15 | body sort -k1,1 -k2,2n | gzip > $BASE_DIR/outputs/7_gxp/$DATASET/\${NAME}.lmm.GxP.txt.gz

# rm \${BASE_NAME}_\${NAME}.bed
# rm \${BASE_NAME}_\${NAME}.bim
# rm \${BASE_NAME}_\${NAME}.log
# rm \${BASE_NAME}_\${NAME}.nosex


# EOA



# ------------------------------------------------------------------------------
# Job 5 phen file preparation for mv-plink

jobfile5=5_phen.tmp # temp file
cat > $jobfile5 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=5_phen.tmp
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/5_phen_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/5_phen_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=22G
#SBATCH --time=04:00:00


MAP=$BASE_DIR/outputs/7_gxp/$DATASET/GxP_plink.map 
PED=$BASE_DIR/outputs/7_gxp/$DATASET/GxP_plink.ped
cp \${MAP} $BASE_DIR/outputs/7_gxp/$DATASET/gwas_multi.map
cp \${PED} $BASE_DIR/outputs/7_gxp/$DATASET/gwas_multi.ped

P=$BASE_DIR/outputs/7_gxp/$DATASET/pheno_table.fam
awk -F " " '{print \$1,\$2,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15}' \${P} > $BASE_DIR/outputs/7_gxp/$DATASET/gwas_multi.phen

var="FID IID PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10"
sed -i "1s/.*/\${var}/" $BASE_DIR/outputs/7_gxp/$DATASET/gwas_multi.phen


EOA



# ------------------------------------------------------------------------------
# Job 6 run mvplink

jobfile6=6_mvplink.tmp # temp file
cat > $jobfile6 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=6_mvplink.tmp
#SBATCH --partition=carl.p
#SBATCH --array=1-55
#SBATCH --output=$BASE_DIR/logs/6_mvplink_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/6_mvplink_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90G
#SBATCH --time=1-00:00:00


INPUT_MV=$BASE_DIR/outputs/lof/mvplink.fofn

#Create a job for all the possible phenotypes and the associated .fam file with just one phenotype at a time
MV=\$(cat \${INPUT_MV} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \${MV}

NAME="\$(cut -d ' ' -f 1 <<<"\${MV}")"
echo \${NAME}

COL="\$(cut -d ' ' -f 2- <<<"\${MV}")"
echo \${COL}

P=\$(echo  /user/doau0129/work/software/plink.multivariate --file gwas_multi --mqfam --mult-pheno gwas_multi.phen --pheno-name \${COL})
echo \${P}

/user/doau0129/work/software/plink.multivariate --file $BASE_DIR/outputs/7_gxp/$DATASET/gwas_multi --mqfam --mult-pheno $BASE_DIR/outputs/7_gxp/$DATASET/gwas_multi.phen --pheno-name \${COL} --out $BASE_DIR/outputs/7_gxp/$DATASET/\${NAME}.plink.mqfam.total


EOA



# ------------------------------------------------------------------------------
# Job 7 format mvplink results

jobfile7=7_slider.tmp # temp file
cat > $jobfile7 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=7_slider.tmp
#SBATCH --partition=carl.p
#SBATCH --array=1-55
#SBATCH --output=$BASE_DIR/logs/7_slider_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/7_slider_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90G
#SBATCH --time=1-00:00:00


body() {
        IFS= read -r header
        printf '%s\n' "\$header"
        "\$@"
}


INPUT_MV=$BASE_DIR/outputs/lof/mvplink.fofn

#Create a job for all the possible phenotypes and the associated .fam file with just one phenotype at a time
MV=\$(cat \${INPUT_MV} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \${MV}

NAME="\$(cut -d ' ' -f 1 <<<"\${MV}")"
echo \${NAME}

FILE=$BASE_DIR/outputs/7_gxp/$DATASET/\${NAME}.plink.mqfam.total.mqfam.total
echo \${FILE}

awk '{print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8}' \${FILE} | cut -d ' ' -f 2,3,6-9 | sed 's/SNP BP /CHROM POS /g' | awk '{sub(/\:.*$/,"",\$1); print \$0}' | awk '{if (\$3!="NA"){ print}}' | body sort -k1,1 -k2,2n | gzip > $BASE_DIR/outputs/7_gxp/$DATASET/\${NAME}.mvplink.txt.gz

win5=50000
step5=5000
echo \${win5}
echo \${step5}

win1=10000
step1=1000
echo \${win1}
echo \${step1}

$BASE_DIR/sh/mvplink_slider.sh $BASE_DIR/outputs/7_gxp/$DATASET/\${NAME}.mvplink.txt.gz \${win5} \${step5}
$BASE_DIR/sh/mvplink_slider.sh $BASE_DIR/outputs/7_gxp/$DATASET/\${NAME}.mvplink.txt.gz \${win1} \${step1}

$BASE_DIR/sh/mvplink_log.sh $BASE_DIR/outputs/7_gxp/$DATASET/\${NAME}.mvplink.txt.gz


EOA



# ------------------------------------------------------------------------------
# Job 8 GxP plots

jobfile8=8_plots.tmp # temp file
cat > $jobfile8 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=8_plots.tmp
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/8_plots_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/8_plots_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=22G
#SBATCH --time=04:00:00


ml hpc-env/8.3
ml R/4.0.2-foss-2019b
ml FriBidi
ml HarfBuzz

ls -1 $BASE_DIR/outputs/7_gxp/$DATASET/*.mvplink.50k.5k.txt.gz > $BASE_DIR/outputs/lof/mvplink_50k.fofn

mkdir $BASE_DIR/figures/7_gxp/
mkdir $BASE_DIR/figures/7_gxp/$DATASET/

Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$DATASET/ $BASE_DIR/figures/7_gxp/$DATASET/ univariate_gemma
Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$DATASET/ $BASE_DIR/figures/7_gxp/$DATASET/ multivariate_plink_PC1
Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$DATASET/ $BASE_DIR/figures/7_gxp/$DATASET/ multivariate_plink_PC2
Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$DATASET/ $BASE_DIR/figures/7_gxp/$DATASET/ multivariate_plink_PC3
Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$DATASET/ $BASE_DIR/figures/7_gxp/$DATASET/ multivariate_plink_PC4
Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$DATASET/ $BASE_DIR/figures/7_gxp/$DATASET/ multivariate_plink_PC5
Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$DATASET/ $BASE_DIR/figures/7_gxp/$DATASET/ multivariate_plink_PC6
Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$DATASET/ $BASE_DIR/figures/7_gxp/$DATASET/ multivariate_plink_PC7
Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$DATASET/ $BASE_DIR/figures/7_gxp/$DATASET/ multivariate_plink_PC8
Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$DATASET/ $BASE_DIR/figures/7_gxp/$DATASET/ multivariate_plink_PC9
Rscript $BASE_DIR/R/gxp_plots.R $BASE_DIR/outputs/7_gxp/$DATASET/ $BASE_DIR/figures/7_gxp/$DATASET/ multivariate_plink_byPCs


EOA

# ------------------------------------------------------------------------------
# Job 9 other plots

jobfile9=9_other_plots.tmp # temp file
cat > $jobfile9 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=9_other_plots.tmp
#SBATCH --partition=carl.p
#SBATCH --array=1-55
#SBATCH --output=$BASE_DIR/logs/9_other_plots_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/9_other_plots_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=22G
#SBATCH --time=04:00:00


ml hpc-env/8.3
ml R-bundle-Bioconductor/3.12-foss-2019b-R-4.0.2
ml FriBidi
ml HarfBuzz


# ls -1 $BASE_DIR/outputs/7_gxp/$DATASET/*.mvplink.50k.5k.txt.gz > $BASE_DIR/outputs/lof/mvplink_50k.fofn

INPUT_AVG=$BASE_DIR/outputs/lof/mvplink_50k.fofn

#Create a job for all the possible phenotypes and the associated .fam file with just one phenotype at a time
AVG=\$(cat \${INPUT_AVG} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \${AVG}

B=\$(basename "\${AVG}")
echo \${B}
NAME=\${B%%.*}
echo \${NAME}
EFF=${DATASET%%_*}
echo \${EFF}

T="1.7"
echo \${T}

mkdir $BASE_DIR/figures/7_gxp/$DATASET/\${T}/
mkdir $BASE_DIR/figures/7_gxp/$DATASET/\${T}/\${NAME}/
mkdir $BASE_DIR/outputs/7_gxp/$DATASET/\${T}/

# Rscript $BASE_DIR/R/gxp_zooms.R $BASE_DIR/outputs/7_gxp/$DATASET/ \${B} $BASE_DIR/figures/7_gxp/$DATASET/\${T}/\${NAME}/ \${NAME} \${T}
python3 $BASE_DIR/python/plot_snp_heatmap.py $BASE_DIR/outputs/7_gxp/$DATASET/\${T}/ \${NAME}.snp.txt $BASE_DIR/images/$DATASET/${DATASET}_modifiedImage.csv $BASE_DIR/ressources/full_mask.tif $BASE_DIR/figures/7_gxp/$DATASET/\${T}/\${NAME}/ \${EFF}


EOA


# ********** Schedule the job launching ***********
# -------------------------------------------------

if [ "$JID_RES" = "jid1" ] || [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ];
then
  echo "*****   0_convert_plink : DONE         **"
else
  jid0=$(sbatch ${jobfile0})
fi


if [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ];
then
  echo "*****   1_plink_binary  : DONE         **"
elif [ "$JID_RES" = jid1 ]
then
  jid1=$(sbatch ${jobfile1})
else
  jid1=$(sbatch --dependency=afterok:${jid0##* } ${jobfile1})
fi


if [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ];
then
  echo "*****   2_prep          : DONE         **"
elif [ "$JID_RES" = "jid2" ]
then
  jid2=$(sbatch ${jobfile2})
else
  jid2=$(sbatch --dependency=afterok:${jid1##* } ${jobfile2})
fi


if [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ];
then
  echo "*****   3_gemma         : DONE         **"
elif [ "$JID_RES" = "jid3" ]
then
  jid3=$(sbatch ${jobfile3})
else
  jid3=$(sbatch --dependency=afterok:${jid2##* } ${jobfile3})
fi


if [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ];
then
  echo "*****   4_windows       : DONE         **"
elif [ "$JID_RES" = "jid4" ]
then
  jid4=$(sbatch ${jobfile4})
else
  jid4=$(sbatch --dependency=afterok:${jid3##* } ${jobfile4})
fi


if [ "$JID_RES" = "jid6" ] || [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ];
then
  echo "*****   5_phen          : DONE         **"
elif [ "$JID_RES" = "jid5" ]
then
  jid5=$(sbatch ${jobfile5})
else
  jid5=$(sbatch --dependency=afterok:${jid4##* } ${jobfile5})
fi
 

if [ "$JID_RES" = "jid7" ] || [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ];
then
  echo "*****   6_mvplink       : DONE         **"
elif [ "$JID_RES" = "jid6" ]
then
  jid6=$(sbatch ${jobfile6})
else
  jid6=$(sbatch --dependency=afterok:${jid5##* } ${jobfile6})
fi
 

if [ "$JID_RES" = "jid8" ] || [ "$JID_RES" = "jid9" ];
then
  echo "*****   7_slider        : DONE         **"
elif [ "$JID_RES" = "jid7" ]
then
  jid7=$(sbatch ${jobfile7})
else
  jid7=$(sbatch --dependency=afterok:${jid6##* } ${jobfile7})
fi


if [ "$JID_RES" = "jid9" ];
then
  echo "*****   8_plots         : DONE         **"
elif [ "$JID_RES" = "jid8" ]
then
  jid8=$(sbatch ${jobfile8})
else
  jid8=$(sbatch --dependency=afterok:${jid7##* } ${jobfile8})
fi


if [ "$JID_RES" = "jid9" ];
then
  jid9=$(sbatch ${jobfile9})
else
  jid9=$(sbatch --dependency=afterok:${jid8##* } ${jobfile9})
fi



# ******** Remove useless files and folders *******
# -------------------------------------------------

rm *tmp
