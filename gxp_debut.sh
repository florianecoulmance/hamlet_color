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
mkdir $BASE_DIR/outputs/gxp/

jobfile1=10_covert_plink.tmp # temp file
cat > $jobfile1 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=10_convert_plink
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/10_convert_plink_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/10_convert_plink_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=02:30:00

INPUT_VCF=$BASE_DIR/outputs/09_1_snpfiltration/filterd_bi-allelic.vcf.gz


smp=(PL17_35puepue PL17_35indpue)
printf "%s " "\${smp[@]}" > $BASE_DIR/outputs/09_1_snpfiltration/change_sample.txt

bcftools reheader -s $BASE_DIR/outputs/09_1_snpfiltration/change_sample.txt -o $BASE_DIR/outputs/09_1_snpfiltration/filterd_bi-allelic_changed.vcf.gz \${INPUT_VCF}

INPUT_CHANGED=$BASE_DIR/outputs/09_1_snpfiltration/filterd_bi-allelic_changed.vcf.gz

vcftools \
      --gzvcf \${INPUT_CHANGED} \
      --plink \
      --out $BASE_DIR/outputs/gxp/GxP_plink

  plink \
      --file $BASE_DIR/outputs/gxp/GxP_plink \
      --recode12 \
      --out $BASE_DIR/outputs/gxp/hapmap

EOA

#Job2
jobfile2=11_plink_binary.tmp # temp file
cat > $jobfile2 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=11_plink_binary
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/11_plink_binary_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/11_plink_binary_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=02:30:00

plink \
    --noweb \
    --file $BASE_DIR/outputs/gxp/GxP_plink \
    --make-bed \
    --out $BASE_DIR/outputs/gxp/GxP_plink_binary

cp $BASE_DIR/outputs/gxp/GxP_plink_binary.fam $BASE_DIR/outputs/gxp/GxP_plink_binary_sauvegarde.fam

EOA

jobfile3=12_prep_gemmad.tmp # temp file
cat > $jobfile3 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=12_prep_gemmad.tmp
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/12_prep_gemmad_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/12_prep_gemmad_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500M
#SBATCH --time=00:15:00

#cp $BASE_DIR/outputs/gxp/GxP_plink_binary_sauvegarde.fam $BASE_DIR/outputs/gxp/trait3/GxP_plink_binary.fam
#cp $BASE_DIR/outputs/gxp/GxP_plink_binary.bed $BASE_DIR/outputs/gxp/trait3/GxP_plink_binary.bed
#cp $BASE_DIR/outputs/gxp/GxP_plink_binary.bim $BASE_DIR/outputs/gxp/trait3/GxP_plink_binary.bim
#cp $BASE_DIR/outputs/gxp/GxP_plink_binary.log $BASE_DIR/outputs/gxp/trait3/GxP_plink_binary.log
#cp $BASE_DIR/outputs/gxp/GxP_plink_binary.nosex $BASE_DIR/outputs/gxp/trait3/GxP_plink_binary.nosex
 
fam=$BASE_DIR/outputs/gxp/testgxp_debut/GxP_plink_binary.fam
pheno=$BASE_DIR/metadata/testgxp_debut

tr=(bhead bbody snout ped pue nig ind combo_spec)
printf "%s\n" "\${tr[@]}" > $BASE_DIR/outputs/listoffiles/testgxp_debut.fofn

#Create joint phenotype and .fam file with all phenotypes
awk -F ";" '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19}' \${pheno} > $BASE_DIR/outputs/gxp/testgxp_debut/pheno_intermediate1
sort -k1 $BASE_DIR/outputs/gxp/testgxp_debut/pheno_intermediate1 > $BASE_DIR/outputs/gxp/testgxp_debut/pheno_intermediate2
sort -k1 \${fam}
join \${fam} $BASE_DIR/outputs/gxp/testgxp_debut/pheno_intermediate2 > $BASE_DIR/outputs/gxp/testgxp_debut/pheno_intermediate3
awk -F " " '{print \$1,\$2,\$3,\$4,\$5,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$20,\$21,\$22,\$23,\$24}' $BASE_DIR/outputs/gxp/testgxp_debut/pheno_intermediate3 > $BASE_DIR/outputs/gxp/testgxp_debut/pheno_intermediate4
echo -e 'label Within_family_ID ID_father ID_mother Sex bhead bbody snout ped gum uni pue nig ind tan chl gut abe may gem flo ran combo_spec' > $BASE_DIR/outputs/gxp/testgxp_debut/pheno_table.fam && cat $BASE_DIR/outputs/gxp/testgxp_debut/pheno_intermediate4 >> $BASE_DIR/outputs/gxp/testgxp_debut/pheno_table.fam


EOA


jobfile4=13_gemmad.tmp # temp file
cat > $jobfile4 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=13_gemmad.tmp
#SBATCH --partition=carl.p
#SBATCH --array=1-8
#SBATCH --output=$BASE_DIR/logs/13_gemmad_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/13_gemmad_%A_%a.err
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


fam=$BASE_DIR/outputs/gxp/testgxp_debut/GxP_plink_binary.fam
INPUT_TR=$BASE_DIR/outputs/listoffiles/testgxp_debut.fofn

#Create a job for all the possible phenotypes and the associated .fam file with just one phenotype at a time
TRAITS=\$(cat \${INPUT_TR} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \${TRAITS}


BASE_NAME=\$(echo  \${fam} | sed 's/.fam//g')
echo \${BASE_NAME}

echo \${BASE_NAME}_\${TRAITS}

mv \${fam} \$BASE_NAME-old.fam
cp \${BASE_NAME}-old.fam \${fam}
cp \${BASE_NAME}.bed \${BASE_NAME}_\${TRAITS}.bed
cp \${BASE_NAME}.bim \${BASE_NAME}_\${TRAITS}.bim
cp \${BASE_NAME}.log \${BASE_NAME}_\${TRAITS}.log
cp \${BASE_NAME}.nosex \${BASE_NAME}_\${TRAITS}.nosex

awk -v t="\${TRAITS}" 'NR==1 {for (i=1; i<=NF; i++) {f[\$i] = i}} {print \$(f["label"]), \$(f["Within_family_ID"]), \$(f["ID_father"]), \$(f["ID_mother"]), \$(f["Sex"]), \$(f[t])}' $BASE_DIR/outputs/gxp/testgxp_debut/pheno_table.fam >  $BASE_DIR/outputs/gxp/testgxp_debut/pheno_header_\${TRAITS}.fam
sed '1d' $BASE_DIR/outputs/gxp/testgxp_debut/pheno_header_\${TRAITS}.fam > \${BASE_NAME}_\${TRAITS}.fam

  # 2) create relatedness matrix of samples using gemma
gemma -bfile \${BASE_NAME}_\${TRAITS} -gk 1 -o /testgxp_debut/\${TRAITS}

  # 3) fit linear model using gemma (-lm)
gemma -bfile \${BASE_NAME}_\${TRAITS} -lm 4 -miss 0.1 -notsnp -o /testgxp_debut/\${TRAITS}.lm

  # 4) fit linear mixed model using gemma (-lmm)
gemma -bfile \${BASE_NAME}_\${TRAITS} -k output/testgxp_debut/\${TRAITS}.cXX.txt -lmm 4 -o /testgxp_debut/\${TRAITS}.lmm

  # 5) reformat output
sed 's/\\trs\\t/\\tCHROM\\tPOS\\t/g; s/\\([0-2][0-9]\\):/\\1\\t/g' output/testgxp_debut/\${TRAITS}.lm.assoc.txt | \
      cut -f 2,3,9-14 | body sort -k1,1 -k2,2n | gzip > $BASE_DIR/outputs/gxp/testgxp_debut/\${TRAITS}.lm.GxP.txt.gz
sed 's/\\trs\\t/\\tCHROM\\tPOS\\t/g; s/\\([0-2][0-9]\\):/\\1\\t/g' output/testgxp_debut/\${TRAITS}.lmm.assoc.txt | \
      cut -f 2,3,8-10,13-15 | body sort -k1,1 -k2,2n | gzip > $BASE_DIR/outputs/gxp/testgxp_debut/\${TRAITS}.lmm.GxP.txt.gz


rm \${BASE_NAME}_\${TRAITS}.bed
rm \${BASE_NAME}_\${TRAITS}.bim
rm \${BASE_NAME}_\${TRAITS}.log

EOA

jobfile5=14_windowsd.tmp # temp file
cat > $jobfile5 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=14_windowsd.tmp
#SBATCH --partition=carl.p
#SBATCH --array=1-8
#SBATCH --output=$BASE_DIR/logs/14_windowsd_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/14_windowsd_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=22G
#SBATCH --time=02:30:00

INPUT_TR=$BASE_DIR/outputs/listoffiles/testgxp_debut.fofn

#Create a job for all the possible phenotypes and the associated .fam file with just one phenotype at a time
TRAITS=\$(cat \${INPUT_TR} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \${TRAITS}

lm=$BASE_DIR/outputs/gxp/testgxp_debut/\${TRAITS}.lm.GxP.txt.gz
echo \${lm}
lmm=$BASE_DIR/outputs/gxp/testgxp_debut/\${TRAITS}.lmm.GxP.txt.gz
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


if [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ];
then
  echo "10_convert_plink DONE                   **"
else
  jid1=$(sbatch ${jobfile1})
fi

if [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ];
then
  echo "*****     11_plink_binary DONE          **"
elif [ "$JID_RES" = jid2 ]
then
  jid2=$(sbatch ${jobfile2})
else
  jid2=$(sbatch --dependency=afterok:${jid1##* } ${jobfile2})
fi

if [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ];
then
  echo "*****     12_prep_gemma DONE              **"
elif [ "$JID_RES" = "jid3" ]
then
  jid3=$(sbatch ${jobfile3})
else
  jid3=$(sbatch --dependency=afterok:${jid2##* } ${jobfile3})
fi

if [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ];
then
  echo "*****     13_gemma DONE              **"
elif [ "$JID_RES" = "jid4" ]
then
  jid4=$(sbatch ${jobfile4})
else
  jid4=$(sbatch --dependency=afterok:${jid3##* } ${jobfile4})
fi

if [ "$JID_RES" = "jid6" ];
then
  echo "*****     14_windows DONE              **"
elif [ "$JID_RES" = "jid5" ]
then
  jid5=$(sbatch ${jobfile5})
else
  jid5=$(sbatch --dependency=afterok:${jid4##* } ${jobfile5})
fi
