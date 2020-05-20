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

vcftools \
      --gzvcf \${INPUT_VCF} \
      --plink \
      --out $BASE_DIR/outputs/gxp/GxP_plink

  plink \
      --file $BASE_DIR/outputs/gxp/GxP_plink \
      --recode12 \
      --out $BASE_DIR/outputs/gxp/hapmap

EOA

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


EOA

jobfile3=12_gemma.tmp # temp file
cat > $jobfile3 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=12_gemma.tmp
#SBATCH --partition=carl.p
#SBATCH --array=1-16
#SBATCH --output=$BASE_DIR/logs/12_gemma_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/12_gemma_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=04:30:00

body() {
	IFS= read -r header
	printf '%s\n' "$header"
	"$@"
}

fam = $BASE_DIR/outputs/gxp/GxP_plink_binary.fam
pheno = $BASE_DIR/metadata/traits

tr=(bars_head bars_body snout peduncle H.gummiguta H.unicolor H.puella H.nigricans H.indigo Tan.hamlet H.chlorurus H.guttavarius H.aberrans H.maya H.gemma H.floridae)
printf "%s\n" "\${tr[@]}" > $BASE_DIR/outputs/listoffiles/traits.fofn

INPUT_TR=$BASE_DIR/outputs/listoffiles/traits.fofn

BASE_NAME=\$(echo  ${fam} | sed 's/.fam//g')

mv ${fam} $BASE_DIR/outputs/gxp/\${BASE_NAME}-old.fam
cp $BASE_DIR/outputs/gxp/\${BASE_NAME}-old.fam ${fam}

#Create joint phenotype and .fam file with all phenotypes
pheno2 = awk -F ";" '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21}' pheno
pheno_fam = join ${fam} ${pheno2}
pheno_fam_table = awk -F " " '{print $1,$2,$3,$4,$5,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26}' ${pheno_fam}
echo -e 'label Within_family_ID ID_father ID_mother Sex bars_head bars_body snout peduncle H.gummiguta H.unicolor H.puella H.nigricans H.indigo Tan.hamlet H.chlorurus H.guttavarius H.aberrans H.maya H.gemma H.floridae' > $BASE_DIR/outputs/gxp/pheno_table.fam && cat ${pheno_fam_table} >> $BASE_DIR/outputs/gxp/pheno_table.fam


#Create a job for all the possible phenotypes and the associated .fam file with just one phenotype at a time
TRAITS=$(cat ${INPUT_TR} | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo $TRAITS

awk -v t="$TRAITS" 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} {print \$(f["label"]), \$(f["Within_family_ID"]), \$(f["ID_father"]), \$(f["ID_mother"]), \$(f["Sex"]), \$(f[t])}' $BASE_DIR/outputs/gxp/pheno_table.fam > ${fam}


  # 2) create relatedness matrix of samples using gemma
gemma -bfile \$BASE_NAME -gk 1 -o ${TRAITS}

  # 3) fit linear model using gemma (-lm)
gemma -bfile \$BASE_NAME -lm 4 -miss 0.1 -notsnp -o ${TRAITS}.lm

  # 4) fit linear mixed model using gemma (-lmm)
gemma -bfile \$BASE_NAME -k output/${TRAITS}.cXX.txt -lmm 4 -o ${TRAITS}.lmm

  # 5) reformat output
sed 's/\\trs\\t/\\tCHROM\\tPOS\\t/g; s/\\([0-2][0-9]\\):/\\1\\t/g' output/${TRAITS}.lm.assoc.txt | \
      cut -f 2,3,9-14 | body sort -k1,1 -k2,2n | gzip > ${TRAITS}.lm.GxP.txt.gz
sed 's/\\trs\\t/\\tCHROM\\tPOS\\t/g; s/\\([0-2][0-9]\\):/\\1\\t/g' output/${TRAITS}.lmm.assoc.txt | \
      cut -f 2,3,8-10,13-15 | body sort -k1,1 -k2,2n | gzip > ${TRAITS}.lmm.GxP.txt.gz

EOA



if [ "$JID_RES" = "jid2" ];
then
  echo "10_convert_plink DONE                   **"
else
  jid1=$(sbatch ${jobfile1})
fi

if [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ];
then
  echo "*****     11_plink_binary DONE          **"
elif [ "$JID_RES" = jid2 ]
then
  jid2=$(sbatch ${jobfile2})
else
  jid2=$(sbatch --dependency=afterok:${jid1##* } ${jobfile2})
fi

if [ "$JID_RES" = "jid4" ];
then
  echo "*****     12_gemma DONE              **"
elif [ "$JID_RES" = "jid3" ]
then
  jid3=$(sbatch ${jobfile3})
else
  jid3=$(sbatch --dependency=afterok:${jid2##* } ${jobfile3})
fi
