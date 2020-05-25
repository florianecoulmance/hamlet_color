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
mkdir $BASE_DIR/outputs/fst/

jobfile1=20_keep_species.tmp # temp file
cat > $jobfile1 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=20_keep_species
#SBATCH --partition=carl.p
#SBATCH --array=1-2
#SBATCH --output=$BASE_DIR/logs/20_keep_species_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/20_keep_species_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=02:30:00

INPUT_VCF=$BASE_DIR/outputs/09_1_snpfiltration/filterd_bi-allelic.vcf.gz
echo \${INPUT_VCF}

tr=(nig pue)
printf "%s\n" "\${tr[@]}" > $BASE_DIR/outputs/listoffiles/fst_species.fofn

INPUT_SP=$BASE_DIR/outputs/listoffiles/fst_species.fofn
SP=\$(cat \${INPUT_SP} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \${SP}

if [ "\${SP}" = "nig" ];
then
  vcfsamplenames \${INPUT_VCF} | \
       grep \${SP} > $BASE_DIR/outputs/fst/\${SP}.pop
fi



if [ "\${SP}" = "pue" ];
then
  vcfsamplenames \${INPUT_VCF} | \
        grep \${SP} | \
        grep -v unipue | \
        grep -v tanpue | \
        grep -v indpue | \
        grep -v gutpue | \
        grep -v chlpue | \
        grep -v abe pue | \
        grep -v nigpue > $BASE_DIR/outputs/fst/\${SP}.pop
fi


vcftools --gzvcf \${INPUT_VCF} \
      --keep $BASE_DIR/outputs/fst/\${SP}.pop \
      --mac 3 \
      --recode \
      --stdout | gzip > $BASE_DIR/outputs/fst/\${SP}.vcf.gz

EOA


jobfile11=21_all_hamlets.tmp # temp file
cat > $jobfile11 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=21_all_hamlets
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/21_all_hamlets_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/21_all_hamlets_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=02:30:00

INPUT_VCF=$BASE_DIR/outputs/09_1_snpfiltration/filterd_bi-allelic.vcf.gz
echo \${INPUT_VCF}

vcfsamplenames \${INPUT_VCF} | \
       grep "abe\\|gum\\|ind\\|may\\|nig\\|pue\\|ran\\|uni\\|gem\\|flo\\|chl\\|tan\\|gut" | \
       awk '{print \$1"\\t"\$1}' | \
       sed 's/\\t.*\\(...\\)\\(...\\)\$/\\t\\1\\t\\2/g' > $BASE_DIR/outputs/fst/hamlets.pop.txt

EOA


jobfile2=22_multi.tmp # temp file
cat > $jobfile2 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=22_multi
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/22_multi_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/22_multi_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=15:00:00

INPUT_VCF=$BASE_DIR/outputs/09_1_snpfiltration/filterd_bi-allelic.vcf.gz
echo \${INPUT_VCF}


INPUT_POP=$BASE_DIR/outputs/fst/hamlets.pop.txt
echo \${INPUT_POP}

awk '{print \$1"\\t"\$2\$3}' \${INPUT_POP} > pop.txt

for k in abebel abeboc abepue chlpue floflo gemflo gumboc gutpue indbel indpue maybel nigbel nigboc nigflo nigpue puebel pueboc pueflo puepue ranbel tanpue uniboc uniflo unipue; do
  grep \$k pop.txt | cut -f 1 > pop.\$k.txt
  done


POP="--weir-fst-pop pop.abebel.txt \
   --weir-fst-pop pop.abeboc.txt \
   --weir-fst-pop pop.abepue.txt \
   --weir-fst-pop pop.chlpue.txt \
   --weir-fst-pop pop.floflo.txt \
   --weir-fst-pop pop.gemflo.txt \
   --weir-fst-pop pop.gumboc.txt \
   --weir-fst-pop pop.gutpue.txt \
   --weir-fst-pop pop.indbel.txt \
   --weir-fst-pop pop.indpue.txt \
   --weir-fst-pop pop.maybel.txt \
   --weir-fst-pop pop.nigbel.txt \
   --weir-fst-pop pop.nigboc.txt \
   --weir-fst-pop pop.nigflo.txt \
  --weir-fst-pop pop.nigpue.txt \
  --weir-fst-pop pop.puebel.txt \
  --weir-fst-pop pop.pueboc.txt \
  --weir-fst-pop pop.pueflo.txt \
  --weir-fst-pop pop.puepue.txt \
  --weir-fst-pop pop.ranbel.txt \
  --weir-fst-pop pop.tanpue.txt \
  --weir-fst-pop pop.uniboc.txt \
  --weir-fst-pop pop.uniflo.txt \
  --weir-fst-pop pop.unipue.txt" 


  # fst by SNP
     # ----------
vcftools --gzvcf \${INPUT_VCF} \
      \$POP \
     --stdout  2> multi_fst_snp.log | \
     gzip > $BASE_DIR/outputs/fst/multi_fst.tsv.gz

     # fst 50kb window
     # ---------------
vcftools --gzvcf \${INPUT_VCF} \
      \$POP \
     --fst-window-step 5000 \
     --fst-window-size 50000 \
     --stdout  2> multi_fst.50k.log | \
     gzip > $BASE_DIR/outputs/fst/multi_fst.50k.tsv.gz

     # fst 10kb window
     # ---------------
vcftools --gzvcf \${INPUT_VCF} \
      \$POP \
     --fst-window-step 1000 \
     --fst-window-size 10000 \
     --stdout  2> multi_fst.10k.log | \
     gzip > $BASE_DIR/outputs/fst/multi_fst_snp.tsv.gz  

EOA



if [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ];
then
  echo "20_keep_species DONE **"
else
  jid1=$(sbatch ${jobfile1})
fi


if [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ];
then
  echo "21_all_hamlets DONE **"
else
  jid1=$(sbatch ${jobfile11})
fi


if [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ];
then
  echo "*****     22_multi DONE          **"
elif [ "$JID_RES" = jid2 ]
then
  jid2=$(sbatch ${jobfile2})
else
  jid2=$(sbatch --dependency=afterok:${jid1##* } ${jobfile2})
fi
