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

awk '{print \$1"\\t"\$2\$3}' \${INPUT_POP} > $BASE_DIR/outputs/fst/pop.txt

for k in abebel abeboc abepue chlpue floflo gemflo gumboc gutpue indbel indpue maybel nigbel nigboc nigflo nigpue puebel pueboc pueflo puepue ranbel tanpue uniboc uniflo unipue; do
  grep \$k pop.txt | cut -f 1 > $BASE_DIR/outputs/fst/pop.\$k.txt
  done


POP="--weir-fst-pop $BASE_DIR/outputs/fst/pop.abebel.txt \
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
vcftools --gzvcf \${INPUT_VCF} \
      \$POP \
     --stdout  2> multi_fst_snp.log | \
     gzip > $BASE_DIR/outputs/fst/multi_fst_snp.tsv.gz

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
     gzip > $BASE_DIR/outputs/fst/multi_fst.10k.tsv.gz

EOA


jobfile3=23_prep_pairwise.tmp # temp file
cat > $jobfile3 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=23_prep_pairwise
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/23_prep_pairwise_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/23_prep_pairwise_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100M
#SBATCH --time=00:00:55


printf 'nig\n%.0s' {1..6} > col1
col2=(flo flo flo boc boc bel)
paste -d " " col1 \${col2} > col2
col3=(boc bel pue bel pue pue)
paste -d " " col2 \${col3} > nig

printf 'pue\n%.0s' {1..6} > col11
col22=(flo flo flo boc boc bel)
paste -d " " col11 \${col22} > col22
col33=(boc bel pue bel pue pue)
paste -d " " col22 \${col33} > pue

awk 'NF' nig pue > $BASE_DIR/outputs/listoffiles/fst_pairwise.fofn
rm nig
rm pue
rm col*


EOA


jobfile4=24_pairwise_fst.tmp # temp file
cat > $jobfile4 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=24_pairwise_fst
#SBATCH --partition=carl.p
#SBATCH --array=1-12
#SBATCH --output=$BASE_DIR/logs/24_pairwise_fst_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/24_pairwise_fst_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100M
#SBATCH --time=00:00:55

INPUT_PW=$BASE_DIR/outputs/listoffiles/fst_pairwise.fofn
PW=\$(cat \${INPUT_PW} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \${PW}

SP=\$(cut -d' ' -f1 <<< \$PW)
LOC1=\$(cut -d' ' -f2 <<< \$PW)
LOC2=\$(cut -d' ' -f3 <<< \$PW)
echo \${SP}
echo \${LOC1}
echo \${LOC2}

POP=$BASE_DIR/outputs/fst/\${SP}.pop
echo \${POP}

VCF=$BASE_DIR/outputs/fst/\${SP}.vcf.gz
echo \${VCF}


grep \${LOC1} \${POP} > $BASE_DIR/outputs/fst/\${SP}_\${LOC1}_\${LOC2}_pop1.txt
grep \${LOC2} \${POP} > $BASE_DIR/outputs/fst/\${SP}_\${LOC1}_\${LOC2}_pop2.txt

vcftools --gzvcf ${VCF} \
      --weir-fst-pop $BASE_DIR/outputs/fst/\${SP}_\${LOC1}_\${LOC2}_pop1.txt \
      --weir-fst-pop $BASE_DIR/outputs/fst/\${SP}_\${LOC1}_\${LOC2}_pop2.txt \
      --fst-window-step 5000 \
      --fst-window-size 50000 \
      --stdout 2> $BASE_DIR/outputs/fst/\${SP}_\${LOC1}_\${LOC2}.50k.log | \
      gzip > $BASE_DIR/outputs/fst/\${SP}_\${LOC1}_\${LOC2}.50k.windowed.weir.fst

  vcftools --gzvcf ${vcf} \
      --weir-fst-pop pop1.txt \
      --weir-fst-pop pop2.txt \
      --fst-window-size 10000 \
      --fst-window-step 1000 \
      --stdout 2> $BASE_DIR/outputs/fst/\${SP}_\${LOC1}_\${LOC2}.10k.log | \
      gzip > $BASE_DIR/outputs/fst/\${SP}_\${LOC1}_\${LOC2}.10k.windowed.weir.fst


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

if [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ];
then
  echo "*****     23_prep_pairwise DONE          **"
elif [ "$JID_RES" = jid3 ]
then
  jid3=$(sbatch ${jobfile3})
else
  jid3=$(sbatch --dependency=afterok:${jid2##* } ${jobfile3})
fi

if [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ];
then
  echo "*****     24_pairwise_fst DONE          **"
elif [ "$JID_RES" = jid4 ]
then
  jid4=$(sbatch ${jobfile4})
else
  jid4=$(sbatch --dependency=afterok:${jid3##* } ${jobfile4})
fi
