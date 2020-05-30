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
#SBATCH --array=1-5
#SBATCH --output=$BASE_DIR/logs/20_keep_species_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/20_keep_species_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=02:30:00

#source activate /user/doau0129/miniconda3/envs/vcfenv

#INPUT_VCF=$BASE_DIR/outputs/09_1_snpfiltration/filterd_bi-allelic.vcf.gz
#echo \${INPUT_VCF}

#smp=(PL17_35puepue PL17_35indpue)
#printf "%s " "\${smp[@]}" > $BASE_DIR/outputs/09_1_snpfiltration/change_sample.txt

#bcftools reheader -s $BASE_DIR/outputs/09_1_snpfiltration/change_sample.txt -o $BASE_DIR/outputs/09_1_snpfiltration/filterd_bi-allelic_changed.vcf.gz \${INPUT_VCF}
#tabix filterd_bi-allelic_changed.vcf.gz

INPUT_CHANGED=$BASE_DIR/outputs/09_1_snpfiltration/filterd_bi-allelic_changed.vcf.gz


tr=(nig pue bel boc puer)
printf "%s\n" "\${tr[@]}" > $BASE_DIR/outputs/listoffiles/fst_species.fofn

INPUT_SP=$BASE_DIR/outputs/listoffiles/fst_species.fofn
SP=\$(cat \${INPUT_SP} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \${SP}

if [ \${SP} = "nig" ];
then
  echo \${SP}
  vcfsamplenames \${INPUT_CHANGED} | grep nig > $BASE_DIR/outputs/fst/\${SP}.pop
fi



if [ \${SP} = "pue" ];
then
  echo \${SP}
  vcfsamplenames \${INPUT_CHANGED} | grep pue | grep -v abe | grep -v chl | grep -v gem | grep -v gum | grep -v gut | grep -v ind | grep -v may | grep -v nig | grep -v ran | grep -v tan | grep -v uni > $BASE_DIR/outputs/fst/\${SP}.pop
fi

if [ \${SP} = "bel" ];
then
  echo \${SP}
  vcfsamplenames \${INPUT_CHANGED} | grep bel | grep -v abe | grep -v chl | grep -v flo | grep -v gem | grep -v gum | grep -v gut | grep -v ran | grep -v tan | grep -v uni > $BASE_DIR/outputs/fst/\${SP}.pop
fi


if [ \${SP} = "boc" ];
then
  echo \${SP}
  vcfsamplenames \${INPUT_CHANGED} | grep boc | grep -v abe | grep -v chl | grep -v flo | grep -v gem | grep -v gum | grep -v gut | grep -v ran | grep -v tan | grep -v ind | grep -v may > $BASE_DIR/outputs/fst/\${SP}.pop
fi

if [ \${SP} = "puer" ];
then
  echo \${SP}
  vcfsamplenames \${INPUT_CHANGED} | grep pue | grep -v abe | grep -v nig | grep -v flo | grep -v gem | grep -v gum | grep -v gut | grep -v ran | grep -v ind | grep -v may | grep -v bel | grep -v boc  > $BASE_DIR/outputs/fst/\${SP}.pop
fi

vcftools --gzvcf \${INPUT_CHANGED} \
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

#source activate /user/doau0129/miniconda3/envs/vcfenv

INPUT_VCF=$BASE_DIR/outputs/09_1_snpfiltration/filterd_bi-allelic_changed.vcf.gz
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

#source activate /user/doau0129/miniconda3/envs/vcfenv

INPUT_VCF=$BASE_DIR/outputs/09_1_snpfiltration/filterd_bi-allelic_changed.vcf.gz
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
printf "%s\n" "\${col2[@]}" > col2
paste -d " " col1 col2 > col12
col3=(boc bel pue bel pue pue)
printf "%s\n" "\${col3[@]}" > col3
paste -d " " col12 col3 > nig

printf 'pue\n%.0s' {1..6} > col11
col22=(flo flo flo boc boc bel)
printf "%s\n" "\${col22[@]}" > col22
paste -d " " col11 col22 > col1122
col33=(boc bel pue bel pue pue)
printf "%s\n" "\${col33[@]}" > col33
paste -d " " col1122 col33 > pue


printf 'bel\n%.0s' {1..6} > col111
col222=(ind ind ind may may nig)
printf "%s\n" "\${col222[@]}" > col222
paste -d " " col111 col222 > col111222
col333=(may nig pue nig pue pue)
printf "%s\n" "\${col333[@]}" > col333
paste -d " " col111222 col333 > bel


printf 'boc\n%.0s' {1..3} > col1111
col2222=(nig nig pue)
printf "%s\n" "\${col2222[@]}" > col2222
paste -d " " col1111 col2222 > col11112222
col3333=(pue uni uni)
printf "%s\n" "\${col3333[@]}" > col3333
paste -d " " col11112222 col3333 > boc


printf 'puer\n%.0s' {1..3} > col11111
col22222=(chl chl pue)
printf "%s\n" "\${col22222[@]}" > col22222
paste -d " " col11111 col22222 > col1111122222
col33333=(pue uni uni)
printf "%s\n" "\${col33333[@]}" > col33333
paste -d " " col1111122222 col33333 > puer

awk 'NF' nig pue bel boc puer > $BASE_DIR/outputs/listoffiles/fst_pairwise.fofn
rm nig
rm pue
rm bel
rm boc
rm puer
rm col*


EOA


jobfile4=24_pairwise_fst.tmp # temp file
cat > $jobfile4 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=24_pairwise_fst
#SBATCH --partition=carl.p
#SBATCH --array=1-24
#SBATCH --output=$BASE_DIR/logs/24_pairwise_fst_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/24_pairwise_fst_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=04:30:00


#source activate /user/doau0129/miniconda3/envs/vcfenv

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

vcftools --gzvcf \${VCF} \
      --weir-fst-pop $BASE_DIR/outputs/fst/\${SP}_\${LOC1}_\${LOC2}_pop1.txt \
      --weir-fst-pop $BASE_DIR/outputs/fst/\${SP}_\${LOC1}_\${LOC2}_pop2.txt \
      --fst-window-step 5000 \
      --fst-window-size 50000 \
      --out $BASE_DIR/outputs/fst/\${SP}_\${LOC1}_\${LOC2}.50k 2> $BASE_DIR/outputs/fst/\${SP}_\${LOC1}_\${LOC2}.50k.log

vcftools --gzvcf \${VCF} \
      --weir-fst-pop $BASE_DIR/outputs/fst/\${SP}_\${LOC1}_\${LOC2}_pop1.txt \
      --weir-fst-pop $BASE_DIR/outputs/fst/\${SP}_\${LOC1}_\${LOC2}_pop2.txt \
      --fst-window-size 10000 \
      --fst-window-step 1000 \
      --out $BASE_DIR/outputs/fst/\${SP}_\${LOC1}_\${LOC2}.10k 2> $BASE_DIR/outputs/fst/\${SP}_\${LOC1}_\${LOC2}.10k.log

gzip $BASE_DIR/outputs/fst/\${SP}_\${LOC1}_\${LOC2}.50k.windowed.weir.fst
gzip $BASE_DIR/outputs/fst/\${SP}_\${LOC1}_\${LOC2}.10k.windowed.weir.fst

EOA


jobfile5=25_global.tmp # temp file
cat > $jobfile5 <<EOA # generate the job file
#!/bin/bash

#SBATCH --job-name=25_global
#SBATCH --partition=carl.p
#SBATCH --output=$BASE_DIR/logs/25_global_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/25_global_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00

ls -1 $BASE_DIR/outputs/fst/*.50k.log > $BASE_DIR/outputs/listoffiles/25_logs50k.fofn 
cd $BASE_DIR/outputs/fst/

cat *.50k.log | \
    grep -E 'Weir and Cockerham|--out' | \
    grep -A 3 50k | \
    sed '/^--/d; s/^.*--out //g; s/.50k//g; /^Output/d; s/Weir and Cockerham //g; s/ Fst estimate: /\t/g' | \
    paste - - - | \
    cut -f 1,3,5 | \

sed 's/^\\(...\\)-/\\1\\t/g' > $BASE_DIR/outputs/fst/fst_globals.txt

EOA


if [ "$JID_RES" = "jid11" ] || [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ];
then
  echo "20_keep_species DONE **"
elif [ "$JID_RES" = "jid1" ]
then
  jid1=$(sbatch ${jobfile1})
else
  jid1=$(sbatch ${jobfile1})
fi


if [ "$JID_RES" = "jid1" ] || [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ];
then
  echo "21_all_hamlets DONE **"
elif [ "$JID_RES" = "jid11" ]
then
  jid11=$(sbatch ${jobfile11})
else
  jid11=$(sbatch ${jobfile11})
fi


if [ "$JID_RES" = "jid1" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ];
then
  echo "*****     22_multi DONE          **"
elif [ "$JID_RES" = jid2 ]
then
  jid2=$(sbatch ${jobfile2})
else
  jid2=$(sbatch --dependency=afterok:${jid11##* } ${jobfile2})
fi

if [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ] || [ "$JID_RES" = "jid6" ];
then
  echo "*****     23_prep_pairwise DONE          **"
elif [ "$JID_RES" = jid3 ]
then
  jid3=$(sbatch ${jobfile3})
else
  jid3=$(sbatch --dependency=afterok:${jid1##* } ${jobfile3})
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

if [ "$JID_RES" = "jid6" ];
then
  echo "*****     25_global DONE          **"
elif [ "$JID_RES" = jid5 ]
then
  jid5=$(sbatch ${jobfile5})
else
  jid5=$(sbatch --dependency=afterok:${jid4##* } ${jobfile5})
fi
