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

tr=(nig, pue)
printf "%s\n" "\${tr[@]}" > $BASE_DIR/outputs/listoffiles/fst_species.fofn

INPUT_SP=$BASE_DIR/outputs/listoffiles/fst_species.fofn
SP=\$(cat \${INPUT_SP} | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
echo \${SP}

if [ "\${SP}" = "nig" ];
then
  vcfsamplenames \${INPUT_VCF} | \
       grep \${SP} > \${SP}.pop
fi



if [ "\${SP}" = "nig" ];
then
  vcfsamplenames \${INPUT_VCF} | \
        grep \${SP} | \
        grep -v unipue | \
        grep -v tanpue | \
        grep -v indpue | \
        grep -v gutpue | \
        grep -v chlpue | \
        grep -v abe pue | \
        grep -v nigpue > \${SP}.pop
fi


vcftools --gzvcf \${INPUT_VCF} \
      --keep \${SP}.pop \
      --mac 3 \
      --recode \
      --stdout | gzip > \${SP}.vcf.gz

EOA
