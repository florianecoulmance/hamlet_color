#!/bin/bash

#SBATCH --job-name=00_readgroups_ubams
#SBATCH --partition=carl.p
#SBATCH --array=2-118
#SBATCH --output=/user/doau0129/work/chapter1_2/logs/00_readgroups_%A_%a.out
#SBATCH --error=/user/doau0129/work/chapter1_2/logs/00_readgroups_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=01:00:00

BASE_DIR=/user/doau0129/work/chapter1_2/
INPUT_META=$BASE_DIR/metadata/metadata_gxp_ben_floridae_complete



LINES=$(cat ${INPUT_META} | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)


IFS=";" read -r -a array <<< "$LINES"
label=${array[1]}
company=${array[7]}
frwdfile=${array[9]}
flowcellidfrwd=${array[15]}
lanefrwd=${array[16]}
revfile=${array[24]}

echo -e "----------------------------"
echo -e "Label:\t${label}\nFwd:\t${frwdfile}\nRev:\t${revfile}"
echo -e "Flowcell:\t${flowcellidfrwd}\nLane:\t${lanefrwd}"
echo -e "Read group:\t${flowcellidfrwd}.${lanefrwd}\nCompany:\t${company}"

gatk --java-options "-Xmx20G" \
    FastqToSam \
    -SM=${label} \
    -F1=$BASE_DIR/data/${frwdfile} \
    -F2=$BASE_DIR/data/${revfile} \
    -O=$BASE_DIR/outputs/00_ubams/${label}.${lanefrwd}.ubam.bam \
    -RG=${label}.${lanefrwd} \
    -LB=${label}".lib1" \
    -PU=${flowcellidfrwd}.${lanefrwd} \
    -PL=Illumina \
    -CN=${company} \
    --TMP_DIR=$BASE_DIR/temp_files
