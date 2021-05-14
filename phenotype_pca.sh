#!/bin/bash
# by: Floriane Coulmance: 11/05/2021
# usage:
# ./phenotype_pca.sh -i /Users/fco/Desktop/PhD/1_CHAPTER1/ -j AB -k fullm -l 54off_59on -m all -n left
# ------------------------------------------------------------------------------
# 
# 
# 
# ------------------------------------------------------------------------------



# ********** Allow to enter bash options **********
# -------------------------------------------------

while getopts i:j:k:l:m:n: option
do
case "${option}"
in
i) BASE_DIR=${OPTARG};; # get the base directory path
j) COLOR_SPACE=${OPTARG};; # get the color space type
k) MASK=${OPTARG};; # get the mask type
l) DATA=${OPTARG};; # get the dataset type
m) SUB_DATA=${OPTARG};; # get the dataset type
n) SIDE=${OPTARG};; # get the dataset type
esac
done



# ********* Create necessary variables *********
# -------------------------------------------------

OVERALL_DAT="${SIDE}_${DATA}"
echo $OVERALL_DAT

ALIGN_FOLDER="left_54off_59on"
echo $ALIGN_FOLDER

DATASET="${COLOR_SPACE}_${MASK}_${SIDE}_${DATA}"
echo $DATASET

FILES_NAME="${COLOR_SPACE}_${MASK}_${DATA}"
echo $FILES_NAME

if [ "$MASK" = "fullm" ];
then
  MASK_FILE=full_mask.tif
elif [ "$MASK" = "bodym" ];
then
  MASK_FILE=body_mask.tif
else
  echo "choose appropriate option"
fi

echo $MASK_FILE



# ********* Create necessary repositories *********
# -------------------------------------------------

# Repo for pca outputs
mkdir $BASE_DIR/1_GENETICS/chapter1/images/
mkdir $BASE_DIR/1_GENETICS/chapter1/images/$DATASET

# Repo for the corresponding dataset
mkdir $BASE_DIR/1_GENETICS/chapter1/figures/7_gxp/$DATASET/



# ********* Run commands *********
# -------------------------------------------------

/Users/fco/miniconda3/bin/python3 python/image_pca_heatmap.py \
         $BASE_DIR/0_IMAGES/convert_png/$ALIGN_FOLDER/3-registred/Modalities/RGB/$SUB_DATA/ \
         $BASE_DIR/0_IMAGES/convert_png/smallDatasetl1/$MASK_FILE \
         $COLOR_SPACE \
         $BASE_DIR/1_GENETICS/chapter1/images/$DATASET/ \
         $BASE_DIR/1_GENETICS/chapter1/figures/7_gxp/$DATASET/ \
         $MASK \
         $OVERALL_DAT \



Rscript R/pcs_plots.R \
         $BASE_DIR/1_GENETICS/chapter1/images/$DATASET/ \
         $COLOR_SPACE \
         $BASE_DIR/1_GENETICS/chapter1/metadata/ \
         ${DATASET}_PCs.csv \
         ${DATASET}_var.csv \
         $BASE_DIR/1_GENETICS/chapter1/metadata/image_metadata.tsv \
         $MASK \
         $OVERALL_DAT \