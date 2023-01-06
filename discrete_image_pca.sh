#!/bin/bash
# by: Floriane Coulmance: 11/05/2021
# usage:
# sbatch discrete_image_pca.sh -i <PATH:/Users/fco/Desktop/PhD/1_CHAPTER1/> -j RGB -k 7 -l all -m fullm
# sbatch discrete_image_pca.sh -i <PATH> -j <COLOR_SPACE> -k <CLUSTER_NUMBER> -l <SUB_DATA> -m <MASK>
# ------------------------------------------------------------------------------
# PATH corresponds to the path to the base directory, all outputs and necessary
# folder will be created by the script 
# COLOR_SPACE : LAB or RGB
# CLUSTER_NUMBER : 7 or any number before
# SUB_DATA : all, on, off, barred, unbarred
# MASK : bodym (for mask without fins), fullm (for mask including all fins)
# ------------------------------------------------------------------------------



# ********** Allow to enter bash options **********
# -------------------------------------------------

while getopts i:j:k:l:m:n: option
do
case "${option}"
in
i) BASE_DIR=${OPTARG};; # get the base directory path
j) COLOR_SPACE=${OPTARG};; # get the color space type
k) CLUSTER_NUMBER=${OPTARG};; # get the mask type
l) SUB_DATA=${OPTARG};; # get the dataset type
m) MASK=${OPTARG};; # get the dataset type
esac
done



# ********* Create necessary variables *********
# -------------------------------------------------
 
# To get the right dataset in the software output folder where all combinations of dataset are found
ALIGN_FOLDER="left_54off_59on"
echo $ALIGN_FOLDER

# Create cluster lable with the appropriate number
DATASET="${COLOR_SPACE}_${MASK}_k${CLUSTER_NUMBER}"
echo $DATASET

# Find the appropriate mask file according to parameter given
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

#Repo for image output
mkdir $BASE_DIR/ressources/images/$ALIGN_FOLDER/3-registred/Modalities/RGB/discrete/
mkdir $BASE_DIR/ressources/images/$ALIGN_FOLDER/3-registred/Modalities/RGB/discrete/$COLOR_SPACE/
mkdir $BASE_DIR/ressources/images/$ALIGN_FOLDER/3-registred/Modalities/RGB/discrete/$COLOR_SPACE/$DATASET/

# Repo for pca outputs
mkdir $BASE_DIR/images/
mkdir $BASE_DIR/images/discrete/
mkdir $BASE_DIR/images/discrete/$COLOR_SPACE/
mkdir $BASE_DIR/images/discrete/$COLOR_SPACE/$DATASET/

# Repo for the corresponding dataset
mkdir $BASE_DIR/figures/
mkdir $BASE_DIR/figures/7_gxp/
mkdir $BASE_DIR/figures/7_gxp/discrete/
mkdir $BASE_DIR/figures/7_gxp/discrete/$COLOR_SPACE/
mkdir $BASE_DIR/figures/7_gxp/discrete/$COLOR_SPACE/$DATASET/



# ********* Run commands *********
# -------------------------------------------------

# Modify images according to color space, perform PCA files and create heatmaps per PCs
/user/doau0129/miniconda3/bin/python3 ./python/phenotype_discrete.py \
         $BASE_DIR/ressources/images/$ALIGN_FOLDER/3-registred/Modalities/RGB/$SUB_DATA/ \
         $BASE_DIR/ressources/images/$MASK_FILE \
         $BASE_DIR/ressources/images/$ALIGN_FOLDER/3-registred/Modalities/RGB/discrete/$COLOR_SPACE/$DATASET/ \
         $BASE_DIR/images/discrete/$COLOR_SPACE/$DATASET/ \
         $BASE_DIR/figures/7_gxp/discrete/$COLOR_SPACE/$DATASET/ \
         $CLUSTER_NUMBER \
         $COLOR_SPACE \
         $MASK \

# Load necessary environment
ml hpc-env/8.3
ml R/4.0.2-foss-2019b
ml FriBidi
ml HarfBuzz

# Create the image metadata table to use in further GWAS, plot the PCA
for folder in $BASE_DIR/images/discrete/$COLOR_SPACE/$DATASET/\*/
do

  echo "$folder"
  echo "$BASE_DIR/images/discrete/$COLOR_SPACE/$DATASET/$folder/"
  LABEL="${DATASET}_${folder}"
  echo $LABEL

  Rscript R/phenotype_discrete_pca.R \
         $BASE_DIR/images/discrete/$COLOR_SPACE/$DATASET/$folder/ \
         $BASE_DIR/metadata/ \
         ${LABEL}_PCs.csv \
         ${LABEL}_var.csv \
         $BASE_DIR/metadata/image_metadata.tsv \
         $LABEL \
         $BASE_DIR/figures/7_gxp/discrete/$COLOR_SPACE/$DATASET/$folder/ \

done