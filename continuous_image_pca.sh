#!/bin/bash
# by: Floriane Coulmance: 11/05/2021
# usage:
# sbatch continuous_image_pca.sh -i <PATH> -j <COLOR_SPACE> -k <MASK> -l <DATA> -m <SUB_DATA>
# sbatch continuous_image_pca.sh -i hamlet_color/ -j LAB -k fullm -l 54off_59on -m all 
# ------------------------------------------------------------------------------
# PATH corresponds to the path to the base directory, all outputs and necessary
# folder will be created by the script 
# COLOR_SPACE : LAB
# MASK : fullm (for mask including all fins)
# DATA : 54off_59on
# SUB_DATA : all
# ------------------------------------------------------------------------------



# Slurm parameters
#SBATCH --job-name=continuous_images.tmp
#SBATCH --partition=carl.p
##SBATCH --output=$BASE_DIR/logs/continuous_images_%A_%a.out
##SBATCH --error=$BASE_DIR/logs/continuous_images_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=22G
#SBATCH --time=02:00:00



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
esac
done



# ********* Create necessary variables *********
# -------------------------------------------------

# To get the right dataset in the software output folder where all combinations of dataset are found
ALIGN_FOLDER="left_54off_59on"
echo $ALIGN_FOLDER
 
# Create appropriate name to give to the output folder name
DATASET="${COLOR_SPACE}_${MASK}_${DATA}"
echo $DATASET

# FILES_NAME="${COLOR_SPACE}_${MASK}_${DATA}"
# echo $FILES_NAME

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

# Repo for pca outputs
mkdir $BASE_DIR/images/
mkdir $BASE_DIR/images/continuous/
mkdir $BASE_DIR/images/continuous/$COLOR_SPACE/
mkdir $BASE_DIR/images/continuous/$COLOR_SPACE/$DATASET/

# Repo for the figures
mkdir $BASE_DIR/figures/
mkdir $BASE_DIR/figures/7_gxp/
mkdir $BASE_DIR/figures/7_gxp/continuous/
mkdir $BASE_DIR/figures/7_gxp/continuous/$COLOR_SPACE/
mkdir $BASE_DIR/figures/7_gxp/continuous/$COLOR_SPACE/$DATASET/



# ********* Run commands *********
# -------------------------------------------------

# Modify images according to color space, perform PCA files and create heatmaps per PCs
/user/doau0129/miniconda3/bin/python3 ./python/phenotype_continuous.py \
         $BASE_DIR/ressources/images/$ALIGN_FOLDER/3-registred/Modalities/RGB/$SUB_DATA/ \
         $BASE_DIR/ressources/images/$MASK_FILE \
         $COLOR_SPACE \
         $BASE_DIR/images/continuous/$COLOR_SPACE/$DATASET/ \
         $BASE_DIR/figures/7_gxp/continuous/$COLOR_SPACE/$DATASET/ \
         $MASK \
         $DATASET \

# Load necessary environment
ml hpc-env/8.3
ml R/4.0.2-foss-2019b
ml FriBidi
ml HarfBuzz

# Create the image metadata table to use in further GWAS, plot the PCA
# Rscript R/phenotype_continuous_pca.R \
#          $BASE_DIR/images/continuous/$COLOR_SPACE/$DATASET/ \
#          $COLOR_SPACE \
#          $BASE_DIR/metadata/ \
#          ${DATASET}_PCs.csv \
#          ${DATASET}_var.csv \
#          $BASE_DIR/metadata/image_metadata.tsv \
#          $MASK \
#          $DATASET \
#          $BASE_DIR/figures/7_gxp/continuous/$COLOR_SPACE/$DATASET/ \
