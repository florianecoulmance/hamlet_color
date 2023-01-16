# Code repository for: Phenotypic and genomic dissection of color pattern variation in a reef fish radiation

This repository contains the original bioinformatic analysis behind the paper *Phenotypic and genomic dissection of color pattern variation in a reef fish radiation* by Coulmance, Akkaynak, Le Poul, Hoppner, McMillan and Puebla.<br> 
It covers all steps from genotyping based on raw sequencing data, over population genetic analysis to the final plotting of the figures used within the publication.<br> 
It covers all image analysis steps from aligned photographs.<br> 

There are four more accompanying repositories for this publication:

    The ENA sequencing repository: contains the raw sequencing data (Accession Nr: ___)
    The photographs repository which contains all raw photographs used in the paper and for the image analysis: ___
    Code for image correction in MatLab: ___
    Software for image alignment: ___


## Setup

In this folder, you will need all necessary files and subfolders to be able to run the different pipelines : 
The folder tree is:<br> 

```
```

All the output files and figures will be created by the various pipelines.<br> 
All figures presented in the paper will be found in the figures/ folder.<br> 

All necessary code for each of these pipelines is found in folders R/ python/ and sh/<br> 
The pipelines are all the files finishing by .sh in the root folder.<br> 


## Run pipelines<br> 

<PATH> is the path to this folder<br> 
<JOB_ID> corresponds to the ID of the job within the pipeline that you want to run if jobfile1 then you should specify -j jid1<br> 
If you want to run the pipeline from the start, do not specify the option -j<br> 
See example for Pipeline 1<br>  

### Pipeline 1: Genotyping<br> 
genotyping.sh<br> 
sbatch genotyping.sh -i <PATH> -j <JOB_ID><br> 
example : sbatch genotyping.sh -i /user/doau0129/work/hamlet_color/ -j jid6 (will start from jobfile6 in genotyping.sh)<br> 

### Pipeline 2: Fst statistics<br> 
fst.sh<br> 
sbatch fst.sh -i <PATH> -j <JOB_ID><br> 

### Pipeline 3: Hybrid and backcross identification<br> 
hybrids.sh<br> 
sbatch hybrids.sh -i <PATH> -j <JOB_ID><br> 

### Pipeline 4: Image analysis<br> 
continuous_image_pca.sh<br> 
sbatch continuous_image_pca.sh -i <PATH> -j <COLOR_SPACE> -k <MASK> -l <DATA> -m <SUB_DATA><br> 
For this script, these parameters should be specified:<br> 
-j LAB<br> 
-k fullm<br> 
-l LAB_54off_59on<br> 
-m all<br> 

### Pipeline 5: GWAS<br> 
gxp.sh<br> 
sbatch gxp.sh -I <PATH> -j <JOB_ID> -k <DATASET><br> 
Be sure to give parameter -k the prefix corresponding to a metadata file.<br> 
If the metadata file is named LAB_fullm_54off_59on_PCs.csv, the prefix will be LAB_fullm_54off_59on<br> 

### Pipeline 6: Alleles analysis<br> 
snp_alleles.sh<br> 
sbatch gxp.sh -I <PATH> -j <JOB_ID> -k <DATASET><br> 
Be sure to give parameter -k the prefix corresponding to a metadata file.<br> 
If the metadata file is named LAB_fullm_54off_59on_PCs.csv, the prefix will be LAB_fullm_54off_59on<br> 


Pipelines should be run in this order:<br> 
Pipeline 2,3 after Pipeline 1<br> 
Pipeline 5 after Pipeline 1 and 4<br> 
Pipeline 6 after Pipeline 5<br> 






