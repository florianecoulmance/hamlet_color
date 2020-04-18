#!/usr/bin/python

import os
import sys
import errno


# set the base directory path of the whole project
base_directory="/user/doau0129/work/chapter1_2/"


# Function to catch error if folders already exists
def try_mkdir(dirname):
    	try:
		os.mkdir(dirname)
	except OSError as exc:
		if exc.errno != errno.EEXIST:
			raise
            return 1
		pass
	return 0



#git0.0
# Split the samples and assign read groups - output as unmapped ubam temp_files
try_mkdir(base_directory+"outputs/")
try_mkdir(base_directory+"outputs/00_ubams/")

if(len(os.lisdir(base_directory+"outputs/00_ubams/"))==0):

    print("----------      split samples running      ----------")

    os.system("sbatch split_samples.sh")
else:
    print("**** split samples DONE ****")




#git0.1
# Mark adapters on the unmappe ubam files - outputs : adapter.bam & adapter.metrics
try_mkdir(base_directory+"outputs/01_adapters/")
try_mkdir(base_directory+"/outputs/01_adapters/adapters/")
try_mkdir(base_directory+"/outputs/01_adapters/metrics/")
try_mkdir(base_directory+"/outputs/listoffiles/")
os.system("ls -1 "+base_directory+"/outputs/01_adapters/adapters/ > "+base_directory+"/outputs/listoffiles/ubams.fofn")

if(len(os.lisdir(base_directory+"outputs/01_adapters/adapters/"))==0 and len(os.lisdir(base_directory+"outputs/01_adapters/metrics/")==0)):

    print("----------      mark adapters running      ----------")

    os.system("sbatch split_samples.sh")
else:
    print("**** mark adapters DONE ****")

#git0.2

try_mkdir(base_directory+"outputs/02_merged_map/")
try_mkdir(base_directory+"outputs/02_merged_map/fatsq_adapters/")
try_mkdir(base_directory+"outputs/02_merged_map/sam_align/")
try_mkdir(base_directory+"outputs/02_merged_map/bam_merge_alignement/")

if(len(os.lisdir(base_directory+"outputs/02_merged_map/fatsq_adapters/"))==0 and len(os.lisdir(base_directory+"outputs/02_merged_map/sam_align/")==0) and len(os.lisdir(base_directory+"outputs/02_merged_map/bam_merge_alignement/")==0)):

    print("----------      map & merge running      ----------")

    os.system("sbatch merge_ubams.sh")
else:
    print("**** map & merge DONE ****")

#git0.3

try_mkdir(base_directory+"outputs/03_mark_duplicates/")
try_mkdir(base_directory+"outputs/03_mark_duplicates/dedup/")
try_mkdir(base_directory+"outputs/03_mark_duplicates/dedup_metrics/")

if(len(os.lisdir(base_directory+"outputs/03_mark_duplicates/dedup/"))==0 and len(os.lisdir(base_directory+"outputs/03_mark_duplicates/dedup_metrics/")==0)):

    print("----------      mark duplicates running      ----------")

    os.system("sbatch mark_duplicates.sh")

else:
    print("**** mark duplicates DONE ****")
