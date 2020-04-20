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
		pass
	return 0



#git0.0
# Split the samples and assign read groups - output as unmapped ubam temp_files
try_mkdir(base_directory+"outputs/")
try_mkdir(base_directory+"outputs/00_ubams/")

if(len(os.listdir(base_directory+"outputs/00_ubams/"))==0):

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

os.system("ls -1 "+base_directory+"outputs/00_ubams/ > "+base_directory+"/outputs/listoffiles/ubams.fofn")

if(len(os.listdir(base_directory+"outputs/01_adapters/adapters/"))==0 and len(os.listdir(base_directory+"outputs/01_adapters/metrics/"))==0):

    print("----------      mark adapters running      ----------")

    os.system("sbatch mark_adapters.sh")
else:
    print("**** mark adapters DONE ****")



#git0.2
#Map and merge
try_mkdir(base_directory+"outputs/02_merged_map/")
try_mkdir(base_directory+"outputs/02_merged_map/fatsq_adapters/")
try_mkdir(base_directory+"outputs/02_merged_map/sam_align/")
try_mkdir(base_directory+"outputs/02_merged_map/bam_merge_alignement/")

os.system("ls -1 "+base_directory+"outputs/01_adapters/adapters/ > "+base_directory+"outputs/listoffiles/adapters.bam.fofn")

if(len(os.listdir(base_directory+"outputs/02_merged_map/fatsq_adapters/"))==0):

	print("----------      map & merge sam to fastq running	----------")
	os.system("sbatch 02_samtofq")
else:
  print("**** map & merge sam to fastq DONE ****")


#os.system("ls -1 "+base_directory+"outputs/02_merged_map/fastq_adapters/ > "+base_directory+"outputs/listoffiles/samtofq.fofn")


#if(len(os.listdir(base_directory+"outputs/02_merged_map/sam_align/"))==0)

#	print("----------      map & merge sam align running ----------")
#	os.system("sbatch 02_samalign")

#else:
#  print("**** map & merge sam align DONE ****")


#os.system("ls -1 "+base_directory+"outputs/02_merged_map/sam_align/ > "+base_directory+"outputs/listoffiles/samalign.fofn")

#if(len(os.listdir(base_directory+"outputs/02_merged_map/bam_merge_alignement/"))==0):
	
#	print("----------      map & merge bam alignment running      ----------")
#	os.system("sbatch 02_merge.sh")
#else:
#    print("**** map & merge bam alignment DONE ****")



#git0.3

#try_mkdir(base_directory+"outputs/03_mark_duplicates/")
#try_mkdir(base_directory+"outputs/03_mark_duplicates/dedup/")
#try_mkdir(base_directory+"outputs/03_mark_duplicates/dedup_metrics/")

#if(len(os.listdir(base_directory+"outputs/03_mark_duplicates/dedup/"))==0 and len(os.listdir(base_directory+"outputs/03_mark_duplicates/dedup_metrics/")==0)):

#    print("----------      mark duplicates running      ----------")

#    os.system("sbatch mark_duplicates.sh")

#else:
#    print("**** mark duplicates DONE ****")
