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

    os.system("sbatch 00_split_samples.sh")
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

    os.system("sbatch 01_mark_adapters.sh")
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


os.system("ls -1 "+base_directory+"outputs/02_merged_map/fatsq_adapters/ > "+base_directory+"outputs/listoffiles/samtofq.fofn")


if(len(os.listdir(base_directory+"outputs/02_merged_map/sam_align/"))==0):

	print("----------      map & merge sam align running ----------")
	os.system("sbatch 02_samalign")
else:
  print("**** map & merge sam align DONE ****")


os.system("ls -1 "+base_directory+"outputs/02_merged_map/sam_align/ > "+base_directory+"outputs/listoffiles/samalign.fofn")

if(len(os.listdir(base_directory+"outputs/02_merged_map/bam_merge_alignement/"))==0):

	print("----------      map & merge bam alignment running      ----------")
	os.system("sbatch 02_merge.sh")
else:
    print("**** map & merge bam alignment DONE ****")



#git0.3
#Mark duplicates
try_mkdir(base_directory+"outputs/03_mark_duplicates/")
try_mkdir(base_directory+"outputs/03_mark_duplicates/sort_sam/")
try_mkdir(base_directory+"outputs/03_mark_duplicates/tags_intermediate/")
try_mkdir(base_directory+"outputs/03_mark_duplicates/mark_duplicates/")
try_mkdir(base_directory+"outputs/03_mark_duplicates/mark_duplicates/duplicates/")
try_mkdir(base_directory+"outputs/03_mark_duplicates/mark_duplicates/metrics/")


os.system("ls -1 "+base_directory+"outputs/02_merged_map/bam_merge_alignement/ > "+base_directory+"outputs/listoffiles/mapped_bams.fofn")
if(len(os.listdir(base_directory+"outputs/03_mark_duplicates/sort_sam/"))==0):
	print("----------       mark duplicates sort_sam running	----------")
	os.system("sbatch 03_sortsam.sh")
else:
	print("**** mark duplicates sort_sam DONE ****")


os.system("ls -1 "+base_directory+"outputs/03_mark_duplicates/sort_sam/ > "+base_directory+"outputs/listoffiles/sort_sam.fofn")
if(len(os.listdir(base_directory+"outputs/03_mark_duplicates/tags_intermediate/"))==0):
	print("----------       mark duplicates tags_intermediate running	----------")
	os.system("sbatch 03_tagsinter.sh")
else:
	print("**** mark duplicates tags_intermediate DONE ****")


os.system("ls -1 "+base_directory+"outputs/03_mark_duplicates/tags_intermediate/*.bam |xargs -n1 basename > "+base_directory+"outputs/listoffiles/tags_intermediate.fofn")
if(len(os.listdir(base_directory+"outputs/03_mark_duplicates/mark_duplicates/duplicates"))==0 and len(os.listdir(base_directory+"outputs/03_mark_duplicates/mark_duplicates/metrics"))==0):
	print("----------       mark duplicates mark_duplicates running	----------")
	os.system("sbatch 03_markdups.sh")
else:
	print("**** mark duplicates mark_duplicates DONE ****")



# git0.4
# Index the mapped and mark duplicates bam files
#
#try_mkdir(base_directory+"outputs/04_indexed/")
#
os.system("ls -1 "+base_directory+"outputs/03_mark_duplicates/mark_duplicates/duplicates/ > "+base_directory+"outputs/listoffiles/duplicates.fofn")
if(len(os.listdir(base_directory+"outputs/03_mark_duplicates/mark_duplicates/duplicates/"))<=117):
	print("----------       index running	----------")
	os.system("sbatch 04_index.sh")
else:
	print("**** index DONE ****")


# git0.5
# Coverage

try_mkdir(base_directory+"outputs/a_coverage/")

os.system("ls -1 "+base_directory+"outputs/03_mark_duplicates/mark_duplicates/duplicates/*.bam |xargs -n1 basename > "+base_directory+"outputs/listoffiles/duplicates_coverage.fofn")
if(len(os.listdir(base_directory+"outputs/03_mark_duplicates/mark_duplicates/duplicates/"))==234 and len(os.listdir(base_directory+"outputs/a_coverage/"))==0):
	print("----------     coverage running	    ----------")
	os.system("sbatch a_coverage.sh")
elif(len(os.listdir(base_directory+"outputs/03_mark_duplicates/mark_duplicates/duplicates/"))==117 and len(os.listdir(base_directory+"outputs/a_coverage/"))==0):
	print("*** coverage NOT RUNNING : step 0.4 indexing NOT COMPLETE ***")
elif(len(os.listdir(base_directory+"outputs/03_mark_duplicates/mark_duplicates/duplicates/"))==0 and len(os.listdir(base_directory+"outputs/a_coverage/"))==0):
	print("*** coverage NOT RUNNING : step 0.3 mark duplicates NOT COMPLETE ***")
elif(len(os.listdir(base_directory+"outputs/a_coverage/"))==117):
	print("**** coverage DONE ****")
else:
	print("**** UNKNOWN ****")


try_mkdir(base_directory+"outputs/b_removed_from_analysis/")
os.system("mv "+base_directory+"outputs/03_mark_duplicates/mark_duplicates/duplicates/AG9RX_47pueboc.1.dedup.bam "+base_directory+"outputs/b_removed_from_analysis/")
os.system("mv "+base_directory+"outputs/03_mark_duplicates/mark_duplicates/duplicates/AG9RX_47pueboc.1.dedup.bai "+base_directory+"outputs/b_removed_from_analysis/")
os.system("mv "+base_directory+"outputs/03_mark_duplicates/mark_duplicates/duplicates/PL17_79abepue.3.dedup.bam "+base_directory+"outputs/b_removed_from_analysis/")
os.system("mv "+base_directory+"outputs/03_mark_duplicates/mark_duplicates/duplicates/PL17_79abepue.3.dedup.bai "+base_directory+"outputs/b_removed_from_analysis/")
os.system("mv "+base_directory+"outputs/03_mark_duplicates/mark_duplicates/duplicates/PL17_101maybel.1.dedup.bam "+base_directory+"outputs/b_removed_from_analysis/")
os.system("mv "+base_directory+"outputs/03_mark_duplicates/mark_duplicates/duplicates/PL17_101maybel.1.dedup.bai "+base_directory+"outputs/b_removed_from_analysis/")
os.system("mv "+base_directory+"outputs/03_mark_duplicates/mark_duplicates/duplicates/PL17_98indbel.3.dedup.bam "+base_directory+"outputs/b_removed_from_analysis/")
os.system("mv "+base_directory+"outputs/03_mark_duplicates/mark_duplicates/duplicates/PL17_98indbel.3.dedup.bai "+base_directory+"outputs/b_removed_from_analysis/")

try_mkdir(base_directory+"outputs/05_genlikely/")
os.system("ls -1 "+base_directory+"outputs/03_mark_duplicates/mark_duplicates/duplicates/*.bam |xargs -n1 basename > "+base_directory+"outputs/listoffiles/duplicates_hapcaller.fofn")
if(len(os.listdir(base_directory+"outputs/05_genlikely/"))==0):
	print("----------     HaplotypeCaller running	    ----------")
	os.system("sbatch 05_genlikely.sh")
else:
	print("*** HaplotypeCaller DONE ***")



#git0.6
if(len(os.listdir(base_directory+"outputs/05_genlikely/"))==226):
	print("----------     CombineGVCF running	    ----------")
	os.system("sbatch 06_combine.sh")
elif(len(os.listdir(base_directory+"outputs/05_genlikely/"))==228):
	print("*** CombineGVCF DONE ***")
else:
	print("either step 05 not done or other problem")
