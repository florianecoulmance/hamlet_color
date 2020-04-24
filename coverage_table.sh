#!/bin/bash

FILES=/user/doau0129/work/chapter1_2/outputs/a_coverage/*

for f in $FILES
do
  echo $f
  name_file=${f##*/}
  echo $name_file
  sample1=${name_file%.*}
  sample2=${sample1%.*}
  sample=${sample2%.*}
  echo $sample
  cov=$(cat $f | awk 'FNR == 8 {print $2}')
  echo $cov

  echo "$sample $cov" >> coverage_table
done
