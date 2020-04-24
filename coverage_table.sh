#!/bin/bash

FILES=/Users/fco/Desktop/BREMEN_OP/chapter1_2/test/*

for f in $FILES
do
  echo $f
  name_file=${f##*/}
  echo $name_file
  sample1=${name_file%.*}
  sample2=${sample1%.*}
  sample=${sample2%.*}
  echo $sample
  cov=$(cat $f | awk '{if(!($2=="")) print $2}')
  echo $cov

  echo "$sample $cov" >> table.txt
done
