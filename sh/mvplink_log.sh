#!/usr/bin/bash
# by: Kosmas Hench: 2019-01-22
# usage gxp_slider <in.GxP.txt.gz> <Windowsize (bp)> <Increment (bp)>
# -------------------------------------------------------------------------------------------------------------------
# !! beware of awk array size maximum - if the sequence lenght is to big or if the window size & step are to small !!
# !!   (hence, if too many windows are created, the array index maximum is early entries are being overwritten     !!
# -------------------------------------------------------------------------------------------------------------------
# $1 = infile; $2 = windowsize; $3 = increment
j=$(echo $1 | sed 's/.txt.gz//g')
RED='\033[0;31m';
NC='\033[0m';

# print the header (the first line of input)
# and then run the specified command on the body (the rest of the input)
# use it in a pipeline, e.g. ps | body grep somepattern
body() {
    IFS= read -r header
    printf '%s\n' "$header"
    "$@"
}


zcat $1 | head -n 1  | cut -d ' ' -f 1,2,3,4,5 > ${j}.logarithm.txt.gz

for k in {01..24};do
  echo -e  "--- ${RED}$j${NC} - $k ---"

	zcat $1 | \
		        awk -v k=LG$k -v OFS='\t' '$1==k {print}' | \
		        cut -d ' ' -f 1,2,3,4,5 | \
		        awk -v OFS=' ' '{$4=-log($4)/log(10);print $0}' | \
            body sort -k2 -n | \
            sed -r '/^\s*$/d' >> ${j}.logarithm.txt.gz

done
