#!/bin/bash
set -e
set -u
set -o pipefail
####################################################################################################################
####################################################################################################################
# a shell script to deal with fna results of DAMe-filter.py to get the files which can be used for UCHIME analyses
####################################################################################################################
####################################################################################################################

# Usage: bash loop_create_fas_forUCHIME.sh

PIPESTART=$(date)

HOMEFOLDER=$(pwd)

echo "Home folder is "${HOMEFOLDER}""

# set variables
INDEX=1

if [ ! -d for_usearch_temp ]
then
	mkdir for_usearch_temp # create a working folder
fi

cat FilteredReads.fna | tr '\t' '_' > for_usearch_temp/FilteredReads_temp.fna
cd for_usearch_temp/
mkdir repeat_reads

awk 'sub(/^>/, "")' FilteredReads_temp.fna > FilteredReads_temp_ID.list # keep the IDs of reads in this list file

read_info=FilteredReads_temp_ID.list
read_names=($(cut -f 1 "$read_info" | uniq))
for read in ${read_names[@]}
do
	cd "${HOMEFOLDER}"
	cd for_usearch_temp/
	mkdir ${read}_working
	cd ${read}_working/
	echo -n "${read}" > "${read}"_ID.txt
	INDEX=$((INDEX+1))
	seqtk subseq ../FilteredReads_temp.fna "${read}"_ID.txt > "${read}".fas
	sequence=$(sed '/^>/d' "${read}".fas)
	cat "${read}"_ID.txt | sed 's/__/_/g' | tr '_' '\t' > "${read}"_ID_tab.txt
	seq 1 $(cat "${read}"_ID_tab.txt | awk '{print $0 "\t" ($4+$5+$6)}' | cut -f 7) > numbers_list.txt
	sed "s/^/>${read}_$((INDEX-1))a/g" numbers_list.txt | sed 's/_.*_/_/g' > IDnumbers_list.txt
	sed "s/$/ ${sequence}/g" IDnumbers_list.txt | tr ' ' '\n' > "${read}"_repeat_reads.fna
	mv "${read}"_repeat_reads.fna ../repeat_reads/
	cd ..
	rm -rvf ${read}_working/

done

find repeat_reads -type f -name "*.fna" -exec cat '{}' ';' > ../FilteredReads_final_forUCHIME.fna
cd ../
rm -rvf for_usearch_temp

echo "Pipeline started at $PIPESTART"
NOW=$(date)
echo "Pipeline ended at   $NOW"
