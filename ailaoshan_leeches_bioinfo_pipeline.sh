#!/bin/bash
set -e
set -u
set -o pipefail
##########################################################################################################
##########################################################################################################
# shell script for ailaoshan leeches' metabarcoding data
##########################################################################################################
##########################################################################################################

# Usage: bash ECEC_ailaishan_leeches_bioinfo_pipeline.sh

# PIPESTART=$(date)

# run script from ~/PATH_TO/ailaishan_leeches_method2019/scripts
# fq files in ailaishan_leeches_method2019/data/raw_data/
# info files about primers, tags and PCR_sets in ailaishan_leeches_method2019/info/
# analysis outputs in ailaishan_leeches_method2019/analysis/

# set variables
geneNUMBER=1
HOMEFOLDER="/Volumes/JYQ/Ailaoshan/ailaoshan_leeches_method2019/"  # do not have a ~ in the path
DAME="/Users/apple/Downloads/DAMe-master-2/bin/"
SWARM="/Users/apple/Downloads/swarm/bin/swarm"

## read in gene folder list and make a bash array
#find data/raw_data/* -maxdepth 0 -type d | sed 's/^.*\///g' > genelist.txt  # find all gene folders in raw_data
#gene_info=genelist.txt # put genelist.txt into variable
#gene_names=($(cut -f 1 "${gene_info}" | uniq)) # convert variable to array this way
#echo "There are" "${#gene_names[@]}" "genes that will be processed:  " "${gene_names[@]}" # echo number of elements in the array
#
#for gene in ${gene_names[@]}  # ${gene_names[@]} is the full bash array
#do
#	echo "Now on gene" ${geneNUMBER} of ${#gene_names[@]}". Moved back to starting directory:"
#	geneNUMBER=$((geneNUMBER+1))
#  INDEX=1
#
#  if [ ! -d ${gene}_DAMe_SORT_outputs ] # if directory ${gene}_DAMe_SORT_outputs does not exist
#  then
#  	mkdir ${gene}_DAMe_SORT_outputs
#  fi
#
#  # read in library folder list and make a bash array
#  find data/raw_data/${gene}/* -maxdepth 0 -type d | sed 's/^.*\///g' > ${gene}_librarylist.txt # find all libraries in this gene folder
#  library_info=${gene}_librarylist.txt # put librarylist.txt into variable
#  library_names=($(cut -f 1 "${library_info}" | uniq)) # convert variable to array this way
#  echo "there are" "${#library_names[@]}" "libraries for ${gene} that will be processed." # echo number of elements in the array
#
#  for library in ${library_names[@]} # ${library_names[@]} is the full bash array
#  do
#    echo "Now on library" ${INDEX} of ${#library_names[@]}
#    INDEX=$((INDEX+1))
#
#    mkdir ${library}_working
#    cp data/raw_data/${gene}/${library}/*.fq.gz ${library}_working
#    cd ${library}_working
#
#    # Use AdapterRemoval to trim Illumina sequencing adapters
#    # brew install adapterremoval
#    # or download from https://github.com/MikkelSchubert/adapterremoval
#    AdapterRemoval --file1 ${library}_R1.fq.gz --file2 ${library}_R2.fq.gz --output1 ${library}_R1_ADtrimed.fq.gz --output2 ${library}_R2_ADtrimed.fq.gz --gzip --trimns
#    rm ${library}_R1.fq.gz; rm -f ${library}_R2.fq.gz
#
#    # Sickle error correction strategies and identified quality trimming
#    # brew install sickle
#    # or download from https://github.com/najoshi/sickle
#    sickle pe -f ${library}_R1_ADtrimed.fq.gz -r ${library}_R2_ADtrimed.fq.gz -t sanger -o ${library}_R1_ADStrimed.fastq.gz -p ${library}_R2_ADStrimed.fastq.gz -s ${library}_trimmed_singles_file.fastq.gz -g >> sickle.out
#    rm -f *_ADtrimed.fq.gz
#
#    # run bfc to do error correction
#    # brew install bfc
#    # or download from https://github.com/lh3/bfc
#    bfc -s 3g -k 25 -t4 ${library}_R1_ADStrimed.fastq.gz | gzip -9 > ${library}_R1_ADStrimed_bfc.fastq.gz
#    bfc -s 3g -k 25 -t4 ${library}_R2_ADStrimed.fastq.gz | gzip -9 > ${library}_R2_ADStrimed_bfc.fastq.gz
#    rm -f *_ADStrimed.fastq.gz
#
#    # use pandaseq to merge pair-end reads
#    # brew install homebrew/science/pandaseq
#    # or download from https://github.com/neufeld/pandaseq/releases
#      # because bfc_outputs can't be used directly by pandaseq, pairfq has to be used to pair reads before running pandaseq
#      # download pairfq from https://github.com/sestaton/Pairfq
#        # because big data need a lot of memory in pairing reads, it's better to divide data into small sets first in case of no enough memory
#        # use fastq-splitter.pl to divide data, download it from http://kirill-kryukov.com/study/tools/fastq-splitter/
#        gunzip -d ${library}_R*_ADStrimed_bfc.fastq.gz
#        cat ${library}_R1_ADStrimed_bfc.fastq | tr "\t" " " > temp.fastq; rm -f ${library}_R1_ADStrimed_bfc.fastq; sed 's/ .*/ 1:N:0:0/g' temp.fastq > ${library}_R1_ADStrimed_bfc.fastq # modify the IDs in bfc_outputs
#        cat ${library}_R2_ADStrimed_bfc.fastq | tr "\t" " " > temp.fastq; rm -f ${library}_R2_ADStrimed_bfc.fastq; sed 's/ .*/ 2:N:0:0/g' temp.fastq > ${library}_R2_ADStrimed_bfc.fastq # modify the IDs in bfc_outputs
#        rm -f temp.fastq
#        fastq-splitter.pl --n-parts 10 ${library}_R1_ADStrimed_bfc.fastq
#        rm -f ${library}_R1_ADStrimed_bfc.fastq
#        fastq-splitter.pl --n-parts 10 ${library}_R2_ADStrimed_bfc.fastq
#        rm -f ${library}_R2_ADStrimed_bfc.fastq
#      partnumber=1
#      find * -maxdepth 0 -name "${library}_R1_ADStrimed_bfc.part-*.fastq" | sed 's/^.*part-//g' | sed 's/.fastq//g' > partnumberlist.txt # find all part bumbers
#      part_info=partnumberlist.txt # put partnumberlist.txt into variable
#      part_names=($(cut -f 1 "${part_info}" | uniq)) # convert variable to array this way
#      for part in "${part_names[@]}"
#      do
#        partnumber=$((partnumber+1))
#        pairfq makepairs -f ${library}_R1_ADStrimed_bfc.part-${part}.fastq -r ${library}_R2_ADStrimed_bfc.part-${part}.fastq -fp R1_paired.part-${part}.fastq -rp R2_paired.part-${part}.fastq -fs r1_sing.part-${part}.fastq -rs r2_sing.part-${part}.fastq --stats
#        rm -f ${library}_R*_ADStrimed_bfc.part-${part}.fastq
#      done
#      cat r1_sing.part-*.fastq > r1_sing.all.fastq # combine all the unpaired reads into one file
#      cat r2_sing.part-*.fastq > r2_sing.all.fastq # combine all the unpaired reads into one file
#      # run pairfq with all the unpaired reads again to make sure that none of the paired reads won't be missed
#      pairfq makepairs -f r1_sing.all.fastq -r r2_sing.all.fastq -fp R1_paired.part-s.fastq -rp R2_paired.part-s.fastq -fs r1_sing.part-s.fastq -rs r2_sing.part-s.fastq --stats
#      rm -f r*.fastq # delete all the unpaired reads data files
#      # combine all the paired reads
#      cat R1_paired.part-*.fastq > R1_paired.fastq; rm -f R1_paired.part-*.fastq
#      cat R2_paired.part-*.fastq > R2_paired.fastq; rm -f R2_paired.part-*.fastq
#    # run pandaseq
#    pandaseq -f R1_paired.fastq -r R2_paired.fastq -A simple_bayesian -B -F -d bfsrk -N -g ${library}_pandaseq_log.txt -w ${library}_ADStrimed_bfc_merged.fastq
#    rm -f R*_paired.fastq
#    gzip -9 ${library}_ADStrimed_bfc_merged.fastq
#
#    #############################################################################################
#    #### DAMe - using PCR replicates and read numbers to filter out bad reads
#    # download DAMe from https://github.com/shyamsg/DAMe
#
#    # DAMe needs each library data in a different folder named "pool*"
#    # prepare a file "${gene}_pool_info.txt" and get pool NUMBER of each library from this file
#    NUMBER=$(grep -E "${library}" "${HOMEFOLDER}"/info/${gene}_pool_info.txt | cut -f 2)
#    mkdir pool${NUMBER}
#    mv ${library}_ADStrimed_bfc_merged.fastq.gz pool${NUMBER}/
#
#    # Sort through each fastq file and determine how many of each tag pair is in each fastq file
#    cd pool${NUMBER}/
#    ${DAME}sort.py -fq ${library}_ADStrimed_bfc_merged.fastq.gz -p "${HOMEFOLDER}"/info/Primers_leeches_${gene}.txt -t "${HOMEFOLDER}"/info/8bp_Tags_leeches.txt
#    # Sort SummaryCounts.txt to SummaryCounts_sorted_${library}_Pool${NUMBER}.txt
#    head -1 SummaryCounts.txt > SummaryCounts_sorted_${library}_Pool${NUMBER}.txt
#    tail -n +2 SummaryCounts.txt | sed "s/Tag//g" | sort -k1,1n -k2,2n | awk 'BEGIN{OFS="\t";}{$1="Tag"$1;$2="Tag"$2; print $0;}' >> SummaryCounts_sorted_${library}_Pool${NUMBER}.txt
#    # Make tag combinations overview
#    ${DAME}splitSummaryByPSInfo.py -p "${HOMEFOLDER}"/info/PCRsetsInfo_${gene}.txt -l ${NUMBER} -s SummaryCounts.txt -o splitSummaryByPSInfo_${library}_Pool${NUMBER}.txt
#    cd ..
#    mv pool${NUMBER} "${HOMEFOLDER}"/${gene}_DAMe_SORT_outputs/
#    cd "${HOMEFOLDER}"
#    rm -rf ${library}_working
#
#  done
#
#  cd "${HOMEFOLDER}"
#  mv ${gene}_DAMe_SORT_outputs analysis/
#
#done
#
## DAMe filter
#  # First fitler.py at -y 1 -t 1 to keep all sequences in the control samples, confirm MINPCR and MINREADS values through the controls' outputs
#  #  -psInfo PSINFO        Text file with the information on the tag combination
#  #                        in each PCR reaction for every sample [Format:
#  #                        sampleName TagNameForwardFromPCR1
#  #                        TagNameReverseFromPCR1 Pool# sampleName
#  #                        TagNameForwardFromPCR2 TagNameReverseFromPCR2 Pool#
#  #                        ...]
#  #  -x X                  Number of PCR rxns performed per sample
#  #  -y Y                  Number of PCR rxns in which the sequence has to be
#  #                        present
#  #  -p P                  The number of pools in which the samples were divided
#  #                        for sequencing (in case of tag combinations repetition
#  #                        due to processing so many samples) [default 1] NOTE:
#  #                        If using pools, each fastq must be in a folder called
#  #                        pool#, in which the sort.py was run for each pool
#  #                        inside the corresponding folder, and this program
#  #                        chimeraCheck.py is run in the parent directory of the
#  #                        pools directories
#  #  -t T                  Number of times a unique sequence has to be present
#  #  -l L                  Minimum sequence length
#
#  # 12S
#  cd "${HOMEFOLDER}"/analysis/12S_DAMe_SORT_outputs
#  mkdir Filter_min1PCRs_min1copies_12S_Controls
#  ${DAME}filter.py -psInfo "${HOMEFOLDER}"/info/PCRsetsInfo_12S_Controls.txt -x 3 -y 1 -p 14 -t 1 -l 81 -o Filter_min1PCRs_min1copies_12S_Controls
#
#  # 16S
#  cd "${HOMEFOLDER}"/analysis/16S_DAMe_SORT_outputs
#  mkdir Filter_min1PCRs_min1copies_16S_Controls
#  ${DAME}filter.py -psInfo "${HOMEFOLDER}"/info/PCRsetsInfo_16S_Controls.txt -x 3 -y 1 -p 13 -t 1 -l 82 -o Filter_min1PCRs_min1copies_16S_Controls
#
#  ####################################################################################################
#  ## After consideration of the controls' results, choose thresholds for filtering
#  # The min PCR and copy number can be informed by looking at negative controls.
#  # After observing the negative controls, set -y and -t higher than the observed values in neg controls
#  ####################################################################################################
#  # Filter - Filtering reads with minPCR and minReads/PCR thresholds
#
#  # 12S
#  cd "${HOMEFOLDER}"/analysis/12S_DAMe_SORT_outputs
#  mkdir Filter_min2PCRs_min20copies_12S
#  ${DAME}filter.py -psInfo "${HOMEFOLDER}"/info/PCRsetsInfo_12S.txt -x 3 -y 2 -p 14 -t 20 -l 81 -o Filter_min2PCRs_min20copies_12S
#
#  # 16S
#  cd "${HOMEFOLDER}"/analysis/16S_DAMe_SORT_outputs
#  mkdir Filter_min2PCRs_min9copies_16S
#  ${DAME}filter.py -psInfo "${HOMEFOLDER}"/info/PCRsetsInfo_16S.txt -x 3 -y 2 -p 13 -t 9 -l 82 -o Filter_min2PCRs_min9copies_16S
#
#####################################################################################################
#################################### remove chimera reads ###########################################
#####################################################################################################
## Prepare the file for removing chimeras
#	# need seqtk in this step
#	# brew install seqtk
#
#	# 12S
#	cd "${HOMEFOLDER}"/analysis/12S_DAMe_SORT_outputs/Filter_min2PCRs_min20copies_12S/
#	bash "${HOMEFOLDER}"/loop_create_fas_forUCHIME.sh
#
#	# 16S
#	cd "${HOMEFOLDER}"/analysis/16S_DAMe_SORT_outputs/Filter_min2PCRs_min9copies_16S/
#	bash "${HOMEFOLDER}"/loop_create_fas_forUCHIME.sh
#
## UCHIME
#	# brew install vsearch
#
#	# 12S
#	cd "${HOMEFOLDER}"/analysis/12S_DAMe_SORT_outputs/Filter_min2PCRs_min20copies_12S/
#	mkdir UCHIME
#	# run -derep_fullseq to dereplicate the same reads, considering that the input file of UCHIME must contain estimated amplicons with abundances specified by size annotations
#	vsearch --derep_fulllength FilteredReads_final_forUCHIME.fna --output UCHIME/uniques.fasta --sizeout
#	cd UCHIME/
#	# run uchime_denovo to find chimeras
#	vsearch --uchime_denovo uniques.fasta --chimeras denovo_chimeras.fasta
#	# remove chimeras
#	vsearch --search_exact ../FilteredReads_final_forUCHIME.fna -db denovo_chimeras.fasta -notmatched seqs_nochimeras.fasta -strand plus
#	gzip -9 ../FilteredReads_final_forUCHIME.fna; mv seqs_nochimeras.fasta ../; cd ../; rm -rf UCHIME
#	# filter reads based on the length
#	vsearch --fastx_filter seqs_nochimeras.fasta -fastaout 12S_seqs_nochimeras_81to117bp.fasta -fastq_minlen 81 -fastq_maxlen 117
#	rm -f seqs_nochimeras.fasta
#	gzip -9 *.fasta; gzip -9 *.fna
#
#	# 16S
#	cd "${HOMEFOLDER}"/analysis/16S_DAMe_SORT_outputs/Filter_min2PCRs_min9copies_16S/
#	mkdir UCHIME
#	# run -derep_fullseq to dereplicate the same reads, considering that the input file of UCHIME must contain estimated amplicons with abundances specified by size annotations
#	vsearch --derep_fulllength FilteredReads_final_forUCHIME.fna --output UCHIME/uniques.fasta --sizeout
#	cd UCHIME/
#	# run uchime_denovo to find chimeras
#	vsearch --uchime_denovo uniques.fasta --chimeras denovo_chimeras.fasta
#	# remove chimeras
#	vsearch --search_exact ../FilteredReads_final_forUCHIME.fna -db denovo_chimeras.fasta -notmatched seqs_nochimeras.fasta -strand plus
#	gzip -9 ../FilteredReads_final_forUCHIME.fna; mv seqs_nochimeras.fasta ../; cd ../; rm -rf UCHIME
#	# filter reads based on the length
#	vsearch --fastx_filter seqs_nochimeras.fasta -fastaout 16S_seqs_nochimeras_82to150bp.fasta -fastq_minlen 82 -fastq_maxlen 150
#	rm -f seqs_nochimeras.fasta
#	gzip -9 *.fasta; gzip -9 *.fna

####################################################################################################
############################################ pick_otu ##############################################
####################################################################################################
# use SWARM + LULU to pick OTUs

# do clustering by SWARM
	# download SWARM from https://github.com/torognes/swarm.git
	# cd swarm/src/; make
	# need USEARCH 5.2.236 and MACQIIME 1.9.1 in this step
	# download USEARCH 5.2.236 from http://www.drive5.com/usearch/; download MACQIIME 1.9.1 from http://www.wernerlab.org/software/macqiime
	# installing USEARCH: sudo mv /Users/apple/Downloads/usearch5.2.236_i86osx32 /usr/local/bin/usearch5.2.236_i86osx32; cd /usr/local/bin/; chmod 755 usearch5.2.236_i86osx32; sudo ln -s usearch5.2.236_i86osx32 usearch5
	# for details of installing MACQIIME 1.9.1, please see http://www.wernerlab.org/software/macqiime/macqiime-installation

#	# 12S
#	cd "${HOMEFOLDER}"/analysis/12S_DAMe_SORT_outputs/Filter_min2PCRs_min20copies_12S/
#	mkdir swarm_lulu; cp 12S_seqs_nochimeras_81to117bp.fasta.gz swarm_lulu/; cd swarm_lulu; gunzip -d 12S_seqs_nochimeras_81to117bp.fasta.gz
#	# SWARM needs the fasta input file which had better be dereplicated First
#	# to obtain the otu-map file of dereplicating step, Pick_cluster_num_V2.pl is used and is written for the output uc file of USEARCH 5.2.236, if you use other versions of USEARCH, you have to rewrite Pick_cluster_num_V2.pl
#	usearch5 -derep_fullseq -cluster 12S_seqs_nochimeras_81to117bp.fasta -seedsout uniques.fasta -sizeout -uc uniques.uc
#	mkdir for_map; mv uniques.uc for_map/; perl "${HOMEFOLDER}"/Pick_cluster_num_V2.pl -id for_map -o usearch_map_temp.txt; sed 's/uniques_/uniques/g' usearch_map_temp.txt > usearch_map.txt
#	sed 's/Cluster/uniques/g' uniques.fasta > temp.txt; sed 's/;size=/_/g' temp.txt > uniques_for_swarm.fasta
#	# run SWARM
#	${SWARM} -t 4 -f -o swarm_OTU_map_temp.txt uniques_for_swarm.fasta > /dev/null
#	# modify swarm_OTU_map.txt to make it available for MACQIIME
#	seq 1 $(grep -c "^" swarm_OTU_map_temp.txt) | sed 's/^/OTU/g' > OTUs_list.txt; paste -d" " OTUs_list.txt swarm_OTU_map_temp.txt | tr ' ' '\t' > swarm_OTU_map.txt
#	# create the files which are needed by LULU
#	pick_rep_set.py -i swarm_OTU_map.txt -f uniques_for_swarm.fasta -o 12S_swarm_otu.fas -m most_abundant
#	sed 's/ .*//g' 12S_swarm_otu.fas > temp.txt; mv temp.txt 12S_swarm_otu.fas
#	cat swarm_OTU_map.txt | tr '\t' '\n' | sed 's/_.*//g' | tr '\n' ' ' | sed 's/ OTU/#OTU/g' | tr '#' '\n' | tr ' ' '\t' > temp.txt; mv temp.txt swarm_OTU_map.txt
#	merge_otu_maps.py -i usearch_map.txt,swarm_OTU_map.txt -o 12S_swarm_otu_merge_map.txt
#	make_otu_table.py -i 12S_swarm_otu_merge_map.txt -o 12S_swarm_merge_otu.biom
#	biom convert -i 12S_swarm_merge_otu.biom -o 12S_swarm_merge_otu.txt --to-tsv
#	sed '1d' 12S_swarm_merge_otu.txt > temp.txt; sed 's/#OTU ID/OTU/g' temp.txt > 12S_swarm_merge_otu.txt
#	vsearch --usearch_global 12S_swarm_otu.fas --db 12S_swarm_otu.fas --self --id .84 --iddef 1 --userout match_list_12S.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10
#	# run LULU and get the final otu-table and otu-sequences files
#	cp "${HOMEFOLDER}"/12S_SWARM_LULU.R .
#	Rscript --vanilla --verbose 12S_SWARM_LULU.R
#	seqtk subseq 12S_swarm_otu.fas <(cut -f 1 12S_otu_table_swarm_lulu.txt) > 12S_otu_table_swarm_lulu.fas
#	mv 12S_otu_table_swarm_lulu.txt ../; mv 12S_otu_table_swarm_lulu.fas ../; cd ..; rm -rf swarm_lulu
	
	# 16S
	cd "${HOMEFOLDER}"/analysis/16S_DAMe_SORT_outputs/Filter_min2PCRs_min9copies_16S/
	mkdir swarm_lulu; cp 16S_seqs_nochimeras_82to150bp.fasta.gz swarm_lulu/; cd swarm_lulu; gunzip -d 16S_seqs_nochimeras_82to150bp.fasta.gz
	# SWARM needs the fasta input file which had better be dereplicated First
	# to obtain the otu-map file of dereplicating step, Pick_cluster_num_V2.pl is used and is written for the output uc file of USEARCH 5.2.236, if you use other versions of USEARCH, you have to rewrite Pick_cluster_num_V2.pl
	usearch5 -derep_fullseq -cluster 16S_seqs_nochimeras_82to150bp.fasta -seedsout uniques.fasta -sizeout -uc uniques.uc
	mkdir for_map; mv uniques.uc for_map/; perl "${HOMEFOLDER}"/Pick_cluster_num_V2.pl -id for_map -o usearch_map_temp.txt; sed 's/uniques_/uniques/g' usearch_map_temp.txt > usearch_map.txt
	sed 's/Cluster/uniques/g' uniques.fasta > temp.txt; sed 's/;size=/_/g' temp.txt > uniques_for_swarm.fasta
	# run SWARM
	${SWARM} -t 4 -f -o swarm_OTU_map_temp.txt uniques_for_swarm.fasta > /dev/null
	# modify swarm_OTU_map.txt to make it available for MACQIIME
	seq 1 $(grep -c "^" swarm_OTU_map_temp.txt) | sed 's/^/OTU/g' > OTUs_list.txt; paste -d" " OTUs_list.txt swarm_OTU_map_temp.txt | tr ' ' '\t' > swarm_OTU_map.txt
	# create the files which are needed by LULU
	pick_rep_set.py -i swarm_OTU_map.txt -f uniques_for_swarm.fasta -o 16S_swarm_otu.fas -m most_abundant
	sed 's/ .*//g' 16S_swarm_otu.fas > temp.txt; mv temp.txt 16S_swarm_otu.fas
	cat swarm_OTU_map.txt | tr '\t' '\n' | sed 's/_.*//g' | tr '\n' ' ' | sed 's/ OTU/#OTU/g' | tr '#' '\n' | tr ' ' '\t' > temp.txt; mv temp.txt swarm_OTU_map.txt
	merge_otu_maps.py -i usearch_map.txt,swarm_OTU_map.txt -o 16S_swarm_otu_merge_map.txt
	make_otu_table.py -i 16S_swarm_otu_merge_map.txt -o 16S_swarm_merge_otu.biom
	biom convert -i 16S_swarm_merge_otu.biom -o 16S_swarm_merge_otu.txt --to-tsv
	sed '1d' 16S_swarm_merge_otu.txt > temp.txt; sed 's/#OTU ID/OTU/g' temp.txt > 16S_swarm_merge_otu.txt
	vsearch --usearch_global 16S_swarm_otu.fas --db 16S_swarm_otu.fas --self --id .84 --iddef 1 --userout match_list_16S.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10
	# run LULU and get the final otu-table and otu-sequences files
	cp "${HOMEFOLDER}"/16S_SWARM_LULU.R .
	Rscript --vanilla --verbose 16S_SWARM_LULU.R
	seqtk subseq 16S_swarm_otu.fas <(cut -f 1 16S_otu_table_swarm_lulu.txt) > 16S_otu_table_swarm_lulu.fas
	mv 16S_otu_table_swarm_lulu.txt ../; mv 16S_otu_table_swarm_lulu.fas ../; cd ..; rm -rf swarm_lulu
