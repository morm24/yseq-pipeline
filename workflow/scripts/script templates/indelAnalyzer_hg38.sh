#!/bin/bash

START=$(date +%s.%N)

YSEQID=${PWD##*/}
# YSEQID="1234" # (the above command simply gets the name of the last segment of the current working directory)

REF="/genomes/0/refseq/hg38/hg38.fa"
BAMFILE_SORTED="${YSEQID}_bwa-mem_hg38_sorted.bam"
AUTO_BED="/genomes/0/refseq/hg38/indelAnalyzer_Autosomal_STRs_hg38.bed"

mkdir -p indelAnalyzer_files
rm -f indelAnalyzer_processBatch.pl
rm -f indelAnalyzer_files/*

cp /genomes/0/script-templates/indelAnalyzer/indelAnalyzer_processBatch.pl .
cp /genomes/0/script-templates/indelAnalyzer/*.js indelAnalyzer_files/
cp /genomes/0/script-templates/indelAnalyzer/test_output_v6 indelAnalyzer_files/
cp /genomes/0/script-templates/indelAnalyzer/*.png indelAnalyzer_files/

perl indelAnalyzer_processBatch.pl -i $YSEQID -r $REF -s $BAMFILE_SORTED -chr chrY -batchLabel Y-STRs -testResultsFile indelAnalyzer_files/test_output_v6 -linkTools
perl indelAnalyzer_processBatch.pl -i $YSEQID -r $REF -s $BAMFILE_SORTED -chr chrM -batchLabel mtdna
perl indelAnalyzer_processBatch.pl -i $YSEQID -r $REF -s $BAMFILE_SORTED -chr all -batchLabel autosomal -definitionsFile $AUTO_BED

END=$(date +%s.%N)

DIFF=$(echo "$END - $START" | bc)

echo ${DIFF}




