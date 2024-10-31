#!/bin/bash
START=$(date +%s.%N)
#clear
# setup parameters

YSEQID=${PWD##*/}
# YSEQID="1234" # (the above command simply gets the name of the last segment of the current working directory)

NUM_THREADS=$(getconf _NPROCESSORS_ONLN)
echo "We can use ${NUM_THREADS} threads."

REF="/genomes/0/refseq/hg38/hg38.fa"

BAMFILE_SORTED="${YSEQID}_bwa-mem_hg38_sorted.bam"


# Get snp.list filename from command line parameter
# Needs to be a tab separated file in the format:
# <chr>\t<position>\t<snp_name>\t<anc>t\<der>\n

SNP_LIST="/genomes/0/refseq/hg38/phenotyping_snp_list.tab.csv"

echo "Using SNP_LIST = ${SNP_LIST}"

SNP_LIST_NAME=${SNP_LIST##*/}
echo "#CHROM	POS	ID	ANC	DER" >${SNP_LIST_NAME}_sorted.tab
sort -t $'\t' -k1,2 -V ${SNP_LIST} >> ${SNP_LIST_NAME}_sorted.tab

bgzip -c ${SNP_LIST_NAME}_sorted.tab > ${SNP_LIST_NAME}_sorted.gz
tabix -s1 -b2 -e2 ${SNP_LIST_NAME}_sorted.gz


bcftools mpileup -C 0 -R ${SNP_LIST_NAME}_sorted.gz -T ${SNP_LIST_NAME}_sorted.gz -f $REF $BAMFILE_SORTED > ${YSEQID}_${SNP_LIST_NAME}_raw.vcf.gz
tabix -p vcf ${YSEQID}_${SNP_LIST_NAME}_raw.vcf.gz

bcftools call -O z -V indels -m -P 0 ${YSEQID}_${SNP_LIST_NAME}_raw.vcf.gz > ${YSEQID}_${SNP_LIST_NAME}_called.vcf.gz
tabix -p vcf ${YSEQID}_${SNP_LIST_NAME}_called.vcf.gz

bcftools annotate -O z -a ${SNP_LIST_NAME}_sorted.gz -c CHROM,POS,ID ${YSEQID}_${SNP_LIST_NAME}_called.vcf.gz > ${YSEQID}_${SNP_LIST_NAME}_annotated.vcf.gz
tabix -p vcf ${YSEQID}_${SNP_LIST_NAME}_annotated.vcf.gz

bcftools query -f "${YSEQID}\t%ID[\t%IUPACGT]\t%CHROM\t%POS\n" ${YSEQID}_${SNP_LIST_NAME}_annotated.vcf.gz > ${YSEQID}_${SNP_LIST_NAME}

INDICATOR=""
rm -f ${YSEQID}_phenotyping.alx
while read -r SAMPLE RSID ALLELE CHROM POS; do
    # Map the SNP information
    #echo "Next Line: $SAMPLE,$RSID,$ALLELE,$CHROM,$POS"
    while read -r M_CHROM M_POS M_RSID M_ANC M_DER; do
    	#echo "$M_CHROM = $CHROM"
    	#echo "$M_POS = $POS"
    	#echo "$ALLELE = $M_ANC OR $M_DER"
    	
		INDICATOR=""
    	if [[ "$ALLELE" == "$M_ANC" ]] 
    		then INDICATOR="-"
    	elif [[ "$ALLELE" == "$M_DER" ]]
    		then INDICATOR="+"
    	elif [[ "ins" == "$M_DER" ]]
    		then INDICATOR=" Check InDel at $CHROM:$POS!"
    		ALLELE="?"
    	elif [[ "del" == "$M_DER" ]]
    		then INDICATOR=" Check InDel at $CHROM:$POS!"
    		ALLELE="?" 		
    	fi
    	
    	echo "$SAMPLE	$RSID	$ALLELE$INDICATOR	${YSEQID}_${SNP_LIST_NAME}" >>${YSEQID}_phenotyping.alx
    done < <(tabix ${SNP_LIST_NAME}_sorted.gz $CHROM:$POS-$POS)
    
done <${YSEQID}_${SNP_LIST_NAME}


echo "Phenotyping alleles have been written to ${YSEQID}_phenotyping.alx"




#cleanup
#rm -f ${SNP_LIST_NAME}_sorted*
#rm -f ${YSEQID}_phenotyping_*raw*
#rm -f ${YSEQID}_phenotyping_*called*
#rm -f ${YSEQID}_phenotyping_*annotated*

END=$(date +%s.%N)

DIFF=$(echo "$END - $START" | bc)
echo ${DIFF}

