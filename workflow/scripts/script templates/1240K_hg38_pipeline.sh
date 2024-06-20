#!/bin/bash
START=$(date +%s.%N)
#clear
# setup parameters

echo "This script needs the plink2 executable"
cp /genomes/0/script-templates/plink2 .

PLINK="./plink2"

YSEQID=${PWD##*/}
# YSEQID="1234" # (the above command simply gets the name of the last segment of the current working directory)

NUM_THREADS=$(getconf _NPROCESSORS_ONLN)
echo "We can use ${NUM_THREADS} threads."

REF_HG19="/genomes/0/refseq/hg19/hg19.fa"
REF_HG38="/genomes/0/refseq/hg38/hg38.fa"

BAMFILE_SORTED="${YSEQID}_bwa-mem_hg38_sorted.bam"
VCF_FILE="${YSEQID}_1240K_hg38.vcf"

# REF_23ANDME="23andMe_all_hg19_ref.tab"
# REF_23ANDME_HG38="23andMe_all_hg38_ref.tab"
REF_1240K="1240K_all_hg19_ref.tab.gz"
REF_1240K_HG38="1240K_hg38.tab.gz"


TEMPLATE_1240K="/genomes/0/refseq/hg19/1240K_hg19.tab.gz"

if [ -f "${TEMPLATE_1240K}" ]; then
    echo "${TEMPLATE_1240K} already stored. Using it. Make sure the file is up to date."
    cp ${TEMPLATE_1240K} .
    zcat ${TEMPLATE_1240K} > 1240K_all_hg19_raw.tab
else 
    echo "${TEMPLATE_1240K} does not exist. BREAK!"
    exit 1
fi
echo "1240k SNP definitions available" 


echo "#CHROM	POS	ID" > 1240K_all_hg19_ref.tab

while IFS=$'\t' read -r chromosome chromosome_position snp
do
	if [[ $index == \#* || $index =~ ^index ]]; then
		echo "skipping $index"
	else 
		CHROM=${chromosome//MT/M}
		echo "${CHROM}	${chromosome_position}	${snp}	--" >> 1240K_all_hg19_unsorted.tab
	fi
done < 1240K_all_hg19_raw.tab


sort -t $'\t' -k1,2 -V 1240K_all_hg19_unsorted.tab > 1240K_all_hg19_sorted.tab
cat 1240K_all_hg19_sorted.tab >> 1240K_all_hg19_ref.tab


# Convert the 1240K TAB file into a VCF
bcftools convert -c CHROM,POS,ID,AA -s SampleName -f ${REF_HG19} --tsv2vcf 1240K_all_hg19_ref.tab -Ov -o 1240K_all_hg19_ref.vcf


bgzip -c 1240K_all_hg19_ref.tab > 1240K_all_hg19_ref.tab.gz
tabix -s1 -b2 -e2 1240K_all_hg19_ref.tab.gz


# CrossMap the 1240K template VCF to hg38
map_chain="/genomes/0/refseq/chain/hg19ToHg38.over.chain"
CrossMap.py vcf $map_chain 1240K_all_hg19_ref.vcf ${REF_HG38} 1240K_all_hg38_ref.vcf

# Sort again (due to inversions during hg19 > hg38 conversion)
bcftools sort -Oz 1240K_all_hg38_ref.vcf -o 1240K_all_hg38_ref_sorted.vcf.gz
tabix -p vcf 1240K_all_hg38_ref_sorted.vcf.gz

# Convert back to TSV
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_all_hg38_ref.tab.gz
tabix -s1 -b2 -e2 1240K_all_hg38_ref.tab.gz
echo "1240K SNP definitions translated to hg38"

# Split in separate chr files for speedup of mpileup
bcftools query -r chr1  -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr1_hg38_ref.tab.gz &
bcftools query -r chr2  -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr2_hg38_ref.tab.gz &
bcftools query -r chr3  -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr3_hg38_ref.tab.gz &
bcftools query -r chr4  -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr4_hg38_ref.tab.gz &
bcftools query -r chr5  -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr5_hg38_ref.tab.gz &
bcftools query -r chr6  -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr6_hg38_ref.tab.gz &
bcftools query -r chr7  -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr7_hg38_ref.tab.gz &
bcftools query -r chr8  -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr8_hg38_ref.tab.gz &
bcftools query -r chr9  -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr9_hg38_ref.tab.gz &
bcftools query -r chr10 -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr10_hg38_ref.tab.gz &
bcftools query -r chr11 -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr11_hg38_ref.tab.gz &
bcftools query -r chr12 -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr12_hg38_ref.tab.gz &
bcftools query -r chr13 -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr13_hg38_ref.tab.gz &
bcftools query -r chr14 -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr14_hg38_ref.tab.gz &
bcftools query -r chr15 -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr15_hg38_ref.tab.gz &
bcftools query -r chr16 -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr16_hg38_ref.tab.gz &
bcftools query -r chr17 -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr17_hg38_ref.tab.gz &
bcftools query -r chr18 -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr18_hg38_ref.tab.gz &
bcftools query -r chr19 -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr19_hg38_ref.tab.gz &
bcftools query -r chr20 -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr20_hg38_ref.tab.gz &
bcftools query -r chr21 -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr21_hg38_ref.tab.gz &
bcftools query -r chr22 -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chr22_hg38_ref.tab.gz &
bcftools query -r chrX -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chrX_hg38_ref.tab.gz &
bcftools query -r chrY -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chrY_hg38_ref.tab.gz &
#bcftools query -r chrM -f'%CHROM\t%POS\t%REF,%ALT\n' 1240K_all_hg38_ref_sorted.vcf.gz | bgzip -c > 1240K_chrM_hg38_ref.tab.gz &
wait

tabix -s1 -b2 -e2 1240K_chr1_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr2_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr3_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr4_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr5_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr6_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr7_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr8_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr9_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr10_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr11_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr12_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr13_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr14_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr15_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr16_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr17_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr18_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr19_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr20_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr21_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chr22_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chrX_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 1240K_chrY_hg38_ref.tab.gz &
#tabix -s1 -b2 -e2 1240K_chrM_hg38_ref.tab.gz &
wait


# Generate 1240K  file

# Parallel SNP calling by chromosome
echo "MpileUp started"

PARAMC=0

bcftools mpileup -R 1240K_chr1_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr1_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr2_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr2_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr3_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr3_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr4_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr4_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr5_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr5_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr6_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr6_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr7_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr7_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr8_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr8_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr9_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr9_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr10_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr10_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr11_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr11_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr12_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr12_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr13_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr13_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr14_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr14_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr15_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr15_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr16_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr16_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr17_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr17_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr18_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr18_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr19_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr19_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr20_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr20_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr21_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr21_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chr22_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr22_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chrX_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chrX_${VCF_FILE}.gz &
bcftools mpileup -R 1240K_chrY_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chrY_${VCF_FILE}.gz &
bcftools mpileup -r chrM                        -Ou -C 50 -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chrM_${VCF_FILE}.gz &
wait

# Concatenate all chromosome VCFs to one big file
bcftools concat -O z chr[1-9]_${VCF_FILE}.gz chr[1-2][0-9]_${VCF_FILE}.gz chr[M,X-Y]_${VCF_FILE}.gz > ${YSEQID}_1240K_called_hg38_vcf.gz
tabix -p vcf ${YSEQID}_1240K_called_hg38_vcf.gz

# Delete no longer needed VCFs
rm -f chr[0-9]_${VCF_FILE}.gz chr[1-2][0-9]_${VCF_FILE}.gz chrX_${VCF_FILE}.gz

echo "mpileup complete"

# Reverse translate hg38 > hg19
back_chain="/genomes/0/refseq/chain/hg38ToHg19.over.chain"
CrossMap.py vcf $back_chain ${YSEQID}_1240K_called_hg38_vcf.gz ${REF_HG19} ${YSEQID}_1240K_translated_hg19.vcf

echo "Translating SNP calls to hg19 complete"

# Sort again (due to inversions during hg38 > hg19 conversion) & 
bcftools sort -Oz ${YSEQID}_1240K_translated_hg19.vcf -o ${YSEQID}_1240K_sorted_hg19.vcf.gz
tabix -p vcf ${YSEQID}_1240K_sorted_hg19.vcf.gz

echo "Sorting hg19 SNPs complete"

# Filter (-R) for tested chip positions only
bcftools view -Oz -R 1240K_all_hg19_ref.tab.gz ${YSEQID}_1240K_sorted_hg19.vcf.gz -o ${YSEQID}_1240K_sorted_filtered_hg19.vcf.gz
tabix -p vcf ${YSEQID}_1240K_sorted_filtered_hg19.vcf.gz

echo "Filtering hg19 SNPs complete"


# Annotate the SNP names (rs numbers)
bcftools annotate -Oz -a 1240K_all_hg19_ref.tab.gz -c CHROM,POS,ID ${YSEQID}_1240K_sorted_filtered_hg19.vcf.gz > ${YSEQID}_1240K_annotated_hg19.vcf.gz
tabix -p vcf ${YSEQID}_1240K_annotated_hg19.vcf.gz

echo "SNP name annotation complete"


if command -v ${PLINK} &> /dev/null
then
  bcftools norm -Ou -m -any ${YSEQID}_1240K_annotated_hg19.vcf.gz | \
  bcftools norm -Ou -f ${REF_HG19} -cs | \
  bcftools annotate -Ob -I +'%CHROM:%POS:%REF:%ALT' > ${YSEQID}_1240K_hg19.bcf.gz
  echo "${YSEQID}_1240K_hg19.bcf.gz created. Running PLINK:" 
  ${PLINK} --bcf ${YSEQID}_1240K_hg19.bcf.gz --const-fid --allow-extra-chr 0 --make-bed --out ${YSEQID}_1240K_plink_hg19
fi

# Cleanup no longer needed files
echo "Cleaning up ..."

rm -f 1240K_chr*ref.tab*
rm -f *${YSEQID}_1240K_hg38.vcf.gz
rm -f ${YSEQID}_1240K_all_hg19_sorted.tab
rm -f ${YSEQID}_1240K_all_hg19.tab
rm -f ${YSEQID}_1240K_all_hg19.txt
#rm -f ${YSEQID}_1240K_annotated_hg19.vcf.gz
#rm -f ${YSEQID}_1240K_annotated_hg19.vcf.gz.tbi
#rm -f ${YSEQID}_1240K_called_hg38_vcf.gz
#rm -f ${YSEQID}_1240K_called_hg38_vcf.gz.tbi
rm -f ${YSEQID}_1240K_sorted_filtered_hg19.vcf.gz
rm -f ${YSEQID}_1240K_sorted_filtered_hg19.vcf.gz.tbi
rm -f ${YSEQID}_1240K_sorted_hg19.vcf.gz
rm -f ${YSEQID}_1240K_sorted_hg19.vcf.gz.tbi
rm -f ${YSEQID}_1240K_translated_hg19.vcf
rm -f ${YSEQID}_1240K_translated_hg19.vcf.unmap
rm -f 1240K_all*

END=$(date +%s.%N)

DIFF=$(echo "$END - $START" | bc)
echo "Finished in ${DIFF} Seconds"





