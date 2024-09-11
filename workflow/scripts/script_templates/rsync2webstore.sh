#!/bin/bash

# Push all customer files to the Strato HiDrive account so that they can be mounted to the 
# http://genomes.yseq.net/WGS/<YSEQID> web server


YSEQID=${PWD##*/}
# YSEQID="1234" # (the above command simply gets the name of the last segment of the current working directory)
WEBSTORE='socotra'
DISK='disk4'
WS_DIR="/webstore/${DISK}"

rsync -e "ssh" --progress -v .htaccess ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress -v .htpasswd ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v ${YSEQID}_1240K* ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v ${YSEQID}_result_summary.txt ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v ${YSEQID}_bwa-mem_hg38_sorted.bam.idxstats.tsv ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v ${YSEQID}_23andMe_all_hg19.zip ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v ${YSEQID}_mtDNA.fasta ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v chrY_called_${YSEQID}_*.vcf.gz ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v chrY_called_${YSEQID}_*.vcf.gz.tbi ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v chrY_cleaned_${YSEQID}_*.vcf.gz ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v chrY_cleaned_${YSEQID}_*.vcf.gz.tbi ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v chrY_derived_${YSEQID}_*.vcf.gz ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v chrY_derived_${YSEQID}_*.vcf.gz.tbi ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v chrY_INDELs_${YSEQID}_*.vcf.gz ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v chrY_INDELs_${YSEQID}_*.vcf.gz.tbi ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v chrY_novel_SNPs_${YSEQID}_*.vcf.gz ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v chrY_novel_SNPs_${YSEQID}_*.vcf.gz.tbi ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v chrY_novel_SNPs_${YSEQID}_*.ods ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v chrY_novel_SNPs_${YSEQID}_*.xls ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v ${YSEQID}_positives.txt ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v ${YSEQID}_negatives.txt ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v ${YSEQID}_cladeFinderOutput.csv ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v ${YSEQID}_bwa-mem_*_sorted.bam.bai ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v ${YSEQID}_*.vcf.gz ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v ${YSEQID}_*.vcf.gz.tbi ${WEBSTORE}:${WS_DIR}/${YSEQID}/
#rsync -e "ssh" --progress --size-only -v ${YSEQID}_bwa-mem_*_sorted.bam ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v ${YSEQID}_bwa-mem_*_chrY.bam ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v ${YSEQID}_bwa-mem_*_chrY.bam.bai ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v ${YSEQID}_bwa-mem_rCRS_chrM.bam ${WEBSTORE}:${WS_DIR}/${YSEQID}/
rsync -e "ssh" --progress --size-only -v ${YSEQID}_bwa-mem_rCRS_chrM.bam.bai ${WEBSTORE}:${WS_DIR}/${YSEQID}/

#rsync -e "ssh" --progress --size-only -v ${YSEQID}_bwa-mem_*_sorted.bam.stats ${WEBSTORE}:${WS_DIR}/${YSEQID}/
#rsync -e "ssh" --progress --size-only -v ${YSEQID}_lobSTR_aSTR.tsv ${WEBSTORE}:${WS_DIR}/${YSEQID}/
#rsync -e "ssh" --progress --size-only -v ${YSEQID}_lobSTR_YSTR.tsv ${WEBSTORE}:${WS_DIR}/${YSEQID}/


if [ -d "hg19" ]; then
  rsync -e "ssh" --progress --size-only --mkpath -v hg19/${YSEQID}_bwa-mem_*_sorted.bam* ${WEBSTORE}:${WS_DIR}/${YSEQID}/hg19/
  rsync -e "ssh" --progress --size-only --mkpath -v hg19/${YSEQID}_23andMe_all_hg19.zip ${WEBSTORE}:${WS_DIR}/${YSEQID}/hg19/
fi

if [ -d "CP086569.1" ]; then
  rsync -e "ssh" --progress --size-only --mkpath -v CP086569.1/${YSEQID}_*.bam* ${WEBSTORE}:${WS_DIR}/${YSEQID}/CP086569.1/
  rsync -e "ssh" --progress --size-only --mkpath -v CP086569.1/CP086569.1_*.vcf* ${WEBSTORE}:${WS_DIR}/${YSEQID}/CP086569.1/
  rsync -e "ssh" --progress --size-only --mkpath -v CP086569.1/CP086569.1_*.ods ${WEBSTORE}:${WS_DIR}/${YSEQID}/CP086569.1/
  #rsync -e "ssh" --progress --size-only --mkpath -v CP086569.1/CP086569.1_*.tsv ${WEBSTORE}:${WS_DIR}/${YSEQID}/CP086569.1/
fi

# Check for membership in the YFull group and set .htpasswd accordingly:
IS_YFULL_MEMBER="`wget -qO- https://www.yseq.net/yfull_group_check.php?cID=${YSEQID}`"
if [[ ${IS_YFULL_MEMBER} == *"1"* ]]; then
  echo "YSEQ ID ${YSEQID} is a member of the YFull group"
  wget -qO- "https://genomes.yseq.net/WGS/yfull_access/join_leave_yfull.py?yseqid=${YSEQID}&joinleave=1"
else 
  if [[ ${IS_YFULL_MEMBER} == *"0"* ]]; then
    echo "YSEQ ID ${YSEQID} is NOT a member of the YFull group"
    wget -qO- "https://genomes.yseq.net/WGS/yfull_access/join_leave_yfull.py?yseqid=${YSEQID}&joinleave=0"
  else  
    echo "ERROR: Couldn't determine if YSEQ ID ${YSEQID} is a member of the YFull group!"
  fi
fi

# Check for membership in the Y-DNA Warehouse group and set .htpasswd accordingly:
IS_YWH_MEMBER="`wget -qO- https://www.yseq.net/ydna_warehouse_group_check.php?cID=${YSEQID}`"
if [[ ${IS_YWH_MEMBER} == *"1"* ]]; then
  echo "YSEQ ID ${YSEQID} is a member of the Y-DNA Warehouse group"
  wget -qO- "https://genomes.yseq.net/WGS/warehouse_access/join_leave_warehouse.py?yseqid=${YSEQID}&joinleave=1"
else 
  if [[ ${IS_YWH_MEMBER} == *"0"* ]]; then
    echo "YSEQ ID ${YSEQID} is NOT a member of the Y-DNA Warehouse group"
    wget -qO- "https://genomes.yseq.net/WGS/warehouse_access/join_leave_warehouse.py?yseqid=${YSEQID}&joinleave=0"
  else  
    echo "ERROR: Couldn't determine if YSEQ ID ${YSEQID} is a member of the Y-DNA warehouse group!"
  fi
fi





