rule extract_all_SNPS:
    input:
        SORTED_BAM =    "results/{YSEQID}_bwa_mem_{REF}_sorted.bam",
        REFSEQ = "resources/refseq/{REF}/{REF}.fa"
    output:
        SNPS = "resources/snps_{REF}.vcf.gz",
        SNPS_TBI = "resources/snps_{REF}.vcf.gz.tbi",

        VCF = "results/chrY_raw_{YSEQID}_{REF}.vcf.gz"
    threads:
        workflow.cores * 1
    shell:
        """
        if [ "{wildcards.REF}" == "hg38" ]; then
            wget http://ybrowse.org/gbrowse2/gff/snps_hg38.vcf.gz
            wget http://ybrowse.org/gbrowse2/gff/snps_hg38.vcf.gz.tbi
        #elif [ "{wildcards.REF}" == "hs1" ]; then
        #   wget http://hs1.ybrowse.org/gbrowse2/gff/snp_hs1.gff3        
        else
            echo "Invalid reference: {wildcards.REF}"
            exit 1
        fi
        
        #get snips out of bam	
        bcftools mpileup -r chrY -C 0 --threads {threads} -O z -f {output.REF} {output.SORTED_BAM} > {output.VCF}
	    tabix {output.VCF}

        #merge all snps (HARRY/ALIEN) and sample
        bcftools merge -m both -O z snps_hg38.vcf.gz chrY_raw_${VCF_FILE}.gz > chrY_merged_${VCF_FILE}.gz
	    tabix chrY_merged_${VCF_FILE}.gz

        #get differences between HARRY/ALIEN and SAMPLE
	    bcftools call -O z -m -P 0 chrY_merged_${VCF_FILE}.gz > chrY_called_${VCF_FILE}.gz
	    tabix chrY_called_${VCF_FILE}.gz


        #remove HARRY/ALIEN from VCF
	    bcftools view -O z -k -s ^HARRY,ALIEN chrY_called_${VCF_FILE}.gz >chrY_cleaned_${VCF_FILE}.gz
	    tabix chrY_cleaned_${VCF_FILE}.gz

        #extraxt ancestral & derived SNPS to seperate files
	    bcftools filter -O z -i '(GT=="1/1" && AA==REF) || (GT=="0/0" && AA==ALT)' chrY_cleaned_${VCF_FILE}.gz >chrY_derived_${VCF_FILE}.gz &
	    bcftools filter -O z -i '(GT=="0/0" && AA==REF) || (GT=="1/1" && AA==ALT)' chrY_cleaned_${VCF_FILE}.gz >chrY_ancestral_${VCF_FILE}.gz &
	    wait
	    tabix chrY_derived_${VCF_FILE}.gz &
	    tabix chrY_ancestral_${VCF_FILE}.gz &
	    wait

        # Find the Y haplogroup
	    bcftools query -f '%ID,' chrY_derived_${VCF_FILE}.gz | sed ':a;N;$!ba;s/\n//g' > ${YSEQID}_positives.txt &
	    bcftools query -f '%ID,' chrY_ancestral_${VCF_FILE}.gz | sed ':a;N;$!ba;s/\n//g' > ${YSEQID}_negatives.txt &
	    wait
        """



