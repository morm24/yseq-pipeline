rule extract_all_SNPS:
    input:
        SORTED_BAM =    "results/{YSEQID}_bwa_mem_{REF}_sorted.bam",
        REFSEQ =        "resources/refseq/{REF}/{REF}.fa",
        SNPS =          "resources/tmp/snps_{REF}.vcf.gz",
        SNPS_TBI =      "resources/tmp/snps_{REF}.vcf.gz.tbi"
    output:        
        RAW_VCF =       temp("results/chrY_raw_{YSEQID}_{REF}.vcf.gz"),
        MERGED_VCF =    temp("results/chrY_merged_{YSEQID}_{REF}.vcf.gz"),
        CALLED_VCF =    temp("results/chrY_called_{YSEQID}_{REF}.vcf.gz"),
        CLEANED_VCF =   temp("results/chrY_cleaned_{YSEQID}_{REF}.vcf.gz"),
        

        DERIVED_VCF =   "results/chrY_derived_{YSEQID}_{REF}.vcf.gz",
        ANCESTRAL_VCF = "results/chrY_ancestral_{YSEQID}_{REF}.vcf.gz",
        POSITIVE_TXT =  "results/{YSEQID}_{REF}_positives.txt",
        NEGATIVE_TXT =  "results/{YSEQID}_{REF}_negatives.txt"



    threads:
        workflow.cores * 1
    shell:
        """

        #get all SNPs from Sample
        bcftools mpileup -r chrY -C 0 --threads {threads} -O z -f {input.REFSEQ} {input.SORTED_BAM} > {output.RAW_VCF}
	    tabix {output.RAW_VCF}

        #merge all snps we have (HARRY/ALIEN) and sample (HARRY ~ No SNPs; ALIAN ~ All SNPs)
        bcftools merge -m both -O z {input.SNPS} {output.RAW_VCF} > {output.MERGED_VCF}
	    tabix {output.MERGED_VCF}

        #get differences between HARRY/ALIEN and SAMPLE
	    bcftools call -O z -m -P 0 {output.MERGED_VCF} > {output.CALLED_VCF} 
	    tabix {output.CALLED_VCF}


        #remove HARRY/ALIEN from VCF
	    bcftools view -O z -k -s ^HARRY,ALIEN {output.CALLED_VCF} > {output.CLEANED_VCF}
	    tabix {output.CLEANED_VCF}

        #extraxt ancestral & derived SNPS to seperate files
	    bcftools filter -O z -i '(GT=="1/1" && AA==REF) || (GT=="0/0" && AA==ALT)' {output.CLEANED_VCF} > {output.DERIVED_VCF} &
	    bcftools filter -O z -i '(GT=="0/0" && AA==REF) || (GT=="1/1" && AA==ALT)' {output.CLEANED_VCF} > {output.ANCESTRAL_VCF} &
	    wait
	    tabix {output.DERIVED_VCF} &
	    tabix {output.ANCESTRAL_VCF} &
	    wait

        # Find the Y haplogroup
	    bcftools query -f '%ID,' {output.DERIVED_VCF} | sed ':a;N;$!ba;s/\n//g' > {output.POSITIVE_TXT} &
	    bcftools query -f '%ID,' {output.ANCESTRAL_VCF} | sed ':a;N;$!ba;s/\n//g' > {output.NEGATIVE_TXT} &
	    wait
        """



