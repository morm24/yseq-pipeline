rule download_snps:
    output:
        SNPS =          temp("resources/tmp/snps_{REF}.vcf.gz"),
        SNPS_TBI =      temp("resources/tmp/snps_{REF}.vcf.gz.tbi")
    shell:
        """
        if [ "{wildcards.REF}" == "hg38" ]; then
            wget -O {output[0]} http://ybrowse.org/gbrowse2/gff/snps_hg38.vcf.gz
            wget -O {output[1]} http://ybrowse.org/gbrowse2/gff/snps_hg38.vcf.gz.tbi
        #elif [ "{wildcards.REF}" == "hs1" ]; then
            #wget -O {output.SNPS} http://ybrowse.org/gbrowse2/gff/snps_hs1.vcf.gz                  still needs to be created
            #wget -O {output.SNPS_TBI} http://ybrowse.org/gbrowse2/gff/snps_hs1.vcf.gz.tbi          still needs to be created
        else
            echo "Invalid reference: {wildcards.REF}"
            exit 1
        fi
        """

rule get_all_SNPs_Sample:
    input:
        SORTED_BAM =    "results/{YSEQID}_bwa_mem_{REF}_sorted.bam",
        REFSEQ =        "resources/refseq/{REF}/{REF}.fa"
    output:
        RAW_VCF =       temp("results/chrY_raw_{YSEQID}_{REF}.vcf.gz")
    threads:
        workflow.cores * 1
    shell:
        """
        bcftools mpileup -r chrY -C 0 --threads {threads} -O z -f {input.REFSEQ} {input.SORTED_BAM} > {output.RAW_VCF}
	    tabix {output.RAW_VCF}
        """

#Human Ancestral Reconstructed Reference of th Y_chromosome; ALIEN ~ ALL Snips are devived allele
rule merge_SNPS_HARRY_ALIEN_SAMPLE:
    input:
        SNPS =          "resources/tmp/snps_{REF}.vcf.gz",
        SNPS_TBI =      "resources/tmp/snps_{REF}.vcf.gz.tbi",
        RAW_VCF =       "results/chrY_raw_{YSEQID}_{REF}.vcf.gz"
    output:
        MERGED_VCF =    temp("results/chrY_merged_{YSEQID}_{REF}.vcf.gz")
    shell:
        """
        bcftools merge -m both -O z {input.SNPS} {input.RAW_VCF} > {output.MERGED_VCF}
	    tabix {output.MERGED_VCF}
        """

rule get_differences_HARRY_ALIEN_SAMPLE:
    input:
        MERGED_VCF =    "results/chrY_merged_{YSEQID}_{REF}.vcf.gz"
    output:
        CALLED_VCF =    temp("results/chrY_called_{YSEQID}_{REF}.vcf.gz")
    shell:
        """
	    bcftools call -O z -m -P 0 {input.MERGED_VCF} > {output.CALLED_VCF} 
	    tabix {output.CALLED_VCF}
        """

rule rm_HARRY_ALIEN_from_VCF:
    input:
        CALLED_VCF =    "results/chrY_called_{YSEQID}_{REF}.vcf.gz"
    output:
        CLEANED_VCF =   temp("results/chrY_cleaned_{YSEQID}_{REF}.vcf.gz")
    shell:
        """
	    bcftools view -O z -k -s ^HARRY,ALIEN {input.CALLED_VCF} > {output.CLEANED_VCF}
	    tabix {output.CLEANED_VCF}
        """

rule extract_ancestral_SNPS:
    input:
        CLEANED_VCF =   "results/chrY_cleaned_{YSEQID}_{REF}.vcf.gz"
    output:
        DERIVED_VCF =   "results/chrY_derived_{YSEQID}_{REF}.vcf.gz"
    threads:
        workflow.cores * 0.5
    shell:
        """
        bcftools filter --threads {threads} -O z -i '(GT=="1/1" && AA==REF) || (GT=="0/0" && AA==ALT)' {input.CLEANED_VCF} > {output.DERIVED_VCF}
        tabix {output.DERIVED_VCF}
        """

rule extract_derived_SNPS:
    input:
        CLEANED_VCF =   "results/chrY_cleaned_{YSEQID}_{REF}.vcf.gz"
    output:
        ANCESTRAL_VCF = "results/chrY_ancestral_{YSEQID}_{REF}.vcf.gz"
    threads:
        workflow.cores * 0.5
    shell:
        """
	    bcftools filter --threads {threads} -O z -i '(GT=="0/0" && AA==REF) || (GT=="1/1" && AA==ALT)' {input.CLEANED_VCF} > {output.ANCESTRAL_VCF} 
        tabix {output.ANCESTRAL_VCF}
        """


