#download all known snps for HS1 or HG38, updated dayly
rule download_snps:
    output:
        SNPS =          results_prefix  / "snp_calling" / "snps_{REF}.vcf.gz",
        SNPS_TBI =      results_prefix  / "snp_calling" / "snps_{REF}.vcf.gz.tbi"
    conda:
        "../envs/get_sample_snps.yaml"
    log:
        results_prefix / "snp_calling" / "logs" / "load_snps_{REF}.log"
    benchmark:
        results_prefix / "snp_calling" / "benchmark" / "load_snps_{REF}.benchmark"
    shell:
        """
        (
        if [ "{wildcards.REF}" == "hg38" ]; then
            wget -O {output[0]} http://ybrowse.org/gbrowse2/gff/snps_hg38.vcf.gz
            wget -O {output[1]} http://ybrowse.org/gbrowse2/gff/snps_hg38.vcf.gz.tbi
        elif [ "{wildcards.REF}" == "hs1" ]; then 
            wget -O {output.SNPS} http://hs1.ybrowse.org/gbrowse2/gff/snps_hs1.vcf.gz         
            wget -O {output.SNPS_TBI} http://hs1.ybrowse.org/gbrowse2/gff/snps_hs1.vcf.gz.tbi 
        else
            echo "Invalid reference: {wildcards.REF}"
            exit 1
        fi
        ) > {log} 2>&1
        """

#extract all SNPs from the sample
rule get_all_SNPs_Sample:
    input:
        SORTED_BAM =    results_prefix / "mapping" / "{YSEQID}_bwa-mem_{REF}_sorted.bam",
        REFSEQ =        ref_prefix / "{REF}/{REF}.fa"
    output:
        VCF_TBI =       temp(results_prefix  / "snp_calling" / "chrY_raw_{YSEQID}_{REF}.vcf.gz.tbi"),
        RAW_VCF =       temp(results_prefix  / "snp_calling" / "chrY_raw_{YSEQID}_{REF}.vcf.gz")
    conda:
        "../envs/get_sample_snps.yaml"

    log:
        results_prefix / "snp_calling" / "logs" / "{YSEQID}_{REF}_bcf_mpileup.log"
    benchmark:
        results_prefix / "snp_calling" / "benchmark" / "{YSEQID}_{REF}_bcf_mpileup.benchmark"
    threads:
        workflow.cores * 1
    shell:
        """
       ( bcftools mpileup -r chrY -C 0 --threads {threads} -O z -f {input.REFSEQ} {input.SORTED_BAM} > {output.RAW_VCF}) > {log} 2>&1
	   ( tabix {output.RAW_VCF} ) >> {log} 2>&1
        """


rule merge_SNPS_HARRY_ALIEN_SAMPLE:
    input:
        SNPS =          results_prefix  / "snp_calling" / "snps_{REF}.vcf.gz",
        SNPS_TBI =      results_prefix  / "snp_calling" / "snps_{REF}.vcf.gz.tbi",
        RAW_VCF =       results_prefix  / "snp_calling" / "chrY_raw_{YSEQID}_{REF}.vcf.gz",
        VCF_TBI =       results_prefix  / "snp_calling" / "chrY_raw_{YSEQID}_{REF}.vcf.gz.tbi"
    output:
        MERGED_VCF =    temp(results_prefix  / "snp_calling" / "chrY_merged_{YSEQID}_{REF}.vcf.gz"),
        VCF_TBI =       temp(results_prefix  / "snp_calling" / "chrY_merged_{YSEQID}_{REF}.vcf.gz.tbi")
    conda:
        "../envs/get_sample_snps.yaml"
    log:
        results_prefix / "snp_calling" / "logs" / "{YSEQID}_{REF}_bcf_merge.log"
    benchmark:
        results_prefix / "snp_calling" / "benchmark" / "{YSEQID}_{REF}_bcf_merge.benchmark"

    shell:
        """
        (bcftools merge -m both -O z {input.SNPS} {input.RAW_VCF} > {output.MERGED_VCF}) > {log} 2>&1
	    tabix {output.MERGED_VCF} >> {log} 2>&1
        """

rule get_differences_HARRY_ALIEN_SAMPLE:
    input:
        MERGED_VCF =    results_prefix  / "snp_calling" / "chrY_merged_{YSEQID}_{REF}.vcf.gz",
        VCF_TBI =       results_prefix  / "snp_calling" / "chrY_merged_{YSEQID}_{REF}.vcf.gz.tbi"
    output:
        CALLED_VCF =    temp(results_prefix  / "snp_calling" / "chrY_called_{YSEQID}_{REF}.vcf.gz"),
        VCF_TBI =       temp(results_prefix  / "snp_calling" / "chrY_called_{YSEQID}_{REF}.vcf.gz.tbi")
    conda:
        "../envs/get_sample_snps.yaml"
    log:
        results_prefix / "snp_calling" / "logs" / "{YSEQID}_{REF}_bcf_call.log"
    benchmark:
        results_prefix / "snp_calling" / "benchmark" / "{YSEQID}_{REF}_bcf_call.benchmark"

    shell:
        """
	    (bcftools call -O z -m -P 0 {input.MERGED_VCF} > {output.CALLED_VCF} ) > {log} 2>&1
	    tabix {output.CALLED_VCF} >> {log} 2>&1
        """

rule rm_HARRY_ALIEN_from_VCF:
    input:
        CALLED_VCF =    results_prefix  / "snp_calling" / "chrY_called_{YSEQID}_{REF}.vcf.gz",
        VCF_TBI =       results_prefix  / "snp_calling" / "chrY_called_{YSEQID}_{REF}.vcf.gz.tbi"
    output:
        #CLEANED_VCF =   temp(results_prefix  / "snp_calling" / "chrY_cleaned_{YSEQID}_{REF}.vcf.gz"),
        CLEANED_VCF =   temp(results_prefix  / "snp_calling" / "chrY_cleaned_{YSEQID}_{REF}.vcf.gz"),
        VCF_TBI =       temp(results_prefix  / "snp_calling" / "chrY_cleaned_{YSEQID}_{REF}.vcf.gz.tbi")
    conda:
        "../envs/get_sample_snps.yaml"
    log:
        results_prefix / "snp_calling" / "logs" / "{YSEQID}_{REF}_bcf__view_rm_HARRY_ALIEN.log"
    benchmark:
        results_prefix / "snp_calling" / "benchmark" / "{YSEQID}_{REF}_bcf__view_rm_HARRY_ALIEN.benchmark"

    shell:
        """
	    (bcftools view -O z -k -s ^HARRY,ALIEN {input.CALLED_VCF} > {output.CLEANED_VCF}) > {log} 2>&1
	    tabix {output.CLEANED_VCF} >> {log} 2>&1
        """

rule extract_derived_SNPS:
    input:
        CLEANED_VCF =   results_prefix  / "snp_calling" / "chrY_cleaned_{YSEQID}_{REF}.vcf.gz",
        VCF_TBI =       results_prefix  / "snp_calling" / "chrY_cleaned_{YSEQID}_{REF}.vcf.gz.tbi"
    output:
        DERIVED_VCF =   results_prefix  / "snp_calling" / "chrY_derived_{YSEQID}_{REF}.vcf.gz",
        VCF_TBI =       results_prefix  / "snp_calling" / "chrY_derived_{YSEQID}_{REF}.vcf.gz.tbi"
    conda:
        "../envs/get_sample_snps.yaml"
    log:
       results_prefix / "snp_calling" / "logs" / "{YSEQID}_{REF}_bcf_filter_ancestral.log"
    benchmark:
        results_prefix / "snp_calling" / "benchmark" / "{YSEQID}_{REF}_bcf_filter_ancestral.benchmark"

    threads:
        workflow.cores * 0.5
    shell:
        """
        (bcftools filter --threads {threads} -O z -i '(GT=="1/1" && AA==REF) || (GT=="0/0" && AA==ALT)' {input.CLEANED_VCF} > {output.DERIVED_VCF}) > {log} 2>&1
        tabix {output.DERIVED_VCF} >> {log} 2>&1
        """

rule extract_ancestral_SNPS:
    input:
        CLEANED_VCF =   results_prefix  / "snp_calling" / "chrY_cleaned_{YSEQID}_{REF}.vcf.gz",
        VCF_TBI =       results_prefix  / "snp_calling" / "chrY_cleaned_{YSEQID}_{REF}.vcf.gz.tbi"

    output:
        ANCESTRAL_VCF = results_prefix  / "snp_calling" / "chrY_ancestral_{YSEQID}_{REF}.vcf.gz",
        VCF_TBI =       results_prefix  / "snp_calling" / "chrY_ancestral_{YSEQID}_{REF}.vcf.gz.tbi"
    conda:
        "../envs/get_sample_snps.yaml"
    log:
        results_prefix / "snp_calling" / "logs" / "{YSEQID}_{REF}_bcf_filter_derived.log"
    benchmark:
        results_prefix / "snp_calling" / "benchmark" / "{YSEQID}_{REF}_bcf_filter_derived.benchmark"

    threads:
        workflow.cores * 0.5
    shell:
        """
	    (bcftools filter --threads {threads} -O z -i '(GT=="0/0" && AA==REF) || (GT=="1/1" && AA==ALT)' {input.CLEANED_VCF} > {output.ANCESTRAL_VCF} )  > {log} 2>&1
        tabix {output.ANCESTRAL_VCF} >> {log} 2>&1
        """

rule get_novel_SNPS:
    input:
        CALLED_VCF =    results_prefix  / "snp_calling" / "chrY_called_{YSEQID}_{REF}.vcf.gz",
        VCF_TBI =       results_prefix  / "snp_calling" / "chrY_called_{YSEQID}_{REF}.vcf.gz.tbi"
    output:
        NOVEL_VCF =     results_prefix  / "snp_calling" / "chrY_novel_SNPs_{YSEQID}_{REF}.vcf.gz",
        NOVEL_TBI =     results_prefix  / "snp_calling" / "chrY_novel_SNPs_{YSEQID}_{REF}.vcf.gz.tbi",
        NOVEL_VCF_TSV =     results_prefix  / "snp_calling" / "chrY_novel_SNPs_{YSEQID}_{REF}.vcf.tsv"
    conda:
        "../envs/get_sample_snps.yaml"
    log:
        results_prefix / "snp_calling" / "logs" / "{YSEQID}_{REF}_bcf_filter_novel.log"
    benchmark:
        results_prefix / "snp_calling" / "benchmark" / "{YSEQID}_{REF}_bcf_filter_novel.benchmark"

    shell:
        """
        (bcftools filter -i 'TYPE="snp" && QUAL>=30 && GT=="1/1" && DP4[2] + DP4[3] >=2 && (DP4[0] + DP4[1]) / DP < 0.1' {input.CALLED_VCF} | bcftools view -n -O z -s ^HARRY,ALIEN > {output.NOVEL_VCF})   > {log} 2>&1
        tabix {output.NOVEL_VCF} >> {log} 2>&1
        (zcat {output.NOVEL_VCF} | grep -v "##" >{output.NOVEL_VCF_TSV} ) >> {log} 2>&1
        """


rule confirm_novel_SNPS:
    input:
        NOVEL_VCF_TSV =     results_prefix  / "snp_calling" / "chrY_novel_SNPs_{YSEQID}_{REF}.vcf.tsv",
        REFSEQ =            ref_prefix / "{REF}/{REF}.fa",
        MM2_INDEX =         ref_prefix / "{REF}/{REF}.fa.mmi"
    output:
        NOVEL_TSV =     results_prefix  / "snp_calling" / "chrY_novel_SNPs_{YSEQID}_{REF}.tsv",
        NOVEL_PASSING_OUT =     results_prefix / "snp_calling" / "novelPassingPositionsForBLATCheck.txt"
    conda:
        "../envs/get_sample_snps.yaml"
    log:
        results_prefix / "snp_calling" / "logs" / "{YSEQID}_{REF}_confirm_novel_snps.log"
    benchmark:
        results_prefix / "snp_calling" / "benchmark" / "{YSEQID}_{REF}_confirm_novel_snps.benchmark"
    params:
        TEMP = tmp_prefix
    threads:
        workflow.cores  
    shell:
        """
        mkdir -p {params}
        (python3 workflow/scripts/identityResolutionTemplateCreator.py -batch {input.NOVEL_VCF_TSV} {output.NOVEL_TSV} {input.REFSEQ} {output.NOVEL_PASSING_OUT} {threads} {params}) > {log} 2>&1
        """
rule indel_calling:
    input:
        CALLED_VCF =    results_prefix  / "snp_calling" / "chrY_called_{YSEQID}_{REF}.vcf.gz"
    output:
        INDEL_VCF =     results_prefix  / "snp_calling" / "chrY_INDELs_{YSEQID}_{REF}.gz"
        INDEL_TBI =     results_prefix  / "snp_calling" / "chrY_INDELs_{YSEQID}_{REF}.gz.tbi"
    conda:
        "../envs/get_sample_snps.yaml"
    log:
        results_prefix / "snp_calling" / "logs" / "{YSEQID}_{REF}_bcf_filter_indels.log"
    benchmark:
        results_prefix / "snp_calling" / "benchmark" / "{YSEQID}_{REF}_bcf_filter_indels.benchmark"

    shell:
        """
        (bcftools filter -i 'TYPE="indel" && QUAL>=30 && GT=="1/1" && DP4[2] + DP4[3] >=2 && (DP4[0] + DP4[1]) / DP < 0.1' {input.CALLED_VCF} | bcftools view -n -O z -s ^HARRY,ALIEN > {output.INDEL_VCF})  > {log} 2>&1
        tabix {output.INDEL_VCF} >> {log} 2>&1
        """