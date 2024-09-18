rule mt_consensus:
    input:
        BAM_SORTED = results_prefix / "{YSEQID}_bwa_mem_{REF}_sorted.bam",
        REFSEQ = "resources/refseq/{REF}/{REF}.fa"

    output:
        VCF = results_prefix / "chrM_{YSEQID}_{REF}.vcf.gz",
        CONSENSUS = results_prefix / "{YSEQID}_{REF}_mtDNA.fasta"

    conda:
        env_path / "bam_process.yaml"
    threads:
        workflow.cores * 1

    shell: 
        """
        bcftools mpileup -r chrM -Ou -C 50 -f {input.REFSEQ} {input.BAM_SORTED} | bcftools call --threads {threads} -O z -v -m -P 0  > {output.VCF}
        tabix {output.VCF}
        samtools faidx {input.REFSEQ} chrM | bcftools consensus {output.VCF} -o {output.CONSENSUS}
        """

rule get_mtdna_differences_process:
    input:
        VCF =           results_prefix / "chrM_{YSEQID}_{REF}.vcf.gz",
        #HAPLOGREP =     "haplogrep/haplogrep-2.1.25.jar"
    output:
        MTDNA_SNPS =    results_prefix / "{YSEQID}_{REF}_MTDNA_SNPS.tsv"
    shell:
        """
        python3 workflow/scripts/script_templates/getMTDNADifferences.py -process {input.VCF} haplogrep/haplogrep-2.1.25.jar {output.MTDNA_SNPS}
        
        """ 
rule get_mtdna_differences_update:
    input:
        MT_VCF =           results_prefix / "chrM_{YSEQID}_{REF}.vcf.gz",
        #HAPLOGREP =     "haplogrep/haplogrep-2.1.25.jar"
        PHYLOEQ_SNPS = results_prefix / "{YSEQID}_{REF}_PHYLOEQ_SNPS.tsv",
        DOWNSTR_SNPS = results_prefix / "{YSEQID}_{REF}_DOWNSTR_SNPS.tsv"
    output:
        TEST_OUTPUT =   results_prefix / "{YSEQID}_{REF}_result_summary.txt"
    shell:
        """
        mkdir -p {output.TEST_OUTPUT}
        python3 workflow/scripts/script_templates/getMTDNADifferences.py -update {output.TEST_OUTPUT}  {output.TEST_OUTPUT}.out {input.MT_VCF} {input.PHYLOEQ_SNPS} {input.DOWNSTR_SNPS} "The differences to the rCRS are" "phylo-equivalent SNPs to" "known downstream SNPs to" "novel SNPs found in sample"        
        """
rule get_mtDifferences_addAlleles:
    input:
        MTDNA_SNPS =    results_prefix / "{YSEQID}_{REF}_MTDNA_SNPS.tsv",
        MT_VCF =        results_prefix / "chrM_{YSEQID}_{REF}.vcf.gz",
        PHYLOEQ_SNPS = results_prefix / "{YSEQID}_{REF}_PHYLOEQ_SNPS.tsv",
        DOWNSTR_SNPS = results_prefix / "{YSEQID}_{REF}_DOWNSTR_SNPS.tsv"
    output:
        ALLELES_TSV = results_prefix / "{YSEQID}_{REF}_addAlleles.tsv"
        #ALLELES_TSV = "addAlleles.tsv"
    shell:
        """
        python3 workflow/scripts/script_templates/getMTDNADifferences.py -addAlleles {output.ALLELES_TSV} {input.MTDNA_SNPS} {wildcards.YSEQID} {input.MT_VCF} {input.PHYLOEQ_SNPS} {input.DOWNSTR_SNPS} chrY_{wildcards.YSEQID}_{wildcards.REF}.vcf.gz
        """

rule get_mtDifferences_create_update_script:
    input:
        TEST_OUTPUT =   results_prefix / "{YSEQID}_{REF}_result_summary.txt",
        #ALLELES_TSV =   results_prefix / "{YSEQID}_{REF}_addAlleles.tsv",
        MTDNA_SNPS =    results_prefix / "{YSEQID}_{REF}_MTDNA_SNPS.tsv",
        MT_VCF =        results_prefix / "chrM_{YSEQID}_{REF}.vcf.gz",
        CHRY_VCF =      results_prefix / "chrY_{YSEQID}_{REF}.vcf.gz",
        PHYLOEQ_SNPS =  results_prefix / "{YSEQID}_{REF}_PHYLOEQ_SNPS.tsv",
        DOWNSTR_SNPS =  results_prefix / "{YSEQID}_{REF}_DOWNSTR_SNPS.tsv"
    output:
        UPDATE_SCRIPT = results_prefix / "{YSEQID}_{REF}_update.sh"
    shell:
        """
        python3 workflow/scripts/script_templates/getMTDNADifferences.py -createUpdateScript \
        addAlleles.tsv {input.TEST_OUTPUT} {input.TEST_OUTPUT}.out {wildcards.YSEQID} {input.MT_VCF} \
        {input.CHRY_VCF} {input.MTDNA_SNPS} {input.PHYLOEQ_SNPS} {input.DOWNSTR_SNPS} \
        "The differences to the rCRS are" "phylo-equivalent SNPs to" "known downstream SNPs to" "novel SNPs found in sample" {output.UPDATE_SCRIPT}
        chmod a+x updateScript.sh
        """
    
rule phenotyping:
    input:

    output:

    shell:
        """
        ./workflow/scripts/script_templates/phenotyping_hg38_pipeline.sh
        """