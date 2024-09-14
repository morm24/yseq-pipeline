rule get_mtdna_differences_process:
    input:
        VCF =           "results/chrM_{YSEQID}_{REF}.vcf.gz",
        #HAPLOGREP =     "haplogrep/haplogrep-2.1.25.jar"
    output:
        MTDNA_SNPS =    "results/{YSEQID}_{REF}_MTDNA_SNPS.tsv"
    shell:
        """
        python3 workflow/scripts/script_templates/getMTDNADifferences.py -process {input.VCF} haplogrep/haplogrep-2.1.25.jar {output.MTDNA_SNPS}
        
        """ 
rule get_mtdna_differences_update:
    input:
        MT_VCF =           "results/chrM_{YSEQID}_{REF}.vcf.gz",
        #HAPLOGREP =     "haplogrep/haplogrep-2.1.25.jar"
        PHYLOEQ_SNPS = "results/{YSEQID}_{REF}_PHYLOEQ_SNPS.tsv",
        DOWNSTR_SNPS = "results/{YSEQID}_{REF}_DOWNSTR_SNPS.tsv"
    output:
        TEST_OUTPUT =   "results/{YSEQID}_{REF}_result_summary.txt"
    shell:
        """
        mkdit -p {output.TEST_OUTPUT}
        python3 workflow/scripts/script_templates/getMTDNADifferences.py -update {output.TEST_OUTPUT}  {output.TEST_OUTPUT}.out {input.MT_VCF} {input.PHYLOEQ_SNPS} {input.DOWNSTR_SNPS} "The differences to the rCRS are" "phylo-equivalent SNPs to" "known downstream SNPs to" "novel SNPs found in sample"        
        """
rule get_mtDifferences_addAlleles:
    input:
        MTDNA_SNPS =    "results/{YSEQID}_{REF}_MTDNA_SNPS.tsv",
        MT_VCF =        "results/chrM_{YSEQID}_{REF}.vcf.gz",
        PHYLOEQ_SNPS = "results/{YSEQID}_{REF}_PHYLOEQ_SNPS.tsv",
        DOWNSTR_SNPS = "results/{YSEQID}_{REF}_DOWNSTR_SNPS.tsv"
    output:
        CHRY_VCF = "results/chrY_{YSEQID}_{REF}.vcf.gz",
        ALLELES_TSV = "results/{YSEQID}_{REF}_addAlleles.tsv"
    shell:
        """
        python3 workflow/scripts/script_templates/getMTDNADifferences.py -addAlleles addAlleles.tsv {input.MTDNA_SNPS} {wildcards.YSEQID} {input.MT_VCF} {input.PHYLOEQ_SNPS} {input.DOWNSTR_SNPS} {output.CHRY_VCF}
        """

rule get_mtDifferences_create_update_script:
    input:
        TEST_OUTPUT =   "results/{YSEQID}_{REF}_result_summary.txt",
        ALLELES_TSV =   "results/{YSEQID}_{REF}_addAlleles.tsv",
        MTDNA_SNPS =    "results/{YSEQID}_{REF}_MTDNA_SNPS.tsv",
        MT_VCF =        "results/chrM_{YSEQID}_{REF}.vcf.gz",
        CHRY_VCF =      "results/chrY_{YSEQID}_{REF}.vcf.gz",
        PHYLOEQ_SNPS =  "results/{YSEQID}_{REF}_PHYLOEQ_SNPS.tsv",
        DOWNSTR_SNPS =  "results/{YSEQID}_{REF}_DOWNSTR_SNPS.tsv"
    output:
        UPDATE_SCRIPT = "results/{YSEQID}_{REF}_update.sh"
    shell:
        """
        python3 workflow/scripts/script_templates/getMTDNADifferences.py -createUpdateScript \
        {input.ALLELES_TSV} {input.TEST_OUTPUT} {input.TEST_OUTPUT}.out {wildcards.YSEQID} {input.MT_VCF} \
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