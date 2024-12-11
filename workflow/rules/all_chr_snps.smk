rule call_chr_snp:
    input: 
        SORTED_BAM =    results_prefix / "mapping" / "{YSEQID}_bwa-mem_{REF}_sorted.bam",
        REFSEQ =        ref_prefix / "{REF}/{REF}.fa"
    output: 
        VCF =           temp("resources/tmp/{chr}_{YSEQID}_{REF}.vcf.gz")
    conda:
        "bam_process"
    threads: 4
    log: results_prefix / "call_chr_snps" /  "log" / "call_chr_snp_{chr}_{YSEQID}_{REF}.log"
    shell: "(bcftools mpileup -r {wildcards.chr} -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.VCF}) > {log} 2>&1" 



chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
rule combine_chr_snps:
    input:
        expand("resources/tmp/{chr}_{{YSEQID}}_{{REF}}.vcf.gz", chr=chromosomes)
    output:
        ALL_CHR_SNPS = results_prefix / "call_chr_snps" / "all_chr_snps_{YSEQID}_{REF}.vcf.gz"
    conda:
        "bam_process" 
    log: results_prefix / "call_chr_snps" / "log" / "combine_chr_snps_{YSEQID}_{REF}.log"
    benchmark:
        results_prefix / "call_chr_snps" / "benchmark" / "combine_chr_snps_{YSEQID}_{REF}.benchmark"
    shell:
        """
        bcftools concat -O z {input} -o {output.ALL_CHR_SNPS}  > {log} 2>&1
	    tabix {output.ALL_CHR_SNPS} > {log} 2>&1
        """