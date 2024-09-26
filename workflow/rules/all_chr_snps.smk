chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
for chr in chromosomes:
    rule:
        name: 
            f"call_all_snps_{chr}"
        input: 
            SORTED_BAM =    results_prefix / "{YSEQID}_bwa_mem_{REF}_sorted.bam",
            REFSEQ =        ref_prefix / "{REF}/{REF}.fa"
        output: 
            VCF =           temp(f"resources/tmp/{chr}_{{YSEQID}}_{{REF}}.vcf.gz")
        conda:
            env_path / "bam_process.yaml"
        threads: 4
        shell: f"bcftools mpileup -r {chr} -Ou -C 0 -f {{input.REFSEQ}} {{input.SORTED_BAM}} | bcftools call -O z --threads {{threads}} -v -V indels -m -P 0 > {{output.VCF}}"


# Define a dictionary to map chromosome names to their respective VCF file paths
#vcf_files = {f"{chr}_VCF": f"resources/tmp/{chr}_{{YSEQID}}_{{REF}}.vcf.gz" for chr in chromosomes}

rule combine_chr_snps:
    input:
        expand("resources/tmp/{chr}_{YSEQID}_{REF}.vcf.gz", chr=chromosomes)
    output:
        ALL_CHR_SNPS = results_prefix / "all_chr_snps_{YSEQID}_{REF}.vcf.gz"
    conda:
        env_path / "bam_process.yaml"
    shell:
        """
        bcftools concat -O z {input} > {output.ALL_CHR_SNPS}
	    tabix {output.ALL_CHR_SNPS}
        """