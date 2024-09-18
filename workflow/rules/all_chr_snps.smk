for chr in ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]:
    rule:
        name: 
            f"call_all_snps_{chr}"
        input: 
            SORTED_BAM =    results_prefix / "{YSEQID}_bwa_mem_{REF}_sorted.bam",
            REFSEQ =        "resources/refseq/{REF}/{REF}.fa"
        output: 
            VCF =           temp(f"resources/tmp/{chr}_{{YSEQID}}_{{REF}}.vcf.gz")
            threads: 4
        shell: f"bcftools mpileup -r {chr} -Ou -C 0 -f {{input.REFSEQ}} {{input.SORTED_BAM}} | bcftools call -O z --threads {{threads}} -v -V indels -m -P 0 > {{output.VCF}}"

rule combine_chr_snps:
    input:
        CHR1_VCF =     "resources/tmp/chr1_{YSEQID}_{REF}.vcf.gz",
        CHR2_VCF =     "resources/tmp/chr2_{YSEQID}_{REF}.vcf.gz",
        CHR3_VCF =     "resources/tmp/chr3_{YSEQID}_{REF}.vcf.gz",
        CHR4_VCF =     "resources/tmp/chr4_{YSEQID}_{REF}.vcf.gz",
        CHR5_VCF =     "resources/tmp/chr5_{YSEQID}_{REF}.vcf.gz",
        CHR6_VCF =     "resources/tmp/chr6_{YSEQID}_{REF}.vcf.gz",
        CHR7_VCF =     "resources/tmp/chr7_{YSEQID}_{REF}.vcf.gz",
        CHR8_VCF =     "resources/tmp/chr8_{YSEQID}_{REF}.vcf.gz",
        CHR9_VCF =     "resources/tmp/chr9_{YSEQID}_{REF}.vcf.gz",
        CHR10_VCF =     "resources/tmp/chr10_{YSEQID}_{REF}.vcf.gz",
        CHR11_VCF =     "resources/tmp/chr11_{YSEQID}_{REF}.vcf.gz",
        CHR12_VCF =     "resources/tmp/chr12_{YSEQID}_{REF}.vcf.gz",
        CHR13_VCF =     "resources/tmp/chr13_{YSEQID}_{REF}.vcf.gz",
        CHR14_VCF =     "resources/tmp/chr14_{YSEQID}_{REF}.vcf.gz",
        CHR15_VCF =     "resources/tmp/chr15_{YSEQID}_{REF}.vcf.gz",
        CHR16_VCF =     "resources/tmp/chr16_{YSEQID}_{REF}.vcf.gz",
        CHR17_VCF =     "resources/tmp/chr17_{YSEQID}_{REF}.vcf.gz",
        CHR18_VCF =     "resources/tmp/chr18_{YSEQID}_{REF}.vcf.gz",
        CHR19_VCF =     "resources/tmp/chr19_{YSEQID}_{REF}.vcf.gz",
        CHR20_VCF =     "resources/tmp/chr20_{YSEQID}_{REF}.vcf.gz",
        CHR21_VCF =     "resources/tmp/chr21_{YSEQID}_{REF}.vcf.gz",
        CHR22_VCF =     "resources/tmp/chr22_{YSEQID}_{REF}.vcf.gz",
        CHRX_VCF =     "resources/tmp/chrX_{YSEQID}_{REF}.vcf.gz",
        CHRY_VCF =     "resources/tmp/chrY_{YSEQID}_{REF}.vcf.gz"
    output:
        ALL_CHR_SNPS = results_prefix / "all_chr_snps_{YSEQID}_{REF}.vcf.gz"
    shell:
        """
        bcftools concat -O z resources/tmp/chr[1-9]_{wildcards.YSEQID}_{wildcards.REF}.vcf.gz resources/tmp/chr[1-2][0-9]_{wildcards.YSEQID}_{wildcards.REF}.vcf.gz resources/tmp/chr[M,X-Y]_{wildcards.YSEQID}_{wildcards.REF}.vcf.gz > {output.ALL_CHR_SNPS}
	    tabix {output.ALL_CHR_SNPS}
        """