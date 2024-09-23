rule call_all_chr_snps:
    input:
        SORTED_BAM =    results_prefix / "{YSEQID}_bwa_mem_{REF}_sorted.bam",
        REFSEQ =        "resources/refseq/{REF}/{REF}.fa"
    output:
        CHR1_VCF =     temp(tmp_prefix / "chr1_{YSEQID}_{REF}.vcf.gz"),
        CHR2_VCF =     temp(tmp_prefix / "chr2_{YSEQID}_{REF}.vcf.gz"),
        CHR3_VCF =     temp(tmp_prefix / "chr3_{YSEQID}_{REF}.vcf.gz"),
        CHR4_VCF =     temp(tmp_prefix / "chr4_{YSEQID}_{REF}.vcf.gz"),
        CHR5_VCF =     temp(tmp_prefix / "chr5_{YSEQID}_{REF}.vcf.gz"),
        CHR6_VCF =     temp(tmp_prefix / "chr6_{YSEQID}_{REF}.vcf.gz"),
        CHR7_VCF =     temp(tmp_prefix / "chr7_{YSEQID}_{REF}.vcf.gz"),
        CHR8_VCF =     temp(tmp_prefix / "chr8_{YSEQID}_{REF}.vcf.gz"),
        CHR9_VCF =     temp(tmp_prefix / "chr9_{YSEQID}_{REF}.vcf.gz"),
        CHR10_VCF =     temp(tmp_prefix / "chr10_{YSEQID}_{REF}.vcf.gz"),
        CHR11_VCF =     temp(tmp_prefix / "chr11_{YSEQID}_{REF}.vcf.gz"),
        CHR12_VCF =     temp(tmp_prefix / "chr12_{YSEQID}_{REF}.vcf.gz"),
        CHR13_VCF =     temp(tmp_prefix / "chr13_{YSEQID}_{REF}.vcf.gz"),
        CHR14_VCF =     temp(tmp_prefix / "chr14_{YSEQID}_{REF}.vcf.gz"),
        CHR15_VCF =     temp(tmp_prefix / "chr15_{YSEQID}_{REF}.vcf.gz"),
        CHR16_VCF =     temp(tmp_prefix / "chr16_{YSEQID}_{REF}.vcf.gz"),
        CHR17_VCF =     temp(tmp_prefix / "chr17_{YSEQID}_{REF}.vcf.gz"),
        CHR18_VCF =     temp(tmp_prefix / "chr18_{YSEQID}_{REF}.vcf.gz"),
        CHR19_VCF =     temp(tmp_prefix / "chr19_{YSEQID}_{REF}.vcf.gz"),
        CHR20_VCF =     temp(tmp_prefix / "chr20_{YSEQID}_{REF}.vcf.gz"),
        CHR21_VCF =     temp(tmp_prefix / "chr21_{YSEQID}_{REF}.vcf.gz"),
        CHR22_VCF =     temp(tmp_prefix / "chr22_{YSEQID}_{REF}.vcf.gz"),
        CHRX_VCF =     temp(tmp_prefix / "chrX_{YSEQID}_{REF}.vcf.gz"),
        CHRY_VCF =     temp(tmp_prefix / "chrY_{YSEQID}_{REF}.vcf.gz")
    conda:
        env_path / "bam_process.yaml"
    threads: 
        workflow.cores * 1
    shell:
        """
	    bcftools mpileup -r chr1 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR1_VCF} &
	    bcftools mpileup -r chr2 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR2_VCF} &
	    bcftools mpileup -r chr3 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR3_VCF} &
	    bcftools mpileup -r chr4 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR4_VCF} &
	    bcftools mpileup -r chr5 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR5_VCF} &
	    bcftools mpileup -r chr6 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR6_VCF} &
	    bcftools mpileup -r chr7 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR7_VCF} &
	    bcftools mpileup -r chr8 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR8_VCF} &
	    bcftools mpileup -r chr9 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR9_VCF} &
	    bcftools mpileup -r chr10 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR10_VCF} &
	    bcftools mpileup -r chr11 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR11_VCF} &
	    bcftools mpileup -r chr12 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR12_VCF} &
	    bcftools mpileup -r chr13 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR13_VCF} &
	    bcftools mpileup -r chr14 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR14_VCF} &
	    bcftools mpileup -r chr15 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR15_VCF} &
	    bcftools mpileup -r chr16 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR16_VCF} &
	    bcftools mpileup -r chr17 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR17_VCF} &
	    bcftools mpileup -r chr18 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR18_VCF} &
	    bcftools mpileup -r chr19 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR19_VCF} &
	    bcftools mpileup -r chr20 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR20_VCF} &
	    bcftools mpileup -r chr21 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR21_VCF} &
	    bcftools mpileup -r chr22 -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0 > {output.CHR22_VCF} &
	    bcftools mpileup -r chrX -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0  > {output.CHRX_VCF} &
	    bcftools mpileup -r chrY -Ou -C 0 -f {input.REFSEQ} {input.SORTED_BAM} | bcftools call -O z --threads {threads} -v -V indels -m -P 0  > {output.CHRY_VCF} &  ##same as in getMTDNADifferences.py addAlleles
	    wait
        """

rule combine_chr_snps:
    input:
        CHR1_VCF =     tmp_prefix / "chr1_{YSEQID}_{REF}.vcf.gz",
        CHR2_VCF =     tmp_prefix / "chr2_{YSEQID}_{REF}.vcf.gz",
        CHR3_VCF =     tmp_prefix / "chr3_{YSEQID}_{REF}.vcf.gz",
        CHR4_VCF =     tmp_prefix / "chr4_{YSEQID}_{REF}.vcf.gz",
        CHR5_VCF =     tmp_prefix / "chr5_{YSEQID}_{REF}.vcf.gz",
        CHR6_VCF =     tmp_prefix / "chr6_{YSEQID}_{REF}.vcf.gz",
        CHR7_VCF =     tmp_prefix / "chr7_{YSEQID}_{REF}.vcf.gz",
        CHR8_VCF =     tmp_prefix / "chr8_{YSEQID}_{REF}.vcf.gz",
        CHR9_VCF =     tmp_prefix / "chr9_{YSEQID}_{REF}.vcf.gz",
        CHR10_VCF =     tmp_prefix / "chr10_{YSEQID}_{REF}.vcf.gz",
        CHR11_VCF =     tmp_prefix / "chr11_{YSEQID}_{REF}.vcf.gz",
        CHR12_VCF =     tmp_prefix / "chr12_{YSEQID}_{REF}.vcf.gz",
        CHR13_VCF =     tmp_prefix / "chr13_{YSEQID}_{REF}.vcf.gz",
        CHR14_VCF =     tmp_prefix / "chr14_{YSEQID}_{REF}.vcf.gz",
        CHR15_VCF =     tmp_prefix / "chr15_{YSEQID}_{REF}.vcf.gz",
        CHR16_VCF =     tmp_prefix / "chr16_{YSEQID}_{REF}.vcf.gz",
        CHR17_VCF =     tmp_prefix / "chr17_{YSEQID}_{REF}.vcf.gz",
        CHR18_VCF =     tmp_prefix / "chr18_{YSEQID}_{REF}.vcf.gz",
        CHR19_VCF =     tmp_prefix / "chr19_{YSEQID}_{REF}.vcf.gz",
        CHR20_VCF =     tmp_prefix / "chr20_{YSEQID}_{REF}.vcf.gz",
        CHR21_VCF =     tmp_prefix / "chr21_{YSEQID}_{REF}.vcf.gz",
        CHR22_VCF =     tmp_prefix / "chr22_{YSEQID}_{REF}.vcf.gz",
        CHRX_VCF =     tmp_prefix / "chrX_{YSEQID}_{REF}.vcf.gz",
        CHRY_VCF =     tmp_prefix / "chrY_{YSEQID}_{REF}.vcf.gz"
    output:
        ALL_CHR_SNPS = results_prefix / "all_chr_snps_{YSEQID}_{REF}.vcf.gz"
    conda:
        env_path / "bam_process.yaml"
    shell:
        """
        bcftools concat -O z resources/tmp/chr[1-9]_{wildcards.YSEQID}_{wildcards.REF}.vcf.gz resources/tmp/chr[1-2][0-9]_{wildcards.YSEQID}_{wildcards.REF}.vcf.gz resources/tmp/chr[M,X-Y]_{wildcards.YSEQID}_{wildcards.REF}.vcf.gz > {output.ALL_CHR_SNPS}
	    tabix {output.ALL_CHR_SNPS}
        """