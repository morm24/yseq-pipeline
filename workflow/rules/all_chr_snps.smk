rule get_all_chr_snps:
    input:
        SORTED_BAM =    "results/{YSEQID}_bwa_mem_{REF}_sorted.bam",
        REFSEQ =        "resources/refseq/{REF}/{REF}.fa"
    output:
        CHR1_VCF =     temp("resources/tmp/chr1_{YSEQID}_{REF}.gz"),
        CHR2_VCF =     temp("resources/tmp/chr2_{YSEQID}_{REF}.gz"),
        CHR3_VCF =     temp("resources/tmp/chr3_{YSEQID}_{REF}.gz"),
        CHR4_VCF =     temp("resources/tmp/chr4_{YSEQID}_{REF}.gz"),
        CHR5_VCF =     temp("resources/tmp/chr5_{YSEQID}_{REF}.gz"),
        CHR6_VCF =     temp("resources/tmp/chr6_{YSEQID}_{REF}.gz"),
        CHR7_VCF =     temp("resources/tmp/chr7_{YSEQID}_{REF}.gz"),
        CHR8_VCF =     temp("resources/tmp/chr8_{YSEQID}_{REF}.gz"),
        CHR9_VCF =     temp("resources/tmp/chr9_{YSEQID}_{REF}.gz"),
        CHR10_VCF =     temp("resources/tmp/chr10_{YSEQID}_{REF}.gz"),
        CHR11_VCF =     temp("resources/tmp/chr11_{YSEQID}_{REF}.gz"),
        CHR12_VCF =     temp("resources/tmp/chr12_{YSEQID}_{REF}.gz"),
        CHR13_VCF =     temp("resources/tmp/chr13_{YSEQID}_{REF}.gz"),
        CHR14_VCF =     temp("resources/tmp/chr14_{YSEQID}_{REF}.gz"),
        CHR15_VCF =     temp("resources/tmp/chr15_{YSEQID}_{REF}.gz"),
        CHR16_VCF =     temp("resources/tmp/chr16_{YSEQID}_{REF}.gz"),
        CHR17_VCF =     temp("resources/tmp/chr17_{YSEQID}_{REF}.gz"),
        CHR18_VCF =     temp("resources/tmp/chr18_{YSEQID}_{REF}.gz"),
        CHR19_VCF =     temp("resources/tmp/chr19_{YSEQID}_{REF}.gz"),
        CHR20_VCF =     temp("resources/tmp/chr20_{YSEQID}_{REF}.gz"),
        CHR21_VCF =     temp("resources/tmp/chr21_{YSEQID}_{REF}.gz"),
        CHR22_VCF =     temp("resources/tmp/chr22_{YSEQID}_{REF}.gz"),
        CHRX_VCF =     temp("resources/tmp/chrX_{YSEQID}_{REF}.gz"),
        CHRY_VCF =     temp("resources/tmp/chrY_{YSEQID}_{REF}.gz")
    shell:
        """
	    bcftools mpileup -r chr1 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR1_VCF} &
	    bcftools mpileup -r chr2 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR2_VCF} &
	    bcftools mpileup -r chr3 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR3_VCF} &
	    bcftools mpileup -r chr4 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR4_VCF} &
	    bcftools mpileup -r chr5 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR5_VCF} &
	    bcftools mpileup -r chr6 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR6_VCF} &
	    bcftools mpileup -r chr7 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR7_VCF} &
	    bcftools mpileup -r chr8 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR8_VCF} &
	    bcftools mpileup -r chr9 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR9_VCF} &
	    bcftools mpileup -r chr10 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR10_VCF} &
	    bcftools mpileup -r chr11 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR11_VCF} &
	    bcftools mpileup -r chr12 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR12_VCF} &
	    bcftools mpileup -r chr13 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR13_VCF} &
	    bcftools mpileup -r chr14 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR14_VCF} &
	    bcftools mpileup -r chr15 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR15_VCF} &
	    bcftools mpileup -r chr16 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR16_VCF} &
	    bcftools mpileup -r chr17 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR17_VCF} &
	    bcftools mpileup -r chr18 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR18_VCF} &
	    bcftools mpileup -r chr19 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR19_VCF} &
	    bcftools mpileup -r chr20 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR20_VCF} &
	    bcftools mpileup -r chr21 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR21_VCF} &
	    bcftools mpileup -r chr22 -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > {output.CHR22_VCF} &
	    bcftools mpileup -r chrX -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0  > {output.CHRX_VCF} &
	    bcftools mpileup -r chrY -Ou -C 0 -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0  > {output.CHRY_VCF} &  ##same as in getMTDNADifferences.py addAlleles
	    wait
        """
rule combine_chr_snps:
    input:
        CHR1_VCF =     "resources/tmp/chr1_{YSEQID}_{REF}.gz",
        CHR2_VCF =     "resources/tmp/chr2_{YSEQID}_{REF}.gz",
        CHR3_VCF =     "resources/tmp/chr3_{YSEQID}_{REF}.gz",
        CHR4_VCF =     "resources/tmp/chr4_{YSEQID}_{REF}.gz",
        CHR5_VCF =     "resources/tmp/chr5_{YSEQID}_{REF}.gz",
        CHR6_VCF =     "resources/tmp/chr6_{YSEQID}_{REF}.gz",
        CHR7_VCF =     "resources/tmp/chr7_{YSEQID}_{REF}.gz",
        CHR8_VCF =     "resources/tmp/chr8_{YSEQID}_{REF}.gz",
        CHR9_VCF =     "resources/tmp/chr9_{YSEQID}_{REF}.gz",
        CHR10_VCF =     "resources/tmp/chr10_{YSEQID}_{REF}.gz",
        CHR11_VCF =     "resources/tmp/chr11_{YSEQID}_{REF}.gz",
        CHR12_VCF =     "resources/tmp/chr12_{YSEQID}_{REF}.gz",
        CHR13_VCF =     "resources/tmp/chr13_{YSEQID}_{REF}.gz",
        CHR14_VCF =     "resources/tmp/chr14_{YSEQID}_{REF}.gz",
        CHR15_VCF =     "resources/tmp/chr15_{YSEQID}_{REF}.gz",
        CHR16_VCF =     "resources/tmp/chr16_{YSEQID}_{REF}.gz",
        CHR17_VCF =     "resources/tmp/chr17_{YSEQID}_{REF}.gz",
        CHR18_VCF =     "resources/tmp/chr18_{YSEQID}_{REF}.gz",
        CHR19_VCF =     "resources/tmp/chr19_{YSEQID}_{REF}.gz",
        CHR20_VCF =     "resources/tmp/chr20_{YSEQID}_{REF}.gz",
        CHR21_VCF =     "resources/tmp/chr21_{YSEQID}_{REF}.gz",
        CHR22_VCF =     "resources/tmp/chr22_{YSEQID}_{REF}.gz",
        CHRX_VCF =     "resources/tmp/chrX_{YSEQID}_{REF}.gz",
        CHRY_VCF =     "resources/tmp/chrY_{YSEQID}_{REF}.gz"
    output:
        ALL_CHR_SNPS = "results/all_chr_snps_{YSEQID}_{REF}.vcf.gz"
    shell:
        """
        bcftools concat -O z resources/tmp/chr[1-9]_{wildcards.YSEQID}_{wildcards.REF}.gz resources/tmp/chr[1-2][0-9]_{wildcards.YSEQID}_{wildcards.REF}.gz resources/tmp/chr[M,X-Y]_{wildcards.YSEQID}_{wildcards.REF}.gz > {output.ALL_CHR_SNPS}
	    tabix {output.ALL_CHR_SNPS}
        """