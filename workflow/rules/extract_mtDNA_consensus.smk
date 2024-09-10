rule mt_consensus:
    input:
        BAM_SORTED = "results/{YSEQID}_bwa_mem_{REF}_sorted.bam",
        REFSEQ = "resources/refseq/{REF}/{REF}.fa"

    output:
        VCF = "results/chrM_{YSEQID}_{REF}.vcf.gz",
        CONSENSUS = "results/{YSEQID}_{REF}_mtDNA.fasta"
    threads:
        workflow.cores * 1

    shell: 
        """
        bcftools mpileup -r chrM -Ou -C 50 -f {input.REFSEQ} {input.BAM_SORTED} | bcftools call --threads {threads} -O z -v -m -P 0  > {output.VCF}
        tabix {output.VCF}
        samtools faidx {input.REFSEQ} chrM | bcftools consensus {output.VCF} -o {output.CONSENSUS}
        """