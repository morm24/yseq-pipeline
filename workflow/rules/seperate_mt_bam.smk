rule seperate_mtBAM:
    input:
        SORTED_BAM =    "results/{YSEQID}_bwa_mem_{REF}_sorted.bam"    
    output:
        MTBAM = "results/{YSEQID}_bwa-mem_{REF}_rCRS_chrM.bam"
    threads: 
        workflow.cores * 1
    shell:
        """
        samtools view -@ {threads} -b -o {output.MTBAM} {input.SORTED_BAM} chrM
        samtools index {output.MTBAM}
        """