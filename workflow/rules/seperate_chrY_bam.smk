rule seperate_yBAM:
    input:
        BAM = "results/{YSEQID}_bwa_mem_{REF}_sorted.bam"
    output:
        YBAM = "results/{YSEQID}_bwa-mem_{REF}_chrY.bam"
    threads: 
        workflow.cores * 1
    shell:
        """
        samtools view -@ {threads} -b -o {output.YBAM} {input.BAM} chrY
        samtools index {output.YBAM}
        """