rule seperate_yBAM:
    input:
        BAM = results_prefix / "{YSEQID}_bwa_mem_{REF}_sorted.bam"
    output:
        YBAM = results_prefix / "{YSEQID}_bwa-mem_{REF}_chrY.bam"
    conda:
        env_path / "bam_process.yaml"
    threads: 
        workflow.cores * 1
    shell:
        """
        samtools view -@ {threads} -b -o {output.YBAM} {input.BAM} chrY
        samtools index {output.YBAM}
        """
rule seperate_mtBAM:
    input:
        SORTED_BAM =    results_prefix / "{YSEQID}_bwa_mem_{REF}_sorted.bam"    
    output:
        MTBAM = results_prefix / "{YSEQID}_bwa-mem_{REF}_rCRS_chrM.bam"
    conda:
        env_path / "bam_process.yaml"
    threads: 
        workflow.cores * 1
    shell:
        """
        samtools view -@ {threads} -b -o {output.MTBAM} {input.SORTED_BAM} chrM
        samtools index {output.MTBAM}
        """