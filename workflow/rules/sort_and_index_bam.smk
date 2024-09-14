rule sort_and_index:
    input:
        BAM = "results/{YSEQID}_bwa-mem_{REF}.bam"
    output:
        SORTED_BAM = "results/{YSEQID}_bwa_mem_{REF}_sorted.bam",
        BAI = "results/{YSEQID}_bwa_mem_{REF}_sorted.bam.bai",
        IDXSTATS = "results/{YSEQID}_bwa_mem_{REF}_sorted.bam.idxstats.tsv"
    params:
        #how do I pass the snakemake num threads to the tool?
        #NUM_THREADS = "-@ 4",
        #SORTDIR = "-T /usr/local/geospiza/var/tmp/"
        SORTDIR = "-T resources/tmp/"
        #OUTPUT =  "-o .SORTED_BAM}"
    threads: workflow.cores * 1
    shell:
        """
        samtools sort -@ {threads} {params.SORTDIR}sorted -o {output.SORTED_BAM} {input.BAM}
	    samtools index -@ {threads} {output.SORTED_BAM}
	    samtools idxstats {output.SORTED_BAM} > {output.IDXSTATS}
        """