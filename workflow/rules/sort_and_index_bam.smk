rule sort_and_index:
    input:
        BAM = "YSEQID_bwa-mem.bam"
    output:
        SORTED_BAM = "YSEQID_sorted.bam",
        BAI = "YSEQID_sorted.bam.bai",
        IDXSTATS = "YSEQID_sorted.bam.idxstats.tsv"
    params:
        NUM_THREADS = "-@ 4",
        SORTDIR = "-T /usr/local/geospiza/var/tmp/'"
        #OUTPUT =  "-o .SORTED_BAM}"
    shell:
    """
        samtools sort {params.NUM_THREADS} {params.SORTDIR}sorted -o {output.SORTED_BAM} {input.BAM}
	    samtools index {params.NUM_THREADS} {output.SORTED_BAM}
	    samtools idxstats {output.SORTED_BAM} > {output.IDXSTATS}
    """