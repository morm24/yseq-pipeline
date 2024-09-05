rule map_bwa_hg38:
    input:
        READS_1 = "resources/sample/YSEQID_R1.fastq.gz",
        READS_2 = "resources/sample/YSEQID_R2.fastq.gz",
        #
        REF     = "./resources/refseq/hg38/hg38.fa"
        #"./resources/refseq/{REF}/{REF}.fa"
    output: 
        BAMFILE = "YSEQID_bwa-mem.bam"
    params:
        BWA = "-M -t 4"
        #"-M -t $NUM_THREADS"
    shell:
        """
        bwa mem {params.BWA} {input.REF} {input.READS_1} {input.READS_2} | 
        samtools view -@ 4 -b -t {input.REF} -o {output.BAMFILE} - 
        """
