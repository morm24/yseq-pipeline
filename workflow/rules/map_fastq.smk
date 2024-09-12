rule index_refseq:
    input:
        REFSEQ = "resources/refseq/{REF}/{REF}.fa"
    output:
        INDEX = "resources/refseq/{REF}/{REF}.fa.fai"
    shell:
        """
        bwa index {input.REFSEQ}
        """
rule map_bwa:
    input:
        READS_1 = "resources/sample/{YSEQID}_R1.fastq.gz",
        READS_2 = "resources/sample/{YSEQID}_R2.fastq.gz",
        #
        REFSEQ  = "resources/refseq/{REF}/{REF}.fa",
        INDEX   = "resources/refseq/{REF}/{REF}.fa.fai"
        #"./resources/refseq/{REF}/{REF}.fa"
    output: 
        BAM = "results/{YSEQID}_bwa-mem_{REF}.bam"
        
    params:
        BWA = "-M -t "
    threads: workflow.cores
    shell:
        """
        bwa mem {params.BWA} {threads} {input.REFSEQ} {input.READS_1} {input.READS_2} | 
        samtools view -@ {threads} -b -t {input.REFSEQ} -o {output.BAM} - 
        """