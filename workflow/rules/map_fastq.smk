rule index_refseq:
    input:
        REFSEQ = "resources/refseq/{REF}/{REF}.fa"
    output:
        BWT = "resources/refseq/{REF}/{REF}.fa.bwt",
        AMB = "resources/refseq/{REF}/{REF}.fa.amb",
        ANN = "resources/refseq/{REF}/{REF}.fa.ann",
        PAC = "resources/refseq/{REF}/{REF}.fa.pac",
        SA  = "resources/refseq/{REF}/{REF}.fa.sa"

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
        BWT = "resources/refseq/{REF}/{REF}.fa.bwt",
        AMB = "resources/refseq/{REF}/{REF}.fa.amb",
        ANN = "resources/refseq/{REF}/{REF}.fa.ann",
        PAC = "resources/refseq/{REF}/{REF}.fa.pac",
        SA  = "resources/refseq/{REF}/{REF}.fa.sa"
        #"./resources/refseq/{REF}/{REF}.fa"
    output: 
        BAM = "results/{YSEQID}_bwa-mem_{REF}.bam"
        
    params:
        BWA = "-M -t " #-t has to be last!
    threads: workflow.cores
    shell:
        """
        bwa mem {params.BWA} {threads} {input.REFSEQ} {input.READS_1} {input.READS_2} | 
        samtools view -@ {threads} -b -t {input.REFSEQ} -o {output.BAM} - 
        """