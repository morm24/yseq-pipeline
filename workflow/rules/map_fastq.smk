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
#swap the code with bamstaistics or sth. like that
rule get_mapping_statistics:
    input:
        BAM = "results/{YSEQID}_bwa_mem_{REF}_sorted.bam",
        IDXSTATS =  "results/{YSEQID}_bwa_mem_{REF}_sorted.bam.idxstats.tsv"

    output:
        STATS = "results/{YSEQID}_{REF}_mapping_stats.txt"
    threads: 
        workflow.cores * 1
    shell:
        """
        bamstats -u -i {input.BAM} -o {output.STATS}
        """