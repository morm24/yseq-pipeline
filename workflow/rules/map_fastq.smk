rule index_refseq_minimap2: 
    input:
        REFSEQ = ref_prefix / "{REF}/{REF}.fa"
    output:
        INDEX = ref_prefix / "{REF}/{REF}.fa.mmi"
    log:
        ref_prefix / "{REF}/{REF}_minimap2_index.log"
    benchmark:
        ref_prefix / "{REF}/{REF}_minimap2_index.benchmark"
    conda:
        "mapping"
    shell:
        """
        minimap2 {input.REFSEQ} -d {output.INDEX} > {log} 2>&1
        """


rule map_minimap2:
    input:
        READS_1 = sample_prefix / "{YSEQID}_R1.fastq.gz",
        READS_2 = sample_prefix / "{YSEQID}_R2.fastq.gz",
        #
        INDEX  = ref_prefix / "{REF}/{REF}.fa.mmi",
        REFSEQ = ref_prefix / "{REF}/{REF}.fa"
        #"./resources/refseq/{REF}/{REF}.fa"
    output: 
        BAM = temp(results_prefix / "mapping" / "{YSEQID}_minimap2_{REF}.bam")
    conda:
        "mapping"
    log:
        results_prefix / "mapping" / "logs" / "{YSEQID}_minimap2_{REF}.log"
    benchmark:
        results_prefix / "mapping" / "benchmark" / "{YSEQID}_minimap2_{REF}.benchmark"
    params:
        "-a -x sr -t " if (config["READS"] == "short") else "-ax map-ont -t" #if (config["READS"] == "nanopore") "-ax map-ont -t" else "-at"#-t has to be last!
    threads: workflow.cores
    shell:
        """
        (minimap2 {params} {threads} {input.INDEX} {input.READS_1} {input.READS_2} | 
        samtools view -@ {threads} -b -t {input.REFSEQ} -o {output.BAM} - ) > {log} 2>&1
        """


rule index_refseq_bwa:
    input:
        REFSEQ = ref_prefix / "{REF}/{REF}.fa"
    output:
        BWT = ref_prefix / "{REF}/{REF}.fa.bwt",
        AMB = ref_prefix / "{REF}/{REF}.fa.amb",
        ANN = ref_prefix / "{REF}/{REF}.fa.ann",
        PAC = ref_prefix / "{REF}/{REF}.fa.pac",
        SA  = ref_prefix / "{REF}/{REF}.fa.sa"
    log:
        ref_prefix / "{REF}/{REF}_bwa_index.log"
    benchmark:
        ref_prefix / "{REF}/{REF}_bwa_index.benchmark"
    conda:
        "mapping"
    shell:
        """
        bwa index {input.REFSEQ} > {log} 2>&1
        """

rule map_bwa:
    input:
        READS_1 = sample_prefix / "{YSEQID}_R1.fastq.gz",
        READS_2 = sample_prefix / "{YSEQID}_R2.fastq.gz",
        #
        REFSEQ  = ref_prefix / "{REF}/{REF}.fa",
        BWT = ref_prefix / "{REF}/{REF}.fa.bwt",
        AMB = ref_prefix / "{REF}/{REF}.fa.amb",
        ANN = ref_prefix / "{REF}/{REF}.fa.ann",
        PAC = ref_prefix / "{REF}/{REF}.fa.pac",
        SA  = ref_prefix / "{REF}/{REF}.fa.sa"
        #"./resources/refseq/{REF}/{REF}.fa"
    output: 
        BAM = temp(results_prefix / "mapping" / "{YSEQID}_bwa-mem_{REF}.bam")
    conda:
        "mapping"
    log:
        results_prefix / "mapping" / "logs" / "{YSEQID}_bwa-mem_{REF}.log"
    benchmark:
        results_prefix / "mapping" / "benchmark" / "{YSEQID}_bwa-mem_{REF}.benchmark"
    params:
        BWA = "-M -t " #-t has to be last!
    threads: workflow.cores
    shell:
        """
        (bwa mem {params.BWA} {threads} {input.REFSEQ} {input.READS_1} {input.READS_2} | 
        samtools view -@ {threads} -b -t {input.REFSEQ} -o {output.BAM} - ) > {log} 2>&1
        """

rule sort_and_index:
    input:
        BAM = results_prefix / "mapping" / "{YSEQID}_bwa-mem_{REF}.bam" if (config["MAPPER"] == "bwa") else results_prefix / "{YSEQID}_minimap2_{REF}.bam" #elif (config["MAPPER"] == "minimap2") results_prefix / "{YSEQID}_minimap2_{REF}.bam"
    output:
        SORTED_BAM = results_prefix / "mapping" / "{YSEQID}_bwa-mem_{REF}_sorted.bam",
        BAI = results_prefix / "mapping" / "{YSEQID}_bwa-mem_{REF}_sorted.bam.bai",
        IDXSTATS = results_prefix / "mapping" / "{YSEQID}_bwa-mem_{REF}_sorted.bam.idxstats.tsv"
    conda:
        "mapping"
    log:
        results_prefix / "mapping" / "logs" / "{YSEQID}_{REF}_samtools_sort.log"
    benchmark:
        results_prefix / "mapping" / "benchmark" / "{YSEQID}_{REF}_samtools_sort.benchmark"
    params:
        #how do I pass the snakemake num threads to the tool?
        #NUM_THREADS = "-@ 4",
        #SORTDIR = "-T /usr/local/geospiza/var/tmp/"
        SORTDIR = "-T resources/tmp/"
        #OUTPUT =  "-o .SORTED_BAM}"
    threads: workflow.cores * 1
    shell:
        """
        samtools sort -@ {threads} {params.SORTDIR}sorted -o {output.SORTED_BAM} {input.BAM}    > {log} 2>&1
	    samtools index -@ {threads} {output.SORTED_BAM}                                         >> {log} 2>&1
	    samtools idxstats {output.SORTED_BAM} | tee {output.IDXSTATS}                           >> {log} 2>&1
        """

#swap the code with bamstaistics or sth. like that
rule get_mapping_statistics:
    input:
        BAM = results_prefix / "mapping" / "{YSEQID}_bwa-mem_{REF}_sorted.bam",
        IDXSTATS =  results_prefix / "mapping" / "{YSEQID}_bwa-mem_{REF}_sorted.bam.idxstats.tsv"

    output:
        STATS = results_prefix / "mapping" / "{YSEQID}_{REF}_mapping_stats.txt"
    conda:
        "mapping"
    log:
        results_prefix / "mapping" / "logs" / "{YSEQID}_{REF}_bamstats.log"
    benchmark:
        results_prefix / "mapping" / "benchmark" / "{YSEQID}_{REF}_bamstats.benchmark"
    threads: 
        workflow.cores * 1
    shell:
        """
        bamstats -u -i {input.BAM} -o {output.STATS} > {log} 2>&1
        """
        
#chromosome = ["chrY", "chrM"]
rule seperate_BAM:
    input:
        BAM = results_prefix / "mapping" / "{YSEQID}_bwa-mem_{REF}_sorted.bam"
        
    output:
        #CHR_BAM = expand(results_prefix / "mapping" / "{{YSEQID}}_bwa-mem_{{REF}}_{chr}.bam", chr=chromosome)
        CHR_BAM = results_prefix / "mapping" / "{YSEQID}_bwa-mem_{REF}_{chr}.bam"

    conda:
        "mapping"
    log:
        results_prefix / "mapping" / "logs" / "{YSEQID}_bwa-mem_{REF}_{chr}.log"
    benchmark:
        results_prefix / "mapping" / "benchmark" / "{YSEQID}_bwa-mem_{REF}_{chr}.benchmark"
    threads: 
        workflow.cores * 1
    shell:
        """
        (samtools view -@ {threads} -b -o {output.CHR_BAM} {input.BAM} {wildcards.chr}) > {log} 2>&1
        samtools index {output.CHR_BAM} >> {log} 2>&1
        """
