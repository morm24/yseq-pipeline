rule index_refseq_minimap2:
    input:
        REFSEQ = ref_prefix / "{REF}/{REF}.fa"
    output:
        INDEX = ref_prefix / "{REF}/{REF}.fa.mmi"
    log:
        log_prefix / "{REF}/{REF}.log"
    conda:
        env_path / "mapping.yaml"
    shell:
        """
        minimap2 {input.REFSEQ} -d {output.INDEX}
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
        BAM = results_prefix / "{YSEQID}_minimap2_{REF}.bam"
    conda:
        env_path / "mapping.yaml"
    params:
        "-a -x sr -t " if (config["READS"] == "short") else "-ax map-ont -t" #if (config["READS"] == "nanopore") "-ax map-ont -t" else "-at"#-t has to be last!
    threads: workflow.cores
    shell:
        """
        minimap2 {params} {threads} {input.INDEX} {input.READS_1} {input.READS_2} | 
        samtools view -@ {threads} -b -t {input.REFSEQ} -o {output.BAM} - 
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
        log_prefix / "{REF}/{REF}.log"
    conda:
        env_path / "mapping.yaml"
    shell:
        """
        bwa index {input.REFSEQ}
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
        BAM = results_prefix / "{YSEQID}_bwa-mem_{REF}.bam"
    conda:
        env_path / "mapping.yaml"
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
        BAM = results_prefix / "{YSEQID}_bwa-mem_{REF}.bam" if (config["MAPPER"] == "bwa") else results_prefix / "{YSEQID}_minimap2_{REF}.bam" #elif (config["MAPPER"] == "minimap2") results_prefix / "{YSEQID}_minimap2_{REF}.bam"
    output:
        SORTED_BAM = results_prefix / "{YSEQID}_bwa_mem_{REF}_sorted.bam",
        BAI = results_prefix / "{YSEQID}_bwa_mem_{REF}_sorted.bam.bai",
        IDXSTATS = results_prefix / "{YSEQID}_bwa_mem_{REF}_sorted.bam.idxstats.tsv"
    conda:
        env_path / "mapping.yaml"
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
        BAM = results_prefix / "{YSEQID}_bwa_mem_{REF}_sorted.bam",
        IDXSTATS =  results_prefix / "{YSEQID}_bwa_mem_{REF}_sorted.bam.idxstats.tsv"

    output:
        STATS = results_prefix / "{YSEQID}_{REF}_mapping_stats.txt"
    conda:
        env_path / "mapping.yaml"
    threads: 
        workflow.cores * 1
    shell:
        """
        bamstats -u -i {input.BAM} -o {output.STATS}
        """
        