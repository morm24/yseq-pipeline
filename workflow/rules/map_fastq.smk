#index teh reference fasta File with minimap2
rule index_refseq_minimap2: 
    input:
        REFSEQ = ref_prefix / "{REF}/{REF}.fa"
    output:
        INDEX = protected(ref_prefix / "{REF}/{REF}.fa.mmi")
    log:
        ref_prefix / "{REF}/{REF}_minimap2_index.log"
    benchmark:
        ref_prefix / "{REF}/{REF}_minimap2_index.benchmark"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        minimap2 {input.REFSEQ} -d {output.INDEX} > {log} 2>&1
        """


#map the Sample reads to the reference fasta file with minimap2
rule map_minimap2:
    input:
        READS = multiext(str(sample_prefix / "{YSEQID}_"), "R2.fastq.gz", "R1.fastq.gz"),
        REFSEQ = multiext(str(ref_prefix / "{REF}/{REF}"),".fa", ".fa.mmi")
    output: 
        BAM = temp(results_prefix / "mapping" / "{YSEQID}_minimap2_{REF}.bam")
    conda:
        "../envs/mapping.yaml"
    log:
        results_prefix / "mapping" / "logs" / "{YSEQID}_minimap2_{REF}.log"
    benchmark:
        results_prefix / "mapping" / "benchmark" / "{YSEQID}_minimap2_{REF}.benchmark"
    params:
        "-a -x sr" if (config["READS"] == "short") else "-ax map-ont" 
    threads: workflow.cores
    shell:
        """
        (minimap2 {params} -t {threads} {input.REFSEQ[1]} {input.READS[0]} {input.READS[1]} | 
        samtools view -@ {threads} -b -t {input.REFSEQ[0]} -o {output.BAM} - ) > {log} 2>&1
        """


#index the reference fasta file with bwa
rule index_refseq_bwa:
    input:
        REFSEQ = ref_prefix / "{REF}/{REF}.fa"
    output:
        protected( multiext(str(ref_prefix / "{REF}/{REF}") ,".fa.bwt", ".fa.amb", ".fa.ann", ".fa.pac", ".fa.sa"))
    log:
        ref_prefix / "{REF}/{REF}_bwa_index.log"
    benchmark:
        ref_prefix / "{REF}/{REF}_bwa_index.benchmark"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        bwa index {input.REFSEQ} > {log} 2>&1
        """


rule map_bwa:
    input:
        READS = multiext(str(sample_prefix / "{YSEQID}_"), "R2.fastq.gz", "R1.fastq.gz"),
        REFSEQ= multiext(str(ref_prefix / "{REF}/{REF}"),".fa", ".fa.bwt", ".fa.amb", ".fa.ann", ".fa.pac", ".fa.sa")
    output: 
        BAM = temp(results_prefix / "mapping" / "{YSEQID}_bwa-mem_{REF}.bam")
    conda:
        "../envs/mapping.yaml"
    log:
        results_prefix / "mapping" / "logs" / "{YSEQID}_bwa-mem_{REF}.log"
    benchmark:
        results_prefix / "mapping" / "benchmark" / "{YSEQID}_bwa-mem_{REF}.benchmark"
    threads: workflow.cores
    shell:
        """
        (bwa mem -M -t {threads} {input.REFSEQ[0]} {input.READS[0]} {input.READS[1]} | 
        samtools view -@ {threads} -b -t {input.REFSEQ[0]} -o {output.BAM} - ) > {log} 2>&1
        """


#sort and index the BAM file (mapped reads) for faster processing in future steps
rule sort_and_index:
    input:
        BAM = str(results_prefix / "mapping" / "{YSEQID}_bwa-mem_{REF}.bam") if (config["MAPPER"] == "bwa") else (results_prefix / "mapping" / "{YSEQID}_minimap2_{REF}.bam" )
    output:
        SORTED_BAM = results_prefix / "mapping" / "{YSEQID}_{REF}_sorted.bam",
        BAI = results_prefix / "mapping" / "{YSEQID}_{REF}_sorted.bam.bai",
        IDXSTATS = results_prefix / "mapping" / "{YSEQID}_{REF}_sorted.bam.idxstats.tsv"
    conda:
        "../envs/mapping.yaml"
    log:
        results_prefix / "mapping" / "logs" / "{YSEQID}_{REF}_samtools_sort.log"
    benchmark:
        results_prefix / "mapping" / "benchmark" / "{YSEQID}_{REF}_samtools_sort.benchmark"
    params:
        SORTDIR = tmp_prefix
    threads: workflow.cores 
    shell:
        """
        mkdir -p {params.SORTDIR}
        samtools sort -@ {threads} -T {params.SORTDIR} -o {output.SORTED_BAM} {input.BAM}    > {log} 2>&1
	    samtools index -@ {threads} {output.SORTED_BAM}                                   >> {log} 2>&1
	    samtools idxstats {output.SORTED_BAM} | tee {output.IDXSTATS}                     >> {log} 2>&1
        """


#Get Mapping statistics of the bam file for quality control
#swap the code with bamstaistics or sth. like that
rule get_mapping_statistics:
    input:
        BAM = results_prefix / "mapping" / "{YSEQID}_{REF}_sorted.bam",
        IDXSTATS =  results_prefix / "mapping" / "{YSEQID}_{REF}_sorted.bam.idxstats.tsv"

    output:
        STATS = results_prefix / "mapping" / "{YSEQID}_{REF}_mapping_stats.txt"
    conda:
        "../envs/mapping.yaml"
    log:
        results_prefix / "mapping" / "logs" / "{YSEQID}_{REF}_bamstats.log"
    benchmark:
        results_prefix / "mapping" / "benchmark" / "{YSEQID}_{REF}_bamstats.benchmark"
    threads: 
        1
    shell:
        """
        bamstats -u -i {input.BAM} -o {output.STATS} > {log} 2>&1
        """


#extract the mtDNA and Y chromosome reads from the BAM file for valifating the results outside the pipeline   
rule seperate_BAM:
    input:
        BAM = results_prefix / "mapping" / "{YSEQID}_{REF}_sorted.bam"
        
    output:
        CHR_BAM = results_prefix / "mapping" / "{YSEQID}_{REF}_{CHR}.bam",
        CHR_BAI = results_prefix / "mapping" / "{YSEQID}_{REF}_{CHR}.bam.bai"

    conda:
        "../envs/mapping.yaml"
    log:
        results_prefix / "mapping" / "logs" / "{YSEQID}_separate_{REF}_{CHR}.log"
    benchmark:
        results_prefix / "mapping" / "benchmark" / "{YSEQID}_separate_{REF}_{CHR}.benchmark"
    threads: 
        workflow.cores * 1
    shell:
        """
        (samtools view -@ {threads} -b -o {output.CHR_BAM} {input.BAM} {wildcards.CHR}) > {log} 2>&1
        samtools index {output.CHR_BAM} >> {log} 2>&1
        """
