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