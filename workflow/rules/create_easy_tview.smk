rule create_easy_tview:
    input:
        SORTED_BAM =    results_prefix / "mapping" / "{YSEQID}_{REF}_sorted.bam",
        REFSEQ = ref_prefix / "{REF}/{REF}.fa"
    output:
        results_prefix / "mapping" / "tview_{YSEQID}_{REF}.sh"
    conda:
        "../envs/bam_process.yaml"
    benchmark:
        results_prefix / "mapping" / "benchmark" / "{YSEQID}_{REF}_create_tview.benchmark"
    shell:
        """
        echo "#!/bin/bash" > {output}
        echo "samtools tview {input.SORTED_BAM} {input.REFSEQ}" >> {output}
        chmod a+x {output}
        """