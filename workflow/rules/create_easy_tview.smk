rule create_easy_tview:
    input:
        SORTED_BAM =    results_prefix / "mapping" / "{YSEQID}_bwa_mem_{REF}_sorted.bam",
        REFSEQ = ref_prefix / "{REF}/{REF}.fa"
    output:
        results_prefix / "mapping" / "tview_{YSEQID}_{REF}.sh"
    conda:
        "bam_process"
    #log:
    #    results_prefix / "{YSEQID}_{REF}_create_tview.log"
    shell:
        """
        echo "#!/bin/bash" > {output}
        echo "samtools tview {input.SORTED_BAM} {input.REFSEQ}" >> {output}
        chmod a+x {output}
        """