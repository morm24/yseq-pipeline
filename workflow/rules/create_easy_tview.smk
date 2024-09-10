rule create_easy_tview:
    input:
        SORTED_BAM =    "results/{YSEQID}_bwa_mem_{REF}_sorted.bam",
        REFSEQ = "resources/refseq/{REF}/{REF}.fa"
    output:
        "results/tview_{YSEQID}_{REF}.sh"
    shell:
        """
        echo "#!/bin/bash" > {output}
        echo "samtools tview {input.SORTED_BAM} {input.REFSEQ}" >> {output}
        chmod a+x {output}
        """