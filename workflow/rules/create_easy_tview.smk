rule create_easy_tview:
    input:
        SORTED_BAM =    results_prefix / "{YSEQID}_bwa_mem_{REF}_sorted.bam",
        REFSEQ = ref_prefix / "{REF}/{REF}.fa"
    output:
        results_prefix / "tview_{YSEQID}_{REF}.sh"
    conda:
        env_path / "bam_process.yaml"
    shell:
        """
        echo "#!/bin/bash" > {output}
        echo "samtools tview {input.SORTED_BAM} {input.REFSEQ}" >> {output}
        chmod a+x {output}
        """