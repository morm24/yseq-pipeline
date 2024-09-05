
rule get_mapping_statistics:
    input:
        BAM = "results/YSEQID_sorted.bam",
        IDXSTATS =  "results/YSEQID_sorted.bam.idxstats.tsv"

    output:
        STATS = "results/YSEQID_mapping_stats.txt"
    shell:
        """
        
        "bamfile size: " > {output.STATS}
        BAM_FILESIZE=`du -kh "{input.BAM}" | cut -f1`
        echo "${{BAM_FILESIZE}}"   >> {output.STATS}
        
        "\navg read length: " >> {output.STATS}
        READLENGTH=$(samtools bam2fq -@ {threads} {input.BAM} | head -n 400000 | awk '{{if(NR%4==2) {{count++; bases += length}} }} END{{print(bases/count)}}')
        echo "${{READLENGTH}}" >> {output.STATS}

        "\nnumber of mapped reads: " >> {output.STATS}
        MAPPED_READS_HG38=$(perl addRows.pl "--column=2" "--filename={input.BAM}.idxstats.tsv")
        echo "${{MAPPED_READS_HG38}}" >> {output.STATS}
        
        "\nnumber of mapped bases: " >> {output.STATS}
        MAPPED_BASES_HG38=$(perl multiply.pl "--multiplicands=$MAPPED_READS_HG38,$READLENGTH" "-round")
        echo "${{MAPPED_BASES_HG38}}" >> {output.STATS}

        "\nnumber of unmapped reads: " >> {output.STATS}
        UNMAPPED_READS=$(perl addRows.pl "--column=3" "--filename={input.BAM}.idxstats.tsv")
        echo "${{UNMAPPED_READS}}" >> {output.STATS}

        "\nnumber of unmapped bases: " >> {output.STATS}
        UNMAPPED_BASES=$(perl multiply.pl "--multiplicands=$UNMAPPED_READS,$READLENGTH" "-round")
        echo "${{UNMAPPED_BASES}}" >> {output.STATS}

        "\nnumber of sequenced bases: " >> {output.STATS}
        SEQUENCED_BASES="$( printf '%s + %s\n' "$MAPPED_BASES_HG38" "$UNMAPPED_BASES" | bc )"
        echo "${{SEQUENCED_BASES}}" >> {output.STATS}

        COVERAGE_HG38=$(awk -v rl="${{READLENGTH}}" '{{x+=$2;m+=$3}}END{{print m*rl/x "x"}}' < {input.IDXSTATS})

        
        
        """