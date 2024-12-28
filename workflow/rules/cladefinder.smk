rule preprocessing_cladefinder_derived:
    input:
        DERIVED_VCF =   results_prefix  / "snp_calling" / "chrY_derived_{YSEQID}_{REF}.vcf.gz",
    output:
        POSITIVE_TXT =  results_prefix / "cladefinder" / "{YSEQID}_{REF}_positives.txt",
    conda:
        "../envs/bam_process.yaml"
    log: 
        results_prefix / "cladefinder" /  "cladefinder.txt.log"
    benchmark:
        results_prefix / "cladefinder" / "benchmark" / "{YSEQID}_{REF}_preprocessing_cladefinder_derived.benchmark"
    shell:
        """
	    (bcftools query -f '%ID,' {input.DERIVED_VCF} | sed ':a;N;$!ba;s/\\n//g' | tee {output.POSITIVE_TXT}) >> {log} 2>&1
        """


rule preprocessing_cladefinder_ancestral:
    input:
        ANCESTRAL_VCF = results_prefix  / "snp_calling" / "chrY_ancestral_{YSEQID}_{REF}.vcf.gz"
    output:
        NEGATIVE_TXT =  results_prefix / "cladefinder" / "{YSEQID}_{REF}_negatives.txt"
    conda:
        "../envs/bam_process.yaml"
    log: 
        results_prefix / "cladefinder" / "cladefinder.txt.log"
    benchmark:
        results_prefix / "cladefinder" / "benchmark" / "{YSEQID}_{REF}_preprocessing_cladefinder_ancestral.benchmark"
    shell:
        """
	    (bcftools query -f '%ID,' {input.ANCESTRAL_VCF} | sed ':a;N;$!ba;s/\\n//g' | tee {output.NEGATIVE_TXT}) >> {log} 2>&1
        """


rule check_HG:
    input:
        POSITIVE_TXT =  results_prefix / "cladefinder" / "{YSEQID}_{REF}_positives.txt",
        NEGATIVE_TXT =  results_prefix / "cladefinder" / "{YSEQID}_{REF}_negatives.txt",
        YFULLTREE =     "resources/tree/latest_YFull_YTree.json"
    output:
        results_prefix / "cladefinder" / "{YSEQID}_{REF}cladeFinderOutput.csv"
    conda:
        "../envs/bam_process.yaml"
    log: 
        results_prefix / "cladefinder" / "log" / "cladefinder.txt.log"
    benchmark:
        results_prefix / "cladefinder" / "benchmark" / "{YSEQID}_{REF}_check_HG.benchmark"
    shell:
        """
        (python workflow/scripts/cladeFinder.py {input.YFULLTREE} {input.POSITIVE_TXT} {input.NEGATIVE_TXT} {output[0]}) >> {log} 2>&1
        """


rule save_HG:
    input:
        CF_CSV =        results_prefix / "cladefinder" / "{YSEQID}_{REF}cladeFinderOutput.csv",
    output:
        HAPLO_DATA =    results_prefix / "cladefinder" / "{YSEQID}_{REF}haploData.txt",
        HAPLO_GROUP =   results_prefix / "cladefinder" / "{YSEQID}_{REF}haploGroup"
    run:
        # Initialize variables
        YFULLHG = "unknown"
        YFULLPATH = "unknown"
        line_counter = 0
        # Input and output file paths
        input_csv = input.CF_CSV
        output_file = output.HAPLO_DATA
        # Read the second line of the input CSV file
        with open(input_csv, 'r') as file:
            lines = file.readlines()
            if len(lines) > 1:
                second_line = lines[1].strip()
                columns = second_line.split('\t')
                YFULLHG = columns[0]
                YFULLPATH = columns[1]

        # Write the output file
        with open(output_file, 'w') as file:
            file.write(f"YFULLHG: {YFULLHG}\n")
            file.write(f"YFULLPATH: {YFULLPATH}\n")
        with open(output.HAPLO_GROUP, 'w') as file:
            file.write(f"{YFULLHG}\n")  


rule get_equivalent_and_downstream_SNPS:
    input:
        YFULLTREE =     "resources/tree/latest_YFull_YTree.json",
        POSITIVES_TXT = results_prefix / "cladefinder" / "{YSEQID}_{REF}_positives.txt",
        NEGATIVES_TXT = results_prefix / "cladefinder" / "{YSEQID}_{REF}_negatives.txt",
        CLEANED_VCF =   results_prefix / "snp_calling" / "chrY_cleaned_{YSEQID}_{REF}.vcf.gz",
        CLEANED_VCF_TBI =   results_prefix / "snp_calling" / "chrY_cleaned_{YSEQID}_{REF}.vcf.gz.tbi",
        IDXSTATS =      results_prefix / "mapping" / "{YSEQID}_bwa-mem_{REF}_sorted.bam.idxstats.tsv",
        HAPLO_GROUP =   results_prefix / "cladefinder" / "{YSEQID}_{REF}haploGroup"  
    output:
        PHYLOEQ_SNPS = results_prefix / "cladefinder" / "{YSEQID}_{REF}_PHYLOEQ_SNPS.tsv",
        DOWNSTR_SNPS = results_prefix / "cladefinder" / "{YSEQID}_{REF}_DOWNSTR_SNPS.tsv"
    conda:
        "../envs/bam_process.yaml" 
    log: 
        results_prefix / "cladefinder" / "log" / "getEQAndDownSNPs.txt.log"
    benchmark:
        results_prefix / "cladefinder" / "benchmark" / "{YSEQID}_{REF}_getEQAndDownSNPs.benchmark"
    shell:
        """
        YFULLHG=$(head -n 1 {input.HAPLO_GROUP})
        python workflow/scripts/getEquivalentAndDownstreamSNPs.py {input.YFULLTREE} "$YFULLHG" {input.CLEANED_VCF} {input.POSITIVES_TXT} {input.NEGATIVES_TXT} {output.PHYLOEQ_SNPS} {output.DOWNSTR_SNPS} {input.IDXSTATS} > {log} 2>&1

        """
    
