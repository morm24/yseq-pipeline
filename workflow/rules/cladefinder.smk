rule preprocessing_cladefinder_derived:
    input:
        DERIVED_VCF =   "results/chrY_derived_{YSEQID}_{REF}.vcf.gz",
    output:
        POSITIVE_TXT =  "results/{YSEQID}_{REF}_positives.txt",
    shell:
        """
	    bcftools query -f '%ID,' {input.DERIVED_VCF} | sed ':a;N;$!ba;s/\\n//g' > {output.POSITIVE_TXT} &
        """

rule preprocessing_cladefinder_ancestral:
    input:
        ANCESTRAL_VCF = "results/chrY_ancestral_{YSEQID}_{REF}.vcf.gz"
    output:
        NEGATIVE_TXT =  "results/{YSEQID}_{REF}_negatives.txt"
    shell:
        """
	    bcftools query -f '%ID,' {input.ANCESTRAL_VCF} | sed ':a;N;$!ba;s/\\n//g' > {output.NEGATIVE_TXT} 
        """

rule check_HG:
    input:
        POSITIVE_TXT =  "results/{YSEQID}_{REF}_positives.txt",
        NEGATIVE_TXT =  "results/{YSEQID}_{REF}_negatives.txt",
        YFULLTREE =     "resources/tree/latest_YFull_YTree.json"
    output:
        temp("results/{YSEQID}_{REF}cladeFinderOutput.csv"),
        
    shell:
        """
        python3 workflow/scripts/script_templates/cladeFinder.py {input.YFULLTREE} {input.POSITIVE_TXT} {input.NEGATIVE_TXT} {output[0]}
        """

rule save_HG:
    input:
        CF_CSV =        "results/{YSEQID}_{REF}cladeFinderOutput.csv",
    output:
        HAPLO_DATA =    "results/{YSEQID}_{REF}haploData.txt",
        HAPLO_GROUP =   temp("results/{YSEQID}_{REF}haploGroup")
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
        POSITIVES_TXT = "results/{YSEQID}_{REF}_positives.txt",
        NEGATIVES_TXT = "results/{YSEQID}_{REF}_negatives.txt",
        CLEANED_VCF =   "results/chrY_cleaned_{YSEQID}_{REF}.vcf.gz",
        IDXSTATS =      "results/{YSEQID}_bwa_mem_{REF}_sorted.bam.idxstats.tsv",
        HAPLO_GROUP =   "results/{YSEQID}_{REF}haploGroup"  


    output:
        PHYLOEQ_SNPS = "results/{YSEQID}_{REF}_PHYLOEQ_SNPS.tsv",
        DOWNSTR_SNPS = "results/{YSEQID}_{REF}_DOWNSTR_SNPS.tsv"
    shell:
        """
        YFULLHG=$(head -n 1 {input.HAPLO_GROUP})
        python3 workflow/scripts/script_templates/getEquivalentAndDownstreamSNPs.py {input.YFULLTREE} "$YFULLHG" {input.CLEANED_VCF} {input.POSITIVES_TXT} {input.NEGATIVES_TXT} {output.PHYLOEQ_SNPS} {output.DOWNSTR_SNPS} {input.IDXSTATS}

        """
    
