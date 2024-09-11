rule download_snps:
    output:
        SNPS =          temp("resources/tmp/snps_{REF}.vcf.gz"),
        SNPS_TBI =      temp("resources/tmp/snps_{REF}.vcf.gz.tbi")
    shell:
        """
        if [ "{wildcards.REF}" == "hg38" ]; then
            wget -O {output[0]} http://ybrowse.org/gbrowse2/gff/snps_hg38.vcf.gz
            wget -O {output[1]} http://ybrowse.org/gbrowse2/gff/snps_hg38.vcf.gz.tbi
        #elif [ "{wildcards.REF}" == "hs1" ]; then
            #wget -O {output.SNPS} http://ybrowse.org/gbrowse2/gff/snps_hs1.vcf.gz                  still needs to be created
            #wget -O {output.SNPS_TBI} http://ybrowse.org/gbrowse2/gff/snps_hs1.vcf.gz.tbi          still needs to be created
        else
            echo "Invalid reference: {wildcards.REF}"
            exit 1
        fi
        """
