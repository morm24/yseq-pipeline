# YSEQ Y-DNA Analysis Pipeline

## Version
**Current Version:** 0.0.0

## Overview
The YSEQ Y-DNA Analysis Pipeline is designed to process paired ended fastq.gz files and extract the Y- and MT-DNA Haplogroups.

## Dependencies
- Python 3.6+
- Snakemake
- conda

- bcftools
- samtools
- samtools
- tabix



## Installation
1. **Clone the repository:**
    ```sh
    git clone https://github.com/morm24/yseq-pipeline.git
    cd yseq-pipeline
    ```

2. **Install Python dependencies:**
    ```sh
    pip install -r requirements.txt
    ```

3. **Ensure the right tools are installed:**
    ```sh
    sudo apt-get install bcftools  samtools python3 snakemake
    ```

## Usage
1. **Prepare the input files:**
    
    ### Add your sample
    Place your samples into the `resources/sample` folder. 
    They have to be Paired end `fastq.gz` files.
    Their name has to be: `{SampleID}_R1.fastq.gz` and `{SampleID}_R2.fastq.gz` 
    
    ### Add the reference sequence
    Place the chosen reference sequence into `resources/refseq/{ref}/{ref}.fa`.
    The Name of the foldername and the fasta seqence have to be the same.

    ### add both to the samples.csv
    Add all samples with their ID to th file `config/samples.csv`. 
    write each sample into a line, separated with a comma.



2. **Run the pipeline:**
    To just run the Pipeline, type the following command with the amount of cores it should use:
    ```sh
    snakemake -pfc {cores}
    ```

3. **Check the output:**
    TBD


## Reference Sequences
Common reference sequences and their Download links are:
- **hg38:** [Download hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)
- **hg19:** [Download hg19](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz)
- **hs1:**  [Download hs1](https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz)

## YFull YTree Updates
To download the latest YFull YTree updates, follow these steps:
1. Visit the [YFull YTree website](https://www.yfull.com/tree/).
2. Download the latest JSON file from the resources section.
3. Place the downloaded JSON file in the `resources/tree` directory and rename it to `latest_YFull_YTree.json`.

## Example
Here is an example of how to run the pipeline with sample data:


## License
This project is licensed under the TBD License. See the [LICENSE](LICENSE) file for details.

## Contact
For any questions or issues, please open an issue on the [GitHub repository](https://github.com/morm24/yseq-pipeline/issues).

---