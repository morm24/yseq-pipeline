# YSEQ Y-DNA Analysis Pipeline

## Version
**Current Version:** 0.0.0

## Overview
The YSEQ Y-DNA Analysis Pipeline is designed to process paired ended fastq.gz files and extract the Y- and MT-DNA Haplogroups.

## Dependencies
- Python 3.6+
- Snakemake
- conda
- pandas
- Path

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

2. **Install Python:**
    ```sh
    apt install python3
    ```

3. **Install snakemake:**
    ```sh
    apt install snakemake
    ```
4. **Install conda & mamba:**
    ### Install and set up miniconda 
    ```sh
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    conda config --set channel_priority strict
    ```
    ### Install mamba:
    ```sh
    conda install -n base -c conda-forge mamba
    ```
5. **Install pandas & Path**
    '''sh
    pip install Pandas Path
    '''
## Usage
1. **Prepare the input files:**
    ### First Methode
    #### Add your sample
    Place your samples into the `resources/sample/{SampleID}` folder.
    The Filename name has to be: `{SampleID}_R1.fastq.gz` if paried end, the second files name must be `{SampleID}_R2.fastq.gz` 
    
    #### Add the reference sequence
    Place the chosen reference sequence into `resources/refseq/{ref}/{ref}.fa`.
    The Name of the foldername and the fasta seqence have to be the same.

    #### add both to the samples.csv
    Add all SampleIDs to the file `config/samples.csv`. 
    Separated with a "," state the name of the reference (file name without ending).    Example: 
    ```
    ID,REF
    63819,hs1
    ```
    #### configurate the settings
    Chose what type of reads you use, the mapping software and result directory in `config/config.yaml`

    ### Second Methode
    #### Add your resource folders
    Open `config/config.yaml`.
    Change the Sample, Reference and Results folders.
    Chose the read type, and mapping software.
    #### Add your sample
    Add all SampleIDs to the file `config/samples.csv`. 
    Separated with a "," note the reference sequence, you want to map and analyze the sample to.
    Example: 
    ```
    ID,REF
    63819,hs1
    ```
    #### configurate the settings



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
```sh 
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```
- **hg19:** [Download hg19](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz)     
```sh 
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
```
- **hs1:**  [Download hs1](https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz)        
```sh  
wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz
```

## YFull YTree Updates  TO DO: change to free version
To download the latest YFull YTree updates, follow these steps:
1. Visit the [YFull YTree website](https://www.yfull.com/tree/).
2. Download the latest JSON file from the resources section.
3. Place the downloaded JSON file in the `resources/tree` directory and rename it to `latest_YFull_YTree.json`.

## Example
When set up corectly, and placed the hs1 sequence into `resources/refseq/hs1/hs1.fa` the following example workflow should run without errors and should take about 5 minutes to finish (faster with more cores): 
`snakemake -pc 1`


## License
This project is licensed under the TBD License. See the [LICENSE](LICENSE) file for details.

## Contact
For any questions or issues, please open an issue on the [GitHub repository](https://github.com/morm24/yseq-pipeline/issues).

---