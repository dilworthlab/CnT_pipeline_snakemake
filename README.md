# Data analysis pipeline for CUT&Tag data


# Summary
This is an analysis pipeline for CUT&Tag data implemented in Snakemake. Cleavage Under Targets and Tagmentation or CUT&Tag is a method developed by the Henikoff lab for epigenomic profiling that employs a Tn5 transposase. A pre-assembled complex consisting of a Tn5 transposase loaded with sequencing adapters and protein A is directed to an antibody-bound, chromatin-associated protein-of-interest where it performs tagmentation. The resulting fragments are PCR amplified and sequenced.

The following pipeline can be used to perform quality-control, alignment as well as Spike-in normalization. The steps are outlined below:

1) MD5 checks of FASTQ reads
2) QC of FASTQ reads
3) Alignment to genome (e.g. mouse, human etc )
4) Alignment to Spike-in genome (e.g. any artificially spiked sequence or carry-over E.coli)
5) Filtering BAMs
6) Spike-in normalization
7) Peak-calling

<img width="1387" alt="Screen Shot 2021-04-08 at 9 02 32 AM" src="https://user-images.githubusercontent.com/20444993/114031286-37bb7a00-9849-11eb-8ddd-397d66b7e3b0.png">




# General Usage and Important considerations

**1) Execute pipeline separately for each group of samples using the same antibody.**

This is because the parameters for Spike-in normalization are dependant on the entire library. Running the pipeline on files from different histone or TF CnT libraries will result in incorrect normalization.

**2) Input Files**
This pipeline starts with raw sequencing FASTQ files, R1 and R2 for each samples. The format **MUST** follow the naming format below:

```
{SAMPLE-NAME}_{READ}.fastq.gz
```

**3) Md5sum hashes**

A text file containing md5sum hashes and corresponding file names **MUST** be provided. The file **Should** be called "md5sum.txt".

**4) Metadata Sheet: Samples.tsv**

The Samples.tsv file allows the user to organize metadata for the library. The minimum information that **MUST** be provided is:

- **Sample name**: this should correspond to the name of the file i.e. ```{SAMPLE-NAME}_{READ}.fastq.gz```

- **Condition**: The two options here are: 'IgG' or 'TargetFile'. This is important because if IgG is specified, it will it use it as a control file when calling peaks.


**5) To download**
Bowtie2 index for either Mus musculus or Homo sapiens

This can be directly downloaded from iGenomes
- Mouse
http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz

- Human
http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz


# STEP 1: Get repository

**1) Set up a Virtual Environment (with or without Conda)**

### With Conda (Make sure anaconda/miniconda is installed)
```
conda env create -f snakemake_CnT_env.yaml
```
```
conda activate CnT_env
```

**2) Clone repo**

``` git clone https://github.com/dilworthlab/CnT_pipeline_snakemake.git ```

**3) Move raw FASTQ files into the directory**

The resulting directory will contain the following:

```
|-- Snakefile
|
|-- config  
|   |-- config.yaml
|   |-- Samples.tsv
|
|-- envs_conda
|   |-- snakemake_CnT_env.yaml    
|
|-- envs_wo_conda
|   |-- Software_Requirements.txt
|
|-- Reference_files
|   |-- Spikein_indices
|       |-- Amp_pbluescript      # (https://www.addgene.org/vector-database/1946/)
|       |-- EcoliK12_index       #(U00096.3 - https://www.ncbi.nlm.nih.gov/nuccore/545778205)
|   |-- cookiecutter
|       |-- default_res_config.yaml
|       |-- config.yaml
|
|-- Scripts
|
|-- ** YOUR FASTQ READS DIR **
|       |-- *.fastq.gz
|       |-- md5sum.txt

```




# STEP 2: Cookiecutter installation and SLURM profile configuration - ONLY NEEDS TO BE DONE ONCE


```
mkdir -p ~/.config/snakemake
```

```
cd  ~/.config/snakemake
```
```
cookiecutter https://github.com/Snakemake-Profiles/slurm.git
```

# In the setup prompt, add your <profile name> (e.g. Slurm_CnT). Others can be left empty.


```
cp /path/to/repo/Reference_files/cookiecutter/default_res_config.yaml ~/.config/snakemake/Slurm_CnT
```

```
cp /path/to/repo/Reference_files/cookiecutter/config.yaml ~/.config/snakemake/Slurm_CnT
```

```
cd ~/.config/snakemake/Slurm_CnT
```

```
# edit slurm-submit.py

add "/path/to/default_res_config.yaml" to:

CLUSTER_CONFIG = <HERE> # (under "cookiecutter arguments")

# Make sure you enter the ABSOLUTE PATH
# wrong = "~/path/to/dir" ('~' is not valid python syntax)
# right = "/home/user/path/to/dir"
```

For more information, refer to this very informative blog post: http://bluegenes.github.io/Using-Snakemake_Profiles/






# STEP 3: Snakefile Configuration



**2) Configure workflow using config.yaml**

- open ./config/config.yaml in a text editor and modify the following lines:
```



# Species (Default: Mus musculus)
"Mus musculus" or "Homo Sapiens"

# Spikein (Default: Amp)
"Amp" or "Bacteria"

# Primary genome index
bowtie2_index: "/path/to/bowtie2/index/prefix" # don't forget prefix

# Picard location/Command
PicardLoc: "java -jar /path/to/picard.jar" # don't forget quotes

# SEACR location
SEACRLoc: /path/to/SEACR-master/SEACR_1.3.sh  # No quotes

# read_lenght (Default = "50")


```





# STEP 4: Run Pipeline with Snakemake

**In the directory containing the Snakefile, call Snakemake with ```--profile <profile.name>```:**


```
snakemake --profile Slurm_CnT  # with example profile name "Slurm_CnT"
```



#### Troubleshooting Snakemake:

- You can do a 'dry-run' with Snakemake. This allows you to check the workflow before running it.
```
snakemake -np --profile Slurm_CnT
```

- Snakemake sometimes locks the working directory. The directory can be unlocked with --unlock.
```
snakemake --unlock --profile Slurm_CnT
```

For more information, refer to the very well documented Snakemake docs - https://snakemake.readthedocs.io/en/stable/project_info/faq.html




## Output files


```
|--Analysis_Results   # Main results directory
|  |-- QC_Rawreads
|  |-- Trimming
|  |-- primary_alignment
|  |-- RPGC_and_Unnormalized_bws  # Non-spikein normalized bigwigs for reference purposes
|  |-- Spikein_alignment
|  |-- Spikein_normalized_bws_bdgs
|  |-- Peaks
|
|-- All_output    # All other output
|  |-- Mapped_reads
|  |-- Processed_reads
|  |-- Spike_mapped_reads
|  |-- Trimmed_reads
|
|-- logs  # Logs of all jobs that were run
```



# References

1- Hatice S. Kaya-Okur, Steven J. Wu, Christine A. Codomo, Erica S. Pledger, Terri D. Bryson, Jorja G. Henikoff, Kami Ahmad, and Steven Henikoff. CUT&Tag for efficient epigenomic profiling of small samples and single cells. Nature Communications, 10, April 2019.

2- Ye Zheng. CUT&Tag Data Processing and Analysis Tutorial.
