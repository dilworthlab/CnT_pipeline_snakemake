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

1) The pipeline should be executed separately for each group of samples using the same antibody. This is because the parameters for Spike-in normalization are dependant on the entire library. Running the pipeline on files from different histone or TF CnT libraries will result in incorrect normalization.

Input files:
This pipeline starts with raw sequencing FASTQ files, R1 and R2 for each sample.

Format:
Fastq.gz files with the following naming format:

```
{SAMPLE-NAME}_{READ}.fastq.gz
```


Output files:

```
|--Analysis_Results   # Main results directory
|  |-- QC_Rawreads
|  |-- Trimming
|  |-- primary_alignment
|  |-- RPGC_and_Unnormalized_bws
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

Metadate file:
Fill in Samples.tsv. This file allows the user to organize metadata for the library. The minimum information that MUST be provided is:

1- Sample name: this should correspond to the name of the file i.e. ```{SAMPLE-NAME}_{READ}.fastq.gz```

2- Condition: The two options here are: 'IgG' or 'TargetFile'. This is important because if IgG is specified, it will it use it as a control file when calling peaks.



# STEP 1: Get repository
- clone repo

``` git clone https://github.com/dilworthlab/CnT_pipeline_snakemake.git ```

- move raw FASTQ files into the directory

The resulting directory will contain the following:

```
|-- Snakefile
|
|-- config  
|   |-- config.yaml
|   |-- Samples.tsv
|
|-- envs_conda
|   |-- snakemake_env.yaml    
|
|-- envs_wo_conda
|   |-- Requirements.txt
|   |-- modules_step2_option1.txt
|
|-- Reference_files
|   |-- Spikein_indices
|       |-- Amp_pbluescript      # (https://www.addgene.org/vector-database/1946/)
|       |-- EcoliK12_index       #(U00096.3 - https://www.ncbi.nlm.nih.gov/nuccore/545778205)
|   |-- cookiecutter
|       |-- default_res_config.yaml
|
|-- Scripts
|
|-- ** YOUR FASTQ READS **

```




# STEP 2: cookiecutter installation and SLURM profile configuration

 - Install cookiecutter - https://cookiecutter.readthedocs.io/en/1.7.2/installation.html

```
 pip install cookiecutter
```
 - Get SLURM profile - https://github.com/Snakemake-Profiles/slurm

```
mkdir ~/.config/snakemake_CnT
cd  ~/.config/snakemake_CnT
cookiecutter https://github.com/Snakemake-Profiles/slurm.git

# In the setup prompt, add your <profile name> (e.g. Slurm_CnT). Others can be left empty.
```

- add cluster details
```
# change dir
cd ~/.config/snakemake_CnT/Slurm_CnT (# profile name)

# open config.yaml and add this line at the end: 'jobs: 100', close file.
vim config.yaml

# mv default_res_config.yaml into this folder
mv  /path/to/repo/Reference_files/cookiecutter/default_res_config.yaml  ~/.config/snakemake_CnT/Slurm_CnT

# edit slurm-submit.py
add "/path/to/default_res_config.yaml" to:
CLUSTER_CONFIG = <HERE> # under "cookiecutter arguments"
```

For more information, refer to this very informative blog post: http://bluegenes.github.io/Using-Snakemake_Profiles/





# STEP 3: Snakefile Configuration

##  Download Bowtie2 index from iGenomes
https://support.illumina.com/sequencing/sequencing_software/igenome.html
```
# change directory to Reference_files
cd Reference_files

# E.g. for Mus musculus (UCSC mm10) genome:
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz

tar -xvf Mus_musculus_UCSC_mm10.tar.gz

```

##  Configure workflow using config.yaml
- Go back to the directory containing the FASTQ files for the Cut and Tag library and contents of the repo.
- open ./config/config.yaml in a text editor and modify the following lines:
```
# Primary genome index location

## e.g. Mus musculus (UCSC mm10)
genome_index: ./Reference_files/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome

## any other index
genome_index: ./Reference_files/path/to/genome/index

# Change "effective genome size" - default is for mouse
bamCov_RPGC:

# Spike-in index location

### Default - sequence for Ampr
Spikein_index: ./Reference_files/Spikein_indices/Amp_pbluescript/Amp_index/Amp_pBlue

### E coli K12
Spikein_index: ./Reference_files/Spikein_indices/EcoliK12_index/EcoliK12Index/EcoliK12

### Any other index
Spikein_index: /path/to/Spikein_index/{index_prefix}



# Picard location/Command
PicardLoc: "java -jar /path/to/picard.jar" # don't forget quotes

# SEACR location
SEACRLoc: /path/to/SEACR-master/SEACR_1.3.sh  # No quotes
```





# STEP 4: Run Pipeline with snakemake

In the directory, call Snakemake with ```--profile <profile.name>```.


```
snakemake --profile Slurm_CnT # with example profile name "Slurm_CnT"
```



### Troubleshooting Snakemake:

- You can do a 'dry-run' with Snakemake. This allows you to check the workflow before running it.
```
snakemake -np --profile Slurm_CnT
```

- Snakemake sometimes locks the working directory. The directory can be unlocked with --unlock.
```
snakemake --unlock --profile Slurm_CnT
```

For more information, refer to the very well documented Snakemake docs - https://snakemake.readthedocs.io/en/stable/project_info/faq.html



# Extra Notes

## Spike-in normalization

Normalization factors for a set of libraries associated with a particular histone or transcription factor
is calculated using reads mapped to the Spike-in sequence. Normalization factors for each sample is
calculated as:

<img src="https://latex.codecogs.com/gif.latex?Normalization&space;Factor&space;=&space;{\frac{Spike_{min}&space;}&space;{Spike_{Sample}&space;}&space;}"/>

where:

1) Spike_min = No. of reads corresponding to the sample with minimum number of mapped
reads within the set,
2) Spike_Sample = No. of total mapped reads corresponding to the sample for which the normalization
factor is being calculated.

# References

1- Hatice S. Kaya-Okur, Steven J. Wu, Christine A. Codomo, Erica S. Pledger, Terri D. Bryson, Jorja G. Henikoff, Kami Ahmad, and Steven Henikoff. CUT&Tag for efficient epigenomic profiling of small samples and single cells. Nature Communications, 10, April 2019.

2- Ye Zheng. CUT&Tag Data Processing and Analysis Tutorial.
