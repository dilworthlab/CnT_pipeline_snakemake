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
1) This pipeline should be run in a High Performance Computing (HPC) environment for best results.

2) The pipeline should be executed separately for each UNIQUE CUT&Tag library. This is because the parameters for Spike-in normalization are dependant on the entire library. Running the pipeline on files from different histone or TF CnT libraries will result in incorrect normalization.

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



# Quickstart

## STEP 1:  Get repository
- clone repo in the directory containing the raw fastq.gz reads. The resulting directory will contain the following:

``` git clone https://github.com/dilworthlab/CnT_pipeline_snakemake.git ```

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
|
|-- Reference_files
|   |-- cookiecutter
|       |-- default_res_config.yaml
|
|-- Scripts
```

## STEP 2: Virtual Environments
- Activate a conda or any other virtual environment.
- Install tools and dependencies based on either the provided enviroment.yaml files (with conda) or requirements.txt file (without conda).

Note: ComputeCanada requests not using conda environments (See: https://docs.computecanada.ca/wiki/Anaconda/en), in which case, the following steps can be followed:

### Option 1: HPC virtual environemnt (without conda)
```
# load all relevant modules
module load python/3.8.2
module load #other-modules

# make and activate virtual environment (e.g. myenv)
virtualenv --no-download ~/myenv
source ~/myenv/bin/activate


# install required software
pip install --no-index -r 2021_requirements.txt
```

### Option 2: Using conda Environment
```
#  create conda environment
conda env create -f envs_conda/snakemake_env.yaml
```
# Note: Consider using Tmux when running the pipeline.  

## STEP 3: cookiecutter installation and SLURM profile configuration

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
add "/path/to/repo/Reference_files/cookiecutter/default_res_config.yaml" to:
CLUSTER_CONFIG = <HERE> # under "cookiecutter arguments"
```

For more information, refer to this very informative blog post: http://bluegenes.github.io/Using-Snakemake_Profiles/


## STEP 4: Configure workflow using config.yaml
- Go back to the directory containing the FASTQ files for the Cut and Tag library and contents of the repo.
- open ./config/config.yaml in a text editor and modify the following lines:
```
PicardLoc: /path/to/Picard/jarfile

Spikein_index: /path/to/Spikein/index

genome_index: /path/to/genome/index

SEACRLoc: /path/to/SEACR-master/SEACR_1.3.sh
```

## STEP 5: Sample set-up
Fill in Samples.tsv. This file allows the user to organize metadata for the library. The minimum information that MUST be provided is:

1- Sample name: this should correspond to the name of the file i.e. ```{SAMPLE-NAME}_{READ}.fastq.gz```

2- Condition: The two options here are: 'IgG' or 'TargetFile'. This is important because if IgG is specified, it will it use it as a control file when calling peaks.


## STEP 6 : Run pipeline with snakemake
Once steps 1 to 5 are completed, you are finally ready to run the pipeline. In the directory, call Snakemake with ```--profile <profile.name>```. E.g.


```
snakemake --profile Slurm_CnT
```



### Troubleshooting Snakemake:

- You can do a 'dry-run' with Snakemake. This allows you to see if the workflow before running it.
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

2- Hatice S. Kaya-Okur, Steven J. Wu, Christine A. Codomo, Erica S. Pledger, Terri D. Bryson, Jorja G. Henikoff, Kami Ahmad, and Steven Henikoff. CUT&Tag for efficient epigenomic profiling of small samples and single cells. Nature Communications, 10, April 2019.

1- Ye Zheng. CUT&Tag Data Processing and Analysis Tutorial.
