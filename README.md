# CutandTag_2021_snakemake


# Summary
Cut and Tag pipeline implemented in snakemake.

# General Usage and Important considerations
1) This pipeline should be run in a High Performance Computing (HPC) environment for best results.

2) The pipeline should be executed separately for each UNIQUE Cut and Tag library. This is because the parameters for Spike-in normalization are dependant on the entire library. Running the pipeline on files from different histone or TF CnT libraries will result in incorrect normalization.

Input files:
This pipeline starts with raw sequencing FASTQ files, R1 and R2 for each sample.

Format:
Fastq.gz files with the following naming format:

```
[SAMPLE-NAME]_[READ].fastq.gz
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



# QUICKSTART

STEP 1
- git clone repo in the directory containing the raw fastq.gz reads. The resulting directory will contain the following:

```
|-- Snakefile
|
|-- config  
|   |-- config.yaml
|   |-- Samples.tsv
|
|-- envs_conda
|   |-- Tool1.yaml
|   |-- Tool2.yaml
|
|-- envs_wo_conda
|   |-- Requirements.txt
|
|-- Scripts
```

STEP 3
- Activate a conda or any other virtual environment.
- Install tools and dependencies based on either the provided enviroment.yaml files (with conda) or requirements.txt file (without conda).

Note: ComputeCanada requests not using conda environments (See: https://docs.computecanada.ca/wiki/Anaconda/en), in which case, the following steps can be followed:

```
module load python/3.8.2
module load #other-modules
```
```
source ~/ENV/bin/activate
```

```
pip install --no-index -r 2021_requirements.txt
```

STEP 4
cookiecutter installation and SLURM profile configuration
 - Install cookiecutter - https://cookiecutter.readthedocs.io/en/1.7.2/installation.html

```
 pip install cookiecutter
```
 - Get SLURM profile - https://github.com/Snakemake-Profiles/slurm

```
mkdir ~/.config/snakemake_CnT
cd  ~/.config/snakemake_CnT
cookiecutter https://github.com/Snakemake-Profiles/slurm.git
```

In the setup prompt, add your <profile name> (e.g. Slurm_CnT). Others can be left empty.

For more information, refer to this very informative blog post: http://bluegenes.github.io/Using-Snakemake_Profiles/

STEP 5
- Go back to the directory containing the FASTQ files for the Cut and Tag library and contents of the repo.
- open ./config/config.yaml in a text editor and modify the following lines:
```
PicardLoc: /path/to/Picard/jarfile

Spikein_index: /path/to/Spikein/index

genome_index: /path/to/genome/index

SEACRLoc: /path/to/SEACR-master/SEACR_1.3.sh
```

STEP 6
Fill in Samples.tsv. This file allows the user to organize metadata for the library. The minimum information that MUST be provided is:
1- Sample name: this should correspond to the name of the file i.e. [SAMPLE-NAME]_[READ].fastq.gz

2- Condition: The two options here are: 'IgG' or 'TargetFile'. This is important because if IgG is specified, it will it use it as a control file when calling peaks.


STEP 7
Once steps 1 to 6 are completed, you are finally ready to run the pipeline. In the directory, call snakemake with --profile <profile.name>. E.g.


```
snakemake --profile Slurm_CnT
```

Troubleshooting Snakemake:

- You can do a 'dry-run' with snakemake. This allows you to see if the workflow is constructed properly.
```
snakemake -np --profile Slurm_CnT
```

- Snakemake sometimes locks the working directory. The directory can be unlocked with --unlock.
```
snakemake --unlock --profile Slurm_CnT
```

For more information, refer to the very well documented Snakemake docs - https://snakemake.readthedocs.io/en/stable/project_info/faq.html
