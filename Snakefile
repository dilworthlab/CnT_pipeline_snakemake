import pandas as pd
import numpy as np
import glob
import os
import fnmatch
import yaml
import math
from pathlib import Path


from snakemake.io import expand, glob_wildcards
from snakemake.utils import min_version
from snakemake.logging import logger

#------------------------------------------------------------------------------------------------------------------------
# Minimum snakemake version
min_version("5.24.1")

# configuration file
configfile:"./config/config.yaml"

# Tabular configuration
samples = pd.read_csv(config["sample_info"], "\t").set_index("TargetFiles")
#------------------------------------------------------------------------------------------------------------------------





#------------------------------------------------------------------------------------------------------------------------
# Multiqc configuration
home = str(Path.home())
datadict = {'log_filesize_limit': 2000000000}

files = glob.glob(os.path.join(home, ".*.yaml"))
if os.path.join(home, ".multiqc_config.yaml") in files:
    print('~/.multiqc_config.yaml is present')
    with open(os.path.join(home, ".multiqc_config.yaml"), 'r') as file:
        values = yaml.safe_load(file)

        if 'log_filesize_limit' in str(values):
            print('log_filesize_limit is already set')
        if values is None:
            with open(os.path.join(home, ".multiqc_config.yaml"), 'w') as file:
                docs = yaml.dump(datadict, file)
                print('added log_filesize_limit to multiqc log file')
        else:
            with open(os.path.join(home, ".multiqc_config.yaml"), 'r') as file:
                new_yaml = yaml.safe_load(file)
                print(type(new_yaml))
                new_yaml.update(datadict)


            with open(os.path.join(home, ".multiqc_config.yaml"),'w') as file:
                yaml.safe_dump(new_yaml, file)
                print('log_filesize_limit: 2000000000 is set')

else:
    with open(os.path.join(home, ".multiqc_config.yaml"), 'w') as file:
        documents = yaml.dump(datadict, file)
        print('made new .multiqc_config.yaml')

#------------------------------------------------------------------------------------------------------------------------





#------------------------------------------------------------------------------------------------------------------------
# Getting directory containing rawreads
wdir = os.getcwd()
for filepath, dirs, allfiles in os.walk(wdir):
    for file_ in allfiles:
        if file_.endswith(".fastq.gz"):
            READS_DIR = filepath
            ROOT_DIR = os.path.relpath(filepath, wdir)


logger.info(f'This is the RawReads dir: {READS_DIR}')




# Describing wildcards
SINGLE_READ = f'{READS_DIR}/{{fastqfile}}_{{read}}.fastq.gz'


READS = set(glob_wildcards(SINGLE_READ).read)
FASTQFILES = set(glob_wildcards(SINGLE_READ).fastqfile)

logger.info(f'Sample file format: {SINGLE_READ}')
logger.info(f'These are the sample names: {FASTQFILES}')
logger.info(f'Each sample has {READS}')
#------------------------------------------------------------------------------------------------------------------------





#------------------------------------------------------------------------------------------------------------------------
# Identify groups
ABGROUPS = samples['AntibodyGroup'].unique().tolist()

# All IGG files
ALL_IGGFILES = samples['IgGFile'].unique().tolist()
if 'None' in ALL_IGGFILES:
    ALL_IGGFILES.remove('None')

ALL_TARGETS = samples.index.tolist()



# Get target reads per group
SampleGroupDict = {}

for ab in ABGROUPS:
    SampleGroupDict_ = {}
    GROUP_TARGETS = list(set (samples.loc[samples['AntibodyGroup'] == ab ].index))
    GROUP_IGGFILES = list(set (samples.loc[samples['AntibodyGroup'] == ab ]['IgGFile']))

    if len(GROUP_IGGFILES) > 1:
        raise SyntaxError('More than one type of IgG file indicated for a group/groups, please edit Samplesheet.tsv')
    else:
        SampleGroupDict_["Targets"] = GROUP_TARGETS
        SampleGroupDict_["IgGfiles"] = GROUP_IGGFILES
        SampleGroupDict[ab] = SampleGroupDict_



for key, item in SampleGroupDict.items():
    T = item['Targets']
    I = item['IgGfiles']
    logger.info(f'Group: {key}, Targets: {T}, IgG: {I}')




#------------------------------------------------------------------------------------------------------------------------





#------------------------------------------------------------------------------------------------------------------------
# Spike
if config['Spikein'] == 'Amp':
    SPIKEINDEX = config['Spikein_index_amp']

if config['Spikein'] == 'Bacteria':
    SPIKEINDEX = config['Spikein_index_Ecoli']
#------------------------------------------------------------------------------------------------------------------------





#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

#localrules: Clean_up

rule all:
    input:
        "logs/Spikein_normalized_bws_bdgs/mockfile.txt",
        expand("Analysis_Results/Peaks/{fastqfile}.stringent.bed",fastqfile=ALL_TARGETS),
        "Analysis_Results/Spikein_alignment/Spike_alignment.html",
        "logs/cleanup.log",
        expand("Analysis_Results/Spikein_normalized_bws_bdgs/ScalingFactors_by_group/{Groupname}_ScalingFactors.csv", Groupname=ABGROUPS),
        "Analysis_Results/Spikein_normalized_bws_bdgs/ScalingFactors.csv",
        expand('Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs/{fastqfile}_Norm.bw', fastqfile=ALL_TARGETS),
        expand('Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/{fastqfile}_Norm.bedgraph', fastqfile=ALL_TARGETS),
        expand('Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs_wDups/{fastqfile}_Norm_wDups.bw', fastqfile=ALL_TARGETS),
        expand('All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.NoDups.bam.bai', fastqfile=FASTQFILES),
        expand('All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.bam.bai', fastqfile=FASTQFILES),
        expand('Analysis_Results/CPM_and_Unnormalized_bws/CPM_normalized_bws/{fastqfile}_CPM.bw', fastqfile=FASTQFILES),
        expand('Analysis_Results/CPM_and_Unnormalized_bws/bws_wo_normalization/{fastqfile}_wo.norm.bw', fastqfile=FASTQFILES),
        expand('Analysis_Results/CPM_and_Unnormalized_bws/bws_wo_norm_wDups/{fastqfile}_wo.norm_wDups.bw', fastqfile=FASTQFILES),
        expand('All_output/Spike_mapped_reads/{fastqfile}.coordsorted.bam',fastqfile=FASTQFILES),
        expand('All_output/Spike_mapped_reads/{fastqfile}.bam', fastqfile=FASTQFILES),
        expand('All_output/Spike_mapped_reads/{fastqfile}.coordsorted.bam.bai', fastqfile=FASTQFILES),
        expand('logs/Spike_Alignment/dupstats/{fastqfile}.dupMarked.bam', fastqfile=FASTQFILES),
        expand('logs/Spike_Alignment/dupstats/{fastqfile}_picard.dupMark.txt', fastqfile=FASTQFILES),
        "Analysis_Results/primary_alignment/filteringbamsStats.html",
        expand('logs/filtered_bams/{fastqfile}_picard.rmDup.txt', fastqfile=FASTQFILES),
        expand('All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.bam', fastqfile=FASTQFILES),
        expand("All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.NoDups.bam", fastqfile=FASTQFILES),
        "Analysis_Results/primary_alignment/Alignment_results.html",
        expand("logs/primary_alignment/PostAlignmentStats/dupstats/{fastqfile}.dupMarked.bam", fastqfile=FASTQFILES),
        expand("logs/primary_alignment/PostAlignmentStats/dupstats/{fastqfile}_picard.dupMark.txt", fastqfile=FASTQFILES),
        expand('All_output/Mapped_reads/{fastqfile}.coordsorted.bam', fastqfile=FASTQFILES),
        expand('All_output/Mapped_reads/{fastqfile}.coordsorted.bam.bai', fastqfile=FASTQFILES),
        expand("All_output/Mapped_reads/{fastqfile}.bam", fastqfile=FASTQFILES),
        "Analysis_Results/Trimming/PostTrimming_QC.html",
        expand("All_output/Trimmed_reads/{fastqfile}_Trimmed_R1.fastq", fastqfile=FASTQFILES),
        expand("All_output/Trimmed_reads/{fastqfile}_Trimmed_R2.fastq", fastqfile=FASTQFILES),
        "Analysis_Results/QC_Rawreads/Rawreads_QC.html",
        expand("Analysis_Results/QC_Rawreads/{fastqfile}_{read}_fastqc.html", fastqfile=FASTQFILES, read=READS),
        expand("Analysis_Results/QC_Rawreads/{fastqfile}_{read}_fastqc.zip", fastqfile=FASTQFILES, read=READS),







#Quality Control FastQC
rule QCrawreads_Fastqc:
    input:
        f'{ROOT_DIR}/{{fastqfile}}_{{read}}.fastq.gz'
    output:
        ('Analysis_Results/QC_Rawreads/{fastqfile}_{read}_fastqc.html'),
        ('Analysis_Results/QC_Rawreads/{fastqfile}_{read}_fastqc.zip')
    log:
        'logs/fastqc_rawreads/{fastqfile}_{read}.log'
    shell:
        """
        fastqc {input[0]} --outdir=./Analysis_Results/QC_Rawreads &>> {log}
        """

rule Compileresults_QC:
    input: expand('Analysis_Results/QC_Rawreads/{fastqfile}_{read}_fastqc.html', read=READS, fastqfile=FASTQFILES)
    output:
     html="Analysis_Results/QC_Rawreads/Rawreads_QC.html"
    log:
        'logs/compileresults/QC.log'
    shell:
        "multiqc ./Analysis_Results/QC_Rawreads --force -v -o ./Analysis_Results/QC_Rawreads -n Rawreads_QC.html &>> {log} "




rule AdapterTrim_Cutadapt:
    input:
        f'{ROOT_DIR}/{{fastqfile}}_R1.fastq.gz',
        f'{ROOT_DIR}/{{fastqfile}}_R2.fastq.gz',
        'Analysis_Results/QC_Rawreads/{fastqfile}_R1_fastqc.html',
        'Analysis_Results/QC_Rawreads/{fastqfile}_R1_fastqc.zip',
        'Analysis_Results/QC_Rawreads/{fastqfile}_R2_fastqc.html',
        'Analysis_Results/QC_Rawreads/{fastqfile}_R2_fastqc.zip'
    output:
        Trim_read1='All_output/Trimmed_reads/{fastqfile}_Trimmed_R1.fastq',
        Trim_read2='All_output/Trimmed_reads/{fastqfile}_Trimmed_R2.fastq'
    log:
        'logs/cutadapt/{fastqfile}.log'
    shell:
        'cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT  -o {output.Trim_read1} -p {output.Trim_read2} {input[0]} {input[1]} &>> {log}'


rule Compileresults_PosttrimQC:
    input: expand('logs/cutadapt/{fastqfile}.log', fastqfile=FASTQFILES)
    output:
        html="Analysis_Results/Trimming/PostTrimming_QC.html"
    log:
        'logs/compileresults/PosttrimQC.log'
    shell:
        "multiqc ./logs/cutadapt --force -v -o ./Analysis_Results/Trimming -n PostTrimming_QC.html &>> {log} "


rule Map_Bowtie2:
    input:
        Trim_read1='All_output/Trimmed_reads/{fastqfile}_Trimmed_R1.fastq',
        Trim_read2='All_output/Trimmed_reads/{fastqfile}_Trimmed_R2.fastq'
    output:
        temp('All_output/Mapped_reads/{fastqfile}.bam'),
        sortedBam='All_output/Mapped_reads/{fastqfile}.coordsorted.bam',
        bai='All_output/Mapped_reads/{fastqfile}.coordsorted.bam.bai'
    params:
        index=config['bowtie2_index'],
        align=config['genome_align'],
        PicardCmd=config['PicardCmd']
    log:
        bowtie2="logs/primary_alignment/bowtie2/{fastqfile}.log",
        picard="logs/primary_alignment/picard_sort/{fastqfile}.log"
    threads: 9
    resources:
        mem_mb=8000,
        runtime=1440
    shell:
        """
        bowtie2 {params.align} -x {params.index} \
        -1 {input.Trim_read1} -2 {input.Trim_read2} 2> {log.bowtie2} | samtools view -bS - > {output[0]}

        {params.PicardCmd} SortSam -I {output[0]} -O {output.sortedBam} -SORT_ORDER coordinate &>> {log.picard}

        samtools index {output.sortedBam}

        """


rule Collect_alignment_stats:
    input:
        sortedbam='All_output/Mapped_reads/{fastqfile}.coordsorted.bam'
    output:
        temp('logs/primary_alignment/PostAlignmentStats/dupstats/{fastqfile}.dupMarked.bam'),
        DupStats='logs/primary_alignment/PostAlignmentStats/dupstats/{fastqfile}_picard.dupMark.txt'
    params:
        PicardCmd=config['PicardCmd']
    log:
        dupstats='logs/primary_alignment/PostAlignmentStats/dupstats/{fastqfile}.log',
        flagstat='logs/primary_alignment/PostAlignmentStats/flagstat/{fastqfile}.log'
    threads: 6
    resources:
        mem_mb=4000,
        runtime=1440
    shell:
        """
        {params.PicardCmd} MarkDuplicates -I {input.sortedbam} -O {output[0]} -METRICS_FILE {output.DupStats} &>> {log.dupstats}

        samtools flagstat {input.sortedbam} &>> {log.flagstat}

        """


rule Compileresults_map:
    input:
        expand(['logs/primary_alignment/PostAlignmentStats/dupstats/{fastqfile}.log',
        'logs/primary_alignment/PostAlignmentStats/flagstat/{fastqfile}.log',
        "logs/primary_alignment/bowtie2/{fastqfile}.log",
        "logs/primary_alignment/picard_sort/{fastqfile}.log",
        'logs/primary_alignment/PostAlignmentStats/dupstats/{fastqfile}.dupMarked.bam',
        'logs/primary_alignment/PostAlignmentStats/dupstats/{fastqfile}_picard.dupMark.txt'], fastqfile=FASTQFILES)
    output:
        html="Analysis_Results/primary_alignment/Alignment_results.html"
    log:
        'logs/compileresults/map.log'
    threads: 6
    resources:
        mem_mb=4000,
        runtime=1440
    shell:
        "multiqc ./logs/primary_alignment --force -v -o ./Analysis_Results/primary_alignment -n Alignment_results.html &>> {log} "



rule Filtering_bams_PicardSamtools:
    input:
        sortedbam='All_output/Mapped_reads/{fastqfile}.coordsorted.bam'
    output:
        'All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.bam',
        'All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.bam.bai',
        NoDupsBam='All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.NoDups.bam',
        rmDups='logs/filtered_bams/{fastqfile}_picard.rmDup.txt',
        nodupsBamindex='All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.NoDups.bam.bai',
    log:
        'logs/filtered_bams/{fastqfile}.log'
    params:
        Prop_paired=config['samtools_proper_paired'],
        Mapq10=config['samtools_mapq'],
        PicardCmd=config['PicardCmd']
    threads: 8
    resources:
        mem_mb=4000,
        runtime=1440
    shell:
        """
        samtools view -bu {params.Prop_paired} {input.sortedbam} | samtools view -b {params.Mapq10} - | samtools sort - -o {output[0]} &>> {log}


        {params.PicardCmd} MarkDuplicates -I {output[0]} \
                                   -O {output.NoDupsBam} \
                                   -REMOVE_DUPLICATES true \
                                   -METRICS_FILE {output.rmDups} &>> {log}

        samtools flagstat {output.NoDupsBam} &>> {log}

        samtools index {output.NoDupsBam}

        samtools index {output[0]}


        """

rule Compileresults_filtering:
    input:
        expand('logs/filtered_bams/{fastqfile}.log', fastqfile=FASTQFILES)
    output:
        html="Analysis_Results/primary_alignment/filteringbamsStats.html"
    log:
        'logs/compileresults/filtering.log'
    shell:
        "multiqc ./logs/filtered_bams --force -v -o ./Analysis_Results/primary_alignment -n filteringbamsStats.html &>> {log} "



rule GetBigwigs_BamCoverage:
    input:
        NoDupsBam='All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.NoDups.bam',
        MAPQfiltBam='All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.bam',
        nodupsBamindex='All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.NoDups.bam.bai',
        MAPQfiltBamindex='All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.bam.bai'
    output:
        bigwig_CPM='Analysis_Results/CPM_and_Unnormalized_bws/CPM_normalized_bws/{fastqfile}_CPM.bw',
        bigwig_WOnorm='Analysis_Results/CPM_and_Unnormalized_bws/bws_wo_normalization/{fastqfile}_wo.norm.bw',
        bigwig_WOnorm_wDups='Analysis_Results/CPM_and_Unnormalized_bws/bws_wo_norm_wDups/{fastqfile}_wo.norm_wDups.bw',
    params:
        bamCov_default=config['bamCov_default'],
        bamCov_min=config['bamCov_min'],
        bamCov_CPM=config['bamCov_CPM']
    log:
        'logs/CPM_and_Unnormalized_bws/{fastqfile}.log'
    threads: 4
    resources:
        mem_mb=2000,
        runtime=1440
    shell:
        """

        # CPM
        bamCoverage --bam {input.NoDupsBam} -o {output.bigwig_CPM} {params.bamCov_CPM}  2> {log}

        # No normalization - bigwig
        bamCoverage --bam {input.NoDupsBam} -o {output.bigwig_WOnorm} {params.bamCov_min} 2> {log}

        # No normalization - with Duplicates
        bamCoverage --bam {input.MAPQfiltBam} -o {output.bigwig_WOnorm_wDups} {params.bamCov_min} 2> {log}

        """



rule Map2Spikein_Bowtie2:
    input:
        Trim_read1='All_output/Trimmed_reads/{fastqfile}_Trimmed_R1.fastq',
        Trim_read2='All_output/Trimmed_reads/{fastqfile}_Trimmed_R2.fastq'
    output:
        temp('All_output/Spike_mapped_reads/{fastqfile}.bam'),
        temp('All_output/Spike_mapped_reads/{fastqfile}.coordsorted.bam'),
        temp('All_output/Spike_mapped_reads/{fastqfile}.coordsorted.bam.bai')
    params:
        alignSpike=config['Spike_align'],
        PicardCmd=config['PicardCmd']
    log:
        Spike_bowtie2="logs/Spike_Alignment/bowtie2/{fastqfile}.log",
        Spike_picard="logs/Spike_Alignment/picard_sort/{fastqfile}.log"
    threads: 9
    resources:
        mem_mb=8000,
        runtime=1440
    shell:
        """
        bowtie2 {params.alignSpike} -x {SPIKEINDEX} \
        -1 {input.Trim_read1} -2 {input.Trim_read2} 2> {log.Spike_bowtie2} | samtools view -Sb - > {output[0]}

        {params.PicardCmd} SortSam -I {output[0]} -O {output[1]} -SORT_ORDER coordinate &>> {log.Spike_picard}

        samtools index {output[1]}

        """


rule Collect_Spikealignment_stats:
    input:
        S_sorted_bam='All_output/Spike_mapped_reads/{fastqfile}.coordsorted.bam'
    output:
        S_DupMarkedBam='logs/Spike_Alignment/dupstats/{fastqfile}.dupMarked.bam',
        S_DupStats='logs/Spike_Alignment/dupstats/{fastqfile}_picard.dupMark.txt'
    params:
        PicardCmd=config['PicardCmd']
    log:
        S_dupstats='logs/Spike_Alignment/dupstats/{fastqfile}.log',
        S_flagstat='logs/Spike_Alignment/flagstat/{fastqfile}.log'
    threads: 4
    resources:
        mem_mb=2000,
        runtime=1440
    shell:
        """
        {params.PicardCmd} MarkDuplicates -I {input.S_sorted_bam} -O {output.S_DupMarkedBam} -METRICS_FILE {output.S_DupStats} &>> {log.S_dupstats}

        samtools flagstat {input.S_sorted_bam} &>> {log.S_flagstat}

        """


rule Compileresults_Spike:
    input: expand('logs/Spike_Alignment/flagstat/{fastqfile}.log', fastqfile=FASTQFILES)
    output:
        html="Analysis_Results/Spikein_alignment/Spike_alignment.html"
    log:
        'logs/compileresults/Spikealign.log'
    shell:
        "multiqc ./logs/Spike_Alignment --force -v -o ./Analysis_Results/Spikein_alignment -n Spike_alignment.html &>> {log} "


rule CalcNormFactors:
    input:
        html="Analysis_Results/Spikein_alignment/Spike_alignment.html"
    output:
        ScalingFactors="Analysis_Results/Spikein_normalized_bws_bdgs/ScalingFactors.csv",
        SFbygroup=expand("Analysis_Results/Spikein_normalized_bws_bdgs/ScalingFactors_by_group/{Groupname}_ScalingFactors.csv", Groupname=ABGROUPS)
    log:
        'logs/compileresults/Scalefactors.log'
    script:
        'Scripts/GetScalingFactorsbyGroup.py'


rule GetNormBwsBdgs_BamCoverage_NonInput:
    input:
        NoDupsBam='All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.NoDups.bam',
        MAPQfiltMappedPairedBam='All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.bam',
        ScalingFactors="Analysis_Results/Spikein_normalized_bws_bdgs/ScalingFactors.csv",
    output:
        bigwig_Spikenorm='Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs/{fastqfile}_Norm.bw',
        bdg_Spikenorm='Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/{fastqfile}_Norm.bedgraph',
        bigwig_Spikenorm_wDups='Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs_wDups/{fastqfile}_Norm_wDups.bw',
    params:
        bamCov_default=config['bamCov_default']
    log:
        "logs/Spikein_normalized_bws_bdgs/{fastqfile}.log",
    run:

        AlignStats = pd.read_csv(input.ScalingFactors).rename(columns={'Unnamed: 0': 'Samplename'}).set_index('Samplename')

        if wildcards.fastqfile in ALL_TARGETS:

            logger.info(f'prefix: {wildcards.fastqfile}')

            groupname_ = samples.loc[samples.index == wildcards.fastqfile]['AntibodyGroup'].tolist()
            groupname = groupname_[0]

            logger.info(f'groupname: {groupname}')

            Sfvalue = AlignStats.loc[wildcards.fastqfile, groupname]

            logger.info(f'Sfvalue: {Sfvalue}')

            #with dups
            shell( "bamCoverage --bam {input.MAPQfiltMappedPairedBam} -o {output.bigwig_Spikenorm_wDups} \
                    --scaleFactor {Sfvalue} {params.bamCov_default} &>> {log} " )

            # without dups
            shell( "bamCoverage --bam {input.NoDupsBam} -o {output.bigwig_Spikenorm} \
                    --scaleFactor {Sfvalue} {params.bamCov_default} &>> {log} " )


            shell( "bamCoverage --bam {input.NoDupsBam} --outFileFormat bedgraph -o {output.bdg_Spikenorm} \
                    --scaleFactor {Sfvalue} {params.bamCov_default} &>> {log} " )


if len(ALL_IGGFILES) != 0:

    rule GetNormBwsBdgs_BamCoverage_Input:
        input:
            NoDupsBam=expand('All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.NoDups.bam', fastqfile=ALL_IGGFILES),
            MAPQfiltMappedPairedBam=expand('All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.bam',fastqfile=ALL_IGGFILES),
            ScalingFactors="Analysis_Results/Spikein_normalized_bws_bdgs/ScalingFactors.csv",
        output:
            logf=expand("logs/Spikein_normalized_bws_bdgs/{fastqfile}.log", fastqfile=ALL_IGGFILES),
            mockfile=temp("logs/Spikein_normalized_bws_bdgs/mockfile.txt")
        params:
            bamCov_default=config['bamCov_default'],
        threads: 9
        resources:
            mem_mb=8000,
            runtime=1440
        run:
            AlignStats = pd.read_csv(input.ScalingFactors).rename(columns={'Unnamed: 0': 'Samplename'}).set_index('Samplename')

            for NoDupsfile in input.NoDupsBam:

                key = os.path.basename(NoDupsfile).split(".")[0]
                logger.info(f"key: {key}")

                for group in AlignStats.columns.tolist():

                    outputfilename = f"{group}_{key}"
                    logger.info(f"outputfilename: {outputfilename}")


                    # manually defining input based on key (derived from NodupsBam):
                    MAPQfiltMappedPairedBam_f= [x for x in input.MAPQfiltMappedPairedBam if key in x]
                    log_f = [x for x in output.logf if key in x]

                    # manually defining outputpaths, outputfiles  based on key (derived from NodupsBam):
                    outputdirprefix1="Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs/Inputs"
                    outputdirprefix2="Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/Inputs"
                    outputdirprefix3="Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs_wDups/Inputs"

                    # manually creating directories

                    shell( "mkdir -p {outputdirprefix1}/{key} " )
                    shell( "mkdir -p {outputdirprefix2}/{key} " )
                    shell( "mkdir -p {outputdirprefix3}/{key} " )

                    Sfvalue = AlignStats.loc[key, group]
                    logger.info(f"Sfvalue: {Sfvalue}")


                    if math.isnan(Sfvalue) == False:


                        ## output
                        bigwig_Spikenorm_f=os.path.join(outputdirprefix1, key, f"{outputfilename}.bw")
                        bdg_Spikenorm_f=os.path.join(outputdirprefix2, key, f"{outputfilename}_Norm.bedgraph")
                        bigwig_Spikenorm_wDups_f=os.path.join(outputdirprefix3, key, f"{outputfilename}_Norm_wDups.bw")


                        #with dups
                        shell( "bamCoverage --bam {MAPQfiltMappedPairedBam_f} -o {bigwig_Spikenorm_wDups_f} \
                                --scaleFactor {Sfvalue} {params.bamCov_default} &>> {log_f} " )

                        # without dups
                        shell( "bamCoverage --bam {NoDupsfile} -o {bigwig_Spikenorm_f} \
                                --scaleFactor {Sfvalue} {params.bamCov_default} &>> {log_f} " )

                        shell( "bamCoverage --bam {NoDupsfile} --outFileFormat bedgraph -o {bdg_Spikenorm_f} \
                                --scaleFactor {Sfvalue} {params.bamCov_default} &>> {log_f} " )

                    with open (output.mockfile, "a") as f:
                        f.write("rule GetNormBwsBdgs_BamCoverage_Input completed!!!!!")
                        f.close()


#optionaloutput="Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs/Inputs/{fastqfile}" if len(ALL_IGGFILES) != 0 else []


if len(ALL_IGGFILES) != 0:

    logger.info("input files present...using rule Peaks_SEACR_ver1")

    rule Peaks_SEACR_ver1:
        input:
            Target='Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/{fastqfile}_Norm.bedgraph',
            mockfile="logs/Spikein_normalized_bws_bdgs/mockfile.txt"
        output:
            Strin='Analysis_Results/Peaks/{fastqfile}.stringent.bed',
        params:
            SEACRLoc=config['SEACRLoc'],
        log:
            "logs/Peaks/{fastqfile}.log"
        threads: 9
        resources:
            mem_mb=8000,
            runtime=1440
        run:
            Control = samples.loc[samples.index.str.contains(wildcards.fastqfile)]['IgGFile'][0]
            Groupname = samples.loc[samples.index.str.contains(wildcards.fastqfile)]['AntibodyGroup'][0]

            logger.info(f"Control: {Control}, Groupname: {Groupname}")

            _IgG_=os.path.join("Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/Inputs/",Control, f"{Groupname}_{Control}_Norm.bedgraph")
            logger.info(f"{params.SEACRLoc} {input.Target} {_IgG_} non stringent Analysis_Results/Peaks/{wildcards.fastqfile} &>> {log}")
            shell ( "bash {params.SEACRLoc} {input.Target} {_IgG_} non stringent Analysis_Results/Peaks/{wildcards.fastqfile} &>> {log} " )
            shell ( "bash {params.SEACRLoc} {input.Target} {_IgG_} non relaxed Analysis_Results/Peaks/{wildcards.fastqfile} &>> {log} " )


else:
    logger.info("no input files...using rule Peaks_SEACR_ver2")
    rule Peaks_SEACR_ver2:
        input:
            Target='Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/{fastqfile}_Norm.bedgraph',
        output:
            Strin='Analysis_Results/Peaks/{fastqfile}.stringent.bed',
        params:
            SEACRLoc=config['SEACRLoc'],
        log:
            "logs/Peaks/{fastqfile}.log"
        threads: 9
        resources:
            mem_mb=8000,
            runtime=1440
        run:
            Control = samples.loc[samples.index.str.contains(wildcards.fastqfile)]['IgGFile'][0]
            Groupname = samples.loc[samples.index.str.contains(wildcards.fastqfile)]['AntibodyGroup'][0]

            logger.info(f"Control: {Control}, Groupname: {Groupname}")


            shell ( "bash {params.SEACRLoc} {input.Target} 0.05 non stringent Analysis_Results/Peaks/{wildcards.fastqfile} &>> {log} " )
            shell ( "bash {params.SEACRLoc} {input.Target} 0.05 non relaxed Analysis_Results/Peaks/{wildcards.fastqfile} &>> {log} " )


#rule Clean_up:
    #input:
        #Strin=expand('Analysis_Results/Peaks/{fastqfile}.stringent.bed', fastqfile=ALL_TARGETS)
    #output:
        #"logs/cleanup.log"
    #threads: 1
    #shell:
        #" rm *.out &>> {output} "
