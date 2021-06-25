import pandas as pd
import numpy as np
import glob
import os
import fnmatch
import yaml
from pathlib import Path


from snakemake.io import expand, glob_wildcards
from snakemake.utils import min_version
from snakemake.logging import logger


# Minimum snakemake version
min_version("5.24.1")

# configuration file
configfile:"./config/config.yaml"

# Tabular configuration
samples = pd.read_csv(config["sample_info"], "\t").set_index("samples")


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




# Getting directory containing rawreads
wdir = os.getcwd()
for filepath, dirs, allfiles in os.walk(wdir):
    for file_ in allfiles:
        if file_.endswith(".fastq.gz"):
            READS_DIR = filepath
            ROOT_DIR = os.path.relpath(filepath, wdir)


logger.info(f'This is the RawReads dir: {READS_DIR}')
logger.info(f'This is the RawReads dir: {ROOT_DIR}')


# Finding md5
def findmd5(ROOT_DIR):
    try:
        for _file_ in os.listdir(READS_DIR):
            if fnmatch.fnmatch(_file_, '*md5*.txt'):
                PATH2MD5 = os.path.join(READS_DIR, _file_)
                logger.info(f'md5sum file: {PATH2MD5}')
                return PATH2MD5

    except:
        logger.info(f'Can not locate md5sum.txt, make sure it is included in the raw-reads directory ')



# Describing wildcards
SINGLE_READ = f'{READS_DIR}/{{fastqfile}}_{{read}}.fastq.gz'


READS = set(glob_wildcards(SINGLE_READ).read)
FASTQFILES = set(glob_wildcards(SINGLE_READ).fastqfile)

# Control vs non-Control files
IGGREADS = set( (samples.loc[samples['condition'].str.contains('IgG', case=False)].index).to_list() )
TARGETS = set( (samples.loc[~samples['condition'].str.contains('IgG', case=False)].index).to_list() )


logger.info(f'Sample file format: {SINGLE_READ}')
logger.info(f'These are the sample names: {FASTQFILES}')
logger.info(f'Each sample has {READS}')
logger.info(f'IgG control: {IGGREADS}')
logger.info(f'Target files: {TARGETS}')


EGS_GRCh38 = {'50': '2308125349', '75': '2747877777', '100': '2805636331', '150': '2862010578', '200': '2887553303'}
EGS_GRCm38 = {'50': '2308125349', '75': '2407883318', '100': '2467481108', '150': '2494787188', '200': '2520869189'}


# Species

READLENGHT = config['read_lenght']

if config['Species'] == 'Mus musculus':
    EFFECTIVEGENOMESIZE = EGS_GRCm38[READLENGHT]


if config['Species'] == 'Homo sapiens':
    EFFECTIVEGENOMESIZE = EGS_GRCh38[READLENGHT]


# Spike
if config['Spikein'] == 'Amp':
    SPIKEINDEX = config['Spikein_index_amp']

if config['Spikein'] == 'Bacteria':
    SPIKEINDEX = config['Spikein_index_Ecoli']


localrules: Clean_up

rule all:
    input:
        expand("Analysis_Results/Peaks/{fastqfile}.stringent.bed",fastqfile=TARGETS),
        "Analysis_Results/Spikein_alignment/Spike_alignment.html",
        "Analysis_Results/Spikein_normalized_bws_bdgs/Spike_align_stats.csv",
        expand("Analysis_Results/Spikein_normalized_bws_bdgs/{fastqfile}_Norm.bedgraph", fastqfile=FASTQFILES),
        expand("Analysis_Results/Spikein_normalized_bws_bdgs/{fastqfile}_Norm_wDups.bw", fastqfile=FASTQFILES),
        expand("Analysis_Results/Spikein_normalized_bws_bdgs/{fastqfile}_Norm.bw", fastqfile=FASTQFILES),
        expand('All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.NoDups.bam.bai', fastqfile=FASTQFILES),
        expand('All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.bam.bai', fastqfile=FASTQFILES),
        expand('Analysis_Results/RPGC_and_Unnormalized_bws/{fastqfile}_RPGC.bw', fastqfile=FASTQFILES),
        expand('Analysis_Results/RPGC_and_Unnormalized_bws/{fastqfile}_wo.norm.bw', fastqfile=FASTQFILES),
        expand('Analysis_Results/RPGC_and_Unnormalized_bws/{fastqfile}_wo.norm_wDups.bw', fastqfile=FASTQFILES),
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






rule check_md5:
    input:
        expand(f'{ROOT_DIR}/{{fastqfile}}_{{read}}.fastq.gz', fastqfile=FASTQFILES, read=READS),
        findmd5
    output:
        directory('logs/md5check')
    log:
        'logs/md5check/md5checks.log'
    run:
        shell( "md5sum --check {input[1]} &>> {log} " )

        shell( " python ./Scripts/md5checks.py " )


#Quality Control FastQC
rule QCrawreads_Fastqc:
    input:
        f'{ROOT_DIR}/{{fastqfile}}_{{read}}.fastq.gz',
        'logs/md5check/md5checks.log'
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
        picardloc=config['PicardLoc']
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

        {params.picardloc} SortSam -I {output[0]} -O {output.sortedBam} -SORT_ORDER coordinate &>> {log.picard}

        samtools index {output.sortedBam}

        """


rule Collect_alignment_stats:
    input:
        sortedbam='All_output/Mapped_reads/{fastqfile}.coordsorted.bam'
    output:
        temp('logs/primary_alignment/PostAlignmentStats/dupstats/{fastqfile}.dupMarked.bam'),
        DupStats='logs/primary_alignment/PostAlignmentStats/dupstats/{fastqfile}_picard.dupMark.txt'
    params:
        picardloc=config['PicardLoc']
    log:
        dupstats='logs/primary_alignment/PostAlignmentStats/dupstats/{fastqfile}.log',
        flagstat='logs/primary_alignment/PostAlignmentStats/flagstat/{fastqfile}.log'
    threads: 6
    resources:
        mem_mb=4000,
        runtime=1440
    shell:
        """
        {params.picardloc} MarkDuplicates -I {input.sortedbam} -O {output[0]} -METRICS_FILE {output.DupStats} &>> {log.dupstats}

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
        picardloc=config['PicardLoc']
    threads: 8
    resources:
        mem_mb=4000,
        runtime=1440
    shell:
        """
        samtools view -bu {params.Prop_paired} {input.sortedbam} | samtools view -b {params.Mapq10} - | samtools sort - -o {output[0]} &>> {log}


        {params.picardloc} MarkDuplicates -I {output[0]} \
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
        bigwig_RPGC='Analysis_Results/RPGC_and_Unnormalized_bws/{fastqfile}_RPGC.bw',
        bigwig_WOnorm='Analysis_Results/RPGC_and_Unnormalized_bws/{fastqfile}_wo.norm.bw',
        bigwig_WOnorm_wDups='Analysis_Results/RPGC_and_Unnormalized_bws/{fastqfile}_wo.norm_wDups.bw',
    params:
        bamCov_default=config['bamCov_default'],
        bamCov_min=config['bamCov_min'],
        bamCov_RPGC=config['bamCov_RPGC']
    log:
        'logs/RPGC_and_Unnormalized_bws/{fastqfile}.log'
    threads: 4
    resources:
        mem_mb=2000,
        runtime=1440
    shell:
        """

        # RPGC
        bamCoverage --bam {input.NoDupsBam} -o {output.bigwig_RPGC} {params.bamCov_RPGC} --effectiveGenomeSize {EFFECTIVEGENOMESIZE} 2> {log}

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
        picardloc=config['PicardLoc']
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

        {params.picardloc} SortSam -I {output[0]} -O {output[1]} -SORT_ORDER coordinate &>> {log.Spike_picard}

        samtools index {output[1]}

        """


rule Collect_Spikealignment_stats:
    input:
        S_sorted_bam='All_output/Spike_mapped_reads/{fastqfile}.coordsorted.bam'
    output:
        S_DupMarkedBam='logs/Spike_Alignment/dupstats/{fastqfile}.dupMarked.bam',
        S_DupStats='logs/Spike_Alignment/dupstats/{fastqfile}_picard.dupMark.txt'
    params:
        picardloc=config['PicardLoc']
    log:
        S_dupstats='logs/Spike_Alignment/dupstats/{fastqfile}.log',
        S_flagstat='logs/Spike_Alignment/flagstat/{fastqfile}.log'
    threads: 4
    resources:
        mem_mb=2000,
        runtime=1440
    shell:
        """
        {params.picardloc} MarkDuplicates -I {input.S_sorted_bam} -O {output.S_DupMarkedBam} -METRICS_FILE {output.S_DupStats} &>> {log.S_dupstats}

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
        SpikeAlignStats="Analysis_Results/Spikein_normalized_bws_bdgs/Spike_align_stats.csv"
    log:
        'logs/compileresults/Scalefacs.log'
    script:
        'Scripts/GetScalingFactors.py'


rule GetNormBwsBdgs_BamCoverage:
    input:
        NoDupsBam='All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.NoDups.bam',
        MAPQfiltMappedPairedBam='All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.bam',
        SpikeAlignStats="Analysis_Results/Spikein_normalized_bws_bdgs/Spike_align_stats.csv"
    output:
        bigwig_Spikenorm='Analysis_Results/Spikein_normalized_bws_bdgs/{fastqfile}_Norm.bw',
        bdg_Spikenorm='Analysis_Results/Spikein_normalized_bws_bdgs/{fastqfile}_Norm.bedgraph',
        bigwig_Spikenorm_wDups='Analysis_Results/Spikein_normalized_bws_bdgs/{fastqfile}_Norm_wDups.bw'
    params:
        prefix=lambda wildcards, output: output[0][45:-8],
        bamCov_default=config['bamCov_default']
    log:
        "logs/Spikein_normalized_bws_bdgs/{fastqfile}.log"
    threads: 6
    resources:
        mem_mb=2000,
        runtime=1440
    run:

        import pandas as pd

        AlignStats = pd.read_csv(input.SpikeAlignStats)
        AlignStats.set_index('Sample', inplace=True)

        Sfvalue = AlignStats.loc[params.prefix, 'ScalingFactors']

        print (Sfvalue)

                #with dups
        shell( "bamCoverage --bam {input.MAPQfiltMappedPairedBam} -o {output.bigwig_Spikenorm_wDups} \
                --scaleFactor {Sfvalue} {params.bamCov_default} 2> {log} " )

        # without dups
        shell( "bamCoverage --bam {input.NoDupsBam} -o {output.bigwig_Spikenorm} \
                --scaleFactor {Sfvalue} {params.bamCov_default} 2> {log} " )


        shell( "bamCoverage --bam {input.NoDupsBam} --outFileFormat bedgraph -o {output.bdg_Spikenorm} \
                --scaleFactor {Sfvalue} {params.bamCov_default} 2> {log} " )


rule Peaks_SEACR:
    input:
        Control=expand( 'Analysis_Results/Spikein_normalized_bws_bdgs/{fastqfile}_Norm.bedgraph', fastqfile=IGGREADS),
        Target=expand('Analysis_Results/Spikein_normalized_bws_bdgs/{fastqfile}_Norm.bedgraph', fastqfile=TARGETS)
    output:
        Strin=expand( 'Analysis_Results/Peaks/{fastqfile}.stringent.bed', fastqfile=TARGETS)
    params:
        SEACRLoc=config['SEACRLoc']
    log:
        expand( "logs/Peaks/{fastqfile}.log", fastqfile=TARGETS)
    threads: 4
    resources:
        mem_mb=2000,
        runtime=1440
    run:
        if len(input.Control) >= 1:
            for IgG in input.Control:
                for Targetfile in input.Target:
                    Name = (os.path.basename(Targetfile).split(".")[0]).replace("_Norm", "")
                    shell ( "bash {params.SEACRLoc} {Targetfile} {IgG} non stringent Analysis_Results/Peaks/{Name} &>> {log} " )
                    shell ( "bash {params.SEACRLoc} {Targetfile} {IgG} non relaxed Analysis_Results/Peaks/{Name} &>> {log} " )

        else:
            for Targetfile in input.Target:
                Name = (os.path.basename(Targetfile).split(".")[0]).replace("_Norm", "")
                shell ( "bash {params.SEACRLoc} {Targetfile} 0.05 non stringent Analysis_Results/Peaks/{Name} &>> {log} " )
                shell ( "bash {params.SEACRLoc} {Targetfile} 0.05 non relaxed Analysis_Results/Peaks/{Name} &>> {log} " )

rule Clean_up:
    input:
        Strin=expand( 'Analysis_Results/Peaks/{fastqfile}.stringent.bed', fastqfile=TARGETS)
    output:
        "logs/cleanup.log"
    shell:
        " rm *.out &>> {output} "
