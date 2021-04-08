import pandas as pd
import numpy as np
import glob
import os


from snakemake.io import expand, glob_wildcards
from snakemake.utils import min_version
from snakemake.logging import logger


# Minimum snakemake version
min_version("5.24.1")

# configuration file
configfile:"./config/config.yaml"

# Tabular configuration
samples = pd.read_csv(config["sample_info"], "\t").set_index("samples")


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

# Control vs non-Control files
IGGREADS = set( (samples.loc[samples['condition'].str.contains('IgG', case=False)].index).to_list() )
TARGETS = set( (samples.loc[~samples['condition'].str.contains('IgG', case=False)].index).to_list() )


logger.info(f'Sample file format: {SINGLE_READ}')
logger.info(f'These are the sample names: {FASTQFILES}')
logger.info(f'Each sample has {READS}')
logger.info(f'IgG control: {IGGREADS}')
logger.info(f'Target files: {TARGETS}')

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
        expand('Analysis_Results/RPGC_and_Unnormalized_bws/{fastqfile}_wo.norm.bedgraph', fastqfile=FASTQFILES),
        expand('All_output/Spike_mapped_reads/{fastqfile}.coordsorted.bam',fastqfile=FASTQFILES),
        expand('All_output/Spike_mapped_reads/{fastqfile}.bam', fastqfile=FASTQFILES),
        expand('All_output/Spike_mapped_reads/{fastqfile}.coordsorted.bam.bai', fastqfile=FASTQFILES),
        expand('logs/Spike_Alignment/dupstats/{fastqfile}.dupMarked.bam', fastqfile=FASTQFILES),
        expand('logs/Spike_Alignment/dupstats/{fastqfile}_picard.dupMark.txt', fastqfile=FASTQFILES),
        "Analysis_Results/primary_alignment/filteringbamsStats.html",
        expand('logs/filtered_bams/{fastqfile}_picard.rmDup.txt', fastqfile=FASTQFILES),
        expand("All_output/Processed_reads/{fastqfile}.MappedPaired.bam", fastqfile=FASTQFILES),
        expand("All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.bam", fastqfile=FASTQFILES),
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
        f'{ROOT_DIR}/hashes'



# Check md5sum
rule check_md5:
    input:
        expand(f'{ROOT_DIR}/{{fastqfile}}_{{read}}.fastq.gz', fastqfile=FASTQFILES, read=READS)
    output:
        f'{ROOT_DIR}/hashes'
    log:
        'logs/md5checks.log'
    run:
        shell( "md5sum {ROOT_DIR}/*.fastq.gz > {ROOT_DIR}/hashes " )
        shell( "md5sum --check {ROOT_DIR}/hashes &>> {log} " )

        shell( " python ./Scripts/md5checks.py " )


#Quality Control FastQC
rule QCrawreads_fastqc:
    input:
        f'{ROOT_DIR}/{{fastqfile}}_{{read}}.fastq.gz',
        f'{ROOT_DIR}/hashes'
    output:
        ('Analysis_Results/QC_Rawreads/{fastqfile}_{read}_fastqc.html'),
        ('Analysis_Results/QC_Rawreads/{fastqfile}_{read}_fastqc.zip')
    log:
        'logs/fastqc_rawreads/{fastqfile}_{read}.log'
    conda:
        'envs_conda/fastqc.yaml'
    shell:
        """
        fastqc {input[0]} --outdir=./Analysis_Results/QC_Rawreads &>> {log}
        """

rule Compileresults_QC:
    input: expand('Analysis_Results/QC_Rawreads/{fastqfile}_{read}_fastqc.html', read=READS, fastqfile=FASTQFILES)
    output:
     html="Analysis_Results/QC_Rawreads/Rawreads_QC.html"
    conda:
        'envs_conda/compile_results.yaml'
    shell:
        "multiqc ./Analysis_Results/QC_Rawreads -o ./Analysis_Results/QC_Rawreads -n Rawreads_QC.html "




rule AdapterTrim_cutadapt:
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
    conda:
        'envs_conda/cutadapt.yaml'
    shell:
        'cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT  -o {output.Trim_read1} -p {output.Trim_read2} {input[0]} {input[1]} &>> {log}'


rule Compileresults_PosttrimQC:
    input: expand('logs/cutadapt/{fastqfile}.log', fastqfile=FASTQFILES)
    output:
        html="Analysis_Results/Trimming/PostTrimming_QC.html"
    conda:
        'envs_conda/compile_results.yaml'
    shell:
        "multiqc ./logs/cutadapt -o ./Analysis_Results/Trimming -n PostTrimming_QC.html "


rule Map_Bowtie2:
    input:
        Trim_read1='All_output/Trimmed_reads/{fastqfile}_Trimmed_R1.fastq',
        Trim_read2='All_output/Trimmed_reads/{fastqfile}_Trimmed_R2.fastq'
    output:
        'All_output/Mapped_reads/{fastqfile}.bam',
        sortedBam='All_output/Mapped_reads/{fastqfile}.coordsorted.bam',
        bai='All_output/Mapped_reads/{fastqfile}.coordsorted.bam.bai'
    params:
        genome=config['genome_index'],
        align=config['genome_align'],
        picardloc=config['PicardLoc']
    log:
        bowtie2="logs/primary_alignment/bowtie2/{fastqfile}.log",
        picard="logs/primary_alignment/picard_sort/{fastqfile}.log"
    conda:
        'envs_conda/bowtie2.yaml'
    threads: 9
    resources:
        mem_mb=8000,
        runtime=1440
    shell:
        """
        bowtie2 {params.align} -x {params.genome} \
        -1 {input.Trim_read1} -2 {input.Trim_read2} 2> {log.bowtie2} | samtools view -bS - > {output[0]}

        {params.picardloc} SortSam -I {output[0]} -O {output.sortedBam} -SORT_ORDER coordinate &>> {log.picard}

        samtools index {output.sortedBam}

        """


rule Collect_alignment_stats:
    input:
        sortedbam='All_output/Mapped_reads/{fastqfile}.coordsorted.bam'
    output:
        ('logs/primary_alignment/PostAlignmentStats/dupstats/{fastqfile}.dupMarked.bam'),
        DupStats='logs/primary_alignment/PostAlignmentStats/dupstats/{fastqfile}_picard.dupMark.txt'
    params:
        picardloc=config['PicardLoc']
    log:
        dupstats='logs/primary_alignment/PostAlignmentStats/dupstats/{fastqfile}.log',
        flagstat='logs/primary_alignment/PostAlignmentStats/flagstat/{fastqfile}.log'
    threads: 4
    resources:
        mem_mb=2000,
        runtime=1440
    shell:
        """
        {params.picardloc} MarkDuplicates -I {input.sortedbam} -O {output[0]} -METRICS_FILE {output.DupStats} &>> {log.dupstats}

        samtools flagstat {input.sortedbam} &>> {log.flagstat}

        """


rule Compileresults_map:
    input:
        expand(['logs/primary_alignment/PostAlignmentStats/dupstats/{fastqfile}.log' ,
        'logs/primary_alignment/PostAlignmentStats/flagstat/{fastqfile}.log',
        "logs/primary_alignment/bowtie2/{fastqfile}.log",
        "logs/primary_alignment/picard_sort/{fastqfile}.log",
        'logs/primary_alignment/PostAlignmentStats/dupstats/{fastqfile}.dupMarked.bam',
        'logs/primary_alignment/PostAlignmentStats/dupstats/{fastqfile}_picard.dupMark.txt'], fastqfile=FASTQFILES)
    output:
        html="Analysis_Results/primary_alignment/Alignment_results.html"
    threads: 4
    resources:
        mem_mb=2000,
        runtime=1440
    shell:
        "multiqc ./logs/primary_alignment -o ./Analysis_Results/primary_alignment -n Alignment_results.html  "



rule filtering_bams:
    input:
        sortedbam='All_output/Mapped_reads/{fastqfile}.coordsorted.bam'
    output:
        MappedPairedBam='All_output/Processed_reads/{fastqfile}.MappedPaired.bam',
        MAPQfiltBam='All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.bam',
        NoDupsBam='All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.NoDups.bam',
        rmDups='logs/filtered_bams/{fastqfile}_picard.rmDup.txt'
    log:
        'logs/filtered_bams/{fastqfile}.log'
    params:
        Prop_paired=config['samtools_proper_paired'],
        Mapq10=config['samtools_mapq'],
        picardloc=config['PicardLoc']
    threads: 4
    resources:
        mem_mb=2000,
        runtime=1440
    shell:
        """
        samtools view -bu {params.Prop_paired} {input.sortedbam} | samtools sort - -o {output.MappedPairedBam} &>> {log}

        samtools view -b {params.Mapq10} {output.MappedPairedBam} | samtools sort - -o {output.MAPQfiltBam} &>> {log}

        {params.picardloc} MarkDuplicates -I {output.MAPQfiltBam} \
                                   -O {output.NoDupsBam} \
                                   -REMOVE_DUPLICATES true \
                                   -METRICS_FILE {output.rmDups} &>> {log}

        samtools flagstat {output.NoDupsBam} &>> {log}

        """

rule Compileresults_filtering:
    input:
        expand('logs/filtered_bams/{fastqfile}.log', fastqfile=FASTQFILES)
    output:
        html="Analysis_Results/primary_alignment/filteringbamsStats.html"
    shell:
        "multiqc ./logs/filtered_bams -o ./Analysis_Results/primary_alignment -n filteringbamsStats.html "


rule BamCoverage_GetBigwigs:
    input:
        NoDupsBam='All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.NoDups.bam',
        MAPQfiltBam='All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.bam'
    output:
        nodupsBamindex='All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.NoDups.bam.bai',
        MAPQfiltBamindex='All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.bam.bai',
        bigwig_RPGC='Analysis_Results/RPGC_and_Unnormalized_bws/{fastqfile}_RPGC.bw',
        bigwig_WOnorm='Analysis_Results/RPGC_and_Unnormalized_bws/{fastqfile}_wo.norm.bw',
        bigwig_WOnorm_wDups='Analysis_Results/RPGC_and_Unnormalized_bws/{fastqfile}_wo.norm_wDups.bw',
        bdg_WOnorm='Analysis_Results/RPGC_and_Unnormalized_bws/{fastqfile}_wo.norm.bedgraph'
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
        samtools index {input.NoDupsBam}

        samtools index {input.MAPQfiltBam}

        # RPGC
        bamCoverage --bam {input.NoDupsBam} -o {output.bigwig_RPGC} {params.bamCov_RPGC} 2> {log}

        # No normalization - bigwig
        bamCoverage --bam {input.NoDupsBam} -o {output.bigwig_WOnorm} {params.bamCov_min} 2> {log}

        # No normalization - bedgraph
        bamCoverage --bam {input.NoDupsBam} --outFileFormat bedgraph -o {output.bdg_WOnorm} {params.bamCov_default} 2> {log}

        # No normalization - with Duplicates
        bamCoverage --bam {input.MAPQfiltBam} -o {output.bigwig_WOnorm_wDups} {params.bamCov_min} 2> {log}

        """



rule Map2Spikein_Bowtie2:
    input:
        Trim_read1='All_output/Trimmed_reads/{fastqfile}_Trimmed_R1.fastq',
        Trim_read2='All_output/Trimmed_reads/{fastqfile}_Trimmed_R2.fastq'
    output:
        S_bam='All_output/Spike_mapped_reads/{fastqfile}.bam',
        S_sortedBam='All_output/Spike_mapped_reads/{fastqfile}.coordsorted.bam',
        S_bai='All_output/Spike_mapped_reads/{fastqfile}.coordsorted.bam.bai'
    params:
        SpikeGenome=config['Spikein_index'],
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
        bowtie2 {params.alignSpike} -x {params.SpikeGenome} \
        -1 {input.Trim_read1} -2 {input.Trim_read2} 2> {log.Spike_bowtie2} | samtools view -Sb - > {output.S_bam}

        {params.picardloc} SortSam -I {output.S_bam} -O {output.S_sortedBam} -SORT_ORDER coordinate &>> {log.Spike_picard}

        samtools index {output.S_sortedBam}

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
    shell:
        "multiqc ./logs/Spike_Alignment -o ./Analysis_Results/Spikein_alignment -n Spike_alignment.html "


rule CalcNormFactors:
    input:
        html="Analysis_Results/Spikein_alignment/Spike_alignment.html"
    output:
        SpikeAlignStats="Analysis_Results/Spikein_normalized_bws_bdgs/Spike_align_stats.csv"
    script:
        'Scripts/GetScalingFactors.py'


rule GetNormBwsBdgs_BamCoverage:
    input:
        NoDupsBam='All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.NoDups.bam',
        MAPQfiltBam='All_output/Processed_reads/{fastqfile}.MappedPaired.MAPQ10.bam',
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
    threads: 4
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
        shell( "bamCoverage --bam {input.MAPQfiltBam} -o {output.bigwig_Spikenorm_wDups} \
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
