#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import pandas as pd


Allsamples = pd.read_csv(snakemake.config['sample_info'], "\t").set_index("samples")

filetoread = "Analysis_Results/Spikein_alignment/Spike_alignment_data/multiqc_bowtie2.txt"

Spike = pd.read_csv(filetoread, sep='\t')
Spike.set_index('Sample', inplace=True)

if len(Spike.index) == len(Allsamples.index):
    Spike["TotalMappedFragments"] = Spike["paired_aligned_multi"] + Spike["paired_aligned_one"]
    Spike.columns = [str(col) + '_spikein' for col in Spike.columns]


    Smin = Spike['TotalMappedFragments_spikein'].min()
    Spike['ScalingFactors'] = Smin / Spike['TotalMappedFragments_spikein']
    Spike.to_csv(snakemake.output.SpikeAlignStats)

else:
    missingfileAll = (set(Allsamples.index)).difference(set(Spike.index))
    print(f'{missingfileAll} is/are missing')
    raise Exception("Sorry, data for some files is missing.....")
