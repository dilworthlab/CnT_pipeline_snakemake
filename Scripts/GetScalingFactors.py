#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import pandas as pd
import yaml

filetoread = "Analysis_Results/Spikein_alignment/Spike_alignment_data/multiqc_bowtie2.txt"

Spike = pd.read_csv(filetoread, sep='\t')
Spike.set_index('Sample', inplace=True)
Spike["TotalMappedFragments"] = Spike["paired_aligned_multi"] + Spike["paired_aligned_one"]
Spike.columns = [str(col) + '_spikein' for col in Spike.columns]

Smin = Spike['TotalMappedFragments_spikein'].min()
Spike['ScalingFactors'] = Smin / Spike['TotalMappedFragments_spikein']
Spike.to_csv(snakemake.output.SpikeAlignStats)










