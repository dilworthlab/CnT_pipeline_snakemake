#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import pandas as pd


samples = pd.read_csv(snakemake.config['sample_info'], "\t").set_index("TargetFiles")

SpikeAlignmentStats_f_ = "Analysis_Results/Spikein_alignment/Spike_alignment_data/multiqc_bowtie2.txt"

SpikeAlignmentStats = pd.read_csv(SpikeAlignmentStats_f_, sep='\t')
SpikeAlignmentStats = SpikeAlignmentStats.set_index('Sample')


Groups = samples['AntibodyGroup'].unique().tolist()


DictsofSamps = {}

# Seperate targetfiles into groups
for AB_group in Groups:
    Targetfiles = list(set (samples.loc[samples['AntibodyGroup'] == AB_group ].index))
    Iggfiles = list(set (samples.loc[samples['AntibodyGroup'] == AB_group ]['IgGFile']))

    
    if ("none" in Iggfiles[0].lower()) or len(Iggfiles) == 0 :

        testlist = []
        for target in Targetfiles:
            f = (SpikeAlignmentStats.loc[SpikeAlignmentStats.index.str.contains(target, case=False, regex=False)].index).tolist()
            testlist.append(f[0])

        

            
        group_df = SpikeAlignmentStats.loc[SpikeAlignmentStats.index.isin(testlist)]
        
        group_df = group_df.copy()
        
        group_df["TotalMappedFragments"] = group_df["paired_aligned_multi"] + group_df["paired_aligned_one"]
        group_df.columns = [str(col) + '_spikein' for col in group_df.columns]


        Smin = group_df['TotalMappedFragments_spikein'].min()
        group_df['ScalingFactors'] = Smin / group_df['TotalMappedFragments_spikein']
        
        NameofFile = f"{AB_group}_ScalingFactors.csv"
        
        group_df.to_csv(os.path.join("Analysis_Results/Spikein_normalized_bws_bdgs/ScalingFactors_by_group", NameofFile))
        
        DictsofSamps.update(dict_df)
            
    else:
        testlist = []
        f1 = (SpikeAlignmentStats.loc[SpikeAlignmentStats.index.str.contains(Iggfiles[0], case=False, regex=False)].index).tolist()
        testlist.append(f1[0])

        
        for target in Targetfiles:
            f = (SpikeAlignmentStats.loc[SpikeAlignmentStats.index.str.contains(target, case=False, regex=False)].index).tolist()
            testlist.append(f[0])
        
        group_df = SpikeAlignmentStats.loc[SpikeAlignmentStats.index.isin(testlist)] 
        
        group_df = group_df.copy()
        
        group_df["TotalMappedFragments"] = group_df["paired_aligned_multi"] + group_df["paired_aligned_one"]
        group_df.columns = [str(col) + '_spikein' for col in group_df.columns]


        Smin = group_df['TotalMappedFragments_spikein'].min()
        group_df['ScalingFactors'] = Smin / group_df['TotalMappedFragments_spikein']
        
        NameofFile = f"{AB_group}_ScalingFactors.csv"
        
        group_df.to_csv(os.path.join("Analysis_Results/Spikein_normalized_bws_bdgs/ScalingFactors_by_group", NameofFile))
        
        DictsofSamps.update(dict_df)

DfofSFs = pd.DataFrame(DictsofSamps.items(), columns=['Sample', 'ScalingFactors'])
DfofSFs = DfofSFs.set_index('Sample')
DfofSFs.to_csv(snakemake.output.ScalingFactors)

