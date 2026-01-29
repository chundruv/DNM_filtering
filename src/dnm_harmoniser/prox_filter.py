import pandas as pd
from tqdm import tqdm
import numpy as np
import sys

x=pd.read_csv(sys.argv[1], sep='\t')
x['END']=x['Start position']
x.loc[x['ALT'].str.len()>1,'END']=x[x['ALT'].str.len()>1]['Start position'] + x[x['ALT'].str.len()>1]['ALT'].str.len() -1
x['proximity_filter']=1
smpls = x['SAMPLE'].unique()
for i in tqdm(range(len(smpls))):
    tmp = x[x['SAMPLE']==smpls[i]].copy()
    tmp = tmp.sort_values('Start position').reset_index(drop=True)

    chroms = tmp['Chromosome'].unique()
    for chrom in chroms:
        chrom_variants = tmp[tmp['Chromosome']==chrom].copy()

        # Keep track of which variants to filter out
        last_kept_pos = -1000  # Initialize to a position far away
        last_kept_end = -1000

        for idx, row in chrom_variants.iterrows():
            current_start = row['Start position']
            current_end = row['END']

            # If this variant is within 20bp of the last kept variant, mark it as 0
            if current_start <= last_kept_end + 20:
                # This variant is too close to the previous kept variant
                x.loc[((x['SAMPLE']==smpls[i]) &
                       (x['Chromosome']==chrom) &
                       (x['Start position']==current_start) &
                       (x['END']==current_end)), 'proximity_filter'] = 0
            else:
                # This variant is far enough, keep it and update the last kept position
                last_kept_pos = current_start
                last_kept_end = current_end
x.to_csv(sys.argv[2], index=False, sep='\t')