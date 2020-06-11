#!/usr/bin/env python
# coding: utf-8

# In[115]:


######## ::: Scan input fasta and prepare data (promoter only, not full genome) ::: ########


# In[116]:


#! conda install biopython


# In[117]:


import re
import os
import sys
import math
import pickle
import pybedtools
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import tensorflow as tf
import multiprocessing as mp
from tensorflow.python.util import deprecation
deprecation._PRINT_DEPRECATION_WARNINGS = False # Supress tensorflow warning
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' # Supress tensorflow warning


dir_path = os.path.dirname(os.path.realpath('prepare_data_prom.ipynb'))

def clean_seq(s):
    ns = s.upper()    
    pattern = re.compile(r'\s+')
    ns = re.sub(pattern, '', ns)
    ns = re.sub(r'[^a-zA-Z]{1}', 'N', ns)
    return ns

## read in fasta
fasta_dir = dir_path + "/../data/promoter/human_epdnew_hg38_noTATA.fa"
fasta = SeqIO.to_dict(SeqIO.parse(open(fasta_dir),'fasta'))

# convert seq to str and clean seq
fasta = {k+' '+v.description.split(' ')[1]: clean_seq(str(v.seq)) for k,v in fasta.items()}

def scan(idx, seq, half_size, key, strand, ref):
    subseq = seq[idx - half_size: idx + half_size]
    
    # make index
    if strand == "-":
        idx1 = len(seq) - (idx + half_size) + 1
        idx2 = len(seq) - (idx - half_size)
    else:
        idx1 = idx - half_size + 1
        idx2 = idx + half_size
    
    chrom_num = "UNKNOWN"
    ref_idx = 0
    if ref is not None: # reference df in bed format, creating index
        gene = key.split(' ')[1]
        ref_idx = ref[ref['name'] == gene]['start'].item() - 5000
        chrom_num = ref[ref['name'] == gene]['chrom'].item()
        
    ck = "{}:{}-{} {} {}".format(chrom_num, idx1+ref_idx, idx2+ref_idx, key, strand)
    
    if not 'N' in subseq:
        return (ck, subseq)
    else:
        return None

    
# read in true TSS bed file
TSS_bed = pybedtools.BedTool(dir_path + "/../data/promoter/human_epdnew_hg38_noTATA.bed")

TSS_bed_df = TSS_bed.to_dataframe()

# multiprocessing
print("Number of processors: ", mp.cpu_count())
pool = mp.Pool(mp.cpu_count())

# globals
scan_step = 100
half_size = 250
batch_size = 128
start = half_size
cutoff = 0.5
proms = {}
scores = {}

i = 1
for key in fasta.keys():
#     print("* Processing " + key)
    print(i, end = " ", flush = True)
    for strand in ["+", "-"]:
#         print("  {} strand".format("Positive" if strand == "+" else "Negative"))
        seq = fasta[key]
        if strand == "-": 
            seq = seq[::-1]
#         print("\t - Scanning...", end = " ")
        putative = pool.starmap(scan, [(i, seq, half_size, key, strand, TSS_bed_df) for i in range(start, len(seq) - half_size, scan_step)]) # multiprocessing
        putative = list(filter(None,putative)) # filter subseq with N
        proms.update({k: v for (k,v) in putative}) # store in dict
        coords = [k for (k,v) in putative] # store coordinates separately
        putative = [v for (k,v) in putative] # use only seq to predict
        probs = []
#         print("Done! Number of subsequences: {}".format(len(putative)))
    i+=1

        
pool.close()    

# create bed file of the scans
bed = []
for i, x in enumerate(list(proms)):
    delim = re.split('[^a-zA-Z0-9+-]',x)
    coords = re.split('-',delim[1])
    bed.append("{}\t{}\t{}\t{}\t1\t{}".format(delim[0],coords[0],coords[1],delim[2]+'_'+delim[3]+'_'+delim[4]+'_'+str(i),delim[5]))
with open(dir_path + "/../data/promoter/human_epdnew_hg38_noTATA_scan_500.bed", 'w+') as f:
    f.write('\n'.join(bed))

    
##### pybedtools #####
# read in true TSS bed file, and make true promoter file -500 to 500 
TSS_bed = pybedtools.BedTool(dir_path + "/../data/promoter/human_epdnew_hg38.bed")
true_prom_bed = TSS_bed.slop(g=dir_path + "/../data/genome/chrom.sizes",l=500,r=499) # -500 to 500

# read in scaned promoter bed files
prom_bed = pybedtools.BedTool(dir_path + "/../data/promoter/human_epdnew_hg38_noTATA_scan_500.bed")

# intersect to get label
# True
overlap_T = prom_bed.intersect(true_prom_bed,
                                s=True, # Require same strandedness.  That is, only report hits in that overlap A on the **same** strand.
                                f=1.0, # Minimum overlap required as a fraction of A
                                u=True, # Write the original A entry **once** if **any** overlaps found in B
#                                 r=True # Require that the fraction overlap be reciprocal for A AND B
                                )

# False
overlap_N = prom_bed.intersect(true_prom_bed, 
                                v=True, # Only report those entries in A that have **no overlaps** with B.
                                s=True, # Require same strandedness.  That is, only report hits in that overlap A on the **same** strand.
                                f=1.0, # Minimum overlap required as a fraction of A
#                                 r=True # Require that the fraction overlap be reciprocal for A AND B
                                )

assert len(prom_bed) == len(overlap_T) + len(overlap_N)

# convert to dataframe, annotate sequences and labels
overlap_T_df = overlap_T.to_dataframe()
overlap_T_df.drop('score', axis=1, inplace=True)
overlap_N_df = overlap_N.to_dataframe()
overlap_N_df.drop('score', axis=1, inplace=True)

# create label
overlap_T_df['label'] = 1
overlap_N_df['label'] = 0

# concatenate
overlap_df = pd.concat([overlap_T_df,overlap_N_df])


# add column of keys in overlap_df, to merge with proms (to get sequences)
overlap_df[['name_1','name_2','name_3','name_4']] = overlap_df['name'].str.split("_", expand = True)
overlap_df['keys'] = overlap_df['chrom']+':'+overlap_df['start'].apply(str)+"-"+overlap_df['end'].apply(str)+" "+ overlap_df['name_1'] +     " " + overlap_df['name_2'] + "_" + overlap_df['name_3'] + " " + overlap_df['strand']
overlap_df = overlap_df.drop(['name_1','name_2','name_3','name_4'], axis=1)


# merge with proms
proms_df = pd.DataFrame.from_dict(proms, orient = 'index', columns=['seq'])
overlap_df = overlap_df.merge(proms_df, left_on = 'keys', right_index = True)


# save
overlap_df.to_csv(dir_path + "/../data/promoter/human_epdnew_hg38_noTATA_scan_500.csv")


