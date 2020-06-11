######## ::: Scan input fasta and make predictions ::: ########
#! conda install biopython

import re
import os
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing as mp

dir_path = os.path.dirname(os.path.realpath('scan.py'))

def clean_seq(s):
    ns = s.upper()    
    pattern = re.compile(r'\s+')
    ns = re.sub(pattern, '', ns)
    ns = re.sub(r'[^a-zA-Z]{1}', 'N', ns)
    return ns

## read in fasta
fasta_dir = dir_path + "/../data/genome/GRCh38.chr22.fa"
fasta = SeqIO.to_dict(SeqIO.parse(open(fasta_dir),'fasta'))

# convert seq to str and clean seq
fasta = {k: clean_seq(str(v.seq)) for k,v in fasta.items()}

# def scan(idx):
def scan(idx, seq, half_size, key, strand):
    subseq = seq[idx - half_size: idx + half_size]
    
    # make index
    if strand == "-":
        idx1 = len(seq) - (idx + half_size) + 1
        idx2 = len(seq) - (idx - half_size)
    else:
        idx1 = idx - half_size + 1
        idx2 = idx + half_size
        
    ck = "{}:{}-{} {}".format(key, idx1, idx2, strand)
    
    if not 'N' in subseq:
        return (ck, subseq)
    else:
        return None

# multiprocessing
print("Number of processors: ", mp.cpu_count())
pool = mp.Pool(mp.cpu_count())

# globals
scan_step = 50
half_size = 250
start = half_size
cutoff = 0.5
proms = {}
scores = {}

for key in fasta.keys():
    print("Scanning " + key)
    for strand in ["+", "-"]:
        print("{} strand".format(strand))
        seq = fasta[key]
        if strand == "-": 
            seq = seq[::-1]
#         putative = pool.map(scan, range(start, len(seq) - half_size, scan_step)) # multiprocessing
        putative = pool.starmap(scan_promid, [(i, seq, half_size, key, strand) for i in range(start, len(seq) -half_size, scan_step)]) # multiprocessing
        putative = list(filter(None,putative)) # filter subseq with N
        proms = {k: v for (k,v) in putative} # store in dict
        coords = [k for (k,v) in putative] # store coordinates separately
        putative = [v for (k,v) in putative] # use only seq to predict
        
        ######### predict here, use putative #########
        
        ######### output is a list: probs #########
        scores = {coords[i]: probs[i] for i in range(len(probs)) if probs[i] > cutoff} # filter
        proms = {k: v for k,v in proms.items() if k in scores.keys()} # filter
        
pool.close()    

