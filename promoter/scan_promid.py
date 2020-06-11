######## ::: Scan input fasta and make predictions ::: ########

#! conda install biopython

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

dir_path = os.path.dirname(os.path.realpath('scan.py'))


def clean_seq(s):
    ns = s.upper()    
    pattern = re.compile(r'\s+')
    ns = re.sub(pattern, '', ns)
    ns = re.sub(r'[^a-zA-Z]{1}', 'N', ns)
    return ns


## read in fasta
fasta_dir = dir_path + "/../data/genome/GRCh38.chr1to22.fa"
fasta = SeqIO.to_dict(SeqIO.parse(open(fasta_dir),'fasta'))

# convert seq to str and clean seq
fasta = {k: clean_seq(str(v.seq)) for k,v in fasta.items()}

def encode(seq, strand):
    enc_mat = np.append(np.eye(4), [[0,0,0,0]], axis=0)
    enc_mat = enc_mat.astype(np.bool)
    mapping_pos = dict(zip("ACGTN", range(5)))
    mapping_neg = dict(zip("TGCAN", range(5)))
    
    if(strand == "+"):
        seq2 = [mapping_pos[i] for i in seq]
    else:
        seq = seq[::-1]
        seq2 = [mapping_neg[i] for i in seq]
    return enc_mat[seq2]

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def scan_promid(idx, seq, half_size, key, strand):
    subseq = seq[idx - half_size: idx + half_size + 1]
    
    # make index
    if strand == "-":
        idx1 = len(seq) - (idx + half_size) + 1
        idx2 = len(seq) - (idx - half_size) + 1
    else:
        idx1 = idx - half_size + 1
        idx2 = idx + half_size + 1
        
    ck = "{}:{}-{} {}".format(key, idx1, idx2, strand)
    
    if not 'N' in subseq:
        return (ck, subseq)
    else:
        return None

    
# multiprocessing
print("Number of processors: ", mp.cpu_count(), flush=True)
pool = mp.Pool(mp.cpu_count())

# globals
scan_step = 100
half_size = 500
batch_size = 128
start = half_size
cutoff = 0.5
proms = {}
scores = {}

for key in fasta.keys():
    print("* Processing " + key, flush=True)
    for strand in ["+", "-"]:
        print("  {} strand".format("Positive" if strand == "+" else "Negative"), flush=True)
        seq = fasta[key]
        if strand == "-": 
            seq = seq[::-1]
            
        print("\t - Scanning...", end = " ", flush=True)
        putative = pool.starmap(scan_promid, [(i, seq, half_size, key, strand) for i in range(start, len(seq) - half_size, scan_step)]) # multiprocessing
        putative = list(filter(None,putative)) # filter subseq with N
        proms.update({k: v for (k,v) in putative}) # store in dict
        coords = [k for (k,v) in putative] # store coordinates separately
        putative = [v for (k,v) in putative] # use only seq to predict
        probs = []
        print("Done! Number of subsequences: {}".format(len(putative)), flush=True)
        ######### predict here, use putative #########
        new_graph = tf.Graph()
        print("\t - Starting Tensorflow...", end = " ", flush=True)
        with tf.Session(graph=new_graph) as sess:
            print("loading saved model...", end = " ", flush=True)
            tf.saved_model.loader.load(sess, [tf.saved_model.tag_constants.SERVING], "/projects/b1017/Jerry/PromID/promid/models/model_scan")
            saver = tf.train.Saver()
            saver.restore(sess, "/projects/b1017/Jerry/PromID/promid/models/model_scan/variables/variables")
            input_x = tf.get_default_graph().get_tensor_by_name("input_prom:0")
            y = tf.get_default_graph().get_tensor_by_name("output_prom:0")
            kr = tf.get_default_graph().get_tensor_by_name("kr:0")
            in_training_mode = tf.get_default_graph().get_tensor_by_name("in_training_mode:0")  
            print("loaded!", flush=True)
            
            # predict
            print("\t - One-hot encoding scanned subsequences...", end = " ", flush=True)
            encoded = [encode(fa, strand) for fa in putative]
            print("Done!", flush=True)
            print("\t - Begin prediction", flush=True)
            i = 1
            for batch in chunks(encoded, batch_size):
                if i % 100 == 0:
                    print("\t\tbatch: ", i, "out of ", math.ceil(len(encoded)/batch_size), flush=True)
                pred = sess.run(y, feed_dict={input_x: batch, kr: 1.0, in_training_mode: False})
                probs.extend([prob[0] for prob in pred])
                i += 1
            print("\t Prediction finished!", flush=True)
                
        ######### output is a list: probs #########
        scores.update({coords[i]: probs[i] for i in range(len(probs)) if probs[i] > cutoff}) # filter

        
# filter
proms = {k: v for k,v in proms.items() if k in scores.keys()}
        
        
pool.close()    


with open(dir_path + "/../results/benchmark/promid/GRCh38.chr1to22.proms.pkl", 'wb') as f:
    pickle.dump(proms, f)

with open(dir_path + "/../results/benchmark/promid/GRCh38.chr1to22.scores.pkl", 'wb') as f:
    pickle.dump(scores, f)

    
# create bed file of the positive predictions
bed = []
for i, x in enumerate(list(proms)):
    delim = re.split('[^a-zA-Z0-9+-]',x)
    coords = re.split('-',delim[1])
    bed.append("{}\t{}\t{}\tprom{}\t0\t{}".format(delim[0],coords[0],coords[1],i,delim[2]))
with open(dir_path + "/../results/benchmark/promid/GRCh38.chr1to22.proms.bed", 'w+') as f:
    f.write('\n'.join(bed))


##### pybedtools #####
# read in true TSS bed file, and make true promoter file -500 to 500 
TSS_bed = pybedtools.BedTool(dir_path + "/../data/promoter/human_epdnew_hg38.bed")
true_prom_bed = TSS_bed.slop(g=dir_path + "/../data/genome/chrom.sizes",l=499,r=500) # -500 to 500

# read in predicted promoter bed files
prom_bed = pybedtools.BedTool(dir_path + "/../results/benchmark/promid/GRCh38.chr22.proms.bed")

# intersect
# TP
overlap_TP = prom_bed.intersect(true_prom_bed,
                                s=True, # Require same strandedness.  That is, only report hits in that overlap A on the **same** strand.
                                f=0.5, # Minimum overlap required as a fraction of A
                                u=True, # Write the original A entry **once** if **any** overlaps found in B
                                r=True # Require that the fraction overlap be reciprocal for A AND B
                                )

# FP
overlap_FP = prom_bed.intersect(true_prom_bed, 
                                v=True, # Only report those entries in A that have **no overlaps** with B.
                                s=True, # Require same strandedness.  That is, only report hits in that overlap A on the **same** strand.
                                f=0.5, # Minimum overlap required as a fraction of A
                                r=True # Require that the fraction overlap be reciprocal for A AND B
                                )

# FN
overlap_FN = true_prom_bed.intersect(prom_bed,
                                     v=True, # Only report those entries in A that have **no overlaps** with B.
                                     s=True, # Require same strandedness.  That is, only report hits in that overlap A on the **same** strand.
                                     f=0.5, # Minimum overlap required as a fraction of A
                                     r=True # Require that the fraction overlap be reciprocal for A AND B
                                    )

# metrics
TP = len(overlap_TP)
FP = len(overlap_FP)
FN = len(overlap_FN)

precision = TP/(TP+FP)
recall = TP/(TP+FN)
F1 = 2*(precision*recall)/(precision+recall)


print("precision = ", precision, flush=True)
print("recall = ", recall, flush=True)
print("F1 = ", F1, flush=True)

