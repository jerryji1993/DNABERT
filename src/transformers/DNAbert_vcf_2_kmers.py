#!/research/bsi/projects/PI/tertiary/Couch_Fergus_coucf/s123456.general_utility/python_virtual/bin/python

"""
This script generates sequence of kmers generated from a sequence
the sequence is created based on the sequence length from the input
 vcf file that contains the mutation
"""


import csv
import argparse
import pyfaidx
import vcf
from pyfaidx import Fasta
import time
import datetime
import logging

def run(vcf_file, kmer, seq_length, reference):
    logger = configure_logger()
    compute_string = parse_vcf(vcf_file, 
                               kmer, 
                               seq_length, 
                               reference, logger)

def configure_logger():
    """
    setting up logging
    """
    logger = logging.getLogger('DNAbert_vcf_to_kmer')
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(time.strftime("DNAbert_vcf_to_kmer-%Y%m%d.log"))
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s'\t'%(name)s'\t'%(levelname)s'\t'%(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


class get_seq(object):

    def __init__(self, seq, ref, alt):
        self.seq =  seq
        self.size = len(self.seq)
        self.index = int(self.size/2) -1
        self.ref =  ref
        self.alt = alt
        self.ref_size = len(self.ref)
        
    def __iter__(self):
        return self

    def get_ref(self):
        ref_index = self.index + self.ref_size
        str_ref = self.seq[self.index:ref_index]
        if str(str_ref) == str(self.ref):
            return True
        else:
            print("ERROR: Ref seq does not match")
            exit()

    def generate_mutant(self):
        check_ref = self.get_ref()
        mut_seq = str(self.seq[:self.index]) + str(self.alt) + str(self.seq[self.index + 1:])
        return mut_seq

    def generate_del(self):
        mut_seq =  str(self.seq[:self.index + 1]) + str(self.seq[self.index + self.ref_size :])
        return mut_seq


def reference_sequence(reference, 
                       seq_length, 
                       chrom, 
                       pos, 
                       deletion_length=False):
    genome = Fasta(reference)
    pad = int(seq_length/2)
    if deletion_length:
        start =  int(pos) - pad
        ex_stop = int(pos) + pad + int(deletion_length)
        stop = int(pos) + pad
        del_seq = genome[chrom][start:ex_stop]
        seq = genome[chrom][start:stop]
        return(del_seq, seq)
    else:
        start =  int(pos) - pad
        stop = int(pos) + pad
        seq = genome[chrom][start:stop]
        return seq


def sliding_windown(seq, kmer):
    seq_str = ''
    for i in range(0,len(seq),1):
        kmer_str = str(seq[i:i + kmer])
        if len(kmer_str) != kmer:
            pass
        else:
            seq_str += (kmer_str) + ' '
    return seq_str


def parse_vcf(vcf_file, kmer, seq_length, reference, logger):
    with open("DNAbert_input.txt", 'w') as fout:
        vcf_reader =  vcf.Reader(open(vcf_file, 'r'))
        for record in vcf_reader:
            if record.is_deletion:
                logger.info("Parsing a deletion chrom {0}, Pos {1}, REF {2}, ALT {3}".format(record.CHROM, 
                                                                                             record.POS,
                                                                                             record.REF,
                                                                                             record.ALT))
                ext_seq, ref_seq = reference_sequence(reference, seq_length, 
                                                      record.CHROM, record.POS, len(record.REF) -1)
                get_results = get_seq(ext_seq, record.REF, record.ALT[0])
                mutant_seq = get_results.generate_del()
                logger.info("{0}, {1}, {2}, {3}, {4}, {5}".format(record.CHROM, 
                                                                  record.POS, 
                                                                  record.REF, 
                                                                  record.ALT[0], 
                                                                  ref_seq, mutant_seq))
                fout.write(sliding_windown(mutant_seq, kmer) + '\n')
            else:
                ref_seq = reference_sequence(reference, seq_length,
                                             record.CHROM, record.POS)
                get_results = get_seq(ref_seq, record.REF, record.ALT[0])
                mutant_seq = get_results.generate_mutant()          
                logger.info("{0}, {1}, {2}, {3}, {4}, {5}".format(record.CHROM, 
                                                                  record.POS, 
                                                                  record.REF, 
                                                                  record.ALT[0], 
                                                                  ref_seq, mutant_seq))
                fout.write(sliding_windown(mutant_seq, kmer) + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', dest='vcf_file',
                        help="vcf file containing the variants to mutate", required=True)
    parser.add_argument('-k', dest='kmer',
                        help="length of kmer to generate a sequence", required=True, type=int)
    parser.add_argument('-s', dest='seq_length',
                        help="length of the string sequence to create kmers", default=512, type=int)
    parser.add_argument('-r', dest='reference',default="/research/bsi/data/refdata/app/gatk_bundle/human/2.8/b37/processed/2015_11_04/allchr.fa",
                        help="reference genome")
    args = parser.parse_args()
    run(args.vcf_file, args.kmer, args.seq_length, args.reference)


 
