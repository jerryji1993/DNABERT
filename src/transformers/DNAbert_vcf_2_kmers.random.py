#!/research/bsi/projects/PI/tertiary/Couch_Fergus_coucf/s123456.general_utility/python_virtual/bin/python

"""
This script generates sequence of kmers generated from a sequence
the sequence is created based on the sequence length from the input
 vcf file that contains the mutation
"""


import argparse
from pysam import FastaFile, VariantFile
import time
import logging
import random
import gzip
from typing import Tuple, Union

def run(vcf_file, kmer, seq_length, reference, outdir, duplicates, random_insert):
    logger = configure_logger()
    if not random_insert:
        compute_string = parse_vcf(vcf_file, 
                                   kmer, 
                                   seq_length, 
                                   reference, 
                                   duplicates, outdir, logger)
    else:
        compute_string = parse_vcf(vcf_file,
                                   kmer,
                                   seq_length,
                                   reference, 
                                   duplicates, outdir, logger)

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

    def __init__(self, seq, ref, alt, index):
        self.seq =  seq
        self.size = len(self.seq)
        self.index = int(index) -1
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
            print("Warning: Ref seq base {0} does not match the ref in vcf {1} for seq {2}".format(str_ref, self.ref, self.seq))
#            exit()

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
    genome = FastaFile(reference)
    pad = int(seq_length/2)
    if deletion_length:
        start =  int(pos) - pad
        ex_stop = int(pos) + pad + int(deletion_length)
        stop = int(pos) + pad
        del_seq = genome.fetch(chrom, start, ex_stop)
        seq = genome.fetch(chrom, start, stop)
        return(del_seq, seq)
    else:
        start =  int(pos) - pad
        stop = int(pos) + pad
        seq = genome[chrom][start:stop]
        return seq

def reference_sequence_random(reference, 
                              seq_length, 
                              chrom, 
                              pos, logger,
                              deletion_length=None) -> Union[Tuple[str, str, int], Tuple[str, int]]:
    genome = FastaFile(reference)
    insert = random.randint(1, seq_length)
    logger.info("position of the insert is {0}".format(insert))
    pad_5 = int(seq_length) 
    pad_3 = int(seq_length) - insert
    if deletion_length:
        start =  int(pos) - insert
        ex_stop = int(pos) + pad_3 + int(deletion_length)
        stop = int(pos) + pad_3
        del_seq = genome.fetch(chrom, start, ex_stop)
        seq = genome.fetch(chrom, start, stop)
        return(del_seq, seq, insert)
    else:
        start =  int(pos) - insert
        stop = int(pos) + pad_3
        seq = genome.fetch(chrom, start, stop)
        return(seq, insert)


def sliding_windown(seq, kmer):
    seq_str = ''
    for i in range(0,len(seq),1):
        kmer_str = str(seq[i:i + kmer])
        if len(kmer_str) != kmer:
            pass
        else:
            seq_str += (kmer_str) + ' '
    return seq_str


def execute_kmers(reference, seq_length,chrom, pos, ref, alt, logger, 
                  deletion=False):
    if deletion:
#        ext_seq, ref_seq = reference_sequence(reference, seq_length, 
#                                              chrom, pos, len(ref) -1)
        ext_seq, ref_seq, index = reference_sequence_random(reference, seq_length,
                                                     chrom, pos, logger, len(ref)-1)
        get_results = get_seq(ext_seq, ref, alt, index)
        mutant_seq = get_results.generate_del()
        # logger.info("{0}, {1}, {2}, {3}, {4}, {5}".format(chrom, 
        #                                                   pos, 
        #                                                   ref, 
        #                                                   alt, 
        #                                                   ref_seq, mutant_seq))
        return ref_seq, mutant_seq, index
    else:
#        ref_seq = reference_sequence(reference, seq_length,
#                                     chrom, pos)
        ref_seq, index = reference_sequence_random(reference, seq_length,
                                            chrom, pos, logger)
        get_results = get_seq(ref_seq, ref, alt, index)
        mutant_seq = get_results.generate_mutant()
        # logger.info("{0}, {1}, {2}, {3}, {4}, {5}".format(chrom, 
        #                                                   pos, 
        #                                                   ref, 
        #                                                   alt, 
        #                                                   ref_seq, mutant_seq))
        return ref_seq, mutant_seq, index

def parse_vcf(vcf_file, kmer, seq_length, 
              reference, duplicates, outdir, logger, random=False):
    out_seq = outdir + "/" + "DNAbert_input_reference_mutant.txt.gz"
    out_vcf = outdir + "/" + "DNAbert_input_filtered.vcf.gz"
    with gzip.open(out_seq, 'wt', compresslevel=6) as fout:
        vcf_in = VariantFile(vcf_file, 'r')
        vcf_in.header.add_meta('INFO', items=[('ID',"INSERT"), ('Number',1), ('Type','Integer'), ('Description','Variant insertion position')])
        vcf_out = VariantFile(out_vcf, 'w', header=vcf_in.header)
        fout.write("seq_ref\tseq_mut\n")
        for rec in vcf_in:
            if rec.alts is None or rec.alts[0] is None:
                pass
            else:
                if len(rec.ref) > 2 and len(rec.alts[0]) > 2:
                    pass
                else:
                    is_deletion = len(rec.alts[0]) < len(rec.ref)
                    for q in range(duplicates):
                        ref_seq, mutant_seq, index = execute_kmers(reference,
                                                            seq_length,
                                                            rec.chrom,
                                                            rec.pos,
                                                            rec.ref,
                                                            rec.alts[0],
                                                            logger,
                                                            is_deletion)
                        fout.write(sliding_windown(ref_seq, kmer) + '\t' + sliding_windown(mutant_seq, kmer) +'\n')
                        rec.info["INSERT"] = index
                        vcf_out.write(rec)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', dest='vcf_file',
                        help="vcf file containing the variants to mutate", required=True)
    parser.add_argument('-k', dest='kmer',
                        help="length of kmer to generate a sequence, default=6", default=6, type=int)
    parser.add_argument('-s', dest='seq_length',
                        help="length of the string sequence to create kmers, default=510", default=510, type=int)
    parser.add_argument('-d', dest='duplicates',
                        help="total number of duplicates, default=1", default=1, type=int)
    parser.add_argument('-r', dest='reference',default="/research/bsi/data/refdata/app/gatk_bundle/human/2.8/b37/processed/2015_11_04/allchr.fa",
                        help="reference genome")
    parser.add_argument('-o', dest='outdir',required=True,
                        help="path to output directory")
    parser.add_argument('-x', dest='random_insert', action='store_false',
                        help="flag to invoke random insert of the mutation, otherwise it would place it in the center of the sequence")
    args = parser.parse_args()
    run(args.vcf_file, args.kmer, 
        args.seq_length, args.reference, args.outdir, args.duplicates, args.random_insert)


 
