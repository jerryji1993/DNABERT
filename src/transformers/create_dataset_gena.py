"""
This script generates sequence of kmers generated from a sequence
the sequence is created based on the sequence length from the input
 vcf file that contains the mutation
"""

import argparse
import pysam
import time
import logging
import gzip
import csv
import tqdm
import random

def run(var_file, seq_length, reference, prefix, is_random, duplicates):
    logger = configure_logger()
    write_seq(var_file,
              seq_length, 
              reference,
              logger,
              prefix,
              is_random,
              duplicates)

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
        self.index = int(index) - 1
        self.ref =  ref
        self.alt = alt
        self.ref_size = len(self.ref)
        #assert seq[self.index] == ref
        
    def __iter__(self):
        return self

    def get_ref(self):
        ref_index = self.index + self.ref_size
        str_ref = self.seq[self.index:ref_index]
        if str(str_ref) == str(self.ref):
            return True
        else:
            print("Warning: Ref seq does not match {0}".format(str_ref))
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
                       deletion_length=False,
                       is_random=False):
    genome = pysam.FastaFile(reference)

    if is_random:
        pad_5 = random.randint(1, seq_length)
        pad_3 = int(seq_length) - pad_5
    else:
        pad_5 = int(seq_length/2)
        pad_3 = int(seq_length/2)
    if deletion_length:
        start =  int(pos) - pad_5
        ex_stop = int(pos) + pad_3 + int(deletion_length)
        stop = int(pos) + pad_3
        del_seq = genome.fetch(chrom, start, ex_stop)
        seq = genome.fetch(chrom, start, stop)
        return(del_seq, seq, pad_5)
    else:
        start =  int(pos) - pad_5
        stop = int(pos) + pad_3
        seq = genome.fetch(chrom, start, stop)
        return (seq, pad_5)

def parse_vcf(vcf_file):
    with pysam.VariantFile(vcf_file, 'r') as vcf_reader:
        for record in vcf_reader.fetch():
            if record.alts[0] is None:
                pass
            
            if len(record.ref) > 3 and len(record.alts[0]) > 3:
                pass
            
            yield (record.chrom,
                   record.pos,
                   record.ref,
                   record.alts[0])

def parse_csv(csv_file):
    if csv_file.endswith("csv.gz"):
        open_func = lambda x: gzip.open(x, "rt")
    else:
        open_func = lambda x: open(x, "rt")
    with open_func(csv_file) as ifile:
        reader = csv.DictReader(ifile)
        for row in tqdm.tqdm(reader, mininterval=0.5):
            yield (row["Chromosome"],
                   row["Start_Position"],
                   row["Tumor_Seq_Allele1"],
                   row["Tumor_Seq_Allele2"])

def write_seq(variant_file: str, seq_length, reference, logger, prefix, is_random, duplicates):
    with gzip.open(f"{prefix}_ref_mut.tsv.gz", 'wt', compresslevel=6) as fout, \
         open(f"{prefix}_filtered.tsv", "w") as tsv:
        tsv.write("chrom\tpos\tref\talt\tinsert\n")
        fout.write("seq_ref\tseq_mut\n")

        if variant_file.endswith("vcf"):
            variant_generator = parse_vcf(variant_file)
        elif variant_file.endswith("csv.gz") or variant_file.endswith("csv"):
            variant_generator = parse_csv(variant_file)

        for chrom, pos, ref, alt in variant_generator:
            for _ in range(duplicates):
                if len(ref) > len(alt):
                    logger.info("Parsing a deletion chrom {0}, Pos {1}, REF {2}, ALT {3}".format(chrom,
                                                                                                pos,
                                                                                                ref,
                                                                                                alt))
                    ext_seq, ref_seq, index = reference_sequence(reference, seq_length, chrom, pos, len(ref) -1, is_random)
                    get_results = get_seq(ext_seq, ref, alt, index)
                    mutant_seq = get_results.generate_del()
                else:
                    ref_seq, index = reference_sequence(reference, seq_length,
                                                chrom, pos, is_random=is_random)
                    get_results = get_seq(ref_seq, ref, alt, index)
                    mutant_seq = get_results.generate_mutant()
                tsv.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{index}\n")
                logger.info("{0}, {1}, {2}, {3}, {4}, {5}".format(chrom, 
                                                                pos, 
                                                                ref, 
                                                                alt, 
                                                                ref_seq, mutant_seq))
                fout.write(ref_seq + '\t' + mutant_seq + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', dest='var_file',
                        help="variant file containing the variants to mutate", required=True)
    parser.add_argument('-p', dest='prefix',
                        help="output prefix", required=True)
    parser.add_argument('-s', dest='seq_length',
                        help="length of the string sequence to create kmers", default=512, type=int)
    parser.add_argument('--random', action="store_true",
                        help="set flag to place mutation at random places", default=False)
    parser.add_argument('-d', dest='duplicates',
                        help="total number of duplicates, default=1", default=1, type=int)
    parser.add_argument('-r', dest='reference',
                        help="reference genome")

    args = parser.parse_args()
    if not args.random:
        args.duplicates = 1
    run(args.var_file, args.seq_length, args.reference, args.prefix, args.random, args.duplicates)
