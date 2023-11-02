#!/research/bsi/tools/biotools/python/3.5.2/bin/python

"""
This scripts annotates the a vcf
from a given input bed file
"""

import os
import pysam.bcftools as bcftools
import collections
import argparse
import sys
import csv
import pandas as pd
import pprint
import subprocess
import base64

BCFTOOLS="/usr/local/biotools/bcftools/bcftools-1.9/bin/bcftools"

def main():
    args = parse_args()
    run(args.vcf_file, args.bed_file, args.annotation_name)


def run(vcf_file, bed_file, annotation_name):
    run_bcftools = annotate_vcf(vcf_file, bed_file, annotation_name)


def parse_args():
    """
    parsing arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v',dest='vcf_file',
                        help='vcf file to annotate',
                        required=True)
    parser.add_argument('-b', dest='bed_file',
                        help='bed file requires it to bgziped and tabixd',
                        required=True)
    parser.add_argument('-a', dest='annotation_name',
                        help='name of the type of annotation',
                        required=True)
    args = parser.parse_args()
    return args


def validate_file(input_file):
    assert os.path.isfile(input_file)
    if input_file.endswith(".gz"):
        pass
    else:
        print("input file needs to be gzip compressed")


def annotate_vcf(vcf_file, bed_file, annotation_name):
    validate_file(vcf_file)
    validate_file(bed_file)
    vcf_header = "<(echo '##INFO=<ID=" + "'" + annotation_name + "'" + ",Number=1,Type=String,Description=" +  '"' + annotation_name + '"' + ">" + "'" + ")"
    cmd = BCFTOOLS + " annotate -a " + bed_file + " -c CHROM,FROM,TO," + annotation_name + " -h " + vcf_header + " " + vcf_file
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE,executable="/bin/bash")
    out, err = p.communicate()
    with open('test.out.vcf', 'w', encoding='utf-8') as fout:
        fout.write(out.decode('utf-8'))

if __name__ == "__main__":
    main()
