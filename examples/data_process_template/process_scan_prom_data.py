import argparse
import os
import csv
import numpy as np
from process_pretrain_data import get_kmer_sentence




def Process(args):

    SCAN_LIST = [int(500/(args.slide-1))*i for i in range(args.slide)]

    old_file = open(args.file_path, "r", encoding="utf-8-sig")
    old_lines = list(csv.reader(old_file, delimiter=",", quotechar=None))[1:]

    if args.output_path:
        root_path = args.output_path + "/"
    else:
        root_path = "/".join(args.file_path.split("/")[:-1]) + "/" + str(args.kmer) + "/"
    if not os.path.exists(root_path):
        os.makedirs(root_path)

    labels = np.array([])
    new_file = open(root_path+"dev.tsv", 'wt')
    tsv_w = csv.writer(new_file, delimiter='\t')
    tsv_w.writerow(["setence", "label"])

    for line in old_lines:
        label = line[6]
        labels = np.append(labels, int(label))

        for index in SCAN_LIST:
            sub_sequence = line[8][index:index+500]
            sub_sentence = get_kmer_sentence(sub_sequence, kmer=args.kmer)
            tsv_w.writerow([sub_sentence, label]) 
    
    np.save(root_path+"label.npy", labels)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--kmer",
        default=1,
        type=int,
        help="K-mer",
    )
    parser.add_argument(
        "--file_path",
        default=None,
        type=str,
        help="The path of the file to be processed",
    )
    parser.add_argument(
        "--output_path",
        default=None,
        type=str,
        help="The path of the processed data",
    )
    parser.add_argument(
        "--slide",
        default=11,
        type=int,
        help="How many 500s to use for the predictes result of 1000",
    )
    args = parser.parse_args()

    Process(args)

    


if __name__ == "__main__":
    main()
