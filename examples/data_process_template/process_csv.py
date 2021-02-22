import csv
import os
import argparse
import random
from process_pretrain_data import get_kmer_sentence


max_length = 0

def Process_pair(args):
    random.seed(24)

    root_path = args.file_path.split('/')[-1]
    train_seq1_file = open(args.file_path+"/"+root_path+"_enhancer.fasta", "r")
    train_seq2_file = open(args.file_path+"/"+root_path+"_promoter.fasta", "r")
    train_label_file = open(args.file_path+"/"+root_path+"_label.txt", "r")
    test_seq1_file = open(args.file_path+"/"+root_path+"_enhancer_test.fasta", "r")
    test_seq2_file = open(args.file_path+"/"+root_path+"_promoter_test.fasta", "r")
    test_label_file = open(args.file_path+"/"+root_path+"_label_test.txt", "r")

    train_seq1 = train_seq1_file.readlines()
    train_seq2 = train_seq2_file.readlines()
    train_label = train_label_file.readlines()
    test_seq1 = test_seq1_file.readlines()
    test_seq2 = test_seq2_file.readlines()
    test_label = test_label_file.readlines()

    train_lines = []
    test_lines = []
    for i in range(len(train_label)):
        train_lines.append([train_seq1[2*i+1], train_seq2[2*i+1], train_label[i]])
    for i in range(len(test_label)):
        test_lines.append([test_seq1[2*i+1], test_seq2[2*i+1], test_label[i]])

    random.shuffle(train_lines)

    if args.dev:
        num_dev = int(len(train_lines)/10)
        dev_lines = train_lines[:num_dev]
        train_lines = train_lines[num_dev:]
    
    output_path = make_path(args)

    f_train = open(os.path.join(output_path, "train.tsv"), 'wt')
    train_w = csv.writer(f_train, delimiter='\t')
    train_w.writerow(["seq1", "seq2", "label"])
    if args.dev:
        f_dev = open(os.path.join(output_path, "dev.tsv"), 'wt')
        dev_w = csv.writer(f_dev, delimiter='\t')
        dev_w.writerow(["seq1", "seq2", "label"])
        os.makedirs(os.path.join(output_path, "test"))
        f_test = open(os.path.join(output_path, "test", "dev.tsv"), 'wt')
        test_w = csv.writer(f_test, delimiter='\t')
        test_w.writerow(["seq1", "seq2", "label"])
    else:
        f_test = open(os.path.join(output_path, "dev.tsv"), 'wt')
        test_w = csv.writer(f_test, delimiter='\t')
        test_w.writerow(["seq1", "seq2", "label"])

    def write_file_pair(lines, writer, seq1_index=0, seq2_index=1, label_index=2):
        for line in lines:
            seq1 = get_kmer_sentence(line[seq1_index],args.kmer)
            seq2 = get_kmer_sentence(line[seq2_index],args.kmer)
            writer.writerow([seq1, seq2, str(int(line[label_index]))])

    write_file_pair(train_lines, train_w)
    write_file_pair(test_lines, test_w)
    
    if args.dev:
        write_file_pair(dev_lines, dev_w)
    

def make_path(args):
    output_path = args.output_path if args.output_path else os.path.join(args.file_path, str(args.kmer))
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    return output_path
    

def Process(args):
    random.seed(24)

    train = os.path.join(args.file_path, "train.csv")
    test = os.path.join(args.file_path, "test.csv")
    train_file =  open(train, "r", encoding="utf-8-sig")
    test_file =  open(test, "r", encoding="utf-8-sig")

    train_lines = list(csv.reader(train_file, delimiter=",", quotechar=None))[1:]
    test_lines = list(csv.reader(test_file, delimiter=",", quotechar=None))[1:]

    random.shuffle(train_lines)
    random.shuffle(test_lines)

    if args.dev:
        num_dev = int(len(train_lines)/9)
        dev_lines = train_lines[:num_dev]
        train_lines = train_lines[num_dev:]

    print(train_lines[0])

    output_path = make_path(args)

    suffix = '.csv' if args.csv else '.tsv'
    delimiter = ',' if args.csv else '\t'


    f_train = open(os.path.join(output_path, "train"+suffix), 'wt')
    train_w = csv.writer(f_train, delimiter=delimiter)
    train_w.writerow(["sentence", "label"])
    if args.dev:
        f_dev = open(os.path.join(output_path, "dev"+suffix), 'wt')
        dev_w = csv.writer(f_dev, delimiter=delimiter)
        dev_w.writerow(["sentence", "label"])
        f_test = open(os.path.join(output_path, "test"+suffix), 'wt')
        test_w = csv.writer(f_test, delimiter=delimiter)
        test_w.writerow(["sentence", "label"])
    else:
        f_test = open(os.path.join(output_path, "dev"+suffix), 'wt')
        test_w = csv.writer(f_test, delimiter=delimiter)
        test_w.writerow(["sentence", "label"])
    


    def write_file(lines, writer, seq_index=2, label_index=3):
        global max_length
        for line in lines:
            sentence = get_kmer_sentence(line[seq_index],args.kmer)
            if len(sentence.split()) > max_length:
                max_length = len(sentence.split())
            writer.writerow([sentence, str(line[label_index])])

    write_file(train_lines, train_w, args.seq_index, args.label_index)
    write_file(test_lines, test_w, args.seq_index, args.label_index)
    
    if args.dev:
        write_file(dev_lines, dev_w)
    

    print("max length: %d" % (max_length))

    




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
        "--dev",
        action="store_true",
        help="Use this flag to split data as (8:1:1), else (9:1)",
    )
    parser.add_argument(
        "--csv",
        action="store_true",
        help="if output csv file or not, if not, output tsv",
    )
    parser.add_argument(
        "--pair",
        action="store_true",
        help="Use this flag to split data as (8:1:1), else (9:1)",
    )
    parser.add_argument(
        "--seq_index",
        default=2,
        type=int,
        help="index of seq in the original csv file",
    )
    parser.add_argument(
        "--label_index",
        default=3,
        type=int,
        help="index of label in the original csv file",
    )
    args = parser.parse_args()
    
    if args.pair:
        Process_pair(args)
    else:
        Process(args)


if __name__ == "__main__":
    main()
