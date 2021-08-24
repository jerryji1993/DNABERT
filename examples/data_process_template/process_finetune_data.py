import argparse
import csv
import os
import random
import numpy as np
from process_pretrain_data import get_kmer_sentence

max_length = 0 

def write_file(lines, path, kmer, head=True, seq_index=0, label_index=1):
        with open(path, 'wt') as f:
            tsv_w = csv.writer(f, delimiter='\t')
            if head:
                tsv_w.writerow(["setence", "label"]) 
            for line in lines:
                if kmer == 0:
                    sentence = str(line[seq_index])
                else:
                    sentence = str(get_kmer_sentence("".join(line[seq_index].split()), kmer))
                if label_index == None:
                    label = "0"
                else:
                    label = str(line[label_index])
                tsv_w.writerow([sentence, label]) 


def Shuffle(args):
    old_file =  open(args.file_path, "r", encoding="utf-8-sig")
    old_lines = list(csv.reader(old_file, delimiter="\t", quotechar=None))[1:]
    random.shuffle(old_lines)

    write_file(old_lines, args.file_path, 0)

def Find_train(args):
    random.seed(args.seed)

    tata = args.file_path + "/TATA_249to50.tsv"
    notata = args.file_path + "/noTATA_249to50.tsv"
    tata_file =  open(tata, "r", encoding="utf-8-sig")
    notata_file =  open(notata, "r", encoding="utf-8-sig")
    tata_lines = list(csv.reader(tata_file, delimiter="\t", quotechar=None))[1:]
    notata_lines = list(csv.reader(notata_file, delimiter="\t", quotechar=None))[1:]

    tata_test = args.file_path + "/tata_test.tsv"
    notata_test = args.file_path + "/notata_test.tsv"
    tata_test_file =  open(tata_test, "r", encoding="utf-8-sig")
    notata_test_file =  open(notata_test, "r", encoding="utf-8-sig")
    tata_test_lines = list(csv.reader(tata_test_file, delimiter="\t", quotechar=None))[1:]
    notata_test_lines = list(csv.reader(notata_test_file, delimiter="\t", quotechar=None))[1:]

    
    train_lines = []

    for line in tata_lines:
        if [line[0], line[1]] not in tata_test_lines:
            train_lines.append([line[0], line[1]])
    

    for line in notata_lines:
        if [line[0], line[1]] not in notata_test_lines:
            train_lines.append([line[0], line[1]])
    
    random.shuffle(train_lines)
    random.shuffle(train_lines)

    # num_dev = int(len(train_lines)/9.0)
    # dev_lines = train_lines[:num_dev]
    # train_lines = train_lines[num_dev:]

    
    write_file(train_lines, args.file_path+"/train.tsv", args.kmer, head=False)
    # write_file(dev_lines, args.file_path+"/dev.tsv", args.kmer)

    for kmer in range(3,7):
        root_path = os.path.join(args.file_path, str(kmer))
        if not os.path.exists(root_path):
            os.makedirs(root_path)

        train_file =  open(os.path.join(args.file_path,"train.tsv"), "r", encoding="utf-8-sig")
        lines = list(csv.reader(train_file, delimiter="\t", quotechar=None))
        train_path = os.path.join(root_path,"train.tsv")

        write_file(lines, train_path, kmer)

        tata_path = os.path.join(root_path, "tata")
        notata_path = os.path.join(root_path, "notata")
        os.makedirs(tata_path)
        os.makedirs(notata_path)

        dev_lines = tata_test_lines+notata_test_lines
        dev_path = os.path.join(root_path,"dev.tsv")

        write_file(tata_test_lines, os.path.join(tata_path, "dev.tsv"), kmer)
        write_file(notata_test_lines, os.path.join(notata_path, "dev.tsv"), kmer)
        write_file(dev_lines, dev_path, kmer)

def Process_1000(args):
    random.seed(args.seed)

    tata_train = args.file_path + "TATA_scan_train.csv"
    notata_train = args.file_path + "noTATA_scan_train.csv"
    tata_train_file =  open(tata_train, "r", encoding="utf-8-sig")
    notata_train_file =  open(notata_train, "r", encoding="utf-8-sig")
    tata_train_lines = list(csv.reader(tata_train_file, delimiter=",", quotechar=None))[1:]
    notata_train_lines = list(csv.reader(notata_train_file, delimiter=",", quotechar=None))[1:]

    tata_test = args.file_path + "/TATA_scan_test.csv"
    notata_test = args.file_path + "/noTATA_scan_test.csv"
    tata_test_file =  open(tata_test, "r", encoding="utf-8-sig")
    notata_test_file =  open(notata_test, "r", encoding="utf-8-sig")
    tata_test_lines = list(csv.reader(tata_test_file, delimiter=",", quotechar=None))[1:]
    notata_test_lines = list(csv.reader(notata_test_file, delimiter=",", quotechar=None))[1:]
    

    print("Original:")
    print("tata train: %d" % (len(tata_train_lines)))
    print("notata train: %d" % (len(notata_train_lines)))
    print("tata test: %d" % (len(tata_test_lines)))
    print("tata test: %d" % (len(notata_test_lines)))

    random.shuffle(tata_train_lines)
    random.shuffle(notata_train_lines)
    random.shuffle(tata_test_lines)
    random.shuffle(notata_test_lines)

    
    notata_train_lines = notata_train_lines[:len(tata_train_lines)]
    notata_test_lines = notata_test_lines[:len(tata_test_lines)]
    with open(os.path.join(args.file_path, "notata_test_id"), "w") as f:
        tsv_w = csv.writer(f, delimiter=',')
        tsv_w.writerow(["index", "chrom", "start", "end", "name", "strand", "keys", "id"]) 
        for line in notata_test_lines:
            tsv_w.writerow([line[0], line[1], line[2], line[3], line[4], line[5], line[7], line[9]]) 



    # print("After:")
    # print("tata train: %d" % (len(tata_train_lines)))
    # print("notata train: %d" % (len(notata_train_lines)))
    # print("tata test: %d" % (len(tata_test_lines)))
    # print("tata test: %d" % (len(notata_test_lines)))

    # train_lines = tata_train_lines + notata_train_lines
    # test_lines = tata_test_lines + notata_test_lines


    # output_path = args.output_path if args.output_path is not None else args.file_path

    # write_file(test_lines, output_path+"/dev.tsv", args.kmer, head=False, seq_index=8, label_index=6)
    # write_file(train_lines, output_path+"/train.tsv", args.kmer, head=False, seq_index=8, label_index=6)
    # write_file(tata_test_lines, output_path+"/tata_dev.tsv", args.kmer, head=False, seq_index=8, label_index=6)
    # write_file(tata_train_lines, output_path+"/tata_train.tsv", args.kmer, head=False, seq_index=8, label_index=6)
    # write_file(notata_test_lines, output_path+"/notata_dev.tsv", args.kmer, head=False, seq_index=8, label_index=6)
    # write_file(notata_train_lines, output_path+"/notata_train.tsv", args.kmer, head=False, seq_index=8, label_index=6)

    # Process_1000_kmer(args, test_lines, train_lines, tata_test_lines, tata_train_lines, notata_test_lines, notata_train_lines)


def Process_1000_kmer(args, test_lines=None, train_lines=None, tata_test_lines=None, tata_train_lines=None, notata_test_lines=None, notata_train_lines=None):
    
    LOAD = True
    output_path = args.output_path if args.output_path is not None else args.file_path

    if test_lines == None:
        path1 = os.path.join(args.file_path,"dev.tsv")
        path2 = os.path.join(args.file_path,"train.tsv")
        path3 = os.path.join(args.file_path,"tata_dev.tsv")
        path4 = os.path.join(args.file_path,"tata_train.tsv")
        path5 = os.path.join(args.file_path,"notata_dev.tsv")
        path6 = os.path.join(args.file_path,"notata_train.tsv")

        file1 =  open(path1, "r", encoding="utf-8-sig")
        file2 =  open(path2, "r", encoding="utf-8-sig")
        file3 =  open(path3, "r", encoding="utf-8-sig")
        file4 =  open(path4, "r", encoding="utf-8-sig")
        file5 =  open(path5, "r", encoding="utf-8-sig")
        file6 =  open(path6, "r", encoding="utf-8-sig")

        test_lines = list(csv.reader(file1, delimiter="\t", quotechar=None))
        train_lines = list(csv.reader(file2, delimiter="\t", quotechar=None))
        tata_test_lines = list(csv.reader(file3, delimiter="\t", quotechar=None))
        tata_train_lines = list(csv.reader(file4, delimiter="\t", quotechar=None))
        notata_test_lines = list(csv.reader(file5, delimiter="\t", quotechar=None))
        notata_train_lines = list(csv.reader(file6, delimiter="\t", quotechar=None))

        LOAD = False



    for kmer in range(3,7):

        print(kmer)
        root_path = os.path.join(output_path, str(kmer))
        if not os.path.exists(root_path):
            os.makedirs(root_path)

        all_path = os.path.join(root_path, "all")
        # tata_path = os.path.join(root_path, "tata")
        notata_path = os.path.join(root_path, "notata")
        os.makedirs(all_path)
        # os.makedirs(tata_path)
        os.makedirs(notata_path)

        if LOAD:
            seq_index=8
            label_index=6
        else:
            seq_index=0
            label_index=1
            
        print("writing dev")
        write_file(test_lines, os.path.join(all_path,"dev.tsv"), kmer, head=False, seq_index=seq_index, label_index=label_index)
        print("writing train")
        write_file(train_lines, os.path.join(all_path,"train.tsv"), kmer, head=False, seq_index=seq_index, label_index=label_index)
        # print("writing tata dev")
        # write_file(tata_test_lines, os.path.join(tata_path,"dev.tsv"), kmer, head=False, seq_index=seq_index, label_index=label_index)
        # print("writing tata train")
        # write_file(tata_train_lines, os.path.join(tata_path,"train.tsv"), kmer, head=False, seq_index=seq_index, label_index=label_index)
        print("writing notata dev")
        write_file(notata_test_lines, os.path.join(notata_path,"dev.tsv"), kmer, head=False, seq_index=seq_index, label_index=label_index)
        print("writing notata train")
        write_file(notata_train_lines, os.path.join(notata_path,"train.tsv"), kmer, head=False, seq_index=seq_index, label_index=label_index)
    

def Process_splice(args):
    # X_train = np.load(os.path.join(args.file_path, "x_train.npy"))
    # X_dev = np.load(os.path.join(args.file_path, "x_dev.npy"))
    # Y_train = np.load(os.path.join(args.file_path, "y_train.npy"))
    # Y_dev = np.load(os.path.join(args.file_path, "y_dev.npy"))

    # assert len(X_train) == len(Y_train)
    # assert len(X_dev) == len(Y_dev)

    # for kmer in range(3,7):
    #     root_path = os.path.join(args.file_path, str(kmer))
    #     os.makedirs(root_path)
    #     f_train = open(os.path.join(root_path, "train.tsv"), "wt")
    #     f_dev = open(os.path.join(root_path, "dev.tsv"), "wt")
    #     tsv_train = csv.writer(f_train, delimiter='\t')
    #     tsv_dev = csv.writer(f_dev, delimiter='\t')
    #     tsv_train.writerow(["seq", "label"])
    #     tsv_dev.writerow(["seq", "label"])

    #     for i, seq in enumerate(X_train):
    #         sequence = get_kmer_sentence(str(seq), kmer)
    #         tsv_train.writerow([sequence, int(Y_train[i])])
        
    #     for j, seq in enumerate(X_dev):
    #         sequence = get_kmer_sentence(str(seq), kmer)
    #         tsv_dev.writerow([sequence, int(Y_dev[j])])
    
    X_test = np.load(os.path.join(args.file_path, "x_test.npy"))
    Y_test = np.load(os.path.join(args.file_path, "y_test.npy"))

    assert len(X_test) == len(Y_test)
   
    for kmer in range(3,7):
        root_path = os.path.join(args.file_path, str(kmer))
        os.makedirs(root_path)
        f_test = open(os.path.join(root_path, "dev.tsv"), "wt")
        tsv_test = csv.writer(f_test, delimiter='\t')
        tsv_test.writerow(["seq", "label"])

        for i, seq in enumerate(X_test):
            sequence = get_kmer_sentence(str(seq), kmer)
            label = int(np.where(Y_test[i]==1)[0])
            tsv_test.writerow([sequence, label])
        
        
def Process_prom_core(args):
    random.seed(args.seed)

    tata = args.file_path + "/TATA.csv"
    notata = args.file_path + "/noTATA.csv"
    tata_file =  open(tata, "r", encoding="utf-8-sig")
    notata_file =  open(notata, "r", encoding="utf-8-sig")
    tata_lines = list(csv.reader(tata_file, delimiter=",", quotechar=None))[1:]
    notata_lines = list(csv.reader(notata_file, delimiter=",", quotechar=None))[1:]

    random.shuffle(tata_lines)
    random.shuffle(notata_lines)

    num_tata_test = int(0.1*len(tata_lines))
    tata_test_lines = tata_lines[:num_tata_test]
    num_notata_test = int(0.1*len(notata_lines))
    notata_test_lines = notata_lines[:num_notata_test]

    train_lines = tata_lines[num_tata_test:] + notata_lines[num_notata_test:]
    if args.dev:
        num_dev = int(len(rest_lines)/9.0)
        dev_lines = train_lines[:num_dev]
        train_lines = train_lines[num_dev:]
    else:
        dev_lines = tata_test_lines + notata_test_lines
    
    print("Number train examples: %d" % (len(train_lines)))
    print("Number dev examples: %d" % (len(dev_lines)))

    for kmer in range(3,7):
        root_path = os.path.join(args.file_path,str(kmer))
        tata_path = os.path.join(root_path, "tata")
        notata_path = os.path.join(root_path, "notata")
        os.makedirs(tata_path)
        os.makedirs(notata_path)

        write_file(tata_test_lines, os.path.join(tata_path,"dev.tsv"), kmer, head=False, seq_index=1, label_index=2)
        write_file(notata_test_lines, os.path.join(notata_path,"dev.tsv"), kmer, head=False, seq_index=1, label_index=2)
        write_file(train_lines, os.path.join(root_path,"train.tsv"), kmer, head=False, seq_index=1, label_index=2)
        write_file(dev_lines, os.path.join(root_path,"dev.tsv"), kmer, head=False, seq_index=1, label_index=2)


def Process_pair(args):
    random.seed(args.seed)

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
    
    output_path = args.output_path if args.output_path else os.path.join(args.file_path, str(args.kmer))
    if not os.path.exists(output_path):
        os.makedirs(output_path)

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


def Process_p53_mut(args):
    random.seed(args.seed)
    
    dev = os.path.join(args.file_path, "dev.csv")
    dev_file =  open(dev, "r", encoding="utf-8-sig")

    lines = list(csv.reader(dev_file, delimiter=",", quotechar=None))[1:]

    print(lines[0])

    for kmer in range(3, 7):
        output_path = args.output_path if args.output_path else os.path.join(args.file_path, str(kmer))
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        
        write_file(lines, os.path.join(output_path, "dev.tsv"), kmer, head=True, seq_index=2, label_index=None)


def Process_p53(args):
    random.seed(args.seed)
    
    train = os.path.join(args.file_path, "train.csv")
    test = os.path.join(args.file_path, "test.csv")
    train_file =  open(train, "r", encoding="utf-8-sig")
    test_file =  open(test, "r", encoding="utf-8-sig")

    train_lines = list(csv.reader(train_file, delimiter=",", quotechar=None))[1:]
    test_lines = list(csv.reader(test_file, delimiter=",", quotechar=None))[1:]
    lines = train_lines + test_lines

    max_length = 0
    for line in lines:
        if len(line[2]) > max_length:
            max_length = len(line[2])

    random.shuffle(train_lines)
    random.shuffle(test_lines)

    if args.dev:
        num_dev = int(len(train_lines)/9)
        dev_lines = train_lines[:num_dev]
        train_lines = train_lines[num_dev:]

    print(train_lines[0])

    for kmer in range(3, 7):
        output_path = args.output_path if args.output_path else os.path.join(args.file_path, str(kmer))
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        
        write_file(train_lines, os.path.join(output_path, "train.tsv"), kmer, head=True, seq_index=2, label_index=3)
        if args.dev:
            write_file(dev_lines, os.path.join(output_path, "dev.tsv"), kmer, head=True, seq_index=2, label_index=3)
            os.makedirs(os.path.join(output_path, "test"))
            write_file(test_lines, os.path.join(output_path, "test", "dev.tsv"), kmer, head=True, seq_index=2, label_index=3)
        else:
            write_file(test_lines, os.path.join(output_path, "dev.tsv"), kmer, head=True, seq_index=2, label_index=3) 

    print("max length: %d" % (max_length))


def Seperate_p53(args):
    random.seed(args.seed)
    
    train = os.path.join(args.file_path, "train.csv")
    test = os.path.join(args.file_path, "test.csv")
    train_file =  open(train, "r", encoding="utf-8-sig")
    test_file =  open(test, "r", encoding="utf-8-sig")

    train_lines = list(csv.reader(train_file, delimiter=",", quotechar=None))[1:]
    test_lines = list(csv.reader(test_file, delimiter=",", quotechar=None))[1:]
    lines = train_lines + test_lines

    POS = []
    NEG = []

    for line in lines:
        if str(line[-1]) == '0':
            NEG.append([line[-2], line[-1]])
        else:
            POS.append([line[-2], line[-1]])

   

    for kmer in range(3,7):
        os.makedirs(os.path.join(args.file_path, "POS", str(kmer)))
        os.makedirs(os.path.join(args.file_path, "NEG", str(kmer)))

        write_file(POS, os.path.join(args.file_path, "POS", str(kmer), "dev.tsv"), kmer=kmer, head=True, seq_index=0, label_index=1)
        write_file(NEG, os.path.join(args.file_path, "NEG", str(kmer), "dev.tsv"), kmer=kmer, head=True, seq_index=0, label_index=1)
    


def Generate_prom_train_dev(args):
    # read TATA and noTATA files
    tata = args.file_path + "/noTATA_249to50.tsv"
    notata = args.file_path + "/TATA_249to50.tsv"
    tata_file =  open(tata, "r", encoding="utf-8-sig")
    notata_file =  open(notata, "r", encoding="utf-8-sig")
    tata_lines = list(csv.reader(tata_file, delimiter="\t", quotechar=None))[1:]
    notata_lines = list(csv.reader(notata_file, delimiter="\t", quotechar=None))[1:]
    
    
    # shuffle all the data and split them
    random.shuffle(tata_lines)
    random.shuffle(notata_lines)
    num_tata_test = int(len(tata_lines)*0.1)
    tata_test_lines = tata_lines[:num_tata_test]
    num_notata_test = int(len(notata_lines)*0.1)
    notata_test_lines = notata_lines[:num_notata_test]
    train_lines = tata_lines[num_tata_test:] + notata_lines[num_notata_test:]
    test_lines = tata_test_lines + notata_test_lines


    write_file(train_lines, args.file_path+"/train.tsv", args.kmer)
    write_file(test_lines, args.file_path+"/dev.tsv", args.kmer)
    write_file(tata_test_lines, args.file_path+"/tata_dev.tsv", args.kmer)
    write_file(notata_test_lines, args.file_path+"/notata_dev.tsv", args.kmer)

def Process_690(args):
    path = args.file_path
    all_folders = os.listdir(path)
    
    count = 0

    for folder in all_folders:
        # load data
        train_seq_path = os.path.join(args.file_path, folder, "train", "sequences_alph.npy")
        test_seq_path = os.path.join(args.file_path, folder, "test", "sequences_alph.npy")
        train_lab_path = os.path.join(args.file_path, folder, "train", "targets.npy")
        test_lab_path = os.path.join(args.file_path, folder, "test", "targets.npy")
        train_sequences = np.load(train_seq_path)
        test_sequences = np.load(test_seq_path)
        train_labels = np.load(train_lab_path)
        test_labels = np.load(test_lab_path)

        train_sequences = train_sequences.reshape(train_sequences.shape[0],1)
        test_sequences = test_sequences.reshape(test_sequences.shape[0],1)
        train_labels = train_labels.reshape(train_labels.shape[0],1)
        test_labels = test_labels.reshape(test_labels.shape[0],1)

        # concat sequence and labels together
        trains = list(np.concatenate((train_sequences, train_labels), axis=1))
        tests = list(np.concatenate((test_sequences, test_labels), axis=1))

        random.seed(args.seed)
        random.shuffle(trains)
        random.shuffle(trains)
        random.shuffle(tests)
        random.shuffle(tests)


        # make output path
        output_path = os.path.join(args.output_path, str(args.kmer), folder)
        if not os.path.exists(output_path):
            os.makedirs(output_path)

        

        # write files 
        f_train = open(os.path.join(output_path, "train.tsv"), 'wt')
        tsv_train = csv.writer(f_train, delimiter='\t')
        tsv_train.writerow(["sequence", "label"])
        for i in range(len(trains)):
            sentence = get_kmer_sentence(trains[i][0].decode("utf-8"), args.kmer)
            tsv_train.writerow([sentence, int(trains[i][1])])

        f_dev = open(os.path.join(output_path, "dev.tsv"), 'wt')
        tsv_dev = csv.writer(f_dev, delimiter='\t')
        tsv_dev.writerow(["sequence", "label"])
        for i in range(len(tests)):
            sentence = get_kmer_sentence(tests[i][0].decode("utf-8"), args.kmer)
            tsv_dev.writerow([sentence, int(tests[i][1])])
        

        count += 1
        print("Finish %s folders" % (count))


def Process_mouse(args):
    random.seed(args.seed)

    files = os.listdir(args.file_path)

    try:
        files.remove("3")
        files.remove("4")
        files.remove("5")
        files.remove("6")
    except ValueError:
        files = files

    files.sort()
    assert len(files) % 2 == 0

    num_task = int(len(files)/2)

    max_length = 0

    for i in range(num_task):
        index = str(i) if i > 9 else "0" + str(i)

        test_name = files[2*i].replace("test", "train")
        train_name = files[2*i+1]
        assert test_name == train_name

        test_file = os.path.join(args.file_path, files[2*i])
        train_file = os.path.join(args.file_path, files[2*i+1])
        train_file =  open(train_file, "r", encoding="utf-8-sig")
        test_file =  open(test_file, "r", encoding="utf-8-sig")
        train_lines = list(csv.reader(train_file, delimiter=",", quotechar=None))[1:]
        test_lines = list(csv.reader(test_file, delimiter=",", quotechar=None))[1:]

        print("dataset %d : %d lines" % (i, len(train_lines)))

        # random.shuffle(train_lines)

        # for kmer in range(3, 7):
        #     os.makedirs(os.path.join(args.file_path, str(kmer), index))
        #     write_file(train_lines, os.path.join(args.file_path, str(kmer), index, "train.tsv"), kmer, head=True, seq_index=2, label_index=3)
        #     write_file(test_lines, os.path.join(args.file_path, str(kmer), index, "dev.tsv"), kmer, head=True, seq_index=2, label_index=3)



def Process(args):
    if args.output_path != None:
        output_path = args.output_path
    else:
        root_path = "/".join(args.file_path.split("/")[:-1]) + "/" + str(args.kmer) + "/"
        output_path = root_path + args.file_path.split("/")[-1]
        if not os.path.exists(root_path):
            os.makedirs(root_path)
 
    old_file =  open(args.file_path, "r", encoding="utf-8-sig")
    lines = list(csv.reader(old_file, delimiter=args.delimiter, quotechar=None))

    write_file(lines, output_path, args.kmer, head=args.head, seq_index=args.seq_index, label_index=args.label_index)
        

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--kmer",
        default=1,
        type=int,
        help="K-mer",
    )
    parser.add_argument(
        "--seed",
        default=24,
        type=int,
        help="Which random seed to use",
    )
    parser.add_argument(
        "--task",
        default="",
        type=str,
        help="which task to do",
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
        "--delimiter",
        default=',',
        type=str,
        help="The path of the processed data",
    )
    parser.add_argument(
        "--head",
        action="store_true",
        help="The path of the processed data",
    )
    parser.add_argument(
        "--dev",
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

    if args.task == "generate_prom":
        Generate_prom_train_dev(args)
    elif  args.task == "shuffle":
        Shuffle(args)
    elif  args.task == "find_train":
        Find_train(args)
    elif  args.task == "prom_1000":
        Process_1000(args)
    elif  args.task == "prom_1000_kmer":
        Process_1000_kmer(args)
    elif  args.task == "splice":
        Process_splice(args)
    elif  args.task == "pair":
        Process_pair(args)
    elif  args.task == "p53":
        Process_p53(args)
    elif  args.task == "p53_mut":
        Process_p53_mut(args)
    elif  args.task == "sep_p53":
        Seperate_p53(args)
    elif  args.task == "690":
        Process_690(args)
    elif  args.task == "mouse":
        Process_mouse(args)
    elif  args.task == "prom-core":
        Process_prom_core(args)
    else:
        Process(args)
   
    

    


if __name__ == "__main__":
    main()
