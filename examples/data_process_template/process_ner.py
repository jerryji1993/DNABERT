import argparse
import csv
import os
import h5py
import numpy as np
import random
from process_pretrain_data import get_kmer_sequence
from multiprocessing import Pool


def generate_example(X, Y, kmer, index):
    # assert X.shape[0] == Y.shape[0]
    lines = []
    for j in range(len(X)):
        if j % 1000 == 0:
            print("%s : %s" % (index, j))

        label = list(np.zeros(200,dtype=int)) + list(np.where(Y[j]==1)[1]) + list(np.zeros(201-kmer,dtype=int))

        sequence = get_kmer_sequence(X[j].decode("utf-8"), kmer)
        lines.append([sequence, label])
    
    return lines


def Process(args):
    filename = args.file_path
    h5 = h5py.File(filename, "r")
    num_chunks = len(h5.keys())//2
    keys = list(h5.keys())[:num_chunks]


    X = []

    for i, key in enumerate(keys):
        x_key = key
        y_key = x_key.replace("X","Y")

        X_l = h5[x_key]
        Y_l = h5[y_key][0]

        X.extend(X_l)
        
        if i == 0:
            Y = Y_l
        else:
            Y = np.concatenate([Y, Y_l], axis=0)

        print("%d : %d, %d, %s" % (i, len(X), Y.shape[0], str(key)))
    
    print(len(X))
    print(len(Y))
    
    n_proc = int(args.n_process)
    print("number of processes for converting feature: " + str(n_proc))
    p = Pool(n_proc)
    indexes = [0]
    len_slice = int(len(X)/n_proc)
    for i in range(1, n_proc+1):
        if i != n_proc:
            indexes.append(len_slice*(i))
        else:
            indexes.append(len(X))
    
    results = []
    
    for i in range(n_proc):
        results.append(p.apply_async(generate_example, args=(X[indexes[i]:indexes[i+1]], Y[indexes[i]:indexes[i+1]], args.kmer, i)))
        print(str(i+1) + ' processor started !')
    
    p.close()
    p.join()

    lines = []
    for result in results:
        lines.extend(result.get())

    
    path = "/".join(args.file_path.split('/')[:-1]) + "/" + str(args.kmer) + "/train.txt"
    print(path)
    file = open(path, "w")
    for line in lines:
        for k, word in enumerate(line[0]):
            file.write(str(word) + " " + str(line[1][k]) + "\n")
        file.write("\n")

    




        

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--kmer",
        default=1,
        type=int,
        help="K-mer",
    )
    parser.add_argument(
        "--n_process",
        default=24,
        type=int,
        help="Number of processes for data processing",
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
    args = parser.parse_args()

    Process(args)
    

    


if __name__ == "__main__":
    main()



