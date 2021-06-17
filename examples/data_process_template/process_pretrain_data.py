import argparse
import random
import numpy as np
     

def cut_no_overlap(length, kmer=1, max_prob=0.5):
    cuts = []
    while length:
        if length <= 509+kmer:
            cuts.append(length)
            break
        else:
            if random.random() > max_prob:
                cut = max(int(random.random()*(509+kmer)), 5)
            else:
                cut = 509+kmer
            cuts.append(cut)
            length -= cut

    return cuts


def sampling(length, kmer=1, sampling_rate=1):
    times = int(length*sampling_rate/256)
    starts = []
    ends = []
    for i in range(times):
        cut = max(int(random.random()*(509+kmer)), 5)
        start = np.random.randint(length-kmer)
        starts.append(start)
        ends.append(start+cut)
    
    return starts, ends


def sampling_fix(length, kmer=1, sampling_rate=1, fix_length=10245):
    times = int(length*sampling_rate/fix_length)
    starts = []
    ends = []
    for i in range(times):
        cut = fix_length
        start = np.random.randint(length-6-fix_length)
        starts.append(start)
        ends.append(start+cut)
    
    return starts, ends


def get_kmer_sentence(original_string, kmer=1, stride=1):
    if kmer == -1:
        return original_string

    sentence = ""
    original_string = original_string.replace("\n", "")
    i = 0
    while i < len(original_string)-kmer:
        sentence += original_string[i:i+kmer] + " "
        i += stride
    
    return sentence[:-1].strip("\"")



def get_kmer_sequence(original_string, kmer=1):
    if kmer == -1:
        return original_string

    sequence = []
    original_string = original_string.replace("\n", "")
    for i in range(len(original_string)-kmer):
        sequence.append(original_string[i:i+kmer])
    
    sequence.append(original_string[-kmer:])
    return sequence

def Process(args):
    old_file = open(args.file_path, "r")
    if args.output_path == None:
        args.output_path = args.file_path

    if args.sampling_rate!=1.0:
        new_file_path = args.output_path + "_sam" + str(args.kmer)
    else:
        new_file_path = args.output_path + "_cut" + str(args.kmer)
    new_file = open(new_file_path, "w")
    line = old_file.readline()
    while line:
        line_length = len(line)
        if args.sampling_rate != 1.0:
            starts, ends = sampling_fix(length=line_length, kmer=args.kmer, sampling_rate=args.sampling_rate, fix_length=args.length)
            for i in range(len(starts)):
                new_line = line[starts[i]:ends[i]]
                sentence = get_kmer_sentence(new_line, kmer=args.kmer)
                new_file.write(sentence + "\n")
            
        else:
            cuts = cut_no_overlap(length=line_length, kmer=args.kmer)
            start = 0
            for cut in cuts:
                new_line = line[start:start+cut]
                sentence = get_kmer_sentence(new_line, kmer=args.kmer)
                start += cut
                new_file.write(sentence + "\n")
                
        line = old_file.readline()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sampling_rate", 
        default=1.0,
        type=float,
        help="We will sample sampling_rate*total_length*2/512 times",
    )
    parser.add_argument(
        "--kmer",
        default=1,
        type=int,
        help="K-mer",
    )
    parser.add_argument(
        "--length",
        default=10000,
        type=int,
        help="Length of the sampled sequence",
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
