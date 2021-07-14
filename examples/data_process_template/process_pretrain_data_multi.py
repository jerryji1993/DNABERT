from multiprocessing import Pool
import copy
import argparse

from process_pretrain_data import Process

# filenames = ['xaa', 'xab', 'xac', 'xad', 'xae', 'xaf', 'xag', 'xah', 'xai', 'xaj', 'xak', 'xal', 'xam', 'xan', 'xao', 'xap', 'xaq', 'xar', 'xas', 'xat', 'xau', 'xav', 'xaw']
# filenames = ['xaa', 'xab']

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
        default="/home/zhihan/dna/data/split/",
        type=str,
        help="The path of the file to be processed",
    )

    args = parser.parse_args()

    # multiprocess
    p = Pool(22)

    for i in range(1,23):
        arg_new = copy.deepcopy(args)
        arg_new.file_path = "/root/data/genome/" + "GRCh38.chr" + str(i) + ".fa"
        arg_new.output_path = "/root/data/sub_001_6140/" + "GRCh38.chr" + str(i) + ".fa"
        # arg_new.file_path = arg_new.output_path + filename
        p.apply_async(Process, args=(arg_new,))
    
    p.close()
    p.join()




if __name__ == "__main__":
  main()
