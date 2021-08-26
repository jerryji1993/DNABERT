#### ::: DNABERT-viz find motifs ::: ####

import os
import pandas as pd
import numpy as np
import argparse
import motif_utils as utils


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data_dir",
        default=None,
        type=str,
        required=True,
        help="The input data dir. Should contain the sequence+label .tsv files (or other data files) for the task.",
    )

    parser.add_argument(
        "--predict_dir",
        default=None,
        type=str,
        required=True,
        help="Path where the attention scores were saved. Should contain both pred_results.npy and atten.npy",
    )

    parser.add_argument(
        "--window_size",
        default=24,
        type=int,
        help="Specified window size to be final motif length",
    )

    parser.add_argument(
        "--min_len",
        default=5,
        type=int,
        help="Specified minimum length threshold for contiguous region",
    )

    parser.add_argument(
        "--pval_cutoff",
        default=0.005,
        type=float,
        help="Cutoff FDR/p-value to declare statistical significance",
    )

    parser.add_argument(
        "--min_n_motif",
        default=3,
        type=int,
        help="Minimum instance inside motif to be filtered",
    )

    parser.add_argument(
        "--align_all_ties",
        action='store_true',
        help="Whether to keep all best alignments when ties encountered",
    )

    parser.add_argument(
        "--save_file_dir",
        default='.',
        type=str,
        help="Path to save outputs",
    )

    parser.add_argument(
        "--verbose",
        action='store_true',
        help="Verbosity controller",
    )

    parser.add_argument(
        "--return_idx",
        action='store_true',
        help="Whether the indices of the motifs are only returned",
    )

    # TODO: add the conditions
    args = parser.parse_args()

    atten_scores = np.load(os.path.join(args.predict_dir,"atten.npy"))
    pred = np.load(os.path.join(args.predict_dir,"pred_results.npy"))
    dev = pd.read_csv(os.path.join(args.data_dir,"dev.tsv"),sep='\t',header=0)
    dev.columns = ['sequence','label']
    dev['seq'] = dev['sequence'].apply(utils.kmer2seq)
    dev_pos = dev[dev['label'] == 1]
    dev_neg = dev[dev['label'] == 0]
    pos_atten_scores = atten_scores[dev_pos.index.values]
    neg_atten_scores = atten_scores[dev_neg.index.values]
    assert len(dev_pos) == len(pos_atten_scores)

    # run motif analysis
    merged_motif_seqs = utils.motif_analysis(dev_pos['seq'],
                                        dev_neg['seq'],
                                        pos_atten_scores,
                                        window_size = args.window_size,
                                        min_len = args.min_len,
                                        pval_cutoff = args.pval_cutoff,
                                        min_n_motif = args.min_n_motif,
                                        align_all_ties = args.align_all_ties,
                                        save_file_dir = args.save_file_dir,
                                        verbose = args.verbose,
                                        return_idx  = args.return_idx
                                    )

if __name__ == "__main__":
    main()


