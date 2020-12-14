import torch
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import numpy as np

from transformers import BertTokenizer, BertModel, DNATokenizer
from process_pretrain_data import get_kmer_sentence


def format_attention(attention):
    squeezed = []
    for layer_attention in attention:
        # 1 x num_heads x seq_len x seq_len
        if len(layer_attention.shape) != 4:
            raise ValueError("The attention tensor does not have the correct number of dimensions. Make sure you set "
                             "output_attentions=True when initializing your model.")
        squeezed.append(layer_attention.squeeze(0))
    # num_layers x num_heads x seq_len x seq_len
    return torch.stack(squeezed)

def get_attention_dna(model, tokenizer, sentence_a, start, end):
    inputs = tokenizer.encode_plus(sentence_a, sentence_b=None, return_tensors='pt', add_special_tokens=True)
    input_ids = inputs['input_ids']
    attention = model(input_ids)[-1]
    input_id_list = input_ids[0].tolist() # Batch index 0
    tokens = tokenizer.convert_ids_to_tokens(input_id_list) 
    attn = format_attention(attention)
    attn_score = []
    for i in range(1, len(tokens)-1):
        attn_score.append(float(attn[start:end+1,:,0,i].sum()))
    return attn_score

def get_real_score(attention_scores, kmer, metric):
    counts = np.zeros([len(attention_scores)+kmer-1])
    real_scores = np.zeros([len(attention_scores)+kmer-1])

    if metric == "mean":
        for i, score in enumerate(attention_scores):
            for j in range(kmer):
                counts[i+j] += 1.0
                real_scores[i+j] += score

        real_scores = real_scores/counts
    else:
        pass

    return real_scores

SEQUENCE = "TGCCTGGCTTTTTGTAATTTTTGAAGAGACGGGGTTTTGCCATGATG"

def Visualize(args):
    if args.kmer == 0:
        KMER_LIST = [3,4,5,6]

        for kmer in KMER_LIST:
            tokenizer_name = 'dna' + str(kmer)
            model_path = os.path.join(args.model_path, str(kmer))
            model = BertModel.from_pretrained(model_path, output_attentions=True)
            tokenizer = DNATokenizer.from_pretrained(tokenizer_name, do_lower_case=False)
            raw_sentence = args.sequence if args.sequence else SEQUENCE
            sentence_a = get_kmer_sentence(raw_sentence, kmer)
            tokens = sentence_a.split()

            attention = get_attention_dna(model, tokenizer, sentence_a, start=args.start_layer, end=args.end_layer)
            attention_scores = np.array(attention).reshape(np.array(attention).shape[0],1)
            # attention_scores[0] = 0
            
            real_scores = get_real_score(attention_scores, kmer, args.metric)
            real_scores = real_scores / np.linalg.norm(real_scores)

            if kmer != KMER_LIST[0]:
                scores += real_scores.reshape(1, real_scores.shape[0])
            else:
                scores = real_scores.reshape(1, real_scores.shape[0])

    else:
        # load model and calculate attention
        tokenizer_name = 'dna' + str(args.kmer)
        model_path = args.model_path
        model = BertModel.from_pretrained(model_path, output_attentions=True)
        tokenizer = DNATokenizer.from_pretrained(tokenizer_name, do_lower_case=False)
        raw_sentence = args.sequence if args.sequence else SEQUENCE
        sentence_a = get_kmer_sentence(raw_sentence, args.kmer)
        tokens = sentence_a.split()

        attention = get_attention_dna(model, tokenizer, sentence_a, start=args.start_layer, end=args.end_layer)
        attention_scores = np.array(attention).reshape(np.array(attention).shape[0],1)
        # attention_scores[0] = 0
        
        real_scores = get_real_score(attention_scores, args.kmer, args.metric)
        scores = real_scores.reshape(1, real_scores.shape[0])
    
    ave = np.sum(scores)/scores.shape[1]
    print(ave)
    print(scores)

    # plot        
    sns.set()
    ax = sns.heatmap(scores, cmap='YlGnBu', vmin=0)
    plt.show()

    


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--kmer",
        default=0,
        type=int,
        help="K-mer",
    )
    parser.add_argument(
        "--model_path",
        default="/home/zhihan/dna/dna-transformers/examples/ft/690/p53-small/TAp73beta/3/",
        type=str,
        help="The path of the finetuned model",
    )
    parser.add_argument(
        "--start_layer",
        default=11,
        type=int,
        help="Which layer to start",
    )
    parser.add_argument(
        "--end_layer",
        default=11,
        type=int,
        help="which layer to end",
    )
    parser.add_argument(
        "--metric",
        default="mean",
        type=str,
        help="the metric used for integrate predicted kmer result to real result",
    )
    parser.add_argument(
        "--sequence",
        default=None,
        type=str,
        help="the sequence for visualize",
    )

    args = parser.parse_args()
    Visualize(args)
        


if __name__ == "__main__":
    main()