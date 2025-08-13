# created by Chuanyi, modified by Jun on infodev server.
import pickle

import torch
from transformers import BertModel, BertConfig, DNATokenizer
import umap, os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Pool, Pipe, Process
import torch


matplotlib.use('Agg')

dir_to_pretrained_model = "/infodev1/non-phi-data/junjiang/DNABERT/model/6-new-12w-0"
path_to_config = 'https://raw.githubusercontent.com/jerryji1993/DNABERT/master/src/transformers/dnabert-config/bert-config-6/config.json'
test_file = "/infodev1/non-phi-data/junjiang/DNABERT/our_data/DNAbert_input_all.txt"  # created using data_intergration.py

output_dir = "/infodev1/non-phi-data/junjiang/DNABERT/our_data/output"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
output_png = os.path.join(output_dir, "test_our_input_all.png")

pickle_fn = os.path.join(output_dir, "embedding_all.pickel")

if os.path.exists(pickle_fn):
    fp = open(pickle_fn, 'rb')
    pooled_encodings, labels = pickle.load(fp)
    fp.close()
else:
    config = BertConfig.from_pretrained(path_to_config)
    tokenizer = DNATokenizer.from_pretrained('dna6')

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(device)
    model = BertModel.from_pretrained(dir_to_pretrained_model, config=config)
    model.to(device)


    lines = open(test_file, 'r').readlines()


    def get_model_input(TCGA_label_line):
        seq, label = TCGA_label_line.strip().split('\t')
        model_input = tokenizer.encode_plus(seq,
                                            add_special_tokens=True,
                                            max_length=512,
                                            pad_to_max_length=True)["input_ids"]
        return [model_input, label]


    pool = Pool(processes=20)
    results_async = pool.map_async(func=get_model_input, iterable=lines)
    results = results_async.get()

    model_input = list(np.array(results)[:, 0])
    labels = list(np.array(results)[:, 1])

    model_input = torch.tensor(model_input, dtype=torch.long)

    batch_size = 8

    pooled_encodings = []
    for i in range(0, model_input.shape[0], batch_size):
        output = model(model_input[i:(i + batch_size), :])
        pooled_encodings.append(output[1].detach().numpy())
        print("Progress: %d/%d" % (i, model_input.shape[0]))

    pooled_encodings = np.vstack(pooled_encodings)
    labels = np.array(labels).astype(int)
    data = (pooled_encodings, labels)
    file = open(pickle_fn, 'wb')
    pickle.dump(data, file)
    file.close()

##################################################
# draw UMAP
##################################################
plt.figure(figsize=(20, 20), dpi=300)
txt_labels = ["ref_ben", "ref_path", "mut_ben", "mut_path"]
colors = ['g', 'y', 'b', 'r']

ax_list = []
txt_list = []
lde = umap.UMAP().fit_transform(pooled_encodings)
for x in np.unique(labels):
    txt_l = txt_labels[int(x)]
    ax = plt.scatter(
         lde[labels == x, 0],
         lde[labels == x, 1],
         s=5,
         alpha=0.3,
         c=colors[int(x)],
         label=txt_l)
    ax_list.append(ax)
    txt_list.append(txt_l)
x_all = lde[:, 0]
y_all = lde[:, 1]
x_lim = (min(x_all), max(x_all))
y_lim = (min(y_all), max(y_all))


from matplotlib.lines import Line2D
legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[0], label=txt_labels[0]),
                Line2D([0], [0], marker='o',  color='w', markerfacecolor=colors[1], label=txt_labels[1]),
                Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[2], label=txt_labels[2]),
                Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[3], label=txt_labels[3])]
plt.legend(handles=legend_elements,  prop={'size': 18})
plt.xlim(x_lim)
plt.ylim(y_lim)
plt.savefig(output_png, dpi=300)

##################################################
# calculate performance
##################################################
from sklearn.model_selection import KFold
from sklearn.utils import shuffle
from sklearn import svm
from sklearn.metrics import f1_score, accuracy_score

labels = np.array(labels).astype(int)
pooled_encodings, labels = shuffle(pooled_encodings, labels)

clf = svm.SVC()
kf = KFold(n_splits=5)
for train_index, test_index in kf.split(pooled_encodings):
    print("TRAIN:", train_index, "TEST:", test_index)
    X_train, X_test = pooled_encodings[train_index, :], pooled_encodings[test_index, :]
    y_train, y_test = labels[train_index], labels[test_index]
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    f1 = f1_score(y_test, y_pred, average='micro')
    acc_score = accuracy_score(y_test, y_pred)
    print("F1 score: %.3f, accuracy: %.3f" %(f1, acc_score))


print("merge labels to two categories")

# reference.benign = 0; reference.pathogenic = 1; mutant.benign = 2; mutant.pathogenic = 3
labels[labels == 2] = 0
labels[labels == 3] = 1

for train_index, test_index in kf.split(pooled_encodings):
    print("TRAIN:", train_index, "TEST:", test_index)
    X_train, X_test = pooled_encodings[train_index, :], pooled_encodings[test_index, :]
    y_train, y_test = labels[train_index], labels[test_index]
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    f1 = f1_score(y_test, y_pred, average='micro')
    acc_score = accuracy_score(y_test, y_pred)
    print("F1 score: %.3f, accuracy: %.3f" %(f1, acc_score))
