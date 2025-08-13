import os

data_dir = "/infodev1/non-phi-data/junjiang/DNABERT/our_data"

reference_data = ["DNAbert_input_reference.benign.txt", "DNAbert_input_reference.pathogenic.txt"]
mutant_data = ["DNAbert_input_mutant.benign.txt", "DNAbert_input_mutant.pathogenic.txt"]

ref_b_lines = open(os.path.join(data_dir, reference_data[0]), 'r').readlines()
print("reference.benign size: %d" % len(ref_b_lines))
ref_p_lines = open(os.path.join(data_dir, reference_data[1]), 'r').readlines()
print("reference.pathogenic size: %d" % len(ref_p_lines))
mut_b_lines = open(os.path.join(data_dir, mutant_data[0]), 'r').readlines()
print("mutant.benign size: %d" % len(mut_b_lines))
mut_p_lines = open(os.path.join(data_dir, mutant_data[1]), 'r').readlines()
print("mutant.pathogenic size: %d" % len(mut_p_lines))

ref_b_lines_labels = len(ref_b_lines) * [0]
ref_p_lines_labels = len(ref_p_lines) * [1]
mut_b_lines_labels = len(mut_b_lines) * [2]
mut_p_lines_labels = len(mut_p_lines) * [3]


def add_label_to_lines(lines, label_list, sep="\t"):
    new_lines = []
    for idx, l in enumerate(lines):
        new_l = l.strip() + sep + str(label_list[idx]) + "\n"
        new_lines.append(new_l)
    return new_lines


def write_lines_to_file(lines, wrt_to):
    f_path = os.path.split(wrt_to)[0]
    if not os.path.exists(f_path):
        os.makedirs(f_path)
    with open(wrt_to, 'w') as fp:
        for l in lines:
            fp.write(l)


all_lines = ref_b_lines + ref_p_lines + mut_b_lines + mut_p_lines
all_labels = ref_b_lines_labels + ref_p_lines_labels + mut_b_lines_labels + mut_p_lines_labels

all_lines_w_labels = add_label_to_lines(all_lines, all_labels)
print(all_lines_w_labels[0])
save_to = os.path.join(data_dir, "DNAbert_input_all.txt")
write_lines_to_file(all_lines_w_labels, save_to)




