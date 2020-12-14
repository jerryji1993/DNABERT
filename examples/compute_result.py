import argparse
import numpy as np
import csv
from copy import deepcopy
from sklearn.metrics import matthews_corrcoef, confusion_matrix, f1_score

def generate_pred(predict_results, i, slide, metric="max"):

    results = predict_results[i*3:(i+1)*3]

    if metric == "max":
        pred = max(results)
    elif metric == "mean":
        pred = np.mean(results)
    elif metric == "second-max":
        pred = np.sort(results)[-2]
    else:
        pass

    return pred

def Compute_scan(args):
    predict_results = np.load(args.pred_path) 
    labels = np.load(args.label_path)
    labels = list(labels.astype(int))
    
    results = []
    for i in range(len(labels)):
        pred = generate_pred(predict_results, i, args.slide, args.metric)
        
        if pred >= args.bound:
            results.append(1)
        else:
            results.append(0)
    a = set(results)
    b = set(labels)
    f1 = f1_score(y_true=labels, y_pred=results)
    mcc = matthews_corrcoef(labels, results)
    tn, fp, fn, tp = confusion_matrix(labels, results).ravel()

    count = 0
    for i in range(len(results)):
        if results[i] == labels[i]:
            count+=1

    print("number of examples: " + str(len(labels)))
    print("number of positive examples: " + str(sum(labels)))
    print("number of negative examples: " + str(len(labels)-sum(labels)))
    print("f1: ", str(f1))
    print("mcc: " + str(mcc))
    print("accuracy: " + str(float(count)/len(results)))
    print("tn:" + str(tn))
    print("fp:" + str(fp))
    print("fn:" + str(fn))
    print("tp:" + str(tp))


def Compute_mouse(args):
    result_file = open(args.pred_path, "r")
    results = result_file.readlines()
    print(len(results))

    all_preds = []
    current_preds = []
    for result in results:
        scores = result.split()
        scores = [scores[0], float(scores[1]), float(scores[2]), float(scores[3]), float(scores[4]), float(scores[5]), float(scores[6]), float(scores[7])]
        if current_preds == [] or scores[0] == current_preds[0][0]:
            current_preds.append(scores)
        else:
            all_preds.append(current_preds)
            current_preds = []
            current_preds.append(scores)
    all_preds.append(current_preds)
    
    print("Number of task: %d" % len(all_preds))

    def get_acc(val):
        return val[1]

    def get_auc(val):
        return val[2]

    tasks = []
    acc = []
    auc = []
    aupr = []
    f1 = []
    mcc = []
    precision = []
    recall = []

    for pred in all_preds:
        if len(pred) < 10 :
            print("Short %s : %d" % (pred[0][0], len(pred)))

        if args.index == "acc":
            pred.sort(key=get_acc)
        elif args.index == "auc":
            pred.sort(key=get_auc)
        else:
            raise ValueError()
        
        BEST = -1
        for i in range(len(pred)):
            if pred[i][1] == pred[-1][1] and pred[i][2] > pred[-1][2]:
                BEST = deepcopy(i)
        tasks.append(pred[0][0])

        best_pred = pred[BEST]
        acc.append(best_pred[1])
        auc.append(best_pred[2])
        aupr.append(best_pred[3])
        f1.append(best_pred[4])
        mcc.append(best_pred[5])
        precision.append(best_pred[6])
        recall.append(best_pred[7])

    acc_ave = np.mean(acc)
    auc_ave = np.mean(auc)
    aupr_ave = np.mean(aupr)
    f1_ave = np.mean(f1)
    mcc_ave = np.mean(mcc)
    precision_ave = np.mean(precision)
    recall_ave = np.mean(recall)


    print("acc: " + str(acc_ave))
    print("auc: " + str(auc_ave))
    print("aupr: " + str(aupr_ave))
    print("f1: ", str(f1_ave))
    print("mcc: " + str(mcc_ave))
    print("precision: ", str(precision_ave))
    print("recall: " + str(recall_ave))

    # find and print the tasks whose results are worst
    ranks = np.argsort(auc)[:args.num_worst]
    print("Top %d worst tasks: " % (args.num_worst))
    for i in ranks:
        print(tasks[i] + "  %3f  %3f" % (acc[i], auc[i]))




def Compute_690(args):
    result_file = open(args.pred_path, "r")
    results = result_file.readlines()

    preds = []

    for result in results:
        scores = result.split()
        preds.append([scores[0], float(scores[1]), float(scores[2]), float(scores[4]), float(scores[5])])
    
    num_results = args.num_results

    num_example = int(len(preds)/num_results)
    print("Num of tasks: %d" % num_example)

    def get_acc(val):
        return val[1]

    def get_auc(val):
        return val[2]
    
    def get_f1(val):
        return val[3]

    def get_mcc(val):
        return val[4]
    
    tasks = []
    acc = []
    auc = []
    f1 = []
    mcc = []

    for i in range(num_example):
        tasks.append(preds[i*num_results][0])

        current_preds = preds[i*num_results:(i+1)*num_results]
        if args.index == "acc":
            current_preds.sort(key=get_acc)
        elif args.index == "auc":
            current_preds.sort(key=get_auc)
        elif args.index == "f1":
            current_preds.sort(key=get_f1)
        elif args.index == "mcc":
            current_preds.sort(key=get_mcc)
        else:
            raise ValueError()
        best_pred = current_preds[-1]
        acc.append(best_pred[1])
        auc.append(best_pred[2])
        f1.append(best_pred[3])
        mcc.append(best_pred[4])

    # calculate and print the average scores
    acc_ave = np.mean(acc)
    auc_ave = np.mean(auc)
    f1_ave = np.mean(f1)
    mcc_ave = np.mean(mcc)


    print("acc: " + str(acc_ave))
    print("auc: " + str(auc_ave))
    print("f1: ", str(f1_ave))
    print("mcc: " + str(mcc_ave))

    # find and print the tasks whose results are worst
    ranks = np.argsort(auc)[:args.num_worst]
    print("Top %d worst tasks: " % (args.num_worst))
    for i in ranks:
        print(tasks[i] + "  %3f  %3f" % (acc[i], auc[i]))
    


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bound",
        default=0.5,
        type=float,
        help="K-mer",
    )
    parser.add_argument(
        "--pred_path",
        default=None,
        type=str,
        help="The path of the predicted result",
    )
    parser.add_argument(
        "--label_path",
        default=None,
        type=str,
        help="The path of the label",
    )
    parser.add_argument(
        "--metric",
        default="max",
        type=str,
        help="The metric of computing predited result (scan)",
    )
    parser.add_argument(
        "--slide",
        default=3,
        type=int,
        help="How many 500s to use for the predictes result of 1000 (scan)",
    )
    parser.add_argument(
        "--task",
        default="scan",
        type=str,
        help="Which task to compute result",
    )
    parser.add_argument(
        "--index",
        default="acc",
        type=str,
        help="Which index to sort result (690)",
    )
    parser.add_argument(
        "--num_results",
        default="10",
        type=int,
        help="Number of results for each task (690)",
    )
    parser.add_argument(
        "--num_worst",
        default="10",
        type=int,
        help="Number of worst tasks to print out (690)",
    )

    args = parser.parse_args()

    if args.task == "scan":
        Compute_scan(args)
    elif args.task == "690":
        Compute_690(args)
    elif args.task == "mouse":
        Compute_mouse(args)
    else:
        raise ValueError()

    


if __name__ == "__main__":
    main()
