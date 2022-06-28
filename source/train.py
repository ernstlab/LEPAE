import argparse, math, sys, scipy.sparse, gzip, random, numpy as np, joblib, linecache, pandas as pd, itertools, copy
from collections import OrderedDict

from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import roc_curve, auc, average_precision_score, mean_squared_error

import torch
import torch.nn as nn

from shared import *

NUM_WORKERS = 4
FIXED_SEED = 413666


def setDecisionTreeHyperparams(args):
    random.seed(args.seed)
    hyperparam = OrderedDict()
    hyperparam["max_depth"] = random.choice([16, 32, 64, 128, 256]) if args.random_search else args.max_depth
    hyperparam["min_samples_split"] = (
        random.choice([0.0005, 0.001, 0.002, 0.005, 0.01]) if args.random_search else args.min_samples_split
    )
    hyperparam["min_samples_leaf"] = (
        random.choice([0.0005, 0.001, 0.002, 0.005, 0.01]) if args.random_search else args.min_samples_leaf
    )
    hyperparam["max_features"] = "auto"
    return hyperparam


def setNeuralNetworkHyperparams(args):

    # Choose a random number of neurons
    def randomNumneuron():
        return random.choice([32, 64, 128])

    random.seed(args.seed)
    hyperparam = OrderedDict()
    hyperparam["batch_size"] = random.choice([16, 32, 64]) if args.random_search else args.batch_size
    hyperparam["learning_rate"] = random.choice([1e-8, 1e-6, 1e-4]) if args.random_search else args.learning_rate
    hyperparam["dropout_rate"] = random.choice([0, 0.25, 0.5]) if args.random_search else args.dropout_rate
    hyperparam["num_neuron_init"] = randomNumneuron() if args.random_search else args.num_neuron_init
    hyperparam["num_neuron_fin"] = randomNumneuron() if args.random_search else args.num_neuron_fin
    return hyperparam


def trainDT(args, training_inputs, training_labels, tuning_inputs, tuning_labels, fout):

    # Define hyperparameters
    hyperparam = setDecisionTreeHyperparams(args)

    # Train a decision tree classifier
    clf = DecisionTreeClassifier(
        random_state=args.seed,
        max_depth=hyperparam["max_depth"],
        min_samples_split=hyperparam["min_samples_split"],
        min_samples_leaf=hyperparam["min_samples_leaf"],
        max_features=hyperparam["max_features"],
    )
    clf.fit(training_inputs, training_labels)
    # Evaluate and print result
    training_result = eval(training_labels, predictDT(clf, training_inputs))
    tuning_result = eval(tuning_labels, predictDT(clf, tuning_inputs))
    cols = [args.seed] + ["%.6f" % s for s in training_result + tuning_result] + list(hyperparam.values())
    fout.write("\t".join([str(c) for c in cols]) + "\n")
    fout.flush()

    return clf, hyperparam, tuning_result[1]


def trainNN(
    args,
    training_upstream_lines_to_read,
    training_downstream_lines_to_read,
    training_labels,
    tuning_inputs,
    tuning_labels,
    fout,
):

    # Define hyperparameters
    hyperparam = setNeuralNetworkHyperparams(args)

    # Set up a Siamese neural network
    net = SiameseNet(args.num_features, hyperparam)
    optimizer = torch.optim.SGD(net.parameters(), lr=hyperparam["learning_rate"], momentum=0.9)
    criterion = nn.BCELoss()  # binary cross entropy loss
    optimizer.zero_grad()

    # Determine which lines in which files to read in advance
    bs = hyperparam["batch_size"]
    num_batch = int(np.floor(len(training_upstream_lines_to_read) / bs) + 1)

    # Start training
    losses = []
    aurocs = []
    best_net = net
    # lv = 0
    for epoch in range(args.max_num_epoch):  # each training epoch
        net.train()  # set to train mode

        line_index = 0
        for j in range(num_batch):  # each batch
            current_batch_size = (
                bs if ((j + 1) * bs) < len(training_upstream_lines_to_read) else len(training_upstream_lines_to_read) % bs
            )
            if current_batch_size == 0:
                break

            x, y = readBatchTorch(
                current_batch_size,
                args.num_features,
                args.frac_feature,
                training_upstream_lines_to_read,
                training_downstream_lines_to_read,
                training_labels,
                line_index,
            )

            line_index += current_batch_size

            # Forward + backward + optimize
            outputs = net(x)  # forward
            loss = criterion(outputs, y)  # calculate BCE loss
            # lv += loss.item()
            # training_loss.append(loss.item())
            loss.backward()  # backward
            optimizer.step()  # optimize

            if args.debug:
                for k in range(current_batch_size):
                    fout.write("#debug#\t%d\t|\t" % y[k])
                    fout.write(" ".join([str(int(s)) for s in x[k, :]]) + "\n")

        net.eval()  # set to evaluation mode

        # Evaluate and print result
        training_scores, _ = predictNNFromFiles(
            args, net, training_upstream_lines_to_read, training_downstream_lines_to_read, training_labels
        )
        training_result = eval(training_labels, training_scores)
        tuning_result = eval(tuning_labels, net(torch.from_numpy(tuning_inputs.todense()).float()).data)
        current_loss = tuning_result[1]
        current_auroc = tuning_result[2]
        losses.append(current_loss)
        aurocs.append(current_auroc)
        cols = [args.seed] + ["%.6f" % s for s in training_result + tuning_result] + [epoch] + list(hyperparam.values())
        fout.write("\t".join([str(c) for c in cols]) + "\n")
        fout.flush()

        # Save the best performing neural network
        if args.save == True:
            if max(aurocs) == current_auroc:
                best_net = copy.deepcopy(net)

        # Check if early stopping is applicable
        if len(aurocs) >= 5:  # allow for five epochs before checking
            m = max(aurocs)  # best performance so far
            if m <= 0.50001:  # never better than random
                print("Early stopping due to poor AUROC performance")
                break
            if m >= aurocs[-1] and m >= aurocs[-2] and m >= aurocs[-3]:
                # no improvement than best in the last three epochs
                print("Early stopping due to no improvement in AUROC")
                break

    net.eval()  # set to evaluation mode

    return best_net, hyperparam, max(aurocs)


def main():
    ### Input arguments ###
    parent_parser = argparse.ArgumentParser(
        prog="base_maker", description="Train a classifier", fromfile_prefix_chars="@", add_help=False
    )

    parent_parser.add_argument("-i", "--training-data-filenames", help="paths to training data files", type=str, nargs="+")
    parent_parser.add_argument("-j", "--tuning-data-filenames", help="paths to tuning data files", type=str, nargs="+")
    parent_parser.add_argument(
        "-k", "--test-data-filenames", help="paths to test data files", type=str, nargs="+", default=[]
    )
    parent_parser.add_argument("-d", "--pairwise-dist", help="distance between windows", type=int, default=1000)
    parent_parser.add_argument("-w", "--window-size", help="genomic window size in bp", type=int, default=1000)
    parent_parser.add_argument("-o", "--output-filename-prefix", type=str, default="tmp")

    parent_parser.add_argument(
        "-a",
        "--positive-training-data-size",
        help="number of samples in positive training data",
        type=int,
        default=10000,
    )
    parent_parser.add_argument(
        "-b", "--positive-tuning-data-size", help="number of samples in positive tuning data", type=int, default=1000
    )
    parent_parser.add_argument(
        "--random-training-data",
        help="if sampling of training data should be dependent on input seed",
        action="store_true",
    )
    parent_parser.add_argument("-n", "--num-features", help="number of features in input vector", type=int, default=6512)
    parent_parser.add_argument(
        "--fixed-batch-size", help="batch size for evaluation (not a hyperparameter)", type=int, default=512
    )

    parent_parser.add_argument(
        "-m", "--random-search", help="if hyperparameters should be randomly set", action="store_true"
    )
    parent_parser.add_argument(
        "-v", "--save", help="if the trained model should be saved after training", action="store_true"
    )
    parent_parser.add_argument("-f", "--frac-feature", help="use fractional features instead of binary", action="store_true")
    parent_parser.add_argument("-u", "--debug", help="print intermediate results to debug", action="store_true")
    parent_parser.add_argument("-s", "--seed", help="random seed", type=int, default=1)

    parent_parser.add_argument("--max-num-epoch", help="maximum number of training epochs", type=int, default=10)

    main_parser = argparse.ArgumentParser()
    subparsers = main_parser.add_subparsers(help="type of classifier to train", dest="classifier", required=True)

    decision_tree_parser = subparsers.add_parser(
        "DT", help="Train a decision tree classifier", fromfile_prefix_chars="@", parents=[parent_parser]
    )
    # decision_tree_parser.add_argument("--hard-voting", help="make trees vote (decision tree)", action="store_true")
    # decision_tree_parser.add_argument("--num-tree", help="number of trees (decision tree)", type=int, default=5)
    decision_tree_parser.add_argument("--max-depth", help="maximum tree depth (decision tree)", type=int, default=32)
    decision_tree_parser.add_argument(
        "--min-samples-split",
        help="minimum number of samples in a node allowed for a split (decision tree)",
        type=float,
        default=2,
    )
    decision_tree_parser.add_argument(
        "--min-samples-leaf", help="minimum number of samples at a leaf (decision tree)", type=float, default=32
    )
    decision_tree_parser.add_argument(
        "--max-features", help="maximum features considered for split (decision tree)", type=str, default="auto"
    )
    # decision_tree_parser.add_argument(
    #     "--bootstrap", help="whether to bootstrap samples (decision tree)", action="store_true"
    # )

    neural_network_parser = subparsers.add_parser(
        "NN", help="Train a neural network classifier", fromfile_prefix_chars="@", parents=[parent_parser]
    )

    neural_network_parser.add_argument("--batch-size", help="batch size", type=int, default=128)
    neural_network_parser.add_argument("--learning-rate", help="learning_rate", type=float, default=0.1)
    neural_network_parser.add_argument("--dropout-rate", help="dropout rate", type=float, default=0.1)

    neural_network_parser.add_argument(
        "--num-layer-init", help="number of hidden layer in initial sub-networks", type=int, default=1
    )
    neural_network_parser.add_argument(
        "--num-layer-fin", help="number of hidden layer in final sub-network", type=int, default=1
    )
    neural_network_parser.add_argument(
        "--num-neuron-init",
        help="number of neurons in the initial sub-network",
        type=int,
        default=8,
    )
    neural_network_parser.add_argument(
        "--num-neuron-fin", help="number of neurons in the final sub-network", type=int, default=8
    )
    # neural_network_parser.add_argument(
    #     "--criteria",
    #     help="criteria for choosing the best classifier or epoch",
    #     type=str,
    #     default="auroc",
    #     choices=["auroc", "loss"],
    # )

    jaccard_parser = subparsers.add_parser(
        "Jaccard", help="Report Jaccard index (no training required)", parents=[parent_parser]
    )

    args = main_parser.parse_args()

    if args.classifier in ["DT", "NN"]:

        # Print out input arguments to progress file
        fout = open(args.output_filename_prefix + ".progress.txt", "w")
        for key, value in vars(args).items():
            if "max_" not in key and "min_" not in key:
                fout.write("# %s: %s\n" % (key, value))

        # Write header
        fout.write("seed\t")
        eval_output = ["loss", "mse", "auroc", "auprc", "pmean", "nmean"]
        fout.write("\t".join([a + b for (a, b) in list(itertools.product(["tr", "va"], eval_output))]) + "\t")
        if args.classifier == "DT":
            fout.write("\t".join(["max_depth", "min_samples_split", "min_samples_leaf", "max_features"]))
        elif args.classifier == "NN":
            fout.write(
                "\t".join(["epoch", "batch_size", "learning_rate", "dropout_rate", "num_neuron_init", "num_neuron_fin"])
            )
        fout.write("\n")

    # Set random seed
    random.seed(args.seed)
    np.random.seed(args.seed)

    training_upstream_lines_to_read, training_downstream_lines_to_read, training_labels = determineLinesToRead(
        args,
        args.training_data_filenames,
        args.positive_training_data_size,
        args.seed if args.random_training_data else FIXED_SEED,
    )
    if args.debug:
        fout.write("\n".join(["#debug#\tp\t" + str(s) for s in training_upstream_lines_to_read]) + "\n")
        fout.write("\n".join(["#debug#\tn\t" + str(s) for s in training_downstream_lines_to_read]) + "\n")
        fout.write("\n".join(["#debug#\tl\t" + str(s) for s in training_labels]) + "\n")

    tuning_inputs, tuning_labels = readDataCSR(
        args, *determineLinesToRead(args, args.tuning_data_filenames, args.positive_tuning_data_size)
    )
    if args.test_data_filenames:
        test_inputs, test_labels = readDataCSR(
            args, *determineLinesToRead(args, args.test_data_filenames, args.positive_tuning_data_size)
        )

    ### TRAIN ###
    if args.classifier == "DT":
        training_inputs, training_labels = readDataCSR(
            args, training_upstream_lines_to_read, training_downstream_lines_to_read, training_labels
        )
        clf, hyperparam, auroc = trainDT(args, training_inputs, training_labels, tuning_inputs, tuning_labels, fout)

        if args.save:
            joblib.dump(clf, args.output_filename_prefix + ".pkl")  # pickle model

            writeScores(
                "%s.training_pred.txt.gz" % args.output_filename_prefix,
                args.frac_feature,
                training_labels,
                predictDT(clf, training_inputs),
            )
            writeScores(
                "%s.tuning_pred.txt.gz" % args.output_filename_prefix,
                args.frac_feature,
                tuning_labels,
                predictDT(clf, tuning_inputs),
            )
            if args.test_data_filenames:
                writeScores(
                    "%s.test_pred.txt.gz" % args.output_filename_prefix,
                    args.frac_feature,
                    test_labels,
                    predictDT(clf, test_inputs),
                )

    elif args.classifier == "NN":
        torch.manual_seed(args.seed)
        torch.use_deterministic_algorithms(True)

        clf, hyperparam, auroc = trainNN(
            args,
            training_upstream_lines_to_read,
            training_downstream_lines_to_read,
            training_labels,
            tuning_inputs,
            tuning_labels,
            fout,
        )

        if args.save:
            joblib.dump(clf, args.output_filename_prefix + ".pkl")  # pickle model

            training_scores, _ = predictNNFromFiles(
                args, clf, training_upstream_lines_to_read, training_downstream_lines_to_read, training_labels, False
            )
            writeScores(
                "%s.training_pred.txt.gz" % args.output_filename_prefix,
                args.frac_feature,
                training_labels,
                training_scores,
            )
            writeScores(
                "%s.tuning_pred.txt.gz" % args.output_filename_prefix,
                args.frac_feature,
                tuning_labels,
                clf(torch.from_numpy(tuning_inputs.todense()).float()).data,
            )
            if args.test_data_filenames:
                writeScores(
                    "%s.test_pred.txt.gz" % args.output_filename_prefix,
                    args.frac_feature,
                    test_labels,
                    clf(torch.from_numpy(test_inputs.todense()).float()).data,
                    test_inputs if args.debug else [],
                )

    elif args.classifier == "Jaccard" and args.save:
        # training_inputs, training_labels = readDataCSR(
        #     args, training_upstream_lines_to_read, training_downstream_lines_to_read, training_labels
        # )

        # writeScores(
        #     "%s.training_pred.txt.gz" % args.output_filename_prefix,
        #     args.frac_feature,
        #     training_labels,
        #     computeJaccardIndex(args, training_inputs),
        # )
        # writeScores(
        #     "%s.tuning_pred.txt.gz" % args.output_filename_prefix,
        #     args.frac_feature,
        #     tuning_labels,
        #     computeJaccardIndex(args, tuning_inputs),
        # )
        if args.test_data_filenames:
            writeScores(
                "%s.test_pred.txt.gz" % args.output_filename_prefix,
                args.frac_feature,
                test_labels,
                computeJaccardIndex(args, test_inputs),
            )

    if args.classifier in ["DT", "NN"]:
        fout.close()


main()
