import argparse, math, sys, scipy.sparse, gzip, random, numpy as np, joblib, linecache, pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc, average_precision_score, mean_squared_error


def setTraininghyperparams(seed, num_features):
    random.seed(seed)
    max_depth = random.choice([32, 64, 128, 256])  # md (n=4)
    min_samples_split = random.choice([0.001, 0.002, 0.005, 0.01])  # mss (n=4)
    min_samples_leaf = random.choice([0.001, 0.002, 0.005, 0.01])  # msl (n=4)
    max_features = random.choice(["auto", "sqrt", "log2"])  # maf (n=3)
    bootstrap = True  # (n=1)
    return max_depth, min_samples_split, min_samples_leaf, max_features, bootstrap


def readData(filenames, pairwise_dist, window_size, num_windows, num_features, frac_feature):

    # 1. Define number of lines to read/skip
    num_lines = np.array([len(open(fn).readlines()) for fn in filenames])  # number of lines in each input file
    num_lines_to_sample = np.floor(
        num_windows * num_lines / np.sum(num_lines)
    )  # number of lines to sample from each input file
    while sum(num_lines_to_sample) < num_windows:  # adjust the number of lines to sample if needed
        num_lines_to_sample[random.randrange(len(filenames))] += 1
    num_lines_between_windows = int(
        pairwise_dist / window_size
    )  # number of lines to skip to define pairs with the specified pairwise distance
    print(
        "Number of lines between %d-kb windows that are %d kb apart: %d"
        % (
            int(window_size / 1000),
            int(pairwise_dist / 1000),
            num_lines_between_windows,
        )
    )

    # 2. Parse selected lines from input files
    num_examples = 0
    dt = float if frac_feature else bool  # data type of input features
    row, col, data = [], [], []  # row indices, column indices, feature values for later defining a sparse matrix
    for i in range(len(filenames)):
        fn = filenames[i]
        current_file_line_index = 0

        # Define randomly selected line indices to read in this file
        random_line_indices = list(
            random.sample(range(num_lines[i] - num_lines_between_windows), int(num_lines_to_sample[i]))
            if (num_lines[i] - num_lines_between_windows) > int(num_lines_to_sample[i])
            else range(num_lines[i] - num_lines_between_windows)
        )

        # Define shuffled line indices to later generate negative pairs
        shuffled_line_indices = np.random.permutation(len(random_line_indices))

        # Read each selected line and generate two positive pairs and two negative pairs
        for j in random_line_indices:
            line = linecache.getline(fn, 1 + j).strip().split("|")
            upstream_nonzero_feature_indices = [int(s) for s in line[1].split()] if len(line) > 1 else []
            upstream_feature_values = (
                [float(s) for s in line[2].split()]
                if (len(line) > 2 and frac_feature)
                else [True] * len(upstream_nonzero_feature_indices)
            )

            line = linecache.getline(fn, 1 + j + num_lines_between_windows).strip().split("|")
            downstream_nonzero_feature_indices = [int(s) for s in line[1].split()] if len(line) > 1 else []
            downstream_feature_values = (
                [float(s) for s in line[2].split()]
                if (len(line) > 2 and frac_feature)
                else [True] * len(downstream_nonzero_feature_indices)
            )

            num_nonzero_feature_indices = len(upstream_nonzero_feature_indices) + len(
                downstream_nonzero_feature_indices
            )

            # Positive pair
            row += [num_examples + current_file_line_index * 4] * num_nonzero_feature_indices
            col += upstream_nonzero_feature_indices + [num_features + s for s in downstream_nonzero_feature_indices]
            data += upstream_feature_values + downstream_feature_values

            # Flipped positive pair
            row += [num_examples + current_file_line_index * 4 + 1] * num_nonzero_feature_indices
            col += downstream_nonzero_feature_indices + [num_features + s for s in upstream_nonzero_feature_indices]
            data += downstream_feature_values + upstream_feature_values

            # Negative pair
            row += [num_examples + current_file_line_index * 4 + 2] * len(upstream_nonzero_feature_indices)
            col += upstream_nonzero_feature_indices
            data += upstream_feature_values

            row += [num_examples + shuffled_line_indices[current_file_line_index] * 4 + 2] * len(
                downstream_nonzero_feature_indices
            )
            col += [num_features + s for s in downstream_nonzero_feature_indices]
            data += downstream_feature_values

            # Flipped negative pair
            row += [num_examples + current_file_line_index * 4 + 3] * len(downstream_nonzero_feature_indices)
            col += downstream_nonzero_feature_indices
            data += downstream_feature_values

            row += [num_examples + shuffled_line_indices[current_file_line_index] * 4 + 3] * len(
                upstream_nonzero_feature_indices
            )
            col += [num_features + s for s in upstream_nonzero_feature_indices]
            data += upstream_feature_values

            current_file_line_index += 1

        num_examples += current_file_line_index * 4

    # 3. Define input data and positive/negative labels based on the data parsed above
    inputs = scipy.sparse.csr_matrix((data, (row, col)), shape=(num_examples, num_features * 2), dtype=dt)
    labels = [True, True, False, False] * int(num_examples / 4)

    return (inputs, labels)


def eval(clf, inputs, labels):
    scores = clf.predict_proba(inputs)
    scores = scores[:, 1]
    mse = mean_squared_error(labels, scores)
    fpr, tpr, _ = roc_curve(labels, scores)
    auroc = auc(fpr, tpr)  # au-ROC
    auprc = average_precision_score(labels, scores)  # au-PRC
    pos_mean_score = np.mean([scores[i] for i in range(len(scores)) if labels[i] == 1])
    neg_mean_score = np.mean([scores[i] for i in range(len(scores)) if labels[i] == 0])
    return mse, auroc, auprc, np.mean(scores), pos_mean_score, neg_mean_score


def main():
    ### Input arguments ###
    parser = argparse.ArgumentParser(
        prog="base_maker",
        description="Train a random forest",
        fromfile_prefix_chars="@",
    )

    # Input data filenames
    parser.add_argument(
        "-i",
        "--training-data-filenames",
        help="paths to training data files",
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "-j",
        "--tuning-data-filenames",
        help="paths to tuning data files",
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "-k",
        "--test-data-filenames",
        help="paths to test data files",
        type=str,
        nargs="+",
        default=[],
    )

    # Pairwise distance
    parser.add_argument(
        "-d",
        "--pairwise-dist",
        help="pairwise distance between windows",
        type=int,
        default=1000,
    )
    parser.add_argument("-w", "--window-size", help="genomic window size in bp", type=int, default=1000)

    # Output filename prefix
    parser.add_argument("-o", "--output-filename-prefix", type=str, default="tmp")

    # Data size
    parser.add_argument(
        "-a",
        "--positive-training-data-size",
        help="number of samples in positive training data",
        type=int,
        default=50000,
    )
    parser.add_argument(
        "-b",
        "--positive-tuning-data-size",
        help="number of samples in positive tuning data",
        type=int,
        default=5000,
    )

    # Options
    parser.add_argument(
        "-m",
        "--random-search",
        help="if hyperparameters should be randomly set",
        action="store_true",
    )
    parser.add_argument(
        "-v",
        "--save",
        help="if the trained model should be saved after training",
        action="store_true",
    )
    parser.add_argument(
        "-f",
        "--frac-feature",
        help="use fraction of overlap as features instead of binary",
        action="store_true",
    )

    # Hyperparameters for training
    # parser.add_argument('-r','--neg-data-ratio',help='ratio of negative samples to positive samples',type=int,default=1)
    parser.add_argument("-s", "--seed", help="random seed", type=int, default=1)
    parser.add_argument(
        "-r",
        "--random-search-n",
        help="number of hyperparameter combinations to try during random search",
        type=int,
        default=1,
    )
    # parser.add_argument('-z','--chunk-size',help='number of samples to read at a time before building a sparse matrix',type=int,default=1000)
    parser.add_argument(
        "-t",
        "--num-trees",
        help="number of trees in a random forest",
        type=int,
        default=10,
    )

    # Feature information --fixed
    parser.add_argument(
        "-n",
        "--num-features",
        help="number of features in input vector",
        type=int,
        default=6819,
    )
    # parser.add_argument('-i','--rnaseq-min',help='minimum expression level in RNA-seq data',type=float,default=8e-05)
    # parser.add_argument('-x','--rnaseq-max',help='maximum expression level in RNA-seq data',type=float,default=1.11729e06)

    # Hyperparameters for random forest
    parser.add_argument("-md", "--max-depth", type=str, default="None")
    parser.add_argument("-mss", "--min-samples-split", type=float, default=2)
    parser.add_argument("-msl", "--min-samples-leaf", type=float, default=32)
    parser.add_argument("-maf", "--max_features", type=str, default="auto")
    parser.add_argument("-bs", "--bootstrap", action="store_true")

    args = parser.parse_args()

    fout = open(args.output_filename_prefix + ".progress.txt", "w")

    # print ('\nReading input parameters...')
    for key, value in vars(args).items():
        if "max_" not in key and "min_" not in key:
            fout.write("# %s: %s\n" % (key, value))

    # Define range for raw RNA-seq values
    # rnaseq_range = [args.rnaseq_min,args.rnaseq_max]

    # Set random seed
    random.seed(args.seed)
    np.random.seed(args.seed)

    # Read training data
    training_inputs, training_labels = readData(
        args.training_data_filenames,
        args.pairwise_dist,
        args.window_size,
        args.positive_training_data_size,
        args.num_features,
        args.frac_feature,
    )
    # print(training_inputs[0, :])
    # print(np.sum(training_inputs, axis=0), np.sum(training_inputs, axis=1))

    # Read tuning data
    tuning_inputs, tuning_labels = readData(
        args.tuning_data_filenames,
        args.pairwise_dist,
        args.window_size,
        args.positive_tuning_data_size,
        args.num_features,
        args.frac_feature,
    )

    test_data_available = args.test_data_filenames
    # Read test data if provided
    if test_data_available:
        test_inputs, test_labels = readData(
            args.test_data_filenames,
            args.pairwise_dist,
            args.window_size,
            args.positive_tuning_data_size,
            args.num_features,
            args.frac_feature,
        )

    fout.write(
        "seed\tmd\tmss\tmsl\tmaf\tbs\ttrmse\ttrauroc\ttrauprc\ttrmean\ttrpmean\ttrnmean\tvamse\tvaauroc\tvaauprc\tvamean\tvapmean\tvanmean"
    )
    if test_data_available:
        fout.write("\ttemse\tteauroc\tteauprc\ttemean\ttepmean\ttenmean")
    fout.write("\n")

    # Train and report performance metrics
    for i in range(args.random_search_n):
        if args.random_search:
            (
                max_depth,
                min_samples_split,
                min_samples_leaf,
                max_features,
                bootstrap,
            ) = setTraininghyperparams(int(i * args.seed), args.num_features)
        else:
            max_depth, min_samples_split, min_samples_leaf, max_features, bootstrap = (
                args.max_depth,
                args.min_samples_split,
                args.min_samples_leaf,
                args.max_features,
                args.bootstrap,
            )
            max_depth = None if max_depth == "None" else int(max_depth)

        seed = i if args.random_search else args.seed
        hyperparams = [
            seed,
            max_depth,
            min_samples_split,
            min_samples_leaf,
            max_features,
            bootstrap,
        ]

        ### Training starts here ###
        clf = RandomForestClassifier(
            random_state=args.seed,
            n_estimators=args.num_trees,
            max_depth=max_depth,
            min_samples_split=min_samples_split,
            min_samples_leaf=min_samples_leaf,
            max_features=max_features,
            bootstrap=bootstrap,
            n_jobs=-1,
        )
        clf.fit(training_inputs, training_labels)
        ### Training ends here ###

        # Evaluate
        training_result = eval(clf, training_inputs, training_labels)
        tuning_result = eval(clf, tuning_inputs, tuning_labels)
        test_result = eval(clf, test_inputs, test_labels) if test_data_available else ()

        # Print out result
        fout.write(
            "\t".join(
                [str(s) for s in hyperparams] + ["%.6f" % s for s in training_result + tuning_result + test_result]
            )
            + "\n"
        )

        # If this is not for hyperparameter search, train one model and quit
        if not args.random_search:
            break

    fout.close()

    if args.random_search:
        df = pd.read_table(args.output_filename_prefix + ".progress.txt", comment="#")
        df = df.sort_values(by="vaauroc", ascending=False)

        with open(args.output_filename_prefix + ".best_hyperparam.txt", "w") as f:
            for p in ["md", "mss", "msl", "maf"]:
                f.write("-%s\n%s\n" % (p, str(df[p].iloc[0])))
            if df["bs"].iloc[0]:
                f.write("-bs\n")

    # Save model and generate predictions for input data if specified
    elif args.save:

        # Pickle model
        joblib.dump(clf, args.output_filename_prefix + ".pkl")

        # Make prediction on training data
        training_y_pred = clf.predict_proba(training_inputs)[:, 1]
        with gzip.open(args.output_filename_prefix + ".training_pred.txt.gz", "wb") as fout:
            for i in range(len(training_labels)):
                l = "%d\t%.6f\n" % (training_labels[i], training_y_pred[i])
                fout.write(l.encode())

        # Make prediction on tuning data
        tuning_y_pred = clf.predict_proba(tuning_inputs)[:, 1]
        with gzip.open(args.output_filename_prefix + ".tuning_pred.txt.gz", "wb") as fout:
            for i in range(len(tuning_labels)):
                l = "%d\t%.6f\n" % (tuning_labels[i], tuning_y_pred[i])
                fout.write(l.encode())

        if test_data_available:
            # Make prediction on test data
            test_y_pred = clf.predict_proba(test_inputs)[:, 1]
            with gzip.open(args.output_filename_prefix + ".test_pred.txt.gz", "wb") as fout:
                for i in range(len(test_labels)):
                    l = "%d\t%.6f\n" % (test_labels[i], test_y_pred[i])
                    fout.write(l.encode())


main()
