import gzip, random, numpy as np, argparse, joblib, linecache
from sklearn.metrics import jaccard_score
from shared import *
import tempfile

NUM_WORKERS = 4


def main():
    parent_parser = argparse.ArgumentParser(
        prog="base_maker", description="Predict score given a trained classifier", add_help=False
    )
    parent_parser.add_argument(
        "-t",
        "--trained-classifier-filenames",
        help="paths to trained classifiers (.pkl or .pkl.gz)",
        type=str,
        nargs="+",
    )
    parent_parser.add_argument(
        "-i",
        "--input-filenames",
        help="paths to input files",
        type=str,
        nargs="+",
    )
    parent_parser.add_argument("-j", "--data-size", help="number of pairs to read", type=int, default=10000)
    parent_parser.add_argument("-d", "--pairwise-dist", help="pairwise distance between windows", type=int, default=1000)
    parent_parser.add_argument("-u", "--subsample", help="subsample pairs", action="store_true")
    parent_parser.add_argument(
        "-a", "--incl-negative", help="generate and make predictions for negative pairs", action="store_true"
    )
    parent_parser.add_argument(
        "-c", "--print-feature", help="print inputs features along with predictions", action="store_true"
    )
    # parent_parser.add_argument(
    #     "-p", "--print-decision-path", help="print decision tree path if using decision tree", action="store_true"
    # )
    # parent_parser.add_argument(
    #     "-l", "--incl-flipped", help="generate and make predictions for flipped pairs", action="store_true"
    # )

    parent_parser.add_argument("-w", "--window-size", help="genomic window size in bp", type=int, default=1000)
    parent_parser.add_argument("-b", "--fixed-batch-size", help="batch size", type=int, default=512)
    parent_parser.add_argument("-s", "--seed", help="random seed", type=int, default=1)
    parent_parser.add_argument("-n", "--num-features", help="number of features in input vector", type=int, default=6512)
    parent_parser.add_argument(
        "-f", "--frac-feature", help="use fraction of overlap as features instead of binary", action="store_true"
    )

    parent_parser.add_argument("-o", "--output-filename", help="path to output file", type=str)

    main_parser = argparse.ArgumentParser()
    subparsers = main_parser.add_subparsers(help="type of classifier to train", dest="classifier", required=True)
    subparsers.add_parser("DT", help="Make predictions using a trained decision tree classifier", parents=[parent_parser])
    subparsers.add_parser("NN", help="Make predictions using a trained neural network classifier", parents=[parent_parser])
    subparsers.add_parser("Jaccard", help="Report Jaccard index (no training required)", parents=[parent_parser])

    args = main_parser.parse_args()

    # Set random seed
    random.seed(args.seed)
    np.random.seed(args.seed)

    # Determine ensemble size from the list of trained classifiers
    ensemble_size = len(args.trained_classifier_filenames) if args.classifier in ["DT", "NN"] else 1

    # Generate a temporary file with input data concatenated
    tmp_fout = tempfile.NamedTemporaryFile(mode="w", delete=True, prefix=args.output_filename, suffix=".tmp", dir='./')
    temp_filename = tmp_fout.name
    print(temp_filename)
    for input_filename in args.input_filenames:
        with open(input_filename, "r") as fin:
            tmp_fout.writelines(fin.readlines())
    tmp_fout.seek(0)

    # Determine which lines to read
    print("Determining which inputs lines to read...")
    upstream_lines_to_read, downstream_lines_to_read, labels = determineLinesToRead(
        args,
        [temp_filename],
        args.data_size,
        subsample=args.subsample,
        incl_negative=args.incl_negative,
        incl_flipped=True,
        shuffle=False,
    )

    # Load previously trained classifier
    print("Generating predictions...")
    if args.classifier in ["DT", "NN"]:
        print("\n".join(args.trained_classifier_filenames))
        classifiers = [joblib.load(c) for c in args.trained_classifier_filenames]

    scores = np.ones((len(upstream_lines_to_read), ensemble_size)) * -1

    if args.classifier in ["DT", "Jaccard"]:
        bs = args.fixed_batch_size
        num_batch = int(np.floor(len(upstream_lines_to_read) / bs) + 1)
        # paths = []
        for i in range(num_batch):
            current_batch_size = bs if ((i + 1) * bs) < len(upstream_lines_to_read) else len(upstream_lines_to_read) % bs
            if current_batch_size == 0:
                continue
            start, end = i * bs, i * bs + current_batch_size

            x, _ = readDataCSR(
                args, upstream_lines_to_read[start:end], downstream_lines_to_read[start:end], labels[start:end]
            )
            if args.classifier == "DT":
                for j in range(ensemble_size):
                    scores[start:end, j] = predictDT(classifiers[j], x)
                # if args.print_decision_path:
                #     indicator, n_nodes_ptr = clf.decision_path(x)
                #     for j in range(current_batch_size):
                #         _, c = indicator[j, n_nodes_ptr[0] : n_nodes_ptr[1]].nonzero()
                #         paths.append(list(c))
            else:
                # print(i, current_batch_size, start, end, x.shape)
                y = computeJaccardIndex(args, x)
                scores[start:end, 0] = y

    elif args.classifier == "NN":
        for i in range(ensemble_size):
            clf = classifiers[i]
            clf.eval()
            scores[:, i], _ = predictNNFromFiles(args, clf, upstream_lines_to_read, downstream_lines_to_read, labels)

    print("Printing predictions to output file...")
    with gzip.open(args.output_filename, "wt") as fout:
        if not args.subsample and not args.incl_negative:
            n = int(len(upstream_lines_to_read) / 2)
            for i in range(n):
                upstream_coord = linecache.getline(*upstream_lines_to_read[i]).strip().split("|")[0].split()[:3]
                downstream_coord = linecache.getline(*downstream_lines_to_read[i]).strip().split("|")[0].split()[:3]
                s1 = np.mean(scores[i, :])
                s2 = np.mean(scores[n + i, :])
                s = (s1 + s2) / 2
                l = "\t".join(
                    upstream_coord
                    + downstream_coord
                    + ["%s-%s-%s" % (upstream_coord[0], upstream_coord[1], downstream_coord[1]), "%d" % (s * 1000), ".", "."]
                    + ["%.5f" % v for v in [s, s1, s2]]
                    + ["%d" % args.window_size, "%d" % args.pairwise_dist, "1"]
                    + [",".join(["%.5f" % v for v in scores[i, :]]), ",".join(["%.5f" % v for v in scores[n + i, :]])]
                )
                fout.write(l + "\n")

        elif args.incl_negative:
            n = int(len(upstream_lines_to_read) / 4)
            for i in range(n):
                upstream_line = linecache.getline(*upstream_lines_to_read[i]).strip().split("|")
                downstream_line = linecache.getline(*downstream_lines_to_read[i]).strip().split("|")
                upstream_coord = upstream_line[0].split()[:3]
                downstream_coord = downstream_line[0].split()[:3]
                s1 = np.mean(scores[i, :])
                s2 = np.mean(scores[n + i, :])
                assert labels[i] == 1
                assert labels[n + i] == 1
                s = (s1 + s2) / 2
                l = "\t".join(
                    upstream_coord
                    + downstream_coord
                    + ["%s-%s-%s" % (upstream_coord[0], upstream_coord[1], downstream_coord[1]), "%d" % (s * 1000), ".", "."]
                    + ["%.5f" % v for v in [s, s1, s2]]
                    + ["%d" % args.window_size, "%d" % args.pairwise_dist, "1"]
                    + [",".join(["%.5f" % v for v in scores[i, :]]), ",".join(["%.5f" % v for v in scores[n + i, :]])]
                )

                # if args.print_decision_path:
                #     l += "\t|\t" + "\t".join([str(s) for s in paths[i]])
                if args.print_feature:
                    l += "\t|\t%s\t|\t%s" % (upstream_line[1], downstream_line[1])
                fout.write(l + "\n")

                upstream_line = linecache.getline(*upstream_lines_to_read[2 * n + i]).strip().split("|")
                downstream_line = linecache.getline(*downstream_lines_to_read[2 * n + i]).strip().split("|")
                upstream_coord = upstream_line[0].split()[:3]
                downstream_coord = downstream_line[0].split()[:3]
                s1 = np.mean(scores[2 * n + i, :])
                s2 = np.mean(scores[3 * n + i, :])
                assert labels[2 * n + 1] == 0
                assert labels[3 * n + i] == 0
                s = (s1 + s2) / 2
                l = "\t".join(
                    upstream_coord
                    + downstream_coord
                    + ["%s-%s-%s" % (upstream_coord[0], upstream_coord[1], downstream_coord[1]), "%d" % (s * 1000), ".", "."]
                    + ["%.5f" % v for v in [s, s1, s2]]
                    + ["%d" % args.window_size, "%d" % args.pairwise_dist, "0"]
                    + [
                        ",".join(["%.5f" % v for v in scores[2 * n + i, :]]),
                        ",".join(["%.5f" % v for v in scores[3 * n + i, :]]),
                    ]
                )
                # if args.print_decision_path:
                #     l += "\t|\t" + "\t".join([str(s) for s in paths[2 * n + i]])
                if args.print_feature:
                    l += "\t|\t%s\t|\t%s" % (upstream_line[1], downstream_line[1])
                fout.write(l + "\n")

    tmp_fout.close()


main()
