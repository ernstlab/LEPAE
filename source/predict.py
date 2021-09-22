import sys, gzip, random, numpy as np, argparse, joblib, linecache


def predict(
    clf,
    input_filename,
    pairwise_dist,
    window_size,
    batch_size,
    num_features,
    frac_feature,
    output_filename,
):
    num_lines = len(open(input_filename).readlines())
    num_lines_between_windows = int(pairwise_dist / window_size)
    print(
        "Number of lines between %d-kb windows that are %d kb apart: %d"
        % (
            int(window_size / 1000),
            int(pairwise_dist / 1000),
            num_lines_between_windows,
        )
    )

    dt = float if frac_feature else bool

    with gzip.open(
        output_filename if output_filename.endswith(".gz") else output_filename + ".gz",
        "wb",
    ) as fout:
        upstream_window_indices = range(num_lines - num_lines_between_windows)
        current_window_index = 0

        for i in range(int(len(upstream_window_indices) / batch_size) + 1):  # iterate through each batch
            current_batch_size = (
                batch_size
                if i < int(len(upstream_window_indices) / batch_size)
                else len(upstream_window_indices) % batch_size
            )  # last batch may have items fewer than N
            if current_batch_size == 0:
                break

            coords = []
            X = np.zeros((current_batch_size, num_features * 2), dtype=dt)
            Y = np.zeros((current_batch_size, num_features * 2), dtype=dt)  # flipped pairs
            for j in range(current_batch_size):  # iterate through each sample within the batch
                line = (
                    linecache.getline(
                        input_filename,
                        1 + upstream_window_indices[current_window_index],
                    )
                    .strip()
                    .split("|")
                )  # read a random line and split by space
                coord1 = line[0].strip()
                nonzero_feature_indices = (
                    np.array([int(s) for s in line[1].split()]) if len(line) > 1 else []
                )  # read indices of nonzero features
                v = [float(s) for s in line[2].split()] if (len(line) > 2 and frac_feature) else True
                if len(nonzero_feature_indices) > 0:
                    X[j, nonzero_feature_indices] = v
                    Y[j, num_features + nonzero_feature_indices] = v

                line = (
                    linecache.getline(
                        input_filename,
                        1 + upstream_window_indices[current_window_index] + num_lines_between_windows,
                    )
                    .strip()
                    .split("|")
                )
                coord2 = line[0].strip()
                nonzero_feature_indices = np.array([int(s) for s in line[1].split()]) if len(line) > 1 else []
                v = [float(s) for s in line[2].split()] if (len(line) > 2 and frac_feature) else True
                if len(nonzero_feature_indices) > 0:
                    X[j, num_features + nonzero_feature_indices] = v
                    Y[j, nonzero_feature_indices] = v

                current_window_index += 1
                coords.append(coord1 + "\t" + coord2)

            # Make prediction on current batch
            X_pred = clf.predict_proba(X)[:, 1]
            Y_pred = clf.predict_proba(Y)[:, 1]
            lines = ""
            for j in range(len(X_pred)):
                lines += "%s\t%.6f\t%.6f\t%.6f\n" % (
                    coords[j],
                    (X_pred[j] + Y_pred[j]) / 2.0,
                    X_pred[j],
                    Y_pred[j],
                )
            fout.write(lines.encode())


def main():

    parser = argparse.ArgumentParser(prog="base_maker", description="Predict score given a trained neural network")
    parser.add_argument(
        "-t",
        "--trained-classifier-filename",
        help="path to trained classifier (.pkl or .pkl.gz)",
        type=str,
    )
    parser.add_argument("-i", "--input-filename", help="path to input file", type=str)

    parser.add_argument(
        "-d",
        "--pairwise-dist",
        help="pairwise distance between windows",
        type=int,
        default=1000,
    )
    parser.add_argument("-w", "--window-size", help="genomic window size in bp", type=int, default=1000)

    parser.add_argument("-b", "--batch-size", help="batch size", type=int, default=1000)

    parser.add_argument("-s", "--seed", help="random seed", type=int, default=1)

    parser.add_argument(
        "-n",
        "--num-features",
        help="number of features in input vector",
        type=int,
        default=6819,
    )
    parser.add_argument(
        "-f",
        "--frac-feature",
        help="use fraction of overlap as features instead of binary",
        action="store_true",
    )
    parser.add_argument("-o", "--output-filename", help="path to output file", type=str)

    args = parser.parse_args()

    # Set random seed
    random.seed(args.seed)
    np.random.seed(args.seed)

    # Load previously trained classifier
    clf = joblib.load(args.trained_classifier_filename)

    # Make predictions
    predict(
        clf,
        args.input_filename,
        args.pairwise_dist,
        args.window_size,
        args.batch_size,
        args.num_features,
        args.frac_feature,
        args.output_filename,
    )


if __name__ == "__main__":
    main()
