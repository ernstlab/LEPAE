import scipy.sparse, numpy as np, joblib, random, linecache, gzip, torch
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc, average_precision_score, mean_squared_error, jaccard_score

# Given a set of files containing input feature data, determine which lines to read to define positive and negative pairs
def determineLinesToRead(
    args, filenames, data_size, seed=413666, subsample=True, incl_negative=True, incl_flipped=True, shuffle=True
):

    # 0. Set seed to ensure the same set of lines are read given the same set of input files
    random.seed(seed)
    np.random.seed(seed)

    # 1. Define number of lines to read/skip
    num_lines = np.array([len(open(fn).readlines()) for fn in filenames])  # number of lines in each file

    if subsample:
        num_lines_to_sample = np.floor(data_size * num_lines / np.sum(num_lines))  # number of lines to sample from each file
        while sum(num_lines_to_sample) < data_size:  # adjust the number of lines to sample if needed
            num_lines_to_sample[random.randrange(len(filenames))] += 1
    else:
        num_lines_to_sample = num_lines

    num_lines_btw_windows = int(
        args.pairwise_dist / args.window_size
    )  # number of lines to skip to define pairs with the specified distance

    # 2. Determine which lines to read
    _upstream_lines_to_read = []
    _downstream_lines_to_read = []
    _labels = []
    for i in range(len(filenames)):
        fn = filenames[i]

        # Add positive pairs (as read from the input file)
        if subsample:
            line_indices = list(
                random.sample(range(num_lines[i] - num_lines_btw_windows), int(num_lines_to_sample[i]))
                if (num_lines[i] - num_lines_btw_windows) > int(num_lines_to_sample[i])
                else range(num_lines[i] - num_lines_btw_windows)
            )
        else:
            line_indices = list(range(num_lines[i] - num_lines_btw_windows))
        _upstream_lines_to_read += [(fn, 1 + j) for j in line_indices]
        _downstream_lines_to_read += [(fn, 1 + j + num_lines_btw_windows) for j in line_indices]
        _labels += [1] * len(line_indices)

        # Add flipped pairs if specified
        if incl_flipped:
            _upstream_lines_to_read += [(fn, 1 + j + num_lines_btw_windows) for j in line_indices]
            _downstream_lines_to_read += [(fn, 1 + j) for j in line_indices]
            _labels += [1] * len(line_indices)

        if incl_negative:
            # Define shuffled line indices to generate negative pairs
            shuffled_line_indices = np.random.permutation(line_indices)
            _upstream_lines_to_read += [(fn, 1 + j) for j in line_indices]
            _downstream_lines_to_read += [(fn, 1 + j + num_lines_btw_windows) for j in shuffled_line_indices]
            _labels += [0] * len(line_indices)

            if incl_flipped:
                _upstream_lines_to_read += [(fn, 1 + j + num_lines_btw_windows) for j in shuffled_line_indices]
                _downstream_lines_to_read += [(fn, 1 + j) for j in line_indices]
                _labels += [0] * len(line_indices)

    # Shuffle the ordering so that each batch has a mix of positive and negative pairs if using batches
    if shuffle:
        s = np.random.permutation(
            len(_upstream_lines_to_read),
        )
        upstream_lines_to_read = [_upstream_lines_to_read[i] for i in s]
        downstream_lines_to_read = [_downstream_lines_to_read[i] for i in s]
        labels = [_labels[i] for i in s]
    else:
        upstream_lines_to_read = _upstream_lines_to_read
        downstream_lines_to_read = _downstream_lines_to_read
        labels = _labels

    return upstream_lines_to_read, downstream_lines_to_read, labels


# Read one batch to use as input to a PyTorch neural network
def readBatchTorch(
    current_batch_size, num_features, frac_feature, upstream_lines_to_read, downstream_lines_to_read, labels, line_index
):
    # Convert a string into a feature index in integer for a downstream window
    def _convert_downstream_feature(s):
        return num_features + int(s)

    x, y = np.zeros((current_batch_size, num_features * 2)), np.zeros((current_batch_size))
    for j in range(current_batch_size):  # each pair in this batch

        # Read a line containing features characterizing the upstream window of the current pair
        upstream_line = linecache.getline(*upstream_lines_to_read[line_index + j]).strip().split("|")
        if len(upstream_line) > 1:

            # Determine which features should be set to a non-zero value
            upstream_nonzero_feature_indices = list(map(int, upstream_line[1].split()))

            # If using binary features, set those features to 1
            # If using fractions as features, set those features to values specified after the second vertical bar
            x[j, upstream_nonzero_feature_indices] = list(map(float, upstream_line[2].split())) if frac_feature else 1

        # Do the same for the downstream window of the current pair
        downstream_line = linecache.getline(*downstream_lines_to_read[line_index + j]).strip().split("|")
        if len(downstream_line) > 1:
            downstream_nonzero_feature_indices = list(map(_convert_downstream_feature, downstream_line[1].split()))
            x[j, downstream_nonzero_feature_indices] = list(map(float, downstream_line[2].split())) if frac_feature else 1

        # Set the label (1 or 0) for the current pair
        y[j] = labels[line_index + j]

    # Convert numpy arrays into tensors accepted by PyTorch
    x, y = torch.autograd.Variable(torch.from_numpy(x).float()), torch.autograd.Variable(torch.from_numpy(y).float())
    return x, y


# Read the entire data and store it in a sparse CSR matrix to use as input to a random forest
def readDataCSR(args, upstream_lines_to_read, downstream_lines_to_read, labels):
    row, col, data = [], [], []  # row indices, column indices, feature values for later defining a sparse matrix
    n = 0
    for i in range(len(upstream_lines_to_read)):
        upstream_line = linecache.getline(*upstream_lines_to_read[i]).strip().split("|")
        upstream_nonzero_feature_indices = [int(s) for s in upstream_line[1].split()] if len(upstream_line) > 1 else []
        downstream_line = linecache.getline(*downstream_lines_to_read[i]).strip().split("|")
        downstream_nonzero_feature_indices = (
            [args.num_features + int(s) for s in downstream_line[1].split()] if len(downstream_line) > 1 else []
        )
        num_nonzero_feature_indices = len(upstream_nonzero_feature_indices) + len(downstream_nonzero_feature_indices)
        row += [i] * num_nonzero_feature_indices
        col += upstream_nonzero_feature_indices + downstream_nonzero_feature_indices

        if args.frac_feature:
            data += [float(s) for s in upstream_line[2].split()] + [float(s) for s in downstream_line[2].split()]
        else:
            data += [1] * num_nonzero_feature_indices
        n += 1

    inputs = scipy.sparse.csr_matrix(
        (data, (row, col)), shape=(n, args.num_features * 2), dtype=float if args.frac_feature else bool
    )
    # label = [True, True, False, False] * data_size

    return inputs, labels


def predictDT(clf, inputs):  # , hard_voting=False, num_workers=4):
    # if not hard_voting:
    #     return clf.predict_proba(inputs)[:, 1]

    # def _predictTree(tree, inputs):
    #     return tree.predict(inputs)

    # votes = joblib.Parallel(n_jobs=num_workers)(
    #     joblib.delayed(_predictTree)(clf.estimators_[i], inputs) for i in range(clf.n_estimators)
    # )
    # scores = np.mean(np.array(votes).T, axis=1)

    return clf.predict_proba(inputs)[:, 1]


def predictNNFromFiles(args, net, upstream_lines_to_read, downstream_lines_to_read, labels, return_inputs=False):
    bs = args.fixed_batch_size
    num_batch = int(np.floor(len(upstream_lines_to_read) / bs) + 1)

    line_index = 0
    scores = np.zeros((len(upstream_lines_to_read)))
    inputs = []
    if return_inputs:
        inputs = np.zeros((len(upstream_lines_to_read), args.num_features * 2))
    for i in range(num_batch):
        current_batch_size = bs if ((i + 1) * bs) < len(upstream_lines_to_read) else len(upstream_lines_to_read) % bs
        # print(i, data_size, num_batch, current_batch_size)

        x, _ = readBatchTorch(
            current_batch_size,
            args.num_features,
            args.frac_feature,
            upstream_lines_to_read,
            downstream_lines_to_read,
            labels,
            line_index,
        )
        line_index += current_batch_size
        scores[i * bs : i * bs + current_batch_size] = net(x).data

        if return_inputs:
            inputs[i * bs : i * bs + current_batch_size, :] = x.numpy()

    return scores, inputs


def computeJaccardIndex(args, inputs):
    inputs = inputs.toarray()

    upstream_inputs = inputs[:, : args.num_features]
    downstream_inputs = inputs[:, args.num_features :]

    score = np.ones(len(inputs)) * -1

    for i in range(len(inputs)):
        x, y = upstream_inputs[i, :], downstream_inputs[i, :]
        if x.sum() > 0 or y.sum() > 0:
            score[i] = jaccard_score(x, y, zero_division=-1)

    return score


def eval(label, scores):
    loss = torch.nn.functional.binary_cross_entropy(
        torch.from_numpy(np.array(scores)).float(), torch.from_numpy(np.array(label)).float()
    ).item()
    mse = mean_squared_error(label, scores)
    fpr, tpr, _ = roc_curve(label, scores)
    auroc = auc(fpr, tpr)  # au-ROC
    auprc = average_precision_score(label, scores)  # au-PRC
    pos_mean_score = np.mean([scores[i] for i in range(len(scores)) if label[i] == 1])
    neg_mean_score = np.mean([scores[i] for i in range(len(scores)) if label[i] == 0])
    return [loss, mse, auroc, auprc, pos_mean_score, neg_mean_score]


def writeScores(output_filename, frac_feature, labels, scores, inputs=[]):
    with gzip.open(output_filename, "wb") as fout:
        for i in range(len(labels)):
            l = "%d\t%.10f" % (labels[i], scores[i])
            if len(inputs) > 0:
                l += "\t"
                l += "\t".join(["%.10f" % float(s) if frac_feature else str(int(s)) for s in list(inputs[i, :])])  # + "\t"
            l += "\n"
            fout.write(l.encode())


class SiameseNet(torch.nn.Module):
    def __init__(self, num_feature, hyperparam):
        super(SiameseNet, self).__init__()
        self.num_feature = num_feature

        self.init_subnetwork = torch.nn.Sequential()

        self.init_subnetwork.add_module("il", torch.nn.Linear(num_feature, int(hyperparam["num_neuron_init"]), bias=False))
        self.init_subnetwork.add_module("id", torch.nn.Dropout(p=hyperparam["dropout_rate"]))
        self.init_subnetwork.add_module("ir", torch.nn.ReLU())

        self.fin_subnetwork = torch.nn.Sequential()
        self.fin_subnetwork.add_module(
            "fl", torch.nn.Linear(int(hyperparam["num_neuron_init"]) * 2, int(hyperparam["num_neuron_fin"]))
        )
        self.fin_subnetwork.add_module("fd", torch.nn.Dropout(p=hyperparam["dropout_rate"]))
        self.fin_subnetwork.add_module("fr", torch.nn.ReLU())
        self.fin_subnetwork.add_module("end", torch.nn.Linear(int(hyperparam["num_neuron_fin"]), 1))
        self.fin_subnetwork.add_module("sigmoid", torch.nn.Sigmoid())

    def forward(self, x):
        a = self.init_subnetwork.forward(x[:, : self.num_feature])
        b = self.init_subnetwork.forward(x[:, self.num_feature :])
        c = torch.cat((a, b), 1)
        y = self.fin_subnetwork.forward(c)
        y = y.view(c.size()[0])
        return y
