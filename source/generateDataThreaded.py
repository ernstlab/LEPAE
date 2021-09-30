import sys, gzip, threading, pandas as pd, argparse, os
from collections import defaultdict

# SOURCE: https://www.tutorialspoint.com/python3/python_multithreading.htm
class myThread(threading.Thread):
    # featureStrings = {key: str() for key in range(1, 3)}  # 5)}
    printInit = 0

    def __init__(
        self,
        thread_id,
        name,
        feature_cat,
        directory,
        input_positions,
        active_indices,
        num_chromatin_states,
        feature_indices,
        frac_feature,
    ):
        threading.Thread.__init__(self)
        self.thread_id = thread_id
        self.name = name
        self.feature_cat = feature_cat
        self.directory = directory
        self.input_positions = input_positions
        self.active_indices = active_indices
        self.num_chromatin_states = num_chromatin_states
        self.feature_index = feature_indices[feature_cat]
        self.frac_feature = frac_feature
        # self.real_values = None

    def run(self):
        if self.feature_cat == 1:
            readDnaseChipFeature(
                self.directory,
                self.input_positions,
                self.active_indices,
                self.frac_feature,
            )
        elif self.feature_cat == 2:
            readChromHmmFeature(
                self.directory,
                self.input_positions,
                self.active_indices,
                self.num_chromatin_states,
                self.feature_index,
                self.frac_feature,
            )
        # elif self.feature == 3:
        #     readCageFeature(self.directory, self.input_positions, self.active_indices,
        #                     self.num_chromatin_states, self.feature_index)
        # elif self.feature == 4:
        #     self.real_values = readRnaSeqFeature(self.directory, self.input_positions,
        #                                        self.active_indices, self.feature_index)

    # def displayFeatureProgress(feature, stringUpdate):
    #     if myThread.printInit:
    #         sys.stdout.write(4 * "\033[F")
    #     else:
    #         myThread.printInit = 1
    #     myThread.featureStrings[feature] = stringUpdate
    #     displayString = str()
    #     for i in range(1, 3):  # 5):
    #         displayString += myThread.featureStrings[i] + "\n"
    #     sys.stdout.write(displayString)
    #     sys.stdout.flush()


def readDnaseChipFeature(dnasechipseq_list_file, input_positions, active_indices, frac_feature):
    # List of files containing position indices with overlapping peak in different DNase-seq and ChIP-seq data
    dnasechipseq_files = open(dnasechipseq_list_file).readlines()

    # Iterate through each file (each experiment in a specific cell-type and with a specific target if it's ChIP-seq)
    for i in range(len(dnasechipseq_files)):
        dnasechipseq_file = dnasechipseq_files[i].strip()
        if os.stat(dnasechipseq_file).st_size == 0:  # skip if empty file
            continue
        try:
            t = pd.read_table(dnasechipseq_file, engine="c", header=None).values
            pos = t[:, 0]
            cnt = t[:, 1]

            valid = list(input_positions.intersection(pos))
            for j in range(len(valid)):
                if frac_feature:
                    active_indices[valid[j]].append((i, cnt[j]))
                else:
                    active_indices[valid[j]].append(i)
        except pd.errors.EmptyDataError as _:  # skip if empty file
            continue

        # Status output
        # p = int((i + 1) / len(dnasechipseq_files) * 100)
        # displayString = (
        #     "\tDNase-seq and ChIP-seq ["
        #     + "=" * p
        #     + " " * (100 - p)
        #     + "] "
        #     + str(p)
        #     + "%\t"
        # )
        # myThread.displayFeatureProgress(1, displayString)


def readChromHmmFeature(
    chromhmm_list_file,
    input_positions,
    active_indices,
    chromhmm_num_states,
    feature_index,
    frac_feature,
):
    num_current_features = feature_index
    #    active_indices = defaultdict(list) # Key: position index, value: non-zero feature index for the given position index

    # List of files containing position indices and their overlapping ChromHMM state for each cell-type
    chromhmm_files = open(chromhmm_list_file).readlines()  # sorted(os.listdir(chromhmm_list_file))
    for i in range(len(chromhmm_files)):  # Iterate through each file (each cell-type)
        chromhmm_file = chromhmm_files[i].strip()
        if os.stat(chromhmm_file).st_size == 0:  # skip if file is empty (unlikely to happen)
            continue
        with gzip.open(chromhmm_file, "rb") as f:
            positions_found = 0
            for line in f:
                l = line.strip().split()
                if frac_feature:
                    position = int(l[0].decode("utf-8"))
                    if position in input_positions:
                        state = num_current_features + int(l[1].decode("utf-8")) - 1
                        cnt = int(l[2].decode("utf-8"))
                        updated = False
                        for (k, v) in active_indices[position]:
                            # print(k, v, state, cnt, k == state)
                            if k == state:  # state annoated the window earlier
                                active_indices[position].remove((k, v))
                                active_indices[position].append((k, v + cnt))
                                updated = True
                                # print("!!!", k, v, v + cnt, active_indices[position])
                                break
                        if not updated:
                            active_indices[position].append((state, cnt))
                else:
                    position = int(l[0].decode("utf-8"))
                    if position in input_positions:
                        positions_found += 1
                        states = l[1].decode("utf-8").split(",")
                        # state = int(state[1:]) if state.startswith('U') else int(state)
                        active_indices[position] += [num_current_features + (int(s) - 1) for s in states]

                    # Optimization added 10 July 2019
                    if positions_found == len(input_positions):
                        break
        num_current_features += chromhmm_num_states

        # Status output
        # p = int((i + 1) / len(chromhmm_files) * 100)
        # displayString = "\tChromHMM [" + "=" * p + " " * (100 - p) + "] " + str(p) + "%"
        # myThread.displayFeatureProgress(2, displayString)


# def readCageFeature(cage_dir, input_positions, active_indices, chromhmm_num_states, feature_index):
# #    num_current_features = 0 # Number of features processed so far
#     num_current_features = feature_index
# #    active_indices = defaultdict(list) # Key: position index, value: non-zero feature index for the given position index
#
#     # A file containing position indices and their CAGE peak data across multiple cell-types
#     cage_file = os.listdir(cage_dir)[0]
#     try:
#         df = pd.read_table(cage_dir+cage_file,engine='c',header=None).as_matrix()
#         positions = df[:,0] # Position indices
#         features = df[:,1:] # Presence of CAGE peak in each position in each cell-type
#         for i in range(len(positions)):
#             if positions[i] in input_positions:
#                 a = num_current_features + np.where(features[i,:]>0)[0] # active CAGE features
#                 for j in a:
#                     active_indices[positions[i]].append(j)
#             # Status output
#             p = int((i+1)/len(positions)*100)
#             displayString = "\tCAGE [" + "=" * p + " "*(100-p) + "] " + str(p)+"% {0}".format(len(features[0]))
#             myThread.displayFeatureProgress(3, displayString)
#
#         num_current_features += len(features[0])
#         #print("DEBUG: CAGE index after update {0}\033[F\033[F".format(num_current_features))
#     except (pd.errors.EmptyDataError,pd.io.common.EmptyDataError) as _:
#         #print('! Empty CAGE data file', cage_file)
#         num_current_features += 1829 if chromhmm_num_states==25 else 1073  ### TODO: THIS WAS HARDCODED. IF USING SPECIES OTHER THAN HUMAN AND MOUSE THIS SHOULD BE FIXED
# #        num_current_features += 1829 if chromhmm_num_states==25 else 13  ### TODO: THIS WAS HARDCODED. IF USING SPECIES OTHER THAN HUMAN AND MOUSE THIS SHOULD BE FIXED
#
# def readRnaSeqFeature(rnaseq_dir, input_positions, active_indices, feature_index):
#     num_current_features = feature_index
#     real_values = defaultdict(list) # Key: position index, value: real value for the non-binary features
#
#     # List of files containing position indices and their RNA-seq level in different cell-types
#     rnaseq_files = sorted(os.listdir(rnaseq_dir))
#     for i in range(len(rnaseq_files)): # Iterate through each file (each cell-type)
#         rnaseq_file = rnaseq_files[i]
#         last_pos = int()
#         try:
#             df = pd.read_table(rnaseq_dir+rnaseq_file,engine='c',header=None).as_matrix()
#             positions = df[:,0] # Position indices
#             signals = df[:,1] # RNA-seq level of each position in the current cell-type
#             for j in range(len(positions)):
#                 if positions[j] in input_positions:
#                     active_indices[positions[j]].append(num_current_features+i)
#                     last_pos = num_current_features+i
#                     real_values[positions[j]].append(signals[j])
#             # Status output
#             p = int((i+1)/len(rnaseq_files)*100)
#             displayString = "\tRNA-seq [" + "=" * p + " "*(100-p) + "] " + str(p)+"% Pos last added {0}".format(last_pos)
#             myThread.displayFeatureProgress(4, displayString)
#         except (pd.errors.EmptyDataError,pd.io.common.EmptyDataError) as _:
#             #print ('! Empty RNA-seq data file',rnaseq_file,i)
#             continue
#
#     return real_values

# Read and save features from separate gzipped files
# def readFeatures(dnasechipseq_list_file,chromhmm_list_file,cage_dir,rnaseq_dir,chromhmm_num_states,input_positions):
def readFeatures(
    dnasechipseq_list_file,
    chromhmm_list_file,
    chromhmm_num_states,
    input_positions,
    frac_feature,
):
    active_indices = defaultdict(
        list
    )  # Key: position index, value: non-zero feature index for the given position index
    # cage_file = os.listdir(cage_dir)[0]
    # df = pd.read_table(cage_dir+cage_file,engine='c',header=None).as_matrix()
    # features = df[:,1:] # Presence of CAGE peak in each position in each cell-type

    # Feature indices are calcualted pre-emptively for multi-threading to append the correct active_indices values at position keys
    # cage_index = chromhmm_index + chromhmm_num_states*len(os.listdir(chromhmm_list_file))
    # rnaseq_index = cage_index + len(features[0]) # len(os.listdir(cage_dir))
    feature_indices = [
        0,
        0,
        len(open(dnasechipseq_list_file).readlines()),
    ]  # , cage_index, rnaseq_index]

    ### Process DNase-seq and ChIP-seq features ###
    # This is the simplest type of feature. We only care about whether there is an overlapping peak at the position
    thread1 = myThread(
        1,
        "Thread-1",
        1,
        dnasechipseq_list_file,
        input_positions,
        active_indices,
        chromhmm_num_states,
        feature_indices,
        frac_feature,
    )

    ### Process ChromHMM features ###
    # We do one-hot encoding for ChromHMM features. For example, if there are 25 states and a position is overlapping
    # with state 5, we have a vector of length 25 with its 5th value set to 1 and the rest set to 0.
    thread2 = myThread(
        2,
        "Thread-2",
        2,
        chromhmm_list_file,
        input_positions,
        active_indices,
        chromhmm_num_states,
        feature_indices,
        frac_feature,
    )

    ### Process CAGE features ###
    # Similarly to DNase-seq and ChIP-seq data, we only care about the presence of a peak for CAGE data.
    # However, the data from different cell-types come in one file so we do not iterate through multiple files here.
    # thread3 = myThread(3, "Thread-3", 3, cage_dir, input_positions, active_indices,
    #                    chromhmm_num_states, feature_indices)

    ### Process RNA-seq features ###
    # thread4 = myThread(4, "Thread-4", 4, rnaseq_dir, input_positions, active_indices,
    #                    chromhmm_num_states, feature_indices)
    # Begin the threads for processing features individualls -- calls run() for each thread instance
    thread1.start()
    thread2.start()
    # thread3.start()
    # thread4.start()

    # Wait for all 4 threads to finish before returning the completed active_indices and real_values from RNAseq feature
    thread1.join()
    thread2.join()
    # thread3.join()
    # thread4.join()

    return active_indices  # ,thread4.real_values


# Write features for each region
def writeFormattedFeatures(f, chrs, starts, ends, position_indices, active_indices, frac_feature):  # ,real_values):
    for i in range(len(position_indices)):  # Iterate through each region to write its features
        ### Output format ###
        # Each lines: chr start end pos_index | <list of all active/non-zero feature indices> | <list real values>
        # Lists are separated by tab
        # Real values correspond to the last active/non-zero feature indices. For example, if there are three real
        # values, they correspond to the three last active/non-zero feature indices. If there is nothing written after
        # the second bar (|), there is no feature with real values. All active features are binary in this case.

        # First, write chromosome, start, end, and position index for the sample
        tp = "\t".join([str(s) for s in [chrs[i], starts[i], ends[i], position_indices[i]]])
        f.write(tp)  # tp.encode())
        window_size = int(ends[i] - starts[i])

        # Second, write active/non-zero feature indices
        if frac_feature:
            d = dict(sorted(active_indices[position_indices[i]]))
            # print(i, tp, d)
            a, b = d.keys(), d.values()
            tp = (
                "\t|\t"
                + "\t".join(["%d" % s for s in a])
                + "\t|\t"
                + "\t".join(["%.6f" % (s / window_size) for s in b])
            )
        else:
            t = sorted(list(set(active_indices[position_indices[i]])))
            tp = "\t|\t" + "\t".join([str(s) for s in t])
        f.write(tp)  # gzf.write(tp.encode())

        # # Third, write real/non-binary values (RNA-seq expression levels)
        # t = real_values[position_indices[i]] # Real values
        # tp = '\t|\t'+'\t'.join([str(round(s,10)) for s in t])
        # gzf.write(tp.encode())

        f.write("\n")  # gzf.write('\n'.encode()) # nextline

        # p = int((i+1)/len(position_indices)*100)
        # sys.stdout.write("\r\t[" + "=" * p + " "*(100-p) + "] " + str(p)+"%")
        # sys.stdout.flush()

    # sys.stdout.write('\n')


def main():

    ### Input arguments ###
    parser = argparse.ArgumentParser(
        prog="base_maker",
        description="Aggregate species-specific features for each training/test sample",
    )
    parser.add_argument("-p", "--position-filename", type=str, help="path to position file")
    # parser.add_argument('-ca','--cage-dir',type=str,help='path to directory containing processed CAGE feature file')
    parser.add_argument(
        "-ch",
        "--chromhmm-list-file",
        type=str,
        help="path to file with list of ChromHMM feature files",
    )
    parser.add_argument(
        "-dn",
        "--dnasechipseq-list-file",
        type=str,
        help="path to file with list of DNase-seq and ChIP-seq feature files",
    )
    # parser.add_argument('-rn','--rnaseq-dir',type=str,help='path to directory containing processed RNA-seq feature file')
    parser.add_argument("-c", "--chrom", type=str, help="chromosome to generate data for", default="all")
    parser.add_argument(
        "-chn",
        "--chromhmm-num-states",
        type=int,
        help="number of ChromHMM states (currently: 25 for human, 15 for mouse)",
    )
    parser.add_argument(
        "-f",
        "--frac-feature",
        help="use fraction of overlap as features instead of binary",
        action="store_true",
    )
    parser.add_argument(
        "-w",
        "--window-size",
        type=int,
        help="size of genomic windows in bp",
    )
    parser.add_argument("-o", "--output-filename", type=str, help="path to output file")
    parser.add_argument(
        "-n",
        "--num-features",
        type=int,
        help="total number of features",
    )
    args = parser.parse_args()
    # print (args)

    ### Read positions ###
    df = pd.read_table(args.position_filename, header=None, names=["chr", "start", "end", "index"])
    if args.chrom != "all":
        df = df[df["chr"] == args.chrom]
        # print('\t%s positions to process in %s' % (len(df), args.chrom))
    df = df.sort_values(by="index").values
    chrs = df[:, 0]  # chromsomes
    starts = df[:, 1]  # start positions
    ends = df[:, 2]  # end positions
    position_indices = df[:, 3]  # position indices
    # print('\t%s positions to process in' % len(position_indices))

    ### Read features from multiple files and write to one output file ###
    # Read features
    active_indices = readFeatures(
        args.dnasechipseq_list_file,
        args.chromhmm_list_file,
        args.chromhmm_num_states,
        set(position_indices),
        args.frac_feature,
    )

    # Write to output
    with open(args.output_filename, "w") as f:
        writeFormattedFeatures(
            f,
            chrs,
            starts,
            ends,
            position_indices,
            active_indices,
            args.frac_feature,
        )


main()
