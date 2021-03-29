import sys,gzip,threading,pandas as pd,argparse,os
from collections import defaultdict

# SOURCE: https://www.tutorialspoint.com/python3/python_multithreading.htm
class myThread (threading.Thread):
    featureStrings = {key: str() for key in range(1,3)} #5)}
    printInit = 0
    
    def __init__(self, threadID, name, feature, directory, input_positions, active_indices, 
                 chrom_states, feature_indices):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.feature = feature
        self.directory = directory
        self.input_positions = input_positions
        self.active_indices = active_indices
        self.chrom_states = chrom_states
        self.featIndex = feature_indices[feature]
        self.real_values = None
      
    def run(self):
        if self.feature == 1:
            readDnaseChipFeature(self.directory, self.input_positions, self.active_indices)
        elif self.feature == 2:
            readChromHmmFeature(self.directory, self.input_positions, self.active_indices, 
                                self.chrom_states, self.featIndex)
        # elif self.feature == 3:
        #     readCageFeature(self.directory, self.input_positions, self.active_indices,
        #                     self.chrom_states, self.featIndex)
        # elif self.feature == 4:
        #     self.real_values = readRnaSeqFeature(self.directory, self.input_positions,
        #                                        self.active_indices, self.featIndex)

    # def displayFeatureProgress(feature, stringUpdate):
    #     if myThread.printInit:
    #         sys.stdout.write(4*"\033[F")
    #     else:
    #         myThread.printInit = 1
    #     myThread.featureStrings[feature] = stringUpdate
    #     displayString = str()
    #     for i in range(1,3): # 5):
    #         displayString += myThread.featureStrings[i] + "\n"
    #     sys.stdout.write(displayString)
    #     sys.stdout.flush()
      
def readDnaseChipFeature(dnase_chipseq_list_file, input_positions, active_indices):
    # List of files containing position indices with overlapping peak in different DNase-seq and ChIP-seq data
    dnase_chipseq_files = open(dnase_chipseq_list_file).readlines() # sorted(os.listdir(dnase_chipseq_list_file))
    # Iterate through each file (each experiment in a specific cell-type and with a specific target if it's ChIP-seq)
    for i in range(len(dnase_chipseq_files)):
        dnase_chipseq_file = dnase_chipseq_files[i].strip()
        if os.stat(dnase_chipseq_file).st_size == 0:
            # print ('! Empty DNase/ChIP-seq data file:',dnase_chipseq_file,i)
            continue
        try:
            positions = pd.read_table(dnase_chipseq_file, engine='c', header=None, squeeze=True).tolist()
            valid = list(input_positions.intersection(positions))
            for j in range(len(valid)):
                active_indices[valid[j]].append(i)
        except pd.errors.EmptyDataError as _: # Some input file may be empty
            # print ('! Empty DNase/ChIP-seq data file:',dnase_chipseq_file,i)
            continue
        # Status output
        # p = int((i+1)/len(dnase_chipseq_files)*100)
        # displayString = "\tDNase-seq and ChIP-seq [" + "=" * p + " "*(100-p) + "] " + str(p)+"%\t"
        # myThread.displayFeatureProgress(1, displayString)

def readChromHmmFeature(chromhmm_list_file, input_positions, active_indices, chromhmm_num_states, featIndex):
    num_current_features = featIndex
    #    active_indices = defaultdict(list) # Key: position index, value: non-zero feature index for the given position index
    
    # List of files containing position indices and their overlapping ChromHMM state for each cell-type
    chromhmm_files = open(chromhmm_list_file).readlines() # sorted(os.listdir(chromhmm_list_file))
    for i in range(len(chromhmm_files)): # Iterate through each file (each cell-type)
        chromhmm_file = chromhmm_files[i].strip()
        if os.stat(chromhmm_file).st_size == 0: # Some input file may be empty
            # print ('! Empty ChromHMM data file',chromhmm_file,i)
            continue
        with gzip.open(chromhmm_file,'rb') as f:
            positions_found = 0
            for line in f:
                l = line.strip().split()
                position = l[0].decode('utf-8')
                if position in input_positions:
                    positions_found += 1
                    state = l[1].decode('utf-8')
                    state = int(state[1:]) if state.startswith('U') else int(state)
                    active_indices[position].append(num_current_features+(state-1))

                # ChromHMM runtime optimization added 10 July 2019
                if positions_found == len(input_positions):
                    break
            num_current_features += chromhmm_num_states
        # Status output
        # p = int((i+1)/len(chromhmm_files)*100)
        # displayString = "\tChromHMM [" + "=" * p + " "*(100-p) + "] " + str(p)+"%"
        # myThread.displayFeatureProgress(2, displayString)
    
# def readCageFeature(cage_dir, input_positions, active_indices, chromhmm_num_states, featIndex):
# #    num_current_features = 0 # Number of features processed so far
#     num_current_features = featIndex
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
# def readRnaSeqFeature(rnaseq_dir, input_positions, active_indices, featIndex):
#     num_current_features = featIndex
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
#def readFeatures(dnase_chipseq_list_file,chromhmm_list_file,cage_dir,rnaseq_dir,chromhmm_num_states,input_positions):
def readFeatures(dnase_chipseq_list_file,chromhmm_list_file,chromhmm_num_states,input_positions):
    active_indices = defaultdict(list) # Key: position index, value: non-zero feature index for the given position index
    # cage_file = os.listdir(cage_dir)[0]
    # df = pd.read_table(cage_dir+cage_file,engine='c',header=None).as_matrix()
    # features = df[:,1:] # Presence of CAGE peak in each position in each cell-type
    
    # Feature indices are calcualted pre-emptively for multi-threading to append the correct active_indices values at position keys
    chromhmm_index = len(open(dnase_chipseq_list_file).readlines())
    # cage_index = chromhmm_index + chromhmm_num_states*len(os.listdir(chromhmm_list_file))
    # rnaseq_index = cage_index + len(features[0]) # len(os.listdir(cage_dir))
    feature_indices = [0, 0, chromhmm_index]#, cage_index, rnaseq_index]

    ### Process DNase-seq and ChIP-seq features ###
    # This is the simplest type of feature. We only care about whether there is an overlapping peak at the position
    thread1 = myThread(1, "Thread-1", 1, dnase_chipseq_list_file, input_positions, active_indices, 
                       chromhmm_num_states, feature_indices)    

    ### Process ChromHMM features ###
    # We do one-hot encoding for ChromHMM features. For example, if there are 25 states and a position is overlapping
    # with state 5, we have a vector of length 25 with its 5th value set to 1 and the rest set to 0.
    thread2 = myThread(2, "Thread-2", 2, chromhmm_list_file, input_positions, active_indices, 
                       chromhmm_num_states, feature_indices)

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

    return active_indices#,thread4.real_values

# Write features for each region
def writeFormattedFeatures(gzf,chrs,starts,ends,position_indices,active_indices): #,real_values):
    for i in range(len(position_indices)): # Iterate through each region to write its features
        ### Output format ###
        # Each lines: chr start end pos_index | <list of all active/non-zero feature indices> | <list real values>
        # Lists are separated by tab
        # Real values correspond to the last active/non-zero feature indices. For example, if there are three real
        # values, they correspond to the three last active/non-zero feature indices. If there is nothing written after
        # the second bar (|), there is no feature with real values. All active features are binary in this case.

        # First, write chromosome, start, end, and position index for the sample
        tp = '\t'.join([str(s) for s in [chrs[i],starts[i],ends[i],position_indices[i]]])
        gzf.write(tp.encode())

        # Second, write active/non-zero feature indices
        t = active_indices[position_indices[i]]
        tp = '\t|\t'+'\t'.join([str(s) for s in t])
        gzf.write(tp.encode())

        # # Third, write real/non-binary values (RNA-seq expression levels)
        # t = real_values[position_indices[i]] # Real values
        # tp = '\t|\t'+'\t'.join([str(round(s,10)) for s in t])
        # gzf.write(tp.encode())

        gzf.write('\n'.encode()) # nextline

        # p = int((i+1)/len(position_indices)*100)
        # sys.stdout.write("\r\t[" + "=" * p + " "*(100-p) + "] " + str(p)+"%")
        # sys.stdout.flush()

    # sys.stdout.write('\n')

# Added 11 July 2019 for multi-threading feature -- each thread can add features at position keys out of sequence; sorting just for output    
def sortActiveIndicesByPosition(active_indices):
    for key in active_indices.keys():
        active_indices[key].sort()

def main():

    ### Input arguments ###
    parser = argparse.ArgumentParser(prog='base_maker',description='Aggregate species-specific features for each training/test sample')
    parser.add_argument('-p','--position-filename',type=str,help='path to position file')
    # parser.add_argument('-ca','--cage-dir',type=str,help='path to directory containing processed CAGE feature file')
    parser.add_argument('-ch','--chromhmm-list-file',type=str,help='path to file with list of ChromHMM feature files')
    parser.add_argument('-dn','--dnase-chipseq-list-file',type=str,help='path to file with list of DNase-seq and ChIP-seq feature files')
    # parser.add_argument('-rn','--rnaseq-dir',type=str,help='path to directory containing processed RNA-seq feature file')
    parser.add_argument('-s','--split',action='store_true',help='whether to split the data (when submitting multiple jobs to parallelize)')
    parser.add_argument('-c','--split-chunk-size',default=1000000,type=int,help='size of each chunk/split if splitting data')
    parser.add_argument('-i','--split-index',type=int,help='split index (starts from 1)')
    parser.add_argument('-chn','--chromhmm-num-states',type=int,help='number of ChromHMM states (currently: 25 for human, 15 for mouse)')
    parser.add_argument('-o','--output-filename',type=str,help='path to output file')
    parser.add_argument('-fn','--num-features',type=int,help='total number of features (currently: 8824 for human, 3313 for mouse)')
    args = parser.parse_args()
    # print (args)


    ### Read positions ###
    # print ('Reading positions...')
    if args.split:
        split_start = (args.split_index-1)*args.split_chunk_size
        split_end = (args.split_index)*args.split_chunk_size
        # print ('\tPosition index ranges from %d to %d' % (split_start,split_end))

    df = pd.read_table(args.position_filename,engine='c',header=None,usecols=[0,1,2,3],names=['chr','start','end','index']).values
    df = df[df[:,3].argsort()]
    if args.split:
        chrs = df[split_start:split_end,0] # chromsomes
        starts = df[split_start:split_end,1] # start positions
        ends = df[split_start:split_end,2] # end positions
        position_indices = df[split_start:split_end,3] # position indices
    else:
        chrs = df[:,0] # chromsomes
        starts = df[:,1] # start positions
        ends = df[:,2] # end positions
        position_indices = df[:,3] # position indices
    # print('\t%s positions read' % len(position_indices))


    ### Read features from multiple files and write to one output file ###
    with gzip.open(args.output_filename if args.output_filename.endswith('.gz') else args.output_filename+'.gz','wb') as gzf:

        # Read features
        # print ('Reading features...')
        active_indices = readFeatures(args.dnase_chipseq_list_file,args.chromhmm_list_file,args.chromhmm_num_states,set(position_indices))
        #active_indices,real_values = readFeatures(args.dnase_chipseq_list_file,args.chromhmm_list_file,args.cage_dir,args.rnaseq_dir,args.chromhmm_num_states,set(position_indices))

        positions_read = set(active_indices.keys())
        overlapping_positions = positions_read.intersection(position_indices)
        # print ('\tOverlapping positions: %d' % len(overlapping_positions))
        
        # Threading causes features to finish out of order, so we need to re-order for consistency
        sortActiveIndicesByPosition(active_indices)

        # Write formatted features
        # print ('Writing formatted features...')
        writeFormattedFeatures(gzf,chrs,starts,ends,position_indices,active_indices) #,real_values)

main()

# time python generateDataThreaded.py -p ../position/all_dist10000_size100000.base.gz -ch ../feature/intersect/dist10000_size100000/hg19_ChromHMM/ -dn ../feature/intersect/dist10000_size100000/hg19_DNaseChIPseq/ -chn 25 -o ../data/all_dist10000_size100000.base.gz -fn 6819

# time python generateDataThreaded.py -p ../position/all_dist20000_size100000.base.gz -ch ../feature/intersect/dist20000_size100000/hg19_ChromHMM/ -dn ../feature/intersect/dist20000_size100000/hg19_DNaseChIPseq/ -chn 25 -o ../data/all_dist20000_size100000.base.gz -fn 6819

# time python generateDataThreaded.py -p ../position/all_dist30000_size100000.base.gz -ch ../feature/intersect/dist30000_size100000/hg19_ChromHMM/ -dn ../feature/intersect/dist30000_size100000/hg19_DNaseChIPseq/ -chn 25 -o ../data/all_dist30000_size100000.base.gz -fn 6819



# time python generateDataThreaded.py -p ../position/chr1_dist10000.base.gz -ch ../feature/intersect/chr1_dist10000/hg19_ChromHMM/ -dn ../feature/intersect/chr1_dist10000/hg19_DNaseChIPseq/ -chn 25 -o ../data/chr1_dist10000.base.gz -fn 6819

# time python generateDataThreaded.py -p ../position/chr1_dist20000.base.gz -ch ../feature/intersect/chr1_dist20000/hg19_ChromHMM/ -dn ../feature/intersect/chr1_dist20000/hg19_DNaseChIPseq/ -chn 25 -o ../data/chr1_dist20000.base.gz -fn 6819

# time python generateDataThreaded.py -p ../position/chr1_dist30000.base.gz -ch ../feature/intersect/chr1_dist30000/hg19_ChromHMM/ -dn ../feature/intersect/chr1_dist30000/hg19_DNaseChIPseq/ -chn 25 -o ../data/chr1_dist30000.base.gz -fn 6819