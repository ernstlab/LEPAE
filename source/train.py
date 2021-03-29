import argparse,math,sys,scipy.sparse,gzip,random,numpy as np
#from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc, average_precision_score, mean_squared_error
import joblib

# Train and save random forest classifier

def printProgress(p):
    sys.stdout.write("\r\t[%s%s] %.3f%%    " % ("=" * int(p), " " * (100 - int(p)), p))
    sys.stdout.flush()

def setTrainingHyperParameters(seed, num_features):
    random.seed(seed)
    np.random.seed(seed)

    max_depth = random.choice([16,32,64,128,256,512,1024,2048,None]) # md (n=9)
    min_samples_split = random.choice([32,64,128,256]) # mss (n=4)
    min_samples_leaf = random.choice([32,64,128,256]) # msl (n=4)
    max_features = random.choice(["auto","sqrt","log2"]) # 32,64,128,256,512,1024,2048,4096,None]) # maf (n=9)
    bootstrap = True #random.choice([True,False]) # (n=2)

    return max_depth,min_samples_split,min_samples_leaf,max_features,bootstrap

def readChunk(files_to_read,current_chunk_size,num_features): #,rnaseq_range):

    # File pointers
    positive_a_data_file = files_to_read[0]
    positive_b_data_file = files_to_read[1]
    negative_a_data_file = files_to_read[2]
    negative_b_data_file = files_to_read[3]

    # Lists storing information needed for building a SciPy CSR matrix for the current chunk
    indptr = [0]
    indices = []
    data = []

    # Positive/negative label
    label = []

    for _ in range(current_chunk_size): # iterate through each pair of one positive sample and one negative sample

        ### Read two positive samples
        l1 = positive_a_data_file.readline().strip().split(b'|')
        l2 = positive_b_data_file.readline().strip().split(b'|')

        # print(current_chunk_size, i, lenO)

        # Read indices of features that should be set to 1
        positive_nonzero_feature_indices_1 = [int(s) for s in l1[1].strip().split()]
        positive_nonzero_feature_indices_2 = [num_features + int(s) for s in l2[1].strip().split()]

        # Same but regions symmetrically flipped
        positive_nonzero_feature_indices_3 = [num_features + int(s) for s in l1[1].strip().split()]
        positive_nonzero_feature_indices_4 = [int(s) for s in l2[1].strip().split()]

        # Normalize RNA-seq values
        # positive_real_valued_features = [(float(s)-rnaseq_range[0])/(rnaseq_range[1]-rnaseq_range[0]) for s in hl[2].strip().split()] if len(hl)>1 else []
        # positive_real_valued_mouse_features = [(float(s)-mouse_rnaseq_range[0])/(mouse_rnaseq_range[1]-mouse_rnaseq_range[0]) for s in ml[2].strip().split()] if len(ml)>1 else []

        # Save data for two positive samples
        indptr += [indptr[-1]+len(positive_nonzero_feature_indices_1)+len(positive_nonzero_feature_indices_2)]
        indptr += [indptr[-1]+len(positive_nonzero_feature_indices_3)+len(positive_nonzero_feature_indices_4)]
        indices += positive_nonzero_feature_indices_1 + positive_nonzero_feature_indices_2 + positive_nonzero_feature_indices_3 + positive_nonzero_feature_indices_4
        data += [1] * (len(positive_nonzero_feature_indices_1)+len(positive_nonzero_feature_indices_2)+len(positive_nonzero_feature_indices_3)+len(positive_nonzero_feature_indices_4))
        #data += [1]*(len(positive_nonzero_feature_indices)-len(positive_real_valued_features)) + positive_real_valued_features\
        #        + [1]*(len(positive_nonzero_mouse_feature_indices)-len(positive_real_valued_mouse_features)) + positive_real_valued_mouse_features
        label += [1,1]



        ### Read two negative samples
        l1 = negative_a_data_file.readline().strip().split(b'|')
        l2 = negative_b_data_file.readline().strip().split(b'|')

        negative_nonzero_feature_indices_1 = [int(s) for s in l1[1].strip().split()]
        negative_nonzero_feature_indices_2 = [num_features + int(s) for s in l2[1].strip().split()]
        negative_nonzero_feature_indices_3 = [num_features + int(s) for s in l1[1].strip().split()]
        negative_nonzero_feature_indices_4 = [int(s) for s in l2[1].strip().split()]

        # Normalize RNA-seq values
        # negative_real_valued_features = [(float(s)-rnaseq_range[0])/(rnaseq_range[1]-rnaseq_range[0]) for s in hl[2].strip().split()] if len(hl)>1 else []
        # negative_real_valued_mouse_features = [(float(s)-mouse_rnaseq_range[0])/(mouse_rnaseq_range[1]-mouse_rnaseq_range[0]) for s in ml[2].strip().split()] if len(ml)>1 else []

        # Save data for two negative samples
        indptr += [indptr[-1]+len(negative_nonzero_feature_indices_1)+len(negative_nonzero_feature_indices_2)]
        indptr += [indptr[-1]+len(negative_nonzero_feature_indices_3)+len(negative_nonzero_feature_indices_4)]
        indices += negative_nonzero_feature_indices_1 + negative_nonzero_feature_indices_2 + negative_nonzero_feature_indices_3 + negative_nonzero_feature_indices_4
        data += [1] * (len(negative_nonzero_feature_indices_1) + len(negative_nonzero_feature_indices_2) + len(negative_nonzero_feature_indices_3) + len(negative_nonzero_feature_indices_4))
        # data += [1]*(len(negative_nonzero_feature_indices)-len(negative_real_valued_features)) + negative_real_valued_features\
        #         + [1]*(len(negative_nonzero_mouse_feature_indices)-len(negative_real_valued_mouse_features)) + negative_real_valued_mouse_features
        label += [0,0]

    return scipy.sparse.csr_matrix((data,indices,indptr),shape=(current_chunk_size*2*2,num_features*2),dtype=bool),label

def csr_vappend(a,b): # vertically combines two Scipy CSR matrices
    a.data = np.hstack((a.data,b.data))
    a.indices = np.hstack((a.indices,b.indices))
    a.indptr = np.hstack((a.indptr,(b.indptr + a.nnz)[1:]))
    a._shape = (a.shape[0]+b.shape[0],b.shape[1])
    return a

def readData(files_to_read,data_size,chunk_size,num_features): #,rnaseq_range):

    # Calculate the number of chucks
    num_chunks = int(math.ceil(data_size/chunk_size))

    # Data array that will eventually store all feature data
    X = None
    Y = []

    for i in range(num_chunks): # iterate through the number of chunks
        current_chunk_size = chunk_size if i <int(data_size/chunk_size) else data_size%chunk_size
        if current_chunk_size==0:
            break
        #print (i,current_chunk_size,i*chunk_size,i*chunk_size+current_chunk_size,lines_to_read[:,i*chunk_size:i*chunk_size+current_chunk_size])
        tmp_X,tmp_Y = readChunk(files_to_read,current_chunk_size,num_features) #,rnaseq_range)

        # Store constructed feature matrix
        if X==None:
            X = tmp_X
            Y = tmp_Y
        else:
            X = csr_vappend(X,tmp_X)
            Y += tmp_Y
        # printProgress((i+1)/num_chunks*100)
    # print (X.shape)

    return X,Y

def eval(clf,inputs,labels):
    scores = clf.predict_proba(inputs)
    scores = scores[:,1]
    mse = mean_squared_error(labels, scores)
    fpr, tpr, _ = roc_curve(labels, scores)
    auroc = auc(fpr, tpr) # au-ROC
    auprc = average_precision_score(labels, scores) # au-PRC
    pos_mean_score = np.mean([scores[i] for i in range(len(scores)) if labels[i]==1])
    neg_mean_score = np.mean([scores[i] for i in range(len(scores)) if labels[i]==0])

    return mse,auroc,auprc,np.mean(scores),pos_mean_score,neg_mean_score

def main():

    ### Input arguments ###
    parser = argparse.ArgumentParser(prog='base_maker',description='Train a random forest',fromfile_prefix_chars='@')

    # Training data filenames
    parser.add_argument('-A','--a-training-data-filename',help='path to region A training data file',type=str)
    parser.add_argument('-B','--b-training-data-filename',help='path to region B training data file',type=str)
    parser.add_argument('-C','--c-training-data-filename',help='path to shuffled region B training data file',type=str)

    # Tuning data filenames
    parser.add_argument('-a','--a-tuning-data-filename',help='path to region A tuning data file',type=str)
    parser.add_argument('-b','--b-tuning-data-filename',help='path to region B tuning data file',type=str)
    parser.add_argument('-c','--c-tuning-data-filename',help='path to shuffled region B tuning data file',type=str)


    # Output prefix
    parser.add_argument('-o','--output-pickle-filename',type=str,default='tmp.pkl')

    # Options
    parser.add_argument('-k','--random-search',help='if hyperparameters should be randomly set',action='store_true')
    parser.add_argument('-v','--save',help='if the trained model should be saved after training',action='store_true')
    parser.add_argument('-d','--sample-data',help='if training data should be sampled',action='store_true')

    # Hyperparameters for training
    parser.add_argument('-r','--neg-data-ratio',help='ratio of negative samples to positive samples',type=int,default=1)
    parser.add_argument('-s', '--seed', help='random seed', type=int, default=1)
    # parser.add_argument('-w', '--random-search-n', help='number of hyperparameter combinations to try during random search', type=int, default=10)
    parser.add_argument('-z','--chunk-size',help='number of samples to read at a time before building a sparse matrix',type=int,default=1000)
    parser.add_argument('-tr','--positive-training-data-size',help='number of samples in positive training data',type=int,default=50000)
    parser.add_argument('-va','--positive-tuning-data-size',help='number of samples in positive tuning data',type=int,default=5000)
    parser.add_argument('-t','--num-trees',help='number of trees in a random forest',type=int,default=10)


    # Feature information --fixed
    parser.add_argument('-n','--num-features',help='number of features in input vector',type=int,default=6819)
    #parser.add_argument('-i','--rnaseq-min',help='minimum expression level in RNA-seq data',type=float,default=8e-05)
    #parser.add_argument('-x','--rnaseq-max',help='maximum expression level in RNA-seq data',type=float,default=1.11729e06)

    # Hyperparameters for random forest
    parser.add_argument('-md','--max-depth',type=str,default="None")
    parser.add_argument('-mss','--min-samples-split',type=int,default=2)
    parser.add_argument('-msl','--min-samples-leaf',type=int, default=32)
    parser.add_argument('-maf','--max_features',type=str,default="auto")
    parser.add_argument('-bs', '--bootstrap', action='store_true')

    args = parser.parse_args()
    # print ('\nReading input parameters...')
    for key, value in vars(args).items():
        if 'max_' not in key and 'min_' not in key:
            print ('# %s:'%key,value)

    # Define range for raw RNA-seq values
    #rnaseq_range = [args.rnaseq_min,args.rnaseq_max]

    # Read training data
    # print ('\nReading training data...')
    files_to_read = [gzip.open(args.a_training_data_filename, 'rb'),
                     gzip.open(args.b_training_data_filename, 'rb'),
                     gzip.open(args.a_training_data_filename, 'rb'),
                     gzip.open(args.c_training_data_filename, 'rb')]
    training_inputs,training_labels = readData(files_to_read,args.positive_training_data_size,args.chunk_size,args.num_features)#, rnaseq_range)

    # print('\nReading tuning data...')
    files_to_read = [gzip.open(args.a_tuning_data_filename, 'rb'),
                     gzip.open(args.b_tuning_data_filename, 'rb'),
                     gzip.open(args.a_tuning_data_filename, 'rb'),
                     gzip.open(args.c_tuning_data_filename, 'rb')]
    tuning_inputs, tuning_labels = readData(files_to_read, args.positive_tuning_data_size, args.chunk_size, args.num_features)# ,rnaseq_range)

    # print('\nResult')
    print('seed\tnr\tmd\tmss\tmsl\tmaf\tbs\ttmse\ttauroc\ttauprc\ttmean\ttpmean\ttnmean\tvmse\tvauroc\tvauprc\tvmean\tvpmean\tvnmean')

    if args.random_search:
        max_depth, min_samples_split, min_samples_leaf, max_features, bootstrap = setTrainingHyperParameters(args.seed, args.num_features)
    else:
        max_depth, min_samples_split, min_samples_leaf, max_features, bootstrap = args.max_depth, args.min_samples_split, args.min_samples_leaf, args.max_features, args.bootstrap
        # max_features = None if max_features=='None' else max_features
        max_depth = None if max_depth=='None' else int(max_depth)

    hyperparameters = [args.seed, args.neg_data_ratio, max_depth, min_samples_split, min_samples_leaf, max_features, bootstrap]

    ####### Training starts here #######
    #print ('\nTraining...')
    clf = RandomForestClassifier(random_state=args.seed,
                                 n_estimators=args.num_trees,
                                 max_depth=max_depth,
                                 class_weight={0:args.neg_data_ratio,1:1},
                                 min_samples_split=min_samples_split,
                                 min_samples_leaf=min_samples_leaf,
                                 max_features=max_features,
                                 bootstrap=bootstrap,
                                 n_jobs=-1)
    clf.fit(training_inputs,training_labels)
    ####### Training ends here #######

    # Evaluate
    training_result = eval(clf,training_inputs,training_labels)
    tuning_result = eval(clf,tuning_inputs,tuning_labels)

    # Print out result
    print('\t'.join([str(s) for s in hyperparameters + [round(s,6) for s in training_result+tuning_result]]))

    # Save model if specified
    if args.save:
        fn = args.output_pickle_filename if args.output_pickle_filename.endswith(".pkl") else args.output_pickle_filename+'.pkl'
        joblib.dump(clf,fn)

main()