import sys,gzip,random,numpy as np,argparse,joblib

def printProgress(p):
    sys.stdout.write("\r[" + "=" * p + " "*(100-p) + "]" + str(p)+"%")
    sys.stdout.flush()

def predict(clf,a_test_filename,b_test_filename,output_filename,test_data_size,batch_size,num_features):

    # Make predictions and write to output
    with gzip.open(a_test_filename,'r') as a_test_file,\
            gzip.open(b_test_filename,'r') as b_test_file,\
            gzip.open(output_filename if output_filename.endswith('.gz') else output_filename+'.gz','w') as fout: # open test data file

        for i in range(int(test_data_size/batch_size)+1): # iterate through each batch
            current_batch_size = batch_size if i<int(test_data_size/batch_size) else test_data_size%batch_size # last batch may have items fewer than N
            if current_batch_size==0:
                break
            X = np.zeros((current_batch_size*2,num_features*2),dtype=bool)

            for j in range(current_batch_size): # iterate through each sample within the batch
                l1 = a_test_file.readline().strip().split(b'|')
                l2 = b_test_file.readline().strip().split(b'|')

                # Read indices of features that should be set to 1
                positive_nonzero_feature_indices_1 = [int(s) for s in l1[1].strip().split()] + [num_features + int(s) for s in l2[1].strip().split()]
                X[2*j, positive_nonzero_feature_indices_1] = True

                # Same but regions symmetrically flipped
                positive_nonzero_feature_indices_2 = [num_features + int(s) for s in l1[1].strip().split()] + [int(s) for s in l2[1].strip().split()]
                X[2*j+1, positive_nonzero_feature_indices_2] = True

            # Make prediction on current batch
            y_pred = clf.predict_proba(X)[:,1]

            # Write and mouse genomic positions and predicted probabilities of the current batch
            sample_output = [str(round(y_pred[j],7)) for j in range(len(y_pred))] # append position data to the front
            l = '\n'.join(sample_output)+'\n'
            fout.write(l.encode())

            # printProgress(int((i+1)/(int(test_data_size/batch_size)+1)*100))

def main():

    parser = argparse.ArgumentParser(prog='base_maker',description='Predict score given a trained neural network')
    parser.add_argument('-t', '--trained-classifier-filename', help='path to trained classifier (.pkl or .pkl.gz)', type=str)
    parser.add_argument('-A','--a-test-data-filename',help='path to region A test data file',type=str)
    parser.add_argument('-B','--b-test-data-filename',help='path to region B test data file',type=str)

    # parser.add_argument('-d', '--test-data-size', help='number of samples in test data', type=int, default=5000)
    parser.add_argument('-s', '--seed', help='random seed', type=int, default=1)
    parser.add_argument('-o', '--output-filename', help='path to output file', type=str)
    parser.add_argument('-b', '--batch-size', help='batch size', type=int, default=1000)

    parser.add_argument('-n','--num-features',help='number of features in input vector',type=int,default=6819)
    # parser.add_argument('-mf','--num-mouse-features',help='number of mouse features in input vector',type=int,default=3113)
    # parser.add_argument('-hrmin','--rnaseq-min',help='minimum expression level in RNA-seq data',type=float,default=8e-05)
    # parser.add_argument('-hrmax','--rnaseq-max',help='maximum expression level in RNA-seq data',type=float,default=1.11729e06)

    args = parser.parse_args()

    # rnaseq_range = [args.rnaseq_min,args.rnaseq_max]
    # mouse_rnaseq_range = [args.mouse_rnaseq_min,args.mouse_rnaseq_max]

    # Set random seed
    random.seed(args.seed)
    np.random.seed(args.seed)

    test_data_size = len(gzip.open(args.a_test_data_filename,'rb').readlines())

    # print ("Loading trained classifier...")
    # Load previously trained classifier
    clf = joblib.load(args.trained_classifier_filename)

    # print ("Making predictions...")
    predict(clf,args.a_test_data_filename,args.b_test_data_filename,args.output_filename,test_data_size,args.batch_size,args.num_features)

if __name__ == "__main__":
    main()

# time python -u predict.py -t ../model/dist10000_size100000/dist10000_size100000_1_1_5000_5000_512_32_16_None_True.pkl -A ../data/all_dist10000_size100000.base.test.a.gz -B ../data/all_dist10000_size100000.base.test.b.gz -o ../model/dist10000_size100000/positive_test_score.txt.gz; time python -u predict.py -t ../model/dist10000_size100000/dist10000_size100000_1_1_5000_5000_512_32_16_None_True.pkl -A ../data/all_dist10000_size100000.base.test.a.gz -B ../data/all_dist10000_size100000.base.test.c.gz -o ../model/dist10000_size100000/negative_test_score.txt.gz


# time python -u predict.py -t ../model/dist20000_size100000/dist20000_size100000_1_1_5000_5000_512_64_8_4096_False.pkl -A ../data/all_dist20000_size100000.base.test.a.gz -B ../data/all_dist20000_size100000.base.test.b.gz -o ../model/dist20000_size100000/dist20000_size100000/positive_test_score.txt.gz; time python -u predict.py -t ../model/dist20000_size100000/dist20000_size100000_1_1_5000_5000_512_64_8_4096_False.pkl -A ../data/all_dist20000_size100000.base.test.a.gz -B ../data/all_dist20000_size100000.base.test.c.gz -o ../model/dist20000_size100000/negative_test_score.txt.gz

# time python -u predict.py -t ../model/dist30000_size100000/dist30000_size100000_1_1_5000_5000_512_64_8_4096_False.pkl -A ../data/all_dist30000_size100000.base.test.a.gz -B ../data/all_dist30000_size100000.base.test.b.gz -o ../model/dist30000_size100000/positive_test_score.txt.gz; time python -u predict.py -t ../model/dist30000_size100000/dist30000_size100000_1_1_5000_5000_512_64_8_4096_False.pkl -A ../data/all_dist30000_size100000.base.test.a.gz -B ../data/all_dist30000_size100000.base.test.c.gz -o ../model/dist30000_size100000/negative_test_score.txt.gz





# time python -u predict.py -t ../model/dist10000_size100000/dist10000_size100000_1_1_5000_5000_512_32_16_None_True.pkl -A ../data/chr1_dist10000.base.a.gz -B ../data/chr1_dist10000.base.b.gz -o ../prediction/chr1_dist10000_score.txt.gz -d 50000

# time python -u predict.py -t ../model/dist20000_size100000/dist20000_size100000_1_1_5000_5000_512_64_8_4096_False.pkl -A ../data/chr1_dist20000.base.a.gz -B ../data/chr1_dist20000.base.b.gz -o ../prediction/chr1_dist20000_score.txt.gz -d 50000

# time python -u predict.py -t ../model/dist30000_size100000/dist30000_size100000_1_1_5000_5000_512_64_8_4096_False.pkl  -A ../data/chr1_dist30000.base.a.gz -B ../data/chr1_dist30000.base.b.gz -o ../prediction/chr1_dist30000_score.txt.gz -d 50000