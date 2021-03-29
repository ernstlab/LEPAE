import os 

### Load config
configfile: "parameters.yaml"
workdir: config["paths"]["working_dir"] 

path_chrom_size_file = config["paths"]["chrom_size_file"]
dist = config["dist"]
seed = config["seed"]

test_chrom = config["test"]["chrom"]
test_start = config["test"]["start"]
test_end = config["test"]["end"]
window_size = config["test"]["window_size"]

tuning_chrom = config["tuning"]["chrom"]

training_data_size = config["training_data_size"]
tuning_data_size = config["tuning_data_size"]
num_trees = config["num_trees"]
random_search_n = config["random_search_n"]

output_dir = config["paths"]["output_dir"]

dnase_chipseq_input_dir = config["paths"]["feature"]["dnase_chipseq_input_dir"]
chromhmm_input_dir = config["paths"]["feature"]["chromhmm_input_dir"]

binary_bed_annot_dir = config["paths"]["binary_bed_annot_dir"]

dnase_chipseq_output_dir = os.path.join(output_dir,"feature/intersect/hg19_DNaseChIPseq")
chromhmm_output_dir = os.path.join(output_dir,"feature/intersect/hg19_ChromHMM")



model_id = "dist%d" % dist

training_prefix = "training_not%s%s" % (tuning_chrom, test_chrom)
tuning_prefix = "tuning_%s" % tuning_chrom
test_prefix = "test_%s" % test_chrom 
prediction_prefix = "prediction_%s" % test_chrom

prefixes = [training_prefix, tuning_prefix, test_prefix, prediction_prefix]

dnase_chipseq_experiments, = glob_wildcards("%s/{experiment}.gz" % dnase_chipseq_input_dir)
chromhmm_experiments, = glob_wildcards("%s/{experiment}.gz" % chromhmm_input_dir)

dnase_chipseq_experiments = sorted(dnase_chipseq_experiments)
chromhmm_experiments = sorted(chromhmm_experiments)

num_features = len(dnase_chipseq_experiments)+len(chromhmm_experiments)*25

gene_features=["genebody","cds","tss","exon","5utr","3utr","intron"]



# Rules to run locally and not as jobs to the cluster 
localrules: all, clean, generate_positions, run_intersect_dnase_chipseq, prepare_input_to_generate_data, split_data, find_best_hyperparam

# Target rule
rule all:
    input:
        os.path.join(output_dir, "prediction/%s.interact" % prediction_prefix),
        expand(os.path.join(output_dir,"figure/{usage}_performance.png"), usage=prefixes[:3]),
        os.path.join(output_dir,"validation/%s_genefeature.txt.gz" % prediction_prefix),
        os.path.join(output_dir,"figure/%s_genefeature.png" % prediction_prefix)

# Rule to run if starting over
rule clean:
    shell: 
        "rm -fr %s" % output_dir



### Pipeline rules

# Sample pairs of genomic bases for training, tuning hyperparamters, testing, and prediction.
# Test data size is assumed to be the same as tuning data size.
# Each line in output files ending with .pair.gz is formatted as follows:
# <chrom> <start> <end> <pair index starting from 1>
# Each line in output files ending with .base.gz is formatted as follows:
# <chrom> <start> <end> <base index --which is pair index followed by "a" for the upstream base and "b" for the downstream base>
rule generate_positions:
    input:
        path_chrom_size_file
    params:
        dist,
        tuning_chrom,
        test_chrom,
        test_start,
        test_end,
        window_size,
        training_data_size,
        tuning_data_size,
        seed, 
        os.path.join(output_dir,"position/%s" % training_prefix),
        os.path.join(output_dir,"position/%s" % tuning_prefix),
        os.path.join(output_dir,"position/%s" % test_prefix),
        os.path.join(output_dir,"position/%s" % prediction_prefix)
    output:
        expand(os.path.join(output_dir,"position/{prefix}.{pb}.gz"), prefix=prefixes, pb=['pair','base'])
    shell:
        """
        . /u/local/Modules/default/init/modules.sh; module load bedtools
        source/generatePositions.sh {input} {params} {output}
        """

# Given sampled bases from generate_positions, run BedTools intersect to determine if they overlap peak calls from DNase-seq and ChIP-seq experiments.
# This is run for each pair of a base set (target) and an experiment.
# Output file consists of indices corresponding to bases that overlap a peak call in the experiment.
rule run_intersect_dnase_chipseq:
    input: 
        position_file = os.path.join(output_dir,"position/{target}.base.gz"),
        experiment_file = os.path.join(dnase_chipseq_input_dir,"{experiment}.gz")     
    output:
        temp(os.path.join(dnase_chipseq_output_dir,"{target}/{experiment}.gz"))
    shell: 
        """
        . /u/local/Modules/default/init/modules.sh; module load bedtools
        bedtools intersect -sorted -a {input.position_file} -b {input.experiment_file} | cut -f 4 | gzip > {output}
        """

# Given sampled bases from generate_positions, run BedTools intersect to determine which chromatin state they overlap with in each epigenome.
# This is run for each pair of a base set (target) and an epigenome (e.g. E001). 
# Output file has two columns, one corresponding to the base index and the other to the chromatin state overlapping the base.
rule run_intersect_chromhmm:
    input: 
        position_file = os.path.join(output_dir,"position/{target}.base.gz"),
        experiment_file = os.path.join(chromhmm_input_dir,"{experiment}.gz")      
    output:
        temp(os.path.join(chromhmm_output_dir,"{target}/{experiment}.gz"))
    shell: 
        """
        . /u/local/Modules/default/init/modules.sh; module load bedtools
        bedtools intersect -sorted -wb -a {input.position_file} -b {input.experiment_file} | cut -f 4,8- | gzip > {output}
        """

# List filenames outputted by BedTools intersect to provide as input to generate_data.
rule prepare_input_to_generate_data:
    input:
        dnase_chipseq_files = sorted(expand(os.path.join(dnase_chipseq_output_dir,"{target}/{experiment}.gz"), experiment=dnase_chipseq_experiments, allow_missing=True)),
        chromhmm_files = sorted(expand(os.path.join(chromhmm_output_dir,"{target}/{experiment}.gz"), experiment=chromhmm_experiments, allow_missing=True))
    output:
        os.path.join(dnase_chipseq_output_dir, "{target}/files.txt"),
        os.path.join(chromhmm_output_dir, "{target}/files.txt")
    run:
        with open(output[0],'w') as f:
            for i in input.dnase_chipseq_files:
                f.write(i+'\n') 
        with open(output[1],'w') as f:
            for i in input.chromhmm_files:
                f.write(i+'\n')

# Generate feature vectors for sampled bases by aggregating output files from BedTools intersect.
# Each line in the .base.gz output file is formatted as follows: 
# <chrom> <start> <end> <base index> | <list of feature indices that should be set to 1>
rule generate_data:
    input:
        position_file = os.path.join(output_dir,"position/{target}.base.gz"),
        dnase_chipseq_list_file = os.path.join(dnase_chipseq_output_dir,"{target}/files.txt"),
        chromhmm_list_file = os.path.join(chromhmm_output_dir,"{target}/files.txt"),
        dnase_chipseq_files = expand(os.path.join(dnase_chipseq_output_dir,"{target}/{experiment}.gz"), experiment=dnase_chipseq_experiments, allow_missing=True),
        chromhmm_files = expand(os.path.join(chromhmm_output_dir,"{target}/{experiment}.gz"), experiment=chromhmm_experiments, allow_missing=True)
    params:
        chromhmm_num_states=25,
        num_features=num_features,
        # data_output_dir=data_output_dir
    threads:
        2
    output: 
        os.path.join(output_dir,"data/{target}.base.gz")
    shell:
        """
        python source/generateDataThreaded.py \
        -p {input.position_file} -ch {input.chromhmm_list_file} -dn {input.dnase_chipseq_list_file} \
        -chn {params.chromhmm_num_states} -o {output} -fn {params.num_features}
        """

# Split feature data into base a vs. base b files such that:
# i-th line in .a.gz file + i-th line in .b.gz file = i-th positive pair
rule split_data:
    input: 
        os.path.join(output_dir,"data/{target}.base.gz")
    params:
        output_prefix = os.path.join(output_dir,"data/{target}"),
        seed = seed
    output:
        expand(os.path.join(output_dir,"data/{target}.{ab}.gz"), ab=['a','b'], allow_missing=True) 
    shell:
        "source/splitData.sh {input} {params}"

# Shuffle lines in .b.gz file and store them in .c.gz file to generate negative pairs such that:
# i-th line in .a.gz file + i-th line in .c.gz file = i-th negative pair
rule shuf_data:
    input: 
        os.path.join(output_dir,"data/{target}.b.gz")
    params:
        seed
    output:
        os.path.join(output_dir,"data/{target}.c.gz")
    shell:
        "source/shufData.sh {input} {params} {output}"

# Random search for random forest classifier hyperparameters.
# Each run with a unique seed trains a random forest with randomly chosen hyperparameters.
# Output file prints out the chosen hyperparameters and performance metrics (MSE, AUROC, AUPRC, etc.).
# For each positive or negative pair, a pair with its upstream base and downstream base are flipped is also provided during training.
# This doubles the number of training and tuning data overall. 
rule hyperparam_search:
    input:
        A = os.path.join(output_dir,"data/%s.a.gz" % training_prefix),
        B = os.path.join(output_dir,"data/%s.b.gz" % training_prefix),
        C = os.path.join(output_dir,"data/%s.c.gz" % training_prefix),
        a = os.path.join(output_dir,"data/%s.a.gz" % tuning_prefix),
        b = os.path.join(output_dir,"data/%s.b.gz" % tuning_prefix),
        c = os.path.join(output_dir,"data/%s.c.gz" % tuning_prefix)
    params:
        # random_search_n = 50,
        num_features = num_features,
        chunk_size = 100,
        training_data_size = training_data_size,
        tuning_data_size = tuning_data_size,
        num_trees = num_trees,
        seed = "{seed}"
    output: 
        os.path.join(output_dir,"classifier/hyperparam_search_%s_seed={seed}_output.txt" % model_id)
    threads:
        12
    shell:
        """
        python -u source/train.py \
        -A {input.A} -B {input.B} -C {input.C} -a {input.a} -b {input.b} -c {input.c} \
        -n {params.num_features} -z {params.chunk_size} -s {params.seed} -t {params.num_trees} \
        -tr {params.training_data_size} -va {params.tuning_data_size} \
        -k \
        > {output}
        """ 

# Find the best combination of hyperparameters given the output files from all runs of hyperparam_search.
rule find_best_hyperparam:
    input:
        expand(os.path.join(output_dir,"classifier/hyperparam_search_%s_seed={seed}_output.txt" % model_id), seed=range(random_search_n))
    output:
       os.path.join(output_dir,"classifier/best_hyperparam_%s.txt" % model_id)
    run:
        import pandas as pd 

        df = pd.read_table(input[0], comment="#")
        for i in range(1,len(input)):
            df = pd.concat([df, pd.read_table(input[i], comment='#')])
        df = df.sort_values(by='vauroc',ascending=False)

        with open(output[0],'w') as f:
            for p in ['md','mss','msl','maf']:
                f.write('-%s\n%s\n' % (p, str(df[p].iloc[0])))
            if df['bs'].iloc[0]:
                f.write('-bs\n')

# Train and save the final random forest classifier with the best combination of hyperparameters.
# As done in hyperparam_search, flipped pairs are added to training and tuning data, doubling the data size.
rule train:
    input:
        hyperparam_file = os.path.join(output_dir,"classifier/best_hyperparam_%s.txt" % model_id),
        A = os.path.join(output_dir,"data/%s.a.gz" % training_prefix),
        B = os.path.join(output_dir,"data/%s.b.gz" % training_prefix),
        C = os.path.join(output_dir,"data/%s.c.gz" % training_prefix),
        a = os.path.join(output_dir,"data/%s.a.gz" % tuning_prefix),
        b = os.path.join(output_dir,"data/%s.b.gz" % tuning_prefix),
        c = os.path.join(output_dir,"data/%s.c.gz" % tuning_prefix)
    params:
        seed = seed,
        num_features = num_features,
        chunk_size = 100,
        training_data_size = training_data_size,
        tuning_data_size = tuning_data_size,
        num_trees = num_trees
    output:
        progress_file = os.path.join(output_dir,"classifier/train_%s_output.txt" % model_id),
        classifier_file = os.path.join(output_dir,"classifier/RF_%s.pkl" % model_id)
    threads:
        12
    shell:
        """
        python -u source/train.py -A {input.A} -B {input.B} -C {input.C} -a {input.a} -b {input.b} -c {input.c} \
        -n {params.num_features} -z {params.chunk_size} -s {params.seed} -t {params.num_trees} \
        -tr {params.training_data_size} -va {params.tuning_data_size} \
        -v -o {output.classifier_file} @{input.hyperparam_file} \
        > {output.progress_file}
        """

# Make predictions for training/tuning/test positive and negative pairs.     
# As with rule predict, it outputs two lines for every pair, one for the input pair and the other for its flipped version.
rule evaluate:
    input:
        classifier_file = os.path.join(output_dir,"classifier/RF_%s.pkl" % model_id),
        a = os.path.join(output_dir,"data/{usage}.a.gz"),
        b = os.path.join(output_dir,"data/{usage}.b.gz"),
        c = os.path.join(output_dir,"data/{usage}.c.gz")
    params:
        num_features = num_features
    output:
        positive_score_file = os.path.join(output_dir,"prediction/{usage}.positive.txt.gz"),
        negative_score_file = os.path.join(output_dir,"prediction/{usage}.negative.txt.gz")
    shell:
        """
        python source/predict.py -t {input.classifier_file} -A {input.a} -B {input.b} -n {params.num_features} -o {output.positive_score_file}
        python source/predict.py -t {input.classifier_file} -A {input.a} -B {input.c} -n {params.num_features} -o {output.negative_score_file}
        """

# Plot ROC curve, precision-recall curve, and scatterplot of predictions made for flipped pairs to evaluate model performance.
# This is based on predictions made for a set of positive and negative pairs for training, tuning, or test.
rule plot_performance:
    input:
        positive_score_file = os.path.join(output_dir,"prediction/{usage}.positive.txt.gz"),
        negative_score_file = os.path.join(output_dir,"prediction/{usage}.negative.txt.gz")
    output:
        os.path.join(output_dir,"figure/{usage}_performance.png")
    run:
        import matplotlib.pyplot as plt, numpy as np, scipy.stats as ss, gzip
        from sklearn.metrics import roc_curve, precision_recall_curve, roc_auc_score, average_precision_score

        p = []
        n = []
        with gzip.open(input.positive_score_file,'rt') as f:
            for line in f:
                p.append(float(line.strip()))
        with gzip.open(input.negative_score_file,'rt') as f:
            for line in f:
                n.append(float(line.strip()))
        y_true = [1]*len(p) + [0]*len(n)
        y_score = p + n

        fig, axarr = plt.subplots(1,3, figsize=(9,3.5), dpi=200)
        
        ax = axarr[0]
        fpr, tpr, _ = roc_curve(y_true, y_score)
        ax.plot(fpr, tpr, color='k')
        ax.text(0.95, 0.05, 'AUROC: %.5f' % roc_auc_score(y_true, y_score), horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)
        ax.set_xlabel('False positive rate')
        ax.set_ylabel('True positive rate')
        ax.set_title('ROC curve')

        ax = axarr[1]
        precision, recall, _ = precision_recall_curve(y_true, y_score)
        ax.plot(recall, precision, color='k')
        ax.text(0.05, 0.05, 'AUPRC: %.5f' % average_precision_score(y_true, y_score), horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
        ax.set_xlabel('Recall')
        ax.set_ylabel('Precision')
        ax.set_title('Precision-recall curve')

        ax = axarr[2]
        a = np.concatenate([p[::2],n[::2]])
        b = np.concatenate([p[1::2],n[1::2]])
        ax.scatter(a,b,color='k',alpha=0.1,s=1)
        ax.text(0.95,0.05,s='PCC: %.5f' % ss.pearsonr(a,b)[0],horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)
        ax.set_xlabel('A-B pair score')
        ax.set_ylabel('B-A pair score')
        ax.set_title('Agreement in flipped pairs')

        plt.tight_layout()
        plt.savefig(output[0])

# Make predictions for the chromosome held out specifically for prediction.
# Odd lines in the output file correspond to predictions made for the pairs provided as input.
# Even lines in the output file correspond to predictions made for the flipped pairs.
# Every two lines correspond to pairs consisting of the same bases, except the second line is for the flipped pair.
rule predict:
    input:
        classifier_file = os.path.join(output_dir,"classifier/RF_%s.pkl" % model_id),
        prediction_a = os.path.join(output_dir,"data/%s.a.gz" % prediction_prefix),
        prediction_b = os.path.join(output_dir,"data/%s.b.gz" % prediction_prefix)
    params:
        num_features = num_features
    output:
        os.path.join(output_dir,"prediction/%s.txt.gz" % prediction_prefix)
    threads:
        12
    shell:
        "python source/predict.py -t {input.classifier_file} -A {input.prediction_a} -B {input.prediction_b} -n {params.num_features} -o {output}"

# Generate an .intersect browser track given the predictions for the chromosome held out specifically for prediction.
# For every pair, we use the average of its two predictions, corresponding to every two lines in the output of predict.
# Pairs with score above the score_threshold is shown in the genome browser as an arc.
# This threshold can be adjusted in the browser. 
rule generate_track:
    input:
        os.path.join(output_dir,"data/%s.a.gz" % prediction_prefix),
        os.path.join(output_dir,"data/%s.b.gz" % prediction_prefix),
        os.path.join(output_dir,"prediction/%s.txt.gz" % prediction_prefix)
    params:
        score_threshold = 950,
        score_name = "pairwise_dist%d" % dist
    output:
        os.path.join(output_dir,"prediction/%s.interact" % prediction_prefix)
    shell:
        """
        . /u/local/Modules/default/init/modules.sh; module load bedtools
        source/generateTrack.sh {input} {params.score_threshold} {params.score_name} {output}
        """

rule run_intersect_gene_annot:
    input: 
        interact_file = os.path.join(output_dir,"prediction/%s.interact" % prediction_prefix),
        gene_annot_bed_files = expand(os.path.join(binary_bed_annot_dir,'hg19.GENCODE_{gene_feature}.bed.gz'), gene_feature=gene_features)
    params:
        gene_features = gene_features
    output:
        os.path.join(output_dir,"validation/%s_genefeature.txt.gz" % prediction_prefix)
    shell:
        """
        . /u/local/Modules/default/init/modules.sh; module load bedtools

        cat {input.interact_file} | tail -n +3 | cut -d"\t" -f 9- > {output}

        for f in {input.gene_annot_bed_files}; do
            bedtools intersect -a {output} -b $f -sorted -c > {output}.tmp
            mv {output}.tmp {output}
        done

        cut -d"\t" -f 6- {output} > {output}.tmp
        mv {output}.tmp {output}

        for f in {input.gene_annot_bed_files}; do
            bedtools intersect -a {output} -b $f -sorted -c > {output}.tmp
            mv {output}.tmp {output}
        done

        paste -d"\t" <(cat {input.interact_file} | tail -n +3 | awk -v OFS="\t" '{{print $9,$10,$11,$14,$15,$16,$6}}') <(cat {output} | cut -d"\t" -f 6-) | gzip > {output}.tmp
        printf "#a.chrom\ta.start\ta.end\tb.chrom\tb.start\tb.end\tSCORE" | gzip > {output}
        for gf in {params.gene_features}; do 
            printf "\ta."$gf | gzip >> {output}
        done
        for gf in {params.gene_features}; do 
            printf "\tb."$gf | gzip >> {output}
        done
        printf "\n" | gzip >> {output}
        cat {output}.tmp >> {output}
        rm {output}.tmp
        """

rule plot_intersect_gene_annot:
    input:
        os.path.join(output_dir,"validation/%s_genefeature.txt.gz" % prediction_prefix)
    params:
        gene_features = gene_features,
        gene_feature_names = ["gene body","CDS","TSS","exon","5' UTR","3' UTR","intron"]
    output:
        os.path.join(output_dir,"figure/%s_genefeature.png" % prediction_prefix)
    run:
        import pandas as pd, matplotlib.pyplot as plt
        
        df = pd.read_table(input[0])

        fig,axarr = plt.subplots(2,2, figsize=(8,6), dpi=200)
        tp = [[df[(df['a.'+gf]>0) | (df['b.'+gf]>0)]['SCORE'] for gf in params.gene_features], 
              [df[df['a.'+gf]>0][df['b.'+gf]>0]['SCORE'] for gf in params.gene_features]]
        gw_mean = df['SCORE'].mean()
        gw_median = df['SCORE'].median()

        t = ['at least\none base','both\nbases']
        for i in range(2):
            ax = axarr[i,0]
            ax.boxplot(tp[i], vert=False, showfliers=False)
            # ax.plot([gw_mean]*2, [-1,len(params.gene_features)+1], '--', color='tab:green', zorder=-99, lw=0.8)
            ax.plot([gw_median]*2, [-1,len(params.gene_features)+1], '--', color='tab:orange', zorder=-99, lw=0.8)
            ax.set_xlim([0,1])
            ax.set_xlabel('score')
            ax.set_ylim([0.5,len(params.gene_features)+0.5])
            ax.set_yticklabels(params.gene_feature_names)
            ax.invert_yaxis()
            ax.set_title('Score distribution of %d-kb pairs with %s overlapping gene annotations' % (int(dist/1000),t[i]), fontsize='medium')

            ax = axarr[i,1]
            for j in range(len(params.gene_features)):
                x = len(tp[i][j])/len(df)*100
                ax.barh(y=j+1, width=x,color='grey')
                ax.text(x=x, y=j+1, s='%.1f' % x if x>0.005 else '%.3f' % x, ha='left', va='center', fontsize='small')
            ax.set_xlim([0,100])
            ax.set_xlabel('percentage')
            ax.set_ylim([0.5,len(params.gene_features)+0.5])
            ax.set_yticks(range(1,len(params.gene_features)+1))
            ax.set_yticklabels(params.gene_feature_names)
            ax.invert_yaxis()
            ax.set_title('Percentage of %d-kb pairs with %s overlapping gene annotations' % (int(dist/1000),t[i]), fontsize='medium')

        plt.tight_layout()
        plt.savefig(output[0])