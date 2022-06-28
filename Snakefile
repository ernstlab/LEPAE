import os, random


### Load config
configfile: "parameters_1kb.yaml"


# Window size
window_size = config["window_size"]

# Pairwise distances between windows
min_dist, max_dist, dist_incr = config["min_dist"], config["max_dist"], window_size
all_dists = list(range(min_dist, max_dist + 1, dist_incr))

# Seed for randomness
seed = config["seed"]
random.seed(seed)

# Chromosomes
all_chroms = ["chr" + str(c) for c in range(1, 23)] + ["chrX"]
odd_chroms = ["chr" + str(c) for c in range(1, 23, 2)]
even_chroms = ["chr" + str(c) for c in range(2, 23, 2)]

num_tuning_chrom = 3
training_chrom = dict()
tuning_chrom = dict()
test_chrom = dict()
pred_chrom = dict()

heldout_chroms = "odd"
sampled_even_chroms = list(random.sample(even_chroms, num_tuning_chrom))
training_chrom[heldout_chroms] = [
    c for c in even_chroms if c not in sampled_even_chroms
]
tuning_chrom[heldout_chroms] = sampled_even_chroms
test_chrom[heldout_chroms] = odd_chroms
pred_chrom[heldout_chroms] = odd_chroms + ["chrX"]

heldout_chroms = "even"
sampled_odd_chroms = list(random.sample(odd_chroms, num_tuning_chrom))
training_chrom[heldout_chroms] = [c for c in odd_chroms if c not in sampled_odd_chroms]
tuning_chrom[heldout_chroms] = sampled_odd_chroms
test_chrom[heldout_chroms] = even_chroms
pred_chrom[heldout_chroms] = even_chroms

# print(pred_chrom)

# for heldout_chrom in all_chroms:
#     sampled_chrom = list(random.sample([c for c in all_chroms if c != heldout_chrom and c != "chrX"], num_tuning_chrom))
#     training_chrom[heldout_chrom] = [c for c in all_chroms if c != heldout_chrom and c != "chrX" and c not in sampled_chrom]
#     tuning_chrom[heldout_chrom] = sampled_chrom
sorted_chroms = [
    "chr" + str(c)
    for c in [1] + list(range(10, 20)) + [2, 20, 21, 22] + list(range(3, 10))
] + ["chrX"]


# Training parameters
training_data_size = config["training_data_size"]
tuning_data_size = config["tuning_data_size"]
num_random_search = config["num_random_search"]
ensemble_size = config["ensemble_size"]

# Decision tree parameters
DT_hyperparam = config["DT"]["hyperparam"]

# Neural network parameters
NN_hyperparam = config["NN"]["hyperparam"]

# Whether to use fraction of feature coverage as feature value (instead of binary)
# frac_feature = config["frac_feature"]


# Paths to input files/directories
path_chrom_size_file = config["paths"]["chrom_size_file"]
dnasechipseq_experiment_paths_file = config["paths"][
    "dnasechipseq_experiment_paths_file"
]
dnasechipseq_data_dir = config["paths"]["dnasechipseq_data_dir"]
chromhmm_experiment_paths_file = config["paths"]["chromhmm_experiment_paths_file"]
chromhmm_data_dir = config["paths"]["dnasechipseq_data_dir"]


# Input features
num_dnasechipseq_experiments = config["num_dnasechipseq_experiments"]
num_chromhmm_experiments = config["num_chromhmm_experiments"]
num_features = (
    num_dnasechipseq_experiments
    + num_chromhmm_experiments * config["num_chromhmm_states"]
)
dnasechipseq_experiments = dict()
with open(dnasechipseq_experiment_paths_file, "r") as f:
    for i in range(num_dnasechipseq_experiments):
        path = f.readline().strip()
        experiment = path.split("/")[-1].replace(".bed.gz", "")
        dnasechipseq_experiments[experiment] = path
chromhmm_experiments = dict()
with open(chromhmm_experiment_paths_file, "r") as f:
    for i in range(num_chromhmm_experiments):
        path = f.readline().strip()
        experiment = path.split("/")[-1].replace(".bed.gz", "")
        chromhmm_experiments[experiment] = path


# Paths to output directories
output_dir = config["paths"]["output_dir"]
dnasechipseq_output_dir = os.path.join(
    output_dir, "feature/intersect/hg38_DNaseChIPseq"
)
chromhmm_output_dir = os.path.join(output_dir, "feature/intersect/hg38_ChromHMM")

# Output file prefixes
pred_prefix = "prediction_min%d_max%d_window%d" % (min_dist, max_dist, window_size)


### Pipeline rules
# Rules to run locally and not as jobs to the cluster
localrules:
    all,
    clean,
    prepare_input_to_generate_data,
    find_best_hyperparam,
    aggregate_bedpe,
    generate_juicer_short,


wildcard_constraints:
    chrom="|".join(sorted_chroms),
    classifier="|".join(["NN", "DT"]),
    method="|".join(["NN", "DT", "Jaccard"]),
    heldout_chroms="|".join(["odd", "even"]),
    norm="|".join(["NONE", "VC_SQRT"]),


# Target rule
rule all:
    input:
        os.path.join(output_dir, "prediction/NN/%s.hic" % pred_prefix),
        expand(
            os.path.join(
                output_dir,
                "classifier/{classifier}/{heldout_chroms}_heldout/train_dist{dist}.test_pred_{ver}.txt.gz",
            ),
            classifier=['DT','NN'],
            heldout_chroms=["odd", "even"],
            dist=all_dists,
            ver=['ensembled_smoothed','ensembled'] + ['seed%d' % s for s in range(ensemble_size)]
        ),
        expand(
            os.path.join(
                output_dir,
                "classifier/Jaccard/{heldout_chroms}_heldout/train_dist{dist}.test_pred.txt.gz",
            ),
            heldout_chroms=["odd", "even"],
            dist=all_dists,
        ),
        os.path.join(output_dir, "prediction/NN/%s.bedpe.gz" % pred_prefix),
        os.path.join(output_dir, "prediction/DT/%s.bedpe.gz" % pred_prefix),
        os.path.join(output_dir, "prediction/Jaccard/%s.bedpe.gz" % pred_prefix),
        os.path.join(output_dir, "prediction/hic/VC_SQRT.txt.gz"),


# Rule to run if starting over
rule clean:
    shell:
        "rm -fr %s" % output_dir


# Define genomic windows across the human genome with a specified window size
# Each line in the output file all_windows.bed.gz is tab delimited and formatted as follows:
# <chrom> <start> <end> <window index>
rule generate_windows:
    input:
        path_chrom_size_file,
    params:
        window_size,
    output:
        os.path.join(output_dir, "window/all_windows.bed.gz"),
    shell:
        """
        . /u/local/Modules/default/init/modules.sh; module load bedtools
        bedtools makewindows -g {input} -w {params} | sort -k1,1 -k2,2n | awk -v OFS="\t" '{{print $1,$2,$3,NR}}' | gzip > {output}
        """


# Given windows from generate_windows, run BedTools intersect to determine if they overlap peak calls from DNase-seq and ChIP-seq experiments.
# This is run for each experiment, outputting a file that contains indices of windows that overlap a peak call in the experiment.
rule run_intersect_dnasechipseq:
    input:
        window_file=os.path.join(output_dir, "window/all_windows.bed.gz"),
    params:
        experiment_path=lambda wildcards: dnasechipseq_experiments[
            wildcards.experiment
        ],
    output:
        os.path.join(dnasechipseq_output_dir, "{experiment}.bed.gz"),
    shell:
        """
        . /u/local/Modules/default/init/modules.sh; module load bedtools
        echo {params.experiment_path}
        curl -L -s {params.experiment_path} | gzip -cd | cut -f 1-3 | sort -k1,1 -k2,2n |\
        bedtools intersect -sorted -wo -a {input.window_file} -b - | cut -f 4,8 | sort -V | uniq | gzip > {output}
        """


# Given windows from generate_windows, run BedTools intersect to determine which chromatin state they overlap with in each epigenome.
# Output file has two columns, one corresponding to the window index and the other to the chromatin state overlapping the window in the epigenome.
rule run_intersect_chromhmm:
    input:
        window_file=os.path.join(output_dir, "window/all_windows.bed.gz"),
    params:
        experiment_path=lambda wildcards: chromhmm_experiments[wildcards.experiment],
    output:
        os.path.join(chromhmm_output_dir, "{experiment}.bed.gz"),
    shell:
        """
        . /u/local/Modules/default/init/modules.sh; module load bedtools
        echo {params.experiment_path}
        curl -L -s {params.experiment_path} | gzip -cd | cut -f 1-4 | sort -k1,1 -k2,2n |\
        bedtools intersect -sorted -wo -a {input.window_file} -b - | cut -f 4,8- | sort -V | uniq | gzip > {output}
        """


# List filenames outputted by BedTools intersect to provide as input to generate_data.
rule prepare_input_to_generate_data:
    input:
        dnasechipseq_files=sorted(
            expand(
                os.path.join(
                    dnasechipseq_output_dir,
                    "{experiment}.bed.gz",
                ),
                experiment=dnasechipseq_experiments.keys(),
                allow_missing=True,
            )
        ),
        chromhmm_files=sorted(
            expand(
                os.path.join(
                    chromhmm_output_dir,
                    "{experiment}.bed.gz",
                ),
                experiment=chromhmm_experiments.keys(),
                allow_missing=True,
            )
        ),
    output:
        os.path.join(dnasechipseq_output_dir, "files.txt"),
        os.path.join(chromhmm_output_dir, "files.txt"),
    run:
        with open(output[0], "w") as f:
            for i in input.dnasechipseq_files:
                f.write(i + "\n")
        with open(output[1], "w") as f:
            for i in input.chromhmm_files:
                f.write(i + "\n")


# Generate feature vectors for sampled windows in each chromosome by aggregating output files from BedTools intersect.
# This is run for each chromosome. Each line in the output file ending with .window.gz is formatted as follows:
# <chrom> <start> <end> <window index> | <list of feature indices that should be set to 1>
rule generate_data:
    input:
        window_file=os.path.join(output_dir, "window/all_windows.bed.gz"),
        dnasechipseq_list_file=os.path.join(dnasechipseq_output_dir, "files.txt"),
        chromhmm_list_file=os.path.join(chromhmm_output_dir, "files.txt"),
        dnasechipseq_files=expand(
            os.path.join(
                dnasechipseq_output_dir,
                "{experiment}.bed.gz",
            ),
            experiment=dnasechipseq_experiments.keys(),
            allow_missing=True,
        ),
        chromhmm_files=expand(
            os.path.join(
                chromhmm_output_dir,
                "{experiment}.bed.gz",
            ),
            experiment=chromhmm_experiments.keys(),
            allow_missing=True,
        ),
    params:
        chromhmm_num_states=25,
        num_features=num_features,
        # frac_feature="-f" if frac_feature else " ",
        window_size=window_size,
    threads: 2
    output:
        os.path.join(output_dir, "data/{chrom}.window"),
    shell:
        """
        time python3 source/generateDataThreaded.py \
        -p {input.window_file} -ch {input.chromhmm_list_file} -dn {input.dnasechipseq_list_file} \
        -chn {params.chromhmm_num_states} -n {params.num_features} -c {wildcards.chrom} \
        -w {params.window_size} \
        -o {output}
        """


# Random search for classifier hyperparameters.
# Output file ending with .progress.txt reports a combination of hyperparameters and its classification performance metric (MSE, AUROC, AUPRC).
# For each positive or negative pair, a pair with its upstream and downstream windows flipped is also provided as an additional pair, doubling the number of training and tuning data overall.
rule hyperparam_search:
    input:
        expand(os.path.join(output_dir, "data/{chrom}.window"), chrom=all_chroms),
    params:
        num_random_search=num_random_search,
        num_features=num_features,
        training_data_size=training_data_size,
        tuning_data_size=tuning_data_size,
        window_size=window_size,
        output_filename_prefix=lambda w: os.path.join(
            output_dir,
            "classifier/%s/%s_heldout/hyperparam_search_dist%s_seed%s"
            % (w.classifier, w.heldout_chroms, w.dist, w.seed),
        ),
        # frac_feature="-f" if frac_feature else " ",
        training_data_filenames=lambda w: expand(
            os.path.join(output_dir, "data/{chrom}.window"),
            chrom=training_chrom[w.heldout_chroms],
        ),
        tuning_data_filenames=lambda w: expand(
            os.path.join(output_dir, "data/{chrom}.window"),
            chrom=tuning_chrom[w.heldout_chroms],
        ),
    output:
        os.path.join(
            output_dir,
            "classifier/{classifier}/{heldout_chroms}_heldout/hyperparam_search_dist{dist}_seed{seed}.progress.txt",
        ),
    threads: 4
    shell:
        """
        time python3 -u source/train.py {wildcards.classifier} \
        --training-data-filenames {params.training_data_filenames} \
        --tuning-data-filenames {params.tuning_data_filenames} \
        --positive-training-data-size {params.training_data_size} \
        --positive-tuning-data-size {params.tuning_data_size} \
        --pairwise-dist {wildcards.dist} \
        --window-size {params.window_size} \
        --seed {wildcards.seed} \
        --random-search \
        --num-features {params.num_features} \
        --output-filename-prefix {params.output_filename_prefix}
        """


# Given classification performance metrics from multiple classifiers, output the hyperparameter values of the best performing one
rule find_best_hyperparam:
    input:
        expand(
            os.path.join(
                output_dir,
                "classifier/{classifier}/{heldout_chroms}_heldout/hyperparam_search_dist{dist}_seed{seed}.progress.txt",
            ),
            seed=range(num_random_search),
            allow_missing=True,
        ),
    params:
        hyperparam=lambda w: config[w.classifier]["hyperparam"],
        criteria="vaauroc",
        ascending=False,
    output:
        os.path.join(
            output_dir,
            "classifier/{classifier}/{heldout_chroms}_heldout/hyperparam_search_dist{dist}.best_hyperparam.txt",
        ),
    run:
        import pandas as pd

        df = pd.concat(
            [pd.read_csv(input[i], comment="#", sep="\t") for i in range(len(input))]
        )
        print(df)
        df = df.sort_values(by=params.criteria, ascending=params.ascending)
        with open(output[0], "w") as f:
            for k in params.hyperparam:
                f.write("--%s\n%s\n" % (k.replace("_", "-"), str(df[k].iloc[0])))


# Train and save a classifier with the best combination of hyperparameters.
# As done in hyperparam_search, flipped pairs are added to training and tuning data, doubling the data size.
rule train:
    input:
        best_hyperparam_file=os.path.join(
            output_dir,
            "classifier/{classifier}/{heldout_chroms}_heldout/hyperparam_search_dist{dist}.best_hyperparam.txt",
        ),
    params:
        num_features=num_features,
        training_data_size=training_data_size,
        tuning_data_size=tuning_data_size,
        window_size=window_size,
        output_filename_prefix=lambda w: os.path.join(
            output_dir,
            "classifier/%s/%s_heldout/train_dist%s_seed%s"
            % (w.classifier, w.heldout_chroms, w.dist, w.seed),
        ),
        # frac_feature="-f" if frac_feature else " ",
        training_data_filenames=lambda w: expand(
            os.path.join(output_dir, "data/{chrom}.window"),
            chrom=training_chrom[w.heldout_chroms],
        ),
        tuning_data_filenames=lambda w: expand(
            os.path.join(output_dir, "data/{chrom}.window"),
            chrom=tuning_chrom[w.heldout_chroms],
        ),
        test_data_filenames=lambda w: expand(
            os.path.join(output_dir, "data/{chrom}.window"),
            chrom=test_chrom[w.heldout_chroms],
        ),
    output:
        expand(
            os.path.join(
                output_dir,
                "classifier/{classifier}/{heldout_chroms}_heldout/train_dist{dist}_seed{seed}.{suffix}",
            ),
            suffix=[
                "progress.txt",
                "pkl",
                "training_pred.txt.gz",
                "tuning_pred.txt.gz",
                "test_pred.txt.gz",
            ],
            allow_missing=True,
        ),
    threads: 4
    shell:
        """
        time python3 -u source/train.py {wildcards.classifier} \
        --training-data-filenames {params.training_data_filenames} \
        --tuning-data-filenames {params.tuning_data_filenames} \
        --test-data-filenames {params.test_data_filenames} \
        --positive-training-data-size {params.training_data_size} \
        --positive-tuning-data-size {params.tuning_data_size} \
        --pairwise-dist {wildcards.dist} \
        --window-size {params.window_size} \
        --seed {wildcards.seed} \
        --random-training-data \
        --num-features {params.num_features} \
        --max-num-epoch 20 \
        --save \
        --output-filename-prefix {params.output_filename_prefix} \
        @{input.best_hyperparam_file}
        """


rule eval_jaccard:
    input:
        expand(os.path.join(output_dir, "data/{chrom}.window"), chrom=all_chroms),
    params:
        num_features=num_features,
        training_data_size=training_data_size,
        tuning_data_size=tuning_data_size,
        window_size=window_size,
        output_filename_prefix=lambda w: os.path.join(
            output_dir,
            "classifier/Jaccard/%s_heldout/train_dist%s" % (w.heldout_chroms, w.dist),
        ),
        seed=seed,
        # frac_feature="-f" if frac_feature else " ",
        training_data_filenames=lambda w: expand(
            os.path.join(output_dir, "data/{chrom}.window"),
            chrom=training_chrom[w.heldout_chroms],
        ),
        tuning_data_filenames=lambda w: expand(
            os.path.join(output_dir, "data/{chrom}.window"),
            chrom=tuning_chrom[w.heldout_chroms],
        ),
        test_data_filenames=lambda w: expand(
            os.path.join(output_dir, "data/{chrom}.window"),
            chrom=test_chrom[w.heldout_chroms],
        ),
    output:
        expand(
            os.path.join(
                output_dir,
                "classifier/Jaccard/{heldout_chroms}_heldout/train_dist{dist}.{suffix}_pred.txt.gz",
            ),
            suffix=["test"],  # ["training", "tuning", "test"],
            allow_missing=True,
        ),
    threads: 4
    shell:
        """
        time python3 -u source/train.py Jaccard \
        --training-data-filenames {params.training_data_filenames} \
        --tuning-data-filenames {params.tuning_data_filenames} \
        --test-data-filenames {params.test_data_filenames} \
        --positive-training-data-size {params.training_data_size} \
        --positive-tuning-data-size {params.tuning_data_size} \
        --pairwise-dist {wildcards.dist} \
        --window-size {params.window_size} \
        --seed {params.seed} \
        --num-features {params.num_features} \
        --save \
        --output-filename-prefix {params.output_filename_prefix}
        """

# Make predictions for the chromosome that was held out specifically for prediction.
rule eval_classifier:
    input:
        classifier_filenames=[os.path.join(
                output_dir,
                "classifier/{classifier}/{heldout_chroms}_heldout/train_dist{dist}_seed{seed}.pkl",
            )
        ],
        test_data_filenames=lambda w: expand(
            os.path.join(output_dir, "data/{chrom}.window"),
            chrom=test_chrom[w.heldout_chroms],
        ),
    params:
        classifier_prefix=os.path.join(
            output_dir, "classifier/{classifier}/{heldout_chroms}_heldout/train_dist"
        ),
        num_features=num_features,
        window_size=window_size,
        seed=seed,
        data_size=tuning_data_size,
    output:
        os.path.join(
            output_dir,
            "classifier/{classifier}/{heldout_chroms}_heldout/train_dist{dist}.test_pred_seed{seed}.txt.gz",
        ),
    threads: 4
    shell:
        """
        time python3 -u source/predict.py {wildcards.classifier} \
        -t {input.classifier_filenames} \
        -i {input.test_data_filenames} \
        -d {wildcards.dist} \
        -w {params.window_size} \
        -n {params.num_features} \
        -s {params.seed} \
        -j {params.data_size} \
        --subsample \
        --incl-negative \
        -o {output}
        """

# Make predictions for the chromosome that was held out specifically for prediction.
rule eval_classifier_ensemble:
    input:
        classifier_filenames=lambda w: expand(
            os.path.join(
                output_dir,
                "classifier/{classifier}/{heldout_chroms}_heldout/train_dist{dist}_seed{seed}.pkl",
            ),
            seed=range(ensemble_size),
            allow_missing=True,
        ),
        test_data_filenames=lambda w: expand(
            os.path.join(output_dir, "data/{chrom}.window"),
            chrom=test_chrom[w.heldout_chroms],
        ),
    params:
        classifier_prefix=os.path.join(
            output_dir, "classifier/{classifier}/{heldout_chroms}_heldout/train_dist"
        ),
        num_features=num_features,
        window_size=window_size,
        seed=seed,
        data_size=tuning_data_size,
    output:
        os.path.join(
            output_dir,
            "classifier/{classifier}/{heldout_chroms}_heldout/train_dist{dist}.test_pred_ensembled.txt.gz",
        ),
    threads: 4
    shell:
        """
        time python3 -u source/predict.py {wildcards.classifier} \
        -t {input.classifier_filenames} \
        -i {input.test_data_filenames} \
        -d {wildcards.dist} \
        -w {params.window_size} \
        -n {params.num_features} \
        -s {params.seed} \
        -j {params.data_size} \
        --subsample \
        --incl-negative \
        -o {output}
        """

# Make predictions for the chromosome that was held out specifically for prediction.
rule eval_classifier_ensemble_smooth:
    input:
        classifier_filenames=lambda w: expand(
            os.path.join(
                output_dir,
                "classifier/{classifier}/{heldout_chroms}_heldout/train_dist{train_dist}_seed{seed}.pkl",
            ),
            train_dist=[
                int(w.dist) - window_size,
                int(w.dist),
                int(w.dist) + window_size,
            ]
            if int(w.dist) > window_size and int(w.dist) < max_dist
            else [int(w.dist) - window_size, int(w.dist)]
            if int(w.dist) == max_dist
            else [int(w.dist), int(w.dist) + window_size],
            seed=range(ensemble_size),
            allow_missing=True,
        ),
        test_data_filenames=lambda w: expand(
            os.path.join(output_dir, "data/{chrom}.window"),
            chrom=test_chrom[w.heldout_chroms],
        ),
    params:
        classifier_prefix=os.path.join(
            output_dir, "classifier/{classifier}/{heldout_chroms}_heldout/train_dist"
        ),
        num_features=num_features,
        window_size=window_size,
        seed=seed,
        data_size=tuning_data_size,
    output:
        os.path.join(
            output_dir,
            "classifier/{classifier}/{heldout_chroms}_heldout/train_dist{dist}.test_pred_ensembled_smoothed.txt.gz",
        ),
    threads: 4
    shell:
        """
        time python3 -u source/predict.py {wildcards.classifier} \
        -t {input.classifier_filenames} \
        -i {input.test_data_filenames} \
        -d {wildcards.dist} \
        -w {params.window_size} \
        -n {params.num_features} \
        -s {params.seed} \
        -j {params.data_size} \
        --subsample \
        --incl-negative \
        -o {output}
        """




# Make predictions for the chromosome that was held out specifically for prediction.
rule predict:
    input:
        classifier_filenames=lambda w: expand(
            os.path.join(
                output_dir,
                "classifier/{classifier}/{heldout_chroms}_heldout/train_dist{train_dist}_seed{seed}.pkl",
            ),
            train_dist=[
                int(w.dist) - window_size,
                int(w.dist),
                int(w.dist) + window_size,
            ]
            if int(w.dist) > window_size and int(w.dist) < max_dist
            else [int(w.dist) - window_size, int(w.dist)]
            if int(w.dist) == max_dist
            else [int(w.dist), int(w.dist) + window_size],
            seed=range(ensemble_size),
            allow_missing=True,
        ),
        input_filename=os.path.join(output_dir, "data/{chrom}.window"),
    params:
        classifier_prefix=os.path.join(
            output_dir, "classifier/{classifier}/{heldout_chroms}_heldout/train_dist"
        ),
        num_features=num_features,
        window_size=window_size,
        # frac_feature="-f" if frac_feature else " ",
        seed=seed,
    output:
        os.path.join(
            output_dir,
            "prediction/{classifier}/{heldout_chroms}_heldout/{chrom}/%s_dist{dist}.bedpe.gz"
            % pred_prefix,
        ),
    threads: 4
    shell:
        """
        time python3 -u source/predict.py {wildcards.classifier} \
        -t {input.classifier_filenames} \
        -i {input.input_filename} \
        -d {wildcards.dist} \
        -w {params.window_size} \
        -n {params.num_features} \
        -s {params.seed} \
        -o {output}
        """


# Make predictions for the chromosome that was held out specifically for prediction.
rule jaccard:
    input:
        input_filename=os.path.join(output_dir, "data/{chrom}.window"),
    params:
        num_features=num_features,
        window_size=window_size,
        seed=seed,
        # frac_feature="-f" if frac_feature else " ",
    output:
        os.path.join(
            output_dir,
            "prediction/Jaccard/{heldout_chroms}_heldout/{chrom}/%s_dist{dist}.bedpe.gz"
            % pred_prefix,
        ),
    threads: 4
    shell:
        """
        time python3 -u source/predict.py Jaccard \
        -i {input.input_filename} \
        -d {wildcards.dist} -w {params.window_size} \
        -n {params.num_features} \
        -s {params.seed} -o {output}
        """


# Sort scores for each chromosome
rule sort_bedpe:
    input:
        expand(
            os.path.join(
                output_dir,
                "prediction/{method}/odd_heldout/{chrom}/%s_dist{dist}.bedpe.gz"
                % pred_prefix,
            ),
            chrom=pred_chrom["odd"],
            dist=all_dists,
            allow_missing=True,
        ),
        expand(
            os.path.join(
                output_dir,
                "prediction/{method}/even_heldout/{chrom}/%s_dist{dist}.bedpe.gz"
                % pred_prefix,
            ),
            chrom=pred_chrom["even"],
            dist=all_dists,
            allow_missing=True,
        ),
    params:
        odd_prefix=lambda w: os.path.join(
            output_dir, "prediction/%s/odd_heldout/" % w.method
        ),
        even_prefix=lambda w: os.path.join(
            output_dir, "prediction/%s/even_heldout/" % w.method
        ),
        odd_chroms=list(pred_chrom["odd"]),
        even_chroms=list(pred_chrom["even"]),
        suffix=pred_prefix,
    output:
        expand(
            os.path.join(
                output_dir,
                "prediction/{method}/odd_heldout/{chrom}/%s.bedpe.gz" % pred_prefix,
            ),
            chrom=pred_chrom["odd"],
            allow_missing=True,
        ),
        expand(
            os.path.join(
                output_dir,
                "prediction/{method}/even_heldout/{chrom}/%s.bedpe.gz" % pred_prefix,
            ),
            chrom=pred_chrom["even"],
            allow_missing=True,
        ),
    threads: 8
    shell:
        """
        echo {params}
        for chrom in {params.odd_chroms}; do
            echo $chrom;
            cat {params.odd_prefix}/"$chrom"/{params.suffix}_dist*.bedpe.gz | gzip -cd | sort -k2,2n -k5,5n | gzip > {params.odd_prefix}/"$chrom"/{params.suffix}.bedpe.gz &
        done
        wait

        for chrom in {params.even_chroms}; do
            echo $chrom;
            cat {params.even_prefix}/"$chrom"/{params.suffix}_dist*.bedpe.gz | gzip -cd | sort -k2,2n -k5,5n | gzip > {params.even_prefix}/"$chrom"/{params.suffix}.bedpe.gz &
        done
        wait
        """

# Aggregate .bedpe.gz files into one .bedpe.gz file
rule aggregate_bedpe:
    input:
        expand(
            os.path.join(
                output_dir,
                "prediction/{method}/odd_heldout/{chrom}/%s.bedpe.gz" % pred_prefix,
            ),
            chrom=pred_chrom["odd"],
            allow_missing=True,
        ),
        expand(
            os.path.join(
                output_dir,
                "prediction/{method}/even_heldout/{chrom}/%s.bedpe.gz" % pred_prefix,
            ),
            chrom=pred_chrom["even"],
            allow_missing=True,
        ),
    params:
        sorted_chroms=sorted_chroms,
        prefix=os.path.join(output_dir, "prediction/{method}"),
        suffix="%s.bedpe.gz" % pred_prefix,
    output:
        os.path.join(output_dir, "prediction/{method}/%s.bedpe.gz" % pred_prefix),
    shell:
        """
        for chrom in {params.sorted_chroms}; do
            cat {params.prefix}/*_heldout/"$chrom"/{params.suffix};
        done > {output}
        """


# Generate a prediction file for each chromosome in "short with score" format of Aiden lab's Juicer software.
# https://github.com/aidenlab/juicer/wiki/Pre#short-with-score-format
rule generate_juicer_short:
    input:
        os.path.join(output_dir, "prediction/{method}/%s.bedpe.gz" % pred_prefix),
    output:
        temp(os.path.join(output_dir, "prediction/{method}/%s.sws" % pred_prefix)),
    shell:
        """
        gzip -cd {input} | awk -v OFS="\t" '{{print 0,$1,$2,0,0,$1,$5,1,$11}}' > {output}
        """


# Run Juicer software tool pre to convert the .sws file outputted by generate_juicer_short to a .hic file that can be viewed on the UCSC Genome Browser.
# https://github.com/aidenlab/juicer/wiki/Pre
# https://genome.ucsc.edu/goldenPath/help/hic.html
rule run_juicer_pre:
    input:
        os.path.join(output_dir, "prediction/{method}/%s.sws" % pred_prefix),
    params:
        resolution=window_size,
        seed=seed,
    output:
        os.path.join(output_dir, "prediction/{method}/%s.hic" % pred_prefix),
    shell:
        """
        . /u/local/Modules/default/init/modules.sh; module load java/jre-1.8.0_281
        time java -Xmx6g -jar source/juicer_tools_1.22.01.jar pre \
        -r {params.resolution} -n --random-seed {params.seed} \
        {input} {output} hg38
        """


# Extract Hi-C interaction data from 4DN that matches our resolution and pairwise distances for each chromosome
rule report_hic:
    input:
        path_chrom_size_file,
    params:
        max_dist=max_dist
    output:
        os.path.join(output_dir, "prediction/hic/{chrom}/{norm}.txt.gz"),
    shell:
        """
        time python3 \
        -u /u/project/ernst/skwon94/pairwise_scoring/source/reportHiC.py \
        -f https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/a98ca64a-861a-4a8c-92e9-586af457b1fb/4DNFI1UEG1HD.hic \
        -c {wildcards.chrom} \
        -l {input} \
        -n {wildcards.norm} \
        -o {output} \
        -m {params.max_dist}
        """


rule aggregate_hic:
    input:
        expand(
            os.path.join(output_dir, "prediction/hic/{chrom}/{norm}.txt.gz"),
            chrom=all_chroms,
            allow_missing=True,
        ),
    output:
        os.path.join(output_dir, "prediction/hic/{norm}.txt.gz"),
    threads: 4
    shell:
        """
        cat {input} | gzip -cd | LC_ALL=C sort -k1,1 --parallel=4 | uniq | gzip > {output}
        """
