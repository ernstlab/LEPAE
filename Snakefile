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
num_tuning_chrom = 3
training_chrom = dict()
tuning_chrom = dict()
for heldout_chrom in all_chroms:
    sampled_chrom = list(
        random.sample(
            [c for c in all_chroms if c != heldout_chrom and c != "chrX"],
            num_tuning_chrom,
        )
    )
    training_chrom[heldout_chrom] = [
        c
        for c in all_chroms
        if c != heldout_chrom and c != "chrX" and c not in sampled_chrom
    ]
    tuning_chrom[heldout_chrom] = sampled_chrom

# Random forest training parameters
training_data_size = config["training_data_size"]
tuning_data_size = config["tuning_data_size"]
num_trees = config["num_trees"]
random_search_n = config["random_search_n"]

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
pred_prefix = "prediction_min%d_max%d_window%d" % (
    min_dist,
    max_dist,
    window_size,
)

# Whether to use fraction of feature coverage as feature value (instead of binary)
frac_feature = config["frac_feature"]


### Pipeline rules
# Rules to run locally and not as jobs to the cluster
localrules:
    all,
    clean,
    prepare_input_to_generate_data,


# Target rule
rule all:
    input:
        os.path.join(output_dir, "prediction/%s.bedpe.gz" % pred_prefix),


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
        bedtools intersect -sorted -wa -a {input.window_file} -b - | cut -f 4 | sort -V | uniq | gzip > {output}
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
        bedtools map -a {input.window_file} -b - -c 4 -o distinct | cut -f 4- | awk '$NF!="."' | gzip > {output}
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
        frac_feature="-f" if frac_feature else " ",
        window_size=window_size,
    threads: 2
    output:
        os.path.join(output_dir, "data/{chrom}.window"),
    shell:
        """
        python3 source/generateDataThreaded.py \
        -p {input.window_file} -ch {input.chromhmm_list_file} -dn {input.dnasechipseq_list_file} \
        -chn {params.chromhmm_num_states} -n {params.num_features} -c {wildcards.chrom} \
        {params.frac_feature} -w {params.window_size} \
        -o {output}
        """


# Random search for random forest classifier hyperparameters.
# Output file ending with .progress.txt reports each combination of hyperparameters and its classification performance metric (MSE, AUROC, AUPRC).
# Output file ending with .best_hyperparam.txt stores the combination of hyperparameters that had the highest tuning AUROC computed from held-out tuning data.
# For each positive or negative pair, a pair with its upstream and downstream windows flipped is also provided as an additional pair, doubling the number of training and tuning data overall.
rule hyperparam_search:
    input:
        expand(os.path.join(output_dir, "data/{chrom}.window"), chrom=all_chroms),
    params:
        random_search_n=random_search_n,
        num_features=num_features,
        training_data_size=training_data_size,
        tuning_data_size=tuning_data_size,
        num_trees=num_trees,
        window_size=window_size,
        output_filename_prefix=os.path.join(
            output_dir, "classifier/{heldout_chrom}/hyperparam_search_dist"
        ),
        seed=seed,
        frac_feature="-f" if frac_feature else " ",
        training_data_filenames=lambda wildcards: expand(
            os.path.join(output_dir, "data/{chrom}.window"),
            chrom=training_chrom[wildcards.heldout_chrom],
        ),
        tuning_data_filenames=lambda wildcards: expand(
            os.path.join(output_dir, "data/{chrom}.window"),
            chrom=tuning_chrom[wildcards.heldout_chrom],
        ),
    output:
        os.path.join(
            output_dir,
            "classifier/{heldout_chrom}/hyperparam_search_dist{dist}.best_hyperparam.txt",
        ),
        os.path.join(
            output_dir,
            "classifier/{heldout_chrom}/hyperparam_search_dist{dist}.progress.txt",
        ),
    threads: 4
    shell:
        """
        python3 -u source/train.py \
        --training-data-filenames {params.training_data_filenames} \
        --tuning-data-filenames {params.tuning_data_filenames} \
        --positive-training-data-size {params.training_data_size} \
        --positive-tuning-data-size {params.tuning_data_size} \
        --pairwise-dist {wildcards.dist} \
        --window-size {params.window_size} \
        --seed {params.seed} \
        --num-trees {params.num_trees} \
        --random-search --random-search-n {params.random_search_n} \
        --num-features {params.num_features} {params.frac_feature} \
        --output-filename-prefix {params.output_filename_prefix}{wildcards.dist}
        """


# Train and save the final random forest classifier with the best combination of hyperparameters.
# As done in hyperparam_search, flipped pairs are added to training and tuning data, doubling the data size.
rule train:
    input:
        best_hyperparam_file=os.path.join(
            output_dir,
            "classifier/{heldout_chrom}/hyperparam_search_dist{dist}.best_hyperparam.txt",
        ),
        all_data_filenames=expand(
            os.path.join(output_dir, "data/{chrom}.window"), chrom=all_chroms
        ),
    params:
        random_search_n=random_search_n,
        num_features=num_features,
        training_data_size=training_data_size,
        tuning_data_size=tuning_data_size,
        num_trees=num_trees,
        window_size=window_size,
        output_filename_prefix=os.path.join(
            output_dir, "classifier/{heldout_chrom}/train_dist"
        ),
        seed=seed,
        frac_feature="-f" if frac_feature else " ",
        training_data_filenames=lambda wildcards: expand(
            os.path.join(output_dir, "data/{chrom}.window"),
            chrom=training_chrom[wildcards.heldout_chrom],
        ),
        tuning_data_filenames=lambda wildcards: expand(
            os.path.join(output_dir, "data/{chrom}.window"),
            chrom=tuning_chrom[wildcards.heldout_chrom],
        ),
        test_data_filenames=lambda wildcards: expand(
            os.path.join(output_dir, "data/{chrom}.window"),
            chrom=wildcards.heldout_chrom,
        ),
    output:
        os.path.join(
            output_dir, "classifier/{heldout_chrom}/train_dist{dist}.progress.txt"
        ),
        os.path.join(output_dir, "classifier/{heldout_chrom}/train_dist{dist}.pkl"),
        os.path.join(
            output_dir,
            "classifier/{heldout_chrom}/train_dist{dist}.training_pred.txt.gz",
        ),
        os.path.join(
            output_dir, "classifier/{heldout_chrom}/train_dist{dist}.tuning_pred.txt.gz"
        ),
        os.path.join(
            output_dir, "classifier/{heldout_chrom}/train_dist{dist}.test_pred.txt.gz"
        ),
    threads: 4
    shell:
        """
        python3 -u source/train.py \
        --training-data-filenames {params.training_data_filenames} \
        --tuning-data-filenames {params.tuning_data_filenames} \
        --test-data-filenames {params.test_data_filenames} \
        --positive-training-data-size {params.training_data_size} \
        --positive-tuning-data-size {params.tuning_data_size} \
        --pairwise-dist {wildcards.dist} \
        --window-size {params.window_size} \
        --seed {params.seed} \
        --num-trees {params.num_trees} \
        --num-features {params.num_features} {params.frac_feature} \
        --save \
        --output-filename-prefix {params.output_filename_prefix}{wildcards.dist} \
        @{input.best_hyperparam_file}
        """


# Make predictions for the chromosome that was held out specifically for prediction.
rule predict:
    input:
        classifier_file=os.path.join(
            output_dir, "classifier/{chrom}/train_dist{dist}.pkl"
        ),
        input_filename=os.path.join(output_dir, "data/{chrom}.window"),
    params:
        num_features=num_features,
        window_size=window_size,
        seed=seed,
        frac_feature="-f" if frac_feature else " ",
    output:
        os.path.join(
            output_dir, "prediction/{chrom}/%s_dist{dist}.txt.gz" % pred_prefix
        ),
    threads: 4
    shell:
        """
        python3 source/predict.py \
        -t {input.classifier_file} -i {input.input_filename} \
        -d {wildcards.dist} -w {params.window_size} \
        -n {params.num_features} {params.frac_feature} \
        -s {params.seed} -o {output}
        """


# Generate a prediction file in "short with score" format of Aiden lab's Juicer software.
# https://github.com/aidenlab/juicer/wiki/Pre#short-with-score-format
rule generate_juicer_short:
    input:
        expand(
            os.path.join(
                output_dir,
                "prediction/{chrom}/%s_dist{dist}.txt.gz" % pred_prefix,
            ),
            chrom=all_chroms,
            dist=all_dists,
        ),
    params:
        input_dirs=expand(
            os.path.join(output_dir, "prediction/{chrom}"), chrom=all_chroms
        ),
        input_prefix=pred_prefix + "_dist",
        input_suffix=".txt.gz",
    output:
        os.path.join(output_dir, "prediction/%s.sws" % pred_prefix),
    shell:
        """
        for d in `\ls -d {params.input_dirs} | sort -k1,1`; do
            cat $d/{params.input_prefix}*{params.input_suffix} | gzip -cd | awk -v OFS="\t" '{{print 0,$1,$2,0,0,$1,$6,1,$9}}' | sort -k3,3n -k7,7n;
        done > {output}
        """


# Generate a prediction file in BEDPE format.
# https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format
rule generate_bed:
    input:
        os.path.join(output_dir, "prediction/%s.sws" % pred_prefix),
    params:
        window_size=window_size,
    output:
        os.path.join(output_dir, "prediction/%s.bedpe.gz" % pred_prefix),
    shell:
        """
        awk -v OFS="\t" '{{ print $2, $3, $3 + {params.window_size}, $2, $7, $7 + {params.window_size}, NR, $NF, ".", ".", {params.window_size}, $7-$3 }}' {input} | gzip > {output}
        """
