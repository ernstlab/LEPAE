# Learning Evidence of Pairwise Association from Epigenomic and TF binding data

**L**earning **E**vidence of **P**airwise **A**ssociation from **E**pigenomic and TF binding data (LEPAE) quantifies evidence of association for pairs of genomic windows based on large-scale epigenomic and TF binding data while considering distance information.

## Genome-wide LEPAE score for hg38

Score tracks are provided in [`.hic`](https://genome.ucsc.edu/goldenPath/help/hic.html) format.

- [1-kb resolution](https://public.hoffman2.idre.ucla.edu/ernst/R0RG6/LEPAE/prediction_min1000_max100000_window1000_070422.hic)
- [10-kb resolution](https://public.hoffman2.idre.ucla.edu/ernst/R0RG6/LEPAE/prediction_min10000_max1000000_window10000_070422.hic)

## Pipeline

This pipeline is written in Snakemake, which enables running on various execution environments. It is highly recommended to run it in a high-performance computing (HPC) environment. We used UCLA's Hoffman2 cluster.

If the following dependencies are available, there is no installation step other than cloning this repository.

### Pre-requisites

- Snakemake 6.8.1
- Python 3.9.6
  - Numpy 1.22.2
  - Pandas 1.3.2
  - Pytorch 1.10.0+cu102
  - Scikit-learn 0.24.2
- Bedtools 2.30.0

### How to run on full dataset on HPC

```
snakemake \
-j 100 \
--configfile parameters_1kb.yaml \
--cluster-config cluster.json \
--cluster "qsub -V -l h_rt={cluster.time},h_data={cluster.memory} -pe shared {cluster.threads} -N {cluster.name} -j y -o {cluster.output}" \
--latency-wait 60 \
--local-cores 1
```

- `-j 100` is to limit the number of parallel jobs to 100. Adjust this number based on resource availability.
- `--configfile parameters_1kb.yaml` specifies the model configuration file
- `--cluster-config cluster.json` specifies the cluster configuration file
- `--cluster "qsub -V -l h_rt={cluster.time},h_data={cluster.memory} -pe shared {cluster.threads} -N {cluster.name} -j y -o {cluster.output}"` specifies the cluster submission command
- `--latency-wait 60` specifies to wait for 60 seconds before checking the status of the submitted jobs. This is optional. Remove or adjust this number accordingly.
- `--local-cores 1` limits the number of local cores used to 1. This is optional

### Required input

| Name                            | Description                                                                                               | Example(s)                                                                                                                 | Where to specify                   |
| ------------------------------- | --------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------- | ---------------------------------- |
| Snakefile                       | Pipeline definition                                                                                       | [Snakefile](Snakefile)                                                                                                     | NA                                 |
| Config YAML                     | Model configuration                                                                                       | [demo/demo.yaml](demo/demo.yaml), [parameters_1kb.yaml](parameters_1kb.yaml), [parameters_10kb.yaml](parameters_10kb.yaml) | `--configfile` or inside Snakefile |
| Cluster JSON                    | Cluster configuration                                                                                     | [cluster.json](cluster.json)                                                                                               | `--cluster-config`                 |
| Chromosome size                 | Chromosome sizes                                                                                          | [hg38.chrom.sizes](hg38.chrom.sizes)                                                                                       | Config YAML                        |
| DNase-seq and ChIP-seq datasets | Text file listing URLs of gzipped BED files containing peak calls from DNase-seq and ChIP-seq experiments | [files_hg38_DNaseChIPseq.txt](files_hg38_DNaseChIPseq.txt)                                                                 | Config YAML                        |
| ChromHMM annoations             | Text file listing URLs of gzipped BED files containing ChromHMM tracks                                    | [files_hg38_ChromHMM.txt](files_hg38_ChromHMM.txt)                                                                         | Config YAML                        |

The 1-kb resolution LEPAE score was generated using [parameters_1kb.yaml](parameters_1kb.yaml) and the 10-kb resolution LEPAE score was generated using [parameters_10kb.yaml](parameters_10kb.yaml). Other files are used for both resolutions. Input and output paths in config YAML and cluster JSON files should be updated accordingly.

#### Config YAML

```yaml
# parameters_1kb.yaml
env: "cluster" # Or "local" to run in a local environment. Used to load bedtools module (See below)

min_dist: 1000 # Minimum pairwise distance in bp. Typically same as window_size
max_dist: 100000 # Maximum pairwise distance in bp
window_size: 1000 # Window size in bp

seed: 413666 # Random seed for reproducibility

training_data_size: 50000 # Number of training data points
tuning_data_size: 5000 # Number of tuning data points
num_random_search: 10 # Number of combinations of hyperparameters to try during random search
ensemble_size: 5 # Number of models to train in an ensemble

num_dnasechipseq_experiments: 3337 # Number of DNase-seq and ChIP-seq experiments
num_chromhmm_experiments: 127 # Number of ChromHMM datasets
num_chromhmm_states: 25 # Number of ChromHMM states

paths: # Paths to input files and output directory
  dnasechipseq_experiment_paths_file: "/u/project/ernst/skwon94/pairwise_scoring/files_hg38_DNaseChIPseq.txt"
  chromhmm_experiment_paths_file: "/u/project/ernst/skwon94/pairwise_scoring/files_hg38_ChromHMM.txt"
  chrom_size_file: "/u/project/ernst/skwon94/pairwise_scoring/hg38.chrom.sizes"
  output_dir: "/u/project/ernst/skwon94/pairwise_scoring/model/gw_1kb"
```

#### Loading `bedtools` module

If `env` is set to `cluster` in the config YAML file, the `bedtools` module, which is available in our HPC environment, is loaded before running any `bedtools` command as shown below. Depending on your environment, this may not be necessary or need to change (e.g. use a conda environment or Docker image instead).

```bash
# Snakefile
...
if [ "{env}" = "cluster" ]; then
    . /u/local/Modules/default/init/modules.sh; module load bedtools
fi
bedtools ...
```

### Demo

This demo is for a local environment and has been tested in a Linux environment using a Windows laptop.

Create a conda environment with dependencies and run the demo inside the environment as follows:

```
conda create -c conda-forge -c bioconda -n lepae snakemake=6.8.1 python=3.9.6 numpy=1.22.2 pandas=1.3.2 pytorch=1.10.0 scikit-learn=0.24.2 bedtools=2.30.0 tabulate=0.8
conda activate lepae
snakemake --configfile demo/demo.yaml -c1
```

In total, 149 steps run in ~15 minutes. Output file `demo/model/prediction/NN/prediction_min10000_max30000_window10000.bedpe.gz` contains the final score predictions for pairwise distances 10, 20, and 30 kb at a 10-kb resolution.

Given the minimal set of input datasets and small training data size, the trained models are not expected to produce meaningful predictions. To take advantage of this pipeline, scale up to use thousands of datasets and larger training data size as done in [parameters_1kb.yaml](parameters_1kb.yaml) and [parameters_10kb.yaml](parameters_10kb.yaml) in a HPC environment.
