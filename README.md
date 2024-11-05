# Pairwise scoring of human genomic windows based on epigenomic and TF binding annotations

## Genome-wide LEPAE score for hg38
Score tracks are provided in [`.hic`](https://genome.ucsc.edu/goldenPath/help/hic.html) format.
- [1-kb resolution](https://public.hoffman2.idre.ucla.edu/ernst/R0RG6/LEPAE/prediction_min1000_max100000_window1000_070422.hic)
- [10-kb resolution](https://public.hoffman2.idre.ucla.edu/ernst/R0RG6/LEPAE/prediction_min10000_max1000000_window10000_070422.hic)

## Software

### Get started
This pipeline is written in Snakemake, which enables the pipeline to run on various execution environments. It is highly recommended to run the pipeline in a high-performance computing (HPC) environment. We used UCLA's Hoffman2 cluster. 

Some minor modifications may be required to run in a different environment. If the following dependencies are available, there is no installation step other than cloning this repository.

#### Software dependencies
- Python 3.9.6
  - Numpy 1.22.2
  - Pandas 1.3.2
  - Pytorch 1.10.0+cu102
  - Scikit-learn 0.24.2
  - StrawC 0.0.9
- Bedtools 2.30.0
- [Juicer](https://github.com/aidenlab/juicer/wiki/Pre) 1.22.01

### Demo --TODO
- Instructions to run on data
- Expected output
- Expected run time for demo on a "normal" desktop computer

### Instruction for use 
How to run the software on your data --TODO: add more details
  ```
  snakemake -j 100 --cluster-config cluster.json --cluster \
  "qsub -V -l h_rt={cluster.time},h_data={cluster.memory} -pe shared {cluster.threads} -N {cluster.name} -j y -o {cluster.output}" \
  --latency-wait 60 --quiet --local-cores 1
  ```
