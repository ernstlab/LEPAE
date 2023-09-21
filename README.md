## Pairwise scoring of human genomic windows based on epigenomic and TF binding annotations

### Genome-wide LEPAE score
Score trakcs are provided in [`.hic`](https://genome.ucsc.edu/goldenPath/help/hic.html) format.
- [1-kb resolution](https://public.hoffman2.idre.ucla.edu/ernst/R0RG6/LEPAE/prediction_min1000_max100000_window1000_070422.hic)
- [10-kb resolution](https://public.hoffman2.idre.ucla.edu/ernst/R0RG6/LEPAE/prediction_min10000_max1000000_window10000_070422.hic)

### Running the pipeline
  It is highly recommended to run the pipeline using high-performance computing (HPC) as follows:
  ```
  snakemake -j 100 --cluster-config cluster.json --cluster \
  "qsub -V -l h_rt={cluster.time},h_data={cluster.memory} -pe shared {cluster.threads} -N {cluster.name} -j y -o {cluster.output}" \
  --latency-wait 60 --quiet --local-cores 1
  ```
