## Pairwise scoring of human genomic windows based on epigenomic and TF binding annotations

### Example of running the pipeline on Hoffman2 cluster
Get an interactive node (```qrsh -l h_rt=8:00:00```) and run: 
  ```
  cd /u/project/ernst/skwon94/pairwise_scoring/
  snakemake -j 100 --cluster-config cluster.json --cluster \
  "qsub -V -l h_rt={cluster.time},h_data={cluster.memory} -pe shared {cluster.threads} -N {cluster.name} -j y -o {cluster.output}" \
  --latency-wait 60 --quiet --local-cores 1
  ```