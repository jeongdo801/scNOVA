#!/bin/bash

snakemake -j 10 --latency-wait 60 --restart-times 1 --use-conda --rerun-incomplete \
--cluster-config Snake.cluster.json --cluster "{cluster.sbatch} -p {cluster.partition} --cpus-per-task {cluster.n} --time {cluster.time} --mem {cluster.mem}" \ 
