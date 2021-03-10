#!/bin/bash

snakemake -j 10 --latency-wait 60 --restart-times 1 --use-conda
