#!/bin/bash

snakemake -j 4 --latency-wait 60 --restart-times 1
