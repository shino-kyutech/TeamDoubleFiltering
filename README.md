# TeamDoubleFiltering

This repository contains the code and Docker environment for participating in the SISAP 2025 Indexing Challenge (Task 1).

## Overview

We propose a double filtering method for efficient nearest neighbor search on large-scale vector datasets.

## How to Build

```bash
docker build -t teamdoublefiltering .
```
## How to Run
```bash
docker run --memory=16g --memory-swap=16g \
-v <path to benchmark-dev-pubmed23.h5>:/app/data/benchmark-dev-pubmed23.h5:ro \
-v <path to directory for result files>:/app/result \
teamdoublefiltering \
/bin/bash /app/autoexec_all.sh
```
## Result Files

The result of search is recorded in files in the folder result: 
result.csv, result1.csv, ... , consisting of lines of the form:

query-ID, distNN, answer-ID[1], dist[1], ... ,
where dataID is renumbered starting from 0.

summary_all.csv contains the summary of search of the form:

trial, file, width, q_bit, ftr_on, nc1, nc2, recall@1, recall@30, filteting, 1st, 2nd, k_NN, ave(ms/q), stdev(ms/q), min(ms/q), max(ms/q),
where, nc1 and nc2 = numbef of candidates of 1st and 2nd filtering,
filtering, 1st, 2nd = cost for overall filtering, 1st filtering, 2nd filtering in seconds,
kNN = cost for reranking to select top 30 from candidates obtained by filtering,
ave, stdev, min, max = average, stdard devience, min, and max of search time (ms/query).

