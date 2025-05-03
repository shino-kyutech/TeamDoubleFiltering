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
/bin/bash /app/autoexec_all_pubmed23.sh
```
## Result Files

The result of search is recorded in files in the folder result: 
result.csv, result1.csv, ... , consisting of lines of the form:

query-ID, distNN, answer-ID[1], dist[1], ... .

