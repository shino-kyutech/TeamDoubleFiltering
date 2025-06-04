# TeamDoubleFiltering

This repository contains the code and Docker environment for participating in the **SISAP 2025 Indexing Challenge (Task 1)**.

## Overview

We propose a **double filtering method** for efficient nearest neighbor search on large-scale vector datasets.  
Our method uses two stages of quantized projections (short and long) to reduce the number of candidate vectors while maintaining high accuracy.

## Repository Structure

```
.
├── Dockerfile
├── autoexec_all.sh               # prepare system and run search
├── autoexec_convert_dataset.sh   # convert vectors of pubmed23 dataset to char vactors
├── autoexec_convert_otest.sh     # convert query vectors to char vectors
├── autoexec_eval.sh              # run search for evaluation using 16 hyper-parameters
├── autoexec_eval_4g.sh           # run search for RAM = 4G environment
├── autoexec_gt_otest_answer.sh   # convert ground truth of otest to answer format (csv)
├── autoexec_make_bucket.sh       # make bucket file for sketch
├── autoexec_prepare.sh           # prepare system (vector conversion, bucket construction, ... )
├── search.sh                     # search using hyper parameters
├── search_by_double_filtering.sh # search using double filtering
├── batch/                        # hyperparameters (nc1 and nc2) for search
├── pivot/                        # pivot files (.csv) for sketch and QSMAP
├── src/                          # Source code for indexing and searching
```

## How to Build the Docker Image

Run the following command in the root directory of the repository:

```
docker build -t doublefiltering .
```

## How to Run the Search

Run the container with the dataset and output directory mounted as volumes:

```
docker run --rm --memory=16g --memory-swap=16g -v <path to benchmark-dev-pubmed23.h5>:/app/data/benchmark-dev-pubmed23.h5:ro   -v <path to directory for result files>:/app/result doublefiltering /bin/bash /app/autoexec_all.sh
```

Replace the paths with your actual dataset path and result output directory.

## How to Run the Search with RAM = 4GB 

Run the container in iterative mode with the dataset and output directory mounted as volumes:

```
docker run --rm --memory=4g --memory-swap=4g -it -v <path to benchmark-dev-pubmed23.h5>:/app/data/benchmark-dev-pubmed23.h5:ro   -v <path to directory for result files>:/app/result doublefiltering /bin/bash
```

After entering container, first run autoexec_prepare.sh and then run autoexec_eval_4g.sh. The results will be stored in directory "temp".

> Note: The script `autoexec_eval_4g.sh` is designed for experiments under strict memory constraints (e.g., 4GB RAM), as discussed in the accompanying short paper.  
> Be sure to run it inside a Docker container launched with the `--memory=4g` and `--memory-swap=4g` options, as shown above, to ensure the memory limit is properly enforced.  
> This evaluation illustrates the effectiveness of our method even under low-resource environments.

## Output Files

The search results will be saved in the `/app/result/` directory, including:

- `result.csv`, `result1.csv`, ..., each line formatted as:

  ```
  query-ID, distNN, answer-ID[1], dist[1], answer-ID[2], dist[2], ...
  ```

  where `query-ID` and `answer-ID` are 0-based indices.

- `summary_all.csv`: a summary of the overall search run, with rows like:

  ```
  trial, file, width, q_bit, ftr_on, nc1, nc2, recall@1, recall@30,
  filtering, 1st, 2nd, k_NN, ave(ms/q), stdev(ms/q), min(ms/q), max(ms/q)
  ```

### Explanation of Metrics

- `nc1`, `nc2`: number of candidates obtained by the 1st and 2nd filtering stages  
- `filtering`, `1st`, `2nd`: filtering costs in seconds  
- `k_NN`: cost for final reranking among filtered candidates  
- `ave`, `stdev`, `min`, `max`: per-query search time statistics (milliseconds)

## Requirements

- Docker (tested with Docker 24.0+)
- Host machine with at least **16GB RAM**

## License

This repository is licensed under the **MIT License**.  
See [LICENSE](./LICENSE) for details.

## Contact

If you have any questions, please contact us at:

- GitHub: [shino-kyutech](https://github.com/shino-kyutech/TeamDoubleFiltering)  
- Email: [shino.kyutech@gmail.com]
