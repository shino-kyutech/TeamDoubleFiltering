#!/bin/bash
# 入力:		PUBMED23 の全体のh5ファイル"benchmark-dev-pubmed23.h5"（2300万個）
#			otest/queries (11000個)
# 出力:		11個の部分質問（各1000），query_20.ftr, ... , query_30.ftr
# Data_id：	40000000 から （オフセット = 1）
# コンパイル：h5cc -O3 ground_trueth_to_csv_PUBMED_v2.c -o ground_trueth_to_csv_PUBMED_v2
# 実行時パラメータ: file.h5 [itest | otest] num_queries k > answer_[i | o].csv
#   1. file.h5       = input file of h5 format
#   2. itest | otest 
#   3. num_queries
#   4. k             = number of NN

h5cc -O3 src/ground_trueth_to_csv_PUBMED_v2.c -o ground_trueth_to_csv_PUBMED_v2

if [ $? == 1 ] ; then
echo compile error
exit 1
fi

h5="data/benchmark-dev-pubmed23.h5"
ds=otest

./ground_trueth_to_csv_PUBMED_v2 $h5 $ds 10000 100 > query/answer_q_otest_r_00_23.csv

