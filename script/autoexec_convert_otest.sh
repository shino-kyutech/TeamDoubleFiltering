#!/bin/bash
# 入力:		PUBMED23 の全体のh5ファイル"benchmark-dev-pubmed23.h5"（2300万個）
#			otest/queries (11000個)
# 出力:		11個の部分質問（各1000），query_20.ftr, ... , query_30.ftr
# Data_id：	40000000 から （オフセット = 1）

pr=convert_h5_to_ftr_PUBMED23_v2

h5cc -O3 -Wall src/$pr.c -o $pr

if [ $? == 1 ] ; then
echo compile error
exit 1
fi

h5="data/benchmark-dev-pubmed23.h5"
ds=otest/queries

for n in {20..30} ; do
	if [ $n -lt 10 ] ; then
		ftr=query/query_0$n.ftr
	else
		ftr=query/query_$n.ftr
	fi
	num=1000
	nn=$(($n - 20))
	from=$(($nn * ${num}))
	offset=40000001
	factor=500
	echo ./$pr "${h5}" ${ds} ${from} ${num} ${offset} ${factor} ${ftr}
	./$pr "${h5}" ${ds} ${from} ${num} ${offset} ${factor} ${ftr}
done

from=0
num=10000
ftr=query/query_otest.ftr
./$pr "${h5}" ${ds} ${from} ${num} ${offset} ${factor} ${ftr}

