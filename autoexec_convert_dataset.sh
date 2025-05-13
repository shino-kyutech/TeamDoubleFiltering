#!/bin/bash
# 入力:		PUBMED23 のメインのデータセット"benchmark-dev-pubmed23.h5"（2300万個）
# 出力:		24個の部分データセット（各100万個），pubmed23_00.ftr, ... , pubmed23_23.ftr
# Data_id：	0 から 99999999（オフセット = 1）

pr=convert_h5_to_ftr_PUBMED23_v2

h5cc -O3 -Wall src/$pr.c -o $pr

if [ $? == 1 ] ; then
echo compile error
exit 1
fi

h5="data/benchmark-dev-pubmed23.h5"
ds=train
ds_dir=ftr

for n in {0..23} ; do
	if [ $n -lt 10 ] ; then
		ftr=$ds_dir/pubmed23_0$n.ftr
	else
		ftr=$ds_dir/pubmed23_$n.ftr
	fi
	num=1000000
	from=$(($n * ${num}))
	offset=1
	factor=500
echo ./$pr "${h5}" ${ds} ${from} ${num} ${offset} ${factor} ${ftr}
     ./$pr "${h5}" ${ds} ${from} ${num} ${offset} ${factor} ${ftr}
done

