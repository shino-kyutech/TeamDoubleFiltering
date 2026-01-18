#!/bin/bash

qpsmap=$1; shift
pivot2=$1; shift
echo $qpsmap
echo $pivot2
dataset=$DATASET

if [ $dataset == "PUBMED23" ] ; then
	qr_prefix=query_
	prefix=pubmed23_
elif [ $dataset == "LAION2B" ] ; then
	qr_prefix=laion2b_query_
	prefix=laion2b_
elif [ $dataset == "DEEP1B" ] ; then
	qr_prefix=query_
	prefix=base_
fi

pr_dir=/app/src
ds_dir=/app/ftr
qr_dir=/app/query
pv_dir=/app/pivot
bk_dir=/app/bkt
sm_dir=/app/smap

p2="$pv_dir/${pivot2}.csv"
if [ ! -e $p2 ]; then
  echo pivot file of qpsmap for 2nd filtering = $p2 does not exist.
  exit
fi

w2=$($pr_dir/pivot_property.sh -w $p2)
echo SMAP_DIM = $w2
d2=$($pr_dir/pivot_property.sh -d $p2)
pt=$($pr_dir/pivot_property.sh -p $p2)
np=$($pr_dir/pivot_property.sh -n $p2)
echo partition type of 2nd pivot = $pt, number of partitioned spaces = $np

sm="$sm_dir/${qpsmap}.sm"

echo sm = $sm

smap_dim=$($pr_dir/qpsmap_header $sm DIM)

echo smap_dim = $smap_dim

if [ $w2 -eq $smap_dim ] ; then
    echo "SMAP_DIM OK"
else
    echo "invalid SMAP_DIM"
fi
