#!/bin/bash
if [ $# -ne 12 ] ; then
echo "Usage> ./search.sh parms:"
echo "1.  width      = width of pivot for 1st filtering"
echo "2.  smap dim.  = dimension of smap"
echo "3.  rs_dir     = directory for result files"
echo "4.  log_file   = search log of each query"
echo "5.  q_bit      = quantization bit of QPSMAP" 
echo "6.  q_range    = quantize range"
echo "7.  plus_half  = plus in recover (PLUS_HALF), or -1 -> use lower bound"
echo "8.  FTR_ON     = ftr data on: 0 -> Second memory, 1 -> RAM, 2 -> 2nd on SSD"
echo "9.  batch_file = parameter file (nc1 and nc2)"
echo "10. program version"
echo "11. #threads   = number of threads"
echo "12. #queries   = number of queries"

exit 1
fi

script=/app/script

width=$1; shift
smap=$1; shift
rs_dir=$1; shift
log=$1; shift
qbit=$1; shift
qrange=$1; shift
phalf=$1; shift
ftr=$1; shift
bt=$1; shift
version=$1; shift
nt=$1; shift
queries=$1; shift

if [ $bt == "NONE" ] ; then
batch=NONE
else
batch=batch/$bt.txt
fi

if [ $ftr -eq 0 ] ; then
summary=summary_w${width}_s${smap}_nc2_${bt}_${qbit}bit_${nt}t_ssd.csv
else
summary=summary_w${width}_s${smap}_nc2_${bt}_${qbit}bit_${nt}t_ram.csv
fi

echo $summary

if [ -z $DATASET ] ; then
echo "$DATASET" is not defined!
exit 1
elif [ $DATASET == "DECAF" ] ; then
range=00_96
query=00_04
nn=1
p2=20
use_pd=0
elif [ $DATASET == "PUBMED23" ] ; then
range=00_23
#query=otest
query=20_29
nn=30
p2=25
use_pd=1
elif [ $DATASET == "LAION2B" ] || [ $DATASET == "LAION100M" ] ; then
range=00_102
query=00_09
nn=1
p2=20
use_pd=0
elif [ $DATASET == "DEEP1B" ] ; then
range=00_99
query=00_09
nn=1
p2=20
use_pd=0
fi

echo RANGE = $range, QUERY = $query

if [ $version == "wsl" ] ; then
script=$SCR
fi
echo $script/search_by_double_filtering_v7.sh sketch_w$width QPSMAP_p$smap $range $query NONE NONE $rs_dir $nn 0 $qbit $qrange $phalf 10 $p2 $nt 7 13 3 1.05 100 0 8 $ftr $use_pd $batch $summary $log $version

if [ $version == "v7" ] ; then
$script/search_by_double_filtering_v7.sh sketch_w$width QPSMAP_p$smap $range $query NONE NONE $rs_dir $nn $queries $qbit $qrange $phalf 10 $p2 $nt 7 13 3 1.05 100 0 8 $ftr $use_pd $batch $summary $log $version
elif [ $version == "v8" ] ; then
$script/search_by_double_filtering_v8.sh sketch_w$width QPSMAP_p$smap $range $query NONE NONE $rs_dir $nn $queries $qbit $qrange $phalf 10 $p2 $nt 7 13 3 1.05 100 0 8 $ftr $use_pd $batch $summary $log $version
elif [ $version == "wsl" ] ; then
$script/search_by_double_filtering_wsl.sh sketch_w$width QPSMAP_p$smap $range $query NONE NONE $rs_dir $nn $queries $qbit $qrange $phalf 10 $p2 $nt 7 13 3 1.05 100 0 8 $ftr $use_pd $batch $summary $log $version
else
$script/search_by_double_filtering.sh sketch_w$width QPSMAP_p$smap $range $query NONE NONE $rs_dir $nn $queries $qbit $qrange $phalf 10 $p2 $nt 7 13 3 1.05 100 0 8 $ftr $use_pd $batch $summary $log $version
fi
