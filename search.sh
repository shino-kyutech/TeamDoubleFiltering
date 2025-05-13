#!/bin/bash
if [ $# -ne 5 ] ; then
echo "Usage>"
echo "1.  width      = width of pivot for 1st filtering"
echo "2.  rs_dir     = directory for result files"
echo "3.  q_bit      = quantization bit of QPSMAP" 
echo "4.  FTR_ON     = ftr data on: 0 -> Second memory, 1 -> RAM"
echo "5.  batch_file = parameter file (nc1 and nc2)"

exit 1
fi

width=$1; shift
rs_dir=$1; shift
qbit=$1; shift
ftr=$1; shift
bt=$1; shift

if [ $bt == "NONE" ] ; then
batch=NONE
else
batch=batch/$bt.txt
fi

if [ $ftr -eq 0 ] ; then
summary=summary_w${width}_nc2_${bt}_${qbit}bit_ssd.csv
else
summary=summary_w${width}_nc2_${bt}_${qbit}bit_ram.csv
fi

echo $summary

./search_by_double_filtering.sh sketch_w$width QPSMAP_p192 00_23 otest NONE NONE $rs_dir 30 0 $qbit 2.0 0.5 10 25 8 7 13 3 1.05 100 0 8 $ftr 1 $batch $summary NONE
