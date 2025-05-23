#!/bin/bash
if [ $# -ne 4 ] ; then
echo "Usage> ./autoexec_eval.sh parms:"
echo "1.  rs_dir     = directory for result files"
echo "2.  log_file   = search log of each query"
echo "3.  program version"
echo "4.  #threads   = number of threads"

exit 1
fi

rs_dir=$1; shift
log=$1; shift
version=$1; shift
nt=$1; shift

# search by double filtering
if [ $log == NONE ] ; then
#./search.sh 16 192 $rs_dir NONE 1 1 batch_16G_w16_1bit_RAM $version $nt
./search.sh 18 192 $rs_dir NONE 1 1 batch_16G_w18_1bit_RAM $version $nt
./search.sh 20 192 $rs_dir NONE 1 1 batch_16G_w20_1bit_RAM $version $nt
else
#./search.sh 16 192 $rs_dir ${log}_w16 1 1 batch_16G_w16_1bit_RAM $version $nt
./search.sh 18 192 $rs_dir ${log}_w18 1 1 batch_16G_w18_1bit_RAM $version $nt
./search.sh 20 192 $rs_dir ${log}_w20 1 1 batch_16G_w20_1bit_RAM $version $nt
fi
