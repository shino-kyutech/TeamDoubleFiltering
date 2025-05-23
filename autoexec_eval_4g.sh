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

echo ./search.sh 16 384 $rs_dir ${log}_w16 1 0 batch_w16_384bit $version $nt

# search by double filtering
if [ $log == NONE ] ; then
./search.sh 16 384 $rs_dir NONE 1 0 batch_w16_384bit $version $nt
./search.sh 18 384 $rs_dir NONE 1 0 batch_w18_384bit $version $nt
./search.sh 20 384 $rs_dir NONE 1 0 batch_w20_384bit $version $nt
else
./search.sh 16 384 $rs_dir ${log}_w16 1 0 batch_w16_384bit $version $nt
./search.sh 18 384 $rs_dir ${log}_w18 1 0 batch_w18_384bit $version $nt
./search.sh 20 384 $rs_dir ${log}_w20 1 0 batch_w20_384bit $version $nt
fi
