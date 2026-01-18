#!/bin/bash
if [ $# -ne 8 ] ; then
echo "Usage> set_app.sh src ftr qry piv bkt smp bch dataset"
echo "1.  src = source  "
echo "2.  ftr = dataset "
echo "3.  qry = query   " 
echo "4.  piv = pivot   " 
echo "5.  bkt = bucket  " 
echo "6.  smp = smap    " 
echo "6.  bch = batch   " 
echo "8.  dataset name  " 
exit
fi

src=$1; shift
ftr=$1; shift
qry=$1; shift
piv=$1; shift
bkt=$1; shift
smp=$1; shift
bch=$1; shift
ds=$1; shift

rm /app/src
rm /app/ftr
rm /app/query
rm /app/pivot
rm /app/bkt
rm /app/batch
rm /app/smap

ln -s $src /app/src
ln -s $ftr /app/ftr
ln -s $qry /app/query
ln -s $piv /app/pivot
ln -s $bkt /app/bkt
ln -s $bch /app/batch
ln -s $smp /app/smap

export DATASET=$ds
