#!/bin/bash
if [ $# -ne 1 ] ; then
echo "Parameter> "
echo "1. pivot       = pivot file (wihtout .csv)"
exit 1
fi

data="PUBMED23"
pr_dir=src
pv_dir=pivot
ds_dir=ftr
bk_dir=bkt

ulimit -Sn 4000


pivot=$1; shift
from=0
to=23
nt=8

pv=${pv_dir}/${pivot}.csv
if [ ! -e $pv ]; then
  echo pivot file = $pv does not exist.
  exit
fi
wd=$($pr_dir/pivot_property.sh -w $pv)
dm=$($pr_dir/pivot_property.sh -d $pv)
pt=$($pr_dir/pivot_property.sh -p $pv)
np=$($pr_dir/pivot_property.sh -n $pv)
rt=$($pr_dir/pivot_property.sh -r $pv)

prefix=pubmed23_

range=$($pr_dir/expand_fn.sh "${ds_dir}/${prefix}" ${from} ${to} ".ftr" 0)

bk="$bk_dir/${pivot}_${range}.bkt"

ds=$($pr_dir/expand_fn.sh "${ds_dir}/${prefix}" ${from} ${to} ".ftr" 1)

cflags0="-O3 -fopenmp "
cflags0="$cflags0 -DNUM_THREADS=$nt"

cflags0="$cflags0 -DCOMPILE_TIME"
cflags0="$cflags0 -DPUBMED23"

cflags0="$cflags0 -DFTR_DIM=$dm"
cflags0="$cflags0 -DPJT_DIM=$wd"
if [ $wd -gt 64 ] ; then
cflags0="$cflags0 -DEXPANDED_SKETCH"
elif [ $wd -gt 32 ] ; then
cflags0="$cflags0 -DWIDE_SKETCH"
else
cflags0="$cflags0 -DNARROW_SKETCH"
fi
cflags0="$cflags0 -DFTR_ON_SECONDARY_MEMORY"
cflags0="$cflags0 -DFTR_ARRANGEMENT_ASIS"
cflags0="$cflags0 -DSEQUENTIAL_FILTERING"

if [ $pt == 3 ] ; then
  cflags0="$cflags0 -DPARTITION_TYPE_PQBP"
  cflags0="$cflags0 -DNUM_PART=$np"
  if [ $rt == ROTATION ] ; then
  cflags0="$cflags0 -DPRE_ROTATION"
  fi
elif [ $pt == 4 ] ; then
  cflags0="$cflags0 -DPARTITION_TYPE_SQBP"
elif [ $pt == 5 ] ; then
  cflags0="$cflags0 -DPARTITION_TYPE_SQBP"
else
  cflags0="$cflags0 -DPARTITION_TYPE_QBP"
fi


cflags0="$cflags0 -DBLOCK_SIZE=200"
cflags="$cflags -DPIVOT_FILE=\"$pv\""
cflags="$cflags -DBUCKET_FILE=\"$bk\""

echo $cflags0
echo $cflags

lb_list="bit_op ftr e_time sketch quick"

lb=""
for l in $lb_list ; do
  lb="$lb $pr_dir/$l.c"
done
echo library = $lb

pr=make_bucket

set -x

gcc $cflags0 $cflags $pr_dir/$pr.c $lb -o $pr -lm

if [ $? == 1 ] ; then
exit 1
fi

time ./$pr $ds

