#!/bin/bash
if [ $# -ne 7 ] ; then
echo "Usage>"
echo "1.  1st pivot  = pivot file for narrow sketch (for 1st filtering)"
echo "2.  2nd pivot  = pivot file for qpsmap (for 2nd filtering)"
echo "3.  range      = range of ftr files" 
echo "4. q_bit      = quantize bit"
echo "5. q_range    = quantize range"
echo "6. use pd     = (0: without, 1: use pd)"
echo "7. #threads   = number of threads"

exit 1
fi

pr=make_qpsmap

dataset=$DATASET
if [ $dataset == "PUBMED23" ] ; then
	prefix=pubmed23_
elif [ $dataset == "LAION2B" ] ; then
	prefix=laion2b_
elif [ $dataset == "DEEP1B" ] ; then
	prefix=base_
elif [ $dataset == "DECAF" ] ; then
	prefix=fc6_
fi

pr_dir=/app/src
ds_dir=/app/ftr
#qr_dir=/app/query
pv_dir=/app/pivot
bk_dir=/app/bkt
sm_dir=/app/smap

#set -x

ulimit -Sn 10000
#ulimit -a

pivot1=$1; shift    # sketch のピボットファイル名（拡張子を除く）バケットファイル名を決めるためだけに使う．
pivot2=$1; shift    # smap のピボットファイル名（拡張子を除く）
range=$1; shift     # QPSMAP を作成する FTR ファイルの範囲（00_23(PUBMED23), 00_99(DEEP1B), ...)
qbit=$1; shift		# QPSMAP の量子化ビット数（1, 2, 3, 4）
qrange=$1; shift	# 量子化の範囲．標準偏差の倍数で指定．±qrange*sigma の範囲を量子化し，それ以外は最大または最小に丸める
use_pd=$1; shift	# PD (Product Decomposition 直積分解（旧：座標分割）を使用するかどうか．PUBMED23では使用（1を指定）LAION2Bでは未使用（0）
nt="$1"; shift		# スレッド数（PUBMED23では，主に8を使用．1のときにはプログラムが未対応の可能性あり）

p1="$pv_dir/${pivot1}.csv"
if [ ! -e $p1 ]; then
  echo pivot file for 1st filtering = $p1 does not exist.
  exit
fi

p2="$pv_dir/${pivot2}.csv"
if [ ! -e $p2 ]; then
  echo pivot file of qpsmap for 2nd filtering = $p2 does not exist.
  exit
fi

b1="$bk_dir/${pivot1}_${range}.bkt"
if [ ! -e $b1 ]; then
  echo bucket file for 1st filtering = $b1 does not exist.
  exit
fi

sm="$sm_dir/${pivot1}_${pivot2}_${range}_${qbit}-bit.sm"

w1=$(./src/pivot_property.sh -w $p1)
echo PJT_DIM = $w1
dm=$(./src/pivot_property.sh -d $p1)
echo FTR_DIM = $dm
pt1=$(./src/pivot_property.sh -p $p1)
np1=$(./src/pivot_property.sh -n $p1)
echo partition type of 1st pivot = $pt1, number of partitioned spaces = $np1

w2=$($pr_dir/pivot_property.sh -w $p2)
echo SMAP_DIM = $w2
d2=$($pr_dir/pivot_property.sh -d $p2)
pt=$($pr_dir/pivot_property.sh -p $p2)
np=$($pr_dir/pivot_property.sh -n $p2)
echo partition type of 2nd pivot = $pt, number of partitioned spaces = $np

if [ $dm -ne $d2 ] ; then
  echo ftr dimensions of 1st and 2nd pivots NOT EQUAL.
  exit 
fi

if [ $nt == 1 ] ; then
    cflags="-O3 -Wall -Wno-strict-overflow"
else
    cflags="-O3 -fopenmp -Wall -Wno-strict-overflow"
fi
cflags="$cflags -DNUM_THREADS=$nt"
cflags="$cflags -DMEMORY_LIMIT=115e9"       # MAX_RAM で決めた方がいい 

cflags="$cflags -DCOMPILE_TIME"				# プログラム編集時に利用する parm.h の定義を無効にするスイッチ
cflags="$cflags -D$dataset"
cflags="$cflags -DFTR_DIM=$dm"
cflags="$cflags -DPJT_DIM=$w1"
cflags="$cflags -DSMAP_DIM=$w2"

cflags="$cflags -DSEQUENTIAL_FILTERING"
cflags="$cflags -DFTR_ON_SECONDARY_MEMORY"
cflags="$cflags -DFTR_ARRANGEMENT_ASIS"

if [ $qbit == 3 ] ; then
    cflags="$cflags -DQUANTIZE_BIT=$qbit"
    cflags="$cflags -DUSE_PACKED_3BIT"
elif [ $qbit == 6 ] ; then
    cflags="$cflags -DQUANTIZE_BIT=3"
    cflags="$cflags -DUSE_PACKED_6BIT"
else
    cflags="$cflags -DQUANTIZE_BIT=$qbit"
    cflags="$cflags -DTINY_IN_CHAR"
fi	

cflags="$cflags -DQUANTIZE_RANGE=$qrange"

if [ $pt1 == 3 ] ; then
	cflags="$cflags -DPARTITION_TYPE_PQBP"
	cflags="$cflags -DNUM_PART=$np1"
else
	cflags="$cflags -DPARTITION_TYPE_QBP"
fi

if [ $pt == 3 ] ; then
	cflags="$cflags -DSMAP_PARTITION_TYPE_PQBP"
	if [ $use_pd -eq 0 ] ; then
		cflags="$cflags -DSMAP_NUM_PART=$np"
	else
		cflags="$cflags -DSMAP_NUM_PART=1"
		cflags="$cflags -DUSE_PD"
		cflags="$cflags -DIGNORE_MED"
	fi
else
	cflags="$cflags -DSMAP_PARTITION_TYPE_QBP"
fi

cflags="$cflags -DBLOCK_SIZE=800"

cflags="$cflags -DPIVOT_FILE=\"$p1\""
cflags="$cflags -DBUCKET_FILE=\"$b1\""
cflags="$cflags -DSMAP_PIVOT_FILE=\"$p2\""
cflags="$cflags -DQPSMAP_FILE=\"$sm\""

echo $cflags

ds=$($pr_dir/expand_filenames.sh $prefix $range .ftr)

files=""
for f in $ds ; do
	if [ ! -e $ds_dir/$f ]; then
	    echo dataset file = $ds_dir/$f does not exist
	    exit 1
	fi
	files="$files $ds_dir/$f"
done

lb_list="bit_op ftr e_time sketch quick smap"
lb=""
for l in $lb_list ; do
  lb="$lb $pr_dir/$l.c"
done
echo library = $lb

gcc $cflags $pr_dir/$pr.c $lb -o $pr -lm

if [ $? == 1 ] ; then
exit 1
fi

time ./$pr $files

