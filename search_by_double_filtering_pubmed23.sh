#!/bin/bash
if [ $# -ne 25 ] ; then
echo "Usage>"
echo "1.  1st pivot  = pivot file for narrow sketch (for 1st filtering)"
echo "2.  2nd pivot  = pivot file for qpsmap (for 2nd filtering)"
echo "3.  range      = range of ftr files" 
echo "4.  1st query  = file number of 1st query file" 
echo "5.  2nd query  = file number of 2nd query file (-1 for skip)" 
echo "6.  3rd query  = file number of 3rd query file (-1 for skip)" 
echo "7.  result.csv = result file" 
echo "8.  num_k      = k (number of neighbors), -2 (compute recall without search)"
echo "9.  num_q      = number of queries (0 -> all)"
echo "10. q_bit      = quantize bit"
echo "11. q_range    = quantize range"
echo "12. plus_half  = plus in recover (PLUS_HALF), or -1 -> use lower bound"
echo "13. p_1st      = 10 * p of D~p for 1st filtering"
echo "14. p_2nd      = 10 * p of D~p for 2nd filtering"
echo "15. #threads   = number of threads"
echo "16. enum_bit   = number of bits to enumerate first"
echo "17. supp_bit   = number of supplementary bits to enumarete"
echo "18. interval   = 1 -> using interval, 0 -> without, 2 -> D~1, 3 -> D~1 using interval, 4 -> D~1 using interval by single-thread"
echo "19. factor_inf = FACTOR_INF"
echo "20. factor2    = FACTOR_INF2"
echo "21. thd_plus   = THREAD_PLUS"
echo "22. #threads_p = number of threads for preparation before search (0 -> same as #threads)"
echo "23. FTR_ON     = 0 -> Second memory, 1 -> RAM"
echo "24. use pd     = (0: without, 1: use pd)"
echo "25. batch file = file name including hyper parameters (nc1 and nc2), or NONE for interactive manner"

exit 1
fi

dataset=PUBMED23
pr_dir=./src
ds_dir=./ftr
qr_dir=./query
pv_dir=./pivot
bk_dir=./bkt

pr=search_by_double_filtering_hamming_smap_v5

#set -x

ulimit -Sn 10000
#ulimit -a

pivot1=$1; shift

p1="$pv_dir/${pivot1}.csv"
if [ ! -e $p1 ]; then
  echo pivot file for 1st filtering = $p1 does not exist.
  exit
fi

pivot2=$1; shift

p2="$pv_dir/${pivot2}.csv"
if [ ! -e $p2 ]; then
  echo pivot file of qpsmap for 2nd filtering = $p2 does not exist.
  exit
fi

range=$1; shift
b1="$bk_dir/${pivot1}_${range}.bkt"
if [ ! -e $b1 ]; then
  echo bucket file for 1st filtering = $b1 does not exist.
  exit
fi

w1=$(./src/pivot_property.sh -w $p1)
echo PJT_DIM = $w1
dm=$(./src/pivot_property.sh -d $p1)
echo FTR_DIM = $dm
pt1=$(./src/pivot_property.sh -p $p1)
np1=$(./src/pivot_property.sh -n $p1)
echo partition type of 1st pivot = $pt1, number of partitioned spaces = $np1

w2=$(./src/pivot_property.sh -w $p2)
echo SMAP_DIM = $w2
d2=$(./src/pivot_property.sh -d $p2)
pt=$(./src/pivot_property.sh -p $p2)
np=$(./src/pivot_property.sh -n $p2)
echo partition type of 2nd pivot = $pt, number of partitioned spaces = $np

if [ $dm -ne $d2 ] ; then
  echo ftr dimensions of 1st and 2nd pivots NOT EQUAL.
  exit 
fi

query=$1; shift
qr_prefix=query_
qr="$qr_dir/${qr_prefix}${query}.ftr"
if [ ! -e $qr ]; then
  echo query file = $qr does not exist.
  exit
fi

an="$qr_dir/answer_q_${query}_r_${range}.csv"
if [ ! -e $an ]; then
  echo answer file = $an does not exist.
  exit
fi

q2nd=$1; shift
if [ $q2nd == NONE ] ; then
  q2=NONE
  a2=NONE
else
  q2="$qr_dir/${qr_prefix}${q2nd}.ftr"
  if [ ! -e $q2 ]; then
    echo query file = $q2 does not exist.
    exit
  fi
  a2="$qr_dir/answer_q_${q2nd}_r_${range}.csv"
  if [ ! -e $a2 ]; then
    echo answer file = $a2 does not exist.
    exit
  fi
fi

q3rd=$1; shift
if [ $q3rd == NONE ] ; then
  q3=NONE
  a3=NONE
else
  q3="$qr_dir/${qr_prefix}${q3rd}.ftr"
  if [ ! -e $q3 ]; then
    echo query file = $q3 does not exist.
    exit
  fi
  a3="$qr_dir/answer_q_${q3rd}_r_${range}.csv"
  if [ ! -e $a3 ]; then
    echo answer file = $a3 does not exist.
    exit
  fi
fi

rs="$1"; shift

nk="$1"; shift
nq="$1"; shift
qbit=$1; shift
qrange=$1; shift
phalf=$1; shift
sp1="$1"; shift
sp2="$1"; shift
nt="$1"; shift
enum=$1; shift
supp=$1; shift
interval=$1; shift
sf=$1; shift
sf2=$1; shift
tp=$1; shift
n0=$1; shift
fo=$1; shift
use_pd=$1; shift
batch=$1; shift

if [ $batch != NONE ] ; then
if [ ! -e $batch ] ; then
echo batch file $batch does not exist.
exit
fi
fi

if [ $interval -le 1 ] ; then
ft="FILTERING_BY_SKETCH_ENUMERATION_HAMMING"
elif [ $interval -le 4 ] ; then
ft="FILTERING_BY_SKETCH_ENUMERATION_C2N"
fi

if [ $fo == 0 ] ; then
fo="SECONDARY_MEMORY"
else
fo="MAIN_MEMORY"
fi

fr="ARRANGEMENT_ASIS"

if [ $nt == 1 ] ; then
	if [ $n0 -lt 1 ] ; then
    cflags="-O3 -Wall -Wno-strict-overflow"
    else
    cflags="-O3 -fopenmp -Wall -Wno-strict-overflow"
    cflags="$cflags -DNUM_THREADS_PREPERATION=$n0"
    fi
    cflags="$cflags -DPARALLEL_ENUM=0"
    cflags="$cflags -DPARA_ENUM_INF=0"
else
    cflags="-O3 -fopenmp -Wall -Wno-strict-overflow"
    if [ $n0 == 0 ] ; then
    cflags="$cflags -DNUM_THREADS_PREPERATION=$nt"
    else
    cflags="$cflags -DNUM_THREADS_PREPERATION=$n0"
    fi
    if [ $interval -ne 4 ] ; then
	    if [ $nt == 2 ] ; then
	        cflags="$cflags -DPARALLEL_ENUM=1"
	        cflags="$cflags -DPARA_ENUM_INF=1"
	    elif [ $nt == 4 ] ; then
	        cflags="$cflags -DPARALLEL_ENUM=2"
	        cflags="$cflags -DPARA_ENUM_INF=2"
	    elif [ $nt == 8 ] ; then
	        cflags="$cflags -DPARALLEL_ENUM=3"
	        cflags="$cflags -DPARA_ENUM_INF=3"
	    elif [ $nt == 16 ] ; then
	        cflags="$cflags -DPARALLEL_ENUM=4"
	        cflags="$cflags -DPARA_ENUM_INF=4"
	    elif [ $nt == 32 ] ; then
	        cflags="$cflags -DPARALLEL_ENUM=5"
	        cflags="$cflags -DPARA_ENUM_INF=5"
	    else
	    echo "Number of threads should be 1, 2, 4, 8 or 16."
	    exit
	    fi
	else
	    if [ $nt == 2 ] ; then
	        cflags="$cflags -DPARALLEL_ENUM=1"
	        cflags="$cflags -DPARA_ENUM_INF=0"
	    elif [ $nt == 4 ] ; then
	        cflags="$cflags -DPARALLEL_ENUM=2"
	        cflags="$cflags -DPARA_ENUM_INF=0"
	    elif [ $nt == 8 ] ; then
	        cflags="$cflags -DPARALLEL_ENUM=3"
	        cflags="$cflags -DPARA_ENUM_INF=0"
	    elif [ $nt == 16 ] ; then
	        cflags="$cflags -DPARALLEL_ENUM=4"
	        cflags="$cflags -DPARA_ENUM_INF=0"
	    elif [ $nt == 32 ] ; then
	        cflags="$cflags -DPARALLEL_ENUM=5"
	        cflags="$cflags -DPARA_ENUM_INF=0"
	    else
	    echo "Number of threads should be 1, 2, 4, 8 or 16."
	    exit
	    fi
	fi
fi
cflags="$cflags -DNUM_THREADS=$nt"
cflags="$cflags -DMEMORY_LIMIT=115e9"

cflags="$cflags -DCOMPILE_TIME"
cflags="$cflags -D$dataset"
cflags="$cflags -DFTR_DIM=$dm"
cflags="$cflags -DPJT_DIM=$w1"
cflags="$cflags -DNARROW_SKETCH"
cflags="$cflags -DSMAP_DIM=$w2"
cflags="$cflags -DENUM_DIM=$enum"
cflags="$cflags -DSPP_BIT=$supp"

cflags="$cflags -DSECOND_FILTERING_KNN_BUFFER"
#cflags="$cflags -DSECOND_FILTERING_SELECT"
#cflags="$cflags -DMIN_LIST=$ml"

cflags="$cflags -DFACTOR_INF=$sf"
cflags="$cflags -DFACTOR_INF2=$sf2"

if [ $interval -eq 1 ] ; then
    cflags="$cflags -DWITH_ENUM_TABLE"
    cflags="$cflags -DSELECT_BY_PRIORITY_AFTER_ENUMERATION"
    cflags="$cflags -DUSE_INTERVAL"
    cflags="$cflags -DINTERVAL_WITH_RUN"
    cflags="$cflags -DINTERVAL_WITH_PRIORITY"
    cflags="$cflags -DLOOP_CONTROL_BY_NUM_SKETCHES"
    if [ $nt -ne 1 ] ; then
        cflags="$cflags -DSELECTION_BY_MULTI_THREAD"
    fi
#    cflags="$cflags -DFACTOR_INF=$sf"
#    cflags="$cflags -DFACTOR_INF2=$sf2"
elif [ $interval -eq 2 ] ; then
    cflags="$cflags -DSELECT_BY_PRIORITY_AFTER_ENUMERATION"
#    cflags="$cflags -DUSE_INTERVAL"
#    cflags="$cflags -DINTERVAL_WITH_RUN"
#    cflags="$cflags -DINTERVAL_WITH_PRIORITY"
#    cflags="$cflags -DUSE_MU_COMMON"
#    cflags="$cflags -DFACTOR_INF=$sf"
#    cflags="$cflags -DFACTOR_INF2=$sf2"
elif [ $interval -eq 3 ] ; then
#    cflags="$cflags -DALTERNATIVE_ENUMERATION"
#    cflags="$cflags -DALTERNATIVE_5"
    cflags="$cflags -DSELECT_BY_PRIORITY_AFTER_ENUMERATION"
    cflags="$cflags -DUSE_INTERVAL"
#    cflags="$cflags -DUSE_COMPACT_INTERVAL"
    cflags="$cflags -DINTERVAL_WITH_RUN"
#    cflags="$cflags -DINTERVAL_WITH_PRIORITY"
    cflags="$cflags -DUSE_MU_COMMON"
#    cflags="$cflags -DFACTOR_INF=$sf"
#    cflags="$cflags -DFACTOR_INF2=$sf2"
elif [ $interval -eq 4 ] ; then
#   cflags="$cflags -DSELECT_BY_PRIORITY_AFTER_ENUMERATION"
    cflags="$cflags -DUSE_INTERVAL"
    cflags="$cflags -DINTERVAL_WITH_RUN"
    cflags="$cflags -DINTERVAL_WITH_PRIORITY"
    cflags="$cflags -DUSE_MU_COMMON"
#    cflags="$cflags -DFACTOR_INF=$sf"
#    cflags="$cflags -DFACTOR_INF2=$sf2"
fi

if [ $qbit -lt 8 ] ; then
	if [ $qbit == 3 ] ; then
		cflags="$cflags -DQUANTIZE_BIT=$qbit"
		cflags="$cflags -DUSE_PACKED_3BIT"
	elif [ $qbit == 6 ] ; then
		cflags="$cflags -DQUANTIZE_BIT=3"
		cflags="$cflags -DUSE_PACKED_6BIT"
	else
		cflags="$cflags -DQUANTIZE_BIT=$qbit"
# cflags="$cflags -DUSE_PACKED_4BIT"
	  cflags="$cflags -DTINY_IN_CHAR"
	fi	
else
	if [ $qbit == 600 ] ; then
	cflags="$cflags -DUSE_PACKED_6BIT"
	qbit=332
	fi
	cflags="$cflags -DQUANTIZE_MIXED_MOD3"
	qb0=$(($qbit / 100))
	qb1=$((($qbit - $qb0 * 100) / 10))
	qb2=$(($qbit - $qb0 * 100 - $qb1 * 10))
	if [ $(($qb0 + $qb1 + $qb2)) -ne 8 ] ; then
	echo invalid quantize mixed bits: sum of bits should be 8.
	exit
	fi
	cflags="$cflags -DQUANTIZE_BIT_0=$qb0"
	cflags="$cflags -DQUANTIZE_BIT_1=$qb1"
	cflags="$cflags -DQUANTIZE_BIT_2=$qb2"
	qbmax=$(($qb0 > $qb1 ? $qb0 : $qb1))
	qbmax=$(($qbmax > $qb2 ? $qbmax : $qb2))
	cflags="$cflags -DQUANTIZE_BIT=$qbmax"
	cflags="$cflags -DTINY_IN_CHAR"
fi

cflags="$cflags -DQUANTIZE_RANGE=$qrange"
#if [ $phalf -lt 0 ] ; then
#cflags="$cflags -DUSE_LOWER_BOUND"
#else
cflags="$cflags -DPLUS_HALF=$phalf"
#fi
cflags="$cflags -DPRIORITY=1"
cflags="$cflags -DSCORE_P_1ST=${sp1}"
cflags="$cflags -DSCORE_P_2ND=${sp2}"
#cflags="$cflags -DTINY_IN_CHAR"

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

#cflags="$cflags -DSMAP_NUM_PART=$np"
else
cflags="$cflags -DSMAP_PARTITION_TYPE_QBP"
fi

cflags="$cflags -DFTR_ON_${fo}"
cflags="$cflags -DFTR_${fr}"

cflags="$cflags -DBLOCK_SIZE=800"
cflags="$cflags -D${ft}"

cflags="$cflags -DPIVOT_FILE=\"$p1\""
cflags="$cflags -DBUCKET_FILE=\"$b1\""
cflags="$cflags -DSMAP_PIVOT_FILE=\"$p2\""
cflags="$cflags -DQUERY_FILE=\"$qr\""
cflags="$cflags -DANSWER_FILE=\"$an\""
cflags="$cflags -DRESULT_FILE=\"result/$rs\""
cflags="$cflags -DQUERY_2ND_FILE=\"$q2\""
cflags="$cflags -DANSWER_2ND_FILE=\"$a2\""
cflags="$cflags -DQUERY_3RD_FILE=\"$q3\""
cflags="$cflags -DANSWER_3RD_FILE=\"$a3\""
if [ $batch != NONE ] ; then
cflags="$cflags -DINPUT_HYPER_PARAMETER=\"$batch\""
fi
cflags="$cflags -DANSWER_DIST_FLOAT"

cflags="$cflags -DNUM_K=$nk"

if [ $tp -gt 0 ] ; then
cflags="$cflags -DTHREAD_PLUS=$tp"
fi

#cflags="$cflags -DSELECT_K_BY_QUICK_SELECT_K"
#cflags="$cflags -DSELECT_K_BY_QUICK_SELECT_K_PARA_WORK"
#cflags="$cflags -DSELECT_K_BY_QUICK_SORT"


if [ $nq > 0 ] ; then
cflags="$cflags -DNUM_Q=$nq"
fi

echo $cflags

prefix=pubmed23_

ds=$(./src/expand_filenames.sh $prefix $range .ftr)

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

