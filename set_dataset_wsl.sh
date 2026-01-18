#!/bin/bash
if [ $# -ne 1 ] ; then
echo "Usage> ./set_dataset_wsl.sh <dataset>"
exit 1
fi

src=./src
script=./script

dataset=$1; shift

if [ $dataset == "PUBMED23" ] ; then
ds_dir=/mnt/ssd/PUBMED23/dataset
qr_dir=/mnt/ssd/PUBMED23/query
pv_dir=/mnt/ssd/PUBMED23/pivot
bk_dir=/mnt/ssd/PUBMED23/bkt
bt_dir=/mnt/u/PUBMED23/batch
sm_dir=/mnt/ssd/PUBMED23/smap
elif [ $dataset == "DECAF" ] || [ $dataset == "YFCC100M" ]; then
ds_dir=/mnt/ssd/DISA_h/dataset
qr_dir=/mnt/ssd/DISA_h/query
pv_dir=/mnt/ssd/DISA_h/pivot
bk_dir=/mnt/ssd/DISA_h/bkt
bt_dir=/mnt/u/DISA_h/batch
sm_dir=/mnt/ssd/DISA_h/smap
elif [ $dataset == "LAION2B" ] || [ $dataset == "LAION100M" ]; then
ds_dir=/mnt/ssd/LAION2B/dataset
qr_dir=/mnt/ssd/LAION2B/query
pv_dir=/mnt/ssd/LAION2B/pivot
bk_dir=/mnt/ssd/LAION2B/bkt
bt_dir=/mnt/u/LAION2B/batch
sm_dir=/mnt/ssd/LAION2B/smap
elif [ $dataset == "DEEP1B" ] ; then
ds_dir=/mnt/ssd/Deep1B/dataset
qr_dir=/mnt/ssd/Deep1B/query
pv_dir=/mnt/ssd/Deep1B/pivot
bk_dir=/mnt/ssd/Deep1B/bkt
bt_dir=/mnt/u/Deep1B/batch
sm_dir=/mnt/ssd/Deep1B/smap
else
echo "Invalid dataset: $dataset"
exit 1
fi

export DATASET=$dataset
export SRC=$src
export SCR=$script
export FTR=$ds_dir
export QUR=$qr_dir
export PIV=$pv_dir
export BKT=$bk_dir
export SMA=$sm_dir

echo DATASET = $DATASET
echo SRC     = $SRC
echo SCR     = $SCR
echo FTR     = $FTR
echo QUR     = $QUR
echo PIV     = $PIV
echo BKT     = $BKT
