#!/bin/bash
if [ $# -ne 3 ] ; then
echo "Usage> ./docker_dun.sh <dataset> <ram> <cpu>"
echo "dataset = PUBMED23, DECAF (YFCC100M), LAION2B (LAION100M), or DEEP1B"
echo "ram     = RAM in GB"
echo "cpu     = number of CPU cores"
exit 1
fi

src=./src
script=./script

dataset=$1; shift
ram=$1; shift
cpu=$1; shift

rs_dir=./${dataset}_result

if [ ! -e $rs_dir ]; then
  echo result directory \"$rs_dir\" dose not exist. mkdir $rs_dir
  mkdir $rs_dir
fi

echo $dataset

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
dataset=LAION2B
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

docker run -it --memory=${ram}g --memory-swap=${ram}g --cpus=${cpu} --name=df_${dataset}_${ram}G_${cpu} --hostname=${dataset}_${ram}G_${cpu} \
-v $src:/app/src \
-v $script:/app/script \
-v $ds_dir:/app/ftr \
-v $qr_dir:/app/query \
-v $pv_dir:/app/pivot \
-v $bk_dir:/app/bkt \
-v $bt_dir:/app/batch \
-v $sm_dir:/app/smap \
-v $rs_dir:/app/result \
-v /etc/localtime:/etc/localtime:ro \
-v /etc/timezone:/etc/timezone:ro \
-e DATASET=${dataset} \
-e MAX_RAM=${ram} \
-e NUM_THREADS=${cpu} \
df_search /bin/bash
